#!/usr/bin/env python3
"""PAR post-step for CNVkit genemetrics in male samples (PLAN.md §14.2a, "After B6", item 1).

CNVkit limitation (fc65941; `CopyNumArray.shift_xx`, unchanged upstream at etal/cnvkit b4509b4): for a sample that
genemetrics treats as male, `shift_xx` adds +1.0 to the log2 of EVERY chrX bin and segment
(`chromosome == chrX`), including the pseudo-autosomal regions, even with `--diploid-parx-genome grch38`. The PARs
are diploid in a male, so every PAR gene (on this panel: CRLF2, PAR1) is reported 1.0 too high.

Post-step (exact by construction):
  1. A call is "treated as male" when genemetrics was given `-x male` (plan column sample_sex = male), or when it
     was given no -x, the arm uses --diploid-parx-genome (s1/s2/s1x arms) and `cnvkit.py sex --diploid-parx-genome
     grch38` on the .cnr reports Male (the same guess_xx call genemetrics makes). s0 arms (no PARX) are skipped.
  2. Copy <sample>.cnr and <sample>.call.cns; subtract 1.0 from log2 of every bin / segment lying wholly inside
     GRCh38 chrX PAR1 (10,001-2,781,479) or PAR2 (155,701,383-156,030,895) (0-based half-open [10000, 2781479),
     [155701382, 156030895)). Bins or segments that overlap a PAR boundary are reported (`span_*` columns) and
     left unchanged: where CBS merged the PAR1 bins into a segment that continues into non-PAR chrX, the PAR gene
     keeps the value genemetrics gives that (mostly haploid non-PAR) segment, +1 included (log: par_check = span).
  3. Re-run genemetrics in the image with the original arguments:
       genemetrics <s>.cnr -s <s>.call.cns -t 0.3 --min-probes 3 -o <s>.genemetrics.tsv --diploid-parx-genome grch38 [-x male]
     (in-process `cnvlib.commands.parse_args(argv).func`, the code path of cnvkit.py).
  4. The original is kept as <sample>.genemetrics.noparfix.tsv; the result replaces <sample>.genemetrics.tsv.
  5. Idempotent: a call with <sample>.genemetrics.noparfix.tsv and <sample>.genemetrics.tsv present is skipped.
Every male call is processed and logged (also when no row changes). The log records, per call: sex source, the
`cnvkit.py sex` result, the observed genemetrics chrX shift (gm log2 - containing call.cns segment log2 on non-PAR
chrX rows; 1 = male, 0 = female), PAR bins/segments changed, rows changed, and whether every non-PAR row is
byte-identical (asserted).

Usage
  genemetrics_parfix.py apply  --log LOG.tsv --plans P.tsv [...] --arm-dirs DIR [...] [--jobs 10]
                               [--samples-file F] [--dry-run OUTDIR] [--rerun-check]
      --arm-dirs     <calls>/<arm> directories (symlinks resolved; duplicates collapsed). Arms whose name ends in
                     _s0 are refused.
      --samples-file restrict to these .cnr paths (validation).
      --dry-run DIR  write the fixed genemetrics to DIR/<run dir key>/ instead of replacing (validation);
                     --rerun-check additionally re-runs genemetrics on the UNMODIFIED inputs and asserts the result
                     equals the original file byte-for-byte (argument check), also for females.
  genemetrics_parfix.py check-s0 --arm-dir <calls>/<arm_s0> --n 10
      verifies that s0 genemetrics applies no chrX shift to LIMS-male samples (gm - segment = 0).
"""
import argparse, collections, concurrent.futures as cf, csv, json, os, subprocess, sys, tempfile

IMG = 'cgp-cnvkit:1.0.0'
W = '/home/wook/Documents/pon_comparison_work'
PARS = [(10000, 2781479), (155701382, 156030895)]   # GRCh38 chrX PAR1, PAR2; 0-based half-open
GM_ARGS = ['-t', '0.3', '--min-probes', '3']
PARX = ['--diploid-parx-genome', 'grch38']


def par_state(chrom, start, end):
    """'in' (wholly inside a PAR), 'span' (overlaps a PAR boundary) or 'out'."""
    if chrom != 'chrX':
        return 'out'
    for a, b in PARS:
        if start >= a and end <= b:
            return 'in'
        if start < b and end > a:
            return 'span'
    return 'out'


# ------------------------------------------------------------------ worker (runs inside the image)
WORKER = r'''
import csv, io, json, os, shutil, subprocess, sys, tempfile
PARS = [(10000, 2781479), (155701382, 156030895)]
def par_state(chrom, start, end):
    if chrom != 'chrX':
        return 'out'
    for a, b in PARS:
        if start >= a and end <= b:
            return 'in'
        if start < b and end > a:
            return 'span'
    return 'out'
from cnvlib import commands
import logging; logging.basicConfig(level=logging.WARNING, format="%(message)s")
def genemetrics(cnr, cns, out, sex):
    argv = ['genemetrics', cnr, '-s', cns, '-t', '0.3', '--min-probes', '3', '-o', out, '--diploid-parx-genome', 'grch38']
    if sex:
        argv += ['-x', sex]
    a = commands.parse_args(argv); a.func(a)
def shift_par(src, dst, start_col, log2_col):
    """copy a .cnr/.call.cns, log2 - 1 for rows wholly inside a PAR; returns (n_in, n_span)."""
    n_in = n_span = 0
    with open(src) as f, open(dst, 'w') as g:
        g.write(f.readline())
        for l in f:
            c = l.rstrip('\n').split('\t')
            st = par_state(c[0], int(c[1]), int(c[2]))
            if st == 'in':
                c[log2_col] = repr(float(c[log2_col]) - 1.0); n_in += 1
                g.write('\t'.join(c) + '\n')
            else:
                n_span += st == 'span'
                g.write(l)
    return n_in, n_span
job = json.load(open(sys.argv[1]))
rd = job['rundir']; out = []
names = [s['name'] for s in job['samples']]
inferred = {}
need = [n for n in names if True]
if need:
    tmp = tempfile.mkdtemp()
    subprocess.run(['cnvkit.py', 'sex'] + [os.path.join(rd, n + '.cnr') for n in need] + ['-o', tmp + '/sex.tsv', '--diploid-parx-genome', 'grch38'],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for r in csv.DictReader(open(tmp + '/sex.tsv'), delimiter='\t'):
        inferred[os.path.basename(r['sample'])[:-4]] = r['sex']
    shutil.rmtree(tmp)
for s in job['samples']:
    n = s['name']; sex = s['sex']; base = os.path.join(rd, n)
    rec = dict(rundir=rd, sample=n, plan_sex=sex or '-', inferred=inferred.get(n, ''))
    male = sex == 'male' or (sex == '' and inferred.get(n) == 'Male')
    rec['treated_male'] = int(male)
    gm, nop = base + '.genemetrics.tsv', base + '.genemetrics.noparfix.tsv'
    dry = job.get('dry')
    if job.get('rerun_check'):
        tmp = tempfile.mkdtemp()
        shutil.copy(base + '.cnr', os.path.join(tmp, n + '.cnr')); shutil.copy(base + '.call.cns', os.path.join(tmp, n + '.call.cns'))
        cwd = os.getcwd(); os.chdir(tmp)
        genemetrics(n + '.cnr', n + '.call.cns', n + '.genemetrics.tsv', sex)
        os.chdir(cwd)
        orig = open(nop if os.path.exists(nop) else gm).read()
        rec['rerun_identical'] = int(open(os.path.join(tmp, n + '.genemetrics.tsv')).read() == orig)
        shutil.rmtree(tmp)
    if not male:
        rec['action'] = 'not_male'; out.append(rec); continue
    if not dry and os.path.exists(nop) and os.path.exists(gm):
        rec['action'] = 'skipped_done'; out.append(rec); continue
    orig_path = nop if os.path.exists(nop) else gm
    tmp = tempfile.mkdtemp()
    nb, nbs = shift_par(base + '.cnr', os.path.join(tmp, n + '.cnr'), 1, 4)
    ns, nss = shift_par(base + '.call.cns', os.path.join(tmp, n + '.call.cns'), 1, 4)
    cwd = os.getcwd(); os.chdir(tmp)
    genemetrics(n + '.cnr', n + '.call.cns', n + '.genemetrics.tsv', sex)
    os.chdir(cwd)
    rec.update(par_bins=nb, span_bins=nbs, par_segs=ns, span_segs=nss)
    new_path = os.path.join(tmp, n + '.genemetrics.tsv')
    if dry:
        os.makedirs(dry, exist_ok=True); shutil.copy(new_path, os.path.join(dry, n + '.genemetrics.tsv'))
        rec['action'] = 'dry'
    else:
        shutil.copy(new_path, gm + '.parfix.tmp')
        if not os.path.exists(nop):
            os.replace(gm, nop)
        os.replace(gm + '.parfix.tmp', gm)
        rec['action'] = 'fixed'
    rec['orig'] = orig_path; rec['new'] = os.path.join(dry, n + '.genemetrics.tsv') if dry else gm
    shutil.rmtree(tmp)
    out.append(rec)
json.dump(out, open(sys.argv[2], 'w'))
'''


def read_rows(p):
    with open(p) as f:
        hdr = f.readline()
        return hdr, [l for l in f]


def seg_lookup(cns):
    segs = collections.defaultdict(list)
    with open(cns) as f:
        next(f)
        for l in f:
            c = l.split('\t'); segs[c[0]].append((int(c[1]), int(c[2]), float(c[4])))
    return segs


def chrx_shift(gm_path, cns_path):
    """observed genemetrics - segment log2 on non-PAR chrX rows: set of rounded deltas."""
    segs = seg_lookup(cns_path); d = set()
    with open(gm_path) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if c[1] != 'chrX' or par_state(c[1], int(c[2]), int(c[3])) != 'out':
                continue
            s, e, lg = int(c[2]), int(c[3]), float(c[4])
            for ss, se, sl in segs['chrX']:
                if ss <= s < se or ss < e <= se:
                    d.add(round(lg - sl, 3)); break
    return d


def par_counts(path):
    n = collections.Counter()
    with open(path) as f:
        next(f)
        for l in f:
            c = l.split('\t'); n[par_state(c[0], int(c[1]), int(c[2]))] += 1
    return n['in'], n['span']


def compare(rec):
    """non-PAR rows byte-identical (asserted); PAR rows before/after. Fixed calls: original = .noparfix."""
    if rec.get('action') not in ('fixed', 'dry', 'skipped_done'):
        return rec
    base = os.path.join(rec['rundir'], rec['sample'])
    if rec['action'] != 'dry':
        rec['orig'] = base + '.genemetrics.noparfix.tsv'; rec['new'] = base + '.genemetrics.tsv'
    rec['par_bins'], rec['span_bins'] = par_counts(base + '.cnr')
    rec['par_segs'], rec['span_segs'] = par_counts(base + '.call.cns')
    h0, r0 = read_rows(rec['orig']); h1, r1 = read_rows(rec['new'])
    # CNVkit writes a 5-column header when genemetrics finds no gene (empty table); otherwise headers must match
    assert h0 == h1 or (not r1 and h1 == 'gene\tchromosome\tstart\tend\tlog2\n'), rec
    rec['empty_after'] = int(not r1)
    isp = lambda l: (lambda c: par_state(c[1], int(c[2]), int(c[3])) != 'out')(l.split('\t'))
    n0 = [l for l in r0 if not isp(l)]; n1 = [l for l in r1 if not isp(l)]
    p0 = [l for l in r0 if isp(l)]; p1 = [l for l in r1 if isp(l)]
    rec['nonpar_identical'] = int(n0 == n1)
    assert n0 == n1, ('non-PAR rows differ', rec['rundir'], rec['sample'])
    rec['par_rows_before'] = len(p0); rec['par_rows_after'] = len(p1)
    rec['rows_changed'] = len(collections.Counter(p0) - collections.Counter(p1)) + len(collections.Counter(p1) - collections.Counter(p0))
    rec['file_changed'] = int(r0 != r1)
    rec['par_before'] = ';'.join(sorted({f"{l.split(chr(9))[0]}:{l.split(chr(9))[4]}" for l in p0}))
    rec['par_after'] = ';'.join(sorted({f"{l.split(chr(9))[0]}:{l.split(chr(9))[4]}" for l in p1}))
    # CRLF2 (PAR) rows after the fix must equal the call.cns segment log2 (no +1)
    segs = seg_lookup(os.path.join(rec['rundir'], rec['sample'] + '.call.cns'))
    ok = True
    for l in p1:
        c = l.split('\t'); s, e = int(c[2]), int(c[3])
        sl = [x[2] for x in segs['chrX'] if x[0] <= s and e <= x[1]]
        ok &= bool(sl) and abs(float(c[4]) - sl[0]) < 5e-6 and abs(sl[0]) >= 0.3
    for l in p0:   # every PAR row that disappeared must be below threshold without the +1
        if l in p1:
            continue
        c = l.split('\t'); s, e = int(c[2]), int(c[3])
        sl = [x[2] for x in segs['chrX'] if x[0] <= s and e <= x[1]]
        ok &= bool(sl) and (abs(sl[0]) < 0.3 or any(m.split('\t')[0] == c[0] for m in p1))
    # a PAR gene whose segment spans the PAR boundary is left as genemetrics reported it (the segment is not
    # modified): par_check = 'span' when the original and fixed PAR rows are then identical
    rec['par_check'] = int(ok) if not rec['span_segs'] else ('span' if p0 == p1 else int(ok))
    return rec


def load_sex(plans):
    sex = {}
    for p in plans:
        for r in csv.DictReader(open(p), delimiter='\t'):
            k = (os.path.realpath(r['outdir']), os.path.basename(r['sample_cnn'])[:-len('.targetcoverage.cnn')])
            v = '' if r['sample_sex'] in ('-', '') else r['sample_sex']
            if k in sex and sex[k] != v:
                sys.exit(f'conflicting sample_sex for {k}: {sex[k]} vs {v}')
            sex[k] = v
    return sex


def run_job(job, workdir):
    jf = os.path.join(workdir, 'job.json'); of = os.path.join(workdir, 'out.json'); wf = os.path.join(workdir, 'worker.py')
    json.dump(job, open(jf, 'w')); open(wf, 'w').write(WORKER)
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-v', f'{W}:{W}', '-v', f'{workdir}:{workdir}',
           '-e', 'HOME=/tmp', '-e', 'MPLCONFIGDIR=/tmp']
    if job.get('dry'):
        os.makedirs(job['dry'], exist_ok=True); cmd += ['-v', f"{job['dry']}:{job['dry']}"]
    cmd += [IMG, 'python3', wf, jf, of]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"{job['rundir']}: {p.stderr[-2000:]}")
    return json.load(open(of))


def apply(a):
    sex = load_sex(a.plans)
    restrict = {os.path.realpath(x.strip()) for x in open(a.samples_file)} if a.samples_file else None
    rundirs = collections.OrderedDict()
    for ad in a.arm_dirs:
        real = os.path.realpath(ad)
        arm = os.path.basename(real)
        if arm.endswith('_s0') or '_s' not in arm:
            sys.exit(f'{ad}: s0 arm (no PARX), refused')
        for dd in sorted(os.listdir(real)):
            for run in sorted(os.listdir(os.path.join(real, dd))):
                rd = os.path.join(real, dd, run)
                rundirs.setdefault(rd, None)
    jobs = []
    for rd in rundirs:
        samples = []
        for fn in sorted(os.listdir(rd)):
            if not fn.endswith('.cnr'):
                continue
            n = fn[:-4]
            if restrict is not None and os.path.join(rd, fn) not in restrict:
                continue
            if not os.path.exists(os.path.join(rd, n + '.genemetrics.tsv')) and not os.path.exists(os.path.join(rd, n + '.genemetrics.noparfix.tsv')):
                continue
            k = (rd, n)
            if k not in sex:
                sys.exit(f'no plan row for {k}')
            samples.append(dict(name=n, sex=sex[k]))
        if samples:
            key = rd[len(W) + 1:].replace('/', '__')
            jobs.append(dict(rundir=rd, samples=samples, rerun_check=a.rerun_check,
                             dry=os.path.join(os.path.realpath(a.dry_run), key) if a.dry_run else None))
    print(f'{len(jobs)} run dirs, {sum(len(j["samples"]) for j in jobs)} calls', file=sys.stderr)
    base = tempfile.mkdtemp(prefix='parfix_')
    recs = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = {}
        for i, j in enumerate(jobs):
            wd = os.path.join(base, f'j{i:04d}'); os.makedirs(wd)
            futs[ex.submit(run_job, j, wd)] = j
        for k, f in enumerate(cf.as_completed(futs)):
            recs.extend(f.result())
            if k % 25 == 0:
                print(f'  {k + 1}/{len(jobs)} run dirs done', file=sys.stderr)
    for r in recs:
        compare(r)
        gm_now = r['orig'] if r.get('action') in ('fixed', 'dry') else os.path.join(r['rundir'], r['sample'] + '.genemetrics.tsv')
        if r.get('action') == 'skipped_done':
            gm_now = os.path.join(r['rundir'], r['sample'] + '.genemetrics.noparfix.tsv')
        r['orig_chrx_shift'] = ','.join(map(str, sorted(chrx_shift(gm_now, os.path.join(r['rundir'], r['sample'] + '.call.cns')))))
        r['arm'] = r['rundir'].split(os.sep)[-3]
    cols = ['arm', 'rundir', 'sample', 'plan_sex', 'inferred', 'treated_male', 'orig_chrx_shift', 'action', 'rerun_identical',
            'par_bins', 'span_bins', 'par_segs', 'span_segs', 'par_rows_before', 'par_rows_after', 'rows_changed',
            'file_changed', 'empty_after', 'nonpar_identical', 'par_check', 'par_before', 'par_after']
    recs.sort(key=lambda r: (r['rundir'], r['sample']))
    with open(a.log, 'w') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore', restval=''); w.writeheader(); w.writerows(recs)
    # consistency of the male determination with the observed original genemetrics chrX shift
    bad = [r for r in recs if r['orig_chrx_shift'] and r['orig_chrx_shift'] != ('1.0' if r['treated_male'] else '0.0')]
    summ = collections.Counter()
    for r in recs:
        summ['calls'] += 1; summ['male'] += r['treated_male']; summ[r['action']] += 1
        summ['files_changed'] += int(r.get('file_changed') or 0); summ['rows_changed'] += int(r.get('rows_changed') or 0)
        summ['span_bins'] += int(r.get('span_bins') or 0); summ['calls_with_span_seg'] += int(r.get('span_segs') or 0) > 0
        summ['par_check_fail'] += r.get('par_check', 1) == 0
        summ['span_par_rows_left'] += r.get('par_check') == 'span' and int(r.get('par_rows_after') or 0) > 0
        summ['rerun_not_identical'] += r.get('rerun_identical', 1) == 0
    summ['sex_vs_observed_shift_mismatch'] = len(bad)
    for k, v in sorted(summ.items()):
        print(f'{k}\t{v}')
    for r in bad[:20]:
        print('MISMATCH', r['rundir'], r['sample'], r['plan_sex'], r['inferred'], r['orig_chrx_shift'], file=sys.stderr)


def check_s0(a):
    arm = os.path.realpath(a.arm_dir); n = 0; out = collections.Counter()
    for dd in sorted(os.listdir(arm)):
        for run in sorted(os.listdir(os.path.join(arm, dd))):
            rd = os.path.join(arm, dd, run)
            for fn in sorted(os.listdir(rd)):
                if not fn.endswith('.genemetrics.tsv') or fn.split('-')[4] != 'M':
                    continue
                s = fn[:-len('.genemetrics.tsv')]
                d = chrx_shift(os.path.join(rd, fn), os.path.join(rd, s + '.call.cns'))
                if not d:
                    continue
                # also PAR rows
                segs = seg_lookup(os.path.join(rd, s + '.call.cns'))
                par = []
                for l in list(open(os.path.join(rd, fn)))[1:]:
                    c = l.split('\t')
                    if par_state(c[1], int(c[2]), int(c[3])) == 'in':
                        sl = [x[2] for x in segs['chrX'] if x[0] <= int(c[2]) and int(c[3]) <= x[1]]
                        par.append(round(float(c[4]) - sl[0], 3))
                print(run, s, 'nonPAR chrX shift', sorted(d), 'PAR shift', par)
                out[tuple(sorted(d))] += 1; n += 1
                if n >= a.n:
                    print(dict(out)); return
    print(dict(out))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('apply')
    p.add_argument('--log', required=True); p.add_argument('--plans', nargs='+', required=True)
    p.add_argument('--arm-dirs', nargs='+', required=True); p.add_argument('--jobs', type=int, default=10)
    p.add_argument('--samples-file'); p.add_argument('--dry-run'); p.add_argument('--rerun-check', action='store_true')
    p = sp.add_parser('check-s0'); p.add_argument('--arm-dir', required=True); p.add_argument('--n', type=int, default=10)
    a = ap.parse_args()
    apply(a) if a.cmd == 'apply' else check_s0(a)
