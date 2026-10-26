#!/usr/bin/env python3
"""PAR post-step for CNVkit genemetrics in samples treated as male (DI-4101). Stdlib only.

Source: `analysis/pon_comparison/scripts/genemetrics_parfix.py` at commit f933f7a of cnv-backbone-purple-atlas
(project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6PpjQ4Yggx4JxXFkQqZxGb, md5 51ada0e597420655c89c86e1a2e13ea6). The batch
driver (plans, arm directories, the in-image worker and the hard-coded work directory) is replaced by three per-sample
subcommands; `par_state`, `shift_par`, the male decision and the `compare` checks are copied unchanged. The app runs
genemetrics in the image between `shift` and `compare` with the original arguments, as the validation did.
Validation: FINAL_RESULTS.md (project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6YBbj4YggvFK194V51V7pX) and DI-4118 (non-UMI).

CNVkit limitation (fc65941, `CopyNumArray.shift_xx`): for a sample that genemetrics treats as male, every chrX bin and
segment is shifted +1.0, the pseudo-autosomal regions included, even with --diploid-parx-genome grch38. The PARs are
diploid in a male, so every PAR gene (on this panel: CRLF2) is reported 1.0 too high. The post-step:
  1. decide  the sample is treated as male when genemetrics gets `-x male`, or no -x with --diploid-parx-genome set
             and `cnvkit.py sex --diploid-parx-genome` on the .cnr reporting Male (the guess genemetrics makes).
             Without --diploid-parx-genome no post-step is run.
  2. shift   copies of <s>.cnr and <s>.call.cns with 1.0 subtracted from log2 of every bin / segment lying wholly
             inside GRCh38 chrX PAR1 [10000, 2781479) or PAR2 [155701382, 156030895) (0-based half-open). Rows that
             overlap a PAR boundary are counted (span) and left unchanged, so a PAR gene whose segment continues into
             non-PAR chrX keeps genemetrics' value (+1 included): par_check = span.
  3. (app)   genemetrics re-run in the image on the copies, with the original arguments.
  4. compare every non-PAR row of the new genemetrics is byte-identical to the original (exit 1 otherwise); PAR rows
             after the fix equal the call.cns segment log2 (par_check 1), or are explained by a span (par_check span).

Usage
  cnvkit_parfix.py decide  --sex-arg male|female|'' --cnr-sex <sex_cnr.tsv> --parx grch38|''
      prints "<1|0>\t<reason>"
  cnvkit_parfix.py shift   <s.cnr> <s.call.cns> <outdir>
      writes <outdir>/<s.cnr>, <outdir>/<s.call.cns>; prints "par_bins=N span_bins=N par_segs=N span_segs=N"
  cnvkit_parfix.py compare <original.genemetrics.tsv> <new.genemetrics.tsv> <s.call.cns>
      prints key=value statistics, ending with par_check; exit 1 if any non-PAR row differs or the headers differ
"""
import argparse
import collections
import csv
import os
import sys

PARS = [(10000, 2781479), (155701382, 156030895)]   # GRCh38 chrX PAR1, PAR2; 0-based half-open
EMPTY_GM_HEADER = 'gene\tchromosome\tstart\tend\tlog2\n'   # CNVkit's header when genemetrics finds no gene


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


def cnr_sex(path):
    """The `sex` call from a one-sample `cnvkit.py sex` table."""
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    if len(rows) != 1:
        raise ValueError(f'{path}: expected one sample row, found {len(rows)}')
    return rows[0]['sex']


def decide(sex_arg, inferred, parx):
    """(treated_male, reason). sex_arg: the -x given to genemetrics ('' = none); inferred: cnvkit.py sex on the .cnr."""
    if sex_arg not in ('male', 'female', ''):
        raise ValueError(f'sex_arg must be male, female or empty, not {sex_arg!r}')
    if not parx:
        return False, 'no --diploid-parx-genome: post-step not applicable'
    if sex_arg == 'male':
        return True, '-x male'
    if sex_arg == 'female':
        return False, '-x female'
    if inferred == 'Male':
        return True, 'no -x; cnvkit.py sex on the .cnr: Male'
    return False, f'no -x; cnvkit.py sex on the .cnr: {inferred}'


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


def par_counts(path):
    n = collections.Counter()
    with open(path) as f:
        next(f)
        for l in f:
            c = l.split('\t'); n[par_state(c[0], int(c[1]), int(c[2]))] += 1
    return n['in'], n['span']


def compare(orig, new, call_cns):
    """genemetrics_parfix.compare for one sample; returns the statistics (raises AssertionError on a failure)."""
    rec = {}
    rec['par_segs'], rec['span_segs'] = par_counts(call_cns)
    h0, r0 = read_rows(orig); h1, r1 = read_rows(new)
    # CNVkit writes a 5-column header when genemetrics finds no gene (empty table); otherwise headers must match
    assert h0 == h1 or (not r1 and h1 == EMPTY_GM_HEADER), ('headers differ', h0, h1)
    rec['empty_after'] = int(not r1)
    isp = lambda l: (lambda c: par_state(c[1], int(c[2]), int(c[3])) != 'out')(l.split('\t'))
    n0 = [l for l in r0 if not isp(l)]; n1 = [l for l in r1 if not isp(l)]
    p0 = [l for l in r0 if isp(l)]; p1 = [l for l in r1 if isp(l)]
    rec['nonpar_identical'] = int(n0 == n1)
    assert n0 == n1, 'non-PAR rows differ'
    rec['par_rows_before'] = len(p0); rec['par_rows_after'] = len(p1)
    rec['rows_changed'] = len(collections.Counter(p0) - collections.Counter(p1)) + len(collections.Counter(p1) - collections.Counter(p0))
    rec['file_changed'] = int(r0 != r1)
    rec['par_before'] = ';'.join(sorted({f"{l.split(chr(9))[0]}:{l.split(chr(9))[4]}" for l in p0}))
    rec['par_after'] = ';'.join(sorted({f"{l.split(chr(9))[0]}:{l.split(chr(9))[4]}" for l in p1}))
    # PAR rows after the fix must equal the call.cns segment log2 (no +1)
    segs = seg_lookup(call_cns)
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('decide')
    p.add_argument('--sex-arg', required=True); p.add_argument('--cnr-sex', required=True)
    p.add_argument('--parx', required=True)
    p = sp.add_parser('shift')
    p.add_argument('cnr'); p.add_argument('call_cns'); p.add_argument('outdir')
    p = sp.add_parser('compare')
    p.add_argument('orig'); p.add_argument('new'); p.add_argument('call_cns')
    a = ap.parse_args(argv)
    try:
        if a.cmd == 'decide':
            male, reason = decide(a.sex_arg, cnr_sex(a.cnr_sex) if a.parx else '', a.parx)
            print(f'{int(male)}\t{reason}')
        elif a.cmd == 'shift':
            os.makedirs(a.outdir, exist_ok=True)
            nb, nbs = shift_par(a.cnr, os.path.join(a.outdir, os.path.basename(a.cnr)), 1, 4)
            ns, nss = shift_par(a.call_cns, os.path.join(a.outdir, os.path.basename(a.call_cns)), 1, 4)
            print(f'par_bins={nb} span_bins={nbs} par_segs={ns} span_segs={nss}')
        else:
            rec = compare(a.orig, a.new, a.call_cns)
            print(' '.join(f'{k}={v}' for k, v in rec.items()))
            if rec['par_check'] == 0:
                print('WARNING: par_check 0: a PAR row does not match its segment after the fix', file=sys.stderr)
    except (AssertionError, ValueError) as e:
        sys.exit(f'cnvkit_parfix.py {a.cmd}: FAILED: {e}')


if __name__ == '__main__':
    main()
