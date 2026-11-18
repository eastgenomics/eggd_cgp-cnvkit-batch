#!/usr/bin/env python3
"""Dead-probe handling before segmentation (../sim2_plan.md §11). Stdlib only; runs in cgp-cnvkit:1.0.0.

  deadprobe_mask.py blacklist <out.tsv> <specimen_map.tsv> <cnr> [<cnr> ...]
      per autosomal bin: specimens (one sequencing each) where the bin is an isolated dead probe (log2 < -3 in a run of
      <= 2 consecutive deep bins; long runs are deletions, e.g. CDKN2A/MTAP, and are not counted) -> site table.
      specimen_map columns: sample_id, specimen, pon_exclusion (sequencings excluded for QC, qc_*, are skipped).
  deadprobe_mask.py mask <arm> <in.cnr> <raw.targetcoverage.cnn> <specimen> <blacklist.tsv|-> <out.cnr> <qc.tsv>
      writes a copy of the .cnr with weight = 0 on the masked bins (CNVkit `segment` drops weight-0 bins before CBS,
      the same operation as --drop-low-coverage) and a QC table of the masked bins.

Arms
  L5 / L3   every bin with log2 < -5 / < -3 (any run length): emulates an adjustable --drop-low-coverage cut-off
            (fc65941 hard-codes log2 < -15 or depth 0)
  M         isolated dead probes: log2 < -3, in a run of <= 2 consecutive deep bins (same chromosome, adjacent
            rows), raw depth < 10 read-equivalents (depth * length / 144) in the .targetcoverage.cnn
  MB        M + recurrent sites: bins deep (log2 < -3) in >= 3 specimens other than this one, masked when log2 < -1
chrX/chrY bins are never masked (sex-dependent depth).
"""
import csv, sys

DEEP, RUN_MAX, RAW_READS, READ_LEN, REC_MIN, REC_LOG2 = -3.0, 2, 10, 144, 3, -1.0


def rows_of(path):
    lines = open(path).read().splitlines()
    return lines[0], [l.split('\t') for l in lines[1:]]


def blacklist(out, spmap, cnrs):
    m = {r['sample_id']: r for r in csv.DictReader(open(spmap), delimiter='\t')}
    seen, sites = set(), {}
    for p in sorted(cnrs):
        s = p.rsplit('/', 1)[-1][:-4]; spec = m[s]['specimen']
        if spec in seen or m[s].get('pon_exclusion', '').startswith('qc'):
            continue
        seen.add(spec)
        hdr, rows = rows_of(p); h = hdr.split('\t'); li = h.index('log2')
        deep = [r[0] not in ('chrX', 'chrY') and float(r[li]) < DEEP for r in rows]
        i = 0
        while i < len(rows):
            if deep[i]:
                j = i
                while j + 1 < len(rows) and rows[j + 1][0] == rows[i][0] and deep[j + 1]:
                    j += 1
                if j - i + 1 <= RUN_MAX:
                    for k in range(i, j + 1):
                        sites.setdefault(f'{rows[k][0]}:{rows[k][1]}-{rows[k][2]}', [rows[k][3], []])[1].append(spec)
                i = j + 1
            else:
                i += 1
    with open(out, 'w') as f:
        f.write('bin\tgene\tn_specimens\tspecimens\n')
        for b, (g, ss) in sorted(sites.items(), key=lambda x: -len(x[1][1])):
            f.write(f'{b}\t{g}\t{len(ss)}\t{",".join(sorted(ss))}\n')
    print(f'{len(seen)} specimens; {len(sites)} deep sites; recurrent (>= {REC_MIN}): '
          f'{sum(1 for v in sites.values() if len(v[1]) >= REC_MIN)}')


def mask(arm, cnr, cnn, spec, bl, out, qc):
    hdr, rows = rows_of(cnr); h = hdr.split('\t'); li, wi = h.index('log2'), h.index('weight')
    _, raw = rows_of(cnn)
    rawd = {(r[0], r[1], r[2]): float(r[4]) * (int(r[2]) - int(r[1])) / READ_LEN for r in raw}
    l2 = [float(r[li]) for r in rows]
    auto = [r[0] not in ('chrX', 'chrY') for r in rows]
    hit, why = set(), {}
    if arm in ('L5', 'L3'):
        t = -5.0 if arm == 'L5' else -3.0
        for i, v in enumerate(l2):
            if auto[i] and v < t:
                hit.add(i); why[i] = f'log2<{t:g}'
    if arm in ('M', 'MB'):
        i = 0
        while i < len(rows):
            if auto[i] and l2[i] < DEEP:
                j = i
                while j + 1 < len(rows) and rows[j + 1][0] == rows[i][0] and l2[j + 1] < DEEP:
                    j += 1
                if j - i + 1 <= RUN_MAX:
                    for k in range(i, j + 1):
                        if rawd.get(tuple(rows[k][:3]), 0.0) < RAW_READS:
                            hit.add(k); why[k] = f'isolated_run{j - i + 1}'
                i = j + 1
            else:
                i += 1
    if arm == 'MB' and bl != '-':
        rec = {}
        for r in csv.DictReader(open(bl), delimiter='\t'):
            others = [x for x in r['specimens'].split(',') if x and x != spec]
            if len(others) >= REC_MIN:
                rec[r['bin']] = len(others)
        for i, r in enumerate(rows):
            k = f'{r[0]}:{r[1]}-{r[2]}'
            if auto[i] and k in rec and l2[i] < REC_LOG2 and i not in hit:
                hit.add(i); why[i] = f'recurrent_site_{rec[k]}'
    with open(out, 'w') as f:
        f.write(hdr + '\n')
        for i, r in enumerate(rows):
            if i in hit:
                r = list(r); r[wi] = '0'
            f.write('\t'.join(r) + '\n')
    with open(qc, 'w') as f:
        f.write('arm\tchromosome\tstart\tend\tgene\tlog2\traw_reads\treason\n')
        for i in sorted(hit):
            r = rows[i]
            f.write(f'{arm}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{l2[i]:g}\t{rawd.get(tuple(r[:3]), 0.0):.1f}\t{why[i]}\n')


if __name__ == '__main__':
    if sys.argv[1] == 'blacklist':
        blacklist(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif sys.argv[1] == 'mask':
        mask(*sys.argv[2:9])
    else:
        sys.exit(__doc__)
