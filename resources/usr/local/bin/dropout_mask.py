#!/usr/bin/env python3
"""Dropout mask (DI-4139): weight 0 on isolated probe dropouts in the segmentation input. Stdlib only.

Source: the `mask` function and constants of `analysis/pon_comparison/scripts/deadprobe_mask.py` at commit ddad79b
of cnv-backbone-purple-atlas (project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JFGy8Vj4Yggp4Zj50GVVzF2Q, md5
382acdecb5be92d6074f8957362bdfd9), arm M only. The `blacklist` mode and the L5 / L3 / MB arms are removed. The
positional CLI and the output bytes are unchanged, so the command is the one validated.
Evaluation: Confluence DV 4850548748 "CNVkit evaluation for Atlas — probe dropouts and the dropout mask" (UMI-consensus
T3 and non-UMI T5); DI-4139 comment 72496 (non-UMI: rule and constants unchanged).

Rule (arm M). An autosomal bin is masked (weight set to 0 in a copy of the .cnr) when all three hold:
  1. .cnr log2 < -3;
  2. it is in a run of <= 2 consecutive bins with log2 < -3 (same chromosome, adjacent rows; longer runs are
     callable deletions and are kept);
  3. it has < 10 raw read-equivalents in the .targetcoverage.cnn (depth * length / 144; a bin absent from the
     .cnn counts as 0).
chrX and chrY are never masked. CNVkit `segment` drops weight-0 bins before CBS, the same operation as
--drop-low-coverage, so only the segmentation input changes.

Usage
  dropout_mask.py mask M <in.cnr> <raw.targetcoverage.cnn> <specimen> - <out.cnr> <qc.tsv>
      writes the masked copy of the .cnr and the QC table of masked bins (header only when none). <specimen> is
      not used by arm M (kept for CLI compatibility); the blacklist argument must be '-'.
  dropout_mask.py qc-header <qc.tsv>
      writes the header-only QC table (the app's output when the mask is switched off).
"""
import sys

DEEP, RUN_MAX, RAW_READS, READ_LEN = -3.0, 2, 10, 144
QC_HEADER = 'arm\tchromosome\tstart\tend\tgene\tlog2\traw_reads\treason\n'


def rows_of(path):
    lines = open(path).read().splitlines()
    return lines[0], [l.split('\t') for l in lines[1:]]


def mask(arm, cnr, cnn, spec, bl, out, qc):
    """Arm M of deadprobe_mask.mask (ddad79b); returns the number of masked bins."""
    if arm != 'M':
        raise ValueError(f'only arm M is supported, not {arm!r}')
    if bl != '-':
        raise ValueError(f"arm M takes no blacklist: expected '-', got {bl!r}")
    hdr, rows = rows_of(cnr); h = hdr.split('\t'); li, wi = h.index('log2'), h.index('weight')
    _, raw = rows_of(cnn)
    rawd = {(r[0], r[1], r[2]): float(r[4]) * (int(r[2]) - int(r[1])) / READ_LEN for r in raw}
    l2 = [float(r[li]) for r in rows]
    auto = [r[0] not in ('chrX', 'chrY') for r in rows]
    hit, why = set(), {}
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
    with open(out, 'w') as f:
        f.write(hdr + '\n')
        for i, r in enumerate(rows):
            if i in hit:
                r = list(r); r[wi] = '0'
            f.write('\t'.join(r) + '\n')
    with open(qc, 'w') as f:
        f.write(QC_HEADER)
        for i in sorted(hit):
            r = rows[i]
            f.write(f'{arm}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{l2[i]:g}\t{rawd.get(tuple(r[:3]), 0.0):.1f}\t{why[i]}\n')
    return len(hit)


def qc_header(qc):
    with open(qc, 'w') as f:
        f.write(QC_HEADER)


def main(argv):
    if len(argv) == 9 and argv[1] == 'mask':
        try:
            n = mask(*argv[2:9])
        except ValueError as e:
            sys.exit(f'dropout_mask.py: {e}')
        print(f'[dropout] masked probes: {n}')
    elif len(argv) == 3 and argv[1] == 'qc-header':
        qc_header(argv[2])
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main(sys.argv)
