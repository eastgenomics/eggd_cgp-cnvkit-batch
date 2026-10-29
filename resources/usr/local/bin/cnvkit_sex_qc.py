#!/usr/bin/env python3
"""Sex resolution and sex QC for one sample (DI-4098, DI-4099, DI-4100, DI-4103; decision D8). Stdlib only.

Source: `analysis/pon_comparison/scripts/sex_check.py` at commit 31a62fc of cnv-backbone-purple-atlas
(project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6X7g04YggpgK1QB942gY64, md5 a6106d50c43cd96ebde9ae3f24360e2b), ported from a
batch over `cnvkit.py sex` tables to one sample. The LIMS-letter parsing, the chrY-presence rule (raw Y log-ratio > -5)
and the sex_mismatch rule are unchanged, and so are the column names, the MultiQC header and the values, except
where D8 extends them:
  - the sex can come from an override input, and the LIMS letter from the BAM name when the sample ID has no 5th
    field; `applied_source` names the source (override / sample_name / bam_name / cnvkit_inference) and
    `applied_sex` gives the inferred sex for samples without -x, e.g. `cnvkit_inferred(Male)`;
  - `sex_mismatch` compares the declared sex (override, else LIMS) with chrY; with no override it is sex_check.py's;
  - appended columns `par_fix_applied`, `par_span_segments`, `dropout_probes`; the MultiQC table adds
    `dropout_probes`.
Validation: the chrY rule and the -5 cut-off, Confluence DV 4854120471 test 4 (non-UMI) and FINAL_RESULTS.md.

Usage
  cnvkit_sex_qc.py resolve --sample-id ID --bam-name NAME --override [male|female|M|F|'']
      prints "<sex_arg>\\t<applied_source>\\t<lims_sex>\\t<declared_sex>": sex_arg is male / female / '' (no -x);
      exit 2 on an override value that is not accepted (D9)
  cnvkit_sex_qc.py qc --sample-id ID --raw-sex RAW.tsv --cnr-sex CNR.tsv --lims-sex M|F|U --override-sex M|F|U
                      --applied-source SRC --par-fix-applied 0|1 --par-span-segments N --dropout-probes N
                      --out PREFIX [--run RUN]
      writes PREFIX.sex_qc.tsv and PREFIX.sex_check_mqc.tsv
"""
import argparse
import csv
import os
import sys

Y_PRESENT = -5.0
SOURCES = ('override', 'sample_name', 'bam_name', 'cnvkit_inference')
OVERRIDE = {'male': 'M', 'm': 'M', 'female': 'F', 'f': 'F', '': 'U'}
TO_ARG = {'M': 'male', 'F': 'female', 'U': ''}
MQC_HEADER = ("# id: 'cnvkit_sex_check'\n# section_name: 'CNVkit sex check'\n"
              "# description: 'LIMS sex vs chrY presence in raw CNVkit coverage (sex_mismatch = 1 needs review); "
              "X log-ratios from raw coverage and from the .cnr for context.'\n# plot_type: 'table'\n")
MQC_COLS = ['lims_sex', 'chrY_sex', 'sex_mismatch', 'raw_Y_logratio', 'raw_X_logratio', 'cnr_X_logratio',
            'applied_sex', 'dropout_probes']


def lims_letter(name):
    """Field 5 of the '-'-split name: M / F, anything else U (sex_check.py); None when the name has < 5 fields."""
    f = name.split('-')
    if len(f) <= 4:
        return None
    return f[4] if f[4] in ('M', 'F') else 'U'


def parse_override(value):
    key = (value or '').strip().lower()
    if key not in OVERRIDE:
        raise ValueError(f"sex override {value!r} not accepted: use male, female, M or F (or leave it empty)")
    return OVERRIDE[key]


def resolve(sample_id, bam_name, override):
    """(sex_arg, applied_source, lims_sex, declared_sex)."""
    ov = parse_override(override)
    lims, src = lims_letter(sample_id), 'sample_name'
    if lims is None:
        stem = os.path.basename(bam_name or '')
        stem = stem[:-4] if stem.endswith('.bam') else stem
        lims, src = lims_letter(stem), 'bam_name'
        if lims is None:
            lims = 'U'
    if ov != 'U':
        return TO_ARG[ov], 'override', lims, ov
    if lims != 'U':
        return TO_ARG[lims], src, lims, lims
    return '', 'cnvkit_inference', lims, 'U'


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def one_row(path, sample):
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    if len(rows) != 1:
        raise ValueError(f'{path}: expected one sample row, found {len(rows)}')
    b = os.path.basename(rows[0]['sample']).rpartition('__')[2]
    for ext in ('.targetcoverage.cnn', '.cnr', '.cnn'):
        if b.endswith(ext):
            b = b[:-len(ext)]
    if b != sample:
        raise ValueError(f'{path}: sample {b!r} is not {sample!r}')
    return rows[0]


def qc_row(a):
    raw, c = one_row(a.raw_sex, a.sample_id), one_row(a.cnr_sex, a.sample_id)
    if a.applied_source not in SOURCES:
        raise ValueError(f'applied source {a.applied_source!r} not one of {SOURCES}')
    ry = fnum(raw['Y_logratio'])
    ypres = ry == ry and ry > Y_PRESENT
    ysex = 'Male' if ypres else 'Female'
    declared = a.override_sex if a.override_sex != 'U' else a.lims_sex
    mm = 'NA' if declared == 'U' else str(int({'M': 'Male', 'F': 'Female'}[declared] != ysex))
    applied = TO_ARG[declared] or f"cnvkit_inferred({c.get('sex', 'NA')})"
    return dict(run=a.run, sample=a.sample_id, lims_sex=a.lims_sex, applied_sex=applied, applied_source=a.applied_source,
                raw_Y_logratio=raw['Y_logratio'], chrY_present=int(ypres), chrY_sex=ysex,
                raw_X_logratio=raw['X_logratio'], raw_cnvkit_sex=raw['sex'],
                cnr_X_logratio=c.get('X_logratio', 'NA'), cnr_cnvkit_sex=c.get('sex', 'NA'),
                sex_mismatch=mm, par_fix_applied=a.par_fix_applied, par_span_segments=a.par_span_segments,
                dropout_probes=a.dropout_probes)


def write_qc(r, out):
    with open(out + '.sex_qc.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(r), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerow(r)
    with open(out + '.sex_check_mqc.tsv', 'w') as fh:
        fh.write(MQC_HEADER)
        fh.write('Sample\t' + '\t'.join(MQC_COLS) + '\n')
        fh.write((r['run'] + '__' if r['run'] else '') + r['sample'] + '\t' + '\t'.join(str(r[c]) for c in MQC_COLS) + '\n')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('resolve')
    p.add_argument('--sample-id', required=True); p.add_argument('--bam-name', default='')
    p.add_argument('--override', default='')
    p = sp.add_parser('qc')
    for opt in ('--sample-id', '--raw-sex', '--cnr-sex', '--applied-source', '--out'):
        p.add_argument(opt, required=True)
    p.add_argument('--lims-sex', required=True, choices=('M', 'F', 'U'))
    p.add_argument('--override-sex', required=True, choices=('M', 'F', 'U'))
    p.add_argument('--par-fix-applied', required=True, choices=('0', '1'))
    p.add_argument('--par-span-segments', required=True, type=int)
    p.add_argument('--dropout-probes', required=True, type=int)
    p.add_argument('--run', default='')
    a = ap.parse_args(argv)
    try:
        if a.cmd == 'resolve':
            print('\t'.join(resolve(a.sample_id, a.bam_name, a.override)))
        else:
            r = qc_row(a)
            write_qc(r, a.out)
            print(f"[sex] chrY {r['chrY_sex']} (raw Y log-ratio {r['raw_Y_logratio']}); "
                  f"lims {r['lims_sex']}; sex_mismatch {r['sex_mismatch']}")
    except ValueError as e:
        print(f'cnvkit_sex_qc.py {a.cmd}: ERROR: {e}', file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
