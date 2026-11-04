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
    `dropout_probes`;
  - the `cnvkit.py sex` tables are validated (port addition): one sample row with every required column, each
    once, and a value for each; a sex call CNVkit writes (Male, Female, or Unknown when the sample has no usable
    chrX bins); finite numeric log-ratios, or CNVkit's literals for no value: `nan` (no usable bins on that
    chromosome; a `nan` raw Y log-ratio is no chrY, as sex_check.py) and `NA` with an Unknown call (read as `nan`,
    as before); and a newline at the end of the last line, as CNVkit writes it (an unterminated line is truncated).
    The row must name the expected file, <sample>.targetcoverage.cnn / <sample>.cnr, optionally <run>__-prefixed;
    the whole sample ID is matched. Anything else is malformed: exit 2.
Validation: the chrY rule and the -5 cut-off, Confluence DV 4854120471 test 4 (non-UMI) and FINAL_RESULTS.md.

Usage
  cnvkit_sex_qc.py resolve --sample-id ID --bam-name NAME --override [male|female|M|F|'']
      prints "<sex_arg>\\t<applied_source>\\t<lims_sex>\\t<declared_sex>": sex_arg is male / female / '' (no -x);
      exit 2 on an override value that is not accepted (D9)
  cnvkit_sex_qc.py qc --sample-id ID --raw-sex RAW.tsv --cnr-sex CNR.tsv --lims-sex M|F|U --override-sex M|F|U
                      --applied-source SRC --par-fix-applied 0|1 --par-span-segments N --dropout-probes N
                      --out PREFIX [--run RUN]
      writes PREFIX.sex_qc.tsv and PREFIX.sex_check_mqc.tsv; exit 2 on a malformed or mismatched sex table
"""
import argparse
import csv
import io
import math
import os
import re
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
    """A validated log-ratio (see check_sex_row) as a float; CNVkit's literals `nan` and `NA` are NaN."""
    return float('nan') if x == 'NA' else float(x)


SEX_COLUMNS = ('sample', 'sex', 'X_logratio', 'Y_logratio')
SEX_VALUES = ('Male', 'Female', 'Unknown')     # CNVkit fc65941 commands.do_sex; Unknown: no usable chrX bins
# a `cnvkit.py sex` log-ratio: a finite '%.3g' number with a sign ('+0.116', '-1.13', '-1.2e+03'), the literal
# 'nan' CNVkit writes when a chromosome has no usable bins (chrY in a .cnr; the policy: no chrY), or the literal 'NA'
# it writes with an Unknown call (no statistics). Anything else, inf and overflow such as -1e999 included, is
# malformed.
NUMBER = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?')


def is_logratio(x, sex):
    return (x == 'nan' or (x == 'NA' and sex == 'Unknown')
            or (NUMBER.fullmatch(x) is not None and math.isfinite(float(x))))


def check_sex_row(path, header, row):
    """Raise ValueError unless a `cnvkit.py sex` row has every required column, each once, and a value for each, a
    sex call CNVkit writes and finite numeric (or CNVkit's `nan` / `NA`) log-ratios."""
    dup = sorted({c for c in header if header.count(c) > 1})
    if dup:     # a later duplicate must not hide an invalid required value
        raise ValueError(f'{path}: duplicate column(s) {", ".join(dup)}')
    missing = [c for c in SEX_COLUMNS if c not in header]
    if missing:
        raise ValueError(f'{path}: missing column(s) {", ".join(missing)} (columns: {", ".join(header)})')
    if len(row) != len(header):
        raise ValueError(f'{path}: the sample row has {len(row)} fields, the header {len(header)}')
    r = dict(zip(header, row, strict=True))
    if r['sex'] not in SEX_VALUES:
        raise ValueError(f"{path}: sex {r['sex']!r} is not one of {', '.join(SEX_VALUES)}")
    for c in ('X_logratio', 'Y_logratio'):
        if not is_logratio(r[c], r['sex']):
            raise ValueError(f"{path}: {c} {r[c]!r} is not a finite number (or CNVkit's literal nan, or NA with an "
                             "Unknown call)")
    return r


def read_sex_table(path):
    """The one sample row of a `cnvkit.py sex` table, validated (check_sex_row). CNVkit (pandas to_csv) ends every
    record with a newline, so a table whose last line is unterminated is truncated and malformed (a raw Y log-ratio
    -5.01e-05 cut to -5.01 is still a number, and would turn the chrY call to Female)."""
    with open(path, newline='') as fh:
        text = fh.read()
    rows = [r for r in csv.reader(io.StringIO(text), delimiter='\t') if r]
    if not rows:
        raise ValueError(f'{path}: empty table')
    if not text.endswith('\n'):
        raise ValueError(f'{path}: the last line is not newline-terminated (a truncated table)')
    if len(rows) != 2:
        raise ValueError(f'{path}: expected one sample row, found {len(rows) - 1}')
    return check_sex_row(path, rows[0], rows[1])


def table_sample(path, name, sample, suffix):
    """Check that a sex-table `sample` value is the file <sample><suffix>, optionally run-prefixed as
    <run>__<sample><suffix> (the validation tables). The whole sample ID is matched, so IDs containing `__` or
    ending in `.cnr` / `.cnn` are not split or stripped."""
    b, want = os.path.basename(name), sample + suffix
    if b != want and not (b.endswith('__' + want) and len(b) > len(want) + 2):
        raise ValueError(f'{path}: sample {b!r} is not {want!r} (or <run>__{want})')


def one_row(path, sample, suffix):
    r = read_sex_table(path)
    table_sample(path, r['sample'], sample, suffix)
    return r


def qc_row(a):
    raw, c = one_row(a.raw_sex, a.sample_id, '.targetcoverage.cnn'), one_row(a.cnr_sex, a.sample_id, '.cnr')
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
