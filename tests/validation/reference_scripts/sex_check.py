#!/usr/bin/env python3
"""CNVkit sex check for the Atlas CGP read-out (reference implementation of the batch-app step).

Per sample:
  LIMS sex       field 5 of the '-'-split sample name: M / F / U (anything else = U)
  raw sex call   `cnvkit.py sex <S>.targetcoverage.cnn --diploid-parx-genome grch38`
                 Raw coverage keeps the non-PAR chrY bins (UTY, RPS4Y2). In the .cnr these bins are always masked,
                 because the PoN depth is the robust centre across mostly female members, i.e. 0.
  cnr sex call   `cnvkit.py sex <S>.cnr --diploid-parx-genome grch38`
                 chrX only (the .cnr has no chrY), so it reflects tumour chrX copy number. Reported, not used for the flag.
  chrY present   raw Y_logratio > -5 (absent chrY sits at about -19 to -29; males at about -0.9 to -2.5 on this panel)
  chrY sex       Male if chrY present, else Female
  applied sex    LIMS M/F (passed as -x to call and genemetrics); LIMS U: CNVkit's own inference (no -x)
  sex_mismatch   1 if LIMS is M/F and chrY sex differs; 0 if they agree; NA if LIMS is U

Writes <out>.sex_qc.tsv (all columns) and <out>.sex_check_mqc.tsv (MultiQC custom content: table plus General
Statistics columns).

Usage: sex_check.py --raw-sex RAW_SEX.tsv --cnr-sex CNR_SEX.tsv --out PREFIX
       (`cnvkit.py sex` tables; the sample is matched on the file basename without .targetcoverage.cnn / .cnr,
        with an optional '<run>__' prefix; samples are keyed on (run, sample))
"""
import argparse, csv, os

Y_PRESENT = -5.0

ap = argparse.ArgumentParser()
ap.add_argument('--raw-sex', required=True)
ap.add_argument('--cnr-sex', required=True)
ap.add_argument('--out', required=True)
a = ap.parse_args()


def load(p):
    d = {}
    for r in csv.DictReader(open(p), delimiter='\t'):
        run, _, b = os.path.basename(r['sample']).rpartition('__')
        for ext in ('.targetcoverage.cnn', '.cnr', '.cnn'):
            if b.endswith(ext):
                b = b[:-len(ext)]
        d[(run, b)] = r          # replicates share a sample name across runs
    return d


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


raw, cnr = load(a.raw_sex), load(a.cnr_sex)
rows = []
for key in sorted(raw):
    run, s = key
    f = s.split('-')
    lims = f[4] if len(f) > 4 and f[4] in ('M', 'F') else 'U'
    ry, rx = fnum(raw[key]['Y_logratio']), fnum(raw[key]['X_logratio'])
    ypres = ry == ry and ry > Y_PRESENT
    ysex = 'Male' if ypres else 'Female'
    c = cnr.get(key, {})
    mm = 'NA' if lims == 'U' else str(int({'M': 'Male', 'F': 'Female'}[lims] != ysex))
    rows.append(dict(run=run, sample=s, lims_sex=lims,
                     applied_sex={'M': 'male', 'F': 'female'}.get(lims, 'cnvkit_inferred'),
                     applied_source='LIMS' if lims != 'U' else 'CNVkit inference (.cnr)',
                     raw_Y_logratio=raw[key]['Y_logratio'], chrY_present=int(ypres), chrY_sex=ysex,
                     raw_X_logratio=raw[key]['X_logratio'], raw_cnvkit_sex=raw[key]['sex'],
                     cnr_X_logratio=c.get('X_logratio', 'NA'), cnr_cnvkit_sex=c.get('sex', 'NA'),
                     sex_mismatch=mm))
with open(a.out + '.sex_qc.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
with open(a.out + '.sex_check_mqc.tsv', 'w') as fh:
    fh.write("# id: 'cnvkit_sex_check'\n# section_name: 'CNVkit sex check'\n"
             "# description: 'LIMS sex vs chrY presence in raw CNVkit coverage (sex_mismatch = 1 needs review); "
             "X log-ratios from raw coverage and from the .cnr for context.'\n# plot_type: 'table'\n")
    cols = ['lims_sex', 'chrY_sex', 'sex_mismatch', 'raw_Y_logratio', 'raw_X_logratio', 'cnr_X_logratio', 'applied_sex']
    fh.write('Sample\t' + '\t'.join(cols) + '\n')
    for r in rows:
        fh.write((r['run'] + '__' if r['run'] else '') + r['sample'] + '\t' + '\t'.join(str(r[c]) for c in cols) + '\n')
n = len(rows); known = [r for r in rows if r['sex_mismatch'] != 'NA']
print(f'{n} samples; LIMS M/F {len(known)}; sex_mismatch=1: {sum(r["sex_mismatch"] == "1" for r in known)}')
