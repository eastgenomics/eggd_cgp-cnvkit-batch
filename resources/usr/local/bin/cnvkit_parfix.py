#!/usr/bin/env python3
"""PAR post-step for CNVkit genemetrics in samples treated as male (DI-4101). Stdlib only.

Source: `analysis/pon_comparison/scripts/genemetrics_parfix.py` at commit f933f7a of cnv-backbone-purple-atlas
(project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6PpjQ4Yggx4JxXFkQqZxGb, md5 51ada0e597420655c89c86e1a2e13ea6). The batch
driver (plans, arm directories, the in-image worker and the hard-coded work directory) is replaced by three per-sample
subcommands; `par_state`, `shift_par`, the male decision and the `compare` checks are copied unchanged, except that the
app fails where the validation only logged: a malformed .cnr sex table is an error, and so is par_check 0, which the
port decides per PAR row against that row's own call.cns segment (step 4; on the 60 validated males per set every
par_check is unchanged). The app runs genemetrics in the image between `shift` and `compare` with the original
arguments, as the validation did.
Validation: FINAL_RESULTS.md (project-JF3p1f04Yggp0vf3PYjy7Q0F:file-JF6YBbj4YggvFK194V51V7pX) and DI-4118 (non-UMI).

CNVkit limitation (fc65941, `CopyNumArray.shift_xx`): for a sample that genemetrics treats as male, every chrX bin and
segment is shifted +1.0, the pseudo-autosomal regions included, even with --diploid-parx-genome grch38. The PARs are
diploid in a male, so every PAR gene (on this panel: CRLF2) is reported 1.0 too high. The post-step:
  1. decide  the sample is treated as male when genemetrics gets `-x male`, or no -x with --diploid-parx-genome set
             and `cnvkit.py sex --diploid-parx-genome` on the .cnr reporting Male (the guess genemetrics makes).
             Without --diploid-parx-genome no post-step is run. Port addition: the sex table is validated
             (cnr_sex), so a malformed one cannot silently skip the post-step.
  2. shift   copies of <s>.cnr and <s>.call.cns with 1.0 subtracted from log2 of every bin / segment lying wholly
             inside GRCh38 chrX PAR1 [10000, 2781479) or PAR2 [155701382, 156030895) (0-based half-open). Rows that
             overlap a PAR boundary are counted (span) and left unchanged, so a PAR gene whose segment continues into
             non-PAR chrX keeps genemetrics' value (+1 included): par_check = span.
  3. (app)   genemetrics re-run in the image on the copies, with the original arguments.
  4. compare every non-PAR row of the new genemetrics is byte-identical to the original (exit 1 otherwise). Port
             additions: both genemetrics tables are read with the checks the comparison needs (read_genemetrics),
             and each PAR row is judged against its own call.cns segment (par_row_verdicts): a row of a segment
             the shift changed must equal the segment log2 L after the fix (par_check 1) and L + 1.0 before it,
             each at a level genemetrics reports (|level| >= 0.3), matched to its before-row by gene, coordinates
             and every field but log2; a row of a segment the shift left unchanged (one spanning a PAR boundary)
             must be unchanged and equal L + 1.0 (par_check span; the gene@segment is printed as span_genes). With
             overlapping segments (overlapping targets) a row is explained when any segment it overlaps or lies
             in explains it (row_segments), as genemetrics reports a gene once per segment it overlaps. Anything else is an unexplained PAR
             discrepancy (par_check 0): exit 1, and the app publishes no outputs. 'Equal' means CNVkit's %.6g
             print of that value, within half a unit in the 6th significant digit at its magnitude (prints_as).

Usage
  cnvkit_parfix.py decide  --sex-arg male|female|'' --cnr-sex <sex_cnr.tsv> --parx grch38|''
      prints "<1|0>\t<reason>"; exit 1 if the sex table is malformed (cnr_sex)
  cnvkit_parfix.py shift   <s.cnr> <s.call.cns> <outdir>
      writes <outdir>/<s.cnr>, <outdir>/<s.call.cns>; prints "par_bins=N span_bins=N par_segs=N span_segs=N"
  cnvkit_parfix.py compare <original.genemetrics.tsv> <new.genemetrics.tsv> <s.call.cns>
      prints key=value statistics, ending with par_check; exit 1 if a genemetrics table cannot be read
      (read_genemetrics), the headers differ, any non-PAR row differs or par_check is 0
"""
import argparse
import collections
import csv
import io
import math
import os
import re
import sys

PARS = [(10000, 2781479), (155701382, 156030895)]   # GRCh38 chrX PAR1, PAR2; 0-based half-open
EMPTY_GM_HEADER = 'gene\tchromosome\tstart\tend\tlog2\n'   # CNVkit's header when genemetrics finds no gene


def par_state(chrom, start, end):
    """'in' (wholly inside PAR boundaries), 'span' (overlaps a PAR boundary) or 'out' (wholly outside PAR boundaries)."""
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


SEX_COLUMNS = ('sample', 'sex', 'X_logratio', 'Y_logratio')
SEX_VALUES = ('Male', 'Female', 'Unknown')     # CNVkit fc65941 commands.do_sex; Unknown: no usable chrX bins
NUMBER = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?')   # finite; CNVkit's literals nan and NA aside


def is_logratio(x, sex):
    """A `cnvkit.py sex` log-ratio: a finite number, CNVkit's literal nan (no usable bins on that chromosome) or NA
    (no statistics: written with an Unknown call). inf and overflow such as -1e999 are not."""
    return (x == 'nan' or (x == 'NA' and sex == 'Unknown')
            or (NUMBER.fullmatch(x) is not None and math.isfinite(float(x))))


def cnr_sex(path):
    """The `sex` call from a one-sample `cnvkit.py sex` table, validated as in cnvkit_sex_qc.py (port addition): every
    required column, each once, and a value for each; a sex call CNVkit writes (Male, Female or Unknown); finite
    numeric log-ratios or CNVkit's nan / NA; a newline at the end of the last line, as CNVkit writes it. Anything else
    is malformed and raises ValueError, so a bad or truncated table cannot silently skip the post-step."""
    with open(path, newline='') as fh:
        text = fh.read()
    rows = [r for r in csv.reader(io.StringIO(text), delimiter='\t') if r]
    if len(rows) != 2:
        raise ValueError(f'{path}: expected one sample row, found {max(len(rows) - 1, 0)}')
    if not text.endswith('\n'):    # a truncated last value can still be a number
        raise ValueError(f'{path}: the last line is not newline-terminated (a truncated table)')
    header, row = rows
    dup = sorted({c for c in header if header.count(c) > 1})
    if dup:
        raise ValueError(f'{path}: duplicate column(s) {", ".join(dup)}')
    missing = [c for c in SEX_COLUMNS if c not in header]
    if missing:
        raise ValueError(f'{path}: not a cnvkit.py sex table: missing column(s) {", ".join(missing)}')
    if len(row) != len(header):
        raise ValueError(f'{path}: the sample row has {len(row)} fields, the header {len(header)}')
    r = dict(zip(header, row, strict=True))
    if r['sex'] not in SEX_VALUES:
        raise ValueError(f"{path}: sex {r['sex']!r} is not one of {', '.join(SEX_VALUES)}")
    for c in ('X_logratio', 'Y_logratio'):
        if not is_logratio(r[c], r['sex']):
            raise ValueError(f"{path}: {c} {r[c]!r} is not a finite number (or CNVkit's literal nan, or NA with an "
                             "Unknown call)")
    return r['sex']


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
    if inferred not in SEX_VALUES:
        raise ValueError(f'cnvkit.py sex on the .cnr: {inferred!r} is not one of {", ".join(SEX_VALUES)}')
    if inferred == 'Male':
        return True, 'no -x; cnvkit.py sex on the .cnr: Male'
    return False, f'no -x; cnvkit.py sex on the .cnr: {inferred}'


class CheckFailed(ValueError):
    """A post-step check failed (the reference used assert, which python -O removes)."""


GM_READ = ('gene', 'chromosome', 'start', 'end', 'log2')   # the genemetrics columns compare reads, by position


def _parses(convert, v):
    try:
        convert(v)
    except ValueError:
        return False
    return True


def read_genemetrics(path):
    """(header line, row lines) of a genemetrics table, with only the checks compare needs to read it (port
    addition): a header whose first columns are gene, chromosome, start, end and log2, and rows with one field per
    header column, integer start and end and a numeric log2. An empty file (no header) raises ValueError, and so
    does any other header or row, so a truncated table cannot pass the comparison."""
    with open(path) as f:
        hdr = f.readline()
        rows = list(f)
    if not hdr:
        raise ValueError(f'{path}: empty file, no header')
    names = hdr.rstrip('\n').split('\t')
    if tuple(names[:len(GM_READ)]) != GM_READ:
        raise ValueError(f'{path}: not a genemetrics header (it must start {", ".join(GM_READ)}): {hdr!r}')
    for i, l in enumerate(rows, 2):
        c = l.rstrip('\n').split('\t')
        if not (len(c) == len(names) and _parses(int, c[2]) and _parses(int, c[3]) and _parses(float, c[4])):
            raise ValueError(f'{path}: line {i} is not a {len(names)}-field genemetrics row with integer start and '
                             f'end and a numeric log2: {l!r}')
    return hdr, rows


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


MALE_SHIFT = 1.0      # CNVkit's chrX shift for a sample treated as male (CopyNumArray.shift_xx: log2 += 1.0)
GM_THRESHOLD = 0.3    # genemetrics -t: a gene is reported when |segment log2| >= 0.3
SIG_DIGITS = 6        # CNVkit writes every float with %.6g (cmdutil.write_dataframe, tabio.write)
FLOAT_SLACK = 1e-12   # float noise: (L - 1.0) + 1.0, and pandas' text parser against Python's on the same digits


def printed_tol(*xs):
    """The largest difference between a value and its %.6g print: half a unit in the 6th significant digit at the
    magnitude of the larger of xs (so a print that rounds up into the next decade, 10.99996 -> 11, is covered), plus
    float noise. A fixed tolerance does not scale: 1e-5 is below half a unit from 10 up. A non-finite value has no
    tolerance (it equals nothing)."""
    m = max(abs(x) for x in xs)
    if not math.isfinite(m):
        return 0.0
    return (0.5 * 10.0 ** (math.floor(math.log10(m)) - (SIG_DIGITS - 1)) if m else 0.0) + FLOAT_SLACK


def prints_as(v, x):
    """v, a log2 read from a genemetrics table, is what CNVkit's %.6g makes of the level x."""
    return abs(v - x) <= printed_tol(v, x)


def before_level(L):
    """The segment level genemetrics reports without the fix: the male shift on the call.cns log2 L."""
    return L + MALE_SHIFT


def after_level(L):
    """The level after the fix, computed as the post-step and CNVkit do: shift writes repr(L - 1.0) (exact), and
    genemetrics adds the male shift back."""
    return (L - MALE_SHIFT) + MALE_SHIFT


def may_be_reported(x):
    """genemetrics can report a segment at level x: |x| >= 0.3, up to float noise at the threshold."""
    return abs(x) >= GM_THRESHOLD - FLOAT_SLACK


def may_be_unreported(x):
    return abs(x) < GM_THRESHOLD + FLOAT_SLACK


def row_segments(segs, s, e):
    """The chrX call.cns segments (start, end, log2) a genemetrics row [s, e) can come from: those it overlaps, and
    those that contain it (the base's lookup; this also covers a zero-length row at a segment's start). genemetrics
    reports a gene once per segment whose range overlaps the gene's bins (by_ranges, mode 'outer'), so with
    non-overlapping segments this is the one segment that contains the row; with overlapping segments (overlapping
    targets) it can be several, and a row may extend past the segment it comes from."""
    return [x for x in segs['chrX'] if (x[0] < e and s < x[1]) or (x[0] <= s and e <= x[1])]


def shifted(seg):
    """shift_par changed the segment: it lies wholly in a PAR."""
    return par_state('chrX', seg[0], seg[1]) == 'in'


def without_log2(c):
    """A genemetrics row's fields other than log2 (index 4): all that the post-step leaves unchanged."""
    return tuple(c[:4] + c[5:])


def par_row_verdicts(p0, p1, segs):
    """Port addition (A1): every PAR row judged against its own call.cns segment, the segment given unshifted. With L
    the segment log2, genemetrics reports a gene at L + 1.0 before the fix and, for a segment the shift changed
    (wholly in a PAR), at L after it, when that value is >= 0.3 in magnitude; a segment the shift left unchanged
    (spanning a PAR boundary) has its genes reported at L + 1.0 both times.

    A row present after the fix is
      spanned    it is present unchanged before the fix (each before-row explains one after-row) and equals L + 1.0,
                 a level genemetrics reports, for an unchanged segment (the documented limitation);
      corrected  it equals L, a level genemetrics reports, for a changed segment;
      otherwise unexplained. A corrected row with no before-row (below) is explained only when |L + 1.0| < 0.3.
    A row present before the fix and not consumed as a spanned row must equal L + 1.0, a reported level, for a
    changed segment, and either be replaced by a corrected after-row with the same gene, coordinates and fields
    other than log2, corrected by that same segment (so a same-name row elsewhere does not count), or have
    |L| < 0.3 (dropped below the threshold once the +1 is gone); otherwise it is unexplained.
    A row's segments are row_segments: one normally, several where segments overlap; the row is explained when any
    of them explains it (Z1: all of them would reject a valid row whose overlapping segments differ in log2, as
    genemetrics reports it once per segment, each at that segment's log2).
    'Equals' is prints_as: the row is CNVkit's %.6g print of that level, within half a unit in its 6th significant
    digit; the levels are computed as CNVkit computes them (before_level, after_level), and a level within float
    noise of the 0.3 threshold may be reported or not.
    Returns ([gene@segment spanned], [gene@coordinates unexplained])."""
    spanned, unexplained = [], []
    before = collections.Counter(p0)           # before-rows not yet consumed by an unchanged spanned after-row
    corrected = collections.defaultdict(list)  # without_log2(row) -> per corrected after-row, the segments fitting it
    for l in p1:
        c = l.split('\t'); s, e, v = int(c[2]), int(c[3]), float(c[4])
        cand = row_segments(segs, s, e)
        unshifted = [x for x in cand if not shifted(x) and prints_as(v, before_level(x[2]))
                     and may_be_reported(before_level(x[2]))]
        fixed = [x for x in cand if shifted(x) and prints_as(v, after_level(x[2]))
                 and may_be_reported(after_level(x[2]))]
        if before[l] and unshifted:
            before[l] -= 1
            spanned.append(f'{c[0]}@chrX:{unshifted[0][0]}-{unshifted[0][1]}')
        elif fixed:
            corrected[without_log2(c)].append(fixed)
        else:
            unexplained.append(f'{c[0]}@chrX:{s}-{e}')
    after = set(p1)
    after_keys = {without_log2(l.split('\t')) for l in p1}
    for l in before.elements():
        c = l.split('\t'); s, e, v = int(c[2]), int(c[3]), float(c[4])
        fits = [x for x in row_segments(segs, s, e) if shifted(x) and prints_as(v, before_level(x[2]))
                and may_be_reported(before_level(x[2]))]
        pending = corrected[without_log2(c)]
        match = next((i for i, fx in enumerate(pending) if set(fx) & set(fits)), None)
        if match is not None:
            pending.pop(match)
            continue
        if any(may_be_unreported(after_level(x[2])) for x in fits):
            continue
        # an unchanged row was judged as an after-row already (reported once); a changed one is named as before
        # the fix, and a row with no counterpart after the fix as dropped
        tag = '' if l in after else ' (before the fix)' if without_log2(c) in after_keys else ' (dropped)'
        if f'{c[0]}@chrX:{s}-{e}{tag}' not in unexplained:
            unexplained.append(f'{c[0]}@chrX:{s}-{e}{tag}')
    for key, pending in corrected.items():
        for fx in pending:
            if not any(may_be_unreported(before_level(x[2])) for x in fx):   # reported before, yet not in the original
                unexplained.append(f'{key[0]}@chrX:{key[2]}-{key[3]} (not in the original)')
    return spanned, unexplained


def compare(orig, new, call_cns):
    """genemetrics_parfix.compare for one sample; returns the statistics (raises CheckFailed when the headers or the
    non-PAR rows differ, and ValueError when a genemetrics table cannot be read)."""
    rec = {}
    rec['par_segs'], rec['span_segs'] = par_counts(call_cns)
    h0, r0 = read_genemetrics(orig); h1, r1 = read_genemetrics(new)
    # CNVkit writes its 5-column header when genemetrics finds no gene (empty table), before or after the fix;
    # otherwise the headers must match (explicit checks, not assert: python -O must not skip them)
    if not (h0 == h1 or (not r1 and h1 == EMPTY_GM_HEADER) or (not r0 and h0 == EMPTY_GM_HEADER)):
        raise CheckFailed(f'headers differ: {h0!r} vs {h1!r}')
    rec['empty_after'] = int(not r1)
    isp = lambda l: (lambda c: par_state(c[1], int(c[2]), int(c[3])) != 'out')(l.split('\t'))
    n0 = [l for l in r0 if not isp(l)]; n1 = [l for l in r1 if not isp(l)]
    p0 = [l for l in r0 if isp(l)]; p1 = [l for l in r1 if isp(l)]
    rec['nonpar_identical'] = int(n0 == n1)
    if n0 != n1:
        raise CheckFailed('non-PAR rows differ')
    rec['par_rows_before'] = len(p0); rec['par_rows_after'] = len(p1)
    rec['rows_changed'] = len(collections.Counter(p0) - collections.Counter(p1)) + len(collections.Counter(p1) - collections.Counter(p0))
    rec['file_changed'] = int(r0 != r1)
    rec['par_before'] = ';'.join(sorted({f"{l.split(chr(9))[0]}:{l.split(chr(9))[4]}" for l in p0}))
    rec['par_after'] = ';'.join(sorted({f"{l.split(chr(9))[0]}:{l.split(chr(9))[4]}" for l in p1}))
    # Each PAR row is judged against its own segment (par_row_verdicts). par_check is 0 when any row is unexplained,
    # 'span' when a row is spanned, or a spanning segment is present and the PAR rows are unchanged (the
    # validation's label: on the 60 validated treated males per set every label is unchanged), and 1 otherwise.
    spanned, unexplained = par_row_verdicts(p0, p1, seg_lookup(call_cns))
    if unexplained:
        rec['unexplained'] = ';'.join(unexplained)
        rec['par_check'] = 0
    elif spanned or (rec['span_segs'] and p0 == p1):
        rec['span_genes'] = ';'.join(spanned) or '-'
        rec['par_check'] = 'span'
    else:
        rec['par_check'] = 1
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
            stats = ' '.join(f'{k}={v}' for k, v in rec.items())
            if rec['par_check'] == 0:
                # an unexplained PAR discrepancy: the corrected genemetrics is not trusted and the job stops
                raise CheckFailed('par_check 0: PAR gene row(s) of the post-step genemetrics do not equal their '
                                  'call.cns segment log2 and are not explained by a segment spanning a PAR '
                                  f"boundary ({rec['unexplained']}); {stats}")
            print(stats)
    except ValueError as e:     # CheckFailed included
        sys.exit(f'cnvkit_parfix.py {a.cmd}: FAILED: {e}')


if __name__ == '__main__':
    main()
