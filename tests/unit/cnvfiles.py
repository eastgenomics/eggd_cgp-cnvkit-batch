"""Builders for small synthetic CNVkit tables (.cnr, .targetcoverage.cnn, .cns) used by the unit tests."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "resources" / "usr" / "local" / "bin"
REFERENCE = REPO / "tests" / "validation" / "reference_scripts"

CNR_HEADER = ["chromosome", "start", "end", "gene", "log2", "depth", "weight"]
CNN_HEADER = ["chromosome", "start", "end", "gene", "depth", "log2"]


def write_table(path, header, rows):
    path = Path(path)
    path.write_text("\n".join(["\t".join(header)] + ["\t".join(str(x) for x in r) for r in rows]) + "\n")
    return path


def bins(spec, chrom="chr1", start=1000, length=144, gap=100, gene="G"):
    """Consecutive bins from a list of log2 values: returns [(chrom, start, end, gene, log2), ...]."""
    out, pos = [], start
    for i, l2 in enumerate(spec):
        out.append((chrom, pos, pos + length, f"{gene}{i}", l2))
        pos += length + gap
    return out


def cnr_rows(bin_list, weight="0.9"):
    return [(c, s, e, g, l2, "100.0", weight) for c, s, e, g, l2 in bin_list]


def cnn_rows(bin_list, depth):
    """depth: a number for every bin, or {(chrom, start): depth} with 500 for the others (raw depth column)."""
    get = (lambda c, s: depth.get((c, s), 500.0)) if isinstance(depth, dict) else (lambda c, s: depth)
    return [(c, s, e, g, get(c, s), "9.0") for c, s, e, g, _ in bin_list]
