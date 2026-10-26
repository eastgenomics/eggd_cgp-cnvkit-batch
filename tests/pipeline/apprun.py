"""Run the app's `main` under PATH stubs for dx, dx-jobutil-add-output, docker and nproc.

The job environment is reproduced as DNAnexus sets it up: each input is an environment variable, file inputs
hold a `{"$dnanexus_link": ...}` JSON value, and inputs not supplied take their `dxapp.json` default (or stay
unset). Nothing touches DNAnexus or Docker.
"""
import dataclasses
import json
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
STUBS = Path(__file__).resolve().parent / "stubs"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
DXAPP = json.loads((REPO / "dxapp.json").read_text())

SAMPLE = "140609591-25328S0023-26CGPV2-8475-M-96770770"
FILE_IDS = {"tumour_bam": "file-BAM", "tumour_bai": "file-BAI", "baits": "file-BAITS", "cn_reference": "file-PON"}


@dataclasses.dataclass
class AppRun:
    returncode: int
    stdout: str
    stderr: str
    workdir: Path
    docker: list      # every `docker run`: {"options", "image", "cmd"}
    dx: list          # every dx argv
    outputs: dict     # job output name -> {"value", "class"}
    uploads: list     # uploaded file names

    @property
    def cnvkit(self):
        """CNVkit argv lists (without the leading `cnvkit.py`), in call order."""
        return [d["cmd"][1:] for d in self.docker if d["cmd"][:1] == ["cnvkit.py"]]

    def subcommand(self, name):
        """The argv lists of every call of one CNVkit subcommand."""
        return [c for c in self.cnvkit if c[0] == name]


def default_inputs(sample_id=SAMPLE, bam_name=None):
    bam_name = bam_name or f"{sample_id}.bam"
    inputs = {name: {"$dnanexus_link": fid} for name, fid in FILE_IDS.items()}
    inputs["sample_id"] = sample_id
    names = {"file-BAM": bam_name, "file-BAI": bam_name + ".bai"}
    return inputs, names


def _env_value(value):
    if isinstance(value, dict):
        return json.dumps(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def run_app(tmp_path, inputs=None, names=None, fail=(), files=None, script="src/code.sh"):
    """Run `main` from `script` in a fresh work directory under tmp_path.

    inputs: input name -> value (None removes an input); names: file ID -> name returned by `dx describe`;
    fail: CNVkit subcommands that should exit 1; files: file ID -> local fixture copied on download.
    """
    base_inputs, base_names = default_inputs()
    if inputs:
        base_inputs.update(inputs)
    if names:
        base_names.update(names)
    for spec in DXAPP["inputSpec"]:   # platform defaults for inputs not supplied
        if spec["name"] not in base_inputs and "default" in spec:
            base_inputs[spec["name"]] = spec["default"]
    base_inputs = {k: v for k, v in base_inputs.items() if v is not None}

    work, log = tmp_path / "work", tmp_path / "log"
    work.mkdir()
    log.mkdir()
    dx_files = {"file-BAITS": str(FIXTURES / "panel.bed")}
    dx_files.update(files or {})
    env = {
        "PATH": f"{STUBS}:{REPO / 'resources' / 'usr' / 'local' / 'bin'}:{os.environ['PATH']}",
        "HOME": str(tmp_path),
        "LC_ALL": "C.UTF-8",
        "STUB_LOG": str(log),
        "STUB_DX_FILES": json.dumps(dx_files),
        "STUB_DX_NAMES": json.dumps(base_names),
        "STUB_FAIL": ",".join(fail),
        # on the worker, resources/ is unpacked to /; here the library is read from the repo
        "CGP_CNVKIT_LIB": str(REPO / "resources" / "usr" / "local" / "lib" / "cgp_cnvkit_batch"),
    }
    env.update({k: _env_value(v) for k, v in base_inputs.items()})
    proc = subprocess.run(
        ["bash", "-c", 'set -eo pipefail; source "$1"; main', "app", str(REPO / script)],
        cwd=work, env=env, capture_output=True, text=True, timeout=120,
    )

    def jsonl(name):
        p = log / name
        return [json.loads(line) for line in p.read_text().splitlines()] if p.exists() else []

    return AppRun(
        returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr, workdir=work,
        docker=jsonl("docker.jsonl"), dx=jsonl("dx.jsonl"),
        outputs={o["name"]: {"value": o["value"], "class": o["class"]} for o in jsonl("outputs.jsonl")},
        uploads=[u["name"] for u in jsonl("uploads.jsonl")],
    )
