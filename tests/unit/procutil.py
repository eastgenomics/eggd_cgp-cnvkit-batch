"""Bounded subprocesses for the tests.

Every child process the tests start has a timeout. A command that can start descendants (bash running the app under
stubs, `docker run`) runs in a new session, so on timeout the whole process group is killed, not only its leader;
`on_timeout` then removes anything the group owned outside it (e.g. a named Docker container).
"""
import os
import signal
import subprocess

SCRIPT_TIMEOUT = 60     # seconds: one bundled script on a small synthetic table
GRACE = 10              # seconds to collect output after the group was killed


def kill_group(proc):
    """SIGKILL the process group led by proc (started with start_new_session=True)."""
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run_group(cmd, timeout, on_timeout=None, check=False, **kwargs):
    """subprocess.run(cmd, timeout=timeout, ...) in a new session. On timeout: kill the process group, call
    on_timeout(), then raise subprocess.TimeoutExpired (with any output collected)."""
    if kwargs.pop("capture_output", False):
        kwargs["stdout"] = kwargs["stderr"] = subprocess.PIPE
    stdin = kwargs.pop("input", None)
    if stdin is not None:
        kwargs["stdin"] = subprocess.PIPE
    with subprocess.Popen(cmd, start_new_session=True, **kwargs) as proc:
        try:
            out, err = proc.communicate(stdin, timeout=timeout)
        except subprocess.TimeoutExpired:
            kill_group(proc)
            try:
                out, err = proc.communicate(timeout=GRACE)
            except subprocess.TimeoutExpired:   # a descendant left the group and holds the pipes open
                out = err = None
            if on_timeout:
                on_timeout()
            raise subprocess.TimeoutExpired(cmd, timeout, output=out, stderr=err) from None
    result = subprocess.CompletedProcess(cmd, proc.returncode, out, err)
    if check:
        result.check_returncode()
    return result


def remove_container(name, timeout=60):
    """`docker rm --force NAME`, bounded; for a named container left by a timed-out `docker run`."""
    try:
        subprocess.run(["docker", "rm", "--force", name], capture_output=True, timeout=timeout)
    except (subprocess.SubprocessError, OSError):
        pass
