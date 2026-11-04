"""Bounded subprocesses for the tests.

Every child process the tests and the test helpers start has a timeout, with one deferred exception: the L4
harness's `--start-from bam` BAM/BAI transfers and their `dx describe` / `dx find` calls (run_local.fetch_bam_inputs)
are unbounded (https://github.com/eastgenomics/eggd_cgp-cnvkit-batch/issues/5). A command that can start descendants
(bash running the app under stubs, `docker`, `dx`, `git`) runs in a new session (run_group; the L4 harness's pipeline
run does the same itself), so on timeout the whole process
group is killed, not only its leader; `cleanup` then removes anything the group owned outside it (e.g. a named Docker
container). The same happens when waiting is interrupted (KeyboardInterrupt or any other exception), which is then
re-raised unchanged. A bundled Python script with no descendants runs under plain subprocess.run(timeout=...).

Every timeout path is bounded and always attempts the cleanup: a failed group kill falls back to killing the leader,
output is collected for at most GRACE seconds, and `cleanup` runs even when the kill failed. The exception raised is
the original one (the TimeoutExpired from the wait, with the output captured so far, or the interrupting one).
Anything that could not be confirmed (kill failure, a process still holding the pipes, a failed wait or output
collection, a failed cleanup) is reported as a note on it (also in its `cleanup_notes`), never in place of it; a
pipe that fails to close after the kill is a note too (in `close_notes`; round 7), but not unconfirmed cleanup. A
caller that reports the failure in its own words keeps all the notes with failure_text(). Any exception from a cleanup
step counts, not only RuntimeError or OSError, and text output is decoded with errors='replace', so invalid UTF-8
from a child (e.g. Docker) cannot raise; remove_containers attempts every name whatever an earlier removal raised
(round 6, V2).
"""
import contextlib
import os
import signal
import subprocess

SCRIPT_TIMEOUT = 60     # seconds: one bundled script on a small synthetic table
GRACE = 10              # seconds to collect output after the group was killed
DOCKER_RM_TIMEOUT = 60  # seconds for `docker rm --force`


def kill_group(proc):
    """SIGKILL the process group led by proc (started with start_new_session=True); if that fails, SIGKILL the
    leader. Returns the problems met (empty when the kill was sent or the group had already exited)."""
    try:
        os.killpg(proc.pid, signal.SIGKILL)
        return []
    except ProcessLookupError:
        return []
    except OSError as e:
        problems = [f"killpg({proc.pid}) failed: {e!r}"]
    try:
        proc.kill()
    except ProcessLookupError:
        pass
    except OSError as e:
        problems.append(f"kill({proc.pid}) failed: {e!r}")
    return problems


def _kill(proc, notes):
    """kill_group(proc), its problems added to notes; an unexpected exception is a note too, never raised in place of
    the exception in flight (round 6, V2)."""
    try:
        notes += kill_group(proc)
    except Exception as e:  # noqa: BLE001 -- reported on the exception in flight, which must not be masked
        notes.append(f"killing process group {proc.pid} failed: {e!r}")


def _close_pipes(proc, exc):
    """Close proc's pipes after the kill (best effort). A failure to close is recorded on exc as a note and in its
    `close_notes` (round 7, W3), never raised in place of it, and is not unconfirmed cleanup (cleanup_unconfirmed
    ignores it): the process is already killed, or a `cleanup_notes` entry says it may not be, so a close failure
    alone must not turn a Docker-check skip into a failure."""
    notes = []
    for f in (proc.stdin, proc.stdout, proc.stderr):
        if f is not None:
            try:
                f.close()
            except Exception as e:  # noqa: BLE001 -- recorded on the exception in flight, which must not be masked
                notes.append(f"closing a pipe of process {proc.pid} failed: {e!r}")
    exc.close_notes = notes
    for note in notes:
        exc.add_note(note)


def _run_cleanup(cleanup, notes):
    """Call cleanup() (if any); a failure is added to notes, never raised in place of the exception in flight."""
    if cleanup is not None:
        try:
            cleanup()
        except Exception as e:  # noqa: BLE001 -- reported on the exception in flight, which must not be masked
            notes.append(f"cleanup unconfirmed: {e if isinstance(e, RuntimeError) else repr(e)}")


def _annotate(exc, notes):
    exc.cleanup_notes = notes
    for note in notes:
        exc.add_note(note)


def details(exc, tail=500):
    """What a failed command's exception carries beyond str(exc), for a caller that reports it in its own words: a
    warning when its cleanup is unconfirmed (cleanup_unconfirmed), the tail of any stderr it captured (TimeoutExpired,
    CalledProcessError) and every note run_group attached (`cleanup_notes`: an unconfirmed kill, wait, output
    collection or cleanup; `close_notes`: a pipe that failed to close)."""
    out = [UNCONFIRMED] if cleanup_unconfirmed(exc) else []
    err = getattr(exc, "stderr", None)
    if isinstance(err, bytes):
        err = err.decode(errors="replace")
    if err and err.strip():
        out.append(f"stderr: {err.strip()[-tail:]}")
    return out + list(getattr(exc, "cleanup_notes", None) or []) + list(getattr(exc, "close_notes", None) or [])


UNCONFIRMED = "cleanup unconfirmed: a child process (or its container) may remain"


def timed_out(what, timeout, exc):
    """The message for a command that timed out: '<what> timed out after <timeout> s', then its details()."""
    return "; ".join([f"{what} timed out after {timeout} s", *details(exc)])


def failure_text(exc):
    """str(exc) with its details(): no diagnostic or cleanup note is dropped."""
    return "; ".join([str(exc), *details(exc)])


def cleanup_unconfirmed(exc):
    """True when run_group could not confirm the cleanup after exc: a process (or a container) may remain. A pipe
    that failed to close (`close_notes`) does not count."""
    return bool(getattr(exc, "cleanup_notes", None))


def run_group(cmd, timeout, cleanup=None, check=False, **kwargs):
    """subprocess.run(cmd, timeout=timeout, ...) in a new session. On timeout: kill the process group, collect output
    for at most GRACE s, call cleanup() whatever happened before, then re-raise the original subprocess.TimeoutExpired
    (its output: everything collected, or what it captured before the timeout when the collection fails; a note for
    each unconfirmed step, the notes also in its `cleanup_notes`).

    On any other exception while waiting (KeyboardInterrupt, ...): kill the process group, wait at most GRACE s for
    the leader, call cleanup() the same way, and re-raise that exception unchanged, with the same notes added.

    Text output is decoded with errors='replace' unless the caller sets `errors`, so output that is not valid in the
    encoding (e.g. invalid UTF-8 from Docker) can never raise UnicodeDecodeError (round 6, V2)."""
    if kwargs.get("text") or kwargs.get("universal_newlines") or kwargs.get("encoding"):
        kwargs.setdefault("errors", "replace")
    if kwargs.pop("capture_output", False):
        kwargs["stdout"] = kwargs["stderr"] = subprocess.PIPE
    stdin = kwargs.pop("input", None)
    if stdin is not None:
        kwargs["stdin"] = subprocess.PIPE
    proc = subprocess.Popen(cmd, start_new_session=True, **kwargs)
    try:
        out, err = proc.communicate(stdin, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        # the original exception is raised, with the output it captured before the timeout (bytes, as communicate
        # gives it) unless the bounded collection after the kill gets more (round 5, U2)
        notes = []
        try:
            _kill(proc, notes)
            try:
                exc.output, exc.stderr = proc.communicate(timeout=GRACE)
            except subprocess.TimeoutExpired as late:   # it holds everything read so far, the earlier part included
                notes.append(f"process {proc.pid} (or a descendant holding its output) still running after the kill "
                             f"and {GRACE} s")
                exc.output = late.output if late.output is not None else exc.output
                exc.stderr = late.stderr if late.stderr is not None else exc.stderr
            except Exception as e:  # noqa: BLE001 -- noted on the TimeoutExpired raised below, never in its place
                notes.append(f"collecting the output of process {proc.pid} after the kill failed: {e!r}")
        finally:
            _run_cleanup(cleanup, notes)
            _annotate(exc, notes)
            _close_pipes(proc, exc)
        raise
    except BaseException as exc:    # e.g. KeyboardInterrupt: do not leave the group, or its container, running
        notes = []
        try:
            _kill(proc, notes)
            try:
                proc.wait(timeout=GRACE)
            except subprocess.TimeoutExpired:
                notes.append(f"process {proc.pid} still running after the kill and {GRACE} s")
            except Exception as e:  # noqa: BLE001 -- noted on the interrupting exception, never in its place
                notes.append(f"waiting for process {proc.pid} after the kill failed: {e!r}")
        finally:
            _run_cleanup(cleanup, notes)
            _annotate(exc, notes)
            _close_pipes(proc, exc)
        raise
    result = subprocess.CompletedProcess(cmd, proc.returncode, out, err)
    if check:
        result.check_returncode()
    return result


def remove_container(name, timeout=DOCKER_RM_TIMEOUT):
    """`docker rm --force NAME`, bounded, for a named container left by a timed-out `docker run`. Returns when the
    container is removed or absent; raises RuntimeError, and only RuntimeError, when that cannot be confirmed (any
    failure to run docker, not only OSError; round 6, V2). Docker's output is decoded with errors='replace'
    (run_group), so invalid UTF-8 in it is reported, not raised."""
    try:   # bounded in its own process group, as the docker CLI (or a wrapper) can start children
        proc = run_group(["docker", "rm", "--force", name], timeout, capture_output=True, text=True)
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(timed_out(f"docker rm --force {name}:", timeout, e)) from None
    except Exception as e:  # noqa: BLE001 -- OSError and anything else: the removal is unconfirmed, reported as such
        raise RuntimeError(f"docker rm --force {name}: {e!r}") from None
    if proc.returncode != 0 and "No such container" not in proc.stderr:
        raise RuntimeError(f"docker rm --force {name}: exit {proc.returncode}: {proc.stderr.strip()[-500:]}")


def remove_containers(names, timeout=DOCKER_RM_TIMEOUT):
    """remove_container for each name, every one attempted; returns the problems (empty when every container is
    confirmed removed or absent). Any exception from a removal, not only RuntimeError, is a problem reported for that
    name, never raised, so the remaining names are still removed (round 6, V2)."""
    problems = []
    for name in names:
        try:
            remove_container(name, timeout)
        except RuntimeError as e:
            problems.append(str(e))
        except Exception as e:  # noqa: BLE001 -- an unconfirmed removal is reported; the other names still run
            problems.append(f"docker rm --force {name}: {e!r}")
    return problems


@contextlib.contextmanager
def removing_containers(*names, timeout=DOCKER_RM_TIMEOUT):
    """Remove the named containers (remove_containers) when the block exits. If the block raised, its exception is
    re-raised unchanged, with each unconfirmed removal as a note; if it did not, an unconfirmed removal raises
    RuntimeError, so the block cannot pass without confirming that its containers are gone."""
    try:
        yield
    except BaseException as exc:
        for p in remove_containers(names, timeout):
            exc.add_note(f"container removal unconfirmed: {p}")
        raise
    problems = remove_containers(names, timeout)
    if problems:
        raise RuntimeError("container removal unconfirmed: " + "; ".join(problems))
