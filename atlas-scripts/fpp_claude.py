#!/usr/bin/python3
"""
fpp_claude.py -- parallel driver for computing the FPP unitary dual with atlas.

Rewrite of fpp.py.  Runs a *fixed* pool of N persistent atlas workers, each
owning one long-lived `atlas` subprocess.  Workers pull (x,lambda) pair indices
from a shared queue, compute the FPP-unitary facets for each pair, and append
the results to their own data file.  A memory governor keeps total RSS under a
cap by asking workers to recycle (exit atlas cleanly between pairs and relaunch
fresh); only an imminent-OOM emergency ever kills a worker mid-computation.

Key design points (see the prompt fpp_prompt.txt):
  * Fixed worker count (E7_s: 500, F4_s test: 10).
  * Global RSS cap (E7_s: 12 TB) enforced by *graceful* recycling.
  * Robust data flow: a unique sentinel terminates every command (no more
    guessing with "END"), big diagnostic output is redirected to a file by
    atlas so the pipe only carries small status lines, and a select()-based
    reader detects hangs.
  * The shared init file (big_unitary_hash) is reset from a reference at start
    and periodically rebuilt by merging every worker's data file back in.
  * Jobs show up in `top` as atlas_<job_id> via per-job symlinks to the atlas
    executable (argv[0]/comm trick).

Typical use:
    ./fpp_claude.py --preset f4s          # 10 workers, tiny caps, test run
    ./fpp_claude.py --preset e7s          # 500 workers, 12 TB, the real thing
    ./fpp_claude.py --preset e7s -n 200 --max-total-gb 8000   # overrides

Stop with Ctrl-C (SIGINT): workers drain gracefully and a final init merge runs.
"""

import argparse
import datetime
import glob
import os
import re
import select
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from queue import Queue, Empty

try:
    import psutil
except ImportError:
    sys.exit("fpp_claude.py requires psutil (pip install --user psutil)")

# --------------------------------------------------------------------------- #
# Fixed locations
# --------------------------------------------------------------------------- #
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
ATLAS_DIR = os.path.dirname(SCRIPTS_DIR)            # .../to_ht_branch_jeff_2
ATLAS_EXE = os.path.join(ATLAS_DIR, "atlas")
DEFAULT_OUTPUT_ROOT = "/scr/jdada11/fpp_claude"
GROUP_SYMBOL = "G_temp"                             # name the init/data files use

GB = 1024.0 ** 3


# --------------------------------------------------------------------------- #
# Presets
# --------------------------------------------------------------------------- #
@dataclass
class Preset:
    name: str
    workers: int
    max_total_gb: float
    max_proc_gb: float
    reference_init: str
    aux_files: list
    # hang detection: a COMPUTING worker pinned near 0% CPU this long is hung
    stall_seconds: float


PRESETS = {
    "f4s": Preset(
        name="f4s", workers=10,
        max_total_gb=3.0, max_proc_gb=0.3,
        reference_init="f4sinitreference.at", aux_files=[],
        stall_seconds=120.0,
    ),
    "e7s": Preset(
        name="e7s", workers=500,
        max_total_gb=12000.0, max_proc_gb=40.0,
        reference_init="e7sinitreference.at",
        aux_files=["edges_F4_E6_E7.at", "coh_ind_E7_centered2PlusE6.at"],
        stall_seconds=1800.0,
    ),
}


@dataclass
class Config:
    name: str
    workers: int
    max_total_gb: float
    max_proc_gb: float
    reference_init: str       # absolute path
    init_file: str            # absolute path (in SCRIPTS_DIR)
    aux_files: list           # absolute paths
    output_root: str
    stall_seconds: float
    reverse: bool = False
    limit: int = -1
    poll_interval: float = 5.0
    merge_interval: float = 120.0
    startup_timeout: float = 900.0      # loading all.at + init + aux (big for E7)
    cmd_timeout: float = 300.0          # short commands (is_finished, write)
    compute_timeout: float = 30 * 86400  # effectively "no timeout"; governor kills hangs
    keep_init: bool = False             # if True, don't reset init from reference
    settings_file: str = ""             # loaded last (overrides defaults); "" = none

    # filled in during setup
    run_dir: str = ""
    logs_dir: str = ""
    symlinks_dir: str = ""
    round: int = 0
    definition_file: str = ""

    @property
    def hard_total_gb(self):
        return self.max_total_gb * 1.10   # emergency ceiling


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def fmt_dur(seconds):
    return re.sub(r"\..*", "", str(datetime.timedelta(seconds=int(seconds))))


class HangError(Exception):
    """A command did not return its sentinel within the timeout."""


class WorkerDied(Exception):
    """The atlas subprocess closed its stdout (died / was killed)."""


class LineReader:
    """Read newline-terminated lines from a raw fd with a per-read timeout.

    Returns the line (bytes, incl. newline) on success, None on timeout, and
    b"" on EOF.  Uses os.read so poll() readiness is accurate (no hidden
    buffering as with TextIOWrapper).

    Uses poll() rather than select(): with hundreds of workers the atlas pipe
    fds exceed select()'s FD_SETSIZE (1024) limit, which raises
    "filedescriptor out of range in select()".  poll() has no such limit.
    Each poll wait is capped (so a multi-day timeout doesn't overflow poll's
    millisecond int) and looped against an overall deadline."""

    _MAX_WAIT_MS = 3600 * 1000   # cap a single poll() wait at 1 h

    def __init__(self, fd):
        self.fd = fd
        self.buf = b""
        self.poller = select.poll()
        self.poller.register(fd, select.POLLIN | select.POLLHUP | select.POLLERR)

    def readline(self, timeout):
        deadline = None if timeout is None else time.time() + timeout
        while b"\n" not in self.buf:
            if deadline is not None:
                remaining = deadline - time.time()
                if remaining <= 0:
                    return None
                wait_ms = min(int(remaining * 1000) + 1, self._MAX_WAIT_MS)
            else:
                wait_ms = self._MAX_WAIT_MS
            if not self.poller.poll(wait_ms):
                continue                     # capped wait elapsed; deadline rechecked
            chunk = os.read(self.fd, 65536)
            if chunk == b"":
                if self.buf:
                    line, self.buf = self.buf, b""
                    return line
                return b""
            self.buf += chunk
        line, self.buf = self.buf.split(b"\n", 1)
        return line + b"\n"


# --------------------------------------------------------------------------- #
# Atlas worker
# --------------------------------------------------------------------------- #
# states
IDLE, COMPUTING, RESTARTING, STOPPED = "idle", "computing", "restarting", "stopped"


class AtlasWorker:
    def __init__(self, job, cfg):
        self.job = job
        self.cfg = cfg
        self.proc = None
        self.reader = None
        self.seq = 0
        self.symlink = os.path.join(cfg.symlinks_dir, "atlas_%d" % job)
        self.data_file = os.path.join(cfg.run_dir, "%d.at" % job)
        # One per-job log holds both our structured lines and the atlas process's
        # (verbose) computation output, which atlas appends directly via redirect.
        # Both writers use O_APPEND and never write at the same instant (the
        # worker is single-threaded), so they interleave cleanly in order.
        self.log_path = os.path.join(cfg.logs_dir, "%d.log" % job)
        self.logfh = open(self.log_path, "a", buffering=1)
        self.stderr_file = os.path.join(cfg.logs_dir, "%d.stderr" % job)
        # cross-thread signals from the governor
        self.recycle_event = threading.Event()    # graceful recycle requested
        # governor bookkeeping (only the governor writes these)
        self.state = STOPPED
        self.compute_started = 0.0
        self.low_cpu_since = 0.0
        # ensure the symlink that gives the atlas_<job> name in top exists
        if not os.path.islink(self.symlink):
            try:
                os.symlink(ATLAS_EXE, self.symlink)
            except FileExistsError:
                pass

    def log(self, msg):
        self.logfh.write("[%s] job %d: %s\n" % (now(), self.job, msg))
        self.logfh.flush()

    # -- process lifecycle ------------------------------------------------- #
    def launch(self):
        args = [self.symlink, "all.at", self.cfg.init_file] + self.cfg.aux_files
        if self.cfg.settings_file:
            args.append(self.cfg.settings_file)
        stderr = open(self.stderr_file, "ab")
        self.proc = subprocess.Popen(
            args, executable=self.symlink, cwd=SCRIPTS_DIR,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=stderr,
            bufsize=0,
        )
        self.reader = LineReader(self.proc.stdout.fileno())
        self.seq = 0
        # Drain the (large) startup banner from loading all.at + init + aux, so
        # later command output isn't mixed with it.
        self._command("", self.cfg.startup_timeout)
        self.log("launched pid=%d (%s)" % (self.proc.pid, os.path.basename(self.symlink)))

    @property
    def pid(self):
        return self.proc.pid if self.proc else None

    def _write(self, s):
        self.proc.stdin.write(s.encode("utf-8"))
        self.proc.stdin.flush()

    def _command(self, atlas_cmd, timeout):
        """Run an atlas command, return the stdout lines emitted before the
        sentinel.  Diagnostic output that the command redirects to a file does
        not appear here -- only what the command prints straight to stdout."""
        if self.proc.poll() is not None:
            raise WorkerDied()
        self.seq += 1
        token = "@@SENTINEL_%d_%d@@" % (self.job, self.seq)
        if atlas_cmd:
            self._write(atlas_cmd + "\n")
        self._write('prints("%s")\n' % token)
        deadline = time.time() + timeout
        out = []
        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                raise HangError("no sentinel for: %s" % (atlas_cmd or "<prints>"))
            line = self.reader.readline(remaining)
            if line is None:
                raise HangError("timeout on: %s" % (atlas_cmd or "<prints>"))
            if line == b"":
                raise WorkerDied()
            s = line.decode("latin-1", "replace").rstrip("\n")
            if token in s:
                return out
            out.append(s)

    def quit(self):
        """Cleanly exit atlas and reap the process."""
        try:
            if self.proc and self.proc.poll() is None:
                self._write("quit\n")
        except (BrokenPipeError, OSError):
            pass
        self._reap()

    def kill(self):
        try:
            if self.proc and self.proc.poll() is None:
                self.proc.kill()
        except OSError:
            pass
        self._reap()

    def _reap(self):
        if not self.proc:
            return
        for stream in (self.proc.stdin, self.proc.stdout):
            try:
                stream.close()
            except OSError:
                pass
        try:
            self.proc.wait(timeout=30)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            self.proc.wait()
        self.proc = None
        self.reader = None

    # -- the actual FPP work ----------------------------------------------- #
    def is_finished(self, k):
        out = self._command("prints(is_finished(%s,%d))" % (GROUP_SYMBOL, k),
                             self.cfg.cmd_timeout)
        for s in out:
            t = s.strip()
            if t == "true":
                return True
            if t == "false":
                return False
        # Unexpected output -- be safe and (re)compute it.
        self.log("is_finished(%d): unexpected output %r" % (k, out))
        return False

    def compute(self, k):
        # The debugging flags (FPP_report_flag, every_lambda_deets_flag, ...)
        # make the computation very verbose.  Append that output straight to the
        # job log via atlas' file redirect (real newlines, and it keeps the
        # control pipe light).  A header line (flushed first) keeps it readable.
        self.log("---- pair %d: FPP_unitary_hash_bottom_layer output ----" % k)
        cmd = '>>"%s" FPP_unitary_hash_bottom_layer(xl_pair(%s,%d))' % (
            self.log_path, GROUP_SYMBOL, k)
        self._command(cmd, self.cfg.compute_timeout)

    def write_pair(self, k, elapsed_s):
        # write_one_pair returns the data string: the unitary block (long_match
        # / finish_num lines) joined by *literal* "\n", then a "void:add(...)"
        # report record, then a trailing "END".  We capture it over stdout and
        # turn the literal "\n" back into real newlines (atlas' file redirect
        # would write them literally, leaving an unloadable file), dropping the
        # "END" marker.  Return bytes written; 0 means atlas errored (the
        # sentinel still returns on failure) so the caller can react.
        out = self._command(
            "prints(write_one_pair(%s,%d,%d,%d))" % (
                GROUP_SYMBOL, k, self.job, int(elapsed_s)),
            self.cfg.cmd_timeout)
        written = 0
        with open(self.data_file, "a") as f:
            for s in out:
                if s.strip() == "END":
                    continue
                s = s.replace("\\n", "\n") + "\n"
                f.write(s)
                written += len(s)
        return written


# --------------------------------------------------------------------------- #
# Init merge (rebuild the shared init file from reference + all worker data)
# --------------------------------------------------------------------------- #
class Merger:
    def __init__(self, cfg, log):
        self.cfg = cfg
        self.log = log
        self.lock = threading.Lock()
        self.request = threading.Event()
        self.last_merge = 0.0
        self.count = 0

    def merge(self, blocking=False):
        """Rebuild cfg.init_file = write(reference_init + grep(all data files)).

        Idempotent: long_match/finish lines only set bits / dedupe in the hash,
        so re-merging every data file each time is safe.  Written to a temp file
        in SCRIPTS_DIR and atomically renamed over the live init.

        blocking=False (periodic merges): skip if one is already running.
        blocking=True (final merge at shutdown): wait for any in-flight merge,
        then run a definitive merge over the now-complete data."""
        if not self.lock.acquire(blocking=blocking):
            return False                     # a merge is already running
        try:
            t0 = time.time()
            self.count += 1
            grep_file = os.path.join(self.cfg.logs_dir, "merge_grep.at")
            data_files = sorted(glob.glob(os.path.join(self.cfg.run_dir, "*.at")))
            with open(grep_file, "wb") as out:
                if data_files:
                    subprocess.run(
                        ["grep", "-hE", r"big_unitary_hash\.(finish|long)"] + data_files,
                        stdout=out, stderr=subprocess.DEVNULL, check=False)
            tmp = self.cfg.init_file + ".tmp.%d" % os.getpid()
            args = [ATLAS_EXE, "all.at", self.cfg.reference_init, grep_file]
            script = '>"%s" big_unitary_hash.write()\nquit\n' % tmp
            proc = subprocess.Popen(args, cwd=SCRIPTS_DIR, stdin=subprocess.PIPE,
                                    stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            _, err = proc.communicate(script.encode("utf-8"), timeout=3600)
            if proc.returncode != 0 or not os.path.exists(tmp) or os.path.getsize(tmp) == 0:
                self.log("MERGE FAILED (rc=%s): %s" %
                         (proc.returncode, err.decode("latin-1", "replace")[:500]))
                if os.path.exists(tmp):
                    os.remove(tmp)
                return False
            os.replace(tmp, self.cfg.init_file)
            self.last_merge = time.time()
            self.log("merge #%d done in %s (%d data files, init=%d bytes)" % (
                self.count, fmt_dur(time.time() - t0), len(data_files),
                os.path.getsize(self.cfg.init_file)))
            return True
        except Exception as e:                       # noqa: BLE001  (log everything)
            self.log("MERGE ERROR: %r" % e)
            return False
        finally:
            self.lock.release()


# --------------------------------------------------------------------------- #
# The run orchestrator
# --------------------------------------------------------------------------- #
class FPPRun:
    def __init__(self, cfg):
        self.cfg = cfg
        self.queue = Queue()
        self.workers = []
        self.threads = []
        self.stop_event = threading.Event()
        self.fail_counts = {}
        self.fail_lock = threading.Lock()
        self.done_count = 0
        self.done_lock = threading.Lock()
        self.main_log_fh = None
        self.merger = None

    # -- logging ----------------------------------------------------------- #
    def log(self, msg):
        line = "[%s] %s\n" % (now(), msg)
        self.main_log_fh.write(line)
        self.main_log_fh.flush()

    # -- setup ------------------------------------------------------------- #
    def setup(self):
        cfg = self.cfg
        os.makedirs(cfg.output_root, exist_ok=True)
        # round = max existing <name>_<k> + 1
        existing = []
        pat = re.compile(r"^%s_(\d+)$" % re.escape(cfg.name))
        for entry in os.listdir(cfg.output_root):
            m = pat.match(entry)
            if m:
                existing.append(int(m.group(1)))
        cfg.round = max(existing) + 1 if existing else 1
        cfg.run_dir = os.path.join(cfg.output_root, "%s_%d" % (cfg.name, cfg.round))
        cfg.logs_dir = os.path.join(cfg.run_dir, "logs")
        cfg.symlinks_dir = os.path.join(cfg.run_dir, "symlinks")
        for d in (cfg.run_dir, cfg.logs_dir, cfg.symlinks_dir):
            os.makedirs(d, exist_ok=True)

        self.main_log_fh = open(os.path.join(cfg.logs_dir, "main.log"), "a", buffering=1)
        self.log("=" * 70)
        self.log("fpp_claude starting; run dir %s" % cfg.run_dir)
        self.log("name=%s workers=%d max_total=%.1f GB max_proc=%.2f GB" % (
            cfg.name, cfg.workers, cfg.max_total_gb, cfg.max_proc_gb))
        self._snapshot()

        # reset the init file from the reference (overwriting), unless --keep-init
        if cfg.keep_init:
            if not os.path.exists(cfg.init_file):
                shutil.copy(cfg.reference_init, cfg.init_file)
                self.log("init %s missing; seeded from reference" % cfg.init_file)
            else:
                self.log("keeping existing init %s" % cfg.init_file)
        else:
            shutil.copy(cfg.reference_init, cfg.init_file)
            self.log("reset init %s <- %s" % (
                os.path.basename(cfg.init_file), os.path.basename(cfg.reference_init)))

        if cfg.settings_file:
            self.log("settings (loaded last): %s" % cfg.settings_file)
        self._write_definition()
        self._record_flags()
        self.merger = Merger(cfg, self.log)
        self._build_queue()

    def _load_args(self):
        """Files each worker loads after all.at, in order: init, aux, settings."""
        args = [self.cfg.init_file] + self.cfg.aux_files
        if self.cfg.settings_file:
            args.append(self.cfg.settings_file)
        return args

    def _record_flags(self):
        """Record the value of every existing *flag* variable (with the settings
        applied) to logs/flags.txt.  The candidate names are scraped from the
        atlas source, then each is queried -- nonexistent ones simply error out
        in atlas and are skipped, so the recorded list is always accurate."""
        cfg = self.cfg
        grep = subprocess.run(
            ["grep", "-rhoE", r"[A-Za-z_][A-Za-z0-9_]*flag[A-Za-z0-9_]* *:?=",
             "--include=*.at", "."],
            cwd=SCRIPTS_DIR, capture_output=True, text=True)
        names = set(re.findall(
            r"[A-Za-z_][A-Za-z0-9_]*flag[A-Za-z0-9_]*", grep.stdout))
        # also record whatever the settings file assigns, even without "flag" in
        # the name (e.g. fund_face_verbose, seat_belt_on)
        if cfg.settings_file and os.path.exists(cfg.settings_file):
            with open(cfg.settings_file) as f:
                for m in re.finditer(r"(?m)^\s*(?:set\s+)?([A-Za-z_]\w*)\s*:?=", f.read()):
                    names.add(m.group(1))
        names = sorted(names)
        if not names:
            return
        script = "".join('prints("@FLAG@|%s|",%s)\n' % (n, n) for n in names) + "quit\n"
        out, _ = self._atlas_oneshot(self._load_args(), script, timeout=900)
        flags = {}
        for line in out.splitlines():
            if line.startswith("@FLAG@|"):
                _, name, val = line.split("|", 2)
                flags[name] = val.strip()
        path = os.path.join(cfg.logs_dir, "flags.txt")
        with open(path, "w") as f:
            for name in sorted(flags):
                f.write("%s = %s\n" % (name, flags[name]))
        self.log("recorded %d flag values (of %d candidates) -> flags.txt" % (
            len(flags), len(names)))

    def _snapshot(self):
        for f in (__file__, "FPP.at", "fpp_settings.at", self.cfg.reference_init):
            src = f if os.path.isabs(f) else os.path.join(SCRIPTS_DIR, f)
            if os.path.exists(src):
                try:
                    shutil.copy(src, self.cfg.logs_dir)
                except OSError as e:
                    self.log("snapshot %s failed: %r" % (src, e))
        try:
            head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=SCRIPTS_DIR,
                                  capture_output=True, text=True).stdout.strip()
            self.log("git HEAD: %s" % head)
        except OSError:
            pass

    def _atlas_oneshot(self, extra_args, script, timeout=600):
        """Run a short atlas job: load all.at + extra_args, feed `script`, return
        (stdout_text, returncode)."""
        args = [ATLAS_EXE, "all.at"] + extra_args
        proc = subprocess.Popen(args, cwd=SCRIPTS_DIR, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(script.encode("utf-8"), timeout=timeout)
        if err:
            self.log("atlas oneshot stderr: %s" %
                     err.decode("latin-1", "replace")[:500])
        return out.decode("latin-1", "replace"), proc.returncode

    def _write_definition(self):
        """Write the standalone group-definition file <name>.at (defines G_temp).
        The init already defines G_temp; this is the canonical group file kept
        alongside the run for reports/tooling, matching the old convention."""
        cfg = self.cfg
        cfg.definition_file = os.path.join(cfg.output_root, "%s.at" % cfg.name)
        tmp = cfg.definition_file + ".tmp"
        script = '>"%s" write_real_form_plus(%s,"%s")\nquit\n' % (
            tmp, GROUP_SYMBOL, GROUP_SYMBOL)
        _, rc = self._atlas_oneshot([cfg.init_file], script)
        if os.path.exists(tmp) and os.path.getsize(tmp) > 0:
            with open(cfg.definition_file, "w") as f:
                f.write("<groups.at\n")
            with open(tmp) as f:
                with open(cfg.definition_file, "a") as g:
                    g.write(f.read())
            os.remove(tmp)
            self.log("wrote group definition %s" % cfg.definition_file)
        else:
            self.log("WARNING: could not write group definition (rc=%s)" % rc)

    def _build_queue(self):
        cfg = self.cfg
        out, rc = self._atlas_oneshot(
            [cfg.init_file],
            "prints(xl_pairs_todo(big_unitary_hash,%s))\nquit\n" % GROUP_SYMBOL)
        # find the list line: "[a,b,c,...]"
        entries = []
        for line in out.splitlines():
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                body = line[1:-1].strip()
                if body:
                    entries = [int(x) for x in body.split(",")]
                break
        if not entries:
            self.log("no (x,lambda) pairs to do (rc=%s). Output tail:\n%s" %
                     (rc, "\n".join(out.splitlines()[-10:])))
            return
        if cfg.reverse:
            entries.reverse()
        if cfg.limit >= 0:
            entries = entries[:cfg.limit]
        for k in entries:
            self.queue.put(k)
        self.log("queued %d (x,lambda) pairs" % len(entries))

    # -- worker thread ----------------------------------------------------- #
    def _record_fail(self, k):
        """Return True if pair k should be quarantined (failed too many times)."""
        with self.fail_lock:
            self.fail_counts[k] = self.fail_counts.get(k, 0) + 1
            n = self.fail_counts[k]
        if n > 2:
            with open(os.path.join(self.cfg.logs_dir, "poison.txt"), "a") as f:
                f.write("%d\n" % k)
            return True
        return False

    def _safe_launch(self, worker, attempts=3):
        """Launch (or relaunch) a worker's atlas, retrying on startup failure.
        Returns False if we should give up (stopping, or repeated failure)."""
        for i in range(attempts):
            if self.stop_event.is_set():
                return False
            try:
                worker.launch()
                return True
            except (HangError, WorkerDied, OSError) as e:
                worker.log("launch attempt %d/%d failed: %r" % (i + 1, attempts, e))
                worker.kill()
                time.sleep(2.0)
        worker.log("giving up after %d launch attempts" % attempts)
        return False

    def worker_thread(self, worker):
        wlog = worker.log
        if not self._safe_launch(worker):
            worker.state = STOPPED
            return
        while not self.stop_event.is_set():
            if worker.recycle_event.is_set():
                wlog("recycle requested; quitting atlas cleanly")
                worker.state = RESTARTING
                worker.quit()
                worker.recycle_event.clear()
                if not self._safe_launch(worker):
                    break
            try:
                k = self.queue.get_nowait()
            except Empty:
                wlog("queue empty; worker exiting")
                break
            try:
                worker.state = IDLE
                if worker.is_finished(k):
                    wlog("pair %d already finished; skipping" % k)
                    continue
                t0 = time.time()
                worker.state = COMPUTING
                worker.compute_started = t0
                worker.low_cpu_since = 0.0
                worker.compute(k)
                worker.state = IDLE
                dt = time.time() - t0
                grew = worker.write_pair(k, dt)
                if grew <= 0:
                    # atlas errored on compute or write (sentinel still returns);
                    # don't count it -- requeue and recycle to reset state.
                    wlog("pair %d produced no output (atlas error?); requeue+recycle" % k)
                    if not self._record_fail(k):
                        self.queue.put(k)
                    else:
                        wlog("pair %d quarantined (poison)" % k)
                    worker.state = RESTARTING
                    worker.quit()
                    if not self._safe_launch(worker):
                        break
                    continue
                with self.done_lock:
                    self.done_count += 1
                    done = self.done_count
                wlog("pair %d done in %s (total done=%d, qsize=%d)" % (
                    k, fmt_dur(dt), done, self.queue.qsize()))
                # trigger an init merge periodically
                if done % 25 == 0:
                    self.merger.request.set()
            except (HangError, WorkerDied) as e:
                worker.state = RESTARTING
                wlog("worker fault on pair %d: %r -> restarting" % (k, e))
                worker.kill()
                if self._record_fail(k):
                    wlog("pair %d quarantined (poison)" % k)
                else:
                    self.queue.put(k)        # retry later
                # always merge after a fault: salvage shared progress
                self.merger.request.set()
                worker.recycle_event.clear()
                if not self._safe_launch(worker):
                    break
            except Exception as e:           # noqa: BLE001
                wlog("UNEXPECTED error on pair %d: %r" % (k, e))
                self.queue.put(k)
                worker.state = RESTARTING
                worker.kill()
                if not self._safe_launch(worker):
                    break
        worker.state = STOPPED
        worker.quit()
        wlog("stopped")

    # -- governor ---------------------------------------------------------- #
    def governor_thread(self):
        cfg = self.cfg
        # prime cpu_percent counters
        procs = {}
        while not self.stop_event.is_set():
            time.sleep(cfg.poll_interval)
            sizes = []          # (rss_gb, worker)
            total = 0.0
            for w in self.workers:
                pid = w.pid
                if pid is None:
                    continue
                try:
                    p = procs.get(pid)
                    if p is None or p.pid != pid:
                        p = psutil.Process(pid)
                        procs[pid] = p
                        p.cpu_percent(None)      # prime
                    rss = p.memory_info().rss / GB
                    cpu = p.cpu_percent(None)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
                total += rss
                sizes.append((rss, w))
                self._check_stall(w, cpu)
            self._enforce_memory(total, sizes)
            self.log("MEM total=%.1f GB across %d workers; qsize=%d done=%d" % (
                total, len(sizes), self.queue.qsize(), self.done_count))

    def _check_stall(self, w, cpu):
        """Flag a COMPUTING worker pinned near 0% CPU for too long as hung.
        The worker thread's compute() is blocked reading stdout; killing the
        proc makes that read hit EOF, so the worker restarts and requeues."""
        if w.state != COMPUTING:
            w.low_cpu_since = 0.0
            return
        if cpu < 2.0:
            if w.low_cpu_since == 0.0:
                w.low_cpu_since = time.time()
            elif time.time() - w.low_cpu_since > self.cfg.stall_seconds:
                self.log("STALL: job %d idle %s while computing -> killing" % (
                    w.job, fmt_dur(time.time() - w.low_cpu_since)))
                w.kill()
                w.low_cpu_since = 0.0
        else:
            w.low_cpu_since = 0.0

    def _enforce_memory(self, total, sizes):
        cfg = self.cfg
        # emergency: above the hard ceiling -> hard-kill the biggest now
        if total > cfg.hard_total_gb:
            sizes.sort(reverse=True)
            self.log("EMERGENCY: total %.1f GB > %.1f GB ceiling; killing biggest" % (
                total, cfg.hard_total_gb))
            for rss, w in sizes[:max(1, len(sizes) // 20)]:
                self.log("emergency kill job %d (%.1f GB)" % (w.job, rss))
                w.kill()
            return
        # per-process soft cap -> graceful recycle
        for rss, w in sizes:
            if rss > cfg.max_proc_gb and not w.recycle_event.is_set():
                self.log("job %d %.2f GB > %.2f GB cap; flag graceful recycle" % (
                    w.job, rss, cfg.max_proc_gb))
                w.recycle_event.set()
        # global soft cap -> recycle the biggest until projected under cap
        if total > cfg.max_total_gb:
            sizes.sort(reverse=True)
            projected = total
            for rss, w in sizes:
                if projected <= cfg.max_total_gb:
                    break
                if not w.recycle_event.is_set():
                    self.log("global %.1f GB > %.1f GB; flag job %d (%.1f GB)" % (
                        total, cfg.max_total_gb, w.job, rss))
                    w.recycle_event.set()
                    projected -= rss          # recycling frees ~its RSS

    # -- merge coordinator ------------------------------------------------- #
    def merge_thread(self):
        cfg = self.cfg
        while not self.stop_event.is_set():
            fired = self.merger.request.wait(timeout=cfg.merge_interval)
            self.merger.request.clear()
            if self.stop_event.is_set():
                break
            need = fired or (time.time() - self.merger.last_merge >= cfg.merge_interval)
            if need:
                self.merger.merge()

    # -- run --------------------------------------------------------------- #
    def run(self):
        cfg = self.cfg
        if self.queue.empty():
            self.log("nothing to do; exiting")
            return
        for job in range(cfg.workers):
            w = AtlasWorker(job, cfg)
            self.workers.append(w)
        gov = threading.Thread(target=self.governor_thread, name="governor", daemon=True)
        mrg = threading.Thread(target=self.merge_thread, name="merger", daemon=True)
        gov.start()
        mrg.start()
        self.log("starting %d workers" % cfg.workers)
        for w in self.workers:
            t = threading.Thread(target=self.worker_thread, args=(w,),
                                 name="job-%d" % w.job)
            t.start()
            self.threads.append(t)
            time.sleep(0.05)    # stagger launches a touch
        # wait for workers to finish (queue drained) or a stop signal
        for t in self.threads:
            while t.is_alive():
                t.join(timeout=1.0)
        self.log("all workers finished; final init merge")
        self.stop_event.set()
        self.merger.merge(blocking=True)
        self.log("run complete: %d pairs computed" % self.done_count)

    def shutdown(self, signum, _frame):
        self.log("signal %d received; draining workers gracefully" % signum)
        self.stop_event.set()
        for w in self.workers:
            w.recycle_event.clear()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_config(argv):
    p = argparse.ArgumentParser(description="Parallel FPP unitary-dual driver.")
    p.add_argument("--preset", choices=sorted(PRESETS), help="f4s (test) or e7s")
    p.add_argument("-G", "--name", help="short name for files/dirs (e.g. f4s)")
    p.add_argument("-n", "--workers", type=int, help="fixed number of workers")
    p.add_argument("--max-total-gb", type=float, help="total RSS cap across workers")
    p.add_argument("--max-proc-gb", type=float, help="per-worker RSS recycle trigger")
    p.add_argument("-i", "--reference-init", help="reference init .at to reset from")
    p.add_argument("--init-file", help="live init file (default <name>_init.at)")
    p.add_argument("-x", "--aux", action="append", default=[],
                   help="auxiliary .at file to load (repeatable)")
    p.add_argument("-o", "--output-root", default=DEFAULT_OUTPUT_ROOT)
    p.add_argument("-r", "--reverse", action="store_true")
    p.add_argument("-q", "--limit", type=int, default=-1,
                   help="cap the queue to this many pairs")
    p.add_argument("--stall-seconds", type=float)
    p.add_argument("--poll-interval", type=float, default=5.0,
                   help="seconds between memory-governor polls")
    p.add_argument("--merge-interval", type=float, default=120.0)
    p.add_argument("--keep-init", action="store_true",
                   help="do NOT reset the init file from the reference")
    p.add_argument("--settings", default="fpp_settings.at",
                   help="settings file loaded last to override defaults "
                        "(debugging flags); default fpp_settings.at")
    p.add_argument("--no-settings", action="store_true",
                   help="do not load any settings/override file")
    a = p.parse_args(argv)

    pre = PRESETS[a.preset] if a.preset else None

    def pick(val, attr, default=None):
        if val is not None:
            return val
        if pre is not None:
            return getattr(pre, attr)
        return default

    name = pick(a.name, "name")
    if not name:
        p.error("need --preset or --name")
    workers = pick(a.workers, "workers")
    reference = pick(a.reference_init, "reference_init")
    if not reference:
        p.error("need --preset or --reference-init")
    reference = reference if os.path.isabs(reference) else os.path.join(SCRIPTS_DIR, reference)
    if not os.path.exists(reference):
        p.error("reference init not found: %s" % reference)
    aux = a.aux if a.aux else (pre.aux_files if pre else [])
    aux = [f if os.path.isabs(f) else os.path.join(SCRIPTS_DIR, f) for f in aux]
    for f in aux:
        if not os.path.exists(f):
            p.error("aux file not found: %s" % f)
    init_file = a.init_file or os.path.join(SCRIPTS_DIR, "%s_init.at" % name)
    if not os.path.isabs(init_file):
        init_file = os.path.join(SCRIPTS_DIR, init_file)

    max_total = pick(a.max_total_gb, "max_total_gb")
    max_proc = pick(a.max_proc_gb, "max_proc_gb")
    if max_total is None or max_proc is None:
        p.error("need --preset, or both --max-total-gb and --max-proc-gb")
    if workers is None:
        p.error("need --preset or --workers (-n)")

    settings = ""
    if not a.no_settings and a.settings:
        settings = a.settings if os.path.isabs(a.settings) else os.path.join(SCRIPTS_DIR, a.settings)
        if not os.path.exists(settings):
            p.error("settings file not found: %s" % settings)

    return Config(
        name=name,
        workers=workers,
        max_total_gb=max_total,
        max_proc_gb=max_proc,
        reference_init=reference,
        init_file=init_file,
        aux_files=aux,
        output_root=a.output_root,
        stall_seconds=pick(a.stall_seconds, "stall_seconds", 900.0),
        reverse=a.reverse,
        limit=a.limit,
        poll_interval=a.poll_interval,
        merge_interval=a.merge_interval,
        keep_init=a.keep_init,
        settings_file=settings,
    )


def main(argv):
    if not os.path.exists(ATLAS_EXE):
        sys.exit("atlas executable not found: %s" % ATLAS_EXE)
    cfg = build_config(argv)
    run = FPPRun(cfg)
    run.setup()
    signal.signal(signal.SIGINT, run.shutdown)
    signal.signal(signal.SIGTERM, run.shutdown)
    run.run()


if __name__ == "__main__":
    main(sys.argv[1:])
