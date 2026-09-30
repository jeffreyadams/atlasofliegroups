#!/usr/bin/python3
"""
fpp_claude.py -- parallel driver for computing the FPP unitary dual with atlas.

Rewrite of fpp.py.  Runs a *fixed* pool of N persistent atlas workers, each
owning one long-lived `atlas` subprocess.  Workers pull (x,lambda) pair indices
from a shared queue, compute the FPP-unitary facets for each pair, and append
the results to their own data file.  A memory governor keeps total RSS under a
cap by asking workers to recycle (exit atlas cleanly between pairs and relaunch
fresh); only an imminent-OOM emergency ever kills a worker mid-computation.

Key design points (see the prompts fpp_prompt.txt / fpp_prompt_2.txt):
  * Fixed worker count (E7_s: 500, F4_s test: 10).
  * Global RSS cap (E7_s: 12 TB) enforced by *graceful* recycling.
  * Robust data flow: a unique sentinel terminates every command (no more
    guessing with "END"), big diagnostic output is redirected to a file by
    atlas so the pipe only carries small status lines, and a poll()-based
    reader detects hangs.
  * The shared init file (big_unitary_hash) is reset from a reference at start
    and periodically rebuilt by merging every worker's data file back in.
  * Jobs show up in `top` as atlas_<job_id> via per-job symlinks to the atlas
    executable (argv[0]/comm trick).
  * Every atlas invocation loads all.at, then the preload files (report.at,
    FPP.at -- the driver-facing API), then the init, aux files and finally the
    settings file, so a worker, a one-shot query and a merge all see the same
    definitions.
  * Provenance: each run records what it actually ran with -- logs/manifest.txt
    (versions, hashes, git HEAD, caps), logs/flags.txt (the value of every
    setting-like atlas variable, with the settings file applied) and
    logs/flags_changed.txt (the diff against the previous run).  These, plus
    snapshots of the inputs, are committed to a git repo at the output root;
    the bulk data is deliberately left untracked.

The driver does not have to live inside atlas-scripts: it finds all.at next to
itself, or under <driver dir>/atlasofliegroups/atlas-scripts, or wherever
--scripts-dir / $FPP_SCRIPTS_DIR points.

Typical use:
    ./fpp_claude.py --preset f4s -n 3     # 3 workers, the small test run
    ./fpp_claude.py --preset e7s          # 500 workers, 12 TB, the real thing
    ./fpp_claude.py --preset e7s -n 200 --max-total-gb 8000   # overrides

Stop with Ctrl-C (SIGINT): workers drain gracefully and a final init merge runs.
"""

import argparse
import datetime
import glob
import hashlib
import json
import os
import platform
import re
import resource
import shlex
import socket
import select
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field, replace as dc_replace
from collections import deque
from queue import Queue, Empty

try:
    import psutil
except ImportError:
    sys.exit("fpp_claude.py requires psutil (pip install --user psutil)")

# --------------------------------------------------------------------------- #
# Fixed locations
# --------------------------------------------------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))


def find_scripts_dir(explicit=None):
    """Locate the atlas-scripts directory (the one holding all.at).

    The driver no longer has to live inside atlas-scripts: it looks next to
    itself, then in the usual checkout layout underneath itself.  $FPP_SCRIPTS_DIR
    and --scripts-dir override."""
    cands = []
    for c in (explicit, os.environ.get("FPP_SCRIPTS_DIR"), HERE,
              os.path.join(HERE, "atlasofliegroups", "atlas-scripts"),
              os.path.join(HERE, "atlas-scripts")):
        if c and c not in cands:
            cands.append(c)
    for d in cands:
        if os.path.exists(os.path.join(d, "all.at")):
            return os.path.realpath(d)
    sys.exit("cannot find atlas-scripts (no all.at in: %s)" % ", ".join(cands))


SCRIPTS_DIR = find_scripts_dir()
ATLAS_DIR = os.path.dirname(SCRIPTS_DIR)
ATLAS_EXE = os.path.join(ATLAS_DIR, "atlas")
DEFAULT_OUTPUT_ROOT = "/scr/jdada11/fpp_claude"
GROUP_SYMBOL = "G_temp"                             # name the init/data files use

# Files loaded right after all.at, before the init file, by *every* atlas
# invocation (workers, one-shots, merges).  These hold the driver-facing API
# (write_one_pair, write_real_form_plus, report_datum, ...) which lives outside
# the branch's all.at manifest.
DEFAULT_PRELOAD = ["report.at", "FPP.at"]

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
    stop_at_remaining: int = 0   # stop with this many pairs outstanding (0=off)


PRESETS = {
    "f4s": Preset(
        name="f4s", workers=10,
        max_total_gb=60.0, max_proc_gb=6.0,
        reference_init="f4sinitreference.at", aux_files=[],
        stall_seconds=300.0,
    ),
    "e7s": Preset(
        name="e7s", workers=500,
        max_total_gb=12000.0, max_proc_gb=40.0,
        reference_init="e7sinitreference.at",
        aux_files=["edges_F4_E6_E7.at", "coh_ind_E7_centered2PlusE6.at"],
        stall_seconds=1800.0,
    ),
    # E8 presets: the reference init files do not exist yet.  Sizes are first
    # guesses to be re-tuned from the E7_s numbers once a real run is possible.
    # Measured 2026-08-25: 5,926,140 pairs (none finished), KGB size 67,110,
    # mean 88 lambda per x.  Sampled cost 2.1 s/pair over 1600 random pairs, so
    # roughly 3500 core-hours.  Reference init loads in 13 s; the edges file
    # adds 154 s and 5.3 GB, so a worker sits at ~5.5 GB before it computes
    # anything -- 784 of them is about 4.3 TB.  Unlike E7_s, cost is
    # uncorrelated with pair index (Spearman 0.09), so --reverse gains nothing.
    "e8q": Preset(
        name="e8q", workers=784,
        max_total_gb=15000.0, max_proc_gb=60.0,
        reference_init="e8qinitreference.at",
        aux_files=["/u02/jdada11/atlasSoftware/to_ht_branch_jeff_2/"
                   "atlas-scripts/edges_F_E.at"],
        stall_seconds=3600.0,
        # The last 200 pairs cost 8h19m of e8q_1's 19h27m -- 43% of the wall
        # clock for 0.003% of the work.  Finish them later with --keep-init.
        stop_at_remaining=200,
    ),
    "e8s": Preset(
        name="e8s", workers=800,
        max_total_gb=12000.0, max_proc_gb=80.0,
        reference_init="e8sinitreference.at",
        aux_files=[],
        stall_seconds=3600.0,
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
    preload_files: list = field(default_factory=list)   # loaded before the init
    git_repo: bool = True               # keep run provenance in git at output_root
    mem_log_interval: float = 60.0      # seconds between governor MEM log lines
    shards: int = 1                     # driver processes to split the workers over
    stagger: float = 0.05               # seconds between worker launches (global)
    persistent_merge: bool = True       # fold worker output as the run proceeds
    fused: bool = True                  # one do_one_pair round trip per pair
    chunk: int = 0                      # consecutive pairs handed out at once (0=auto)
    giveback_after: float = 60.0        # hand a slow run back after this long
    share: bool = False                 # broadcast new parameters to live workers
    share_cap_mb: float = 4.0           # most a worker applies in one go
    share_concurrency: int = 32         # workers absorbing at once, machine-wide
    stop_at_remaining: int = 0          # stop once this few pairs are outstanding (0=off)

    # filled in during setup
    run_dir: str = ""
    logs_dir: str = ""
    snapshot_dir: str = ""
    symlinks_dir: str = ""
    round: int = 0
    definition_file: str = ""
    bcast_file: str = ""                # append-only stream of new parameters

    @property
    def merge_base(self):
        """What a merge starts from before folding in this run's data files.

        A fresh run rebuilds the init as reference + everything the workers
        wrote, which is complete because the run computed every pair.  A
        --keep-init run did not: its run directory holds only the pairs that
        were still outstanding, so starting from the reference would write an
        init containing the reference plus a handful of pairs and silently
        discard everything earlier runs accumulated.  On E8_q that would have
        replaced 211,017 parameters with about 49,890."""
        return self.init_file if self.keep_init else self.reference_init

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


def fmt_gb(x):
    return "%.2f GB" % x if x < 10 else "%.1f GB" % x


def pair_times(logs_dir):
    """Every per-pair elapsed time (seconds) recorded in logs/<job>.pairs."""
    out = []
    for fn in sorted(glob.glob(os.path.join(logs_dir, "*.pairs"))):
        try:
            with open(fn, errors="replace") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split(",")
                    if len(parts) >= 2:
                        try:
                            out.append(float(parts[1]))
                        except ValueError:
                            pass
        except OSError:
            continue
    return out


def percentile(sorted_vals, q):
    if not sorted_vals:
        return 0.0
    i = int(round(q * (len(sorted_vals) - 1)))
    return sorted_vals[max(0, min(len(sorted_vals) - 1, i))]


def build_summary(run_dir):
    """Render <run_dir>/logs/summary.txt from what the run left on disk.

    Reads logs/stats.json (published every governor tick, so this works on an
    interrupted run too) and the per-worker logs/<job>.pairs timing records."""
    run_dir = os.path.abspath(run_dir)
    logs = os.path.join(run_dir, "logs")
    if not os.path.isdir(logs):
        return "no logs directory in %s\n" % run_dir
    st = {}
    sp = os.path.join(logs, "stats.json")
    if os.path.exists(sp):
        try:
            with open(sp) as f:
                st = json.load(f)
        except (OSError, ValueError):
            st = {}
    times = sorted(pair_times(logs))
    det = st.get("workers_detail", [])
    active = [w for w in det if w.get("first_launch")]

    started = st.get("started", 0.0)
    updated = st.get("updated", 0.0)
    finished = st.get("finished") or 0.0
    end = finished or updated
    wall = max(0.0, end - started) if started else 0.0
    comp0 = st.get("compute_started") or started
    comp_wall = max(0.0, end - comp0) if comp0 else 0.0

    worker_cpu = sum(w.get("cpu_seconds", 0.0) for w in det)
    worker_wall = sum(max(0.0, (w.get("last_stop") or end) - w["first_launch"])
                      for w in active)
    compute_secs = sum(w.get("compute_seconds", 0.0) for w in det)
    child = st.get("child_cpu", {})
    child_cpu = child.get("user", 0.0) + child.get("system", 0.0)
    n = len(active) or st.get("workers", 0)
    done = st.get("done", len(times))

    def line(label, value, note=""):
        return "  %-32s %s%s" % (label, value, ("   " + note) if note else "")

    L = []
    A = L.append
    A("FPP run summary: %s" % st.get("run", os.path.basename(run_dir)))
    A("=" * 72)
    A(line("run directory", run_dir))
    if started:
        A(line("started", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(started))))
        A(line("finished" if finished else "last update",
               time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end)),
               "" if finished else "(run not finished)"))
    A("")
    A("WORK")
    A(line("(x,lambda) pairs queued", st.get("queue_total", "?")))
    A(line("pairs computed", done))
    if st.get("skipped_already_finished"):
        A(line("pairs skipped (already finished)", st["skipped_already_finished"]))
    if st.get("queue_remaining"):
        A(line("pairs still queued", st["queue_remaining"]))
    if st.get("requeued"):
        A(line("pairs requeued after a fault", st["requeued"]))
    if st.get("quarantined"):
        A(line("pairs quarantined (poison)", st["quarantined"]))
    up = st.get("unitary_params")
    if up is not None:
        A(line("unitary parameters in shared hash", up,
               "(long_match lines in the rebuilt init)"))
    upw = st.get("unitary_params_written")
    if upw is not None:
        note = "(written by the workers)" if upw == up else \
               "(written by the workers -- DIFFERS from the hash)"
        A(line("unitary parameters written", upw, note))
    A("")
    A("TIME")
    A(line("wall time, whole run", fmt_dur(wall)))
    A(line("wall time, computing", fmt_dur(comp_wall), "(after setup)"))
    A(line("worker wall time, summed", fmt_dur(worker_wall),
           "= %d workers alive" % n if n else ""))
    A(line("CPU time, workers", fmt_dur(worker_cpu)))
    A(line("CPU time, all atlas children", fmt_dur(child_cpu),
           "(workers + merges + one-shot queries)"))
    if n:
        A(line("CPU time per worker, mean", fmt_dur(worker_cpu / n)))
        A(line("wall time per worker, mean", fmt_dur(worker_wall / n)))
    if comp_wall > 0 and n:
        A(line("parallel efficiency", "%.0f%%" % (100.0 * worker_cpu /
                                                  (comp_wall * n)),
               "worker CPU / (computing wall x workers)"))
        A(line("mean concurrency", "%.2f" % (worker_cpu / comp_wall),
               "cpus busy on average"))
    A("")
    A("PER (x,lambda) PAIR" + ("" if times else "  (no per-pair records found)"))
    if times:
        tot = sum(times)
        A(line("pairs measured", len(times)))
        A(line("summed compute time", fmt_dur(tot)))
        A(line("mean", "%.3f s" % (tot / len(times))))
        A(line("median", "%.3f s" % percentile(times, 0.50)))
        A(line("90th / 99th percentile", "%.3f s / %.3f s" % (
            percentile(times, 0.90), percentile(times, 0.99))))
        A(line("max", "%.3f s" % times[-1]))
        if comp_wall > 0:
            A(line("throughput", "%.1f pairs/s" % (len(times) / comp_wall),
                   "= %.0f pairs per worker-hour" % (
                       3600.0 * len(times) / worker_wall) if worker_wall else ""))
        if compute_secs > 0 and worker_cpu > 0:
            A(line("compute share of worker CPU", "%.0f%%" % (
                100.0 * compute_secs / worker_cpu),
                "rest is atlas startup, is_finished, write_one_pair"))
    A("")
    A("MEMORY")
    caps = st.get("caps", {})
    A(line("peak total RSS", fmt_gb(st.get("peak_total_rss_gb", 0.0)),
           "(cap %.1f GB)" % caps.get("max_total_gb", 0.0)))
    if det:
        A(line("peak single-worker RSS",
               fmt_gb(max(w.get("peak_rss_gb", 0.0) for w in det)),
               "(recycle trigger %.2f GB)" % caps.get("max_proc_gb", 0.0)))
    if child.get("maxrss_gb"):
        A(line("largest atlas process seen", fmt_gb(child["maxrss_gb"]),
               "(getrusage ru_maxrss)"))
    mg = st.get("merges", {})
    if mg.get("count"):
        A("")
        persistent = mg.get("persistent")
        A("INIT MERGES (persistent)" if persistent else "INIT MERGES")
        A(line("ingest passes" if persistent else "merges run", mg["count"]))
        # A persistent merger folds while the workers run, so its time overlaps
        # the compute phase and must not be charged against wall clock; a batch
        # merge blocks, so there the share of wall is the point.
        A(line("time folding" if persistent else "time in merges",
               fmt_dur(mg.get("seconds", 0.0)),
               "overlapped with the run" if persistent else
               ("%.0f%% of wall" % (100.0 * mg.get("seconds", 0.0) / wall)
                if wall else "")))
        A(line("mean pass" if persistent else "mean merge",
               "%.1f s" % (mg.get("seconds", 0.0) / mg["count"])))
        A(line("data folded" if persistent else "last merge input",
               "%.1f MB" % (mg.get("last_input_bytes", 0) / 1e6)))
        A(line("rebuilt init size", "%.1f MB" % (mg.get("last_init_bytes", 0) / 1e6)))
    if det:
        A("")
        A("PER WORKER")
        A("  %5s %8s %8s %10s %10s %6s %9s %7s %7s" % (
            "job", "pairs", "skipped", "wall", "cpu", "cpu%", "peak RSS",
            "recycle", "faults"))
        for w in sorted(det, key=lambda d: -d.get("pairs", 0))[:40]:
            ww = max(0.0, (w.get("last_stop") or end) - (w.get("first_launch") or end))
            cs = w.get("cpu_seconds", 0.0)
            A("  %5d %8d %8d %10s %10s %5.0f%% %8.2fG %7d %7d" % (
                w.get("job", -1), w.get("pairs", 0), w.get("skipped", 0),
                fmt_dur(ww), fmt_dur(cs), (100.0 * cs / ww) if ww > 0 else 0.0,
                w.get("peak_rss_gb", 0.0), w.get("recycles", 0),
                w.get("faults", 0)))
        if len(det) > 40:
            A("  ... %d more workers (busiest 40 shown)" % (len(det) - 40))
    A("")
    return "\n".join(L) + "\n"


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
        # One CSV line per completed pair, at full precision: atlas only ever
        # sees whole seconds (reportDatum.x_lambda_total_time is an int), which
        # rounds every small group to zero.  Columns:
        #   pair, elapsed_seconds, start_epoch, bytes_written
        self.pairs_path = os.path.join(cfg.logs_dir, "%d.pairs" % job)
        new_pairs = not os.path.exists(self.pairs_path)
        self.pairsfh = open(self.pairs_path, "a", buffering=1)
        if new_pairs:
            self.pairsfh.write("# pair,elapsed_seconds,start_epoch,bytes\n")
        # parameter sharing: how far this worker has read the broadcast stream
        self.bcast_off = 0
        self.bcast_bytes = 0        # bytes of shared parameters applied
        self.bcast_seconds = 0.0    # time spent applying them
        self.bcast_applies = 0
        # cross-thread signals from the governor
        self.recycle_event = threading.Event()    # graceful recycle requested
        # governor bookkeeping (only the governor writes these)
        self.state = STOPPED
        self.compute_started = 0.0
        self.low_cpu_since = 0.0
        self.peak_rss = 0.0
        self.cpu_seconds = 0.0        # CPU of atlas processes already replaced
        self._cpu_pid = None          # pid the running total below belongs to
        self._cpu_cur = 0.0           # latest CPU total for that pid
        # run accounting (worker thread writes these)
        self.first_launch = 0.0
        self.last_stop = 0.0
        self.pairs_done = 0
        self.pairs_skipped = 0
        self.compute_seconds = 0.0
        self.recycles = 0
        self.faults = 0
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
        # A relaunched worker reads init_file again, which during a run is still
        # the reference: everything shared so far is gone, so replay it.
        self.bcast_off = 0
        args = ([self.symlink, "all.at"] + self.cfg.preload_files
                + [self.cfg.init_file] + self.cfg.aux_files)
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
        if not self.first_launch:
            self.first_launch = time.time()
        # Drain the (large) startup banner from loading all.at + init + aux, so
        # later command output isn't mixed with it.
        self._command("", self.cfg.startup_timeout)
        if self.cfg.fused:
            # do_one_pair's result is parked here so the compute's (verbose)
            # output can be redirected to the job log while only the data
            # string comes back over the pipe -- both in one round trip.
            self._command('set res_cur = ""', self.cfg.cmd_timeout)
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

    @property
    def cpu_total(self):
        return self.cpu_seconds + self._cpu_cur

    def note_cpu(self, pid, total):
        """Record the CPU total observed for atlas process `pid`.  A worker gets
        a fresh process on every recycle, so the previous process's last known
        total is rolled into cpu_seconds when the pid changes."""
        if self._cpu_pid != pid:
            self.cpu_seconds += self._cpu_cur
            self._cpu_pid, self._cpu_cur = pid, 0.0
        self._cpu_cur = total

    def sample_cpu(self):
        """Read this atlas process's CPU total.  Called by the governor each
        poll and once more just before the process is reaped, so the accounting
        loses at most the last few milliseconds of a worker's life."""
        pid = self.pid
        if pid is None:
            return
        try:
            t = psutil.Process(pid).cpu_times()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return
        self.note_cpu(pid, t.user + t.system)

    def _reap(self):
        if not self.proc:
            return
        self.sample_cpu()
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

    def record_pair(self, k, elapsed_s, started_at, written):
        self.pairs_done += 1
        self.compute_seconds += elapsed_s
        self.pairsfh.write("%d,%.4f,%.3f,%d\n" % (k, elapsed_s, started_at, written))

    def apply_broadcast(self):
        """Fold parameters other workers have found into this worker's hash.

        Each worker starts from the reference init and, without this, only ever
        learns what it discovers itself: on E8_q every worker ended holding
        about 50,280 FPP-unitary faces while 211,017 were known across the
        machine, and 93% of its hash lookups missed.  The broadcast file is the
        parameter-bearing blocks the merger has already isolated, so replaying
        it here costs only the parse -- no reload of the edges file, which is
        what made the earlier restart-based version of this idea uneconomic.

        Reads at most share_cap_mb per call and always stops on a line boundary;
        whatever is left is picked up after the next pair.  Returns bytes
        applied."""
        path = self.cfg.bcast_file
        if not path:
            return 0
        try:
            size = os.path.getsize(path)
        except OSError:
            return 0
        if size <= self.bcast_off:
            return 0
        cap = int(self.cfg.share_cap_mb * 1e6)
        with open(path, "rb") as f:
            f.seek(self.bcast_off)
            buf = f.read(min(cap, size - self.bcast_off))
        end = buf.rfind(b"\n")
        if end < 0:
            return 0                       # no complete line yet
        payload = buf[:end + 1]
        t0 = time.time()
        # The payload is megabytes and atlas answers it with megabytes of its
        # own (every `set` in a block echoes).  Writing it all before reading
        # deadlocks: the pipes are ~64 KB each way, so atlas blocks in
        # pipe_write on stdout while we block writing stdin.  Push from a
        # second thread and drain stdout here, concurrently.
        self.seq += 1
        token = "@@SENTINEL_%d_%d@@" % (self.job, self.seq)
        data = payload + ('prints("%s")\n' % token).encode("utf-8")
        werr = []

        def _push():
            try:
                self.proc.stdin.write(data)
                self.proc.stdin.flush()
            except Exception as e:                   # noqa: BLE001
                werr.append(e)

        th = threading.Thread(target=_push, name="bcast-%d" % self.job,
                              daemon=True)
        th.start()
        deadline = time.time() + self.cfg.cmd_timeout + 60.0 * len(payload) / 1e6
        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                raise HangError("broadcast apply stalled (%d bytes)" % len(payload))
            line = self.reader.readline(remaining)
            if line is None:
                raise HangError("broadcast apply timed out")
            if line == b"":
                raise WorkerDied()
            if token in line.decode("latin-1", "replace"):
                break
        th.join(timeout=5.0)
        if werr:
            raise WorkerDied()
        self.bcast_off += len(payload)
        self.bcast_bytes += len(payload)
        self.bcast_seconds += time.time() - t0
        self.bcast_applies += 1
        return len(payload)

    def do_pair(self, k):
        """is_finished + compute + write in a single round trip.

        The old path sent three commands per pair and each re-derived
        (x,lambda) from the flat index -- four times in all, counting the one
        inside report_datum.  do_one_pair derives it once.  The redirect puts
        the computation's verbose output in the job log, exactly as before,
        while `prints(res_cur)` returns just the data.

        Returns (status, compute_seconds, bytes_written); status is one of
        "done", "skip" (already finished) or "fail" (atlas errored)."""
        self.log("---- pair %d: FPP_unitary_hash_bottom_layer output ----" % k)
        cmd = ('>>"%s" res_cur := do_one_pair(%s,%d,%d)\nprints(res_cur)'
               % (self.log_path, GROUP_SYMBOL, k, self.job))
        out = self._command(cmd, self.cfg.compute_timeout)
        ms, data, started = None, [], False
        for s in out:
            t = s.strip()
            if not started:
                if t == "SKIP":
                    return ("skip", 0.0, 0)
                if t.startswith("TIME="):
                    started = True
                    try:
                        ms = int(t[5:])
                    except ValueError:
                        ms = 0
                continue
            if t == "END":
                break
            data.append(s)
        if ms is None:
            return ("fail", 0.0, 0)
        written = 0
        with open(self.data_file, "a") as f:
            for line in data:
                line = line.replace("\\n", "\n") + "\n"
                f.write(line)
                written += len(line)
        # atlas timed the computation itself, so this excludes the round trip
        return ("done", ms / 1000.0, written)

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
        self.seconds = 0.0
        self.last_duration = 0.0
        self.last_input_bytes = 0
        self.last_init_bytes = 0

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
            merge_input = os.path.join(self.cfg.logs_dir, "merge_input.at")
            data_files = sorted(glob.glob(os.path.join(self.cfg.run_dir, "[0-9]*.at")))
            with open(merge_input, "wb") as out:
                if data_files:
                    # Feed back everything the workers wrote except the per-pair
                    # report records, which need rA/add (defined only for
                    # analysis).  Do NOT try to select the interesting lines:
                    # on this branch a pair's unitary parameters arrive as a
                    # multi-line block (set tempParams / tempParams[k]:=... /
                    # uhash(G).fill(tempParams)) that only works read in order.
                    subprocess.run(["grep", "-hv", "^void:add("] + data_files,
                                   stdout=out, stderr=subprocess.DEVNULL,
                                   check=False)
            tmp = self.cfg.init_file + ".tmp.%d" % os.getpid()
            args = ([ATLAS_EXE, "all.at"] + self.cfg.preload_files
                    + [self.cfg.merge_base, merge_input])
            # jeff_sizes_flag makes write() emit set_xl_sizes with the real
            # sizes.  It defaults to false on this branch, and without it the
            # rebuilt init has no xl_sizes at all: every xl_pair / is_finished /
            # xl_pairs_todo lookup against it fails.
            script = ('jeff_sizes_flag:=true\n'
                      '>"%s" big_unitary_hash.write()\nquit\n' % tmp)
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
            self.last_duration = time.time() - t0
            self.seconds += self.last_duration
            self.last_input_bytes = os.path.getsize(merge_input)
            self.last_init_bytes = os.path.getsize(self.cfg.init_file)
            self.log("merge #%d done in %s (%d data files, %.1f MB in, "
                     "init=%.1f MB)" % (
                         self.count, fmt_dur(self.last_duration), len(data_files),
                         os.path.getsize(merge_input) / 1e6,
                         os.path.getsize(self.cfg.init_file) / 1e6))
            return True
        except Exception as e:                       # noqa: BLE001  (log everything)
            self.log("MERGE ERROR: %r" % e)
            return False
        finally:
            self.lock.release()


# --------------------------------------------------------------------------- #
# Persistent merge (fold worker output into one live hash as the run proceeds)
# --------------------------------------------------------------------------- #
class PersistentMerger:
    """One long-lived atlas process that accumulates every worker's findings.

    The batch Merger re-parses the whole corpus on every merge.  Measured on the
    e7s_3 corpus: parsing costs ~2.6 s/MB and the write costs ~3 s, so the final
    merge of 397 MB took 23m36s -- 30% of that run's wall clock, all of it
    blocking, with 783 of 784 CPUs idle.  Nothing about that work needs to
    happen at the end: this process stays up for the whole run, is fed only the
    bytes each data file has grown by since last time, and when the workers stop
    it only has to write.  Every byte is parsed exactly once, concurrently with
    the computation, and the blocking tail becomes the ~3 s write.

    Worker data files are appended one block per pair -- several lines
    (set Nullp / set tempParams / tempParams[k]:=parameter(...) /
    uhash(G).fill(tempParams) / finish_num(...)) terminated by that pair's
    `void:add(...)` report record.  Only whole blocks are fed, so a file caught
    mid-append is simply consumed up to its last complete block and the rest
    waits for the next pass.  The void:add records are dropped: reading them
    needs rA/add, which only read_reports.at defines.

    Memory is modest -- ingesting the raw form peaked at 0.73 GB, against the
    13 GB it costs to load the written (long_match) form back in."""

    CHUNK = 8 << 20            # bytes feed per sentinel round trip

    def __init__(self, cfg, log):
        self.cfg = cfg
        self.log = log
        self.lock = threading.Lock()
        self.request = threading.Event()
        self.proc = None
        self.reader = None
        self.seq = 0
        self.offsets = {}          # data file -> bytes consumed so far
        self.count = 0             # ingest passes that fed something
        self.seconds = 0.0
        self.bytes_fed = 0
        self.bytes_parsed = 0
        self.bytes_shared = 0      # parameter blocks written to the broadcast
        # Completion state, accumulated across every worker.  A shared
        # parameter is only actionable if the receiving worker also knows the
        # pair is COMPLETE: finalizing (x,lambda,nu) yields (x',lambda',nu'),
        # and "not in the hash" only means "not unitary" when (x',lambda') is
        # finished.  finish_num ORs bits, so re-sending an x's accumulated mask
        # is idempotent -- one line per x beats one per completed pair.
        self.finish_masks = {}     # x -> OR of every mask seen for it
        self.finish_dirty = set()  # x whose mask changed since the last flush
        self.finish_lines = 0      # completion lines written to the broadcast
        self.bcast_fh = (open(cfg.bcast_file, "ab", buffering=0)
                         if cfg.share and cfg.bcast_file else None)
        self.last_duration = 0.0
        self.last_merge = 0.0
        self.last_input_bytes = 0
        self.last_init_bytes = 0
        self.failed = False
        self.wrote_init = False
        self.stderr_path = os.path.join(cfg.logs_dir, "merger.stderr")

    # -- process ----------------------------------------------------------- #
    def start(self):
        """Launch the merger on cfg.merge_base: the reference for a fresh run,
        the live init for a --keep-init continuation (see Config.merge_base)."""
        args = ([ATLAS_EXE, "all.at"] + self.cfg.preload_files
                + [self.cfg.merge_base])
        if self.cfg.settings_file:
            args.append(self.cfg.settings_file)
        try:
            self.proc = subprocess.Popen(
                args, cwd=SCRIPTS_DIR, stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=open(self.stderr_path, "ab"),
                bufsize=0)
            self.reader = LineReader(self.proc.stdout.fileno())
            self._round_trip(b"", self.cfg.startup_timeout)
            self.log("persistent merger up (pid %d) on %s%s" % (
                self.proc.pid, os.path.basename(self.cfg.merge_base),
                " (continuing)" if self.cfg.keep_init else ""))
            return True
        except Exception as e:                       # noqa: BLE001
            self.log("persistent merger failed to start: %r" % e)
            self.failed = True
            return False

    def _round_trip(self, payload, timeout):
        """Feed bytes, then wait for a sentinel so we know atlas consumed them."""
        if self.proc is None or self.proc.poll() is not None:
            raise WorkerDied()
        self.seq += 1
        token = "@@MERGE_%d@@" % self.seq
        data = (payload or b"") + ('prints("%s")\n' % token).encode("utf-8")
        werr = []

        def _push():
            try:
                self.proc.stdin.write(data)
                self.proc.stdin.flush()
            except Exception as e:                   # noqa: BLE001
                werr.append(e)

        # Same pipe deadlock as AtlasWorker.apply_broadcast: feed from a second
        # thread so this one can drain atlas's stdout while it is being fed.
        th = threading.Thread(target=_push, name="merge-feed", daemon=True)
        th.start()
        deadline = time.time() + timeout
        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                raise HangError("merger: no sentinel after %.0f s" % timeout)
            line = self.reader.readline(remaining)
            if line is None:
                raise HangError("merger: timed out")
            if line == b"":
                raise WorkerDied()
            if token in line.decode("latin-1", "replace"):
                th.join(timeout=5.0)
                if werr:
                    raise WorkerDied()
                return

    def _data_files(self):
        return sorted(glob.glob(os.path.join(self.cfg.run_dir, "[0-9]*.at")))

    def _ingest_file(self, path):
        """Feed the whole blocks this file has grown by.  Returns bytes consumed."""
        off = self.offsets.get(path, 0)
        try:
            size = os.path.getsize(path)
        except OSError:
            return 0
        if size <= off:
            return 0
        with open(path, "rb") as f:
            f.seek(off)
            buf = f.read(size - off)
        # a block ends at its void:add record; consume up to the last complete one
        cut = buf.rfind(b"\nvoid:add((")
        if cut < 0:
            if not buf.startswith(b"void:add(("):
                return 0                      # no complete block yet
            end = buf.find(b"\n")
        else:
            end = buf.find(b"\n", cut + 1)
        if end < 0:
            return 0
        chunk = buf[:end + 1]
        # Only 1.59% of pairs yield a parameter; the other 98.4% write a block
        # whose only content is its finish mark, wrapped in three lines of
        # boilerplate (set Nullp / set tempParams for i:0 / fill(tempParams)).
        # Feeding that boilerplate costs 29% of the corpus to no effect -- an
        # empty block needs only its finish_num.  Blocks that do carry
        # parameters are fed whole and in order, which is what they require.
        kept, block, shared = [], [], []
        for ln in chunk.split(b"\n"):
            if not ln:
                continue
            if ln.startswith(b"void:add(("):
                if any(x.startswith(b"void:tempParams[") for x in block):
                    kept.extend(block)
                    # A parameter-bearing block is self-contained (it sets Nullp
                    # and tempParams before filling), so it can be replayed into
                    # any live worker.  Only 1.59% of blocks carry parameters,
                    # which is what makes broadcasting them affordable; the
                    # empty 98.4% would cost far more than they could save.
                    shared.extend(block)
                else:
                    kept.extend(x for x in block
                                if x.startswith(b"big_unitary_hash.finish_num"))
                block = []
            else:
                block.append(ln)
        for ln in kept:
            if ln.startswith(b"big_unitary_hash.finish_num("):
                try:
                    body = ln.split(b"(", 1)[1].rstrip(b")\n")
                    parts = body.split(b",")
                    xnum, mask = int(parts[1]), int(parts[2])
                except (IndexError, ValueError):
                    continue
                prev = self.finish_masks.get(xnum, 0)
                new = prev | mask
                if new != prev:
                    self.finish_masks[xnum] = new
                    self.finish_dirty.add(xnum)
        if shared and self.bcast_fh is not None:
            payload = b"\n".join(shared) + b"\n"
            self.bcast_fh.write(payload)
            self.bcast_fh.flush()
            self.bytes_shared += len(payload)
        parts, batch, size_so_far = [], [], 0
        for ln in kept:
            batch.append(ln)
            size_so_far += len(ln) + 1
            if size_so_far >= self.CHUNK:      # split on line boundaries only
                parts.append(b"\n".join(batch) + b"\n")
                batch, size_so_far = [], 0
        if batch:
            parts.append(b"\n".join(batch) + b"\n")
        for part in parts:
            self._round_trip(part, 300.0 + 40.0 * len(part) / 1e6)
            self.bytes_parsed += len(part)
        self.offsets[path] = off + len(chunk)
        return len(chunk)

    def _flush_finish_masks(self):
        """Publish the accumulated completion masks for every x that changed.

        One line per x, carrying the OR of every mask any worker has reported
        for it, so a receiving worker learns about pairs completed by all the
        others.  Without this a broadcast parameter cannot be acted on: the
        algorithm may only conclude "not unitary" from a hash miss when the
        pair is known finished."""
        if not self.cfg.share or self.bcast_fh is None or not self.finish_dirty:
            return
        out = []
        for xnum in sorted(self.finish_dirty):
            out.append(b"big_unitary_hash.finish_num(%s,%d,%d)"
                       % (GROUP_SYMBOL.encode(), xnum, self.finish_masks[xnum]))
        self.finish_dirty.clear()
        payload = b"\n".join(out) + b"\n"
        self.bcast_fh.write(payload)
        self.bcast_fh.flush()
        self.bytes_shared += len(payload)
        self.finish_lines += len(out)

    def ingest(self):
        """One pass over the data files.  Skipped if a pass is already running."""
        if self.failed or self.proc is None:
            return False
        if not self.lock.acquire(blocking=False):
            return False
        try:
            t0 = time.time()
            fed = sum(self._ingest_file(p) for p in self._data_files())
            self._flush_finish_masks()
            if fed:
                self.count += 1
                self.bytes_fed += fed
                self.last_duration = time.time() - t0
                self.seconds += self.last_duration
                self.last_merge = time.time()
                self.last_input_bytes = self.bytes_fed
                self.log("merger ingested %.1f MB in %s (%.1f MB folded, "
                         "%.1f MB actually parsed)" % (
                             fed / 1e6, fmt_dur(self.last_duration),
                             self.bytes_fed / 1e6, self.bytes_parsed / 1e6))
            return True
        except (HangError, WorkerDied, OSError, ValueError) as e:
            self.log("PERSISTENT MERGER FAILED: %r -- will fall back to a batch "
                     "merge at the end" % e)
            self.failed = True
            return False
        finally:
            self.lock.release()

    def finish(self):
        """Final ingest, then write the init.  True if the init was written."""
        if self.failed or self.proc is None:
            return False
        with self.lock:
            try:
                t0 = time.time()
                fed = sum(self._ingest_file(p) for p in self._data_files())
                self.bytes_fed += fed
                self.last_input_bytes = self.bytes_fed
                self.log("merger final ingest: %.1f MB (%.1f MB total) in %s" % (
                    fed / 1e6, self.bytes_fed / 1e6, fmt_dur(time.time() - t0)))
                tmp = self.cfg.init_file + ".tmp.%d" % os.getpid()
                t1 = time.time()
                # jeff_sizes_flag makes write() emit set_xl_sizes with the real
                # sizes; without it every xl_pair / is_finished lookup against
                # the rebuilt init fails.
                self._round_trip(
                    ('jeff_sizes_flag:=true\n>"%s" big_unitary_hash.write()\n'
                     % tmp).encode("utf-8"), 3600.0)
                if not os.path.exists(tmp) or os.path.getsize(tmp) == 0:
                    raise ValueError("merger produced no init")
                os.replace(tmp, self.cfg.init_file)
                self.last_init_bytes = os.path.getsize(self.cfg.init_file)
                self.seconds += time.time() - t0
                self.count += 1
                self.wrote_init = True
                self.log("persistent merger wrote init (%.1f MB) in %s "
                         "(write alone %s)" % (
                             self.last_init_bytes / 1e6,
                             fmt_dur(time.time() - t0), fmt_dur(time.time() - t1)))
                return True
            except (HangError, WorkerDied, OSError, ValueError) as e:
                self.log("PERSISTENT MERGER FAILED at write: %r" % e)
                self.failed = True
                return False

    def quit(self):
        try:
            if self.proc and self.proc.poll() is None:
                self.proc.stdin.write(b"quit\n")
                self.proc.stdin.flush()
                self.proc.wait(timeout=120)
        except (OSError, subprocess.TimeoutExpired):
            try:
                self.proc.kill()
            except OSError:
                pass
        self.proc = None

    def detach(self):
        """Called in a forked child: drop the parent's merger pipes."""
        for stream in ("stdin", "stdout"):
            try:
                getattr(self.proc, stream).close()
            except (OSError, AttributeError):
                pass
        self.proc = None
        self.reader = None


# --------------------------------------------------------------------------- #
# The run orchestrator
# --------------------------------------------------------------------------- #
class FPPRun:
    def __init__(self, cfg):
        self.cfg = cfg
        self.queue = Queue()
        # Pairs handed back from a slow run.  They are known to be expensive --
        # that is why they were handed back -- so they must be taken before the
        # ordinary queue.  Queue is FIFO: on e7s_8 the returned pairs went to
        # the back of a queue holding ~980 chunks and waited 14 minutes anyway.
        self.priority = deque()
        self.priority_lock = threading.Lock()
        # Absorbing the broadcast costs 44 s in one process but 1430 s when all
        # 784 do it at once: 256 concurrent absorbers already run 18.6x slower
        # on a 1568-CPU machine, so it is memory bandwidth, not CPU.  Workers
        # synchronise naturally (they all apply just after finishing a pair),
        # which is the worst case.  Cap how many absorb at a time.  Each shard
        # is a separate process and inherits its own copy of this semaphore, so
        # the machine-wide figure is divided across them.
        self._stop_floor_hit = False
        self._share_slots = max(1, cfg.share_concurrency // max(1, cfg.shards))
        self.share_sem = threading.Semaphore(self._share_slots)
        self.workers = []
        self.threads = []
        self.stop_event = threading.Event()
        self.fail_counts = {}
        self.fail_lock = threading.Lock()
        self.done_count = 0
        self.done_lock = threading.Lock()
        self.queue_total = 0
        self.skipped_count = 0
        self.requeued_count = 0
        self.quarantined_count = 0
        self.peak_total_rss = 0.0
        self.started_at = time.time()
        self.compute_started_at = 0.0
        self.finished_at = 0.0
        self.main_log_fh = None
        self.merger = None
        self.merger_persistent = None
        self.shard = None                 # set in a forked shard child
        self.chunk_size = 0
        self.child_pids = []
        self._git_identity = []
        self._last_mem_log = 0.0
        self.unitary_params = None
        self.unitary_params_written = None

    # -- logging ----------------------------------------------------------- #
    def log(self, msg):
        line = "[%s] %s\n" % (now(), msg)
        self.main_log_fh.write(line)
        self.main_log_fh.flush()

    # -- setup ------------------------------------------------------------- #
    def setup(self, argv=()):
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
        cfg.snapshot_dir = os.path.join(cfg.logs_dir, "snapshot")
        cfg.symlinks_dir = os.path.join(cfg.run_dir, "symlinks")
        cfg.bcast_file = os.path.join(cfg.run_dir, "broadcast.at")
        for d in (cfg.run_dir, cfg.logs_dir, cfg.snapshot_dir, cfg.symlinks_dir):
            os.makedirs(d, exist_ok=True)

        self.main_log_fh = open(os.path.join(cfg.logs_dir, "main.log"), "a", buffering=1)
        self.log("=" * 70)
        self.log("fpp_claude starting; run dir %s" % cfg.run_dir)
        self.log("name=%s workers=%d max_total=%.1f GB max_proc=%.2f GB" % (
            cfg.name, cfg.workers, cfg.max_total_gb, cfg.max_proc_gb))
        self.log("scripts dir %s" % SCRIPTS_DIR)
        if cfg.preload_files:
            self.log("preloaded before init: %s" %
                     " ".join(os.path.basename(f) for f in cfg.preload_files))
        if cfg.share:
            # created empty so workers can stat it from the first pair on
            open(cfg.bcast_file, "ab").close()
            self.log("parameter sharing ON: new parameters broadcast via %s "
                     "(cap %.1f MB per apply)" % (
                         os.path.basename(cfg.bcast_file), cfg.share_cap_mb))
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
        self._git_setup()
        self._write_definition()
        self._write_report_reader()
        self._write_manifest(argv)
        self._record_flags()
        self._git_commit("%s_%d start: %d workers, %.0f GB cap" % (
            cfg.name, cfg.round, cfg.workers, cfg.max_total_gb))
        self.merger = Merger(cfg, self.log)
        if cfg.persistent_merge:
            self.merger_persistent = PersistentMerger(cfg, self.log)
        self._build_queue()

    def _load_args(self):
        """Files each worker loads after all.at, in order: init, aux, settings."""
        args = [self.cfg.init_file] + self.cfg.aux_files
        if self.cfg.settings_file:
            args.append(self.cfg.settings_file)
        return args

    # Names matching any of these (case-insensitively) are treated as settings
    # that steer the atlas computation and get recorded.  On top of these we
    # record every identifier bound to a bare literal (true/false/number) at
    # top level in the .at sources, and everything the settings file assigns.
    SETTING_PATTERNS = (
        "flag", "verbose", "deets", "cutoff", "seat_belt",
        r"_mult\b", r"_frac\b", r"_factor\b", r"_factors\b", r"_level\b",
        r"_limit\b", r"_bound\b", r"_max\b", r"_min\b", r"_size\b",
        r"_depth\b", r"_threshold\b", r"_on\b", r"_off\b",
    )

    # A top-level binding of a plain literal: "set foo = true", "bar:=1/5", ...
    _ASSIGN_RE = re.compile(
        r"(?m)^[ \t]*(?:set[ \t]+)?!?([A-Za-z_][A-Za-z0-9_]*)[ \t]*\+?:?=[ \t]*(.*)$")
    _LITERAL_RE = re.compile(r"^(true|false|-?\d+(?:/\d+)?)[ \t]*(?:\{.*)?$")

    def _flag_candidates(self):
        """Scrape candidate setting names from the .at sources + settings file.

        Only source-sized files are scanned; the multi-megabyte data files
        (edges_*.at, coh_ind_*.at, *initreference.at) hold no settings and
        scanning them is pure cost."""
        interesting = re.compile("|".join(self.SETTING_PATTERNS), re.I)
        names = set()
        for fn in sorted(glob.glob(os.path.join(SCRIPTS_DIR, "*.at"))):
            try:
                if os.path.getsize(fn) > 2 * 1024 * 1024:
                    continue
                text = open(fn, errors="replace").read()
            except OSError:
                continue
            for m in self._ASSIGN_RE.finditer(text):
                name, rhs = m.group(1), m.group(2).strip()
                if interesting.search(name) or self._LITERAL_RE.match(rhs):
                    names.add(name)
        # whatever the settings file assigns counts, whatever it is called
        sf = self.cfg.settings_file
        if sf and os.path.exists(sf):
            with open(sf, errors="replace") as f:
                for m in self._ASSIGN_RE.finditer(f.read()):
                    names.add(m.group(1))
        # atlas keywords / obvious non-values that would just error out
        names -= {"set", "let", "in", "then", "if", "fi", "do", "od", "true",
                  "false", "void", "int", "rat", "bool", "string"}
        return sorted(names)

    def _record_flags(self):
        """Record the value of every setting-like atlas variable, with the run's
        settings file applied, to logs/flags.txt.

        Candidates are scraped from the sources and then queried one per atlas
        command: a name that does not exist (or is not printable) simply errors
        out on its own line and is skipped, so what lands in flags.txt is always
        what this run actually saw.  logs/flags_changed.txt records the diff
        against the previous run in the same output root, which is usually the
        thing you actually want to look at."""
        cfg = self.cfg
        names = self._flag_candidates()
        if not names:
            self.log("no flag candidates found; skipping flags.txt")
            return
        script = "".join('prints("@FLAG@|%s|",%s)\n' % (n, n) for n in names) + "quit\n"
        out, _ = self._atlas_oneshot(self._load_args(), script, timeout=1800,
                                     quiet=True)
        flags = {}
        for line in out.splitlines():
            if line.startswith("@FLAG@|"):
                _, name, val = line.split("|", 2)
                val = val.strip()
                if len(val) > 200:
                    val = val[:200] + " ...[truncated]"
                flags[name] = val
        # names the settings file sets, flagged in the output for quick reading
        overridden = set()
        if cfg.settings_file and os.path.exists(cfg.settings_file):
            with open(cfg.settings_file, errors="replace") as f:
                for m in self._ASSIGN_RE.finditer(f.read()):
                    overridden.add(m.group(1))
        path = os.path.join(cfg.logs_dir, "flags.txt")
        with open(path, "w") as f:
            f.write("# atlas settings in effect for run %s_%d\n" % (cfg.name, cfg.round))
            f.write("# settings file: %s\n" % (cfg.settings_file or "(none)"))
            f.write("# '*' marks a value the settings file assigns\n")
            for name in sorted(flags):
                f.write("%s%s = %s\n" %
                        ("*" if name in overridden else " ", name, flags[name]))
        self.log("recorded %d settings (of %d candidates) -> flags.txt" %
                 (len(flags), len(names)))
        self._diff_flags(flags)

    def _diff_flags(self, flags):
        """Write logs/flags_changed.txt: what differs from the previous run."""
        cfg = self.cfg
        prev_path = None
        for r in range(cfg.round - 1, 0, -1):
            cand = os.path.join(cfg.output_root, "%s_%d" % (cfg.name, r),
                                "logs", "flags.txt")
            if os.path.exists(cand):
                prev_path = cand
                break
        if not prev_path:
            return
        prev = {}
        with open(prev_path, errors="replace") as f:
            for line in f:
                if line.startswith("#") or "=" not in line:
                    continue
                name, val = line.split("=", 1)
                prev[name.strip().lstrip("*").strip()] = val.strip()
        lines = []
        for name in sorted(set(prev) | set(flags)):
            a, b = prev.get(name), flags.get(name)
            if a != b:
                lines.append("%s: %s -> %s" % (name, "(absent)" if a is None else a,
                                               "(absent)" if b is None else b))
        path = os.path.join(cfg.logs_dir, "flags_changed.txt")
        with open(path, "w") as f:
            f.write("# %s_%d vs %s\n" % (cfg.name, cfg.round, prev_path))
            f.write("\n".join(lines) + ("\n" if lines else ""))
        self.log("%d settings differ from %s -> flags_changed.txt" %
                 (len(lines), os.path.relpath(prev_path, cfg.output_root)))

    # -- provenance -------------------------------------------------------- #
    @staticmethod
    def _sha256(path):
        h = hashlib.sha256()
        try:
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1 << 20), b""):
                    h.update(chunk)
        except OSError as e:
            return "unreadable (%r)" % e
        return h.hexdigest()

    def _file_line(self, path):
        try:
            st = os.stat(path)
            when = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime))
            return "%-24s %12d  %s  %s" % (os.path.basename(path), st.st_size,
                                           when, self._sha256(path)[:16])
        except OSError as e:
            return "%-24s  MISSING (%r)" % (os.path.basename(path), e)

    def _snapshot(self):
        """Copy the small inputs that define this run into logs/snapshot/, so the
        run can be reproduced (or explained) later from the run directory alone.
        Big unchanging inputs (aux data files) are recorded by hash only."""
        cfg = self.cfg
        for src in self._provenance_files():
            if os.path.exists(src) and os.path.getsize(src) <= 8 * 1024 * 1024:
                try:
                    shutil.copy(src, cfg.snapshot_dir)
                except OSError as e:
                    self.log("snapshot %s failed: %r" % (src, e))

    def _provenance_files(self):
        """Every file whose contents affect this run, in load order."""
        cfg = self.cfg
        files = [os.path.abspath(__file__), os.path.join(SCRIPTS_DIR, "all.at")]
        files += cfg.preload_files
        files += [cfg.reference_init] + list(cfg.aux_files)
        if cfg.settings_file:
            files.append(cfg.settings_file)
        return files

    def _git_info(self):
        def git(*a):
            try:
                r = subprocess.run(["git"] + list(a), cwd=SCRIPTS_DIR,
                                   capture_output=True, text=True)
                return r.stdout.strip()
            except OSError as e:
                return "(git failed: %r)" % e
        return (git("rev-parse", "HEAD"), git("rev-parse", "--abbrev-ref", "HEAD"),
                git("status", "--short"))

    def _write_manifest(self, argv):
        """logs/manifest.txt: everything needed to say what this run *was*."""
        cfg = self.cfg
        head, branch, dirty = self._git_info()
        vm = psutil.virtual_memory()
        lines = [
            "run             %s_%d" % (cfg.name, cfg.round),
            "started         %s" % now(),
            "host            %s" % socket.gethostname(),
            "user            %s" % (os.environ.get("USER") or "?"),
            "command         %s" % " ".join(shlex.quote(a) for a in
                                            [os.path.abspath(__file__)] + list(argv)),
            "cwd             %s" % os.getcwd(),
            "python          %s (psutil %s)" % (platform.python_version(),
                                                getattr(psutil, "__version__", "?")),
            "machine         %d cpus, %.0f GB RAM" % (psutil.cpu_count(), vm.total / GB),
            "",
            "scripts_dir     %s" % SCRIPTS_DIR,
            "atlas           %s" % ATLAS_EXE,
            "                %s" % self._file_line(ATLAS_EXE),
            "git HEAD        %s (branch %s)" % (head, branch),
            "git status      %s" % (dirty.replace("\n", "\n                ")
                                    if dirty else "(clean)"),
            "",
            "workers         %d" % cfg.workers,
            "max_total_gb    %.1f   (hard ceiling %.1f)" % (cfg.max_total_gb,
                                                            cfg.hard_total_gb),
            "max_proc_gb     %.2f" % cfg.max_proc_gb,
            "stall_seconds   %.0f" % cfg.stall_seconds,
            "poll_interval   %.1f" % cfg.poll_interval,
            "merge_interval  %.1f" % cfg.merge_interval,
            "shards          %d" % cfg.shards,
            "stagger         %.3f" % cfg.stagger,
            "persistent_merge %s" % cfg.persistent_merge,
            "fused           %s" % cfg.fused,
            "chunk           %d" % cfg.chunk,
            "giveback_after  %.1f" % cfg.giveback_after,
            "reverse         %s" % cfg.reverse,
            "limit           %s" % (cfg.limit if cfg.limit >= 0 else "(none)"),
            "keep_init       %s" % cfg.keep_init,
            "output_root     %s" % cfg.output_root,
            "run_dir         %s" % cfg.run_dir,
            "",
            "atlas load order (all.at first, settings last):",
        ]
        for f in self._provenance_files()[1:]:
            lines.append("  " + self._file_line(f))
        lines.append("")
        lines.append("driver:")
        lines.append("  " + self._file_line(os.path.abspath(__file__)))
        path = os.path.join(cfg.logs_dir, "manifest.txt")
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        self.log("wrote manifest.txt (git HEAD %s%s)" %
                 (head[:12], ", dirty" if dirty else ""))

    # -- provenance git repo ----------------------------------------------- #
    GITIGNORE = """\
# fpp_claude provenance repo.
#
# Only the small files that describe a run are tracked: the manifest, the
# recorded atlas settings, the main log, and snapshots of the inputs.  The
# bulk output (per-worker .at data, per-job logs, merge inputs, rebuilt init
# files) stays on disk and out of git -- it is far too big for a repository.
*
!.gitignore
!*/
!/*.at
!*/logs/main.log
!*/logs/manifest.txt
!*/logs/flags.txt
!*/logs/flags_changed.txt
!*/logs/poison.txt
!*/logs/snapshot/*
!*/logs/summary.txt
!*/logs/stats.json
!*/read_reports.at
"""

    def _git(self, *args, **kw):
        cmd = ["git", "-C", self.cfg.output_root]
        if self._git_identity:
            cmd += self._git_identity
        cmd += list(args)
        try:
            return subprocess.run(cmd, capture_output=True, text=True, **kw)
        except OSError as e:
            self.log("git %s failed: %r" % (args[0], e))
            return None

    def _git_setup(self):
        cfg = self.cfg
        self._git_identity = []
        if not cfg.git_repo:
            return
        who = subprocess.run(["git", "config", "user.email"], cwd=cfg.output_root,
                             capture_output=True, text=True)
        if who.returncode != 0 or not who.stdout.strip():
            self._git_identity = ["-c", "user.name=fpp_claude",
                                  "-c", "user.email=fpp_claude@localhost"]
        if not os.path.isdir(os.path.join(cfg.output_root, ".git")):
            r = self._git("init", "-q")
            if r is None or r.returncode != 0:
                self.log("git init failed; provenance repo disabled")
                cfg.git_repo = False
                return
            self.log("initialised provenance git repo at %s" % cfg.output_root)
        gi = os.path.join(cfg.output_root, ".gitignore")
        if not os.path.exists(gi):
            with open(gi, "w") as f:
                f.write(self.GITIGNORE)

    def _git_commit(self, message):
        """Commit the provenance files for this run.  Paths are added with -f so
        the commit contents never depend on .gitignore guessing right."""
        cfg = self.cfg
        if not cfg.git_repo:
            return
        paths = [os.path.join(cfg.output_root, ".gitignore"), cfg.snapshot_dir,
                 os.path.join(cfg.run_dir, "read_reports.at")]
        for name in ("manifest.txt", "flags.txt", "flags_changed.txt",
                     "main.log", "poison.txt", "summary.txt", "stats.json"):
            paths.append(os.path.join(cfg.logs_dir, name))
        if cfg.definition_file:
            paths.append(cfg.definition_file)
        paths = [p for p in paths if os.path.exists(p)]
        if not paths:
            return
        r = self._git("add", "-f", "--", *paths)
        if r is None or r.returncode != 0:
            self.log("git add failed: %s" % (r.stderr.strip()[:300] if r else "?"))
            return
        r = self._git("commit", "-q", "-m", message)
        if r is not None and r.returncode not in (0, 1):    # 1 = nothing to commit
            self.log("git commit failed: %s" % r.stderr.strip()[:300])
        else:
            self.log("git: %s" % message)

    def _atlas_oneshot(self, extra_args, script, timeout=None, quiet=False):
        """Run a short atlas job: load all.at + preload + extra_args, feed
        `script`, return (stdout_text, returncode).

        quiet=True suppresses the stderr report: the settings probe deliberately
        asks for names that may not exist, so its "Undefined identifier" errors
        are the expected outcome, not something worth putting in main.log."""
        args = [ATLAS_EXE, "all.at"] + self.cfg.preload_files + extra_args
        if timeout is None:
            # A --keep-init restart loads an init that has grown with the run:
            # e8q_init.at reached 31.4 MB (211,017 parameters) against the
            # 7.96 MB reference, and a flat 600 s was not enough to parse it,
            # which killed the restart at setup.  Allow 600 s plus 120 s per MB
            # of everything being loaded.
            mb = sum(os.path.getsize(a) for a in extra_args
                     if isinstance(a, str) and os.path.exists(a)) / 1e6
            timeout = 600.0 + 120.0 * mb
        proc = subprocess.Popen(args, cwd=SCRIPTS_DIR, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(script.encode("utf-8"), timeout=timeout)
        if err and not quiet:
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

    def _count_unitary(self):
        """(params in the rebuilt init, params the workers wrote).

        The first is the size of the shared big_unitary_hash -- one long_match
        line per parameter in the init that write() produced.  The second counts
        what the workers emitted (a tempParams[k]:=parameter(...) line each, or
        long_match on older data).  They should agree; if they do not, a merge
        lost something."""
        cfg = self.cfg
        in_init = written = None
        try:
            with open(cfg.init_file, errors="replace") as f:
                in_init = sum(1 for line in f if "long_match" in line)
        except OSError:
            pass
        data_files = sorted(glob.glob(os.path.join(cfg.run_dir, "[0-9]*.at")))
        if data_files:
            r = subprocess.run(["grep", "-c", "-e", r"tempParams\[", "-e",
                                "long_match"] + data_files,
                               capture_output=True, text=True)
            total = 0
            for ln in r.stdout.splitlines():
                num = ln.rsplit(":", 1)[-1] if ":" in ln else ln
                try:
                    total += int(num)
                except ValueError:
                    pass
            written = total
        return in_init, written

    def write_stats(self, final=False):
        """Publish logs/stats.json -- the machine-readable state of this run.

        Written on every governor tick as well as at the end, so a run that is
        interrupted or dies can still be summarised with --report."""
        if self.shard is not None:
            return self._write_shard_stats(final)
        cfg = self.cfg
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        if final:
            for w in self.workers:
                w.sample_cpu()
            self.unitary_params, self.unitary_params_written = self._count_unitary()
        data = {
            "run": "%s_%d" % (cfg.name, cfg.round),
            "name": cfg.name, "round": cfg.round,
            "started": self.started_at,
            "compute_started": self.compute_started_at or None,
            "finished": self.finished_at or None,
            "updated": time.time(),
            "workers": cfg.workers,
            "queue_total": self.queue_total,
            "done": self.done_count,
            "queue_remaining": self._pairs_remaining(),
            "skipped_already_finished": self.skipped_count,
            "requeued": self.requeued_count,
            "quarantined": self.quarantined_count,
            "unitary_params": self.unitary_params,
            "unitary_params_written": self.unitary_params_written,
            "peak_total_rss_gb": self.peak_total_rss,
            "caps": {"max_total_gb": cfg.max_total_gb,
                     "max_proc_gb": cfg.max_proc_gb},
            # RUSAGE_CHILDREN covers every process the driver reaped: the
            # workers, the merge jobs and the one-shot queries.
            "child_cpu": {"user": ru.ru_utime, "system": ru.ru_stime,
                          "maxrss_gb": ru.ru_maxrss / (1024.0 ** 2)},
            "merges": self._merge_stats(),
            "workers_detail": self._worker_detail(),
        }
        path = os.path.join(cfg.logs_dir, "stats.json")
        tmp = path + ".tmp"
        try:
            with open(tmp, "w") as f:
                json.dump(data, f, indent=1)
            os.replace(tmp, path)
        except OSError as e:
            self.log("stats.json write failed: %r" % e)

    def _write_summary(self):
        self.write_stats(final=True)
        text = build_summary(self.cfg.run_dir)
        path = os.path.join(self.cfg.logs_dir, "summary.txt")
        with open(path, "w") as f:
            f.write(text)
        for ln in text.splitlines():
            self.log("| " + ln)
        return path

    def _write_report_reader(self):
        """Write <run_dir>/read_reports.at.

        Each worker appends a `void:add(<reportDatum>)` record per pair -- the
        timings and counter values for that pair.  Reading them back needs rA
        and add, which the rebuilt init does not carry.  Load this file after
        the init and then the data files, e.g.

            cd atlas-scripts
            ../atlas all.at report.at FPP.at <name>_init.at \
                     <run_dir>/read_reports.at <run_dir>/0.at ...
            summary(rA)
            show_sorted(rA, sorted_xl_pair_numbers(rA), 0, 20)
        """
        cfg = self.cfg
        path = os.path.join(cfg.run_dir, "read_reports.at")
        with open(path, "w") as f:
            f.write(
                "{load after the init file, before the per-worker data files,\n"
                " to read the void:add(...) timing records they contain}\n"
                "set rA=make_report(number_xl_pairs(%s))\n"
                "set add(reportDatum rd)=void:rA[rd.xl_pair_number]:=rd\n"
                % GROUP_SYMBOL)
        self.log("wrote %s" % os.path.basename(path))

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
        # Hand out runs of consecutive indices, not single pairs.  xl_pair is
        # memoized per x and E7_s averages 97 lambdas per x, so consecutive
        # indices almost always hit the memo: measured 0.0011 s against 0.0860 s
        # for scattered indices, i.e. the whole per-pair overhead.  A shared
        # queue of single pairs sends consecutive indices to 784 different
        # workers and every one of them pays the cold cost.
        chunk = cfg.chunk
        if chunk <= 0:                   # ~20 chunks per worker: locality, but
            chunk = min(256, max(1,      # still fine-grained enough to balance
                                 len(entries) // max(1, cfg.workers * 20)))
        chunks = [entries[i:i + chunk] for i in range(0, len(entries), chunk)]
        for c in chunks:
            self.queue.put(c)
        self.queue_total = len(entries)
        self.chunk_size = chunk
        self.log("queued %d (x,lambda) pairs in %d chunks of %d consecutive "
                 "indices" % (len(entries), len(chunks), chunk))

    # -- worker thread ----------------------------------------------------- #
    def _pairs_remaining(self):
        """Pairs left to do.  The queue holds chunks now, so qsize() would
        count runs rather than pairs."""
        return max(0, self.queue_total - self.done_count - self.skipped_count
                   - self.quarantined_count)

    def _record_fail(self, k):
        """Return True if pair k should be quarantined (failed too many times)."""
        with self.fail_lock:
            self.fail_counts[k] = self.fail_counts.get(k, 0) + 1
            n = self.fail_counts[k]
        if n > 2:
            with open(os.path.join(self.cfg.logs_dir, "poison.txt"), "a") as f:
                f.write("%d\n" % k)
            self.quarantined_count += 1
            return True
        self.requeued_count += 1
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
        pending = []          # the rest of the consecutive run this worker holds
        idle_logged = False
        chunk_taken = time.time()
        if not self._safe_launch(worker):
            worker.state = STOPPED
            return
        while not self.stop_event.is_set():
            if worker.recycle_event.is_set():
                wlog("recycle requested; quitting atlas cleanly")
                worker.state = RESTARTING
                worker.recycles += 1
                worker.quit()
                worker.recycle_event.clear()
                if not self._safe_launch(worker):
                    break
            if not pending:
                with self.priority_lock:
                    got = self.priority.popleft() if self.priority else None
                if got is not None:
                    pending = [got]
                    chunk_taken = time.time()
                    idle_logged = False
                try:
                    if not pending:
                        pending = list(self.queue.get_nowait())
                except Empty:
                    # An empty queue does not mean the work is done: another
                    # worker may still be holding a run it is about to hand
                    # back.  On e7s_6 every worker went home the moment the
                    # queue drained, so when the straggler returned its tail
                    # there was nobody left to take it and one process ground
                    # through 18 pairs alone.  Wait until every pair is
                    # accounted for.
                    if self._pairs_remaining() <= 0 and not self.priority:
                        wlog("all pairs accounted for; worker exiting")
                        break
                    if not idle_logged:
                        wlog("queue empty but %d pairs outstanding; standing by"
                             % self._pairs_remaining())
                        idle_logged = True
                    worker.state = IDLE
                    time.sleep(0.25)
                    continue
                if got is None:
                    idle_logged = False
                    chunk_taken = time.time()
            k = pending.pop(0)
            held = time.time() - chunk_taken
            slow = held > self.cfg.giveback_after
            if pending and (slow or (self.queue.empty() and not self.priority)):
                # Whoever holds an expensive run would otherwise grind through
                # it alone while the rest of the machine drains.  On e7s_5 that
                # cost 70 minutes with one worker busy and 783 idle; the worst
                # run of 129 consecutive indices is 220 minutes of serial work.
                #
                # Waiting for the queue to drain is too late: on e7s_7 that
                # happened 13 minutes in, and the last pair to finish (697 s of
                # work) had been sitting in a held run until 15.4 minutes.  So
                # also give back on elapsed time.  A typical run of 129 pairs
                # takes about 26 s, well under the threshold, so only genuinely
                # expensive runs are broken up and the locality that chunking
                # buys is kept for the bulk of the queue.
                giving = len(pending)
                with self.priority_lock:
                    self.priority.extend(pending)
                pending = []
                wlog("returned %d pairs of my run (%s after %s)" % (
                    giving, "slow" if slow else "queue drained", fmt_dur(held)))
            try:
                worker.state = IDLE
                if self.cfg.fused:
                    t0 = time.time()
                    worker.state = COMPUTING
                    worker.compute_started = t0
                    worker.low_cpu_since = 0.0
                    status, dt, grew = worker.do_pair(k)
                    worker.state = IDLE
                    if status == "skip":
                        wlog("pair %d already finished; skipping" % k)
                        worker.pairs_skipped += 1
                        with self.done_lock:
                            self.skipped_count += 1
                        continue
                else:
                    if worker.is_finished(k):
                        wlog("pair %d already finished; skipping" % k)
                        worker.pairs_skipped += 1
                        with self.done_lock:
                            self.skipped_count += 1
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
                    worker.faults += 1
                    if not self._record_fail(k):
                        self.queue.put([k])
                    else:
                        wlog("pair %d quarantined (poison)" % k)
                    worker.state = RESTARTING
                    worker.quit()
                    if not self._safe_launch(worker):
                        break
                    continue
                worker.record_pair(k, dt, t0, grew)
                with self.done_lock:
                    self.done_count += 1
                    done = self.done_count
                wlog("pair %d done in %s (total done=%d, left=%d)" % (
                    k, fmt_dur(dt), done, self._pairs_remaining()))
                # trigger an init merge periodically
                if done % 25 == 0:
                    self._request_merge()
                # fold in what other workers have found since the last pair.
                # Done between pairs, never during one, so a worker's hash only
                # changes at a point where nothing is part-computed.
                if self.cfg.share and self.share_sem.acquire(blocking=False):
                    # non-blocking: a worker that cannot get a slot keeps
                    # computing and tries again after its next pair, so the cap
                    # never leaves anyone idle waiting to absorb
                    try:
                        n = worker.apply_broadcast()
                    finally:
                        self.share_sem.release()
                    if n and worker.bcast_applies % 20 == 1:
                        wlog("applied %.1f MB of shared parameters "
                             "(%.1f MB total in %s)" % (
                                 n / 1e6, worker.bcast_bytes / 1e6,
                                 fmt_dur(worker.bcast_seconds)))
            except (HangError, WorkerDied) as e:
                worker.state = RESTARTING
                worker.faults += 1
                wlog("worker fault on pair %d: %r -> restarting" % (k, e))
                worker.kill()
                if self._record_fail(k):
                    wlog("pair %d quarantined (poison)" % k)
                else:
                    self.queue.put([k])      # retry later
                # always merge after a fault: salvage shared progress
                self._request_merge()
                worker.recycle_event.clear()
                if not self._safe_launch(worker):
                    break
            except Exception as e:           # noqa: BLE001
                wlog("UNEXPECTED error on pair %d: %r" % (k, e))
                worker.faults += 1
                self.requeued_count += 1
                self.queue.put([k])
                worker.state = RESTARTING
                worker.kill()
                if not self._safe_launch(worker):
                    break
        worker.state = STOPPED
        worker.quit()
        worker.last_stop = time.time()
        wlog("stopped after %d pairs, %s computing, %s cpu" % (
            worker.pairs_done, fmt_dur(worker.compute_seconds),
            fmt_dur(worker.cpu_total)))

    # -- governor ---------------------------------------------------------- #
    def governor_thread(self):
        cfg = self.cfg
        # prime cpu_percent counters
        procs = {}
        while not self.stop_event.is_set():
            time.sleep(cfg.poll_interval)
            sizes = []          # (rss_gb, worker)
            total = 0.0
            live = set()
            for w in self.workers:
                pid = w.pid
                if pid is None:
                    continue
                live.add(pid)
                try:
                    p = procs.get(pid)
                    if p is None or p.pid != pid:
                        p = psutil.Process(pid)
                        procs[pid] = p
                        p.cpu_percent(None)      # prime
                    rss = p.memory_info().rss / GB
                    cpu = p.cpu_percent(None)
                    times = p.cpu_times()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
                w.note_cpu(pid, times.user + times.system)
                if rss > w.peak_rss:
                    w.peak_rss = rss
                total += rss
                sizes.append((rss, w))
                self._check_stall(w, cpu)
            self._check_stop_floor()
            # workers are recycled a lot on a long run; don't hoard dead pids
            if len(procs) > len(live):
                procs = {k: v for k, v in procs.items() if k in live}
            if total > self.peak_total_rss:
                self.peak_total_rss = total
            self._enforce_memory(total, sizes)
            # The governor polls often (to catch memory spikes early) but logging
            # every poll buries everything else; log on a slower cadence.
            if time.time() - self._last_mem_log >= cfg.mem_log_interval:
                self._last_mem_log = time.time()
                biggest = max(s[0] for s in sizes) if sizes else 0.0
                self.log("MEM total=%.1f GB across %d workers (max %.2f GB); "
                         "left=%d done=%d" % (total, len(sizes), biggest,
                                               self._pairs_remaining(), self.done_count))
                # keep stats.json current so an interrupted or killed run can
                # still be summarised afterwards
                self.write_stats()

    def _check_stop_floor(self):
        """Stop once few enough pairs are left MACHINE-WIDE, not per shard.

        The last pairs of a run are its most expensive: on E8_q the final 200
        cost 8h19m of a 19h27m run.  When the outstanding count reaches the
        floor we stop dispatching and kill whatever is still computing, so the
        shards unwind and the merger writes the init exactly as at a normal
        finish; the rest can be completed later with --keep-init.

        An earlier version applied the floor per shard, which was wrong: each
        shard then ran until IT was nearly done, so the run ended long after the
        machine as a whole was, with the last stretch on a fraction of the
        workers.  On e8q_5 that left 3 shards grinding with 98 workers up
        instead of 784.  Shards are separate processes, so each publishes its
        own count and reads everyone's; no parent involvement is needed."""
        floor = self.cfg.stop_at_remaining
        if floor <= 0 or self._stop_floor_hit:
            return
        mine = self._pairs_remaining()
        me = self.shard if self.shard is not None else 0
        try:
            tmp = os.path.join(self.cfg.logs_dir, "remaining_%d.tmp" % me)
            with open(tmp, "w") as f:
                f.write("%d\n" % mine)
            os.replace(tmp, os.path.join(self.cfg.logs_dir,
                                         "remaining_%d.txt" % me))
        except OSError:
            pass
        # A shard whose governor has died stops updating its file, leaving a
        # stale -- and possibly very large -- value.  Believing it would disable
        # the floor for the whole run, which is exactly what happened on e8q_6:
        # four governors died and 1,291,708 phantom pairs kept the sum above the
        # threshold forever.  Treat anything not updated recently as finished.
        stale_after = max(60.0, 10.0 * self.cfg.poll_interval)
        now = time.time()
        total = 0
        seen = 0
        stale = 0
        for f in glob.glob(os.path.join(self.cfg.logs_dir, "remaining_*.txt")):
            try:
                if now - os.path.getmtime(f) > stale_after:
                    stale += 1
                    seen += 1
                    continue    # count it as zero: that shard is not reporting
                total += int(open(f).read().strip())
                seen += 1
            except (OSError, ValueError):
                return          # a partial read: try again next tick
        if seen < self.cfg.shards:
            return              # not every shard has reported yet
        if total > floor:
            return
        if stale:
            self.log("stop-at-remaining: %d shard file(s) stale (>%.0f s); "
                     "counting them as done" % (stale, stale_after))
        self._stop_floor_hit = True
        self.log("stop-at-remaining: %d pairs left machine-wide across %d "
                 "shards (floor %d, this shard holds %d); stopping and "
                 "abandoning what is still running"
                 % (total, seen, floor, mine))
        self.stop_event.set()
        for w in self.workers:
            if w.state == COMPUTING:
                self.log("  abandoning job %d (computing)" % w.job)
                w.kill()

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
            sizes.sort(key=lambda t: t[0], reverse=True)
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
            sizes.sort(key=lambda t: t[0], reverse=True)
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
        """Rebuild the shared init periodically.

        Workers ask for a merge often (every 25 pairs, and after every fault),
        but a merge reloads the whole data corpus, so honouring each request
        means merging back to back: on the F4_s test that was 77% of wall time,
        and on E7_s the corpus is orders of magnitude bigger.  Requests are
        therefore only wake-up hints -- the gap between merges is at least
        merge_interval, and at least 4x the duration of the last merge, so
        merging can never take more than ~20% of the run no matter how big the
        data gets.  Nothing is lost by merging less often: each worker's results
        are already durable in its own data file, and the merge only affects how
        soon other workers get to share them."""
        cfg = self.cfg
        while not self.stop_event.is_set():
            pm = self.merger_persistent
            if pm is not None and not pm.failed:
                # incremental: only what the data files have grown by
                if self.stop_event.wait(timeout=cfg.merge_interval):
                    break
                pm.ingest()
                continue
            self.merger.request.wait(timeout=cfg.merge_interval)
            self.merger.request.clear()
            if self.stop_event.is_set():
                break
            gap = max(cfg.merge_interval, 4.0 * self.merger.last_duration)
            if time.time() - self.merger.last_merge >= gap:
                self.merger.merge()

    # -- run --------------------------------------------------------------- #
    def _worker_detail(self):
        return [
            {"job": w.job, "pairs": w.pairs_done, "skipped": w.pairs_skipped,
             "compute_seconds": round(w.compute_seconds, 3),
             "cpu_seconds": round(w.cpu_total, 3),
             "peak_rss_gb": round(w.peak_rss, 4),
             "first_launch": w.first_launch or None,
             "last_stop": w.last_stop or None,
             "recycles": w.recycles, "faults": w.faults}
            for w in self.workers]

    def _request_merge(self):
        """Worker hint that there is new data.  A shard child has no merger of
        its own -- the parent process drives the single persistent one."""
        if self.merger is not None:
            self.merger.request.set()

    def _merge_stats(self):
        """Report whichever merger actually produced the init."""
        pm, bm = self.merger_persistent, self.merger
        m = pm if (pm is not None and pm.wrote_init) else None
        if m is None and bm is not None and bm.count:
            m = bm
        if m is None:
            m = pm or bm
        if m is None:
            return {"count": 0, "seconds": 0.0,
                    "last_input_bytes": 0, "last_init_bytes": 0}
        return {"count": m.count, "seconds": m.seconds,
                "last_input_bytes": m.last_input_bytes,
                "last_init_bytes": m.last_init_bytes,
                "persistent": m is pm}

    def _finalize_merge(self):
        """Write the run's init.  The persistent merger only has the tail of the
        data left to fold, so this is seconds rather than the tens of minutes a
        batch merge over the whole corpus costs."""
        pm = self.merger_persistent
        if pm is not None and not pm.failed:
            ok = pm.finish()
            pm.quit()
            if ok:
                return
        self.log("falling back to a batch merge over the whole corpus")
        if self.merger is None:
            self.merger = Merger(self.cfg, self.log)
        self.merger.merge(blocking=True)

    def _write_shard_stats(self, final=False):
        """A forked shard writes its own slice of the stats; the parent folds
        the shard files into logs/stats.json."""
        cfg = self.cfg
        if final:
            for w in self.workers:
                w.sample_cpu()
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        data = {
            "shard": self.shard,
            "workers": cfg.workers,
            "queue_total": self.queue_total,
            "done": self.done_count,
            "queue_remaining": self._pairs_remaining(),
            "skipped_already_finished": self.skipped_count,
            "requeued": self.requeued_count,
            "quarantined": self.quarantined_count,
            "peak_total_rss_gb": self.peak_total_rss,
            "compute_started": self.compute_started_at or None,
            "updated": time.time(),
            "child_cpu": {"user": ru.ru_utime, "system": ru.ru_stime,
                          "maxrss_gb": ru.ru_maxrss / (1024.0 ** 2)},
            "workers_detail": self._worker_detail(),
        }
        path = os.path.join(cfg.logs_dir, "stats_%d.json" % self.shard)
        tmp = path + ".tmp"
        try:
            with open(tmp, "w") as f:
                json.dump(data, f, indent=1)
            os.replace(tmp, path)
        except OSError as e:
            self.log("stats_%d.json write failed: %r" % (self.shard, e))

    def _aggregate_stats(self, final=False):
        """Fold the shard stats files into logs/stats.json (same schema as the
        unsharded run, so --report and build_summary are unchanged)."""
        cfg = self.cfg
        details = []
        done = skipped = requeued = quarantined = qrem = 0
        peak = cu = cs = mx = 0.0
        for path in sorted(glob.glob(os.path.join(cfg.logs_dir, "stats_*.json"))):
            try:
                with open(path) as f:
                    d = json.load(f)
            except (OSError, ValueError):
                continue
            details.extend(d.get("workers_detail", []))
            done += d.get("done", 0)
            skipped += d.get("skipped_already_finished", 0)
            requeued += d.get("requeued", 0)
            quarantined += d.get("quarantined", 0)
            qrem += d.get("queue_remaining", 0)
            # shards peak at about the same moment, so the sum is a fair estimate
            peak += d.get("peak_total_rss_gb", 0.0)
            c = d.get("child_cpu", {})
            cu += c.get("user", 0.0)
            cs += c.get("system", 0.0)
            mx = max(mx, c.get("maxrss_gb", 0.0))
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        if final:
            # reaped shards fold their whole subtree into our RUSAGE_CHILDREN
            cu, cs = max(cu, ru.ru_utime), max(cs, ru.ru_stime)
            mx = max(mx, ru.ru_maxrss / (1024.0 ** 2))
            self.unitary_params, self.unitary_params_written = self._count_unitary()
        self.done_count = done
        self.skipped_count = skipped
        self.requeued_count = requeued
        self.quarantined_count = quarantined
        self.peak_total_rss = max(self.peak_total_rss, peak)
        data = {
            "run": "%s_%d" % (cfg.name, cfg.round),
            "name": cfg.name, "round": cfg.round,
            "started": self.started_at,
            "compute_started": self.compute_started_at or None,
            "finished": self.finished_at or None,
            "updated": time.time(),
            "workers": cfg.workers,
            "shards": cfg.shards,
            "queue_total": self.queue_total,
            "done": done,
            "queue_remaining": qrem,
            "skipped_already_finished": skipped,
            "requeued": requeued,
            "quarantined": quarantined,
            "unitary_params": self.unitary_params,
            "unitary_params_written": self.unitary_params_written,
            "peak_total_rss_gb": self.peak_total_rss,
            "caps": {"max_total_gb": cfg.max_total_gb,
                     "max_proc_gb": cfg.max_proc_gb},
            "child_cpu": {"user": cu, "system": cs, "maxrss_gb": mx},
            "merges": self._merge_stats(),
            "workers_detail": details,
        }
        path = os.path.join(cfg.logs_dir, "stats.json")
        tmp = path + ".tmp"
        try:
            with open(tmp, "w") as f:
                json.dump(data, f, indent=1)
            os.replace(tmp, path)
        except OSError as e:
            self.log("stats.json write failed: %r" % e)

    # -- run --------------------------------------------------------------- #
    def run(self):
        if self.queue.empty():
            self.log("nothing to do; exiting")
            return
        if self.cfg.shards > 1:
            self._run_sharded()
        else:
            self._run_single()

    def _finish_run(self):
        self.finished_at = time.time()
        self.log("run complete: %d pairs computed" % self.done_count)
        summary = self._write_summary()
        self.log("wrote %s" % summary)
        self._git_commit("%s_%d done: %d pairs computed" % (
            self.cfg.name, self.cfg.round, self.done_count))

    def _run_single(self):
        cfg = self.cfg
        for job in range(cfg.workers):
            self.workers.append(AtlasWorker(job, cfg))
        if self.merger_persistent is not None and not self.merger_persistent.start():
            self.merger_persistent = None
        gov = threading.Thread(target=self.governor_thread, name="governor", daemon=True)
        mrg = threading.Thread(target=self.merge_thread, name="merger", daemon=True)
        gov.start()
        mrg.start()
        self.log("starting %d workers" % cfg.workers)
        self.compute_started_at = time.time()
        for w in self.workers:
            t = threading.Thread(target=self.worker_thread, args=(w,),
                                 name="job-%d" % w.job)
            t.start()
            self.threads.append(t)
            time.sleep(cfg.stagger)
        for t in self.threads:
            while t.is_alive():
                t.join(timeout=1.0)
        self.log("all workers finished; folding the tail of the data")
        self.stop_event.set()
        self._finalize_merge()
        self._finish_run()

    def _run_sharded(self):
        """Split the workers over several driver processes.

        One driver process cannot feed 784 workers: measured on e7s_3 the gap
        between a worker finishing one pair and being given the next was 0.884 s
        against 0.407 s of actual compute, and workers were busy only 34% of the
        compute phase.  The same atlas processes driven one-per-python-process
        show a 0.0086 s gap even at 784-way concurrency, so the cost is GIL
        contention among 784 driver threads, not the pairs and not the hardware
        (which only charges 1.32x at that concurrency)."""
        cfg = self.cfg
        entries = []
        while True:
            try:
                entries.append(self.queue.get_nowait())
            except Empty:
                break
        P = max(1, min(cfg.shards, cfg.workers))
        slices = [entries[i::P] for i in range(P)]
        counts = [cfg.workers // P + (1 if i < cfg.workers % P else 0)
                  for i in range(P)]
        bases, b = [], 0
        for c in counts:
            bases.append(b)
            b += c
        if self.merger_persistent is not None and not self.merger_persistent.start():
            self.merger_persistent = None
        self.log("sharding %d pairs (%d chunks) over %d driver processes "
                 "(%d-%d workers each, %.0f GB cap each)" % (
                     sum(len(c) for c in entries), len(entries), P,
                     min(counts), max(counts), cfg.max_total_gb / P))
        self.compute_started_at = time.time()
        pids = {}
        for i in range(P):
            self.main_log_fh.flush()
            pid = os.fork()
            if pid == 0:
                os._exit(self._child_main(i, P, slices[i], counts[i], bases[i]))
            pids[pid] = i
        self.child_pids = list(pids)
        last_ingest = time.time()
        last_stats = 0.0
        while pids:
            try:
                pid, status = os.waitpid(-1, os.WNOHANG)
            except ChildProcessError:
                break
            if pid:
                i = pids.pop(pid, None)
                self.log("shard %s (pid %d) exited, status %d" % (i, pid, status))
                continue
            time.sleep(1.0)
            t = time.time()
            pm = self.merger_persistent
            if (pm is not None and not pm.failed
                    and t - last_ingest >= cfg.merge_interval):
                pm.ingest()
                last_ingest = time.time()
            if t - last_stats >= cfg.poll_interval:
                self._aggregate_stats()
                last_stats = t
        self.stop_event.set()
        self.log("all shards finished; folding the tail of the data")
        self._finalize_merge()
        self.finished_at = time.time()
        self._aggregate_stats(final=True)
        text = build_summary(cfg.run_dir)
        with open(os.path.join(cfg.logs_dir, "summary.txt"), "w") as f:
            f.write(text)
        for ln in text.splitlines():
            self.log("| " + ln)
        self.log("run complete: %d pairs computed" % self.done_count)
        self._git_commit("%s_%d done: %d pairs computed" % (
            cfg.name, cfg.round, self.done_count))

    def _child_main(self, shard, nshards, entries, nworkers, job_base):
        """Body of a forked shard: its own workers, governor and stats file."""
        try:
            if self.merger_persistent is not None:
                self.merger_persistent.detach()   # the parent owns that process
                self.merger_persistent = None
            self.merger = None
            self.shard = shard
            self.cfg = dc_replace(self.cfg, workers=nworkers,
                                  max_total_gb=self.cfg.max_total_gb / nshards)
            cfg = self.cfg
            self.queue = Queue()
            self.priority = deque()
            for c in entries:
                self.queue.put(c)
            self.queue_total = sum(len(c) for c in entries)
            self.workers = [AtlasWorker(job_base + j, cfg) for j in range(nworkers)]
            self.threads = []
            self.log("shard %d: %d pairs in %d chunks, jobs %d-%d, %.0f GB cap" % (
                shard, self.queue_total, len(entries), job_base,
                job_base + nworkers - 1, cfg.max_total_gb))
            gov = threading.Thread(target=self.governor_thread,
                                   name="governor", daemon=True)
            gov.start()
            self.compute_started_at = time.time()
            # keep the *global* launch rate at cfg.stagger even though the shards
            # launch in parallel: pairs landing in the first minute of e7s_3 cost
            # 2.87x, decaying to the 1.32x steady state only after ~5 minutes
            time.sleep(shard * cfg.stagger)
            for w in self.workers:
                t = threading.Thread(target=self.worker_thread, args=(w,),
                                     name="job-%d" % w.job)
                t.start()
                self.threads.append(t)
                time.sleep(cfg.stagger * nshards)
            for t in self.threads:
                while t.is_alive():
                    t.join(timeout=1.0)
            self.stop_event.set()
            self._write_shard_stats(final=True)
            self.log("shard %d done: %d pairs computed" % (shard, self.done_count))
            return 0
        except BaseException as e:                   # noqa: BLE001
            try:
                self.log("shard %s CRASHED: %r" % (shard, e))
            except Exception:                        # noqa: BLE001
                pass
            return 1

    def shutdown(self, signum, _frame):
        self.log("signal %d received; draining workers gracefully" % signum)
        self.stop_event.set()
        for pid in self.child_pids:
            try:
                os.kill(pid, signum)
            except OSError:
                pass
        for w in self.workers:
            w.recycle_event.clear()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_config(argv):
    p = argparse.ArgumentParser(
        description="Parallel FPP unitary-dual driver.",
        epilog="Reporting: %(prog)s --report <run_dir> rebuilds and prints "
               "<run_dir>/logs/summary.txt for any run, finished or not.")
    p.add_argument("--preset", choices=sorted(PRESETS), help="parameter bundle: f4s (small test), e7s, e8q, e8s")
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
    p.add_argument("--shards", type=int, default=None,
                   help="driver processes to split the workers over "
                        "(default: one per ~50 workers; 1 disables sharding)")
    p.add_argument("--stagger", type=float, default=0.05,
                   help="seconds between worker launches, globally")
    p.add_argument("--giveback-after", type=float, default=60.0,
                   help="seconds before a worker hands the rest of its run "
                        "back to the queue so idle workers can help")
    p.add_argument("--share", action="store_true",
                   help="broadcast newly found parameters to live workers, so "
                        "each worker's hash grows with the run instead of "
                        "holding only its own discoveries")
    p.add_argument("--share-cap-mb", type=float, default=4.0,
                   help="most a worker applies in one go (default 4 MB)")
    p.add_argument("--share-concurrency", type=int, default=32,
                   help="how many workers may absorb shared parameters at once, "
                        "machine-wide (default 32).  Absorbing costs 44 s alone "
                        "but 1430 s when all 784 do it together, so this cap is "
                        "what makes sharing affordable")
    p.add_argument("--stop-at-remaining", type=int, default=None,
                   help="stop gracefully once this many pairs are still "
                        "outstanding, instead of waiting for the last few.  On "
                        "E8_q the final 200 pairs cost 8h19m of a 19h27m run "
                        "-- 43%% of the wall clock for 0.003%% of the work.  The "
                        "abandoned pairs can be finished later with "
                        "--keep-init.  0 disables.")
    p.add_argument("--chunk", type=int, default=0,
                   help="consecutive pair indices handed to a worker at once "
                        "(0 = auto, ~20 chunks per worker)")
    p.add_argument("--no-fused", action="store_true",
                   help="send is_finished/compute/write as three separate "
                        "commands instead of one do_one_pair call")
    p.add_argument("--no-persistent-merge", action="store_true",
                   help="use the old batch merge instead of folding worker "
                        "output into a live hash as the run proceeds")
    p.add_argument("--keep-init", action="store_true",
                   help="do NOT reset the init file from the reference")
    p.add_argument("--settings", default="fpp_settings.at",
                   help="settings file loaded last to override defaults "
                        "(debugging flags); default fpp_settings.at")
    p.add_argument("--no-settings", action="store_true",
                   help="do not load any settings/override file")
    p.add_argument("--scripts-dir",
                   help="atlas-scripts directory (default: found next to this "
                        "script, or $FPP_SCRIPTS_DIR)")
    p.add_argument("--preload", action="append",
                   help="file loaded after all.at and before the init "
                        "(repeatable; default: %s)" % " ".join(DEFAULT_PRELOAD))
    p.add_argument("--no-preload", action="store_true",
                   help="load nothing between all.at and the init file")
    p.add_argument("--no-git", action="store_true",
                   help="do not keep run provenance in a git repo at the "
                        "output root")
    p.add_argument("--mem-log-interval", type=float, default=60.0,
                   help="seconds between governor memory lines in main.log")
    a = p.parse_args(argv)

    global SCRIPTS_DIR, ATLAS_DIR, ATLAS_EXE
    if a.scripts_dir:
        SCRIPTS_DIR = find_scripts_dir(a.scripts_dir)
        ATLAS_DIR = os.path.dirname(SCRIPTS_DIR)
        ATLAS_EXE = os.path.join(ATLAS_DIR, "atlas")
    if not os.path.exists(ATLAS_EXE):
        p.error("atlas executable not found: %s" % ATLAS_EXE)

    def resolve(name):
        """Absolute paths as given; relative ones looked for beside the driver,
        in the current directory, then in atlas-scripts."""
        if os.path.isabs(name):
            return name
        for base in (os.getcwd(), HERE, SCRIPTS_DIR):
            cand = os.path.join(base, name)
            if os.path.exists(cand):
                return os.path.abspath(cand)
        return os.path.join(SCRIPTS_DIR, name)

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
    reference = resolve(reference)
    if not os.path.exists(reference):
        p.error("reference init not found: %s" % reference)
    aux = a.aux if a.aux else (pre.aux_files if pre else [])
    aux = [resolve(f) for f in aux]
    for f in aux:
        if not os.path.exists(f):
            p.error("aux file not found: %s" % f)
    preload = [] if a.no_preload else [resolve(f) for f in (a.preload or DEFAULT_PRELOAD)]
    for f in preload:
        if not os.path.exists(f):
            p.error("preload file not found: %s" % f)
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
        settings = resolve(a.settings)
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
        preload_files=preload,
        git_repo=not a.no_git,
        mem_log_interval=a.mem_log_interval,
        shards=(a.shards if a.shards is not None
                else max(1, (workers + 49) // 50)),
        stagger=a.stagger,
        persistent_merge=not a.no_persistent_merge,
        fused=not a.no_fused,
        chunk=a.chunk,
        giveback_after=a.giveback_after,
        share=a.share,
        share_cap_mb=a.share_cap_mb,
        share_concurrency=a.share_concurrency,
        # explicit flag wins, including an explicit 0 to disable a preset's
        # floor -- which is what a --keep-init run finishing a handful of
        # stragglers needs, since the preset's 200 would fire immediately
        stop_at_remaining=(a.stop_at_remaining if a.stop_at_remaining is not None
                           else pick(None, "stop_at_remaining", 0)),
    )


def main(argv):
    if "--report" in argv:
        i = argv.index("--report")
        if i + 1 >= len(argv):
            sys.exit("--report needs a run directory")
        run_dir = argv[i + 1]
        text = build_summary(run_dir)
        try:
            with open(os.path.join(run_dir, "logs", "summary.txt"), "w") as f:
                f.write(text)
        except OSError as e:
            sys.stderr.write("could not write summary.txt: %r\n" % e)
        sys.stdout.write(text)
        return
    cfg = build_config(argv)
    run = FPPRun(cfg)
    run.setup(argv)
    signal.signal(signal.SIGINT, run.shutdown)
    signal.signal(signal.SIGTERM, run.shutdown)
    run.run()


if __name__ == "__main__":
    main(sys.argv[1:])
