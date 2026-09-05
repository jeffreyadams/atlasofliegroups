# E8_q run plan

> **STATUS 2026-09-05: COMPLETE.** All 5,926,140 pairs are finished and the
> answer is **211,017 unitary parameters**, held in
> `atlas-scripts/e8q_init.at` (31.4 MB, written 2026-09-05 16:56 by e8q_9).
> See [Outcome](#outcome) at the end. The plan below is kept as written on
> 2026-08-25. **Do not launch it again**: a bare `--preset e8q` run without
> `--keep-init` resets `e8q_init.at` from the reference and overwrites the
> answer. A byte-identical copy is `e8q_init.at.FINAL-20260905-1656`.
Prepared 2026-08-25 from direct measurement. Everything below was measured on
euphrates today, not estimated from E7_s.

## Launch command

```
cd /u02/jdada11/claudeFPP
./fpp_claude.py --preset e8q --shards 16 --stagger 0.08 --merge-interval 60 --chunk 1
```

Run it inside `screen` (or tmux) — this is an overnight job:

```
screen -S e8q
cd /u02/jdada11/claudeFPP
./fpp_claude.py --preset e8q --shards 16 --stagger 0.08 --merge-interval 60 --chunk 1
```
Detach with `Ctrl-A D`; reattach later with `screen -r e8q`.

The preset now supplies the rest: 784 workers, 15000 GB total cap, 60 GB
per-worker recycle trigger, the reference init and the edges file.

Note there is no `--reverse`. See "Ordering" below.

## The work

| | E8_q | E7_s |
|---|---|---|
| KGB size | 67,110 | 20,926 |
| (x,lambda) pairs | 5,926,140 | 2,025,524 |
| mean lambda per x | 88.3 | 96.8 |
| parameters already known | 49,881 | 33,160 |
| pairs already finished | **0** | 0 |

All 5,926,140 pairs must be computed. The 49,881 known representations seed
`big_unitary_hash` but mark no pair complete (`TODO=5926140`, and all 67,110
entries of the `G_temp_XLK` map are zero).

## Files loaded

In order, as the driver passes them to every atlas process:

| file | size | why |
|---|---|---|
| `all.at` | — | the usual prelude |
| `report.at` | 12 KB | `report_datum`, `make_report` |
| `FPP.at` | 8 KB | `write_one_pair`, `do_one_pair` |
| `e8qinitreference.at` | 7.96 MB | the reference hash + `set_xl_sizes` + XLK map |
| `edges_F_E.at` | 984.87 MB | edge data; the E8 block is 99% of it |
| `fpp_settings.at` | 687 B | flags, incl. `jeff_sizes_flag:=true` |

`edges_F_E.at` is taken from
`/u02/jdada11/atlasSoftware/to_ht_branch_jeff_2/atlas-scripts/` by absolute path
(it is not in the working atlas-scripts, and at 985 MB is not worth copying).

### Deliberately NOT loaded

`E8qCohIndPlus.at` and `E8qCohIndUnipPlus.at`. Comparing every `parameter(...)`
argument:

```
E8qCohIndUnipPlus.at   49,708 params  -> a strict subset of E8qCohIndPlus.at
E8qCohIndPlus.at       49,881 params  -> identical set to e8qinitreference.at
union of all three     49,881 params
```

`e8qinitreference.at` already contains everything in both, in a different order
(the signature of having been produced by `big_unitary_hash.write()`). Loading
either CohInd file would cost ~13 MB of parsing to insert duplicates the hash
discards. The 173 parameters in CohIndPlus but not in CohIndUnipPlus are the
ones not reachable from unipotent; they are in the reference too.

## Processors and memory

**784 workers, 16 driver shards (49 workers each).** 784 is half the machine's
1568 CPUs, the limit set for E7_s. Sixteen shards keeps ~49 threads per driver
process; one process driving 784 workers was the bottleneck in e7s_3, costing
0.884 s of idle per pair against 0.0086 s for a dedicated driver.

**Memory: ~4.6 TB expected, cap 15 TB.**

| | measured |
|---|---|
| atlas + all.at + reference init | 0.29 GB, 12.6 s |
| + `edges_F_E.at` | 6.02 GB peak, 166 s total |
| 784 workers at that level | 4.61 TB |

Headroom is deliberate: the probe processes computed 400 pairs each, whereas a
real worker will do ~7,560 (5,926,140 / 784). E7_s workers grew about 1.5 GB
over 2,600 pairs; if E8_q grows similarly per pair, expect 8-11 GB per worker,
so 6-9 TB. Still inside the cap, and `--max-proc-gb 60` plus the global
governor will recycle workers if it goes further.

**Disk:** the run directory goes to `/scr/jdada11/fpp_claude/e8q_1` (94 TB
free). Expect 20-30 GB: ~2.2 GB of parameter data plus per-job verbose logs,
which ran ~2.5 KB per pair on E7_s.

## Expected duration

From 1,600 random pairs computed cold in four parallel processes:

```
mean 4.086 s   median 1.445 s   p90 7.62 s   p99 26.98 s   max 715.9 s
bootstrap 95% CI on the mean: 3.19 - 5.31 s
```

| | core-hours | wall on 784 workers at 88% |
|---|---|---|
| low  | 5,254 | 7.6 h |
| **point** | **6,727** | **9.8 h** |
| high | 8,743 | 12.7 h |

**Plan for about 10 hours, and do not be alarmed by 13.** For comparison E7_s
was 254 core-hours and 24 minutes, so this is roughly 26x the computation.

The uncertainty is irreducible at this sample size because the cost is
concentrated in the tail: **the top 1% of pairs carry 32% of the total, the top
10% carry 63%.** The estimated mean moved from 3.81 s (40 samples) to 2.13 s
(713) to 4.09 s (1600) as single expensive pairs entered and were diluted.

## Ordering: no `--reverse`

For E7_s, cost rose steeply with pair index — the top decile held 27% of the
compute, the bottom 2.9% — so `--reverse` was an effective way to start the
expensive work first. **That does not hold for E8_q:**

```
Pearson r(index, cost) = 0.056        Spearman rank r = 0.091
index band share of cost: 7.9%, 22.2%, 20.1%, 22.1%, 27.6%
```

Cost is spread almost evenly across the index range, so reversing gains
nothing. Expensive pairs will surface throughout the run, including near the
end.

## Dispatch: `--chunk 1`

Single-pair dispatch, i.e. the e7s_4 configuration, which remains the best
wall-clock result on E7_s (24:09). Chunking consecutive indices saves 31% of
worker CPU through `xl_pair` memoisation, but on E7_s it cost a factor of five
in wall time and four scheduling fixes never fully recovered it.

There is a real argument that chunking would behave better here — E7_s's
failure came from all the expensive pairs clustering in a few chunks, and E8_q's
cost is index-independent, so chunks would have near-average cost. If this first
run succeeds, a second with `--chunk 129` is worth about 30% and the driver now
carries the protections (60-second give-back, priority requeue, workers waiting
for the tail). But a first run should be the configuration known to work.

## What to watch

1. **The tail.** The worst sampled pair is 715.9 s, from 0.027% of the pairs.
   The true worst is likely much larger — possibly hours — and with no
   index/cost correlation it may start late. If the run sits at 99.9% with one
   worker busy, that is why, and it is not a bug.
2. **Memory growth per worker**, the least-constrained extrapolation here.
   `logs/main.log` prints a `MEM total=...` line per shard per minute.
3. **The merger keeping up.** `grep "merger ingested" logs/main.log`. It folded
   757.8 MB during E7_s with the final ingest at 0.0 MB; E8_q is roughly 3x
   that and should still keep pace across a 10-hour run.

## Monitoring

```
D=/scr/jdada11/fpp_claude/e8q_1
tail -f $D/logs/main.log
/usr/bin/python3 -c "import json;d=json.load(open('$D/logs/stats.json'));\
print('%d/%d  %.1f%%'%(d['done'],d['queue_total'],100*d['done']/d['queue_total']))"
```

`logs/stats.json` is refreshed every governor tick, so progress is available at
any moment and an interrupted run can still be summarised with
`./fpp_claude.py --report $D`.

## If it has to be stopped

`Ctrl-C` (or SIGTERM to the driver) drains workers gracefully and the persistent
merger still writes the init. To continue afterwards:

```
./fpp_claude.py --preset e8q --keep-init --shards 16 --chunk 1
```

`--keep-init` keeps the existing `e8q_init.at` instead of resetting from the
reference, and `xl_pairs_todo` then returns only the pairs still outstanding.
Without that flag the run starts from the reference again — which is the correct
default, and is how every E7_s run began.

## Verifying when it finishes

```
cd /u02/jdada11/claudeFPP/atlasofliegroups/atlas-scripts
printf 'prints("TODO_LEFT=",#xl_pairs_todo(big_unitary_hash,G_temp))\n\
prints("UHASH=",big_unitary_hash.uhash(G_temp).size())\nquit\n' | \
  ../atlas all.at report.at FPP.at e8q_init.at fpp_settings.at
```

Expect `TODO_LEFT=0`. `UHASH` is the answer — unknown in advance; E7_s found
237,641 unitary parameters from 2,025,524 pairs.

**Result, run 2026-09-05 17:35 on the final init:** `TODO_LEFT=0`, `UHASH=211017`.
The query took 13 min 12 s and 129 GB, almost all of it parsing the init.

## Outcome

E8_q finished on 2026-09-05 at 16:56. Two runs produced the answer:

| run | when | workers | pairs | pair compute | worker CPU | wall |
|---|---|---|---|---|---|---|
| e8q_1 | 08-25 18:32 to 08-26 14:00 | 784 | 5,926,131 | 6,882 core-h | 8,453 core-h | 19 h 28 m |
| e8q_9 | 09-03 12:48 to 09-05 16:56 | 16 | 9 | 218 core-h | 221 core-h | 52 h 08 m |
| **total** | | | **5,926,140** | **7,100 core-h** | **8,674 core-h** | **71.6 h** |

**The answer is 211,017 unitary parameters** (E7_s: 237,641 from 2,025,524
pairs; E8_q has nearly three times the pairs and fewer unitary parameters).
e8q_1 found all of them. The nine pairs e8q_9 completed contributed none: the
only change the merger made to `e8q_init.at` was the completion mask of six
`G_temp_XLK` entries, checked by diffing against the pre-e8q_9 copy.

Where it lives:

| file | what |
|---|---|
| `atlas-scripts/e8q_init.at` | **the answer**: 31,432,187 bytes, mtime 2026-09-05 16:56 |
| `atlas-scripts/e8q_init.at.FINAL-20260905-1656` | byte-identical safety copy, made 2026-09-05 |
| `e8q_init.at.SAFE-20260828-0437`, `.SAFE-20260903-1248` | the pre-e8q_9 state: same 211,017 parameters, nine pairs still open; identical to each other, superseded |
| `e8q_summary_20260826.pdf` | write-up of e8q_1 |
| `fpp_sharing_20260826.pdf` | the parameter-sharing A/B (e8q_7 vs e8q_8) |
| `fpp_report_e8q_20260826.pdf` | the E8_s measurement made from the E8_q machinery |

### What happened between the plan and the answer

- **e8q_1 (production).** All but 38 pairs were done at 12 h 16 m; those 38,
  every one from the top of the KGB, held 784 workers for seven more hours
  while 29 of them finished. Stopped deliberately at 19 h 28 m with nine pairs
  abandoned after 8-9 h each. Pair compute was 6,882 core-h against the plan's
  6,727 point estimate (+2.3%); wall was 19.5 h against the 9.8 h planned,
  entirely because of the tail: the top 1% of pairs carried 61.5% of the
  compute, not the 32% the 1,600-pair probe suggested. The worst pair that
  finished inside e8q_1 took 9 h 10 m.
- **e8q_2 (08-27) and e8q_4 (08-28)** were the first two `--keep-init`
  attempts at the nine, and each exposed a driver bug. e8q_2 died in setup: the
  one-shot's flat 600 s timeout could not parse the 31 MB init. e8q_4's merger
  was based on `e8qinitreference.at`, so finishing it would have written 49,881
  parameters over the 211,017; it was killed (its governor had also crashed at
  start with the tuple-comparison TypeError). All three bugs were fixed in
  `fpp_claude.py` before e8q_9.
- **e8q_3, e8q_5, e8q_6, e8q_7, e8q_8** were experiments (parameter sharing,
  the writexlam index, the sharing A/B) run from the reference into separate
  init files (`e8q_share_init.at`, `e8q_idx_init.at`, `e8q_ab_init.at`). They
  never touched `e8q_init.at`. Their hash counts of 211,014 and 211,009 differ
  from the answer only by the pairs each abandoned at its `--stop-at-remaining`
  floor.
- **e8q_9** ran the nine with the fixed driver: 16 workers,
  `--max-proc-gb 400`, `--stop-at-remaining 0`, the current `fpp_settings.at`
  (`logs/flags_changed.txt` lists four flags that differ from e8q_8:
  `face_verts_hash_flag`, `local_vertices_hash_flag`, `jag_KGB_frac` added,
  `khash_log_flag` removed).

### The nine hardest pairs

All nine had been killed after 8-9 hours in e8q_1. In e8q_9, against the
completed hash:

| pair | x | lambda | e8q_9 time | worker peak RSS |
|---|---|---|---|---|
| 5,925,863 | 67,078 | [3, 1, 1, 0, 1, 1, 2, 1] | 51.0 h | 165 GB |
| 5,924,456 | 66,989 | [4, 0, 0, 1, 1, 1, 1, 2] | 49.4 h | 160 GB |
| 5,925,409 | 67,045 | [3, -1, 1, 1, 1, 1, 2, 1] | 49.0 h | 164 GB |
| 5,925,416 | 67,046 | [4, 1, -1, 1, 1, 1, 2, 1] | 48.1 h | 162 GB |
| 5,924,461 | 66,989 | [4, 0, 0, 1, 1, 1, 2, 1] | 7.0 h | 121 GB |
| 5,922,202 | 66,885 | [-3, 0, 3, 1, 1, 1, 2, 2] | 6.9 h | 123 GB |
| 5,918,863 | 66,763 | [-3, 1, 4, 0, 1, 1, 1, 3] | 6.8 h | 113 GB |
| 5,913,723 | 66,601 | [-2, 1, 3, 1, 0, 1, 1, 3] | 83 s | 122 GB |
| 5,918,864 | 66,763 | [-3, 1, 4, 0, 1, 1, 2, 2] | 77 s | 114 GB |

The four two-day pairs each spent essentially all their time inside one or two
`is_unitary` tests (the 51 h pair: 30 tests, 18 of them on one K-type, 44 h in
two interrupted tests). Two pairs that had run over nine hours in e8q_1
finished in under 90 s here; whether that is the completed hash, the changed
settings or the fixed index was not determined. The idle workers' 113-123 GB
is startup memory (see below), so a hard pair adds roughly 40 GB on top.

### Cost of `--keep-init` at this init size

Restarting from the 31 MB init is expensive and should be budgeted:

| | from the reference (e8q_1) | `--keep-init` on the final init (e8q_9) |
|---|---|---|
| setup one-shot (group definition) | seconds | 12 min |
| worker launch | 166 s | 15 min; one worker's first attempt timed out (`HangError('timeout on: <prints>')`) and was retried |
| worker RSS before computing | 6 GB | 113-123 GB |
| peak worker RSS on a pair | 25 GB | 165 GB |

So the preset's 60 GB recycle trigger cannot be used with `--keep-init` here;
e8q_9 ran with `--max-proc-gb 400` and 16 workers.

### If E8_q ever has to be recomputed

Only if the code changes. The measured best configuration is the sharing arm
of the A/B (e8q_8): `--share --share-concurrency 32` with the preset's
`--stop-at-remaining 200`, which did the bulk in 6 h 50 m and 3,177 core-h
against e8q_7's 10 h 50 m and 6,526 core-h. Then a straggler pass:

```
./fpp_claude.py --preset e8q --keep-init --shards 16 --chunk 1 \
    --stop-at-remaining 0 --max-proc-gb 400 --max-total-gb 6000
```

with about two days and 165 GB per worker budgeted for it. Back up
`e8q_init.at` first.
