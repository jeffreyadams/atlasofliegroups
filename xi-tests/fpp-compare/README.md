# FPP unitary-dual cross-branch comparison

Compares `FPP_unitary_hash_bottom_layer(G)` on this (covering) branch vs the
reference branch `~/atlasSoftware/to_ht_branch_alt` (origin of the FPP code).

- `batch.sh <tree> <outdir> <tag>` — runs all 9 groups on a tree, dumping each
  group's sorted unitary-param set + count + elapsed_ms.
- `compare.sh <outdir>` — diffs COVER vs REF sorted param sets per group.
- `gen_table.sh <outdir>` — emits the LaTeX table.
- `cover_batch.out` / `ref_batch.out` — the parallel-pass run logs.
- `reclean.out` — uncontended (sequential) re-measurement of F4_s and Sp(8,R).
- `compare-results.txt` — pipe-delimited summary (KEY|GROUP|NP_cov|NP_ref|ms_cov|ms_ref|setmatch|ndiff).

Result: all 9 groups (SL3/4/5, Sp4/6/8, G2_s, F4_s, E6_s) produce IDENTICAL
unitary-param sets on both branches. Timing comparable except F4_s (~2.9x
slower on covering; traced to master vs to_ht_branch_alt C++ cores, not the
port or xi=0 cover code). Full writeup: coverNotes/fpp-comparison.tex.
