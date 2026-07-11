# all_parameters_gamma at rho and rho/2 — covering vs stock master

For each of 9 groups (SL3/4/5, Sp4/6/8, G2_s, F4_s, E6_s), enumerate
all_parameters_gamma(G, rho(G)) and all_parameters_gamma(G, rho(G)/2) on the
covering branch and on pristine ~/atlasSoftware/master, and compare.

- apg_batch.sh <tree> <outdir> <tag> : run all 9 groups, dump params at rho (A)
  and rho/2 (B) plus counts.
- apg_compare.sh <outdir> : per group, covering-vs-master set diff at rho and rho/2.
- apg-compare-results.txt : pipe-delimited summary
  (KEY|GROUP|Nrho_cov|Nrho_mas|match_rho|Nhalf_cov|Nhalf_mas|match_half).
- SL3.example.log : full parameter listing for SL(3,R) (worked example).

Result: covering == master EXACTLY (parameter-set level) for all 9 groups at
BOTH infinitesimal characters. rho/2 (half-integral: <rho/2,alpha^v>=1/2 for
all simple alpha) gives ~1/4 to ~1/2 as many params as rho. Writeup:
coverNotes/apg-comparison.tex.
