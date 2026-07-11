#!/bin/bash
# usage: apg_batch.sh <tree_dir> <outdir> <tag>
TREE="$1"; OUT="$2"; TAG="$3"
KEYS=(SL3 SL4 SL5 Sp4 Sp6 G2s Sp8 F4s E6s)
EXPR=("SL(3,R)" "SL(4,R)" "SL(5,R)" "Sp(4,R)" "Sp(6,R)" "G2_s" "Sp(8,R)" "F4_s" "E6_s")
for idx in "${!KEYS[@]}"; do
  K="${KEYS[$idx]}"; G="${EXPR[$idx]}"
  IN="$OUT/${K}.in"; LOG="$OUT/${K}.${TAG}.log"
  cat > "$IN" << EOF
<K_highest_weights.at
<groups.at
set G=$G
prints("RHO ", rho(G))
set A=all_parameters_gamma(G, rho(G))
set B=all_parameters_gamma(G, rho(G)/2)
prints("N_RHO ", #A)
prints("N_RHOHALF ", #B)
for p in A do prints("A ", p) od
for p in B do prints("B ", p) od
EOF
  cd "$TREE"
  timeout 1200 ./atlas --path=atlas-scripts < "$IN" > "$LOG" 2>&1
  RC=$?
  NR=$(grep -m1 "^N_RHO " "$LOG" | awk '{print $2}')
  NH=$(grep -m1 "^N_RHOHALF " "$LOG" | awk '{print $2}')
  ERR=$(grep -ciE "error|undefined|assert.*fail" "$LOG")
  echo "[$TAG] $K ($G): rc=$RC N_rho=${NR:-NA} N_rhohalf=${NH:-NA} errlines=$ERR"
done
echo "[$TAG] APG_DONE"
