#!/bin/bash
# usage: batch.sh <tree_dir> <outdir> <tag>
TREE="$1"; OUT="$2"; TAG="$3"
KEYS=(SL3 SL4 SL5 Sp4 Sp6 G2s Sp8 F4s E6s)
EXPR=("SL(3,R)" "SL(4,R)" "SL(5,R)" "Sp(4,R)" "Sp(6,R)" "G2_s" "Sp(8,R)" "F4_s" "E6_s")
for idx in "${!KEYS[@]}"; do
  K="${KEYS[$idx]}"; G="${EXPR[$idx]}"
  IN="$OUT/${K}.in"
  LOG="$OUT/${K}.${TAG}.log"
  cat > "$IN" << EOF
<FPP_globalDirac.at
<groups.at
set G=$G
set t0=elapsed_ms()
FPP_unitary_hash_bottom_layer(G)
set t1=elapsed_ms()
set HH=big_unitary_hash.uhash(G)
prints("NPARAMS ", HH.size())
prints("ELAPSED_MS ", t1-t0)
for i:HH.size() do prints("P ", HH.index(i)) od
EOF
  cd "$TREE"
  START=$(date +%s)
  timeout 5400 ./atlas --path=atlas-scripts < "$IN" > "$LOG" 2>&1
  RC=$?
  WALL=$(( $(date +%s) - START ))
  NP=$(grep -m1 "^NPARAMS " "$LOG" | awk '{print $2}')
  EM=$(grep -m1 "^ELAPSED_MS " "$LOG" | awk '{print $2}')
  echo "[$TAG] $K ($G): rc=$RC nparams=${NP:-NONE} elapsed_ms=${EM:-NA} wall_s=$WALL"
done
echo "[$TAG] BATCH_DONE"
