#!/bin/bash
# Compare COVER vs REF logs per group; emit a data table + per-group diff files.
SP="$1"
KEYS=(SL3 SL4 SL5 Sp4 Sp6 G2s Sp8 F4s E6s)
EXPR=("SL(3,R)" "SL(4,R)" "SL(5,R)" "Sp(4,R)" "Sp(6,R)" "G2_s" "Sp(8,R)" "F4_s" "E6_s")
echo "KEY|GROUP|NP_COVER|NP_REF|MS_COVER|MS_REF|SETMATCH|NDIFF"
for idx in "${!KEYS[@]}"; do
  K="${KEYS[$idx]}"; G="${EXPR[$idx]}"
  CL="$SP/${K}.COVER.log"; RL="$SP/${K}.REF.log"
  NPC=$(grep -m1 "^NPARAMS " "$CL" 2>/dev/null | awk '{print $2}')
  NPR=$(grep -m1 "^NPARAMS " "$RL" 2>/dev/null | awk '{print $2}')
  MSC=$(grep -m1 "^ELAPSED_MS " "$CL" 2>/dev/null | awk '{print $2}')
  MSR=$(grep -m1 "^ELAPSED_MS " "$RL" 2>/dev/null | awk '{print $2}')
  # sorted param sets (strip "P " prefix)
  grep "^P " "$CL" 2>/dev/null | sed 's/^P //' | sort > "$SP/${K}.COVER.set"
  grep "^P " "$RL" 2>/dev/null | sed 's/^P //' | sort > "$SP/${K}.REF.set"
  if [ -s "$SP/${K}.COVER.set" ] && [ -s "$SP/${K}.REF.set" ]; then
    NDIFF=$(comm -3 "$SP/${K}.COVER.set" "$SP/${K}.REF.set" | grep -c .)
    if [ "$NDIFF" -eq 0 ]; then SM="YES"; else SM="NO"; fi
    comm -3 "$SP/${K}.COVER.set" "$SP/${K}.REF.set" > "$SP/${K}.diff"
  else
    NDIFF="NA"; SM="MISSING"
  fi
  echo "$K|$G|${NPC:-NA}|${NPR:-NA}|${MSC:-NA}|${MSR:-NA}|$SM|$NDIFF"
done
