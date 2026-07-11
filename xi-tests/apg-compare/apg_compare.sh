#!/bin/bash
SP="$1"
KEYS=(SL3 SL4 SL5 Sp4 Sp6 G2s Sp8 F4s E6s)
EXPR=("SL(3,R)" "SL(4,R)" "SL(5,R)" "Sp(4,R)" "Sp(6,R)" "G2_s" "Sp(8,R)" "F4_s" "E6_s")
echo "KEY|GROUP|NRHO_C|NRHO_M|MATCH_RHO|NHALF_C|NHALF_M|MATCH_HALF"
for idx in "${!KEYS[@]}"; do
  K="${KEYS[$idx]}"; G="${EXPR[$idx]}"
  CL="$SP/${K}.COVER.log"; ML="$SP/${K}.MASTER.log"
  NRC=$(grep -m1 "^N_RHO " "$CL" 2>/dev/null | awk '{print $2}')
  NRM=$(grep -m1 "^N_RHO " "$ML" 2>/dev/null | awk '{print $2}')
  NHC=$(grep -m1 "^N_RHOHALF " "$CL" 2>/dev/null | awk '{print $2}')
  NHM=$(grep -m1 "^N_RHOHALF " "$ML" 2>/dev/null | awk '{print $2}')
  # sorted param sets
  grep "^A " "$CL" 2>/dev/null | sort > "$SP/${K}.A.C"; grep "^A " "$ML" 2>/dev/null | sort > "$SP/${K}.A.M"
  grep "^B " "$CL" 2>/dev/null | sort > "$SP/${K}.B.C"; grep "^B " "$ML" 2>/dev/null | sort > "$SP/${K}.B.M"
  if [ -s "$SP/${K}.A.C" ] && [ -s "$SP/${K}.A.M" ]; then
    DA=$(comm -3 "$SP/${K}.A.C" "$SP/${K}.A.M" | grep -c .); [ "$DA" -eq 0 ] && MA=YES || MA="NO($DA)"
  else MA="MISSING"; fi
  if [ -s "$SP/${K}.B.C" ] && [ -s "$SP/${K}.B.M" ]; then
    DB=$(comm -3 "$SP/${K}.B.C" "$SP/${K}.B.M" | grep -c .); [ "$DB" -eq 0 ] && MB=YES || MB="NO($DB)"
  else MB="MISSING"; fi
  echo "$K|$G|${NRC:-NA}|${NRM:-NA}|$MA|${NHC:-NA}|${NHM:-NA}|$MB"
done
