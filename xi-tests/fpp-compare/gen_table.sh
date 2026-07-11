#!/bin/bash
# Emit LaTeX table rows and a conclusion line from compare.sh output.
SP="$1"
DATA=$($SP/compare.sh "$SP" | tail -n +2)   # skip header
ALLMATCH=1
{
echo '\begin{center}'
echo '\begin{tabular}{lrrrrl}'
echo '\toprule'
echo 'Group & \#unit.\ (cov) & \#unit.\ (ref) & ms (cov) & ms (ref) & Set match \\'
echo '\midrule'
while IFS='|' read -r KEY GROUP NPC NPR MSC MSR SM NDIFF; do
  [ -z "$KEY" ] && continue
  # LaTeX-escape underscores in group name
  GL=$(echo "$GROUP" | sed 's/_/\\_/g')
  case "$SM" in
    YES) MATCH='\yes' ;;
    NO)  MATCH="\\no ($NDIFF)"; ALLMATCH=0 ;;
    *)   MATCH='\textit{missing}'; ALLMATCH=0 ;;
  esac
  echo "\\code{$GL} & $NPC & $NPR & $MSC & $MSR & $MATCH \\\\"
done <<< "$DATA"
echo '\bottomrule'
echo '\end{tabular}'
echo '\end{center}'
} > $SP/table.tex
if [ "$ALLMATCH" -eq 1 ]; then
  echo 'ALLMATCH' > $SP/verdict.txt
else
  echo 'MISMATCH' > $SP/verdict.txt
fi
echo "table.tex written; verdict=$(cat $SP/verdict.txt)"
