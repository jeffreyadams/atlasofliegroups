#! /bin/sh
# As check-scripts.sh, but files that took at least $SLOW seconds to load on the last run,
# and are unchanged since then, are not checked. (Slowness is not size: the files over 1 MB
# load in under a second; warmE7s_to20000.at is 1.6 KB and takes 100 s.)
# Not absolutely safe: a slow file is skipped even if a file it includes has changed.
# check-scripts-fast.last-run lists the slow files; its modification time is the last-run
# date. Delete it to check everything. A slow file that fails is not listed, so it is
# checked again next time.

STAMP=$(pwd)/check-scripts-fast.last-run
NEW=$STAMP.new
SLOW=5

: > "$NEW"   # this run's start time: files changed during the run are checked next time
if [ -f "$STAMP" ]
then echo "Last run: $(date -r "$STAMP")"
else echo "No last run: checking everything"
fi

cd atlas-scripts >/dev/null
skipped=0
for i in $(git ls-files | grep '\.at$')
do
if [ -f "$STAMP" ] && grep -qxF "$i" "$STAMP" && ! [ "$i" -nt "$STAMP" ]
then skipped=$((skipped+1)); echo "$i" >> "$NEW"; continue
fi
#echo $i;   #uncomment this to get a running report on files checked
start=$(date +%s)
if ../atlas <$i >/dev/null 2>/dev/null
then if [ $(( $(date +%s) - start )) -ge $SLOW ]; then echo "$i" >> "$NEW"; fi
else echo Problem loading $i
fi
done
cd .. >/dev/null

mv "$NEW" "$STAMP"
echo "Done ($skipped slow unchanged files skipped)"
