#! /bin/sh

cd atlas-scripts >/dev/null
for i in $(git ls-files | grep '\.at$')
do if ! ../atlas <$i >/dev/null 2>/dev/null; then echo Problem loading $i; fi
done
echo Done
