#! /bin/sh
#.ax files are not checked, use this for large files, for example bigblocke8.ax

cd atlas-scripts >/dev/null
for i in $(git ls-files | grep '\.at$')
do
#echo $i;   #uncomment this to get a running report on files checked
if ! ../atlas <$i >/dev/null 2>/dev/null; then echo Problem loading $i; fi
done
echo Done
