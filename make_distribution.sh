#!/bin/bash

here=$(pwd)
version=$1
distr=atlas_$version
distr_subdirs="messages sources rx-scripts"
distr_files="COPYRIGHT LICENSE README CHANGES"
distr_all_src="error interface io utilities structure gkmod test stand-alone"
dirs_with_cweb="interpreter stand-alone io"  # directories containing *.w files
rx_files="basic groups misc kl iterate_deform hermitian unitary translate \
 lietypes det sp4 my"


echo Building $distr

mkdir $distr
for sd in $distr_subdirs ; do mkdir $distr/$sd; done
for f in $distr_files ; do ln -s $here/$f $distr; done
ln -s $here/distr_INSTALL $distr/INSTALL
ln -s $here/distr_Makefile $distr/Makefile
ln -s $here/messages/* $distr/messages/
for f in $rx_files ; do ln -s $here/rx-scripts/$f.rx $distr/rx-scripts ; done
ln -s $here/rx-scripts/examples $distr/rx-scripts/
ln -s $here/doc/modules $distr/
ln -s $here/sources/interpreter/*.help $distr/
for sd in sources/*; do mkdir $distr/$sd; done
for sd in $distr_all_src
  do ln -s $here/sources/$sd/*.{h,cpp} $distr/sources/$sd/
  done
ln -s $here/sources/interface/input*readline.c $distr/sources/interface/
rm $distr/sources/stand-alone/'*'.h # a dangling link of this name was created!

for sd in $dirs_with_cweb
do pushd $distr/sources/$sd >/dev/null
  for f in $here/sources/$sd/*.w
  do ctanglex ++ -lbhp $f - $(basename $f .w)
  done
  popd >/dev/null
done
(cd $distr/sources/interpreter >/dev/null; bison -d --report=state \
  $here/sources/interpreter/parser.y)

cat sources/*/*.d >$distr/dependencies

tar cvz --file=$distr.tgz --dereference $distr
rm -r $distr
