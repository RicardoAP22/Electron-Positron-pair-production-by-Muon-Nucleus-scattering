#/bin/csh

cd $VFEM/photoeffect
mkdir -p $REFERENCE
cd $REFERENCE

mkdir -p test
mkdir -p plotXS
mkdir -p plot22

rm -rf *.gz

source $EMSRC/Photoeffect/runA.csh >! res.out

gzip *.out 
#
