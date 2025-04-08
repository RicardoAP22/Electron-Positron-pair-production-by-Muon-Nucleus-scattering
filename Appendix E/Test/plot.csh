#/bin/csh

cd $VFEM/photoeffect/$REFERENCE

gunzip -fq *.gz

root -b -q $EMSRC/Photoeffect/plotXS.C >! plot.out
root -b -q $EMSRC/Photoeffect/plotXS_ratio.C >>& plot.out
root -b -q $EMSRC/Photoeffect/plot22A.C >>& plot.out
root -b -q $EMSRC/Photoeffect/plot22E.C >>& plot.out

gzip -fq *.out 

source $EMSRC/Photoeffect/summary.csh $REFERENCE
cd $VFEM/photoeffect/$REFERENCE

#
