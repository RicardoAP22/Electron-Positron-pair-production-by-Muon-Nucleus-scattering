#!/bin/csh -f

cd $VFEM/photoeffect/plotXS
mkdir -p $1
cd $1
rm -rf *

set iL = 0;
while ($iL < 98)
    @ iL++;
    if(-r $VFEM/photoeffect/${1}/plotXS/ApicXS_$iL.png) then
	ln -s ../../${1}/plotXS/ApicXS_$iL.png
    endif
end

cd ../../plotXSr
mkdir -p $1
cd $1
rm -rf *

set iL = 0;
while ($iL < 98)
    @ iL++;
    if(-r $VFEM/photoeffect/${1}/plotXS/ApicXSr_$iL.gif) then
	ln -s ../../${1}/plotXS/ApicXSr_$iL.gif
    endif
end

cd ../../plotEnergy
mkdir -p $1
cd $1
rm -rf *

set iL = 0;
while ($iL < 98)
    @ iL++;
    if(-r $VFEM/photoeffect/${1}/plot22/ApicE_$iL.gif) then
	ln -s ../../${1}/plot22/ApicE_$iL.gif
    endif
end

cd ../../plotAngle
mkdir -p $1
cd $1
rm -rf *

set iL = 0;
while ($iL < 98)
    @ iL++;
    if(-r $VFEM/photoeffect/${1}/plot22/ApicA_$iL.gif) then
	ln -s ../../${1}/plot22/ApicA_$iL.gif
    endif
end

cd $VFEM/photoeffect

#
