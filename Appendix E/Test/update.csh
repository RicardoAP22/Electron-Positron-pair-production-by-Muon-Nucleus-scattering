#!/bin/csh -f

set ir = 4;
set listR = (geant4-09-05-ref-09 geant4-10-00-ref-00 geant4-10-01-ref-00 geant4-10-02-ref-00);

set il = 4;
set listL = (plotXS  plotXSr plotAngle plotEnergy);
set listB = (plotXS plotXS plot22  plot22);
set listA = (ApicXS_ ApicXSr_ ApicA_ ApicE_);

set in = 98;

echo "step1"

#*** Main loop ***
set iL = 0;
while ($iL < $il)
    @ iL++;
    mkdir -p $listL[$iL]
    set iR = 0;
    echo "step3 " $iL 
    while ($iR < $ir)
      @ iR++;
      echo "step4 " $iR 
      mkdir -p $listL[$iL]/$listR[$iR]
      cd $listL[$iL]/$listR[$iR]
      set iN = 0; 
      while ($iN < $in)
        @ iN++;
        set x = ../../$listR[$iR]/$listB[$iL]/$listA[$iL]$iN;
        if(-r $x.gif) then
          ln -s $x.gif
        endif
        if(-r $x.png) then
          ln -s $x.png
        endif
      end
      cd ../../
    end
end

