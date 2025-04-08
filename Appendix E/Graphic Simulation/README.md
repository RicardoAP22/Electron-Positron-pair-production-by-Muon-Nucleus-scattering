This code is for running a graphic GEANT4 simulation. World its a 1m3 cube of "galatic" material with a thin plate of carbon in the middle of 0.75 cm of thickness.

The objetive its view the new G4RiGeMuPairProductionModel on the G4MuPairProduction process.

Beam its circular (1cm radius) and compose of muons with 160 GeV.

This code its an overall modification of example Basic B3a with the DetectorConstruction and PhysicsList extracted from extendend  electromagnetic example Em17 with a few tweaks on it.

To compile, make a folder named build and inside write "Cmake .." on your terminal. Then "make" . Finally run the executable MuonE.

For windows, insert CmakeList.txt on CMAKE and then compile with Visual Studio. Open the executable.

This test its for a GEANT4 GUI.

1. Recent Developments in Geant4, J. Allison et al., Nucl. Instrum. Meth. A 835 (2016) 186-225
2. Geant4 Developments and Applications, J. Allison et al., IEEE Trans. Nucl. Sci. 53 (2006) 270-278
3. Geant4 - A Simulation Toolkit, S. Agostinelli et al., Nucl. Instrum. Meth. A 506 (2003) 250-303