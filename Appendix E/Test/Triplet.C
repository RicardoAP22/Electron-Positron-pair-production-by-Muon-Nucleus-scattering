//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     photo.cc
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 September 2011 from test30.cc
//
//      Modifications:
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4Region.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4TouchableHistory.hh"

#include "G4LivermorePhotoElectricModel.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PEEffectFluoModel.hh"

#include "G4ParticleChangeForGamma.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include "G4UImanager.hh"

int main(int argc, char** argv)
{
  G4UImanager::GetUIpointer();
  G4cout << "========================================================" << G4endl;
  G4cout << "======        Photoeffect Test Starts           ========" << G4endl;
  G4cout << "========================================================" << G4endl;
  //-----------------------------------------------------------------------------
  // ------- Initialisation 

  G4int     verbose  = 0;
  G4int     stat     = 100000;
  G4int     Z = 13;
  G4double  energy   = 1.*MeV;
  G4Material* mat = 0;

  // convert string to Z
  G4String sz = "";
  if(argc >= 2) { 
    sz = argv[1];
    std::istringstream is(sz);
    is >> Z;
  }

  // convert string to energy
  G4String se = "";
  if(argc >= 3) { 
    se = argv[2];
    std::istringstream is(se);
    is >> energy;
    energy *= MeV;
  }

  // convert string to statistics
  G4String ss = "";
  if(argc >= 4) { 
    ss = argv[3];
    std::istringstream is(ss);
    is >> stat;
  }

  // -------- Material
  G4NistManager::Instance()->SetVerbose(2);
  G4String mname = "G4_";
  if(Z < 99) {
    const G4Element* elm = G4NistManager::Instance()->FindOrBuildElement(Z);
    mname += elm->GetName();
  } else {
    mname += "WATER";
  }
  mat = G4NistManager::Instance()->FindOrBuildMaterial(mname);
  G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

  G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(mat);
  couple->SetIndex(0);

  G4double coeff = mat->GetDensity()*cm2/g;

  if(!mat) { exit(1); }
  /*
  G4cout << "========================================================" << G4endl;
  G4cout << "======        Photoeffect Test Starts           ========" << G4endl;
  G4cout << "========================================================" << G4endl;
  */
  G4cout << "E(MeV)= " << energy/MeV << "  Z= " << Z 
	 << "  " << mat->GetName() << G4endl;

  // Track
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  const G4ParticleDefinition* part = G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  
  G4DynamicParticle dParticle(part,aDirection,energy);
  dParticle.SetPolarization(0.0,1.0,0.0);

  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  //--------- Geometry definition

  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,mat,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);

  G4ProductionCutsTable* theCoupleTable = 
    G4ProductionCutsTable::GetProductionCutsTable();
  theCoupleTable->UpdateCoupleTable(pFrame);

  G4Region* reg = new G4Region("DefaultRegionForTheWorld");
  reg->AddRootLogicalVolume(lFrame);

  // -------------------------------------------------------------------
  // -------- Start run processing

  G4StateManager* g4State=G4StateManager::GetStateManager();
  if (! g4State->SetNewState(G4State_Init)) {
    G4cout << "error changing G4state"<< G4endl;;   
  }

  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  de->SetFluo(true);
  de->InitialiseAtomicDeexcitation();
  de->SetFluo(false);

  // Instantiate models
  G4PenelopePhotoElectricModel penpe;
  G4LivermorePhotoElectricModel liv;
  G4LivermorePhotoElectricModel livp;
  G4PEEffectFluoModel pee;

  livp.SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized());

  // Set particle change object
  G4ParticleChangeForGamma* fParticleChange = new G4ParticleChangeForGamma();
  liv.SetParticleChange(fParticleChange, 0);
  livp.SetParticleChange(fParticleChange, 0);
  penpe.SetParticleChange(fParticleChange, 0);
  pee.SetParticleChange(fParticleChange, 0);
   
  G4DataVector cuts;
  cuts.push_back(keV);

  // Initilise models
  penpe.Initialise(part, cuts);
  liv.Initialise(part, cuts);
  livp.Initialise(part, cuts);
  pee.Initialise(part, cuts);

  penpe.SetCurrentCouple(couple);
  liv.SetCurrentCouple(couple);
  livp.SetCurrentCouple(couple);
  pee.SetCurrentCouple(couple);

  // ------- Histograms name 
  Histo    histo;
  G4String hname = "test/photo_" + sz + "_" + se;
  
  G4int nbins = 500;
  G4double xmin = -6.0;
  G4double xmax =  4.0;
  if(99 == Z) {
    nbins = 650;
    xmin = -9.0;
  }
  G4double dx = (xmax - xmin)/G4double(nbins);
  G4double x0 = xmin - 0.5*dx;

  histo.Add1D("0","PhotoEffect cross section Penelope",nbins,xmin,xmax);
  histo.Add1D("1","PhotoEffect cross section Livermore",nbins,xmin,xmax);
  histo.Add1D("2","PhotoEffect cross section Polarized Model",nbins,xmin,xmax);
  histo.Add1D("3","PhotoEffect cross section Standard Model",nbins,xmin,xmax);
  histo.Add1D("4","empty",nbins,xmin,xmax);

  histo.Add1D("5","PhotoEffect cross section ratio Penelope",nbins,xmin,xmax);
  histo.Add1D("6","PhotoEffect cross section ratio Livermore",nbins,xmin,xmax);
  histo.Add1D("7","PhotoEffect cross section ratio Polarized Model",nbins,xmin,xmax);
  histo.Add1D("8","PhotoEffect cross section ratio Standard Model",nbins,xmin,xmax);
  histo.Add1D("9","empty",nbins,xmin,xmax);

  histo.Add1D("10","Electron Energy Spectra Penelope respective to beam energy",nbins,-3,3);
  histo.Add1D("11","Electron Energy Spectra Livermore respective to beam energy",nbins,-3,3);
  histo.Add1D("12","Electron Energy Spectra Polarized Model respective to beam energy",nbins,-3,3);
  histo.Add1D("13","Electron Energy Spectra Standard Model respective to beam energy",nbins,-3,3);
  histo.Add1D("14","empty",nbins,-3,3);

  histo.Add1D("15","Electron cosTheta Penelope",nbins,-1,1);
  histo.Add1D("16","Electron cosTheta Livermore",nbins,-1,1);
  histo.Add1D("17","Electron cosTheta Polarized Model",nbins,-1,1);
  histo.Add1D("18","Electron cosTheta Standard Model",nbins,-1,1);


  histo.SetFileName(hname);
  histo.Book();
  G4cout << "Histograms are booked output file <" << hname << "> "
	 << G4endl;
  
  // -------- Track
  G4Track* track = new G4Track(&dParticle,0.0,aPosition);
  G4TouchableHandle fpTouchable(new G4TouchableHistory());
  track->SetTouchableHandle(fpTouchable);

  // -------- Step
  if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
    { G4cout << "G4StateManager PROBLEM! " << G4endl; }

  // -------- Event loop
  for (G4int iter=0; iter<nbins; ++iter) {

    x0 += dx;
    G4double e = std::pow(10., x0)*MeV;
    dParticle.SetKineticEnergy(e);

    G4double sig1 = penpe.CrossSectionPerVolume(mat,part,e)/coeff;
    G4double sig2 = liv.CrossSectionPerVolume(mat,part,e)/coeff;
    G4double sig3 = livp.CrossSectionPerVolume(mat,part,e)/coeff;
    G4double sig4 = pee.CrossSectionPerVolume(mat,part,e,0.,e)/coeff;

    if(verbose>=2) { 
      G4cout << "### " << iter << "-th event start " 
	     << " E(MeV)= " << e/MeV
	     << " pen(g/cm2)= " << sig1 
	     << " liv(g/cm2)= " << sig2 
             << " Polarized Model (g/cm2)= " << sig3
	     << " Standard Model (g/cm2)= " << sig4
	     << G4endl;

    }
    
    histo.Fill(0,x0,sig1);
    histo.Fill(1,x0,sig2);
    histo.Fill(2,x0,sig3);
    histo.Fill(3,x0,sig4);

    if(sig1 > 0.0) {
      sig2 /= sig1;
      if(sig2 > 100.) { sig2 = 0.0; }
      sig3 /= sig1;
      if(sig3 > 100.) { sig3 = 0.0; }
      sig4 /= sig1;
      if(sig4 > 100.) { sig4 = 0.0; }
    } else {
      sig2 = 0.0;
      sig3 = 0.0;
      sig4 = 0.0;
    }
    histo.Fill(5,x0,1.0);
    histo.Fill(6,x0,sig2);
    histo.Fill(7,x0,sig3);
    histo.Fill(8,x0,sig4);
  }

  std::vector<G4DynamicParticle*> vdp;
  dParticle.SetKineticEnergy(energy);
  G4double cost, e1;
  
  G4Timer* timer = new G4Timer();
  timer->Start();
  
  for (G4int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    penpe.SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    G4double E1 = vdp[0]->GetKineticEnergy();
    e1 = log10((energy - E1)/keV);
    histo.Fill(10,e1,1.0); 
    cost = vdp[0]->GetMomentumDirection().z();
    histo.Fill(15,cost,1.0);  
    delete vdp[0]; 
    vdp.clear();
  }
  
  timer->Stop();
  G4cout << "Penelope:  "  << *timer << G4endl;
  delete timer;
  // Penelope end
  G4cout << "  "   << G4endl;

  // Livermore begin
  timer = new G4Timer();
  timer->Start();

  for (G4int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    liv.SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    G4double E1 = vdp[0]->GetKineticEnergy();
    e1 = log10((energy - E1)/keV);
    histo.Fill(11,e1,1.0); 
    cost = vdp[0]->GetMomentumDirection().z();
    histo.Fill(16,cost,1.0);  
    delete vdp[0]; 
    vdp.clear();
  }
  timer->Stop();
  G4cout << "Livermore:  "  << *timer << G4endl;
  delete timer;
  // Livermore end
  G4cout << "  "   << G4endl;

  // Polarized begin
  timer = new G4Timer();
  timer->Start();

   for (G4int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    livp.SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    G4double E1 = vdp[0]->GetKineticEnergy();
    e1 = log10((energy - E1)/keV);
    histo.Fill(12,e1,1.0); 
    cost = vdp[0]->GetMomentumDirection().z();
    histo.Fill(17,cost,1.0);  
    delete vdp[0]; 
    vdp.clear();
  }
  timer->Stop();
  G4cout << "LivermorePolarized:  "  << *timer << G4endl;
  delete timer;
  // Polarized end
  G4cout << "  "   << G4endl;

  // Stand Model begin
  timer = new G4Timer();
  timer->Start();

  for (G4int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    pee.SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    G4double E1 = vdp[0]->GetKineticEnergy();
    e1 = log10((energy - E1)/keV);
    //    e1 = vdp[0]->GetKineticEnergy()/energy;
    histo.Fill(13,e1,1.0); 
    cost = vdp[0]->GetMomentumDirection().z();
    histo.Fill(18,cost,1.0);  
    delete vdp[0]; 
    vdp.clear();
  }
  timer->Stop();
  G4cout << "Stand Model:  "  << *timer << G4endl;
  delete timer;
  // Stand model end
  G4cout << "  "   << G4endl;
  
  // -------- Committing the transaction with the tree
  G4double norm = 1.0;
  if(stat > 0.0) { norm /= G4double(stat); }
  for(G4int i=10; i<20; ++i) { histo.ScaleH1(i, norm); }

  if(verbose > 0) { G4cout << "###### Save histograms" << G4endl; }
  histo.Save();

  if(verbose > 0) {
    G4cout << "###### End of run # " << G4endl;
  }

  delete pFrame;
  delete lFrame;
  delete sFrame;
  partTable->DeleteAllParticles();

  G4cout << "###### End of test #####" << G4endl;
}
