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
//      File name:     MuonE.cc
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 September 2024
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
#include "G4ParticleChangeForLoss.hh"
#include "G4MuIonisation.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4Region.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4TouchableHistory.hh"

#include "G4MuPairProductionModel.hh"
#include "G4RiGeMuPairProductionModel.hh"

#include "G4ParticleChangeForGamma.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include "G4UImanager.hh"

#include "G4ThreeVector.hh"

int main(int argc, char** argv)
{
  G4UImanager::GetUIpointer();
  G4cout << "========================================================" << G4endl;
  G4cout << "==========         Muon Test Starts           ==========" << G4endl;
  G4cout << "========================================================" << G4endl;
  //-----------------------------------------------------------------------------
  // ------- Initialisation

  G4int         verbose = 0;
  G4int         stat = 1e6;
  G4int         Z = 6;
  G4double      energy = 160.*GeV;
  G4Material*   mat = 0;

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
    energy *= GeV;
  }

  auto param = G4EmParameters::Instance();

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

  G4cout << "E(MeV)= " << energy/MeV << "  Z= " << Z
	 << "  " << mat->GetName() << G4endl;

  // Track
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  const G4ParticleDefinition* part = G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
  G4Positron::Positron();
  G4Electron::Electron();
  G4Gamma::Gamma();

  G4double particleMass = part->GetPDGMass();
  G4DynamicParticle dParticle(part, aDirection, energy - particleMass);
  //dParticle.SetPolarization(0.0,1.0,0.0);
  G4double eMass = CLHEP::electron_mass_c2;

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
  G4RiGeMuPairProductionModel mtmpp;
  G4MuPairProductionModel mpp;

  // Set particle change object
  G4ParticleChangeForLoss* fParticleChange = new G4ParticleChangeForLoss();
  mtmpp.SetParticleChange(fParticleChange, 0);
  mpp.SetParticleChange(fParticleChange, 0);

  G4DataVector cuts;
  cuts.push_back(keV);

  // Initilise models
  mtmpp.Initialise(part, cuts);
  mpp.Initialise(part, cuts);


  mtmpp.SetCurrentCouple(couple);
  mpp.SetCurrentCouple(couple);

  // ------- Histograms name
  Histo    histo;
  G4String hname = "muon_models";

  G4double pi = CLHEP::pi;

  G4int nbins = 100;
  G4double xmin =  0.0;
  G4double xmax =  10.0;
  if(99 == Z) {
    nbins = 650;
    xmin = -9.0;
  }
  G4double dx = (xmax - xmin)/G4double(nbins);
  G4double x0 = xmin - 0.5*dx;

  histo.Add1D("0","Energy of e+",nbins,0.2,160);
  histo.Add1D("1","Energy of e+",nbins,0,160);

  histo.Add1D("2","Energy of final mu",nbins,10.23,160);
  histo.Add1D("3","Energy of final mu",nbins,0,160);

  histo.Add1D("4","Energy of e-",nbins,0.2,160);
  histo.Add1D("5","Energy of e-",nbins,0,160);

  histo.Add1D("6","Polar angle of final mu",nbins,0.2e-3,4.8e-3);
  histo.Add1D("7","Polar angle of final mu",nbins,0,pi);

  histo.Add1D("8","Polar angle of e-",nbins,0,32e-3);
  histo.Add1D("9","Polar angle of e-",nbins,0,pi);

  histo.Add1D("10","Polar angle of e+",nbins,0,32e-3);
  histo.Add1D("11","Polar angle of e+",nbins,0,pi);

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

  std::vector<G4DynamicParticle*> vdp;

  G4double coste, costp, cost;
  G4double te, tp, t;
  G4double E, Ee , Ep;

  G4ThreeVector v1, v2, v3;

  G4Timer* timer = new G4Timer();
  timer->Start();

  for (G4int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    vdp.clear();
    mtmpp.SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    G4int N = (G4int) vdp.size();

    if (N == 3) {
      v1 = vdp[0]->GetMomentumDirection();
      v2 = vdp[1]->GetMomentumDirection();
      v3 = vdp[2]->GetMomentumDirection();

      E = (vdp[2]->GetKineticEnergy() + particleMass)/1e3;
      histo.Fill(2,E,1.0);
      Ee = (vdp[0]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(4,Ee,1.0);
      Ep = (vdp[1]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(0,Ep,1.0);
      cost = vdp[2]->GetMomentumDirection().z()/v3.mag();
      t = std::acos(cost);
      histo.Fill(6,t,1.0);
      coste = vdp[0]->GetMomentumDirection().z()/v1.mag();
      te = std::acos(coste);
      histo.Fill(8,te,1.0);
      costp = vdp[1]->GetMomentumDirection().z()/v2.mag();
      tp = std::acos(costp);
      histo.Fill(10,tp,1.0);

      delete vdp[0];
      delete vdp[1];
      delete vdp[2];
    }
    else{
      v1 = vdp[0]->GetMomentumDirection();
      v2 = vdp[1]->GetMomentumDirection();

      Ee = (vdp[0]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(4,Ee,1.0);
      Ep = (vdp[1]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(0,Ep,1.0);
      coste = vdp[0]->GetMomentumDirection().z()/v1.mag();
      te = std::acos(coste);
      histo.Fill(8,te,1.0);
      costp = vdp[1]->GetMomentumDirection().z()/v2.mag();
      tp = std::acos(costp);
      histo.Fill(10,tp,1.0);

      delete vdp[0];
      delete vdp[1];
    }
  }

  timer->Stop();
  G4cout << "MTMPP:  "  << *timer << G4endl;
  delete timer;
  // 2Gamma end
  G4cout << "  "   << G4endl;

  // Livermore begin
  timer = new G4Timer();
  timer->Start();

  for (G4int iter=0; iter<stat; ++iter) {
    fParticleChange->InitializeForPostStep(*track);
    vdp.clear();
    mpp.SampleSecondaries(&vdp,couple,&dParticle,0.0,energy);
    G4int N = (G4int) vdp.size();

    if (N == 3) {
      v1 = vdp[0]->GetMomentumDirection();
      v2 = vdp[1]->GetMomentumDirection();
      v3 = vdp[2]->GetMomentumDirection();

      E = (vdp[2]->GetKineticEnergy() + particleMass)/1e3;
      histo.Fill(3,E,1.0);
      Ee = (vdp[0]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(5,Ee,1.0);
      Ep = (vdp[1]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(1,Ep,1.0);
      cost = vdp[2]->GetMomentumDirection().z()/v3.mag();
      t = std::acos(cost);
      histo.Fill(7,t,1.0);
      coste = vdp[0]->GetMomentumDirection().z()/v1.mag();
      te = std::acos(coste);
      histo.Fill(9,te,1.0);
      costp = vdp[1]->GetMomentumDirection().z()/v2.mag();
      tp = std::acos(costp);
      histo.Fill(11,tp,1.0);
      
      delete vdp[0];
      delete vdp[1];
      delete vdp[2];
    }
    else{
      v1 = vdp[0]->GetMomentumDirection();
      v2 = vdp[1]->GetMomentumDirection();

      Ee = (vdp[0]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(5,Ee,1.0);
      Ep = (vdp[1]->GetKineticEnergy() + eMass)/1e3;
      histo.Fill(1,Ep,1.0);
      coste = vdp[0]->GetMomentumDirection().z()/v1.mag();
      te = std::acos(coste);
      histo.Fill(9,te,1.0);
      costp = vdp[1]->GetMomentumDirection().z()/v2.mag();
      tp = std::acos(costp);
      histo.Fill(11,tp,1.0);

      delete vdp[0];
      delete vdp[1];
    }
  }
  timer->Stop();
  G4cout << "MPP:  "  << *timer << G4endl;
  delete timer;
  // Livermore end
  G4cout << "  "   << G4endl;

  // -------- Committing the transaction with the tree
  //G4double norm = 1.0;
  //if(stat > 0.0) { norm /= G4double(stat); }
  //for(G4int i=10; i<20; ++i) { histo.ScaleH1(i, norm); }

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
