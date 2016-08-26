//
#include "RunAction.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
//
#include "DetectorConstruction.hh"
#include "G4VSensitiveDetector.hh"
#include <math.h>
#include "TH1D.h"
#include "TFile.h"
#include <sstream>
#include "TTree.h"
#include "TROOT.h"

RunAction::RunAction(DetectorConstruction* DC){
  //
  
  //take the DetectorConstruction pointer given when this object is created (in main) and copy to local member
  myDC = DC;
  //
  
}

RunAction::~RunAction(){

}

void RunAction::BeginOfRunAction(const G4Run* aRun){

  //Arguments:Name of file,type of rewrite option(UPDATE means append),comments
  RootOP = new TFile("Edep.root","RECREATE","LS Sim Output");
  
  //arguments: Name, Title, number of bins, lower level range, upper
  //level range.
  //
  //All Energies are passed in MeV.
  //
  EdepNoQuench = new TH1D("Edep with no Quench","Energy Deposition in Scintillant per Event, no quench correction",210000,0.,210000.);
  //
  //EdepQuenchPerEvt = new TH1D("Edep with Quench per Evt","Energy Deposition in Scintillant per Event, quench correction applied per event",100000,0.,1000.);
  //
  //EdepQuenchPerInt = new TH1D("Edep with Quench per Int","Energy Deposition in Scintillant per Event, quench correction applied per Interaction",100000,0.,1000.);
  EdepRear = new TH1D("Edep with no Quench, rear vial","Energy Deposition in rear Scintillant per Event, no quench correction",210000,0.,210000.);

  EdepTree = new TTree("Results", "Particles causing Edep in Sensitive Vol");

  EdepTree->Branch("Type", &IntType);
  EdepTree->Branch("KinEn", &KinEnIn);
  EdepTree->Branch("Edep", &totalEdep);
  EdepTree->Branch("TrLenPrior", &TrackLenPrior);
  EdepTree->Branch("TrLenInVol", &TrackLenInVol);
  EdepTree->Branch("VertX", &VertX);
  EdepTree->Branch("VertY", &VertY);
  EdepTree->Branch("VertR", &VertR);
  EdepTree->Branch("FrontHodoEdep", &HFrontEn);
  EdepTree->Branch("RearHodoEdep", &HRearEn);
  EdepTree->Branch("PrimaryFrontHodo?", &FrontHazPrimary);
  EdepTree->Branch("PrimaryRearHodo?", &RearHazPrimary);
  EdepTree->Branch("NumberOfPhotons", &NumOptPhots);
  EdepTree->Branch("OptPhotonEnergies", &PhotEnVec);
  EdepTree->Branch("OptPhotonParentID", &PhotParentIDVec);
  EdepTree->Branch("OptPhotonProcess", &PhotProcVec);
  EdepTree->Branch("OptPhotonEnergyAtWindowBoth", &PhotWinEn);
  //EdepTree->Branch("UpWinPhotonDistFromPrim", &WinUpDistance);
  //EdepTree->Branch("DownWinPhotonDistFromPrim", &WinDownDistance);
  EdepTree->Branch("UpstreamPMTphotonEn", &PMTUpEn);
  EdepTree->Branch("UpstreamPMTparentIDs", &PMTUpParID);
  EdepTree->Branch("DownstreamPMTphotonEn", &PMTDownEn);
  EdepTree->Branch("DownstreamPMTparentIDs", &PMTDownParID);
  EdepTree->Branch("UpstreamPMTNumPhotons", &PMTUpNumPhotons);
  EdepTree->Branch("DownstreamPMTNumPhotons", &PMTDownNumPhotons);
  //EdepTree->Branch("UpstreamPMTPhotoCathodePhotons", &PMTUpPCPhotons);
  //EdepTree->Branch("DownstreamPMTPhotoCathodePhotons", &PMTDownPCPhotons);
  EdepTree->Branch("UpPMTMeasPhotonEn", &PMTUpMeasEn);
  EdepTree->Branch("UpPMTMeasParentIDs", &PMTUpMeasParID);
  EdepTree->Branch("DownPMTMeasPhotonEn", &PMTDownMeasEn);
  EdepTree->Branch("DownPMTMeasParentIDs", &PMTDownMeasParID);
  EdepTree->Branch("UpPMTMeasNumPhotons", &PMTUpMeasPhotons);
  EdepTree->Branch("DownPMTMeasNumPhotons", &PMTDownMeasPhotons);

}



void RunAction::TallyEdepNoQuench(const G4double thisEdep){

  //G4cout << "Inside TallyEdepNoQuench, Edep = " << thisEdep << G4endl;

  //G4double  thisEdepkeV = thisEdep*1000;

  //Fill ROOT histogram
  EdepNoQuench->Fill(thisEdep);

}

void RunAction::TallyEdepRear(const G4double thisEdep){
  
  //G4cout << "Inside TallyEdepRear, Edep = " << thisEdep << G4endl;

  G4double  thisEdepkeV = thisEdep*1000;

  //Fill ROOT histogram
  EdepRear->Fill(thisEdepkeV);

}


void RunAction::TallyEvtData(G4int type, G4double KinEn, G4double Edep,
			     G4double TrLenPrior, G4double TrLenInVol,
			     G4double VxX, G4double VxY, G4double VxR,
			     G4double HFrEn, G4double HReEn,
			     G4bool FrHazPrim, G4bool ReHazPrim, 
			     G4int NumPhotons,
			     vector< vector< double > > &PhotEn,
			     vector< vector< double > >& PhotParentID,
			     vector< string >& PhotProc,
			     vector< vector< double > >& WinEn,
			     //vector< vector< double > >& DistUp,
			     //vector< vector< double > >& DistDown,
			     vector< vector< double > >& UPMTEn,
			     vector< vector< double > >& UParID,
			     vector< vector< double > >& DPMTEn,
			     vector< vector< double > >& DParID,
			     G4double UpPhots, G4double DownPhots,
			     // G4double UpPCPhots, G4double DownPCPhots,
			     vector< vector< double > >& MeasUPMTEn,
			     vector< vector< double > >& MeasUParID,
			     vector< vector< double > >& MeasDPMTEn,
			     vector< vector< double > >& MeasDParID,
			     G4double MeasUpPhots, G4double MeasDownPhots){
  //Set the branch values and fill the tree.
  IntType = type;
  KinEnIn = KinEn;
  totalEdep = Edep;
  TrackLenPrior = TrLenPrior;
  TrackLenInVol = TrLenInVol;
  VertX = VxX;
  VertY = VxY;
  VertR = VxR;
  HFrontEn = HFrEn;
  HRearEn = HReEn;
  FrontHazPrimary = FrHazPrim;
  RearHazPrimary = ReHazPrim;
  NumOptPhots = NumPhotons;
  //G4cout << "Number of Photons = " << NumPhotons << G4endl;
  PhotEnVec = PhotEn;
  PhotParentIDVec = PhotParentID;
  PhotProcVec = PhotProc;
  PhotWinEn = WinEn;
  //WinUpDistance = DistUp;
  //WinDownDistance = DistDown;
  PMTUpEn = UPMTEn;
  PMTUpParID = UParID;
  PMTDownEn = DPMTEn;
  PMTDownParID = DParID;
  PMTUpNumPhotons = UpPhots;
  PMTDownNumPhotons = DownPhots;
  //PMTUpPCPhotons = UpPCPhots;
  //PMTDownPCPhotons = DownPCPhots;
  PMTUpMeasEn = MeasUPMTEn;
  PMTUpMeasParID = MeasUParID;
  PMTDownMeasEn = MeasDPMTEn;
  PMTDownMeasParID = MeasDParID;
  PMTUpMeasPhotons = MeasUpPhots;
  PMTDownMeasPhotons = MeasDownPhots;

  EdepTree->Fill();

}

//task to be carried out at the end of the run
void RunAction::EndOfRunAction(const G4Run* aRun){
  //get the number of primary particles being simulated for this run
  G4double NumberOfEvents = aRun->GetNumberOfEventToBeProcessed();
 
  // Name the histograms
  G4String SpectName;
  G4String Beginning = "EdepNoQuench";
  //stringstream s;
  //G4double RunNumber = aRun->GetRunID();
  //s << RunNumber;
  //G4String RunNum;
  //s >> RunNum;
  //SpectName = Beginning+RunNum;
  
  //EdepNoQuench->SetName(Beginning);
  //EdepNoQuench->Scale(1/NumberOfEvents);
  //EdepNoQuench->Write();

  //Beginning = "EdepNoQuenchRear";
  //SpectName = Beginning+RunNum;
  //EdepRear->SetName(Beginning);
  //EdepRear->Scale(1/NumberOfEvents);
  //EdepRear->Write();
 
  //Beginning = "EdepQuenchPerEvtFromRunNum";
  //SpectName = Beginning+RunNum;
  //EdepQuenchPerEvt->SetName(SpectName);
  //EdepQuenchPerEvt->Scale(1/NumberOfEvents);
  //EdepQuenchPerEvt->Write();

  // Beginning = "EdepQuenchPerIntFromRunNum";
  //SpectName = Beginning+RunNum;
  //EdepQuenchPerInt->SetName(SpectName);
  //EdepQuenchPerInt->Scale(1/NumberOfEvents);
  //EdepQuenchPerInt->Write();

  RootOP->Write();
  RootOP->Close();
  

}
