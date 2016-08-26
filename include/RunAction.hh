#ifndef RunAction_hh
#define RunAction_hh 1
//
#include "G4UserRunAction.hh"
#include "G4UnitsTable.hh"
#include "DetectorConstruction.hh"

#include "TH1D.h"
#include <fstream>

//declare the DetectorConstruction class as we will define a pointer later
class DetectorConstruction;
class TTree;
class TFile;

//needed for using standard libraries
using namespace std;

//run action class, carries out tasks at the begin and end of each run
//the concept of a run incorporates a fixed geometry, fixed beam conditions, simulation of number of primaries
//begins with /run/beamOn command and finishes with tracking of last secondary to zero energy 

class RunAction : public G4UserRunAction {
//
public:
//run action class needs pointer ot the detector construction class in order to get details of the readout geometry
//accepts pointer to detector construction class

  RunAction(DetectorConstruction*);
  ~RunAction();

    //

private:
//an ofstream to access the output file
  ofstream outfile;
  TFile* RootOP;
  TH1D *EdepNoQuench; TH1D* EdepRear;

  TH1D *EdepQuenchPerEvt;
  TH1D *EdepQuenchPerInt;
  TTree* EdepTree;

  G4int IntType;
  G4double KinEnIn;
  G4double totalEdep;
  G4double TrackLenPrior;
  G4double TrackLenInVol;
  G4double VertX;
  G4double VertY;
  G4double VertR;
  G4double HFrontEn;
  G4double HRearEn;
  G4bool FrontHazPrimary;
  G4bool RearHazPrimary;
  G4int NumOptPhots;
  vector< vector< G4double > > PhotEnVec;
  vector< vector< double > > PhotParentIDVec;
  vector< string > PhotProcVec;
  vector< vector< double > > PhotWinEn;
  vector< vector< double > > WinUpDistance;
  vector< vector< double > > WinDownDistance;
  vector< vector< double > > PMTUpEn;
  vector< vector< double > > PMTUpParID;
  vector< vector< double > > PMTDownEn;
  vector< vector< double > > PMTDownParID;
  vector< vector< double > > PMTUpMeasEn;
  vector< vector< double > > PMTUpMeasParID;
  vector< vector< double > > PMTDownMeasEn;
  vector< vector< double > > PMTDownMeasParID;
  G4double PMTUpNumPhotons;
  G4double PMTDownNumPhotons;
  G4double PMTUpMeasPhotons;
  G4double PMTDownMeasPhotons;
  //G4double PMTUpPCPhotons;
  //G4double PMTDownPCPhotons;

//local pointer for detector construction class
    DetectorConstruction* myDC;

public:

//note argument of these methods is a pointer to a G4Run object
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void TallyEdepNoQuench(const G4double);
  //void TallyEdepQuenchPerEvt(const G4double);
  //void TallyEdepQuenchPerInt(const G4double);
  void TallyEdepRear(const G4double);
  //Argument list: Interaction type, Kinetic Energy in, Total Edep, Track length
  //prior to entering liquid, track length in liquid, vertex X, vertex Y, vertR,
  //Hfront Energy, Hrear Energy, Hfront Primary?, Hrear Primary?,
  //Tot # of opt photons, Opt Phot Init En Vector,
  //Opt Phot ParentID, Opt Phot Proc Vector,
  //Opt Phot En at Window,
  //COMMENTED OUT Distance from inferred primary track at upstream window,
  //COMMENTED OUT Distance from inferred primary track at downstream window,
  //Opt Phot En at upstream PMT,
  //ParentIDs at upstream PMT,
  //Opt Phot En at downstream PMT,
  //ParentIDs at downstream PMT,
  //# of photons received at upstream PMT,
  //# of photons received at downstream PMT,
  //Note: 'Measured' below means hitting Photocathode and after QE has been
  //accounted for.
  //Measured Opt Phot En at upstream PMT,
  //Measured ParentIDs at upstream PMT,
  //Measured Opt Phot En at downstream PMT,
  //Measured ParentIDs at downstream PMT,
  //Measued # Photons at upstream PMT,
  //Measured # photons at downstream PMT.
  void TallyEvtData(G4int, G4double, G4double, G4double, G4double,
		    G4double, G4double, G4double,
		    G4double, G4double, G4bool, G4bool,
		    G4int, vector< vector< G4double > >&,
		    vector< vector< double > >&, vector< string >&,
		    vector< vector< double > >&,
		    //vector< vector< double > >&,
		    //vector< vector< double > >&,
		    vector< vector< double > >&,
		    vector< vector< double > >&,
		    vector< vector< double > >&,
		    vector< vector< double > >&,
		    G4double,
		    G4double,
		    vector< vector< double > >&,
		    vector< vector< double > >&,
		    vector< vector< double > >&,
		    vector< vector< double > >&,
		    G4double,
		    G4double);

};

#endif
