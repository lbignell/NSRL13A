/*G4VSensitiveDetector is an abstract base class which represents a detector. The principal mandate of a sensitive detector is the construction of hit objects using information from steps along a particle track. The ProcessHits() method of G4VSensitiveDetector  performs this task using G4Step objects as input*/
#ifndef SensitiveDetector_h
#define SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Run.hh"
//#include "boost/variant.hpp"

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SensitiveDetector : public G4VSensitiveDetector
{

public:
    SensitiveDetector(G4String name);
    ~SensitiveDetector();

    /*This method is invoked at the beginning of each event. The argument of this method is an object of the G4HCofThisEvent class. Hits collections, where hits produced in this particular event are stored, can be associated to the G4HCofThisEvent object in this method. The hits collections associated with theG4HCofThisEvent  object during this method can be used for ``during the event processing'' digitization.G4Event has a G4HCofThisEvent class object, that is a container class of collections of hits. Hits collections are stored by their pointers, of the type of the base class.*/
    void Initialize(G4HCofThisEvent*HCE);
    /*This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. The firstargument of this method is a G4Step  object of the current step. The second argument is a G4TouchableHistory object for the ``Readout geometry'' describedin the next section. The second argument is NULL for the case ``Readout geometry'' is not assigned to this sensitive detector. In this method, one or moreG4VHit objects should be constructed if the current step is meaningful for your detector.*/
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    /*This method is invoked at the end of each event. The argument of this methd is the same object as the previous method. Hits collections occasionally created in your sensitive detector can be associated to the G4HCofThisEvent object.*/
    void EndOfEvent(G4HCofThisEvent*HCE);

  //This function will accept an energy, apply Birk's correction for Ionisation
  //quenching, and return the quenched energy.
  G4double ApplyBirksQuench(G4double Energy);
  G4double EvaluateBirks(G4double Energy);
  G4bool SampleQY(G4double PhotEn);

private:
  //G4double Dose;
  G4double runnum;

  G4bool MultipleIn;
  G4double KinEnIn;
  G4int EvtType;
  G4double TrackLenInVol;
  G4double TrackLenPrior;
  G4double VertexX;
  G4double VertexY;
  G4double VertexR;


  //Indexed by track ID.
  //1st bin is Vertex X position, 2nd is Vertex Y, 3rd vertex Z, 4th Energy.
  G4double SecondaryArray[1000][4];

  //Family tree.
  //Logged by track ID
  //0th argument is the track ID of the parent.
  //1st argument is the particle type (0=Gamma, 1=electron, 2=positron).
  G4int FamilyTree[1000][2];

  //Events Array.
  //The argument is the trackID of the ancestor that started the interaction.
  //The value stored is the energy deposited by that track and any electron
  //secondaries.
  G4double Events[1000];
  G4double DirectGammaInts[1000];

  //Different types of Edep Collections
  G4double EdepThisEventUnquenched;
  G4double EdepQuenchPerEvent;
  G4double EdepQuenchPerInteraction;

  //Primary Range. Pass this data on the primary beta tracklength for binning
  //into the primary beta range.
  G4double ThisPrimaryRange;
  
  //Number of Bremm interactions. Use this as a simple counter for the # of
  //Brem interactions occuring in the scintillant. Pass to RunAction at the end
  //of the event.
  G4int numInts;
  G4int numCorrs;
  G4int Counter;

  G4int TrackID;
  G4int ParentID;
  G4double TrackLength;
  G4double VertexEnergy;
  G4bool First;

  //Range fit parameters (9th order polynomial coeffs, lowest to highest order)
  G4double Const;
  G4double SpdOfLight;
  G4double ZonA;
  G4double MeanExEn;
  G4double ElecRestEnergy;
  G4double DensityEffectCorr;
  G4double ScintDensity;
  G4double kB;

  //typedef boost::variant<G4double, string, G4int> InitData;
  //vector< vector< InitData > > InitVec;

  G4int TotOptPhotons;
  vector< vector< double > > InitEnVec;
  //vector< string > InitProcVec;
  vector< vector< double > > ParentIDVec;
  vector< string > ProcVec;
  vector< double > dblvec;
  vector< int > intvec;
  vector< vector< double > > WinEnVec;
  vector< vector< double > > distVecUp;
  vector< vector< double > > distVecDown;

  G4bool HitWindow;
  G4int PhotonID;

  G4ThreeVector PrimEntryLocn;
  G4ThreeVector PrimMomDirn;

};

#endif
