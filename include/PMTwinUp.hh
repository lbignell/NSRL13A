/*G4VSensitiveDetector is an abstract base class which represents a detector. The principal mandate of a sensitive detector is the construction of hit objects using information from steps along a particle track. The ProcessHits() method of G4VSensitiveDetector  performs this task using G4Step objects as input*/
#ifndef PMTwinUp_h
#define PMTwinUp_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Run.hh"

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PMTwinUp : public G4VSensitiveDetector
{

public:
    PMTwinUp(G4String name);
    ~PMTwinUp();

    /*This method is invoked at the beginning of each event. The argument of this method is an object of the G4HCofThisEvent class. Hits collections, where hits produced in this particular event are stored, can be associated to the G4HCofThisEvent object in this method. The hits collections associated with theG4HCofThisEvent  object during this method can be used for ``during the event processing'' digitization.G4Event has a G4HCofThisEvent class object, that is a container class of collections of hits. Hits collections are stored by their pointers, of the type of the base class.*/
    void Initialize(G4HCofThisEvent*HCE);
    /*This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. The firstargument of this method is a G4Step  object of the current step. The second argument is a G4TouchableHistory object for the ``Readout geometry'' describedin the next section. The second argument is NULL for the case ``Readout geometry'' is not assigned to this sensitive detector. In this method, one or moreG4VHit objects should be constructed if the current step is meaningful for your detector.*/
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    /*This method is invoked at the end of each event. The argument of this methd is the same object as the previous method. Hits collections occasionally created in your sensitive detector can be associated to the G4HCofThisEvent object.*/
    void EndOfEvent(G4HCofThisEvent*HCE);

  vector< vector< double > >& GetEnSpec();
  vector< vector< double > >& GetParentIDs();
  vector< vector< double > >& GetMeasEnSpec();
  vector< vector< double > >& GetMeasParentIDs();
  G4double GetTotalPhotons();
  G4double GetMeasPhotons();
  G4double GetPhotsOnPhotoCathode();

  G4bool PMTHazFired(G4double);//Argument is Opt photon wavelength in MeV.

private:
  //G4double Dose;
  G4double runnum;

  G4bool HazPrimary;

  //1st column cerenkov, 2nd column scintillation.
  vector< vector< double > > PhotEn;
  vector< vector< double > > ParentID;
  vector< double > dblvec;
  vector< vector< double > > PMTFirePhotEnIn_nm;
  vector< vector< double > > PMTFireParentID;

  G4double TotPhotons;
  G4double MeasPhotons;
  G4double PCphotons;
  G4int TrackID;

};

#endif
