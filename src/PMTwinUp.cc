#include "PMTwinUp.hh"
#include "G4Step.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4SteppingManager.hh"
#include <iterator>
#include "G4TrackVector.hh"
#include "G4TrackStatus.hh"

#include "DetectorConstruction.hh"

PMTwinUp::PMTwinUp(G4String name) : G4VSensitiveDetector(name){
  //runnum = 0;
  name = "PMTwinUp_log";
}

PMTwinUp::~PMTwinUp(){;}

/*This method is invoked at the beginning of each event. The argument of this method is an object of the G4HCofThisEvent class. Hits collections, where hits produced in this particular event are stored, can be associated to the G4HCofThisEvent object in this method. The hits collections associated with the G4HCofThisEvent  object during this method can be used for ``during the event processing'' digitization.*/
void PMTwinUp::Initialize(G4HCofThisEvent* HCE){

  runnum = 0;
  HazPrimary = false;
  TotPhotons = 0;
  MeasPhotons = 0;
  PCphotons = 0;
  TrackID = 0;
  dblvec.clear();
  PhotEn.clear();
  ParentID.clear();
  PMTFireParentID.clear();
  PMTFirePhotEnIn_nm.clear();

  for(int i = 0; i<3; i++){
    PhotEn.push_back(dblvec);
    ParentID.push_back(dblvec);
    PMTFirePhotEnIn_nm.push_back(dblvec);
    PMTFireParentID.push_back(dblvec);
  }

}

vector< vector< double > >& PMTwinUp::GetEnSpec(){return PhotEn;}
vector< vector< double > >& PMTwinUp::GetParentIDs(){return ParentID;}
vector< vector< double > >& PMTwinUp::GetMeasEnSpec()
{return PMTFirePhotEnIn_nm;}
vector< vector< double > >& PMTwinUp::GetMeasParentIDs()
{return PMTFireParentID;}
G4double PMTwinUp::GetTotalPhotons(){return TotPhotons;}
G4double PMTwinUp::GetMeasPhotons(){return MeasPhotons;}
G4double PMTwinUp::GetPhotsOnPhotoCathode(){return PCphotons;}

/*This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. The first argument of this method is a G4Step  object of the current step. The second argument is a G4TouchableHistory object for the ``Readout geometry'' described in the next section. The second argument is NULL for the case ``Readout geometry'' is not assigned to this sensitive detector. In this method, one or more G4VHit objects should be constructed if the current step is meaningful for your detector.*/
G4bool PMTwinUp::ProcessHits(G4Step* theStep, G4TouchableHistory*){

  //Get optical photons
  if(theStep->GetTrack()->GetParticleDefinition()->GetParticleName()
     =="opticalphoton"){
    
    //Check the track ID to avoid re-counting particles that don't hit the
    //effective area on their first step.
    if(theStep->GetTrack()->GetTrackID()!=TrackID){
      //Increment counter
      TotPhotons++;
      //Update TrackID
      TrackID = theStep->GetTrack()->GetTrackID();
    }

    DetectorConstruction* myDC = 
      (DetectorConstruction*)(G4RunManager::GetRunManager()->
			      GetUserDetectorConstruction());

    G4double thePMTGap = myDC->GetPMTGap();
    //G4cout << "thePMTGap = " << thePMTGap << G4endl;
    //Calculate photocathode location.
    G4double PMTglassOffset = (136.36*mm/2) + 4.7625*mm + (5*mm);
    G4double PClocation = -(PMTglassOffset+thePMTGap);

    //G4cout << "PClocation = " << PClocation << G4endl;
    //    G4cout << "Photon Z position = "
    //	   << theStep->GetPostStepPoint()->GetPosition().z() << G4endl;

    //Get photons hitting photocathode edge.
    //if(theStep->GetPostStepPoint()->GetPosition().z()==PClocation){
      //G4cout << "Photon at Photocathode!" << G4endl;

      //Test to see if photons are hitting the 'effective area' of the photo-
      //cathode. The R7723 used has an effective area diameter of 46 mm. The
      //model used is a uniform area characterised entirely by the PMT QE (i.e.
      //doesn't depend on photon angle of incidence), with zero response outside
      //the effective area.
      //This is the location of the centre of the PMT glass along the Y-axis;
      G4double CentreY = (-(139.48*mm)/2 + 34.5*mm) -25.76*mm;
      //Centre along X-axis is the origin.
      G4double PhotX = theStep->GetPostStepPoint()->GetPosition().x();
      G4double AbsPhotY = theStep->GetPostStepPoint()->GetPosition().y();
      G4double PhotY = AbsPhotY-CentreY;
      G4double PhotR = sqrt(PhotX*PhotX+PhotY*PhotY);

      //Add a semi-arbitrary inefficiency factor to the PMT response, to
      //optimise the output to that measured.
      G4double FFactor = G4UniformRand();

      //if((PhotR<=(46.*mm/2))&&(FFactor<0.88)){
      if((PhotR<=(44.*mm/2))&&(FFactor<1)){
      //Photon hitting effective area.
	//G4cout << "Photon is hitting effective area!" << G4endl;

	G4bool IsDetected =
	  PMTHazFired(theStep->GetPostStepPoint()->GetKineticEnergy());
	
	if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	   =="Cerenkov"){
	  //Log Cerenkov photons
	  PhotEn.at(0).push_back(theStep->GetPostStepPoint()->GetKineticEnergy());
	  ParentID.at(0).push_back(theStep->GetTrack()->GetParentID());
	  if(IsDetected){
	    MeasPhotons++;
	    //First, convert En (MeV) to wavelength (nm).
	    G4double Wavelength = 1239.84187/
	      (theStep->GetPostStepPoint()->GetKineticEnergy()/eV);
	    PMTFirePhotEnIn_nm.at(0).push_back(Wavelength);
	    PMTFireParentID.at(0).push_back(theStep->GetTrack()->GetParentID());
	  }
	}
	else if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
		=="OpWLS"){
	  //Log scintillation photons.
	  PhotEn.at(1).push_back(theStep->GetPostStepPoint()->GetKineticEnergy());
	  ParentID.at(1).push_back(theStep->GetTrack()->GetParentID());
	  if(IsDetected){
	    MeasPhotons++;
	    //First, convert En (MeV) to wavelength (nm).
	    G4double Wavelength = 1239.84187/
	      (theStep->GetPostStepPoint()->GetKineticEnergy()/eV);
	    PMTFirePhotEnIn_nm.at(1).push_back(Wavelength);
	    PMTFireParentID.at(1).push_back(theStep->GetTrack()->GetParentID());
	  }
	}
	else if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
		=="Scintillation"){
	  //Log scintillation photons.
	  PhotEn.at(2).push_back(theStep->GetPostStepPoint()
				 ->GetKineticEnergy());
	  ParentID.at(2).push_back(theStep->GetTrack()->GetParentID());
	  if(IsDetected){
	    MeasPhotons++;
	    //First, convert En (MeV) to wavelength (nm).
	    G4double Wavelength = 1239.84187/
	      (theStep->GetPostStepPoint()->GetKineticEnergy()/eV);
	    PMTFirePhotEnIn_nm.at(2).push_back(Wavelength);
	    PMTFireParentID.at(2).push_back(theStep->GetTrack()->GetParentID());
	  }
	}
	else{
	  //Print a warning.
	  G4cout << "WARNING: a photon was incident on the upstream PMT from an"
		 << " unknown process! The process name was: "
		 << theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
		 << G4endl;
	}

	theStep->GetTrack()->SetTrackStatus(fStopAndKill);

      }
      //else{
      //G4cout << "Photon is NOT hitting effective area!" << G4endl;
      //}
      // }
  }

  
  return true;  
}


/*This method is invoked at the end of each event. The argument of this method is the same object as the previous method. Hits collections occasionally created in your sensitive detector can be associated to the G4HCofThisEvent object.*/
void PMTwinUp::EndOfEvent(G4HCofThisEvent*)
{
  
  //get run action pointer
  //RunAction* myRunAction = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
  
  //if(myRunAction){
    
  //}
  
}

G4bool PMTwinUp::PMTHazFired(G4double Lambda){

  //Energy units in QE file are MeV.

  //Get the QE vector.
  DetectorConstruction* myDC = 
    (DetectorConstruction*)(G4RunManager::GetRunManager()->
			    GetUserDetectorConstruction());

  vector< vector< double > > QEdata = myDC->GetQEdata();
  //first column is lambda, second is QE.

  //get a vector iterator, starting at beginning of vector.
  std::vector<double>::reverse_iterator itWL = QEdata.at(0).rbegin();
  std::vector<double>::reverse_iterator itQE = QEdata.at(1).rbegin();

  //G4cout << "QEdata.at(0).size() = " << QEdata.at(0).size() << G4endl;

  G4double prevWL = *itWL;
  G4double prevQE = *itQE;
  G4double thisWL = 0.;
  G4double thisQE = 0.;
  G4double theQE = 0;
  while(itWL != QEdata.at(0).rend()){//Loop shouldn't break with this condition
    // unless I have an error in my code.
    itWL++;
    itQE++;
    thisWL = *itWL;
    thisQE = *itQE;
    if((prevWL<Lambda)&&(Lambda<thisWL)){
      //Interpolate btw QE points and do MC trial of QE.
      theQE = prevQE + (Lambda-prevWL)*(thisQE-prevQE)/(thisWL-prevWL);
      //Now do MC trial to say whether photon is detected.
      //G4cout << "prevQE = " << prevQE << ", thisQE = " << thisQE 
      //     << ", theQE = " << theQE << G4endl;
      G4double trial = G4UniformRand();
      if(trial<=theQE){
	return true;
      }
      else{
	return false;
      }
    }
    else{//Initialise for next iteration.
      prevWL = thisWL;
      prevQE = thisQE;
    }
  }

  //If we've made it this far, the loop has finished, so return a warning.
  //G4cout<<"WARNING: a photon has hit PMTdown that was out of range of the QE "
  //	 << "data! Photon wavelength = " << Lambda << " MeV" << G4endl;
  //After looking at the events that make it through the while loop using the
  //output above; it is only the very long wavelength photons that make it
  //though (<~1.5 eV). These will have zero QE, so return false.
  return false;

}
