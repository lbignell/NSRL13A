#include "SensitiveDetector.hh"
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
#include "G4SDManager.hh"
#include "HodoRear.hh"
#include "HodoFront.hh"
#include "PMTwinUp.hh"
#include "PMTwinDown.hh"

SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name){
  //runnum = 0;
  name = "liquid_log";

}

SensitiveDetector::~SensitiveDetector(){;}

/*This method is invoked at the beginning of each event. The argument of this method is an object of the G4HCofThisEvent class. Hits collections, where hits produced in this particular event are stored, can be associated to the G4HCofThisEvent object in this method. The hits collections associated with the G4HCofThisEvent  object during this method can be used for ``during the event processing'' digitization.*/
void SensitiveDetector::Initialize(G4HCofThisEvent* HCE){
  runnum = 0;
  TrackID = 0;
  ParentID = 0;
  Counter = 0;
  EdepThisEventUnquenched = 0;

  //initialise process array at start of each event.
  for(int i=0; i<2; i++){
    for(int j = 0; j<1000; j++){
      SecondaryArray[j][i] = 0;
      FamilyTree[j][i] = 0;
      Events[j] = 0;
      DirectGammaInts[j] = 0;
    }
  }

  //Set constants for Stopping Power calculation
  Const = (0.15353745)*((MeV*cm2)/g);
  SpdOfLight = (2.99792458*pow(10.,10.))*(cm/s);
  ElecRestEnergy = 0.511003*MeV;
  DensityEffectCorr = 0.;

  MultipleIn = false;
  KinEnIn = 0.;
  EvtType = 0;
  TrackLenInVol = 0.;
  TrackLenPrior = 0.;
  VertexX = 0.;
  VertexY = 0.;
  VertexR = 0.;

  TotOptPhotons = 0;
  dblvec.clear();
  intvec.clear();
  InitEnVec.clear();
  ParentIDVec.clear();
  WinEnVec.clear();
  ProcVec.clear();
  //InitProcVec.clear();
  for(int i = 0; i<3; i++){
    //G4cout << "Initialising vector dimensions" << G4endl;
    InitEnVec.push_back(dblvec);
    ParentIDVec.push_back(dblvec);
    WinEnVec.push_back(dblvec);
    distVecUp.push_back(dblvec);
    distVecDown.push_back(dblvec);
  }
  //G4cout << "InitEnVec.size() = " << InitEnVec.size() << G4endl;
  //G4cout << "WinEnVec.size() = " << WinEnVec.size() << G4endl;


  HitWindow = false;

  PrimEntryLocn.setX(0);
  PrimEntryLocn.setY(0);
  PrimEntryLocn.setZ(0);
  PrimMomDirn.setX(0);
  PrimMomDirn.setY(0);
  PrimMomDirn.setZ(0);

  PrimaryIsOpPhot = false;

}

/*This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. The first argument of this method is a G4Step  object of the current step. The second argument is a G4TouchableHistory object for the ``Readout geometry'' described in the next section. The second argument is NULL for the case ``Readout geometry'' is not assigned to this sensitive detector. In this method, one or more G4VHit objects should be constructed if the current step is meaningful for your detector.*/
G4bool SensitiveDetector::ProcessHits(G4Step* theStep, G4TouchableHistory*){

  //A bit of code to check for reflections:
  if((HitWindow)&&(theStep->GetTrack()->GetTrackID()==PhotonID)){
    //G4cout << "Photon reflected from UVT Acrylic, removing hit!" << G4endl;
    if(ProcVec.size()!=0){
      ProcVec.pop_back();
    }
    HitWindow = false;
  } 

  //Here's the plan: I'll classify events into different categories:
  //1. Multiple particles entered the vial and deposited energy.
  //2. Only the primary proton entered the vial and deposited energy.
  //3. Only a gamma ray entered the vial and deposited energy.
  //4. Only a secondary proton entered the vial and deposited energy.
  //5. Only an electron entered the vial and deposited energy.
  //6. Only a neutron entered the vial and deposited energy.
  //7. Only an some other particle entered the vial and deposited energy.

  //For each of these, I'll save the total energy deposit in the vial for that
  //event and the kinetic energy of the incident particle (or the sum if there
  //were multiple).

  //Check whether the particle is stepping into the volume from outside...
  if((theStep->GetTrack()->GetTrackID()!=TrackID)&&
     (theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()!="liquid_log")
     ){
    //stepping into material for the first time
    if(TrackID!=0){
      MultipleIn = true;
      EvtType = 1;
    }      
    
    TrackID = theStep->GetTrack()->GetTrackID();
    KinEnIn += theStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy();

    //Collect the vertex locations of the entering particle.
    VertexX = theStep->GetTrack()->GetVertexPosition().getX();
    VertexY = theStep->GetTrack()->GetVertexPosition().getY();
    VertexR = sqrt(VertexX*VertexX + VertexY*VertexY);

    
    TrackLenPrior += theStep->GetTrack()->GetTrackLength();
    
    string thename = (theStep->GetPreStepPoint()->GetMaterial()->GetName());

    if(theStep->GetTrack()->GetTrackID() != 1){
      thename = (theStep->GetTrack()->GetDefinition()->GetParticleName());
      if(thename == "gamma"){
	EvtType = 3;
      }
      else if(thename == "proton"){
	EvtType = 4;
      }
      else if(thename == "e-"){
	EvtType = 5;
      }
      else if(thename == "neutron"){
	EvtType = 6;
      }
      else{
	EvtType = 7;
      }
    }
    else{
      //G4cout << "PRIMARY PARTICLE" << endl;
      EvtType = 2;
      //This is called when the primary particle enters the liquid, so I can use
      //this statement to grab the initial primary particle location.
      PrimEntryLocn = theStep->GetPreStepPoint()->GetPosition();
      PrimMomDirn = theStep->GetPreStepPoint()->GetMomentumDirection();

    }
    
  }


  //Now accumulate the path length in the scintillator.
  if((theStep->GetTrack()->GetTrackID()==TrackID)&&
     (theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()!="liquid_log")
     ){
    //The particle that we're tracking, get the step length.
    TrackLenInVol+=theStep->GetStepLength();
  }

  //Need to add in alternate Edep collection, simple Edep per step.
  EdepThisEventUnquenched += theStep->GetTotalEnergyDeposit()*MeV;

  //Flag for special case where I'm doing a validation measurement using
  //optical photons.
  if((theStep->GetTrack()->GetTrackID()==1)&&
    (theStep->GetTrack()->GetParticleDefinition()->GetParticleName()
     =="opticalphoton")){ PrimaryIsOpPhot = true;}


  ////////////////////////////////////////////////////////////////////////////  
  //Get optical photons
  if((theStep->GetTrack()->GetParticleDefinition()->GetParticleName()
      =="opticalphoton")&&(theStep->GetTrack()->GetTrackID()!=1)){

    ///////First step stuff///////
    //First thing that I need to do is identify optical photons upon creation
    //within the volume.
    G4int StepNum = theStep->GetTrack()->GetCurrentStepNumber();

    //As a single photon won’t have the first step number more than once, then
    //it should be ok to use these events to increment the tally of the number
    //of optical photons generated.
    if( StepNum == 1 ){      

      //Initialization data will be read out using vectors at the end of
      //each event (don’t forget to clear it after readout). The first column
      //of the vector pertains to Cerenkov photons, the second, Scintillation.
      
      if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	 =="Cerenkov"){
	//Energy upon production.
	InitEnVec.at(0)
	  .push_back(theStep->GetTrack()->GetVertexKineticEnergy());
	ParentIDVec.at(0).push_back(theStep->GetTrack()->GetParentID());
      }
      else if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	      =="OpWLS"){
	//First step of a WLS photon.
	InitEnVec.at(1)
	  .push_back(theStep->GetTrack()->GetVertexKineticEnergy());
	//Get the particle ID of the parent.
	ParentIDVec.at(1).push_back(theStep->GetTrack()->GetParentID());
      }
      else if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	      =="Scintillation"){
	//Scintillation photon
	InitEnVec.at(2)
	  .push_back(theStep->GetTrack()->GetVertexKineticEnergy());
	//Get the particle ID of the parent.
	ParentIDVec.at(2).push_back(theStep->GetTrack()->GetParentID());
      }
      else{
	//Some other process for generating optical photons, issue warning.
	G4cout << "WARNING: an optical photon was generated by a process other"
	       << " than 'Cerenkov' or 'OpWLS'. The process was: "
	       << theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	       << G4endl;
      }

      //Increment number of opt photons.
      TotOptPhotons++;

    }
    
    
    ///////Absorption/process stuff///////
    //What I want from this section is to get an idea of the number of
    //particles Rayleigh scattered, absorbed (by the water, not by the walls),
    //and that escape through the acrylic window.
    //I don’t want to record photon transport events, or events where the
    //photons are absorbed on the black walls.
    
    //Again, I think I’ll use a vector to contain the relevant info.
    //ProcVec: | Process Type | Energy | Creator Process |
    
    if( (theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	 != "Transportation") ){
      //I'm not sure, but i think "Transportation" covers the step being limited
      //by the maximum step size or encountering a boundary.
      
      //Something happened. Log the data.
      ProcVec.push_back(theStep->GetPostStepPoint()->GetProcessDefinedStep()
			->GetProcessName());


      //The code below is obselete; the QY is now implemented in MyOpWLS.
      //if(theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
      // == "OpWLS"){
	//Apply QY to WLS photons.
	//Rather than fixed QY, implement QY as function of exciting photon En.
	//if(SampleQY(theStep->GetTrack()->GetVertexKineticEnergy())){
	  //do nothing?

	//Or, turn off process
	//if(false){

	//}
	//else{
      //theStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
      //  return true;
      //}
      //}
      
    }
    //Also need to cover the case where the post step point lands in the
    //UVT Acrylic
    else if(theStep->GetPostStepPoint()->GetMaterial()->GetName()
	    == "G4_POLYACRYLONITRILE"){
      //Log as a "Window" event.
      string windowD = "WindowDownstream";
      string windowU = "WindowUpstream";

      if(theStep->GetPostStepPoint()->GetPosition().z()<0){
	//Upstream window
	ProcVec.push_back(windowU);
      }
      else{
	//Downstream window
	ProcVec.push_back(windowD);
      }

      //Get track ID. If the very next time this method is called it is still 
      //the same track ID, then the photon was reflected. (this is doubly true
      //as all photons that are incident on the PMT window are killed).
      PhotonID = theStep->GetTrack()->GetTrackID();
      HitWindow = true;

      //This is where I can implement the 'distance from primary trajectory'
      //First, get the momentum direction of the photon and where on the window
      //it has hit. This can define the line.
      G4ThreeVector HitLocn = theStep->GetPostStepPoint()->GetPosition();
      G4ThreeVector HitDirn = theStep->GetPostStepPoint()
	->GetMomentumDirection();

      //The primary proton trajectory is modeled as follows:
      //1. the entry location into the liquid volume is logged.
      //2. the momentum direction of the first step in the liquid volume is
      //   logged.
      //3. I assume that the proton remains undeflected along its path through
      //   the liquid. This should be a reasonable approximation for most events

      //Now to calculate the distance between the two 3D lines.
      //formula d = |<PQ>.< <u> x <v> >|/ |< <u> x <v> >|
      // where P and Q are points on the lines, and <u> and <v> are vectors
      // along the lines. I'll keep this notation.
      G4ThreeVector PQ;
      PQ.setX(HitLocn.x()-PrimEntryLocn.x());
      PQ.setY(HitLocn.y()-PrimEntryLocn.y());
      PQ.setZ(HitLocn.z()-PrimEntryLocn.z());

      G4ThreeVector uXv = PrimMomDirn.cross(HitDirn);

      G4double distance = abs(PQ.dot(uXv))/(uXv.mag());

      if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	 =="Cerenkov"){
	//log photon energy for Cerenkov photons.
	//G4cout << "Logging a Cerenkov photon at Window... ";
	WinEnVec.at(0).push_back(theStep->GetPostStepPoint()->
				 GetKineticEnergy());
	if(HitLocn.z()<0){
	  distVecUp.at(0).push_back(distance);
	}
	else{
	  distVecDown.at(0).push_back(distance);
	}
      }
      else if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	      =="OpWLS"){
	WinEnVec.at(1).push_back(theStep->GetPostStepPoint()->
				 GetKineticEnergy());
	if(HitLocn.z()<0){
	  distVecUp.at(1).push_back(distance);
	}
	else{
	  distVecDown.at(1).push_back(distance);
	}
      }
      else if(theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	      =="Scintillation"){
	WinEnVec.at(2).push_back(theStep->GetPostStepPoint()->
				 GetKineticEnergy());
	if(HitLocn.z()<0){
	  distVecUp.at(2).push_back(distance);
	}
	else{
	  distVecDown.at(2).push_back(distance);
	}
      }
      else{
	//Print an error (the energies are now out of sync).
	G4cout << "ERROR: unexpected optical photon event!. "
	       << "The creator process was: "
	       << theStep->GetTrack()->GetCreatorProcess()->GetProcessName()
	       << "Don't trust the opt photon energy spectrum at the window "
	       << "for this run"
	       << G4endl;
      }


    }
    else if(theStep->GetPostStepPoint()->GetMaterial()->GetName()
	    == "G4_POLYSTYRENE"){
      //Implement the reflection probability.
      G4double ReflProb = 0;// = 0.065;
      if(G4UniformRand()<ReflProb){
	//Allow photon to live (do nothing).
      }
      else{
	//kill the photon.
	theStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }

    }
  }
  
  return true;  
}


/*This method is invoked at the end of each event. The argument of this method is the same object as the previous method. Hits collections occasionally created in your sensitive detector can be associated to the G4HCofThisEvent object.*/
void SensitiveDetector::EndOfEvent(G4HCofThisEvent*){
  //only need to do anything if there were some interactions.
  if((EdepThisEventUnquenched!=0)||(PrimaryIsOpPhot)){
    
    //get run action pointer
    RunAction* myRunAction = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
    
    if(myRunAction){
      //Tally unquenched Edeps.
      //G4cout << "Calling RunAction->TallyEdepNoQuench" << G4endl;
      myRunAction->TallyEdepNoQuench(EdepThisEventUnquenched);

      //Get the H2/H1 Edep and whether primary was incident.
      G4SDManager* SDman = G4SDManager::GetSDMpointer();

      HodoFront* pH2 =(HodoFront*)SDman->FindSensitiveDetector("hodoFront_log");
      G4double EnFront = pH2->GetEdep();
      G4bool HPFront = pH2->HadPrimary();

      HodoRear* pH1 = (HodoRear*)SDman->FindSensitiveDetector("hodoRear_log");
      G4double EnRear = pH1->GetEdep();
      G4bool HPRear = pH1->HadPrimary();

      PMTwinUp* PMTup = (PMTwinUp*)SDman->FindSensitiveDetector("PMTwinUp_log");
      vector< vector< double > > PMTEnUp = PMTup->GetEnSpec();
      vector< vector< double > > ParIDUp = PMTup->GetParentIDs();
      vector< vector< double > > MeasPMTEnUp = PMTup->GetMeasEnSpec();
      vector< vector< double > > MeasParIDUp = PMTup->GetMeasParentIDs();
      G4double PhotonsUp = PMTup->GetTotalPhotons();
      G4double MeasPhotonsUp = PMTup->GetMeasPhotons();
      //G4double PhotsOnPCUp = PMTup->GetPhotsOnPhotoCathode();

      PMTwinDown* PMTdown = 
	(PMTwinDown*)SDman->FindSensitiveDetector("PMTwinDown_log");
      vector< vector< double > > PMTEnDown = PMTdown->GetEnSpec();
      vector< vector< double > > ParIDDown = PMTdown->GetParentIDs();
      vector< vector< double > > MeasPMTEnDown = PMTdown->GetMeasEnSpec();
      vector< vector< double > > MeasParIDDown = PMTdown->GetMeasParentIDs();
      G4double PhotonsDown = PMTdown->GetTotalPhotons();
      G4double MeasPhotonsDown = PMTdown->GetMeasPhotons();
      //G4double PhotsOnPCDown = PMTdown->GetPhotsOnPhotoCathode();

      //G4cout << "Num Opt Photons (SD) = " << TotOptPhotons << G4endl;
      myRunAction->TallyEvtData(EvtType, KinEnIn, EdepThisEventUnquenched,
				TrackLenPrior, TrackLenInVol,
				VertexX, VertexY, VertexR,
				EnFront, EnRear, HPFront, HPRear,
				TotOptPhotons, InitEnVec, ParentIDVec,
				ProcVec, WinEnVec, //distVecUp, distVecDown,
				PMTEnUp, ParIDUp, PMTEnDown, ParIDDown,
				PhotonsUp, PhotonsDown,
				MeasPMTEnUp, MeasParIDUp,
				MeasPMTEnDown, MeasParIDDown,
				MeasPhotonsUp, MeasPhotonsDown);

    }
  }
}


G4double SensitiveDetector::ApplyBirksQuench(G4double En){
  //method to get the ionisation quench correction.

  //I can get the Mean Ionisation Potential in the following ways:
  //1.If I'm happy with the way G4 calculates this (see web) then I can get it
  //  via theStep()->GetTrack()->GetMaterial()->GetIonisation()
  //  ->GetMeanExcitationEnergy() method.
  //2.Define the Mean Excitation Energy as an entry in the Materials Property
  //  Table which can be accessed via theStep->GetTrack()->GetMaterial()
  //  ->GetMaterialPropertiesTable()->GetConstProperty(ref) method.
  //  This would also be a convenient way of getting kB.
  //3.Calculate it myself here.
  //4.Define it here as a number which will need to be changed depending on the
  //  scintillant used.

  //Going to go with Option 2. on advice from Li.
  //Code is implemented in ProcessHits.

  //Need to evaluate integral by either Simpson's rule or Trapezoidal rule
  //Take approximation over 10000 points.

  if(En!=0){
    G4int numPoints = 1000;
    G4double SumPart = 0;
    for(int i = 1; i<(numPoints); i++){
      SumPart += EvaluateBirks((En*i)/numPoints);
    }

    G4double BraketedBit = (((EvaluateBirks(En)+EvaluateBirks(0.))/2)+SumPart);
    G4double QuenchIntegral = BraketedBit/numPoints;
    G4double QuenchedEnergy = En*QuenchIntegral;

    //Then return the quenched energy.
    return QuenchedEnergy;
  }
  else{return 0;}
}

G4double SensitiveDetector::EvaluateBirks(G4double En){
  if(En!=0){
    G4bool Interp = false;
    G4double EnOrig;
    if(En<0.1*keV){
      EnOrig = En;
      En = 0.1*keV;
      Interp = true;
    }
    
    //Get necessary constants sorted out.
    G4double BetaSquared = 1 - pow((ElecRestEnergy/(ElecRestEnergy+En)),2);
    G4double Tau = (En/ElecRestEnergy);
    G4double Fminus = (1-BetaSquared)*(1+((Tau*Tau)/8)-((2*Tau+1)*log(2)));
    
    //Calculate the stopping power
    G4double dEdX = ScintDensity*Const*(1/BetaSquared)*ZonA*(log(pow((En/MeanExEn),2))+log(1+Tau/2)+Fminus-DensityEffectCorr);
    
    if(Interp==true){
      //Interpolate to zero for < 100eV particles.
      dEdX = dEdX*(EnOrig/(0.1*keV));
    }

    //Apply the basic Birk's formula
    G4double BirksCorr = (1/(1+kB*dEdX));

    return BirksCorr;  
  }
  else{return 0;}
}

G4bool SensitiveDetector::SampleQY(G4double PhotEn){
  //This is where I get the Daya bay QY as a function of energy. Eventually
  //I should implement this as a MaterialsPropertyVector.


  //Get the QY vector.
  DetectorConstruction* myDC = 
    (DetectorConstruction*)(G4RunManager::GetRunManager()->
			    GetUserDetectorConstruction());

  vector< vector< double > > QYdata = myDC->GetQYdata();
  //first column is lambda, second is QY.

  //get a vector iterator, starting at beginning of vector.
  std::vector<double>::reverse_iterator itWL = QYdata.at(0).rbegin();
  std::vector<double>::reverse_iterator itQY = QYdata.at(1).rbegin();

  //G4cout << "QYdata.at(0).size() = " << QYdata.at(0).size() << G4endl;

  G4double prevWL = *itWL;
  G4double prevQY = *itQY;
  G4double thisWL = 0.;
  G4double thisQY = 0.;
  G4double theQY = 0;
  while(itWL != QYdata.at(0).rend()){//Loop shouldn't break with this condition
    // unless I have an error in my code.
    itWL++;
    itQY++;
    thisWL = *itWL;
    thisQY = *itQY;
    if((prevWL<PhotEn)&&(PhotEn<thisWL)){
      //Interpolate btw QY points and do MC trial of QY.
      theQY = prevQY + (PhotEn-prevWL)*(thisQY-prevQY)/(thisWL-prevWL);
      //Now do MC trial to say whether photon is detected.
      //G4cout << "prevQY = " << prevQY << ", thisQY = " << thisQY 
      //     << ", theQY = " << theQY << G4endl;
      G4double trial = G4UniformRand();
      if(trial<=theQY){
      //Uniformly reduce the QY by a factor!
      //if(trial<=theQY*0.25){
	return true;
      }
      else{
	return false;
      }
    }
    else{//Initialise for next iteration.
      prevWL = thisWL;
      prevQY = thisQY;
    }
  }

  //If we've made it this far, the loop has finished, so return a warning.
  //G4cout<<"WARNING: a photon initiated OpWLS that was out of range of the QY "
  //	<< "data! Photon wavelength = " << 1240./(PhotEn/eV) << " nm" << G4endl;
  //After looking at the events that make it through the while loop using the
  //output above; it is only the very long wavelength photons that make it
  //though (<~1.5 eV). These will have zero QY, so return false.
  return false;


}
