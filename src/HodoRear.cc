#include "HodoRear.hh"
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


HodoRear::HodoRear(G4String name) : G4VSensitiveDetector(name){
  //runnum = 0;
  name = "HodoRear_log";
}

HodoRear::~HodoRear(){;}

/*This method is invoked at the beginning of each event. The argument of this method is an object of the G4HCofThisEvent class. Hits collections, where hits produced in this particular event are stored, can be associated to the G4HCofThisEvent object in this method. The hits collections associated with the G4HCofThisEvent  object during this method can be used for ``during the event processing'' digitization.*/
void HodoRear::Initialize(G4HCofThisEvent* HCE){
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

  HazPrimary = false;
}

G4double HodoRear::GetEdep(){return EdepThisEventUnquenched;}
G4bool HodoRear::HadPrimary(){return HazPrimary;}

/*This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. The first argument of this method is a G4Step  object of the current step. The second argument is a G4TouchableHistory object for the ``Readout geometry'' described in the next section. The second argument is NULL for the case ``Readout geometry'' is not assigned to this sensitive detector. In this method, one or more G4VHit objects should be constructed if the current step is meaningful for your detector.*/
G4bool HodoRear::ProcessHits(G4Step* theStep, G4TouchableHistory*){
  //G4cout << "Particle Stepping in Rear Vial!" << G4endl;
  
  if(theStep->GetTrack()->GetTrackID()==1){
    HazPrimary = true;
  }

  //Get the material constants.
  //MeanExEn = theStep->GetTrack()->GetMaterial()->GetMaterialPropertiesTable()->GetConstProperty("MeanExcitationEnergy");
  //ZonA = theStep->GetTrack()->GetMaterial()->GetMaterialPropertiesTable()->GetConstProperty("ZonA");
  //kB = theStep->GetTrack()->GetMaterial()->GetMaterialPropertiesTable()->GetConstProperty("kB");
  //ScintDensity = (theStep->GetTrack()->GetMaterial()->GetDensity());//(g/cm3);

  //First things first, log the particle type and parent ID.
  //G4int thisTrackID = theStep->GetTrack()->GetTrackID();
  
  //Assign Parent ID
  //FamilyTree[thisTrackID][0] = theStep->GetTrack()->GetParentID();

  //get particle type and assign numerical value (0=gamma, 1=e-, 2=e+)
  //G4String thisTypeString = theStep->GetTrack()->GetDefinition()->GetParticleName();
  // G4int thisTypeInt;
  //if(thisTypeString == "gamma"){
    //FamilyTree[thisTrackID][1]=0;
    //}
  //else if(thisTypeString == "e-"){
    //FamilyTree[thisTrackID][1]=1;
    //}
  //else if(thisTypeString == "e+"){
    //FamilyTree[thisTrackID][1]=2;
    //}
  //else{
    //G4cout << "Particle not gamma, e- or e+!" << G4endl;
    //}
  

  //Log energy from gammas that interact but don't produce secondaries (below
  //the range cut).  
  //if(theStep->GetTrack()->GetDefinition()->GetParticleName() == "gamma"){
    //if(theStep->GetTotalEnergyDeposit()!=0){
      //Log the energy as an isolated, separate Edep event.
      //DirectGammaInts[Counter] += (theStep->GetTotalEnergyDeposit())*MeV;  
      //Counter++;
      //}
    //}


  //Look for secondary electrons.
  //if(theStep->GetTrack()->GetDefinition()->GetParticleName() == "e-"){
    
    //G4int m = 0;
    
    //need to find TrackID of last gamma Parent.
    //G4int TrackID = thisTrackID;
    //G4int ParentID = 0;
    //G4cout << "Secondary Electron Stepping in the Scintillant..." <<G4endl;
    //G4cout << "This TrackID is: " << TrackID << G4endl;
    //G4cout << "Entering Loop" << G4endl;
    //while(true){
      //G4cout << "m = " << m << G4endl;
      //G4cout << "Track ID = " << TrackID << G4endl;
      //G4cout << "Particle Type (enum) = " << FamilyTree[TrackID][1] << G4endl;
      //increment counter
      //m++;
      
      //ParentID = FamilyTree[TrackID][0];
      //G4cout << "ParentID = " << ParentID << G4endl;
      //Get the ID of the parent. The if statement here protects against
      //trying to follow tracks outside the scintillant.
      //Need to double-check this doesn't exclude the primary gammas...
      //G4cout << "Parent Type (enum) = " << FamilyTree[ParentID][1] << G4endl;
      //if(FamilyTree[ParentID][1]==0){
	//Parent was a gamma, we've found 1st interaction!
	//break;
	//}
      //else{
	//if(ParentID!=0){
	  //do the loop again, one up the chain
	  //TrackID = ParentID;
	  //}
	//}

      //if(TrackID==0){
	//Either we've stepped past the primary (we shouldn't get there), or
	//this track wasn't part of the scintillant (we shouldn't get there).
	//G4cout << "thisAncestorID = 0!!! This code isn't behaving as expected! THIS MUST BE FIXED!!!" <<G4endl;
	//}
      
      //For the case of electrons entering the scintillant (usually
      //backscattered) with no gamma parent, we need an exception.
      //Amazingly, this is the only change in the code that I need to account
      //for this. Guess I'm just a genius.
      //if(m>1000){
	//G4cout << "Breaking the AncestorID loop - m too large" << G4endl;
	//break;
	//}
      //}	
    //No matter what type of event it was, need to add Edep this step to Events
    //array for processing later. AncestorID should be chosen correctly by the
    //while loop above.
    //G4cout << "Logging Edep in Scintillant (per Event)" << G4endl;
    //Events[TrackID] += (theStep->GetTotalEnergyDeposit())*MeV;
    //G4cout << "TrackID logged to: " << TrackID << G4endl;
    //G4cout << "Energy Logged (MeV) = " << theStep->GetTotalEnergyDeposit() << G4endl << G4endl;
    
    //}
  
  //Need to add in alternate Edep collection, simple Edep per step.
  EdepThisEventUnquenched += theStep->GetTotalEnergyDeposit()*MeV;

  //G4cout << "Energy deposit so far this event = " << EdepThisEventUnquenched
  //	 << G4endl;

  return true;  
}


/*This method is invoked at the end of each event. The argument of this method is the same object as the previous method. Hits collections occasionally created in your sensitive detector can be associated to the G4HCofThisEvent object.*/
void HodoRear::EndOfEvent(G4HCofThisEvent*)
{
  //only need to do anything if there were some interactions.
  if(EdepThisEventUnquenched!=0){
    
    //get run action pointer
    RunAction* myRunAction = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
    
    if(myRunAction){
      //Tally unquenched Edeps.
      //G4cout << "Calling RunAction->TallyEdepRear" << G4endl;
      //myRunAction->TallyEdepRear(EdepThisEventUnquenched);
      
      //Quench the Edep of this event.
      //G4double EdepQuenchedPerEvent=ApplyBirksQuench(EdepThisEventUnquenched);
      //Tally this.
      //myRunAction->TallyEdepQuenchPerEvt(EdepQuenchedPerEvent);

      //G4double EdepQuenchedPerInt = 0;
    
      //G4double EventsSum = 0;
      //for(int i=0; i<1000; i++){
      //EventsSum += Events[i];
      //}

      //Quench the Edeps per interaction.
      //for(int i=0; i<1000; i++){
      //if(Events[i]!=0){
	  //Apply quench to every member of the Events array (Q(0)=0).
      //  EdepQuenchedPerInt += ApplyBirksQuench(Events[i]);
      //}

      //if(DirectGammaInts[i]!=0){
	  //Apply quench to gamma ints.
      //  EdepQuenchedPerInt += ApplyBirksQuench(DirectGammaInts[i]);
      ///}

      //}
      //Tally this.
      // myRunAction->TallyEdepQuenchPerInt(EdepQuenchedPerInt);

    }
  }
}


G4double HodoRear::ApplyBirksQuench(G4double En){
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

G4double HodoRear::EvaluateBirks(G4double En){
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
