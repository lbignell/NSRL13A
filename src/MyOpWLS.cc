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
// $Id: MyOpWLS.cc 71487 2013-06-17 08:19:40Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon WaveLength Shifting (WLS) Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        MyOpWLS.cc
// Description: Discrete Process -- Wavelength Shifting of Optical Photons
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     2005-07-28 - add G4ProcessType to constructor
//              2006-05-07 - add G4VWLSTimeGeneratorProfile
// mail:        gum@triumf.ca
//              jparcham@phys.ualberta.ca
//
////////////////////////////////////////////////////////////////////////

#include "MyOpWLS.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpProcessSubType.hh"

#include "G4WLSTimeGeneratorProfileDelta.hh"
#include "G4WLSTimeGeneratorProfileExponential.hh"

//Include ROOT.
#include "TROOT.h"
#include "TLeaf.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <cmath>

/////////////////////////
// Class Implementation
/////////////////////////

/////////////////
// Constructors
/////////////////

using namespace std;

MyOpWLS::MyOpWLS(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  SetProcessSubType(fOpWLS);

  theIntegralTable = NULL;
  theQYTable = NULL;
 
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }

  WLSTimeGeneratorProfile = 
       new G4WLSTimeGeneratorProfileDelta("WLSTimeGeneratorProfileDelta");

}

////////////////
// Destructors
////////////////

MyOpWLS::~MyOpWLS()
{
  if (theIntegralTable != 0) {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
  }
  delete WLSTimeGeneratorProfile;
}

////////////
// Methods
////////////

void MyOpWLS::BuildPhysicsTable(const G4ParticleDefinition&)
{
  //This method is called during initialization.
    if (!theIntegralTable) BuildThePhysicsTable();
    if(!theQYTable) BuildTheQYTable();
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
MyOpWLS::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if (verboseLevel>0) {
    G4cout << "\n** Photon absorbed! **" << G4endl;
  }
  
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  const G4MaterialPropertyVector* WLS_Intensity = 
    aMaterialPropertiesTable->GetProperty("WLSCOMPONENT"); 

  if (!WLS_Intensity)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  G4double primaryEnergy = aTrack.GetDynamicParticle()->GetKineticEnergy();

  //////////////////////////////////////////////////////////////////////////////
  //This block of code is a good place to implement my QY sampling, provided I
  //can get the incident photon energy.
  G4int NumPhotons = 0;

  G4PhysicsOrderedFreeVector* QYValues =
    (G4PhysicsOrderedFreeVector*)((*theQYTable)(aMaterial->GetIndex()));


  double theQY = 0;
  //cout << "Calculating QY... QYValues = " << QYValues;

  if(QYValues){
    //Case where energy is lower than the min energy; set to min value.
    if(primaryEnergy<QYValues->GetMinLowEdgeEnergy())
      {
	theQY = QYValues->GetMinValue();
	//	cout << "Energy = " << 1240./(primaryEnergy*1000000) << " nm, "
	//<< "QY = " << theQY << endl;
      }

    //Case where energy is higher than the max energy; set to max value.
    else if(QYValues->GetMaxLowEdgeEnergy()<primaryEnergy)
      {
	theQY = QYValues->GetMaxValue();
	//cout << "Energy = " << 1240./(primaryEnergy*1000000) << " nm, "
	//<< "QY = " << theQY << endl;
      }

    //Set to the nearest energy bin.
    else{
    theQY = QYValues->Value(primaryEnergy);
    //cout << "Energy = " << 1240./(primaryEnergy*1000000) << " nm, "
    //	 << "QY = " << theQY << endl;
    }
  }


  //Do a monte carlo trial. Later on I could implement multiple photon emission
  //if I believe that can happen.
  if(G4UniformRand()<theQY){
    NumPhotons = 1;
    //cout << "WLS photon!" << endl;
  }
  else{
    //return unchanged primary and no secondaries.
    aParticleChange.SetNumberOfSecondaries(0);
    //cout << "QY MC returned 0 photons" << endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  //The original implementation 
  //  if (aMaterialPropertiesTable->ConstPropertyExists("WLSMEANNUMBERPHOTONS")) {
  //
  // G4double MeanNumberOfPhotons = aMaterialPropertiesTable->
  //                                GetConstProperty("WLSMEANNUMBERPHOTONS");
  //
  // NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));

  // if (NumPhotons <= 0) {

        // return unchanged particle and no secondaries

  //    aParticleChange.SetNumberOfSecondaries(0);

  //    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  // }

  //}

  aParticleChange.SetNumberOfSecondaries(NumPhotons);

  //cout << "Number of Secondaries = " << NumPhotons <<  endl;

  /////////////////////////////////////////////////////////////////////////////

  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the WLS Integral for this material
  // new G4PhysicsOrderedFreeVector allocated to hold CII's

  G4double WLSTime = 0.*ns;
  G4PhysicsOrderedFreeVector* WLSIntegral = 0;

  WLSTime   = aMaterialPropertiesTable->
    GetConstProperty("WLSTIMECONSTANT");
  //WLSIntegral =
  //(G4PhysicsOrderedFreeVector*)((*theIntegralTable)(materialIndex));
   
  // Max WLS Integral
  
  //G4double CIImax = WLSIntegral->GetMaxValue();
 
  G4int NumberOfPhotons = NumPhotons;
 
  for (G4int i = 0; i < NumPhotons; i++) {

    G4double sampledEnergy;
    
    // Don't
    // Make sure the energy of the secondary is less than that of the primary

    //for (G4int j = 1; j <= 100; j++) {

        // Determine photon energy

        //G4double CIIvalue = G4UniformRand()*CIImax;
        
    //Call function to look up my Em/Ex matrix data.

    //sampledEnergy = WLSIntegral->GetEnergy(CIIvalue);

    sampledEnergy = MyOpWLS::GetEmEnergy(primaryEnergy);

    //cout << "Wavelength Shifted! Primary wavelength = "
    //	 << 1240./(primaryEnergy/eV)
    //	 << ", Sampled Wavelength = "
    //	 << 1240./(sampledEnergy/eV) << endl; 

    if (verboseLevel>1) {
      G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
      //G4cout << "CIIvalue =      " << CIIvalue << G4endl;
    }

        //if (sampledEnergy <= primaryEnergy) break;
	//}
    /////////////////////////////////////////////////////////////////////
    //I can remove this protection in my code as I'm using realistic data.
    // If no such energy can be sampled, return one less secondary, or none

    //if (sampledEnergy > primaryEnergy) {
    // if (verboseLevel>1)
    // G4cout << " *** One less WLS photon will be returned ***" << G4endl;
    // NumberOfPhotons--;
    // aParticleChange.SetNumberOfSecondaries(NumberOfPhotons);
    // if (NumberOfPhotons == 0) {
    //    if (verboseLevel>1)
    //    G4cout << " *** No WLS photon can be sampled for this primary ***"
    //           << G4endl;
    //    // return unchanged particle and no secondaries
    //    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    // }
    // continue;
    //}

    // Generate random photon direction
    
    G4double cost = 1. - 2.*G4UniformRand();
    G4double sint = std::sqrt((1.-cost)*(1.+cost));

    G4double phi = twopi*G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    
    G4double px = sint*cosp;
    G4double py = sint*sinp;
    G4double pz = cost;
    
    // Create photon momentum direction vector
    
    G4ParticleMomentum photonMomentum(px, py, pz);
    
    // Determine polarization of new photon
    
    G4double sx = cost*cosp;
    G4double sy = cost*sinp;
    G4double sz = -sint;
    
    G4ThreeVector photonPolarization(sx, sy, sz);
    
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);
    
    phi = twopi*G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);
    
    photonPolarization = cosp * photonPolarization + sinp * perp;
    
    photonPolarization = photonPolarization.unit();
    
    // Generate a new photon:
    
    G4DynamicParticle* aWLSPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
			    photonMomentum);
    aWLSPhoton->SetPolarization
      (photonPolarization.x(),
       photonPolarization.y(),
       photonPolarization.z());
    
    aWLSPhoton->SetKineticEnergy(sampledEnergy);
    
    // Generate new G4Track object:
    
    // Must give position of WLS optical photon

    G4double TimeDelay = WLSTimeGeneratorProfile->GenerateTime(WLSTime);
    G4double aSecondaryTime = (pPostStepPoint->GetGlobalTime()) + TimeDelay;

    G4ThreeVector aSecondaryPosition = pPostStepPoint->GetPosition();

    G4Track* aSecondaryTrack = 
      new G4Track(aWLSPhoton,aSecondaryTime,aSecondaryPosition);
   
    aSecondaryTrack->SetTouchableHandle(aTrack.GetTouchableHandle()); 
    // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);
    
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    
    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel>0) {
    G4cout << "\n Exiting from MyOpWLS::DoIt -- NumberOfSecondaries = " 
	   << aParticleChange.GetNumberOfSecondaries() << G4endl;  
  }
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the wavelength shifting process
// --------------------------------------------------
//

void MyOpWLS::BuildThePhysicsTable()
{
  if (theIntegralTable) return;
  
  const G4MaterialTable* theMaterialTable = 
    G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  // create new physics table
  
  if(!theIntegralTable)theIntegralTable = new G4PhysicsTable(numOfMaterials);
  
  // loop for materials
  
  for (G4int i=0 ; i < numOfMaterials; i++)
    {
      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
	new G4PhysicsOrderedFreeVector();
      
      // Retrieve vector of WLS wavelength intensity for
      // the material from the material's optical properties table.
      
      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
	aMaterial->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {

	G4MaterialPropertyVector* theWLSVector = 
	  aMaterialPropertiesTable->GetProperty("WLSCOMPONENT");

	if (theWLSVector) {
	  
	  // Retrieve the first intensity point in vector
	  // of (photon energy, intensity) pairs
	  
	  G4double currentIN = (*theWLSVector)[0];
	  
	  if (currentIN >= 0.0) {

	    // Create first (photon energy) 
	   
	    G4double currentPM = theWLSVector->Energy(0);
	    
	    G4double currentCII = 0.0;
	    
	    aPhysicsOrderedFreeVector->
	      InsertValues(currentPM , currentCII);
	    
	    // Set previous values to current ones prior to loop
	    
	    G4double prevPM  = currentPM;
	    G4double prevCII = currentCII;
	    G4double prevIN  = currentIN;
	    
	    // loop over all (photon energy, intensity)
	    // pairs stored for this material

            for (size_t j = 1;
                 j < theWLSVector->GetVectorLength();
                 j++)	    
	      {
		currentPM = theWLSVector->Energy(j);
		currentIN = (*theWLSVector)[j];
		
		currentCII = 0.5 * (prevIN + currentIN);
		
		currentCII = prevCII +
		  (currentPM - prevPM) * currentCII;
		
		aPhysicsOrderedFreeVector->
		  InsertValues(currentPM, currentCII);
		
		prevPM  = currentPM;
		prevCII = currentCII;
		prevIN  = currentIN;
	      }
	  }
	}
      }
	// The WLS integral for a given material
	// will be inserted in the table according to the
	// position of the material in the material table.

	theIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
    }
}


//This is essentially a duplicate of BuildThePhysicsTable, but rather than get
//the CDF for the "WLSCOMPONENT" variable, it gets it for "QUANTUMYIELD".
void MyOpWLS::BuildTheQYTable()
{
  if (theQYTable) return;

  cout << "////////////////////////////////////////////////////////////////////"
       << endl;
  cout << "Building the QY Table" << endl;
  cout << "////////////////////////////////////////////////////////////////////"
       << endl;

  const G4MaterialTable* theMaterialTable = 
    G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  // create new physics table
  theQYTable = new G4PhysicsTable(numOfMaterials);
  
  // loop for materials
  for (G4int i=0 ; i < numOfMaterials; i++)
    {      
      // Retrieve vector of WLS wavelength intensity for
      // the material from the material's optical properties table.
      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
	aMaterial->GetMaterialPropertiesTable();

      G4MaterialPropertyVector* theQYVector = new G4MaterialPropertyVector();

      if (aMaterialPropertiesTable) {
	theQYVector = aMaterialPropertiesTable->GetProperty("QUANTUMYIELD");
	//cout << "TheQYVector address: " << theQYVector << endl;
	//if(theQYVector){
	//cout << "Found QY property; length = " << 
	//theQYVector->GetVectorLength() << endl;
	//theQYVector->DumpValues();
	//}
	//else{cout << "theQYVector not set for this material." << endl;}
	//theQYVector->GetLowEdge();
      }
	// The WLS integral for a given material
	// will be inserted in the table according to the
	// position of the material in the material table.

	theQYTable->insertAt(i,theQYVector);
    }
}


// GetMeanFreePath
// ---------------
//
G4double MyOpWLS::GetMeanFreePath(const G4Track& aTrack,
 				         G4double ,
				         G4ForceCondition* )
{
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4double thePhotonEnergy = aParticle->GetTotalEnergy();

  G4MaterialPropertiesTable* aMaterialPropertyTable;
  G4MaterialPropertyVector* AttenuationLengthVector;
	
  G4double AttenuationLength = DBL_MAX;

  aMaterialPropertyTable = aMaterial->GetMaterialPropertiesTable();

  if ( aMaterialPropertyTable ) {
    AttenuationLengthVector = aMaterialPropertyTable->
      GetProperty("WLSABSLENGTH");
    if ( AttenuationLengthVector ){
      AttenuationLength = AttenuationLengthVector->
	Value(thePhotonEnergy);
      //G4cout << "WLS Abs length found! AttenuationLength = "
      //     << AttenuationLength << " mm"
      //     << ". Photon Energy = " << (1240./(thePhotonEnergy/eV)) << "nm"
      //     << G4endl;
    }
    else {
      //G4cout << "No WLS absorption length specified" << G4endl;
      //G4cout << "The Photon Energy = " << (1240./(thePhotonEnergy/eV)) << "nm"
      //     << G4endl;
    }
  }
  else {
    //G4cout << "No WLS absortion length specified" << G4endl;
  }
  
  return AttenuationLength;
}

void MyOpWLS::UseTimeProfile(const G4String name)
{
  if (name == "delta")
    {
      delete WLSTimeGeneratorProfile;
      WLSTimeGeneratorProfile = 
             new G4WLSTimeGeneratorProfileDelta("delta");
    }
  else if (name == "exponential")
    {
      delete WLSTimeGeneratorProfile;
      WLSTimeGeneratorProfile =
             new G4WLSTimeGeneratorProfileExponential("exponential");
    }
  else
    {
      G4Exception("MyOpWLS::UseTimeProfile", "em0202",
                  FatalException,
                  "generator does not exist");
    }
}

G4double MyOpWLS::GetEmEnergy(G4double ExEn){
  //Convert ExEn to LambEx in nm
  G4double LambEx = 1239.84187/(ExEn/eV);  
  //if(LambEx>450){
  //cout << "Sampling the Ex/Em matrix!" << endl;
  //}

  //Loop the ExEmData events until a longer wavelength event is found.
  int theEvent = 0;
  int UpperEvt = 0;
  bool interp = true;
  bool FoundEvt = false;

  //if(LambEx>450){
  //cout << "Finding the next-highest wavelength, LambEx = " <<LambEx << "... ";
  //}
  while(theEvent<ExEmData.size()){
    if(ExEmData.at(theEvent).at(0).at(0)>LambEx){
      UpperEvt = theEvent;
      FoundEvt = true;
      break;
    }
    theEvent++;
  }
  //if(LambEx>450){
  //cout << "Done! (UpperEvt = " << UpperEvt << ")." << endl;
  //cout << "Initializing some parameters...";
  //}

  int LowerEvt = 0;
  if(!FoundEvt){return ExEn;}//Wavelength too long, just return same.
  else if(UpperEvt!=0){LowerEvt = UpperEvt-1;}//interpolate
  else{interp = false;}//Don't interpolate beyond bounds.
  
  double UpperEx = ExEmData.at(UpperEvt).at(0).at(0);
  vector<double> UpperInten = ExEmData.at(UpperEvt).at(2);
  double LowerEx = 0;
  vector<double> LowerInten;
  if(interp){
    LowerEx = ExEmData.at(LowerEvt).at(0).at(0);
    LowerInten = ExEmData.at(LowerEvt).at(2);
  }
  vector<double> InterpDist;
  double RunningSum = 0;
  vector<double> InterpCDF;
  //if(LambEx>450){
    //cout << "Done! (LowerEvt = " << LowerEvt << ", LowerEx = " << LowerEx
    //   << "Calculating interpolated spectrum...";
  //}

  double theUpperInten = 0;

  for(int i = 0; i<UpperInten.size(); i++){
    if(interp){
      //This check is because my input data is a bit stupid.
      if(isnan(UpperInten.at(i))){theUpperInten = 0;}
      else{theUpperInten = UpperInten.at(i);}
      
      InterpDist.push_back(LowerInten.at(i) + (LambEx-LowerEx)*
			   (theUpperInten-LowerInten.at(i))
			   /(UpperEx-LowerEx));
      //  if(isnan(InterpDist.back())){
      //cout << "NaN found!" << endl
      //     << "i = " << i << endl
      //     << "LowerInten.at(i) = " << LowerInten.at(i) << endl
      //     << "UpperInten.at(i) = " << UpperInten.at(i) << endl
      //     << "LambEx = " << LambEx << endl
      //     << "LowerEx = " << LowerEx << endl
      //     << "UpperEx = " << UpperEx << endl;
      //}
    }
    else{InterpDist.push_back(UpperInten.at(i));}

    RunningSum += InterpDist.back();
    InterpCDF.push_back(RunningSum);
  }

  //if(LambEx>450){
  //cout << " Done!" << endl
  //   << "RunningSum = " << RunningSum << endl
  //   << "Sampling optical photon wavelength to return... ";
  //}

  //Perform a monte carlo sample
  G4double sample = G4UniformRand()*RunningSum;
  //if(LambEx>450){
  //cout << "Sample = " << sample << endl;
  //}

  G4int theIndex = 0;
  while(theIndex<InterpCDF.size()){
    if(InterpCDF.at(theIndex)>sample){
      //Get the value of wavelength at theIndex
      //  if(LambEx>450){
      //cout << "Done! Returning..." << endl;
      //}
      //Data are in nm.
      //cout << "theIndex = " << theIndex << endl;
      //cout << "returning this value: " << 1240./(ExEmData.at(LowerEvt).at(1).at(theIndex)/eV) << endl;
      return 1240./(ExEmData.at(LowerEvt).at(1).at(theIndex)/eV);
    }
    theIndex++;
  }

  cout << "ERROR: Ex/Em matrix sampling finished with no result" << endl;
  cout << "LambdaEx = " << LambEx << " nm." << endl;
  return 0;

}


void* MyOpWLS::GetPointerToValue(TBranch* theBranch, int entry,
				const char* name){
  theBranch->GetEntry(entry);
  TLeaf* theLeaf = theBranch->GetLeaf(name);
  return theLeaf->GetValuePointer();
}


//Arguments: File name, Tree name containing Ex/Em data, Branch names for
//excitation wavelength and emission intensity. Bins are assumed to be nm.
void MyOpWLS::SetExEmData(string fname){
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //Currently, the data are stored in the ExEmData vector, with the exciting
  //wavelength for event i stored in ExEmData.at(i).at(0).at(0), the emitting
  //wavelengths stored in ExEmData.at(i).at(1).at(:), and the normalized 
  //stored at ExEmData.at(i).at(2).at(:).
  //
  //Perhaps it is also better to assume that the data are in a formatted text
  //file, rather than a .root file, for portability and so that the user doesn't
  //need to have root installed.

  //I assume that the events are stored in such a way that they are
  //monotonically increasing with wavelength
  //cout << "Opening the Ex/Em root file and getting branches... ";
  TFile* f = new TFile(fname.c_str());
  TTree* theTree = (TTree*)f->Get("FluorSpec");
  TBranch* ExBranch = (TBranch*)theTree->GetBranch("Lambda_ex");
  TBranch* EmBranch = (TBranch*)theTree->GetBranch("Wavelength");
  TBranch* IntenBranch = (TBranch*)theTree->GetBranch("Intensity");
  //cout << "Done!" << endl;

  int nEntries = theTree->GetEntries();
  //cout << "Number of Entries in the tree = " << nEntries << endl;

  ExEmData.clear();

  vector<double> ExWavelength;
  vector<double> EmWavelengths;
  vector<double> EmIntensities;
  vector< vector<double> > theData;
  double theIntegral;

  for(int i = 0; i<nEntries; i++){
    ExWavelength.clear();
    theData.clear();
    theIntegral = 0;

    ExWavelength.push_back(*(double*)
			   (MyOpWLS::GetPointerToValue(ExBranch, i,
						       ExBranch->GetName())));

    EmWavelengths = *(vector<double>*)
      (MyOpWLS::GetPointerToValue(EmBranch, i, EmBranch->GetName()));

    EmIntensities = *(vector<double>*)
      (MyOpWLS::GetPointerToValue(IntenBranch, i, IntenBranch->GetName()));

    //Load the data into slots 1, 2, and 3.
    theData.push_back(ExWavelength);
    theData.push_back(EmWavelengths);

    //Normalize the emission intensities.
    for(size_t j = 0; j<EmIntensities.size(); j++){
      theIntegral += EmIntensities.at(j);
    }
    for(size_t j = 0; j<EmIntensities.size(); j++){
      EmIntensities.at(j) = EmIntensities.at(j)/theIntegral;
    }
    theData.push_back(EmIntensities);

    ExEmData.push_back(theData);

  }

  return;
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //The block of code below is the old way to read in the data

  //initialize
  //G4int EvtLower = 0;
  //G4int EvtUpper = 0;
  //G4double UpperEx; 
  //G4double LowerEx = 0;

  //Loop until a longer lambdaEx is found.
  //while(EvtUpper<nEntries){
  //UpperEx = *(G4double*)(MyOpWLS::GetPointerToValue(ExBranch, EvtUpper,
  //						     ExBranch->GetName())); 

//if(UpperEx>LambEx){
//    break;
//  }

//  EvtUpper++;
//  LowerEx = UpperEx;
//}

//if(EvtUpper!=0){
//  EvtLower = EvtUpper-1;
//}

  //  vector<double> wavelength;
  //  vector<double> IntenLower;
  //vector<double> IntenUpper;

  //TBranch* IntenBranch = (TBranch*)theTree->GetBranch("Intensity");
  //TBranch* WLBranch = (TBranch*)theTree->GetBranch("Wavelength");

  //if(IntenBranch){
  // cout << "IntenBranch returned successfully." << endl;
  //}

  //if(WLBranch){
  // wavelength = *(vector<double>*)
  // (MyOpWLS::GetPointerToValue(WLBranch, EvtUpper, WLBranch->GetName()));
  //}
  //else{
  //cout << "WLBranch not returned!" << endl;
  //return 0;
  //}

  //IntenLower = *(vector<double>*)
  //(MyOpWLS::GetPointerToValue(IntenBranch, EvtLower, IntenBranch->GetName()));

  //IntenUpper = *(vector<double>*)
  //(MyOpWLS::GetPointerToValue(IntenBranch, EvtUpper, IntenBranch->GetName()));

  //G4double IntegLower = 0;
  //G4double IntegUpper = 0;

  //cout << "IntenLower.size() = " << IntenLower.size() << endl;
  //cout << "IntenUpper.size() = " << IntenUpper.size() << endl;

  //Get integral for normalization.
  //if(IntenUpper.size()==IntenLower.size()){
  //for(int i = 0; i<IntenUpper.size(); i++){

  //	IntegLower += IntenLower.at(i);
  //	IntegUpper += IntenUpper.at(i);

	//cout << "IntenLower.at(" << i << ") = " << IntenLower.at(i) << endl;
	//cout << "IntenUpper.at(" << i << ") = " << IntenUpper.at(i) << endl;

  //  }
  //}
  //else{
  //cout << "ERROR: The Em/Ex Matrix vectors for intensity and wavelength are"
  //	 << " not the same size! Returning photon energy = 0!" << endl;
  //cout << "IntenUpper.size() = " << IntenUpper.size() << endl
  //	 << "IntenLower.size() = " << IntenLower.size() << endl;
  //  //<< "wavelength.size() = " << wavelength.size() << endl;
  //f->Close();
  //delete f;
  //return 0;
  //}

  //  cout << "IntegLower = " << IntegLower << endl 
  //	 << "IntegUpper = " << IntegUpper << endl;
  
  //  vector<double> InterpDist;
  //vector<double> InterpCDF;
  //G4double RunningSum = 0;

  //  G4double NormLowVal;
  //G4double NormUpVal;

  //Normalize the distributions, make a weighted distribution (interpolate),
  //make a cumulative distribution fn for sampling.
  //for(int i = 0; i<IntenUpper.size(); i++){

  //NormLowVal = IntenLower.at(i)/IntegLower;
  //NormUpVal = IntenUpper.at(i)/IntegUpper;

  //InterpDist.push_back(NormLowVal +  
  //			 (LambEx-LowerEx)*
  //			 (NormUpVal-NormLowVal)/(UpperEx-LowerEx));

  //RunningSum += InterpDist.back();

  //InterpCDF.push_back(RunningSum);

  //}

  //if (verboseLevel>1){
  // cout << "Ex/Em normalisation integral = " << RunningSum << endl;
  //}

  //Don't need to check normalisation, as next line takes care of it.
  //G4double sample = G4UniformRand()*RunningSum;
  // G4int theIndex = 0;

  //cout << "InterpCDF.size() = " << InterpCDF.size() << endl;

  //End of aformentioned block of code...
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
}
