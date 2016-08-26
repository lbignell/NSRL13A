//This will be a simple goemetry just to check and see if I can get it working

#include "DetectorConstruction.hh"
//include header files for all classes used in the methods
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4RotationMatrix.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Isotope.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>

//for sensitive detector definition
#include "SensitiveDetector.hh"
#include "HodoRear.hh"
#include "HodoFront.hh"
#include "PMTwinUp.hh"
#include "PMTwinDown.hh"
#include "G4SDManager.hh"

#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "DetectorMessenger.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4LogicalBorderSurface.hh"

using namespace std;

//constructor / destructor do nothing
DetectorConstruction::DetectorConstruction(){
  AbsThickness = 233*mm;
  PMTGap = 0.1*mm;
  DetMess = new DetectorMessenger(this);
  //QEdata.clear();
}

DetectorConstruction::~DetectorConstruction(){ 

}

void DetectorConstruction::SetAbsThickness(G4double value){
  AbsThickness = value;
}

G4double DetectorConstruction::GetAbsThickness(){
  return AbsThickness;
}

void DetectorConstruction::SetPMTGap(G4double value){
  PMTGap = value;
}

G4double DetectorConstruction::GetPMTGap(){
  return PMTGap;
}

vector< vector< double > >& DetectorConstruction::GetQEdata(){
  return QEdata;
}

vector< vector< double > >& DetectorConstruction::GetQYdata(){
  return QYdata;
}

void DetectorConstruction::UpdateGeometry(){
  //My original code.
  //G4RunManager::GetRunManager()->DefineWorldVolume(Construct(), true);
  //G4RunManager::GetRunManager()->GeometryHasBeenModified();


  //Code taken from:
  //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/eventtrackmanage/970/1/2.html?inline=-1
  //Delete existing geom.
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
 
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  //G4RunManager::GetRunManager->PhysicsHasBeenModified();
  //G4RegionStore::GetInstance()->UpdateMaterialList(experimentalHall_phys);

}

G4VPhysicalVolume* DetectorConstruction::Construct(){
/* materials definition */

/*define the elements that will be used in our materials*/
//define hydrogen 
  G4double A = 1.01 * g/mole;
  G4double Z = 1;
  G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

  //define oxygen
  A = 16.0 * g/mole;
   Z = 8;
  G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

  //define nitrogen
  A = 14.0 * g/mole;
  Z = 7;
  G4Element* elN = new G4Element("Nitrogen", "N", Z, A);

  //define carbon
  A = 12.0107 * g/mole;
  Z = 6;
  G4Element* elC = new G4Element("Carbon", "C", Z, A);

  //define Silicon
  A = 28.086 * g/mole;
  Z = 14;
  G4Element* elSi = new G4Element("Silicon", "Si", Z, A);

  //define sodium
  A = 22.990 * g/mole;
  Z = 11;
  G4Element* elNa = new G4Element("Sodium", "Na", Z, A);

  //define Calcium
  A = 40.08 * g/mole;
  Z = 20;
  G4Element* elCa = new G4Element("Calcium", "Ca", Z, A);

  //define Phosphorus
  A = 30.973761*g/mole;
  Z = 15;
  G4Element* elP = new G4Element("Phosphorus", "P", Z, A);

  //define Sulphur
  A = 32.065*g/mole;
  Z = 16;
  G4Element* elS = new G4Element("Sulphur", "S", Z, A);

  //define Aluminium
  A = 26.961538*g/mole;
  Z = 13;
  G4Element* elAl = new G4Element("Aluminium", "Al",Z,A);

  //define Iron
  A = 55.845*g/mole;
  Z = 26;
  G4Element* elFe = new G4Element("Iron", "Fe",Z,A);
 
  //Define Copper
  A = 63.546*g/mole;
  Z = 29;
  G4Element* elCu = new G4Element("Copper", "Cu",Z,A);

  //Define Manganese
  A = 54.938045*g/mole;
  Z = 25;
  G4Element* elMn = new G4Element("Manganese", "Mn", Z, A);

  //Define Magnesium
  A = 24.3050*g/mole;
  Z = 12;
  G4Element* elMg = new G4Element("Magnesium", "Mg", Z, A);

  //Define Titanium
  A = 47.867*g/mole;
  Z = 22;
  G4Element* elTi = new G4Element("Titanium", "Ti", Z, A);

  //Define Chromium
  A = 51.9961*g/mole;
  Z = 24;
  G4Element* elCr = new G4Element("Chromium", "Cr",Z,A);

  //Define Zinc
  A = 65.409*g/mole;
  Z = 30;
  G4Element* elZn = new G4Element("Zinc", "Zn",Z,A);

  //Define Boron
  A = 10.811*g/mole;
  Z = 5;
  G4Element* elB = new G4Element("Boron", "B",Z,A);

  //Define Potassium
  A = 39.0983*g/mole;
  Z = 19;
  G4Element* elK = new G4Element("Potassium", "K",Z,A);

  //constructor of the G4Material class requires arguments: string 
  //conG4OpticalSurfacetaining name of material, density, number of elements
  G4Material* water = new G4Material("water", 1.0 * g/cm3, 2);
  water->AddElement(elH,2);
  water->AddElement(elO,1);

  /*now we define air for the world volume*/
  G4Material* air = new G4Material("dry air", 0.01*mg/cm3, 2, kStateGas, 293*kelvin, 1*atmosphere);
  //we can also specify the percentage (by mass) composition
  air->AddElement(elN, 75*perCent);
  air->AddElement(elO, 25*perCent);

  //High Density Polyethylene
  //G4Material* HDPE = new G4Material("HD Poly-Ethylene", 0.941*g/cm3, 2);
  //G4Material* HDPE = new G4Material("HD Poly-Ethylene", 0.93*g/cm3, 2);
  G4Material* HDPE = new G4Material("HD Poly-Ethylene", 0.97*g/cm3, 2);
  HDPE->AddElement(elC, 33.33333*perCent);
  HDPE->AddElement(elH, 66.66667*perCent);

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* water_nist = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* air_nist = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* HDPE_nist = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* Polystyrene_nist = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material* Acrylic_nist = nist->FindOrBuildMaterial("G4_POLYACRYLONITRILE");
  G4Material* glass_nist = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Element* C_nist = nist->FindOrBuildElement("C");
  G4Element* H_nist = nist->FindOrBuildElement("H");

  G4Material* Scint = new G4Material("LABscintillator", 0.86*g/cm3, 2);
  //Scint->AddElement(C_nist, (18/(18+30)));
  //Scint->AddElement(H_nist, (30/(18+30)));
  Scint->AddElement(C_nist, 18*12);
  Scint->AddElement(H_nist, 30*1);


  //Make WbLS as used in expt.
  G4double WbLSfraction = 0.999999;
  G4double WbLSdensity = (water_nist->GetDensity())*(1-WbLSfraction) + 
    (Scint->GetDensity())*WbLSfraction;
  G4Material* WbLS = new G4Material("WbLS", WbLSdensity, 2);
  //kStateLiquid, 293*kelvin, 1*atmosphere);
  WbLS->AddMaterial(water_nist, (1-WbLSfraction));
  WbLS->AddMaterial(Scint, WbLSfraction);

  //Aluminium Alloy used in the walls Alloy #6063
  //G4Material* AlAlloy=new G4Material("Al Alloy", 2680.*kg/m3, 9, kStateSolid);
  //Adding elements in this way seems pointless now but is useful for later
  //AlAlloy->AddElement(elAl, 98.5*perCent);
  //AlAlloy->AddElement(elSi, 0.40*perCent);
  //AlAlloy->AddElement(elFe, 0.175*perCent);
  //AlAlloy->AddElement(elCu, 0.05*perCent);
  //AlAlloy->AddElement(elMn, 0.05*perCent);
  //AlAlloy->AddElement(elMg, 0.675*perCent);
  //AlAlloy->AddElement(elCr, 0.05*perCent);
  //AlAlloy->AddElement(elZn, 0.05*perCent);
  //AlAlloy->AddElement(elTi, 0.05*perCent);



  //Getting data to fill the acrylic MPT.
  G4MaterialPropertiesTable* MPTacrylic = new G4MaterialPropertiesTable();

  //Test out the capability of the function to read out the data from the opt
  //property files.
  //cout << "Reading data from optical property file" << G4endl;
  FILE* FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/acrylic_abslength.vec", "r");
  //FILE* FilePtr = fopen("r7723.vec", "r");
  cout << "File open! Pointer = " << FilePtr << G4endl;
  //vector < vector < double > >& thedata;
  if(FilePtr!=0){
    //call the function
    GetOptInfo(FilePtr, 1*mm);
  }
  else{
    G4cout << "ERROR: could not open file!" << G4endl;
  }

  G4cout << "Successfully read file!" << G4endl;
  G4int cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }


  //Writing data to material properties table. Arguments: variable name, photon
  //energy, property values, number of entries.
  //Geant4 requires that the optical properties be passed as arrays, so I'll 
  //copy the vectors to arrays.
  const G4int size1 = thedata.at(0).size();
  G4double En1[size1];
  //std::copy(thedata.at(0).begin(), thedata.at(0).end(), En1);
  G4double AbsLen1[size1];//can assume vectors are same size
  //std::copy(thedata.at(1).begin(), thedata.at(1).end(), AbsLen1);
  for(int i = 0; i<size1; i++){
    En1[i] = thedata.at(0).at(i);
    AbsLen1[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex2[size2] << G4endl;
  }
  MPTacrylic->AddProperty("ABSLENGTH", En1, AbsLen1, size1);

  thedata.clear();

  //Read out the acryling refractive index data
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/acrylic_rindex.vec", "r");
  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }

  const G4int size2 = thedata.at(0).size();
  G4double En2[size2];
  //std::copy(thedata.at(0).begin(), thedata.at(0).end(), En2);
  G4double Rindex2[size2];
  for(int i = 0; i<size2; i++){
    En2[i] = thedata.at(0).at(i);
    Rindex2[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex2[size2] << G4endl;
  }
  //std::copy(thedata.at(1).begin(), thedata.at(1).end(), Rindex2);
  MPTacrylic->AddProperty("RINDEX", En2, Rindex2, size2);

  //MPTacrylic->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
  //MPTacrylic->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries);
  
  //MPTacrylic->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  //MPTacrylic->AddConstProperty("RESOLUTIONSCALE",1.0);
  //MPTacrylic->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  //MPTacrylic->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  //MPTacrylic->AddConstProperty("YIELDRATIO",0.8);
  
  Acrylic_nist->SetMaterialPropertiesTable(MPTacrylic);

  G4MaterialPropertiesTable* MPTwater = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* MPTWbLS = new G4MaterialPropertiesTable();

  thedata.clear();
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/water_abslength.vec", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WaterAbsSeg_FULL.csv", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WaterAbsSeg_PARTIAL2.csv", "r");
  //essentially de-activating this process for LS (it is handled by the QY) so
  //set crazy long abs length.
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/LongAbsLength.csv", "r");

  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1*mm);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }

  G4cout << "thedata.at(0).at(0) = " << thedata.at(0).at(0) << G4endl;
  G4cout << "thedata.at(0).at(1) = " << thedata.at(0).at(1) << G4endl;
  G4cout << "thedata.at(1).at(0) = " << thedata.at(1).at(0) << G4endl;
  G4cout << "thedata.at(1).at(1) = " << thedata.at(1).at(1) << G4endl;

  const G4int size3 = thedata.at(0).size();
  G4double En3[size3];
  //std::copy(thedata.at(0).begin(), thedata.at(0).end(), En3);
  G4double AbsLen3[size3];
  //std::copy(thedata.at(1).begin(), thedata.at(1).end(), AbsLen3);
  for(int i = 0; i<size3; i++){
    En3[i] = thedata.at(0).at(i);
    AbsLen3[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex2[size2] << G4endl;
  }
  MPTwater->AddProperty("ABSLENGTH", En3, AbsLen3, size3);
  MPTWbLS->AddProperty("ABSLENGTH", En3, AbsLen3, size3);

  thedata.clear();
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/water_rindex.vec", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/water_rindex_edit.vec", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WaterRindexSeg_FULL.csv", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WaterRindexSeg_PARTIAL2.csv", "r");
  //DB LS Rindex, while monotonically increasing
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/DBScintRindex.csv", "r");
  //DB LS Rindex, to max value
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/DBScintRindex_toMaxVal.csv", "r");


  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }

  G4cout << "thedata.at(0).at(0) = " << thedata.at(0).at(0) << G4endl;
  G4cout << "thedata.at(0).at(1) = " << thedata.at(0).at(1) << G4endl;
  G4cout << "thedata.at(1).at(0) = " << thedata.at(1).at(0) << G4endl;
  G4cout << "thedata.at(1).at(1) = " << thedata.at(1).at(1) << G4endl;


  const G4int size4 = thedata.at(0).size();
  G4double En4[size4];
  //std::copy(thedata.at(0).begin(), thedata.at(0).end(), En4);
  G4double Rindex4[size4];
  //std::copy(thedata.at(1).begin(), thedata.at(1).end(), Rindex4);
  for(int i = 0; i<size4; i++){
    En4[i] = thedata.at(0).at(i);
    Rindex4[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex4[i] << G4endl;
  }
  MPTwater->AddProperty("RINDEX", En4, Rindex4, size4);
  MPTWbLS->AddProperty("RINDEX", En4, Rindex4, size4);

  G4MaterialPropertyVector* theRindex = MPTwater->GetProperty("RINDEX");
  G4cout << "Rindex max value = " << theRindex->GetMaxValue() << G4endl;
  G4cout << "Rindex min value = " << theRindex->GetMinValue() << G4endl;
  G4cout << "Rindex at 200 nm (6.1992 eV): "
	 << theRindex->Value(6.1992*eV) << G4endl;
  G4cout << "MinLowEdgeEnergy = " << theRindex->GetMinLowEdgeEnergy() << G4endl;
  G4cout << "MaxLowEdgeEnergy = " << theRindex->GetMaxLowEdgeEnergy() << G4endl;

  //This is the WLS photon energy spectrum.
  //It has not been normalised.
  thedata.clear();
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSem.csv", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSemExtrap.csv", "r");
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/scintillator_fastcomponent.vec", "r");
  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }
  
  const G4int size5 = thedata.at(0).size();
  G4double En5[size5];
  //std::copy(thedata.at(0).begin(), thedata.at(0).end(), En5);
  G4double WLSEm5[size5];
  //std::copy(thedata.at(1).begin(), thedata.at(1).end(), WLSEm5);
  for(int i = 0; i<size5; i++){
    En5[i] = thedata.at(0).at(i);
    WLSEm5[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex2[size2] << G4endl;
  }
  MPTWbLS->AddProperty("WLSCOMPONENT", En5, WLSEm5, size5);


  //This is the Scint photon energy spectrum.
  //It has not been normalised.
  thedata.clear();
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSem.csv", "r");
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/scintillator_fastcomponent.vec", "r");
  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }


  //The relevant property name is "FASTCOMPONENT".
  const G4int size6 = thedata.at(0).size();
  G4double En6[thedata.at(0).size()];
  std::copy(thedata.at(0).begin(), thedata.at(0).end(), En6);
  G4double ScintEm6[thedata.at(1).size()];
  std::copy(thedata.at(1).begin(), thedata.at(1).end(), ScintEm6);
  for(int i = 0; i<size6; i++){
    En6[i] = thedata.at(0).at(i);
    ScintEm6[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex2[size2] << G4endl;
  }

  //G4MaterialPropertiesTable* Scnt_MPT = new G4MaterialPropertiesTable();

  MPTWbLS->AddProperty("FASTCOMPONENT", En6, ScintEm6, size6);
  //Scnt_MPT->AddProperty("SLOWCOMPONENT", Scnt_PP, Scnt_SLOW, NUMENTRIES);

  MPTWbLS->AddConstProperty("SCINTILLATIONYIELD", 11522./MeV);
  //The resolution yield is indeterminate. It will affect the breadth of the
  //photon number distribution, and can therefore be tuned to the measured value
  //if that makes sense.
  MPTWbLS->AddConstProperty("RESOLUTIONSCALE", 1.0);
  MPTWbLS->AddConstProperty("FASTTIMECONSTANT",  1.*ns);
  //Scnt_MPT->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  MPTWbLS->AddConstProperty("YIELDRATIO", 1.);

  //Scnt->SetMaterialPropertiesTable(Scnt_MPT);



  water_nist->SetMaterialPropertiesTable(MPTwater);
  
  //Now set the properties for air...
  const G4int nEntries = 32 + 26;
  G4double PhotonEnergy[nEntries] =
    { 1.5*eV, 1.6*eV, 1.7*eV, 1.8*eV, 1.9*eV,
      2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV,
      4.2*eV, 4.3*eV, 4.4*eV, 4.5*eV, 4.6*eV,
      4.7*eV, 4.8*eV, 4.9*eV, 5.0*eV, 5.1*eV,
      5.2*eV, 5.3*eV, 5.4*eV, 5.5*eV, 5.6*eV,
      5.7*eV, 5.8*eV, 5.9*eV, 6.0*eV, 6.1*eV,
      6.2*eV};

  //Now set the properties for air...
  //extra values have been added arbritrarily
  G4double RindexAir[nEntries] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00};

  G4double AbsLenAir[nEntries] =
    { 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km, 1.00*km,
      1.00*km, 1.00*km};

  G4MaterialPropertiesTable* MPTair = new G4MaterialPropertiesTable();
  MPTair->AddProperty("RINDEX", PhotonEnergy, RindexAir, nEntries);
  MPTair->AddProperty("ABSLENGTH", PhotonEnergy, AbsLenAir, nEntries);

  air_nist->SetMaterialPropertiesTable(MPTair);


  //Calculate the refractive index of fused silica (from 
  //http://refractiveindex.info/legacy/?group=GLASSES&material=F_SILICA
  //).
  G4double RindexSilica[nEntries];
  for(int i = 0; i<nEntries; i++){
    //calculate the wavelength in micrometers.
    G4double WLum = (1239.84187/PhotonEnergy[i])/1000;
    RindexSilica[i] =
      sqrt( 1 +
	    (0.6961663*WLum*WLum)/(WLum*WLum-0.0684043*0.0684043) +
	    (0.4079426*WLum*WLum)/(WLum*WLum-0.1162414*0.1162414) +
	    (0.8974794*WLum*WLum)/(WLum*WLum-9.896161*9.896161) );
  }
  
  G4double AbsLenSilica[nEntries] =
    { 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m, 1.00*m,
      1.00*m, 1.00*m};

  G4MaterialPropertiesTable* MPTsilica = new G4MaterialPropertiesTable();
  MPTsilica->AddProperty("RINDEX", PhotonEnergy, RindexSilica, nEntries);
  MPTsilica->AddProperty("ABSLENGTH", PhotonEnergy, AbsLenSilica, nEntries);

  glass_nist->SetMaterialPropertiesTable(MPTsilica);

  G4double RindexBox[nEntries];
  G4double AbsLenBox[nEntries];
  G4double ReflBox[nEntries];
  G4double TransmitEffBox[nEntries];
  for(int i = 0; i<nEntries; i++){
    RindexBox[i] = 1.6;
    AbsLenBox[i] = 1*nm;
    ReflBox[nEntries] = 0.05;//0.5% reflectance
    TransmitEffBox[nEntries] = 0.5;
  }

  G4MaterialPropertiesTable* MPTbox = new G4MaterialPropertiesTable();
  MPTbox->AddProperty("REFLECTIVITY", PhotonEnergy, ReflBox, nEntries);
  //MPTbox->AddProperty("EFFICIENCY", PhotonEnergy, TransmitEffBox, nEntries);
  //MPTbox->AddProperty("RINDEX", PhotonEnergy, RindexBox, nEntries);
  //MPTbox->AddProperty("ABSLENGTH", PhotonEnergy, AbsLenBox, nEntries);

  Polystyrene_nist->SetMaterialPropertiesTable(MPTbox);

  //Implement QE data grab here.
  thedata.clear();
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/r7723.vec", "r");
  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }
  QEdata = thedata;

  //Implement QY data grab here
  thedata.clear();
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/QYDayaBay.csv", "r");
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/QYDayaBayZeroBelow200nm.csv", "r");
  if(FilePtr!=0){
    GetOptInfo(FilePtr, 1);
  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }
  QYdata = thedata;

  //WbLS implementation:
  // -- Absorption  = water absorption.
  // -- Rindex = water Rindex.
  // -- WLS absorption = WbLSExtraAbs.
  // -- WLS emission = scint spec.
  // -- WLS time constant = 1*ns (arbitrary for now).

  //Now get the WbLS absorption
  thedata.clear();
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs.csv", "r");
  //Extrap to benzene abs data, then approximate benzene trend.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs_ExtrapExtrapDown.csv", "r");
  //Extrap to scaled benzene data, then follow water abs.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs_ExtrapFollowWater.csv", "r");
  //Extrap short wavelength absorption, then assume constant for shorter.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs_ExtrapFlat.csv", "r");
  //Extrap to unscaled benzene data, then remain constant.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs_ExtrapFlatScaleOnlyDilution.csv", "r");
  //Assume 4x concentration of abs data, benzene @ short w/L then flat
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLS4xAbs_BenzeneShortWLFlat.csv", "r");
  //Assume 8x concentration of abs data, benzene @ short w/L then flat
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLS8xAbs_BenzeneShortWLFlat.csv", "r");
  //Assume 100x concentration of abs data, benzene @ short w/L then flat
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLS100xAbs_BenzeneShortWLFlat.csv", "r");
  //Assume 1000x concentration of abs data, benzene @ short w/L then flat
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLS1000xAbs_BenzeneShortWLFlat.csv", "r");
  //Assume 10000x concentration of abs data, benzene @ short w/L then flat
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLS10000xAbs_BenzeneShortWLFlat.csv", "r");
  //Normalise to DAYA BAY water, extrap to Benzene@short w/L then flat.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbsDBwater_ExtrapFlatBenzene.csv", "r");
  //Normalise to 10x DAYA BAY water, extrap to Benzene@short w/L then flat.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs10xDBwater_ExtrapFlatBenzene.csv", "r");
  //Normalise to 100x DAYA BAY water, extrap to Benzene@short w/L then flat.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs100xDBwater_ExtrapFlatBenzene.csv", "r");
  //Normalise to 500x DAYA BAY water, extrap to Benzene@short w/L then flat.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs500xDBwater_ExtrapFlatBenzene.csv", "r");
  //Normalise to 500x DAYA BAY water, extrap to Benzene@short w/L then flat.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs700xDBwater_ExtrapFlatBenzene.csv", "r");
  //Normalise to 1000x DAYA BAY water, extrap to Benzene@short w/L then flat.
  //FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/WbLSExtraAbs1000xDBwater_ExtrapFlatBenzene.csv", "r");
  //Daya bay pure LS absorption; units are in METRES! 
  FilePtr = fopen("/mnt/hgfs/share/ScintData/nsrl13a_wavelength_data_20140417_new/DBScintAbsorption.csv", "r");

  if(FilePtr!=0){
    //GetOptInfo(FilePtr, 1*mm);
    GetOptInfo(FilePtr, 1*m);

  }
  else{
    G4cout << "ERROR: Could not open file!" << G4endl;
  }
  G4cout << "Successfully read file!" << G4endl;
  cols = thedata.size();
  G4cout << "Number of columns = " << cols << G4endl;
  for(int i=0; i < cols; i++){
    G4cout << "Number of entries in column # " << i << " = "
	   << thedata.at(i).size() << G4endl;
  }

  const G4int size7 = thedata.at(0).size();
  G4double En7[size7];
  //std::copy(thedata.at(0).begin(), thedata.at(0).end(), En7);
  G4double Rindex7[size7];
  //std::copy(thedata.at(1).begin(), thedata.at(1).end(), Rindex7);
  for(int i = 0; i<size7; i++){
    En7[i] = thedata.at(0).at(i);
    Rindex7[i] = thedata.at(1).at(i);
    //G4cout << "RINDEX[" << i << "] = " << Rindex2[size2] << G4endl;
  }
  MPTWbLS->AddProperty("WLSABSLENGTH", En7, Rindex7, size7);

  G4MaterialPropertyVector* theWLSabs = MPTWbLS->GetProperty("WLSABSLENGTH");
  G4cout << "WLSAbs max value = " << theWLSabs->GetMaxValue() << G4endl;
  G4cout << "WLSAbs min value = " << theWLSabs->GetMinValue() << G4endl;
  G4cout << "WLSAbs at 200 nm (6.1992 eV): "
	 << theWLSabs->Value(6.1992*eV) << G4endl;
  G4cout << "WLSAbs at 165 nm (7.51419315151515 eV): "
	 << theWLSabs->Value(7.51419315151515*eV) << G4endl;
  G4cout << "WLSAbs at 700 nm (1.77120267142857 eV): "
	<< theWLSabs->Value(1.77120267142857*eV) << G4endl;
  G4cout << "MinLowEdgeEnergy = " << theWLSabs->GetMinLowEdgeEnergy() << G4endl;
  G4cout << "MaxLowEdgeEnergy = " << theWLSabs->GetMaxLowEdgeEnergy() << G4endl;


  MPTWbLS->AddConstProperty("WLSTIMECONSTANT", 1.*ns);

  WbLS->SetMaterialPropertiesTable(MPTWbLS);
  WbLS->GetIonisation()->SetBirksConstant(0.1*mm/MeV);


/*the volumes: */

  G4double worldx = 1 * m;
  G4double worldy = 1 * m;
  G4double worldz = 1 * m;

//the whole simulation must be contained within a "world" volume
//defining a volume requires definition of solid, logical, and physical volume
//solid volume is the shape, has dimensions
  G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

//logical volume: has a solid volume as a member, a material, last 3???
  G4LogicalVolume* logical_world = new G4LogicalVolume(world, air_nist, "world_log", 0,0,0);

  //make the world invisible!
  logical_world->SetVisAttributes(G4VisAttributes::Invisible);

  //physical volume: G4PVPlacement class, represents a logical volume that is placed somewhere in the mother volumeG4OpticalSurface
  G4VPhysicalVolume* physical_world = new G4PVPlacement(0,G4ThreeVector(),logical_world, "world_phys", 0, false, 0);


  //Add the first hodoscope (H2) as a rectangular prism of plastic scintillator.
  //I'm unsure of the dimensions; for now I'll assume that it is a 25 mm cube.
  G4Box* hodo = new G4Box("hodoscope", 20*mm, 20*mm, 0.5*mm);

  G4LogicalVolume* hodo_front_log =
    new G4LogicalVolume(hodo, Polystyrene_nist, "hodoFront_log", 0,0,0);

  //G4RunManager to see which sensitive detectors there are
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //create SensitiveDetector object
  HodoFront* SDHoFront = new HodoFront("hodoFront_log");
  //pass new sensitive detector to manager
  SDman->AddNewDetector(SDHoFront);
  hodo_front_log->SetSensitiveDetector(SDHoFront);

  G4VPhysicalVolume* hodoFront_phys =
    new G4PVPlacement(0, G4ThreeVector(0, 0, -508*mm), hodo_front_log,
		      "hodoFront_phys", logical_world, false, 0);

  G4LogicalVolume* hodo_rear_log =
    new G4LogicalVolume(hodo, Polystyrene_nist, "hodoRear_log", 0,0,0);
  

  //create SensitiveDetector object
  HodoRear* SDHoRear = new HodoRear("hodoRear_log");
  //pass new sensitive detector to manager
  SDman->AddNewDetector(SDHoRear);
  hodo_rear_log->SetSensitiveDetector(SDHoRear);

  G4VPhysicalVolume* hodoRear_phys = 
    new G4PVPlacement(0, G4ThreeVector(0, 0, 508*mm), hodo_rear_log,
		      "hodoRear_phys", logical_world, false, 0);

  //Create the scintillating volume box.
  //To make this I need to make a box equal in dimensions to the outer walls of
  //ABS plastic box, then make another box equal in dimensions to the inner
  //walls of the ABS plastic box. Then I need to subtract some circular holes
  //from the wall of this solid, where the UV acrylic goes.
  //from the technical drawings, the thickness of the ABS is different in
  //different locations, making my life more difficult...
  //Outer dimensions:
  G4double widthOuter = 88.9*mm;
  G4double heightOuter = 139.48*mm;// = 126.78+12.70 mm
  G4double depthOuter = 136.36*mm;// = 130 mm (measured) + 2*3.18 mm
  //Inner dimensions:
  G4double widthInner =63.5*mm;// = 88.9 - 2*12.7 mm (both walls same thickness)
  G4double heightInner = 117.25*mm;// = 139.48 - 12.7 - 7.94 - 1.59 mm
  //(bottom is 12.7 mm, top is 9.53 mm thick, offset by +3.17 mm along y).
  G4double heightInnerOffset = 3.17*mm;
  G4double depthInner = 130*mm;//(measured, both walls same thickness).
  G4double PMTHoleDiam = 57.15*mm;
  //PMT centre is offset from BOTTOM of vessel by +43.98 mm. This corresponds to
  //being offset from the centre of the vessel by -25.76 mm.
  G4double PMTHoleOffset = -25.76*mm;

  G4Box* OuterBox =
    new G4Box("OuterBox", widthOuter/2, heightOuter/2, depthOuter/2);

  G4Box* InnerBox =
    new G4Box("InnerBox", widthInner/2, heightInner/2, depthInner/2);

  G4VSolid* MTBox =
    new G4SubtractionSolid("Empty Box", OuterBox, InnerBox,
			   0, G4ThreeVector(0, heightInnerOffset, 0));

  G4Tubs* PMTHole = new G4Tubs("PMT Hole", 0, PMTHoleDiam/2, 100*cm,0,360*deg);

  G4VSolid* theBox =
    new G4SubtractionSolid("The detector box", MTBox, PMTHole,
			   0, G4ThreeVector(0, PMTHoleOffset, 0));

  G4LogicalVolume* theBox_log = new G4LogicalVolume(theBox, Polystyrene_nist,
						   "theBox_log", 0,0,0);

  G4VisAttributes* VA_theBox = new G4VisAttributes();
  VA_theBox->SetColour(1,0,1);
  theBox_log->SetVisAttributes(VA_theBox);

  //Need to set optical properties of black ABS here.

  //theBoxOffset determines the displacement from the centre of the LS box to 
  //the beam centre. The technical drawings indicate that this is 34.5 mm below
  //the top of the outside of the box. Therefore, the beam offset from the 
  //centre of the box is:
  G4double theBoxOffset = - heightOuter/2 + 34.5*mm;
  G4VPhysicalVolume* theBox_phys =
    new G4PVPlacement(0, G4ThreeVector(0, theBoxOffset, 0), theBox_log,
		    "theBox_phys", logical_world, false, 0);


  //Model the inner UV acrylic window. Diameter = PMTHoleDiam,
  //thickness on tech drawing = 0.125" = 3.175 mm. To make life easier, I'll
  //set the thickness to be the wall thickness = 3.18*mm.
  G4double InnerWinThick = 3.18*mm;
  G4Tubs* WinHole = new G4Tubs("PMT Hole", 0,
			       PMTHoleDiam/2, InnerWinThick/2, 0, 360*deg);

  G4Tubs* WinHoleCut = new G4Tubs("SubSolid", PMTHoleDiam/2, PMTHoleDiam,
  				  100*cm, 0, 360*deg);

  //Make a plane with thickness equal to window thickness.
  //G4Box* PlaneCut = new G4Box("SubSolid", 100*cm, 100*cm, InnerWinThick/2);

  //Make another WinHole that has been rotated about Y-axis by 10*deg, then
  //cut off edge bits.
  //G4RotationMatrix* RotHole = new G4RotationMatrix();
  //G4double Ang = 10.*deg;
  //RotHole->rotateY(Ang);
  //G4double RotOffset = (PMTHoleDiam/2)*tan(Ang);
  //G4SubtractionSolid* WinHolRot =
  //new G4SubtractionSolid("PMT Hole, Rotated", WinHole, WinHoleCut, RotHole,
  //			   G4ThreeVector(0,0,0));

  //Now model the UV acrylic
  G4LogicalVolume* window_log = new G4LogicalVolume(WinHole, Acrylic_nist,
  						    "window_log", 0,0,0);

  //G4LogicalVolume* winRot_log = new G4LogicalVolume(WinHolRot, Acrylic_nist,
  //						    "window_log", 0,0,0);

  //G4VisAttributes* VA_winRot = new G4VisAttributes();
  //VA_winRot->SetColour(0,1,0);
  //VA_winRot->SetForceSolid(true);
  //winRot_log->SetVisAttributes(VA_winRot);

  //Place one hole upstream and one downstream.
  //Offset relative to beam is equal to theBoxOffset + PMTHoleOffset. Offset
  //along z axis is +/- (depthOuter/2 - InnerWinThick/2)
  //(=offset - thickness/2).
  G4VPhysicalVolume* windowUp_phys =
    new G4PVPlacement(0,G4ThreeVector(0, (theBoxOffset + PMTHoleOffset), 
				      -(depthOuter/2)+(InnerWinThick/2)),
		      window_log, "windowUp_phys", logical_world, false, 0);

  G4VPhysicalVolume* windowDown_phys = 
    new G4PVPlacement(0, G4ThreeVector(0, (theBoxOffset+PMTHoleOffset),
				       (depthOuter/2)-(InnerWinThick/2)),
		      window_log, "windowDown_phys", logical_world, false, 0);
  
  //Now need to make the UV transparent glass on the windows to the tubes
  //holding the PMTs.
  //Diameter = 2.8750" = 73.025 mm, thickness = 0.1875" = 4.7625 mm.
  G4double OuterWinDiam = 73.025*mm;
  G4double OuterWinThick = 4.7625*mm;
  G4Tubs* OuterWin = new G4Tubs("Outer Window", 0,
				OuterWinDiam/2, OuterWinThick/2, 0, 360*deg);

  G4LogicalVolume* OuterWin_log = new G4LogicalVolume(OuterWin, Acrylic_nist,
						      "OuterWin_log", 0,0,0);

  //Place the outer windows. The offset along y is the same as inner window.
  //The z-axis location is +/- (depthOuter/2 + OuterWinThick/2)
  G4VPhysicalVolume* outerWinUp_phys = 
    new G4PVPlacement(0, G4ThreeVector(0, (theBoxOffset + PMTHoleOffset),
				       -(depthOuter/2)-(OuterWinThick/2)),
		      OuterWin_log, "OuterWinUp_phys", logical_world, false, 0);
  
  G4VPhysicalVolume* outerWinDown_phys = 
    new G4PVPlacement(0, G4ThreeVector(0, (theBoxOffset + PMTHoleOffset),
				       (depthOuter/2)+(OuterWinThick/2)),
		      OuterWin_log, "OuterWinDown_phys",
		      logical_world, false, 0);

  //Now to make the PMT glass. According to drawings, diam = 53 mm. Thickness is
  //indeterminate, however I'll make it 5 mm. Actually, accroding to Hamamatsu
  //datasheet, the PMT diameter is 52 mm. I've modified accordingly.
  G4double PMTglassDiam = 51*mm;
  G4double PMTglassThick = 5*mm;
  G4Tubs* PMTglass = new G4Tubs("PMT glass", 0, PMTglassDiam/2, PMTglassThick/2,
				0, 360*deg);

  G4LogicalVolume* PMTglass_log_1 = new G4LogicalVolume(PMTglass, glass_nist,
							"PMTglassUp_log",
							0,0,0);

  G4LogicalVolume* PMTglass_log_2 = new G4LogicalVolume(PMTglass, glass_nist,
							"PMTglassDown_log",
							0,0,0);

  //Create rotation matrix to tilt the PMT glass.
  //G4RotationMatrix* RM = new G4RotationMatrix();
  //RM->rotateX(0.018*rad);

  //create SensitiveDetector object
  //Upstream PMT
  PMTwinUp* SDPMTup = new PMTwinUp("PMTwinUp_log");
  //pass new sensitive detector to manager
  SDman->AddNewDetector(SDPMTup);
  PMTglass_log_1->SetSensitiveDetector(SDPMTup);

  //Downstream PMT
  PMTwinDown* SDPMTdown = new PMTwinDown("PMTwinDown_log");
  //pass new sensitive detector to manager
  SDman->AddNewDetector(SDPMTdown);
  PMTglass_log_2->SetSensitiveDetector(SDPMTdown);

  //PMT glass has same y-axis offset as acrylic window.
  //Z offset = +/- (depthOuter/2 + OuterWinThick + PMTglassThick/2)
  G4double PMTglassOffset = (depthOuter/2) + OuterWinThick + (PMTglassThick/2);

  G4VPhysicalVolume* PMTglassUp_phys =
    new G4PVPlacement(0, G4ThreeVector(0, (theBoxOffset + PMTHoleOffset),
				       -(PMTglassOffset+PMTGap)),
		      PMTglass_log_1, "PMTglassUp_phys",
		      logical_world, false, 0);

  G4VPhysicalVolume* PMTglassDown_phys = 
    new G4PVPlacement(0, G4ThreeVector(0, (theBoxOffset + PMTHoleOffset),
				       (PMTglassOffset+PMTGap)), PMTglass_log_2,
		      "PMTglassDown_phys", logical_world, false, 0);
  

  //Now to add the WbLS/water to the mix...
  //The dimensions of the liquid are identical to those of the inner box volume.
  //Water case
  //G4LogicalVolume* liquid_log = new G4LogicalVolume(InnerBox, water_nist,
  //						    "liquid_log", 0,0,0);
  //WbLS case
  G4LogicalVolume* liquid_log = new G4LogicalVolume(InnerBox, WbLS,
  						    "liquid_log", 0,0,0);

  //Need to register liquid as sensitive detector...
  SensitiveDetector* SDLiquid = new SensitiveDetector("liquid_log");
  //get pointer to the sensitive detector manager: this class is used by
  //pass to manager
  SDman->AddNewDetector(SDLiquid);

  //now we pass the sensitive detector pointer to the logical volume of our 
  //scoring volume
  liquid_log->SetSensitiveDetector(SDLiquid);


  //The liquid placement vector will be the same as the subtraction solid used
  //to create the space for it in the first place.
  G4VPhysicalVolume* liquid_phys =
    new G4PVPlacement(0, G4ThreeVector(0, theBoxOffset + heightInnerOffset, 0),
		      liquid_log, "liquid_phys", logical_world, false, 0);


  //////////Do optical surfaces...////////////
  G4OpticalSurface* OptSurf_WaterBox =
    new G4OpticalSurface("Optical surface, Water-Box", glisur, polished, 
  			 dielectric_metal);//I don't want G4 calculating
  //refractions or reflections at this boundary, so treat box as 'metal'.
  G4LogicalBorderSurface* WaterToBox =
    new G4LogicalBorderSurface("WaterToBox", liquid_phys, theBox_phys,
  			       OptSurf_WaterBox);
 

  //G4OpticalSurface* OptSurf_WaterAcrylic =
  //new G4OpticalSurface("Optical surface, Water-Acrylic");
  //OptSurf_WaterAcrylic->SetModel(unified);
  //OptSurf_WaterAcrylic->SetType(dielectric_dielectric);
  //OptSurf_WaterAcrylic->SetFinish(ground);
  //OptSurf_WaterAcrylic->SetPolish(0.1);
  //OptSurf_WaterAcrylic->SetSigmaAlpha(1.);  

  //G4LogicalBorderSurface* WaterToAcrylicUp =
  //new G4LogicalBorderSurface("WaterToAcrylicUp", liquid_phys, windowUp_phys,
  //			       OptSurf_WaterAcrylic);
  
  //G4LogicalBorderSurface* WaterToAcrylicDown =
  //new G4LogicalBorderSurface("WaterToAcrylicDown", liquid_phys, 
  //			       windowDown_phys, OptSurf_WaterAcrylic);
  
  
//G4OpticalSurface* OptSurf_AcrylicAir =
//new G4OpticalSurface("Optical surface, Acrylic-Air");
//OptSurf_AcrylicAir->SetModel(glisur);
//OptSurf_AcrylicAir->SetType(dielectric_dielectric);
//OptSurf_AcrylicAir->SetFinish(ground);
//OptSurf_AcrylicAir->SetPolish(0.5);
 
//G4LogicalBorderSurface* AcrylicUpToAir =
//new G4LogicalBorderSurface("AcrylicUpToAir", windowUp_phys, physical_world,
//			       OptSurf_AcrylicAir);
  
//G4LogicalBorderSurface* AcrylicDownToAir =
//new G4LogicalBorderSurface("AcrylicDownToAir", windowDown_phys,
//			       physical_world, OptSurf_AcrylicAir);
  
  
//G4OpticalSurface* OptSurf_AirPMT =
//new G4OpticalSurface("Optical surface, Air-PMT");
//OptSurf_AirPMT->SetModel(glisur);
//OptSurf_AirPMT->SetType(dielectric_dielectric);
//OptSurf_AirPMT->SetFinish(ground);
//OptSurf_AirPMT->SetPolish(0.5);
  
//G4LogicalBorderSurface* AirToPMTUp =
//new G4LogicalBorderSurface("AirToPMTUp", physical_world, PMTglassUp_phys,
//			       OptSurf_AirPMT);
  
//G4LogicalBorderSurface* AirToPMTDown =
//new G4LogicalBorderSurface("AirToPMTDown", physical_world,
//			       PMTglassDown_phys, OptSurf_AirPMT);
  
    
  return physical_world;
}


void DetectorConstruction::GetOptInfo(FILE* pfile, G4double unit){
  //Create the container for the data.
  vector < double > wl;
  vector < double > prop;
  int nvals = 0;
  float wavelength;
  float property;  
  //G4cout << "Entering while loop" << G4endl;
  //now loop through the file until EOF
  while(1){
    //file format is, wavelength [whitespace] property. 
    //G4cout << "Calling fscanf" << G4endl;
    nvals = fscanf(pfile, "%f %f ", &wavelength, &property);
    //G4cout << "fscanf called; nvals = " << nvals << G4endl;
    if(nvals == 2){
      //store in vectors.
      wl.push_back((1239.84187/wavelength)*eV);//convert from nm to eV
      prop.push_back(property*unit);
      //G4cout << "wavelength = " << wavelength << ", property = " << property
      //     << G4endl;
    }
    else if(nvals == -1){
      break;
    }
    else{
      printf("FILE READ ERROR.\n");
    }
  }
  //G4cout << "Finished reading file, outputting" << G4endl;
  //finished file, return the vector
  thedata.push_back(wl);
  thedata.push_back(prop);
}
