//Detector Messenger Class controls UI Commands for adjusting detector params
#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
//#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"


DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:Detector(Det)
{ 
  //Set up Directories
  Dir = new G4UIdirectory("/CustomCommands");
  Dir->SetGuidance("UI commands to modify this simulation");
  
  detDir = new G4UIdirectory("/CustomCommands/det/");
  detDir->SetGuidance("Detector Geometry Commands");
       
  AbsThicknessCmd = new G4UIcmdWithADoubleAndUnit("/CustomCommands/det/setAbsThickness",this);  
  AbsThicknessCmd->SetGuidance("Set HDPE Absorber thickness");
  AbsThicknessCmd->SetParameterName("Size",false);
  AbsThicknessCmd->SetRange("Size>0.");
  AbsThicknessCmd->SetUnitCategory("Length");
  AbsThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);  

  PMTGapCmd = new G4UIcmdWithADoubleAndUnit("/CustomCommands/det/setPMTGap",this);  
  PMTGapCmd->SetGuidance("Set gap between PMT and UVT Acrylic");
  PMTGapCmd->SetParameterName("Size",false);
  PMTGapCmd->SetRange("Size>=0.");
  PMTGapCmd->SetUnitCategory("Length");
  PMTGapCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);  

  UpdateCmd = new G4UIcmdWithoutParameter("/CustomCommands/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  LightYieldCmd = new G4UIcmdWithADouble("/CustomCommands/det/setLightYield",
					 this);
  LightYieldCmd->SetGuidance("Set the scintillator light yield, in ph/MeV");
  LightYieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle,
				    G4State_GeomClosed,G4State_EventProc);
      
}


DetectorMessenger::~DetectorMessenger()
{
  delete LightYieldCmd;
  delete UpdateCmd;
  delete AbsThicknessCmd;
  delete PMTGapCmd;
  delete detDir;
  delete Dir;
}


void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == UpdateCmd)
    { Detector->UpdateGeometry(); }

  if( command == AbsThicknessCmd)
    { Detector->SetAbsThickness(AbsThicknessCmd->GetNewDoubleValue(newValue));}

  if( command == PMTGapCmd)
    { Detector->SetPMTGap(PMTGapCmd->GetNewDoubleValue(newValue));}

  if( command == LightYieldCmd)
    { Detector->SetLightYield(LightYieldCmd->GetNewDoubleValue(newValue));}

}
