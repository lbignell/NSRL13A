/*a Geant4 applicatoin requires the use of in built GEANT4 classes, and user defined classes which inherit behaviour from GEANT4 base classes. 
Anyways, for every class that is used in the main function, we must include the defitinion which is contained in the corresponding header file: className.hh
See below for an explanation of what each class does
*/
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"//Using PhysicsFactory instead...
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
//
#include "G4VisExecutive.hh"
//
#include <time.h>
#include "Randomize.hh"
#include "WebRandomizer.hh"

//User actions
#include "RunAction.hh"

using namespace std;

/*the main function, every C++ program has one*/
int main(int argc, char** argv){
 
  //apparently has a resolution of 1 sec
  //time_t clocktime = time(NULL);
  //The following command grabs an integer random number between 1 and 10^8
  //from www.random.org
  int seed = random_int(NULL, 1);
  G4cout << "RANDOM SEED IS: " << seed << G4endl;
  if(seed==-1){
    //WebRandomizer issued an error (we're not connected, server down, etc.)
    //Use clock for seed instead
    seed = time(NULL);
  }
  //CLHEP::HepRandom::setTheSeed(seed);

  //need to make the random number seed different for each RUN

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(seed));
  G4cout << "Random engine = " << CLHEP::HepRandom::getTheEngine()->name()
	 << G4endl;
  G4cout << "Random seed (from engine) = " << CLHEP::HepRandom::getTheSeed()
	 << G4endl;


  //get the run manager pointer
  G4RunManager* rm = new G4RunManager();

  //detector construction pointer
  DetectorConstruction* detector = new DetectorConstruction();  
  //register detector
  rm->SetUserInitialization(detector);

  /////////////////Not using cusomised physics list////////////////////////
  //point to physics list
  G4VUserPhysicsList* physics = new PhysicsList();
  //register it
  rm->SetUserInitialization(physics);
  /////////////////////////////////////////////////////////////////////////

  //Using custom physics list instead of phys list factory.
  //G4PhysListFactory factory;
  //G4VModularPhysicsList* physlist =
  //factory.GetReferencePhysList("QGSP_BIC_PEN");
  //physlist->SetVerboseLevel(1);
  //rm->SetUserInitialization(physlist);

  //INITIALISE!!
  rm->Initialize();
    
  //get pointer to pga
  G4VUserPrimaryGeneratorAction* pga = new PrimaryGeneratorAction();
  //now give the run manager access to the primary generator action object
  rm->SetUserAction(pga);    

  G4UserRunAction* myRA = new RunAction(detector);

  rm->SetUserAction(myRA);


  // Get the pointer to the User Interface manager
  G4UImanager * UI = G4UImanager::GetUIpointer();
  G4VisManager* visManager = new G4VisExecutive();
  G4UIsession* session = new G4UIterminal(new G4UItcsh); 

  //if only 1 command line argument, we run in interactive mode
  if(argc == 1){

  /*now create an instance of the G4VisManager class and intialise the object*/
  visManager->Initialize();

/*this sets up the user interface to run in interactive mode */
	UI->ApplyCommand("/control/execute vis.mac");//now, we run in interactive mode so tell the UI manager to read the vis.mac macro file and 
	session->SessionStart();
	//delete session;
	//delete visManager;
    } else { 
//otherwise we run in batch mode
	G4String command = "/control/execute ";//create first part of command
	G4String fileName = argv[1];//second part is the file name that was typed at the command line 
	UI->ApplyCommand(command+fileName);//join the two and pass to the UI manager for interpretation
    }

/*the memory that was dynamically allocated for the run manager must be freed*/
 
  delete session;
  delete visManager;
  //delete rm;
  rm->~G4RunManager();

  return 0;
}
