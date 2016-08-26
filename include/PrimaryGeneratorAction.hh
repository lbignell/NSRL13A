#ifndef PrimaryGeneratorAction_hh 
#define PrimaryGeneratorAction_hh 1

//as this user defined class is derived from a geant4 base class, we need to include the class definition
#include "G4VUserPrimaryGeneratorAction.hh"

//this class is used to specify details of each particle, referred to as a 'generator'
class G4ParticleGun;
class G4Event;//geant4 class which represents an event, ie. particle history

//the primary generator action class is used to specify how each primary particle is generated

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

//you must define this method, it is called by the G4RunManager
//run manager passes the pointer to an event object, it will be given attributes from the Particle Gun
    void GeneratePrimaries(G4Event*);

private:
//private member of this class, a pointer to an object of another class
    G4ParticleGun* gun;

};

#endif
