#!/bin/bash
#The default behavior will be to start with setLightYield 101 and 475 MeV.
#This is located in RunDefault.mac
cp RunDefault.mac RunSimple.mac
################
#0th run
./NSRL13Av2 RunSimple.mac > log0.txt
echo Finished 0th run
mv Edep.root /mnt/hgfs/share/beamrun2/Simulations/FinalSims_150729/WbLS_475MeV_116phMeV.root
################
################
#Ninth run
perl -pi -e 's/116./117./' RunSimple.mac
./NSRL13Av2 RunSimple.mac > log9.txt
echo Finished Ninth run!
mv Edep.root /mnt/hgfs/share/beamrun2/Simulations/FinalSims_150729/WbLS_475MeV_116phMeV.root
################
perl -pi -e 's/117./118./' RunSimple.mac
./NSRL13Av2 RunSimple.mac >> log9.txt
mv Edep.root /mnt/hgfs/share/beamrun2/Simulations/FinalSims_150729/WbLS_475MeV_118phMeV.root
################
perl -pi -e 's/118./119./' RunSimple.mac
./NSRL13Av2 RunSimple.mac >> log9.txt
mv Edep.root /mnt/hgfs/share/beamrun2/Simulations/FinalSims_150729/WbLS_475MeV_119phMeV.root
################
perl -pi -e 's/119./120./' RunSimple.mac
./NSRL13Av2 RunSimple.mac >> log9.txt
mv Edep.root /mnt/hgfs/share/beamrun2/Simulations/FinalSims_150729/WbLS_475MeV_120phMeV.root
################
##All Done!##