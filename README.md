# Geant4_Simulation_Silas

Commands to run geant4 simulation code

-------###### To compilate ####------

Inside MPD directory:

mkdir MPD-build

sudo cmake -DGeant4_DIR={directory where geant4 is installed}/geant4.{version}-install/lib/Geant4-{version} {directory where is MPD code}

*Example:

*sudo cmake -DGeant4_DIR=/home/silas/Geant4/geant4.10.07.p03-install/lib/Geant4-10.7.3 /data/MPD

sudo make    or sudo make -jN   N is number of core


-------###### To Run ####------


./mpd   



