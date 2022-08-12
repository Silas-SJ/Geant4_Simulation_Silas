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
/// \file mpdDetectorConstruction.cc
/// \brief Implementation of the mpdDetectorConstruction class

#include "mpdDetectorConstruction.hh"
#include "mpdCalorimeterSD.hh"
#include "mpdPMTSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpBoundaryProcess.hh"
//#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* mpdDetectorConstruction::fMagFieldMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdDetectorConstruction::mpdDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdDetectorConstruction::~mpdDetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* mpdDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_CONCRETE");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* mpdDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 12;
 // G4double absoThickness = 10.*mm;
 // G4double gapThickness =  5.*mm;
 // G4double calorSizeXY  = 10.*cm;

  G4double absoThickness = 300.*mm;
  G4double gapThickness =  3300.*mm;
  calorSizeX  = 230.*m;
  calorSizeY  = 25.*m;
 
  G4double scintThickness = 20.*mm;
  G4double scintSizeXY  = 40.*cm;
    
  G4double detabsThickness = 5.*mm;
  G4double detabsSizeXY  = 40.*cm;
 
  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = fNofLayers * layerThickness;
  auto worldSizeX = 1.2 * calorSizeX;
  auto worldSizeY = 1.2 * calorSizeY;
  auto worldSizeZ  = 1.2 * calorThickness;

//Scintillator + waveguide
  G4double scinti_x = 40.*cm/2.;
  G4double scinti_y = 40.*cm/2.;
  G4double scinti_z = 2.5*cm/2.;
  G4double pDz = 30.*cm/2.;
  G4double pDy1 = 40.*cm/2.;
  G4double pDy2 = 7.5*cm/2.;
  G4double pDx1 = 2.5*cm/2.;
  G4double pDx2 = 2.5*cm/2.;
  G4double pPhi = 90.*deg;
  G4double PMT_pRMin = 0*cm;
  G4double PMT_pRMax = 2.6*cm;
  G4double vidro_Dz = 3.*mm/2.;
  G4double PMT_SPhi = 0.*deg;
  G4double PMT_DPhi = 360.*deg;
  G4double lamina2_Dz = 1.*mm/2.;
//ghost volume for the detector (scintillator + waveguide +PMT)
  G4double caixa_x = 72.*cm/2;
  G4double caixa_y = 40.1*cm/2;
  G4double caixa_z = 5.2*cm/2;
  
//================== Get materials ==================

  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
  auto absorberMaterial = G4Material::GetMaterial("G4_CONCRETE");
  auto gapMaterial = G4Material::GetMaterial("G4_AIR");
  auto scintMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto detabsMaterial = G4Material::GetMaterial("G4_Pb");

//Vidro
  G4Element* C = new G4Element("Carbono" , "C", 6 , 12.01*g/mole);
  G4Element* H = new G4Element("Hidrogenio", "H", 1 , 1.01*g/mole);
  G4Material* Al = new G4Material("Aluminum", 13.,26.98*g/mole,2.7*g/cm3);
  G4Material* Glass = new G4Material("Glass", 1.032*g/cm3,2);
  Glass->AddElement(C,91.533*perCent);
  Glass->AddElement(H,8.467*perCent);
    
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("mpdDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }
    
//============== Materials Properties ================

      const G4int nEntries = 105; //new
      
   G4double PhotonEnergy[nEntries] = {2.501*eV, 2.638*eV, 2.625*eV, 2.499*eV, 2.508*eV, 2.521*eV, 2.533*eV, 2.545*eV, 2.557*eV, 2.569*eV, 2.581*eV, 2.593*eV, 2.604*eV, 2.616*eV, 2.633*eV, 2.645*eV, 2.656*eV, 2.666*eV, 2.677*eV, 2.687*eV, 2.707*eV, 2.716*eV, 2.726*eV, 2.734*eV, 2.751*eV, 2.759*eV, 2.767*eV, 2.774*eV, 2.781*eV, 2.787*eV, 2.793*eV, 2.797*eV, 2.801*eV, 2.806*eV, 2.811*eV, 2.815*eV, 2.821*eV, 2.827*eV, 2.834*eV, 2.851*eV, 2.862*eV, 2.874*eV, 2.886*eV, 2.899*eV, 2.909*eV, 2.918*eV, 2.923*eV, 2.927*eV, 2.931*eV, 2.934*eV, 2.938*eV, 2.941*eV, 2.945*eV, 2.948*eV, 2.951*eV,2.953*eV, 2.956*eV, 2.974*eV, 2.962*eV, 2.980*eV, 2.983*eV, 2.989*eV, 2.994*eV, 2.998*eV,  3.003*eV, 3.007*eV, 3.011*eV, 3.013*eV, 3.016*eV, 3.018*eV, 3.020*eV, 3.021*eV, 3.023*eV,  3.025*eV, 3.026*eV, 3.027*eV, 3.028*eV, 3.030*eV, 3.031*eV, 3.032*eV, 3.0327*eV, 3.033*eV, 3.034*eV, 3.035*eV, 3.036*eV, 3.037*eV, 3.038*eV, 3.040*eV, 3.041*eV, 3.043*eV, 3.046*eV, 3.049*eV, 3.051*eV, 3.054*eV, 3.057*eV, 3.061*eV, 3.065*eV, 3.070*eV, 3.081*eV, 3.087*eV, 3.094*eV, 3.102*eV, 3.111*eV, 3.112*eV, 3.115*eV}; //105
   
    const G4int nEntries2 = 32;
      
    G4double PhotonEnergy2[nEntries2] = {2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV, 2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV, 2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV, 2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV, 2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV, 3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV, 3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV, 3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };//32
    
    G4double RefractiveIndex1[nEntries] ={1.333, 1.334, 1.335, 1.337, 1.338, 1.34, 1.341, 1.342, 1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,  1.3475, 1.348, 1.3485, 1.3492, 1.35, 1.3505, 1.351, 1.3518, 1.3522, 1.3530, 1.3535, 1.354, 1.3545, 1.355,  1.3555, 1.356,  1.3568, 1.3572, 1.358, 1.3585, 1.359, 1.3595, 1.36, 1.3608, 1.3609, 1.3610, 1.3611, 1.3612, 1.3613, 1.3614, 1.3615, 1.3616, 1.3617, 1.3618, 1.3619, 1.3620, 1.3621, 1.3622, 1.3623, 1.3624, 1.3625, 1.3626, 1.3627, 1.3628, 1.3629, 1.3630, 1.3631, 1.3632, 1.3633, 1.3634, 1.3635, 1.3636, 1.3637, 1.3638, 1.3639, 1.3640, 1.3641, 1.3642, 1.3643, 1.3644, 1.3645, 1.3646, 1.3647, 1.3648, 1.3649, 1.3650, 1.3651, 1.3652, 1.3653, 1.3654, 1.3655, 1.3656, 1.3657, 1.3658, 1.3659, 1.3660, 1.3661, 1.3662, 1.3663, 1.3664, 1.3665, 1.3666, 1.3667, 1.3668, 1.3669, 1.3670, 1.3671, 1.3674, 1.3675};//105
    
    G4double RefractiveIndex2[nEntries] = { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};//105
    
    G4double ScintilFast[nEntries] = { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};//105
    
  G4double ScintilSlow[nEntries] = {0.1137, 0.2203, 0.2140, 0.1116, 0.1235, 0.1307, 0.1370, 0.1442, 0.1524, 0.1611, 0.1699,0.1796, 0.1896, 0.1993, 0.2169, 0.2254, 0.2376, 0.2503, 0.2625, 0.2753, 0.2893, 0.3037, 0.3183, 0.3336, 0.3498, 0.3663, 0.3831, 0.3999, 0.4174, 0.4357, 0.4544, 0.4734, 0.4933, 0.5142, 0.5351, 0.5557, 0.5762, 0.5968, 0.6170, 0.6366, 0.6547, 0.6715, 0.6865, 0.6986, 0.7071, 0.7125, 0.7185, 0.7288, 0.7450, 0.7649, 0.7861, 0.8073, 0.8285, 0.8497, 0.8712, 0.8927, 0.9139, 0.9354, 0.9572, 0.9787, 0.9985, 0.9989, 1.001, 0.9799, 0.9607, 0.9413, 0.9206, 0.8996, 0.8789, 0.8579, 0.8364, 0.8148, 0.7929, 0.7710, 0.7491, 0.7273, 0.7054, 0.6832, 0.6610, 0.6389, 0.6167, 0.5945, 0.5723, 0.5502, 0.5280, 0.5058, 0.4837, 0.4615, 0.4393, 0.4174, 0.3956, 0.3737, 0.3518, 0.3302, 0.3086, 0.2868, 0.2652, 0.2439, 0.2226, 0.2016, 0.1812, 0.1609, 0.1405, 0.1210, 0.1024};//105 //0.08514, 0.06932, 0.06654, 0.06373}; //109
  
    
  G4double Absorption1[nEntries2] = {3.448*m, 4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m, 15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m, 45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m, 52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m, 30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m, 17.500*m, 14.500*m};//32
   
    // Air
    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
    gapMaterial->SetMaterialPropertiesTable(myMPT2);
   
    // Scintillator
    G4MaterialPropertiesTable* myMPT5 = new G4MaterialPropertiesTable();
    myMPT5->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
    myMPT5->AddProperty("ABSLENGTH",    PhotonEnergy2, Absorption1,     nEntries2);
    myMPT5->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
    myMPT5->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries);

    myMPT5->AddConstProperty("SCINTILLATIONYIELD",10./keV);
 //   myMPT5->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
    myMPT5->AddConstProperty("RESOLUTIONSCALE",1.0);
    myMPT5->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
    myMPT5->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
    myMPT5->AddConstProperty("YIELDRATIO",0.8);

    scintMaterial->SetMaterialPropertiesTable(myMPT5);
    
    G4double Glass_RIND[nEntries]= { 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52, 1.52};//105

    G4double Glass_AbsLength[nEntries]= { 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm};//105

    
    G4MaterialPropertiesTable *Glass_mt = new G4MaterialPropertiesTable();
    Glass_mt->AddProperty("ABSLENGTH",PhotonEnergy,Glass_AbsLength,nEntries);
    Glass_mt->AddProperty("RINDEX",PhotonEnergy,Glass_RIND,nEntries);

    Glass->SetMaterialPropertiesTable(Glass_mt);
    
    
//================= Geometry Construction =====================
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeX/2, calorSizeY/2, calorThickness/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
  //                                 
  // Layer
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeX/2, calorSizeY/2, layerThickness/2); //its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name
  auto layerLV2
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer2");         // its name
  for (int i=0;i<fNofLayers;i++){
      if(i==8){
          new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -calorThickness/2+layerThickness*(i+1./2)), // its position
                  layerLV2,       // its logical volume
                  "Layer2",           // its name
                  calorLV,          // its mother  volume
                  false,            // no boolean operation
                  i,                // copy number
                  fCheckOverlaps);  // checking overlaps
      }else{
          new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(0., 0., -calorThickness/2+layerThickness*(i+1./2)), // its position
                        layerLV,       // its logical volume
                        "Layer",           // its name
                        calorLV,          // its mother  volume
                        false,            // no boolean operation
                        i,                // copy number
                        fCheckOverlaps);  // checking overlaps
      }
      
  }
 
  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeX/2, calorSizeY/2, absoThickness/2); // its size
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

    new G4PVPlacement(
                     0,                // no rotation
                     G4ThreeVector(0., 0., -gapThickness/2), // its position
                     absorberLV,       // its logical volume
                     "Abso",           // its name
                     layerLV2,          // its mother  volume
                     false,            // no boolean operation
                     1,                // copy number
                     fCheckOverlaps);  // checking overlaps
  //
  // Gap
  //
  auto gapS 
    = new G4Box("Gap",             // its name
                 calorSizeX/2, calorSizeY/2, gapThickness/2); // its size
                         
  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name
  auto gapLV2
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV2");         // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
    G4VPhysicalVolume* gap2_phys =
    new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV2,            // its logical volume
                   "Gap2",            // its name
                   layerLV2,          // its mother  volume
                   false,            // no boolean operation
                   1,                // copy number. "detabsorber"
                   fCheckOverlaps);  // checking overlaps
    //
    // Scintillator
    //
    
    
//Volume "Fantasma" (caixa) para colocar o Cintilador dentro e poder mexer com ele agrupado

  G4Box* caixa = new G4Box("ghostbox",caixa_x,caixa_y,caixa_z);

  G4LogicalVolume* caixa_log = new G4LogicalVolume(caixa,defaultMaterial,"caixa",0,0,0);
  G4double detposX = 5*m-calorSizeX/2;
  G4double detposY = 2*m-calorSizeY/2;
  G4VPhysicalVolume* caixa_phys[12];
  caixa_phys[0] = new G4PVPlacement(0,G4ThreeVector(detposX+pDz+vidro_Dz, detposY,scinti_z+gapThickness/2-68*cm),caixa_log,"caixa",gapLV2,false,1); // copy 1 - bottom detector
  caixa_phys[1] = new G4PVPlacement(0,G4ThreeVector(detposX+pDz+vidro_Dz, detposY,scinti_z+gapThickness/2-detabsThickness/2-77*cm),caixa_log,"caixa",gapLV2,false,0); // copy 0 - top detector
    
//Cintilador (paralelepípedo do detector)

G4Box* Scinti_Box = new G4Box("Cintilador", scinti_x, scinti_y, scinti_z);

G4LogicalVolume* Scinti_log = new G4LogicalVolume(Scinti_Box,scintMaterial,"Scinti_log",0,0,0);

G4VPhysicalVolume* Scinti_phys = new G4PVPlacement(0,G4ThreeVector(-15.9*cm,0,0),Scinti_log,"Scinti",caixa_log,false,0);

//Cintilador (trapézio do detector)

G4RotationMatrix rm;
rm.rotateY(pPhi);

G4Trap* trapezio_Trap = new G4Trap("Trapezio",pDx1, pDx2, pDy1, pDy2, pDz);

G4LogicalVolume* trapezio_log = new G4LogicalVolume(trapezio_Trap,Glass,"trapezio_log",0,0,0);

G4VPhysicalVolume* trapezio_phys = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(35*cm-15.9*cm,0,0)),trapezio_log,"Trapezio",caixa_log,false,0);

//Vidro 2 (fotomultiplicadora)

G4RotationMatrix rm2;
rm2.rotateY(pPhi);

G4Tubs* vidro_Tubs = new G4Tubs("Foto_raquete",PMT_pRMin,PMT_pRMax,vidro_Dz,PMT_SPhi,PMT_DPhi);

G4LogicalVolume* vidro_log = new G4LogicalVolume(vidro_Tubs,Glass,"vidro",0,0,0);

G4VPhysicalVolume* vidro_phys = new G4PVPlacement(G4Transform3D(rm2,G4ThreeVector(50*cm+vidro_Dz-15.9*cm,0,0)),vidro_log,"Vidro",caixa_log,false,0);
    
//Photocatodo (fotocátodo da fotomultiplicadora, cujo volume vai "dentro" do vidro)

G4Tubs* photocathode = new G4Tubs("photocathode_tube",PMT_pRMin,PMT_pRMax,lamina2_Dz,PMT_SPhi,PMT_DPhi);

G4LogicalVolume* photocathode_log = new G4LogicalVolume(photocathode,Al,"photocathode_log",0,0,0);

G4VPhysicalVolume* photocathode_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),photocathode_log,"photocathode",vidro_log,false,0);

    
       //
       // Detector Absorber
       //
       auto detabsS
         = new G4Box("Detabs",             // its name
                      detabsSizeXY/2, detabsSizeXY/2, detabsThickness/2); // its size
                              
       auto detabsLV
         = new G4LogicalVolume(
                      detabsS,             // its solid
                      detabsMaterial,      // its material
                      "DetabsLV");         // its name
                                        
       new G4PVPlacement(
                      0,                // no rotation
                      G4ThreeVector(detposX, detposY, scintThickness/2+gapThickness/2-detabsThickness/2-72*cm), // its position
                      detabsLV,            // its logical volume
                      "Detabs",            // its name
                      gapLV2,          // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps

    //--------------------------
    // Interfaces do Cintilador
    //--------------------------

    //Interface Cintilador-Ar

      G4OpticalSurface* Op_Scinti_Ar = new G4OpticalSurface("Scinti_Ar");
      Op_Scinti_Ar -> SetType(dielectric_metal);
      Op_Scinti_Ar -> SetFinish(polished);
      Op_Scinti_Ar -> SetModel(unified);

    //  G4LogicalBorderSurface* Int_Scinti_Ar = new G4LogicalBorderSurface("Scinti_Ar",Scinti_phys,caixa_phys[1],Op_Scinti_Ar);

    //Interface Vidro(trapézio)-Ar

      G4OpticalSurface* Op_Glass_Ar = new G4OpticalSurface("Glass_Ar");
      Op_Glass_Ar -> SetType(dielectric_metal);
      Op_Glass_Ar -> SetFinish(polished);
      Op_Glass_Ar -> SetModel(unified);

     // G4LogicalBorderSurface* Int_Glass_Ar = new G4LogicalBorderSurface("Glass_Ar",trapezio_phys,caixa_phys[1],Op_Glass_Ar);

    //Interface Vidro(fotomultiplicadora)-Photocatodo

      G4OpticalSurface* photocathode_opsurf= new G4OpticalSurface("photocathode_opsurf",glisur,polished,dielectric_metal);
      G4LogicalSkinSurface* Int_vidro_photocathode = new G4LogicalSkinSurface("photocathode_surf",photocathode_log,photocathode_opsurf);

    //Interface Vidro(fotomultiplicadora)-Ar

      G4OpticalSurface* Op_Vidro_Ar = new G4OpticalSurface("Vidro_Ar");
      Op_Vidro_Ar -> SetType(dielectric_metal);
      Op_Vidro_Ar -> SetFinish(polished);
      Op_Vidro_Ar -> SetModel(unified);

     for (int i=0; i<2;i++) {
        G4LogicalBorderSurface* Int_Scinti_Ar = new G4LogicalBorderSurface("Scinti_Ar",Scinti_phys,caixa_phys[i],Op_Scinti_Ar);
        G4LogicalBorderSurface* Int_Glass_Ar = new G4LogicalBorderSurface("Glass_Ar",trapezio_phys,caixa_phys[i],Op_Glass_Ar);
        G4LogicalBorderSurface* Int_Vidro_Ar = new G4LogicalBorderSurface("Vidro_Ar",photocathode_phys,caixa_phys[i],Op_Vidro_Ar);
     }

//================================== Propriedades opticas das interfaces ===================================

    //-------------
    // Cintilador
    //-------------
    const G4int num = 2;
  
      G4double Ephoton[num] = {1.367*eV, 4.136*eV};
    //  G4double Ephoton2[num] = {2.034*eV,4.136*eV};

      G4double Reflectivity2[num] = {0.2, 0.3};
      G4double Refractive2[num] = {1.333,1.3675};
    //  G4double Efficiency2[num] = {90.,95.};
      G4double Efficiency2[num] = {0.2,0.25};


      G4MaterialPropertiesTable* scinti_ar = new G4MaterialPropertiesTable();
      scinti_ar -> AddProperty("REFLECTIVITY", Ephoton, Reflectivity2, num);
      scinti_ar -> AddProperty("REFRACTIVE", Ephoton, Refractive2, num);
      Op_Scinti_Ar -> SetMaterialPropertiesTable(scinti_ar);
      Op_Glass_Ar -> SetMaterialPropertiesTable(scinti_ar);


      G4MaterialPropertiesTable* photocathode_mt = new G4MaterialPropertiesTable();
      photocathode_mt -> AddProperty("EFFICIENCY", Ephoton, Efficiency2, num);
      photocathode_mt -> AddProperty("REFLECTIVITY", Ephoton, Reflectivity2, num);
      photocathode_opsurf -> SetMaterialPropertiesTable(photocathode_mt);

      G4MaterialPropertiesTable* vidro_ar = new G4MaterialPropertiesTable();
      vidro_ar -> AddProperty("REFLECTIVITY", Ephoton, Reflectivity2, num);
      vidro_ar -> AddProperty("REFRACTIVE", Ephoton, Refractive2, num);
      Op_Vidro_Ar -> SetMaterialPropertiesTable(vidro_ar);
      
rtiesTable(scinti_ar);
*/
  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  caixa_log->SetVisAttributes (G4VisAttributes::GetInvisible());
  auto visAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  detabsLV->SetVisAttributes(visAttributes);
  auto visAttributes2 = new G4VisAttributes(G4Colour(0.0, 0.8, 0.0));
  Scinti_log->SetVisAttributes(visAttributes2);
  auto visAttributes3 = new G4VisAttributes(G4Colour(0.8, 0.8, 0.0));
  photocathode_log->SetVisAttributes(visAttributes3);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //

  auto absoSD
    = new mpdCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto gapSD 
    = new mpdCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);

  auto scintSD
    = new mpdCalorimeterSD("ScintSD", "ScintHitsCollection", 2); // 2 scintillator
  G4SDManager::GetSDMpointer()->AddNewDetector(scintSD);
  SetSensitiveDetector("Scinti_log",scintSD);

  auto detabsSD
    = new mpdCalorimeterSD("DetabsSD", "DetabsHitsCollection", 1); // 1 detector absorber
  G4SDManager::GetSDMpointer()->AddNewDetector(detabsSD);
  SetSensitiveDetector("DetabsLV",detabsSD);
    
  auto pmtSD
    = new mpdPMTSD("mpdPMTSD");
    G4SDManager::GetSDMpointer()->AddNewDetector( pmtSD );
    SetSensitiveDetector( "vidro",pmtSD );
    
  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
