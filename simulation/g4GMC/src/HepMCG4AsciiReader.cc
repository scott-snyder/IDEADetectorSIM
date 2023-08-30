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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMCG4AsciiReader.cc
/// \brief Implementation of the HepMCG4AsciiReader class
//
//

#include "HepMCG4AsciiReader.hh"
#include "HepMCG4AsciiReaderMessenger.hh"
#include "HepMC3/Print.h"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4AsciiReader::HepMCG4AsciiReader()
  :  filename("xxx.dat"), verbose(0)
{
  asciiInput= std::make_unique<HepMC3::ReaderAscii>(filename);

  messenger= std::make_unique<HepMCG4AsciiReaderMessenger>(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4AsciiReader::~HepMCG4AsciiReader()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4AsciiReader::Initialize()
{
  asciiInput= std::make_unique<HepMC3::ReaderAscii>(filename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::unique_ptr<HepMC3::GenEvent> HepMCG4AsciiReader::GenerateHepMCEvent()
{
  auto evt = std::make_unique<HepMC3::GenEvent>();
  if (!asciiInput->read_event(*evt)) return 0; // no more event

  if(verbose>0) HepMC3::Print::listing (*evt);

  return evt;
}