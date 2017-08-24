// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemRPOrganization
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  
//         Created:  Tue May 16 10:14:34 CEST 2006
//

// system include files

// user include files
#include "SimG4CMS/Forward/interface/TotemRPOrganization.h"
#include "SimG4CMS/Forward/interface/TotemNumberMerger.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4Step.hh"

#include <iostream>

//******************************************************************** Constructor and destructor

TotemRPOrganization :: ~TotemRPOrganization() {
}

//
// member functions
//

uint32_t TotemRPOrganization :: GetUnitID(const G4Step* aStep) const {
  return const_cast<TotemRPOrganization *>(this)->GetUnitID(aStep);
}

uint32_t TotemRPOrganization :: GetUnitID(const G4Step* aStep) {

G4VPhysicalVolume* physVol;
  unsigned int arm = 0;
  unsigned int station = 0;
  unsigned int roman_pot = 0;
  unsigned int detector = 0;

  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();

  for(int ii = 0; ii < touch->GetHistoryDepth(); ii++ )
  {
    physVol = touch->GetVolume(ii);
    if(physVol->GetName() == "RP_Silicon_Detector")
    {
      detector = physVol->GetCopyNo();
    }
    else if(physVol->GetName() == "RP_box_primary_vacuum")
    {
      int cpy_no = physVol->GetCopyNo();
      arm = (cpy_no/100)%10;
      station = (cpy_no/10)%10;
      roman_pot = cpy_no%10;
    }

#ifdef SCRIVI
    edm::LogInfo("TotemRP") << "physVol=" << physVol->GetName()<< ", level=" << ii
	  << ", physVol->GetCopyNo()=" << physVol->GetCopyNo()<< endl;
#endif
  }

  return TotemRPDetId(arm, station, roman_pot, detector).rawId();
}
