#ifndef Forward_TotemSD_h
#define Forward_TotemSD_h
// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemSD
//
/**\class TotemSD TotemSD.h SimG4CMS/Forward/interface/TotemSD.h
 
 Description: Stores hits of Totem in appropriate  container
 
 Usage:
    Used in sensitive detector builder 
 
*/
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
//
 
// system include files

// user include files

#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveTkDetector.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"

#include "SimG4CMS/Forward/interface/TotemG4Hit.h"
#include "SimG4CMS/Forward/interface/TotemG4HitCollection.h"
#include "SimG4CMS/Forward/interface/TotemVDetectorOrganization.h"
 
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
 
#include <string>


class G4Step;
class G4HCofThisEvent;
class TrackingSlaveSD;
class TotemTestHitHBNtuple;

class TotemSD : public SensitiveTkDetector,
		public Observer<const BeginOfEvent*>,
		public Observer<const EndOfEvent*> {

public:

  TotemSD(std::string, const DDCompactView &, const SensitiveDetectorCatalog &,
	  edm::ParameterSet const &, const SimTrackManager*);
  virtual ~TotemSD();


  void Print_Hit_Info();

  void Initialize(G4HCofThisEvent * HCE);
  virtual void EndOfEvent(G4HCofThisEvent * eventHC);
  void clear();
  void DrawAll();
  void PrintAll();

  void fillHits(edm::PSimHitContainer&, std::string use);

 private:
  virtual void clearHits();
  virtual G4bool ProcessHits(G4Step * step, G4TouchableHistory * tHistory);
  virtual uint32_t setDetUnitId(G4Step * step);
  virtual void update(const BeginOfEvent *);
  virtual void update (const ::EndOfEvent*);

  void SetNumberingScheme(TotemVDetectorOrganization* scheme);

  TrackingSlaveSD* slave;
  TotemVDetectorOrganization * numberingScheme;

private:
  int verbosity_;

  G4ThreeVector SetToLocal(G4ThreeVector globalPoint);
  void           GetStepInfo(G4Step* aStep);
  void           CreateNewHit();
  G4ThreeVector  PosizioEvo(const G4ThreeVector&,double ,double ,double, double,int&);
  void           StoreHit(TotemG4Hit*);
  void           ResetForNewPrimary();
  void           Summarize();
  bool IsPrimary(const G4Track * track);

private:

  // Data relative to primary particle (the one which triggers a shower)
  // These data are common to all Hits of a given shower.
  // One shower is made of several hits which differ by the
  // unit ID (cristal/fiber/scintillator) and the Time slice ID.

  G4ThreeVector               entrancePoint;
  double                       incidentEnergy;

  G4String name;
  G4int                       hcID;
  TotemG4HitCollection*       theHC;
  const SimTrackManager*      theManager;

  int                         tsID; 
  TotemG4Hit*                 currentHit;
  G4Track*                    theTrack;
  G4VPhysicalVolume*          currentPV;
  unsigned int                    unitID, previousUnitID;
  int                         primaryID, tSliceID;  
  double                      tSlice;

  G4StepPoint*                preStepPoint;
  G4StepPoint*                postStepPoint;
  G4ThreeVector hitPoint;
  G4ThreeVector exitPoint;
  G4ThreeVector theLocalEntryPoint;
  G4ThreeVector theLocalExitPoint;

  double Pabs;
  double p_x, p_y, p_z;
  double Tof;
  double Eloss;
  short ParticleType;

  double ThetaAtEntry;
  double PhiAtEntry;

  int ParentId;
  double Vx,Vy,Vz;


  //
  // Hist
  //
  static TotemTestHitHBNtuple* theNtuple;
  int                         eventno;
};

#endif
