#ifndef SimG4CMS_Forward_TotemRP_H
#define SimG4CMS_Forward_TotemRP_H

#include "SimG4CMS/Forward/interface/TotemG4Hit.h"
#include "SimG4CMS/Forward/interface/TotemG4HitCollection.h"
#include "SimG4CMS/Forward/interface/TotemRPHisto.h"
#include "SimG4CMS/Forward/interface/TotemRPHistoClass.h"
#include "SimG4CMS/Forward/interface/TotemRPHistoManager.h"

#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfJob.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

// system include files
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <map>

class G4Step;
class BeginOfEvent;
class BeginOfJob;
class BeginOfTrack;
class EndOfEvent;
class EndOfTrack;

class TotemRP : public SimWatcher,
    public Observer<const BeginOfJob *>,
    public Observer<const BeginOfEvent *>, 
    public Observer<const EndOfEvent *>, 
    public Observer<const G4Step *>,
    public Observer<const EndOfTrack *>, 
    public Observer<const BeginOfTrack *>
{
 public:
  TotemRP(const edm::ParameterSet &p);
  virtual ~TotemRP();
  
 private:
  void update(const BeginOfEvent * evt);
  void update(const EndOfEvent * evt);
  void update(const G4Step * step);
  void update(const EndOfTrack * end_of_track);
  void update(const BeginOfTrack * beg_of_track);
  void update(const BeginOfJob * job);
  
  inline bool IsThroughPlane(double z1, double z2, double z0);
  void FillIfLeavesRP220Station(const G4Step * step);
  void KillParticlesBetweenStations(const G4Step * step);
  void FillIfParticleEntersRP(const G4Step * step);
  void FillIfParticleLeavesRP(const G4Step * step);
  void FillIfParticleLeavesFrontWallOfRP(const G4Step * aStep);
  void InitializePhysicalDetMap();
  void Indent(int indent);
  void PrintParticleTreeNode(G4PrimaryParticle* particle, int indent);
  void PrintPrimaryVertex(G4PrimaryVertex* primaryVertex, int indent);
  bool IsGoodForTrack(G4PrimaryParticle* pp);
  bool IsPrimary(const G4Track * track);


 private:
  const int primary_proton_id_code;
  const int particle_leaving_220_right_station_id_code;
  const int particle_leaving_220_left_station_id_code;
  const int particle_leaving_147_right_station_id_code;
  const int particle_leaving_147_left_station_id_code;
  const int particle_entering_RP_id_code;
  const int particle_leaving_RP_id_code;
  const int particle_leaving_front_wall_of_RP_id_code;
  const int primary_proton_inelastic_event_in_RP_station;
  const int particle_entering_station_id_code;
  
  TotemRPHisto * histos;
  int event_no;
  std::map<std::string, int> PhysicalDetMap;
  
  //Keep parameters to instantiate TotemRPHistoManager later
  std::string                             fileName;
  std::vector<std::string>                names;

  // Private Tuples
  std::string RP_debugfileName;
  std::auto_ptr<TotemRPHistoManager>    tuplesManager;
  TotemRPHistoClass *                   tuples;

  std::string nomeFile;
  bool verbosity_;
};

#endif  //TotemRP_TotemRP_H
