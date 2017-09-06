#ifndef SimG4Core_PrintGeomInfoAction_H
#define SimG4Core_PrintGeomInfoAction_H

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
    
#include "G4NavigationHistory.hh"

#include <iostream>
#include <vector>
#include <map>
#include <string>

// TotemRP specific START

#include "SimG4CMS/Forward/interface/TotemRPHisto.h"
#include "SimG4CMS/Forward/interface/TotemRPHistoClass.h"
#include "SimG4CMS/Forward/interface/TotemRPHistoManager.h"

#include "SimG4Core/Notification/interface/BeginOfJob.h"
#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"

// TotemRP specific END

class BeginOfJob;
class BeginOfRun;
class EndOfRun;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

typedef std::map< G4VPhysicalVolume*, G4VPhysicalVolume*, std::less<G4VPhysicalVolume*> > mpvpv;
typedef std::multimap< G4LogicalVolume*, G4VPhysicalVolume*, std::less<G4LogicalVolume*> > mmlvpv;

class PrintGeomInfoAction : public SimWatcher,
			    public Observer<const BeginOfJob *>,
//			    public Observer<const BeginOfRun *>,    -- not used in TotemRP
				public Observer<const EndOfRun *>,
				public Observer<const BeginOfEvent *>,
                public Observer<const EndOfEvent *>,
                public Observer<const G4Step *>,
                public Observer<const EndOfTrack *>,
                public Observer<const BeginOfTrack *>{

public:
  PrintGeomInfoAction(edm::ParameterSet const & p);
  ~PrintGeomInfoAction();
private:

//  void update(const BeginOfJob * job); // - duplicated in TotemRPs
  void update(const BeginOfRun * run);
  void dumpSummary(std::ostream& out = std::cout);
  void dumpG4LVList(std::ostream& out = std::cout);
  void dumpG4LVTree(std::ostream& out = std::cout);
  void dumpMaterialList(std::ostream& out = std::cout);
  void dumpG4LVLeaf(G4LogicalVolume * lv, unsigned int leafDepth, unsigned int count, std::ostream & out = std::cout);
  int countNoTouchables();
  void add1touchable(G4LogicalVolume * lv, int & nTouch);
  void dumpHierarchyTreePVLV(std::ostream& out = std::cout);
  void dumpHierarchyLeafPVLV(G4LogicalVolume * lv, unsigned int leafDepth, std::ostream & out = std::cout);
  void dumpLV(G4LogicalVolume * lv, unsigned int leafDepth, std::ostream & out = std::cout);
  void dumpPV(G4VPhysicalVolume * pv, unsigned int leafDepth, std::ostream & out = std::cout);
  void dumpTouch(G4VPhysicalVolume * pv, unsigned int leafDepth, std::ostream & out = std::cout);
  std::string spacesFromLeafDepth(unsigned int leafDepth);
  void dumpSolid(G4VSolid * sol, unsigned int leafDepth, std::ostream & out = std::cout);
  G4VPhysicalVolume * getTopPV();
  G4LogicalVolume * getTopLV();

	// TotemRP specific START
	void InitializePhysicalDetMap();

    void update(const BeginOfEvent * evt);
    void update(const EndOfEvent * evt);
    void update(const G4Step * step);
    void update(const EndOfTrack * end_of_track);
    void update(const BeginOfTrack * beg_of_track);
    void update(const BeginOfJob * job);
	void update(const EndOfRun * end_of_run);

    void PrintPrimaryVertex(G4PrimaryVertex* primaryVertex, int indent);
    void PrintParticleTreeNode(G4PrimaryParticle* particle, int indent);
    void Indent(int indent);
    bool IsGoodForTrack(G4PrimaryParticle* pp);

    // update G4Step required:
    void FillIfLeavesRP220Station(const G4Step * step);
    void FillIfParticleEntersRP(const G4Step * step);
    void FillIfParticleLeavesRP(const G4Step * step);
    void FillIfParticleLeavesFrontWallOfRP(const G4Step * aStep);
    bool IsPrimary(const G4Track * track);


    // TotemRP specific END


private:
  bool                     _dumpSummary, _dumpLVTree, _dumpLVList;
  bool                     _dumpMaterial;
  bool                     _dumpLV, _dumpSolid, _dumpAtts, _dumpSense;
  bool                     _dumpPV, _dumpRotation, _dumpReplica, _dumpTouch;
  std::string              name;
  int                      nchar;
  // std::vector<std::string> names;	- duplicated in TotemRPs
  mpvpv                    thePVTree;
  G4VPhysicalVolume *      theTopPV; 
  G4NavigationHistory      fHistory;

	// TotemRP properties START
	std::vector<std::string>            names;
	std::string                     	fileName;
	std::string 						RP_debugfileName;
	std::string 						nomeFile;
	bool verbosity_;

	TotemRPHisto * histos;
    std::auto_ptr<TotemRPHistoManager>    tuplesManager;
    TotemRPHistoClass *                   tuples;


    std::map<std::string, int> PhysicalDetMap;
	int event_no;


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
	// TotemRP properties END
};

#endif
