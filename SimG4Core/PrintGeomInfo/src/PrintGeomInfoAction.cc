#include "SimG4Core/PrintGeomInfo/interface/PrintGeomInfoAction.h"

//#include "SimG4Core/Notification/interface/BeginOfJob.h" // duplicated in TotemRPs
//#include "SimG4Core/Notification/interface/BeginOfRun.h" // duplicated in TotemRPs

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDLogicalPart.h"
#include "DetectorDescription/Core/interface/DDSplit.h"
#include "DetectorDescription/Core/interface/DDValue.h"

#include "G4Run.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4Track.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4TransportationManager.hh"

// TotemRP specific START

#include "SimG4CMS/Forward/interface/TotemG4HitCollection.h"

#include "SimG4Core/Notification/interface/BeginOfJob.h"
#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"

#include "G4SDManager.hh"
#include "G4ParticleTable.hh"

// for EndOfEvent
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"

// G4Step
#include "SimG4Core/Notification/interface/TrackInformation.h"

// BeginOfTrack
#include "G4VProcess.hh"


// TotemRP specific END

#include <set>
#include <map>

using namespace CLHEP;

PrintGeomInfoAction::PrintGeomInfoAction(const edm::ParameterSet &p)
// TotemRP specific START
:primary_proton_id_code(-1),
particle_leaving_220_right_station_id_code(-2),
particle_leaving_220_left_station_id_code(-6),
particle_leaving_147_right_station_id_code(-7),
particle_leaving_147_left_station_id_code(-8),
particle_entering_RP_id_code(-3),
particle_leaving_RP_id_code(-4),
particle_leaving_front_wall_of_RP_id_code(-9),
primary_proton_inelastic_event_in_RP_station(-5),
particle_entering_station_id_code(-10)
// TotemRP specific END
{

  // JUST TAKING ARGUMENTS FROM CONFIG FILE
  _dumpSummary = p.getUntrackedParameter<bool>("DumpSummary", true);
  _dumpLVTree  = p.getUntrackedParameter<bool>("DumpLVTree",  true);
  _dumpMaterial= p.getUntrackedParameter<bool>("DumpMaterial",false);
  _dumpLVList  = p.getUntrackedParameter<bool>("DumpLVList",  false);
  _dumpLV      = p.getUntrackedParameter<bool>("DumpLV",      false);
  _dumpSolid   = p.getUntrackedParameter<bool>("DumpSolid",   false);
  _dumpAtts    = p.getUntrackedParameter<bool>("DumpAttributes", false);
  _dumpPV      = p.getUntrackedParameter<bool>("DumpPV",      false);
  _dumpRotation= p.getUntrackedParameter<bool>("DumpRotation",false);
  _dumpReplica = p.getUntrackedParameter<bool>("DumpReplica", false);
  _dumpTouch   = p.getUntrackedParameter<bool>("DumpTouch",   false);
  _dumpSense   = p.getUntrackedParameter<bool>("DumpSense",   false);
  name  = p.getUntrackedParameter<std::string>("Name","*");
  nchar = name.find("*");
  name.assign(name,0,nchar);
  // names = p.getUntrackedParameter<std::vector<std::string> >("Names");

    // TotemRP specific START

    edm::ParameterSet m_Anal = p.getParameter<edm::ParameterSet>("TotemRP");
    verbosity_ = m_Anal.getParameter<bool>("Verbosity");
    fileName = m_Anal.getParameter<std::string>("FileName");
    nomeFile = m_Anal.getParameter<std::string>("FileNameOLD");
    names = m_Anal.getParameter<std::vector<std::string> >("Names");

    G4cout << "TotemRP specific START "
           << "\tfileName: " << fileName << "\n"
           << "\tnomeFile: " << nomeFile << "\n"
           << "\tverbosity: " << verbosity_ << "\n"
           << "\tnames: \n";
    for (unsigned int i=0; i<names.size(); i++) {
        G4cout << " " << names[i];
    }
    G4cout << G4endl;

    histos = new TotemRPHisto(nomeFile);
    event_no = -1;
    InitializePhysicalDetMap();


    // TotemRP specific END


    // PRINTING THESE ARGUMENTS...
  G4cout << "PrintGeomInfoAction:: <HELLO WE SHOULD SEE IT> initialised with verbosity levels:"
	 << " Summary   " << _dumpSummary << " LVTree   " << _dumpLVTree
	 << " LVList    " << _dumpLVList  << " Material " << _dumpMaterial
	 << "\n                                                        "
	 << " LV        " << _dumpLV      << " Solid    " << _dumpSolid
	 << " Attribs   " << _dumpAtts
	 << "\n                                                        "
	 << " PV        " << _dumpPV      << " Rotation " << _dumpRotation
	 << " Replica   " << _dumpReplica
	 << "\n                                                        "
	 << " Touchable " << _dumpTouch << " for names (0-" << nchar
	 << ") = " << name
	 << "\n                                                        "
	 << " Sensitive " << _dumpSense << " for " << names.size()
	 << " namess";

}


// TotemRP specific START
void PrintGeomInfoAction::InitializePhysicalDetMap()
{
    G4cout << "DUPA :) JESTEM W INITALIZEPHYSICALDETMAP" << G4endl;
    PhysicalDetMap["Other"] = 0;
    PhysicalDetMap["RP_Silicon_Detector"] = 1;
    PhysicalDetMap["RP_Separ_Spacer"] = 2;
    PhysicalDetMap["RP_Separ_Frame_5"] = 3;
    PhysicalDetMap["RP_PCB"] = 4;
    PhysicalDetMap["RP_Left_Right_Wall"] = 5;
    PhysicalDetMap["RP_front_wall_6"] = 6;
    PhysicalDetMap["RP_Front_Frame_3"] = 7;
    PhysicalDetMap["RP_Device_Vert_Corp_3"] = 8;
    PhysicalDetMap["RP_Device_Hor_Corp_3"] = 9;
    PhysicalDetMap["RP_box_secondary_vacuum"] = 10;
    PhysicalDetMap["RP_box_primary_vacuum"] = 11;
    PhysicalDetMap["RP_bottom_wall_6"] = 12;
    PhysicalDetMap["RP_bottom_foil"] = 13;
    PhysicalDetMap["RP_220_Right_Station_Vacuum_5"] = 14;
    PhysicalDetMap["RP_220_Right_Station_Tube_5"] = 15;
    PhysicalDetMap["RP_220_Right_Station_Tube_4"] = 16;
    PhysicalDetMap["RP_220_Right_Station_Tube_3"] = 17;
    PhysicalDetMap["RP_220_Right_Station_Tube_2"] = 18;
    PhysicalDetMap["RP_220_Right_Station_Tube_1"] = 19;
    PhysicalDetMap["RP_220_Left_Station_Vacuum_5"] = 20;
    PhysicalDetMap["RP_220_Left_Station_Tube_5"] = 21;
    PhysicalDetMap["RP_220_Left_Station_Tube_4"] = 22;
    PhysicalDetMap["RP_220_Left_Station_Tube_3"] = 23;
    PhysicalDetMap["RP_220_Left_Station_Tube_2"] = 24;
    PhysicalDetMap["RP_220_Left_Station_Tube_1"] = 25;
    PhysicalDetMap["RP_147_Right_Station_Vacuum_5"] = 26;
    PhysicalDetMap["RP_147_Right_Station_Tube_5"] = 27;
    PhysicalDetMap["RP_147_Right_Station_Tube_4"] = 28;
    PhysicalDetMap["RP_147_Right_Station_Tube_3"] = 29;
    PhysicalDetMap["RP_147_Right_Station_Tube_2"] = 30;
    PhysicalDetMap["RP_147_Right_Station_Tube_1"] = 31;
    PhysicalDetMap["RP_147_Left_Station_Vacuum_5"] = 32;
    PhysicalDetMap["RP_147_Left_Station_Tube_5"] = 33;
    PhysicalDetMap["RP_147_Left_Station_Tube_4"] = 34;
    PhysicalDetMap["RP_147_Left_Station_Tube_3"] = 35;
    PhysicalDetMap["RP_147_Left_Station_Tube_2"] = 36;
    PhysicalDetMap["RP_147_Left_Station_Tube_1"] = 37;
}

void PrintGeomInfoAction::Indent(int indent)
{
    for(int i=0; i<indent; i++)
    {
        edm::LogInfo("TotemRP")<<" ";
    }
}

bool PrintGeomInfoAction::IsGoodForTrack(G4PrimaryParticle* pp)
{
    G4ParticleDefinition* pd = pp->GetG4code();

    if(!pd)
    {
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        pd = particleTable->FindParticle(pp->GetPDGcode());
    }

    if(!pd)
    {
        return false;
    }
    else if(!(pd->IsShortLived()))
    {
        return true;
    }
// Following two lines should be removed if the user does not want to make shortlived
// primary particle with proper decay table to be converted into a track.
    else if(pd->GetDecayTable())
    {
        return true;
    }
    return false;
}

void PrintGeomInfoAction::PrintParticleTreeNode(G4PrimaryParticle* particle, int indent)
{
    while(particle!=0)
    {
        Indent(indent);

        if(particle->GetG4code())
            edm::LogInfo("TotemRP")<<particle->GetG4code()->GetParticleName();

        edm::LogInfo("TotemRP")<<" "<<particle->GetPDGcode()
                               <<" "<<particle->GetMomentum()<<" "
                               <<particle->GetProperTime()<<", "<<IsGoodForTrack(particle)<<std::endl;
        PrintParticleTreeNode(particle->GetDaughter(), indent+2);
        particle = particle->GetNext();
    }
}

void PrintGeomInfoAction::PrintPrimaryVertex(G4PrimaryVertex* primaryVertex, int indent)
{
    if(primaryVertex!=0)
    {
        G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
        Indent(indent);
        edm::LogInfo("TotemRP")<<"Primary vertex: ("<<primaryVertex->GetX0()<<", "<<primaryVertex->GetY0()
                               <<", "<<primaryVertex->GetZ0()<<"), "<<primaryVertex->GetT0()<<", "
                               <<primaryVertex->GetWeight()<<std::endl;
        PrintParticleTreeNode(primaryParticle, indent+2);
    }
}


bool PrintGeomInfoAction::IsPrimary(const G4Track * track)
{
    TrackInformation* info
            = dynamic_cast<TrackInformation*>( track->GetUserInformation() );
    return info && info->isPrimary();
}

//=================================================================== job

void PrintGeomInfoAction::update(const BeginOfJob * job){
    G4cout << "Hello from: PrintGeomInfoAction.cc BeginOfJob\n";

    // Ntuples
    tuplesManager.reset(new TotemRPHistoManager(fileName));

}

//=================================================================== event

void PrintGeomInfoAction::update(const BeginOfEvent * evt){
    G4cout << "Hello from: PrintGeomInfoAction.cc BeginOfEvent\n";

    // create tuple object
    tuples = new TotemRPHistoClass();
    //debug_event = new RPDebugEvent();

    int iev = (*evt)()->GetEventID();

    if(verbosity_)
    {
        LogDebug("TotemRP") << "TotemRP: Begin of event = " << iev;
        edm::LogInfo("TotemRP") << " Begin event !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! event " << (*evt)()->GetEventID() << std::endl;
    }

    // access to the G4 hit collections
    G4HCofThisEvent* allHC = (*evt)()->GetHCofThisEvent();
    int HCRPid = G4SDManager::GetSDMpointer()->GetCollectionID("TotemHitsRP");
    TotemG4HitCollection* theHC_RP = (TotemG4HitCollection*) allHC->GetHC(HCRPid);
    LogDebug("ForwardSim") << "TotemTestGem :: Hit Collection for TotemHitsRP of ID " << HCRPid << " is obtained at " << theHC_RP;
    G4int nentriesRP = 0;
    if (HCRPid >= 0 && theHC_RP > 0) {
        nentriesRP = theHC_RP->entries();
        if(verbosity_)
        {
            edm::LogInfo("TotemRP")<<nentriesRP<<"  entries in the collection"<<std::endl;
        }
    }

    int PrimNum = (*evt)()->GetNumberOfPrimaryVertex();
    int evtnum = (*evt)()->GetEventID();
    event_no = evtnum;

    for(int i=0; i<PrimNum; i++)
    {
        G4PrimaryVertex* PrimVert = (*evt)()->GetPrimaryVertex(i);
//    edm::LogInfo("TotemRP")<<std::endl<<"Vx="<<PrimVert->GetX0()<<"  Vy="<<PrimVert->GetY0()<<"  Vz="<<PrimVert->GetZ0()<<std::endl;
        int num_of_particles = PrimVert->GetNumberOfParticle();
        for(int j=0; j<num_of_particles; j++)
        {
            G4PrimaryParticle* primary = PrimVert->GetPrimary(j);
            if(primary->GetPDGcode()==2212)
            {
                G4cout << "Hello from loop in BeginOfEvent - filling histos\n";

                histos->set_EVT(evtnum);
                histos->set_UID(primary_proton_id_code);   //for primary particles
                histos->set_Ptype(primary->GetPDGcode());
                histos->set_TID(0);
                histos->set_PID(0);
                histos->set_ELoss(0);
                histos->set_PABS(primary->GetMomentum().mag()/GeV);
                histos->set_p_x(primary->GetPx()/GeV);
                histos->set_p_y(primary->GetPy()/GeV);
                histos->set_p_z(primary->GetPz()/GeV);
                histos->set_VX(PrimVert->GetX0()/mm);
                histos->set_VY(PrimVert->GetY0()/mm);
                histos->set_VZ(PrimVert->GetZ0()/mm);

                histos->set_X(0);
                histos->set_Y(0);
                histos->set_Z(0);
                histos->set_Loc_X(0);
                histos->set_Loc_Y(0);
                histos->set_Loc_Z(0);
                histos->set_X_Exit(0);
                histos->set_Y_Exit(0);
                histos->set_Z_Exit(0);
                histos->set_Loc_X_Exit(0);
                histos->set_Loc_Y_Exit(0);
                histos->set_Loc_Z_Exit(0);

                histos->fillNtuple();
            }
        }
    }

    for(int i=0; i<PrimNum; i++)
    {
        G4PrimaryVertex* PrimVert = (*evt)()->GetPrimaryVertex(i);
        if(verbosity_)
            PrintPrimaryVertex(PrimVert, 0);
    }


}

void PrintGeomInfoAction::update(const EndOfEvent * evt){
    G4cout << "Hello from: PrintGeomInfoAction.cc EndOfEvent\n";


    if(verbosity_)
        edm::LogInfo("TotemRP") << " Fill event " << (*evt)()->GetEventID() << std::endl;

    // access to the G4 hit collections
    G4HCofThisEvent* allHC = (*evt)()->GetHCofThisEvent();

    int ihit = 0;
    for (unsigned int in=0; in<names.size(); in++)
    {
        int HCRPid = G4SDManager::GetSDMpointer()->GetCollectionID(names[in]);
        if(HCRPid == -1)
        {


        }
        else
        {
            TotemG4HitCollection* theHC_RP = (TotemG4HitCollection*) allHC->GetHC(HCRPid);
            if(verbosity_)
                LogDebug("TotemRP") << "TotemRP :: Hit Collection for " <<names[in]
                                    << " of ID " << HCRPid << " is obtained at " << theHC_RP;

            G4int nentriesRP = 0;
            nentriesRP = theHC_RP->entries();

            if(HCRPid >= 0 && theHC_RP > 0 && nentriesRP>0)
            {
                for(ihit = 0; ihit <nentriesRP; ihit++)
                {
                    TotemG4Hit* aHit = (*theHC_RP)[ihit];

                    int evtnum = (*evt)()->GetEventID();
                    //edm::LogInfo("TotemRP")<<"event no: "<<evtnum<<std::endl;

                    int UID = aHit->getUnitID();
                    int Ptype = aHit->getParticleType();
                    int TID = aHit->getTrackID();
                    int PID = aHit->getParentId();
                    double ELoss =  aHit->getEnergyLoss();
                    double PABS =  aHit->getPabs();
                    double x = aHit->getEntry().x();
                    double y = aHit->getEntry().y();
                    double z = aHit->getEntry().z();
                    double lx = aHit->getLocalEntry().x();
                    double ly = aHit->getLocalEntry().y();
                    double lz = aHit->getLocalEntry().z();

                    double x_ex = aHit->getExit().x();
                    double y_ex = aHit->getExit().y();
                    double z_ex = aHit->getExit().z();
                    double lx_ex = aHit->getLocalExit().x();
                    double ly_ex = aHit->getLocalExit().y();
                    double lz_ex = aHit->getLocalExit().z();

                    double vx = aHit->getVx();
                    double vy = aHit->getVy();
                    double vz = aHit->getVz();

                    double p_x = aHit->get_p_x();
                    double p_y = aHit->get_p_y();
                    double p_z = aHit->get_p_z();

                    int prim_vert_id = -1;

                    if(PID!=0)
                    {
                        G4ThreeVector myPoint(vx, vy, vz);
                        G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
                        G4VPhysicalVolume* myVolume = theNavigator->LocateGlobalPointAndSetup(myPoint);
                        prim_vert_id = PhysicalDetMap[myVolume->GetName()];
                    }


                    G4cout << "Hello from loop in EndOfEvent - filling histos\n";


                    histos->set_EVT(evtnum);

                    TotemRPDetId det_id((uint32_t)UID);
                    histos->set_UID(det_id.plane()+det_id.getRPDecimalId());

                    histos->set_Ptype(Ptype) ;
                    histos->set_TID(TID) ;
                    histos->set_PID(PID);
                    histos->set_ELoss(ELoss) ;
                    histos->set_PABS(PABS) ;
                    histos->set_VX(vx);
                    histos->set_VY(vy) ;
                    histos->set_VZ(vz);

                    histos->set_X(x);
                    histos->set_Y(y);
                    histos->set_Z(z);
                    histos->set_Loc_X(lx);
                    histos->set_Loc_Y(ly);
                    histos->set_Loc_Z(lz);
                    histos->set_X_Exit(x_ex);
                    histos->set_Y_Exit(y_ex);
                    histos->set_Z_Exit(z_ex);
                    histos->set_Loc_X_Exit(lx_ex);
                    histos->set_Loc_Y_Exit(ly_ex);
                    histos->set_Loc_Z_Exit(lz_ex);
                    histos->set_p_x(p_x);
                    histos->set_p_y(p_y);
                    histos->set_p_z(p_z);
                    histos->set_prim_ver_id(prim_vert_id);
                    histos->fillNtuple();

                    tuples->fillHit(UID, Ptype, TID, PID, ELoss, PABS, p_x, p_y, p_z, vx, vy, vz,
                                    x, y, z, lx, ly, lz, x_ex, y_ex, z_ex, lx_ex, ly_ex, lz_ex, prim_vert_id);

                }
            }
        }
    }

    // tu juz sa tuple do tuplesManager
    tuplesManager->fillTree(tuples);

    if(verbosity_)
        LogDebug("TotemRP") << "TotemRP:: --- after fillTree";

}

//=================================================================== each STEP

void PrintGeomInfoAction::FillIfLeavesRP220Station(const G4Step * aStep)
{
    const G4StepPoint* thePreStepPoint = aStep->GetPreStepPoint();
    const G4StepPoint* thePostStepPoint = aStep->GetPostStepPoint();
    const G4ThreeVector & hitPoint = thePreStepPoint->GetPosition();
    double x = hitPoint.x();
    double y = hitPoint.y();
    double z = hitPoint.z();

    const double pipe_inner_radius = 80*mm; //mm

    if( thePostStepPoint && thePostStepPoint->GetPhysicalVolume()
        && aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()== "DDDWorld"
        && ( (x*x+y*y)<=pipe_inner_radius*pipe_inner_radius) )
    {
        if(thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_220_Right_Station_Vacuum_5")
        {
            //edm::LogInfo("TotemRP")<<"!!read"<<std::endl;
            //std::cout<<"220 right station left"<<std::endl;
            const G4Track *theTrack = aStep->GetTrack();

            G4cout << "Hello from G4Step - FillIfLeavesRP220Station - filling\n";

            histos->set_EVT(event_no);
            histos->set_UID(particle_leaving_220_right_station_id_code);   //for particles going from RP station 220m
            histos->set_Ptype( theTrack->GetDefinition()->GetPDGEncoding() );
            if(IsPrimary(theTrack))
                histos->set_PID(0);
            else histos->set_PID(theTrack->GetParentID());
            histos->set_ELoss(0);
            histos->set_PABS(thePreStepPoint->GetMomentum().mag()/GeV);
            histos->set_p_x(thePreStepPoint->GetMomentum().x()/GeV);
            histos->set_p_y(thePreStepPoint->GetMomentum().y()/GeV);
            histos->set_p_z(thePreStepPoint->GetMomentum().z()/GeV);
            histos->set_VX(theTrack->GetVertexPosition().x()/mm);
            histos->set_VY(theTrack->GetVertexPosition().y()/mm);
            histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

            histos->set_X(x);
            histos->set_Y(y);
            histos->set_Z(z);
            histos->set_Loc_X(0);
            histos->set_Loc_Y(0);
            histos->set_Loc_Z(0);
            histos->set_X_Exit(0);
            histos->set_Y_Exit(0);
            histos->set_Z_Exit(0);
            histos->set_Loc_X_Exit(0);
            histos->set_Loc_Y_Exit(0);
            histos->set_Loc_Z_Exit(0);

            histos->fillNtuple();
            //    edm::LogInfo("TotemRP")<<"leaving!!!! "<<theTrack->GetDefinition()->GetParticleName()<<std::endl;
        }
        else if(thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_220_Left_Station_Vacuum_5")
        {
            //edm::LogInfo("TotemRP")<<"!!read"<<std::endl;
            //std::cout<<"220 left station left"<<std::endl;
            const G4Track *theTrack = aStep->GetTrack();

            histos->set_EVT(event_no);
            histos->set_UID(particle_leaving_220_left_station_id_code);   //for particles going from RP station 220m
            histos->set_Ptype( theTrack->GetDefinition()->GetPDGEncoding() );
            if(IsPrimary(theTrack))
                histos->set_PID(0);
            else histos->set_PID(theTrack->GetParentID());
            histos->set_ELoss(0);
            histos->set_PABS(thePreStepPoint->GetMomentum().mag()/GeV);
            histos->set_p_x(thePreStepPoint->GetMomentum().x()/GeV);
            histos->set_p_y(thePreStepPoint->GetMomentum().y()/GeV);
            histos->set_p_z(thePreStepPoint->GetMomentum().z()/GeV);
            histos->set_VX(theTrack->GetVertexPosition().x()/mm);
            histos->set_VY(theTrack->GetVertexPosition().y()/mm);
            histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

            histos->set_X(x);
            histos->set_Y(y);
            histos->set_Z(z);
            histos->set_Loc_X(0);
            histos->set_Loc_Y(0);
            histos->set_Loc_Z(0);
            histos->set_X_Exit(0);
            histos->set_Y_Exit(0);
            histos->set_Z_Exit(0);
            histos->set_Loc_X_Exit(0);
            histos->set_Loc_Y_Exit(0);
            histos->set_Loc_Z_Exit(0);

            histos->fillNtuple();
            //    edm::LogInfo("TotemRP")<<"leaving!!!! "<<theTrack->GetDefinition()->GetParticleName()<<std::endl;
        }
    }
    else if( thePostStepPoint && thePostStepPoint->GetPhysicalVolume()
             && aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()== "Beam_150_220_L"
             && ( (x*x+y*y)<=pipe_inner_radius*pipe_inner_radius) &&
             thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_147_Left_Station_Vacuum_5")
    {
        //edm::LogInfo("TotemRP")<<"!!read"<<std::endl;
        //std::cout<<"150 left station left"<<std::endl;
        const G4Track *theTrack = aStep->GetTrack();

        histos->set_EVT(event_no);
        histos->set_UID(particle_leaving_147_left_station_id_code);   //for particles going from RP station 220m
        histos->set_Ptype( theTrack->GetDefinition()->GetPDGEncoding() );
        if(IsPrimary(theTrack))
            histos->set_PID(0);
        else histos->set_PID(theTrack->GetParentID());
        histos->set_ELoss(0);
        histos->set_PABS(thePreStepPoint->GetMomentum().mag()/GeV);
        histos->set_p_x(thePreStepPoint->GetMomentum().x()/GeV);
        histos->set_p_y(thePreStepPoint->GetMomentum().y()/GeV);
        histos->set_p_z(thePreStepPoint->GetMomentum().z()/GeV);
        histos->set_VX(theTrack->GetVertexPosition().x()/mm);
        histos->set_VY(theTrack->GetVertexPosition().y()/mm);
        histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

        histos->set_X(x);
        histos->set_Y(y);
        histos->set_Z(z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        histos->fillNtuple();
//    edm::LogInfo("TotemRP")<<"leaving!!!! "<<theTrack->GetDefinition()->GetParticleName()<<std::endl;
    }
    else if( thePostStepPoint && thePostStepPoint->GetPhysicalVolume()
             && aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()== "Beam_150_220_R"
             && ( (x*x+y*y)<=pipe_inner_radius*pipe_inner_radius) &&
             thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_147_Right_Station_Vacuum_5")
    {
        //edm::LogInfo("TotemRP")<<"!!read"<<std::endl;
        //std::cout<<"150 right station left"<<std::endl;
        const G4Track *theTrack = aStep->GetTrack();

        histos->set_EVT(event_no);
        histos->set_UID(particle_leaving_147_right_station_id_code);   //for particles going from RP station 220m
        histos->set_Ptype( theTrack->GetDefinition()->GetPDGEncoding() );
        if(IsPrimary(theTrack))
            histos->set_PID(0);
        else histos->set_PID(theTrack->GetParentID());
        histos->set_ELoss(0);
        histos->set_PABS(thePreStepPoint->GetMomentum().mag()/GeV);
        histos->set_p_x(thePreStepPoint->GetMomentum().x()/GeV);
        histos->set_p_y(thePreStepPoint->GetMomentum().y()/GeV);
        histos->set_p_z(thePreStepPoint->GetMomentum().z()/GeV);
        histos->set_VX(theTrack->GetVertexPosition().x()/mm);
        histos->set_VY(theTrack->GetVertexPosition().y()/mm);
        histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

        histos->set_X(x);
        histos->set_Y(y);
        histos->set_Z(z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        histos->fillNtuple();
//    edm::LogInfo("TotemRP")<<"leaving!!!! "<<theTrack->GetDefinition()->GetParticleName()<<std::endl;
    }
}


void PrintGeomInfoAction::FillIfParticleEntersRP(const G4Step * aStep)
{
    const G4StepPoint* thePreStepPoint = aStep->GetPreStepPoint();
    const G4StepPoint* thePostStepPoint = aStep->GetPostStepPoint();

    if( thePostStepPoint && thePostStepPoint->GetPhysicalVolume()
        && thePostStepPoint->GetPhysicalVolume()->GetName()== "RP_box_primary_vacuum"
        && (thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_220_Right_Station_Vacuum_5"
            || thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_220_Left_Station_Vacuum_5"
            || thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_147_Right_Station_Vacuum_5"
            || thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_147_Left_Station_Vacuum_5"))
    {
        const G4Track *theTrack = aStep->GetTrack();
        const G4ThreeVector & hitPoint = thePostStepPoint->GetPosition();
        double x = hitPoint.x();
        double y = hitPoint.y();
        double z = hitPoint.z();
        double rp_id = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();

        G4cout << "Hello from G4Step - FillIfParticleEntersRP - filling\n";


        histos->set_EVT(event_no);
        histos->set_UID(particle_entering_RP_id_code);   //for particles entering RP
        histos->set_Ptype(theTrack->GetDefinition()->GetPDGEncoding());
        if(IsPrimary(theTrack))
            histos->set_PID(0);
        else histos->set_PID(theTrack->GetParentID());
        histos->set_ELoss(rp_id);  //rp_id in eloss
        histos->set_PABS(thePostStepPoint->GetMomentum().mag()/GeV);
        histos->set_p_x(thePostStepPoint->GetMomentum().x()/GeV);
        histos->set_p_y(thePostStepPoint->GetMomentum().y()/GeV);
        histos->set_p_z(thePostStepPoint->GetMomentum().z()/GeV);
        histos->set_VX(theTrack->GetVertexPosition().x()/mm);
        histos->set_VY(theTrack->GetVertexPosition().y()/mm);
        histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

        histos->set_X(x);
        histos->set_Y(y);
        histos->set_Z(z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        histos->fillNtuple();
//    edm::LogInfo("TotemRP")<<"entering RP!!!! "<<theTrack->GetDefinition()->GetParticleName()<<" "<<rp_id<<std::endl;
    }
}


void PrintGeomInfoAction::FillIfParticleLeavesFrontWallOfRP(const G4Step * aStep)
{
    const G4StepPoint* thePreStepPoint = aStep->GetPreStepPoint();
    const G4StepPoint* thePostStepPoint = aStep->GetPostStepPoint();

    if( thePostStepPoint && thePostStepPoint->GetPhysicalVolume()
        && thePostStepPoint->GetPhysicalVolume()->GetName()== "RP_box_secondary_vacuum"
        && thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_front_wall_6"
        && thePreStepPoint->GetPhysicalVolume()->GetCopyNo()==0)
    {
        const G4Track *theTrack = aStep->GetTrack();
        const G4ThreeVector & hitPoint = thePostStepPoint->GetPosition();
        double x = hitPoint.x();
        double y = hitPoint.y();
        double z = hitPoint.z();

        const G4VTouchable* touchable = theTrack->GetTouchable();
        const G4VPhysicalVolume *ph_vol = touchable->GetVolume(1);
        double rp_id = ph_vol->GetCopyNo();
//    edm::LogInfo("TotemRP")<<"rp id:"<<rp_id<<" name:"<<theTrack->GetDefinition()->GetParticleName()<<std::endl;

        G4cout << "Hello from G4Step - FillIfParticleLeavesFrontWallOfRP - filling\n";


        histos->set_EVT(event_no);
        histos->set_UID(particle_leaving_front_wall_of_RP_id_code);   //for particles leaving the front wall of RP
        histos->set_Ptype(theTrack->GetDefinition()->GetPDGEncoding());
        if(IsPrimary(theTrack))
            histos->set_PID(0);
        else histos->set_PID(theTrack->GetParentID());
        histos->set_ELoss(rp_id);  //rp_id in eloss
        histos->set_PABS(thePostStepPoint->GetMomentum().mag()/GeV);
        histos->set_p_x(thePostStepPoint->GetMomentum().x()/GeV);
        histos->set_p_y(thePostStepPoint->GetMomentum().y()/GeV);
        histos->set_p_z(thePostStepPoint->GetMomentum().z()/GeV);
        histos->set_VX(theTrack->GetVertexPosition().x()/mm);
        histos->set_VY(theTrack->GetVertexPosition().y()/mm);
        histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

        histos->set_X(x);
        histos->set_Y(y);
        histos->set_Z(z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        histos->fillNtuple();
//    edm::LogInfo("TotemRP")<<"leaving front wall of RP!!!! "<<theTrack->GetDefinition()->GetParticleName()<<" "<<rp_id<<std::endl;
    }
}


void PrintGeomInfoAction::FillIfParticleLeavesRP(const G4Step * aStep)
{
    const G4StepPoint* thePreStepPoint = aStep->GetPreStepPoint();
    const G4StepPoint* thePostStepPoint = aStep->GetPostStepPoint();

    if( thePostStepPoint && thePostStepPoint->GetPhysicalVolume()
        && thePreStepPoint->GetPhysicalVolume()->GetName()== "RP_box_primary_vacuum"
        && (thePostStepPoint->GetPhysicalVolume()->GetName()== "RP_220_Right_Station_Vacuum_5"
            || thePostStepPoint->GetPhysicalVolume()->GetName()== "RP_220_Left_Station_Vacuum_5"
            || thePostStepPoint->GetPhysicalVolume()->GetName()== "RP_147_Right_Station_Vacuum_5"
            || thePostStepPoint->GetPhysicalVolume()->GetName()== "RP_147_Left_Station_Vacuum_5"))
    {
        const G4Track *theTrack = aStep->GetTrack();
        const G4ThreeVector & hitPoint = thePreStepPoint->GetPosition();
        double x = hitPoint.x();
        double y = hitPoint.y();
        double z = hitPoint.z();
        double rp_id = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();

        G4cout << "Hello from G4Step - FillIfParticleLeavesRP - filling\n";

        histos->set_EVT(event_no);
        histos->set_UID(particle_leaving_RP_id_code);   //for particles going from RP station 220m
        histos->set_Ptype( theTrack->GetDefinition()->GetPDGEncoding() );
        if(IsPrimary(theTrack))
            histos->set_PID(0);
        else histos->set_PID(theTrack->GetParentID());
        histos->set_ELoss(rp_id);  //rp_id in eloss
        histos->set_PABS(thePreStepPoint->GetMomentum().mag()/GeV);
        histos->set_p_x(thePreStepPoint->GetMomentum().x()/GeV);
        histos->set_p_y(thePreStepPoint->GetMomentum().y()/GeV);
        histos->set_p_z(thePreStepPoint->GetMomentum().z()/GeV);
        histos->set_VX(theTrack->GetVertexPosition().x()/mm);
        histos->set_VY(theTrack->GetVertexPosition().y()/mm);
        histos->set_VZ(theTrack->GetVertexPosition().z()/mm);

        histos->set_X(x);
        histos->set_Y(y);
        histos->set_Z(z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        histos->fillNtuple();
//    edm::LogInfo("TotemRP")<<"leaving RP!!!! "<<theTrack->GetDefinition()->GetParticleName()<<" "<<rp_id<<std::endl;
    }
}


void PrintGeomInfoAction::update(const G4Step * aStep){
    G4cout << "Hello from: PrintGeomInfoAction.cc G4Step\n";

    const G4ThreeVector& pre_step_pos = aStep->GetPreStepPoint()->GetPosition();
    const G4ThreeVector& post_step_pos = aStep->GetPostStepPoint()->GetPosition();

    if(verbosity_ && aStep->GetTrack()->GetDefinition()->GetPDGEncoding()==2212)
    {
        edm::LogInfo("TotemRP") << aStep->GetTrack()->GetDefinition()->GetParticleName() << " "
                                << "trackID:" << aStep->GetTrack()->GetTrackID()
                                << " parentID:" << aStep->GetTrack()->GetParentID() << " "
                                << aStep->GetPreStepPoint()->GetPosition()<<" "
                                << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()<< " "
                                << aStep->GetPostStepPoint()->GetPosition()<<" "
                                <<"primary_vert="<<aStep->GetTrack()->GetVertexPosition()<<std::endl;
    }

    FillIfLeavesRP220Station(aStep);

//do not uncomment ever in cmssw with parameterization
    //KillParticlesBetweenStations(aStep);

    FillIfParticleEntersRP(aStep);
    FillIfParticleLeavesRP(aStep);
    FillIfParticleLeavesFrontWallOfRP(aStep);

    if(verbosity_)
    {
        edm::LogInfo("TotemRP")<<"step info, "<<aStep->GetTrack()->GetTrackID()<<", pre_step:";
        if(aStep && aStep->GetPreStepPoint() && aStep->GetPreStepPoint()->GetPhysicalVolume())
            edm::LogInfo("TotemRP")<<aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()<<" "
                                   <<pre_step_pos<<" "<<aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
        edm::LogInfo("TotemRP")<<", post_step:";
        if(aStep && aStep->GetPostStepPoint() && aStep->GetPostStepPoint()->GetPhysicalVolume())
            edm::LogInfo("TotemRP")<<aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()<<" "
                                   <<post_step_pos<<" "<<aStep->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo();
        edm::LogInfo("TotemRP")<<std::endl;
    }
}

//=================================================================== track

void PrintGeomInfoAction::update(const BeginOfTrack * beg_of_track){
    G4cout << "Hello from: PrintGeomInfoAction.cc BeginOfTrack\n";

    const G4Track *track = (*beg_of_track)();

    TrackInformation* info
            = dynamic_cast<TrackInformation*>( (*beg_of_track)()->GetUserInformation() );


    if(track->GetDynamicParticle()->GetPDGcode()==2212
       && track->GetCreatorProcess()
       && track->GetCreatorProcess()->GetProcessType()==fParameterisation
       && track->GetCreatorProcess()->GetProcessName()=="TotemRPParameterisationProcess"
       && info && info->isPrimary())
    {
        G4ThreeVector vert = track->GetVertexPosition();

        G4ThreeVector position = track->GetPosition();
        double pos_x = position.x()/mm;
        double pos_y = position.y()/mm;
        double pos_z = position.z()/mm;

        G4cout << "Hello from if in BeginOfTrack - filling histos\n";


        histos->set_EVT(event_no);
        histos->set_UID(particle_entering_station_id_code);   //for particles going from RP station 220m
        histos->set_Ptype(2212);

        histos->set_ELoss(0);
        histos->set_PABS(track->GetMomentum().mag()/GeV);
        histos->set_p_x(track->GetMomentum().x()/GeV);
        histos->set_p_y(track->GetMomentum().y()/GeV);
        histos->set_p_z(track->GetMomentum().z()/GeV);
        histos->set_VX(vert.x()/mm);
        histos->set_VY(vert.y()/mm);
        histos->set_VZ(vert.z()/mm);

        histos->set_X(pos_x);
        histos->set_Y(pos_y);
        histos->set_Z(pos_z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        histos->fillNtuple();
    }
}

void PrintGeomInfoAction::update(const EndOfTrack * end_of_track){
    const G4Track *track = (*end_of_track)();

    G4cout << "Hello from: PrintGeomInfoAction.cc EndOfTrack\n"
           << "IsPrimary(track): " << IsPrimary(track) << "\n"
           << "track->GetDefinition()->GetPDGEncoding(): " << track->GetDefinition()->GetPDGEncoding() <<"\n"
           << "track->GetVolume()->GetName(): " << track->GetVolume()->GetName() << "\n";

    if(IsPrimary(track) && track->GetDefinition()->GetPDGEncoding()==2212 &&
       track->GetVolume()->GetName().substr(0, 3)=="RP_")
    {
        const G4VTouchable* touchable = track->GetTouchable();
        int hist_depth = touchable->GetHistoryDepth();
        double rp_id=-1;
        for(int i=0; i<hist_depth; ++i)
        {
            const G4VPhysicalVolume *ph_vol = touchable->GetVolume(i);
            if(ph_vol->GetName()=="RP_box_primary_vacuum")
            {
                rp_id = ph_vol->GetCopyNo();
                break;
            }
        }

        const G4ThreeVector & position = track->GetPosition();
        double x = position.x()/mm;
        double y = position.y()/mm;
        double z = position.z()/mm;

        G4cout << "Hello from if in EndOfTrack - filling histos\n";

        histos->set_EVT(event_no);
        histos->set_UID(primary_proton_inelastic_event_in_RP_station);   //for particles going from RP station 220m
        histos->set_Ptype(2212);

        histos->set_ELoss(rp_id);  //rp_id in eloss
        histos->set_PABS(track->GetMomentum().mag()/GeV);
        histos->set_p_x(track->GetMomentum().x()/GeV);
        histos->set_p_y(track->GetMomentum().y()/GeV);
        histos->set_p_z(track->GetMomentum().z()/GeV);
        histos->set_VX(track->GetVertexPosition().x()/mm);
        histos->set_VY(track->GetVertexPosition().y()/mm);
        histos->set_VZ(track->GetVertexPosition().z()/mm);

        histos->set_X(x);
        histos->set_Y(y);
        histos->set_Z(z);
        histos->set_Loc_X(0);
        histos->set_Loc_Y(0);
        histos->set_Loc_Z(0);
        histos->set_X_Exit(0);
        histos->set_Y_Exit(0);
        histos->set_Z_Exit(0);
        histos->set_Loc_X_Exit(0);
        histos->set_Loc_Y_Exit(0);
        histos->set_Loc_Z_Exit(0);

        G4VPhysicalVolume* myVolume = track->GetVolume();
        int vol_copy_no = myVolume->GetCopyNo();
        int prim_vert_id = PhysicalDetMap[myVolume->GetName()];
        histos->set_prim_ver_id(prim_vert_id);
        histos->set_PID(vol_copy_no);

        histos->fillNtuple();
//    edm::LogInfo("TotemRP")<<"prim proton died "<<prim_vert_id<<std::endl;
    } else {
        G4cout << "Hello from if in EndOfTrack - NO!!! filling!!!! histos\n";
    }


    // SHOULD BE REMOVED START
    histos->rt_hf->cd();
    histos->ntuple->Write();

    histos->rt_hf->Close();
    //delete rt_hf;

    edm::LogInfo("TotemRP") << std::endl << "TotemRPHisto - from my simwatcher:P: End writing user histograms " << std::endl;
    // SHOULD BE REMOBED END
}

// TotemRP specific END


PrintGeomInfoAction::~PrintGeomInfoAction() {
   // G4cout << "Hello from: Destructive destructor\n";
    delete histos;
}

// duplicated in TotemRPs
//void PrintGeomInfoAction::update(const BeginOfJob * job) {
//
////  if (_dumpSense) {
////    edm::ESTransientHandle<DDCompactView> pDD;
////    (*job)()->get<IdealGeometryRecord>().get(pDD);
////
////    G4cout << "PrintGeomInfoAction::Get Printout of Sensitive Volumes "
////	   << "for " << names.size() << " Readout Units" << G4endl;
////    for (unsigned int i=0; i<names.size(); i++) {
////      std::string attribute = "ReadOutName";
////      std::string sd        = names[i];
////      DDSpecificsMatchesValueFilter filter{DDValue(attribute,sd,0)};
////      DDFilteredView fv(*pDD,filter);
////      G4cout << "PrintGeomInfoAction:: Get Filtered view for "
////	     << attribute << " = " << sd << G4endl;
////      bool dodet = fv.firstChild();
////
////      std::string spaces = spacesFromLeafDepth(1);
////
////      while (dodet) {
////	const DDLogicalPart & log = fv.logicalPart();
////	std::string lvname = log.name().name();
////	DDTranslation tran = fv.translation();
////	std::vector<int> copy = fv.copyNumbers();
////
////	unsigned int leafDepth = copy.size();
////	G4cout << leafDepth << spaces << "### VOLUME = " << lvname
////	       << " Copy No";
////	for (int k=leafDepth-1; k>=0; k--) G4cout << " " << copy[k];
////	G4cout << " Centre at " << tran << " (r = " << tran.Rho()
////	       << ", phi = " << tran.phi()/deg << ")" << G4endl;
////	dodet = fv.next();
////      }
////    }
////  }
//}
  
void PrintGeomInfoAction::update(const BeginOfRun * run) {

//  theTopPV = getTopPV();
//
//  if (_dumpSummary)  dumpSummary(G4cout);
//  if (_dumpLVTree)   dumpG4LVTree(G4cout);
//
//  //---------- Dump list of objects of each class with detail of parameters
//  if (_dumpMaterial) dumpMaterialList(G4cout);
//  if (_dumpLVList)   dumpG4LVList(G4cout);
//
//  //---------- Dump LV and PV information
//  if (_dumpLV || _dumpPV || _dumpTouch) dumpHierarchyTreePVLV(G4cout);
}

void PrintGeomInfoAction::dumpSummary(std::ostream & out) {

//  //---------- Dump number of objects of each class
//  out << " @@@@@@@@@@@@@@@@@@ Dumping G4 geometry objects Summary " << G4endl;
//  if (theTopPV == 0) {
//    out << " No volume created " << G4endl;
//    return;
//  }
//  out << " @@@ Geometry built inside world volume: " << theTopPV->GetName() << G4endl;
//  // Get number of solids (< # LV if several LV share a solid)
//  const G4LogicalVolumeStore * lvs = G4LogicalVolumeStore::GetInstance();
//  std::vector<G4LogicalVolume *>::const_iterator lvcite;
//  std::set<G4VSolid *> theSolids;
//  for (lvcite = lvs->begin(); lvcite != lvs->end(); lvcite++)
//    theSolids.insert((*lvcite)->GetSolid());
//  out << " Number of G4VSolid's: " << theSolids.size() << G4endl;
//  out << " Number of G4LogicalVolume's: " << lvs->size() << G4endl;
//  const G4PhysicalVolumeStore * pvs = G4PhysicalVolumeStore::GetInstance();
//  out << " Number of G4VPhysicalVolume's: " << pvs->size() << G4endl;
//  out << " Number of Touchable's: " << countNoTouchables() << G4endl;
//  const G4MaterialTable * matTab = G4Material::GetMaterialTable();
//  out << " Number of G4Material's: " << matTab->size() << G4endl;
}

void PrintGeomInfoAction::dumpG4LVList(std::ostream & out) {

//  out << " @@@@@@@@@@@@@@@@ DUMPING G4LogicalVolume's List  " << G4endl;
//  const G4LogicalVolumeStore * lvs = G4LogicalVolumeStore::GetInstance();
//  std::vector<G4LogicalVolume*>::const_iterator lvcite;
//  for (lvcite = lvs->begin(); lvcite != lvs->end(); lvcite++)
//    out << "LV:" << (*lvcite)->GetName() << "\tMaterial: " << (*lvcite)->GetMaterial()->GetName() << G4endl;
}

void PrintGeomInfoAction::dumpG4LVTree(std::ostream & out) {

//  out << " @@@@@@@@@@@@@@@@ DUMPING G4LogicalVolume's Tree  " << G4endl;
//  G4LogicalVolume * lv = getTopLV();
//  dumpG4LVLeaf(lv,0,1,out);
}

void PrintGeomInfoAction::dumpMaterialList(std::ostream & out) {
//
//  out << " @@@@@@@@@@@@@@@@ DUMPING G4Material List ";
//  const G4MaterialTable * matTab = G4Material::GetMaterialTable();
//  out << " with " << matTab->size() << " materials " << G4endl;
//  std::vector<G4Material*>::const_iterator matite;
//  for (matite = matTab->begin(); matite != matTab->end(); matite++)
//    out << "Material: " << (*matite) << G4endl;
}

void PrintGeomInfoAction::dumpG4LVLeaf(G4LogicalVolume * lv, unsigned int leafDepth, unsigned int count, std::ostream & out) {

//  for (unsigned int ii=0; ii < leafDepth; ii++) out << "  ";
//  out << " LV:(" << leafDepth << ") " << lv->GetName() << " (" << count
//      << ")" << G4endl;
//  //--- If a volume is placed n types as daughter of this LV, it should only be counted once
//  std::map<G4LogicalVolume*, unsigned int> lvCount;
//  std::map<G4LogicalVolume*, unsigned int>::const_iterator cite;
//  for (int ii = 0; ii < lv->GetNoDaughters(); ii++) {
//    cite = lvCount.find(lv->GetDaughter(ii)->GetLogicalVolume());
//    if (cite != lvCount.end()) lvCount[cite->first] = (cite->second) + 1;
//    else lvCount.insert(std::pair< G4LogicalVolume*,unsigned int>(lv->GetDaughter(ii)->GetLogicalVolume(),1));
//  }
//  for (cite = lvCount.begin(); cite != lvCount.end(); cite++)
//    dumpG4LVLeaf((cite->first), leafDepth+1, (cite->second), out);
}

int PrintGeomInfoAction::countNoTouchables() {

  int nTouch = 0;
  G4LogicalVolume * lv = getTopLV();
  add1touchable(lv, nTouch);
  return nTouch;
}

void PrintGeomInfoAction::add1touchable(G4LogicalVolume * lv, int & nTouch) {
//
//  int siz = lv->GetNoDaughters();
//  for(int ii = 0; ii < siz; ii++)
//    add1touchable(lv->GetDaughter(ii)->GetLogicalVolume(), ++nTouch);
}
 
void PrintGeomInfoAction::dumpHierarchyTreePVLV(std::ostream & out) {

//  //dumps in the following order:
//  //    1) a LV with details
//  //    2) list of PVs daughters of this LV with details
//  //    3) list of LVs daughters of this LV and for each go to 1)
//
//  //----- Get top PV
//  G4LogicalVolume*  topLV = getTopLV();
//
//  //----- Dump this leaf (it will recursively dump all the tree)
//  dumpHierarchyLeafPVLV(topLV, 0, out);
//  dumpPV(theTopPV, 0, out);
//
//  //----- Dump the touchables (it will recursively dump all the tree)
//  if (_dumpTouch) dumpTouch(theTopPV, 0, out);
}

void PrintGeomInfoAction::dumpHierarchyLeafPVLV(G4LogicalVolume * lv, unsigned int leafDepth, std::ostream & out) {
//
//  //----- Dump this LV
//  dumpLV(lv, leafDepth, out);
//
//  //----- Get LV daughters from list of PV daughters
//  mmlvpv lvpvDaughters;
//  std::set< G4LogicalVolume * > lvDaughters;
//  int NoDaughters = lv->GetNoDaughters();
//  while ((NoDaughters--)>0) {
//    G4VPhysicalVolume * pvD = lv->GetDaughter(NoDaughters);
//    lvpvDaughters.insert(mmlvpv::value_type(pvD->GetLogicalVolume(), pvD));
//    lvDaughters.insert(pvD->GetLogicalVolume());
//  }
//
//  std::set< G4LogicalVolume * >::const_iterator scite;
//  mmlvpv::const_iterator mmcite;
//
//  //----- Dump daughters PV and LV
//  for (scite = lvDaughters.begin(); scite != lvDaughters.end(); scite++) {
//    std::pair< mmlvpv::iterator, mmlvpv::iterator > mmER = lvpvDaughters.equal_range(*scite);
//    //----- Dump daughters PV of this LV
//    for (mmcite = mmER.first ; mmcite != mmER.second; mmcite++)
//      dumpPV((*mmcite).second, leafDepth+1, out);
//    //----- Dump daughters LV
//    dumpHierarchyLeafPVLV(*scite, leafDepth+1, out );
//  }
}
 
void PrintGeomInfoAction::dumpLV(G4LogicalVolume * lv, unsigned int leafDepth, std::ostream & out) {
//
//  std::string spaces = spacesFromLeafDepth(leafDepth);
//
//  //----- dump name
//  if (_dumpLV) {
//    out << leafDepth << spaces << "$$$ VOLUME = " << lv->GetName()
//	<< "  Solid: " << lv->GetSolid()->GetName() << "  MATERIAL: "
//	<< lv->GetMaterial()->GetName() << G4endl;
//    if (_dumpSolid)
//      dumpSolid(lv->GetSolid(), leafDepth, out); //----- dump solid
//
//    //----- dump LV info
//    //--- material
//    if (_dumpAtts) {
//      //--- Visualisation attributes
//      const G4VisAttributes * fVA = lv->GetVisAttributes();
//      if (fVA!=0) {
//	out <<  spaces << "  VISUALISATION ATTRIBUTES: " << G4endl;
//	out <<  spaces << "    IsVisible " << fVA->IsVisible() << G4endl;
//	out <<  spaces << "    IsDaughtersInvisible " << fVA->IsDaughtersInvisible() << G4endl;
//	out <<  spaces << "    Colour " << fVA->GetColour() << G4endl;
//	out <<  spaces << "    LineStyle " << fVA->GetLineStyle() << G4endl;
//	out <<  spaces << "    LineWidth " << fVA->GetLineWidth() << G4endl;
//	out <<  spaces << "    IsForceDrawingStyle " << fVA->IsForceDrawingStyle() << G4endl;
//	out <<  spaces << "    ForcedDrawingStyle " << fVA->GetForcedDrawingStyle() << G4endl;
//      }
//
//      //--- User Limits
//      G4UserLimits * fUL = lv->GetUserLimits();
//      G4Track dummy;
//      if (fUL!=0) {
//	out <<  spaces << "    MaxAllowedStep " << fUL->GetMaxAllowedStep(dummy) << G4endl;
//	out <<  spaces << "    UserMaxTrackLength " << fUL->GetUserMaxTrackLength(dummy) << G4endl;
//	out <<  spaces << "    UserMaxTime " << fUL->GetUserMaxTime(dummy) << G4endl;
//	out <<  spaces << "    UserMinEkine " << fUL->GetUserMinEkine(dummy) << G4endl;
//	out <<  spaces << "    UserMinRange " << fUL->GetUserMinRange(dummy) << G4endl;
//      }
//
//      //--- other LV info
//      if (lv->GetSensitiveDetector())
//	out << spaces << "  IS SENSITIVE DETECTOR " << G4endl;
//      if (lv->GetFieldManager())
//	out << spaces << "  FIELD ON " << G4endl;
//
//      // Pointer (possibly NULL) to optimisation info objects.
//      out <<  spaces
//	  << "        Quality for optimisation, average number of voxels to be spent per content "
//	  << lv->GetSmartless() << G4endl;
//
//      // Pointer (possibly NULL) to G4FastSimulationManager object.
//      if (lv->GetFastSimulationManager())
//	out << spaces << "     Logical Volume is an envelope for a FastSimulationManager "
//	    << G4endl;
//      out << spaces << "     Weight used in the event biasing technique = "
//	  << lv->GetBiasWeight() << G4endl;
//    }
//  }
}	

void PrintGeomInfoAction::dumpPV(G4VPhysicalVolume * pv, unsigned int leafDepth, std::ostream & out) {
//
//  std::string spaces = spacesFromLeafDepth(leafDepth);
//
//  //----- PV info
//  if (_dumpPV) {
//    std::string mother = "World";
//    if (pv->GetMotherLogical()) mother = pv->GetMotherLogical()->GetName();
//    out << leafDepth << spaces << "### VOLUME = " << pv->GetName()
//	<< " Copy No " << pv->GetCopyNo() << " in " << mother
//	<< " at " << pv->GetTranslation();
//  }
//  if (!pv->IsReplicated()) {
//    if (_dumpPV) {
//      if(pv->GetRotation() == 0) out << " with no rotation" << G4endl;
//      else  if(!_dumpRotation)   out << " with rotation" << G4endl; //just rotation name
//      else                       out << " with rotation " << *(pv->GetRotation()) << G4endl;
//    }
//  } else {
//    if (_dumpReplica ) {
//      out << spaces << "    It is replica: " << G4endl;
//      EAxis axis;
//      int nReplicas;
//      double width;
//      double offset;
//      bool consuming;
//      pv->GetReplicationData(axis, nReplicas, width, offset, consuming);
//      out << spaces << "     axis " << axis << G4endl
//	  << spaces << "     nReplicas " << nReplicas << G4endl;
//      if (pv->GetParameterisation() != 0)
//	out << spaces << "    It is parameterisation " << G4endl;
//      else
//	out << spaces << "     width " << width << G4endl
//	    << spaces << "     offset " << offset << G4endl
//	    << spaces << "     consuming" <<  consuming << G4endl;
//      if (pv->GetParameterisation() != 0)
//	out << spaces << "    It is parameterisation " << G4endl;
//    }
//  }
}

void PrintGeomInfoAction::dumpTouch(G4VPhysicalVolume * pv, unsigned int leafDepth, std::ostream & out) {
//
//  std::string spaces = spacesFromLeafDepth(leafDepth);
//  if (leafDepth == 0) fHistory.SetFirstEntry(pv);
//  else fHistory.NewLevel(pv, kNormal, pv->GetCopyNo());
//
//  G4ThreeVector globalpoint = fHistory.GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,0));
//  G4LogicalVolume * lv = pv->GetLogicalVolume();
//
//  std::string mother = "World";
//  if (pv->GetMotherLogical()) mother = pv->GetMotherLogical()->GetName();
//  std::string lvname = lv->GetName();
//  lvname.assign(lvname,0,nchar);
//  if (lvname == name)
//    out << leafDepth << spaces << "### VOLUME = " << lv->GetName()
//	<< " Copy No " << pv->GetCopyNo() << " in " << mother
//	<< " global position of centre " << globalpoint << " (r = "
//	<<  globalpoint.perp() << ", phi = " <<  globalpoint.phi()/deg
//	<< ")" << G4endl;
//
//  int NoDaughters = lv->GetNoDaughters();
//  while ((NoDaughters--)>0) {
//    G4VPhysicalVolume * pvD = lv->GetDaughter(NoDaughters);
//    if (!pvD->IsReplicated()) dumpTouch(pvD, leafDepth+1, out);
//  }
//
//  if (leafDepth > 0) fHistory.BackLevel();
}

std::string PrintGeomInfoAction::spacesFromLeafDepth(unsigned int leafDepth) {

  std::string spaces;
  unsigned int ii;
  for(ii = 0; ii < leafDepth; ii++) { spaces += "  "; }

  return spaces;
}

void PrintGeomInfoAction::dumpSolid(G4VSolid * sol, unsigned int leafDepth, std::ostream & out) {
//
//  std::string spaces = spacesFromLeafDepth(leafDepth);
//  out << spaces << *(sol) << G4endl;
}

G4VPhysicalVolume * PrintGeomInfoAction::getTopPV() {

  return G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
}

G4LogicalVolume * PrintGeomInfoAction::getTopLV() { 
  return theTopPV->GetLogicalVolume();
}


