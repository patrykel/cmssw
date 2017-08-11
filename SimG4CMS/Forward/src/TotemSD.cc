// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemSD
//
// Implementation:
//     <Notes on implementation>
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
//

// system include files

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Notification/interface/G4TrackToParticleID.h"
#include "SimG4Core/Geometry/interface/SensitiveDetectorCatalog.h"
#include "SimG4Core/Physics/interface/G4ProcessTypeEnumerator.h"

#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/SimHitMaker/interface/TrackingSlaveSD.h"

#include "SimG4CMS/Forward/interface/TotemSD.h"
#include "SimG4CMS/Forward/interface/TotemNumberMerger.h"
#include "SimG4CMS/Forward/interface/TotemT1NumberingScheme.h"
#include "SimG4CMS/Forward/interface/TotemT2NumberingSchemeGem.h"
#include "SimG4CMS/Forward/interface/TotemRPNumberingScheme.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <iostream>
#include <vector>
#include <string>
#include "FWCore/MessageLogger/interface/MessageLogger.h"


//
// constructors and destructor
//
TotemSD::TotemSD(std::string name, const DDCompactView & cpv,
		 const SensitiveDetectorCatalog & clg,
		 edm::ParameterSet const & p, const SimTrackManager* manager) :
  SensitiveTkDetector(name, cpv, clg, p), numberingScheme(0), name(name),
  hcID(-1), theHC(0), theManager(manager), currentHit(0), theTrack(0), 
  currentPV(0), unitID(0),  previousUnitID(0), preStepPoint(0), 
  postStepPoint(0), eventno(0){

  //Add Totem Sentitive Detector Names
  collectionName.insert(name);

  //Parameters
  edm::ParameterSet m_p = p.getParameter<edm::ParameterSet>("TotemSD");
  int verbn = m_p.getUntrackedParameter<int>("Verbosity");
 
  SetVerboseLevel(verbn);
  LogDebug("ForwardSim") 
    << "*******************************************************\n"
    << "*                                                     *\n"
    << "* Constructing a TotemSD  with name " << name << "\n"
    << "*                                                     *\n"
    << "*******************************************************";

  slave  = new TrackingSlaveSD(name);

  //
  // Now attach the right detectors (LogicalVolumes) to me
  //
  const std::vector<std::string>& lvNames = clg.logicalNames(name);
  this->Register();
  for (std::vector<std::string>::const_iterator it=lvNames.begin();
       it !=lvNames.end(); it++) {
    this->AssignSD(*it);
    edm::LogInfo("ForwardSim") << "TotemSD : Assigns SD to LV " << (*it);
  }

  if      (name == "TotemHitsT1") {
    numberingScheme = dynamic_cast<TotemVDetectorOrganization*>(new TotemT1NumberingScheme(1));
  } else if (name == "TotemHitsT2Si") {
    numberingScheme = dynamic_cast<TotemVDetectorOrganization*>(new TotemT2NumberingSchemeGem(3));
  } else if (name == "TotemHitsT2Gem") {
    numberingScheme = dynamic_cast<TotemVDetectorOrganization*>(new TotemT2NumberingSchemeGem(4));
  } else if (name == "TotemHitsRP") {
    numberingScheme = dynamic_cast<TotemVDetectorOrganization*>(new TotemRPNumberingScheme(3));
  } else {
    edm::LogWarning("ForwardSim") << "TotemSD: ReadoutName not supported\n";
  }
  
  edm::LogInfo("ForwardSim") << "TotemSD: Instantiation completed";
} 

TotemSD::~TotemSD() { 
  if (slave)           delete slave; 
  if (numberingScheme) delete numberingScheme;
}

G4bool TotemSD::ProcessHits(G4Step * aStep, G4TouchableHistory * )
{
  if (aStep == NULL)
  {
  	return true;
  }
  else
  {
    GetStepInfo(aStep);
      CreateNewHit();
      //LogDebug("TotemRP")<<"New hit created"<<std::endl;
	    return true;
	}
}

void TotemSD::Print_Hit_Info()
{
  LogDebug("TotemRP") << theTrack->GetDefinition()->GetParticleName()
       << " TotemSD CreateNewHit for"
       << " PV "     << currentPV->GetName()
       << " PVid = " << currentPV->GetCopyNo()
       << " Unit "   << unitID;
  LogDebug("TotemRP") << " primary "    << primaryID
       << " time slice " << tSliceID
       << " of energy " << theTrack->GetTotalEnergy()
       << " Eloss " << Eloss
       << " positions ";
       printf("(%10f,%10f,%10f)",preStepPoint->GetPosition().x(),preStepPoint->GetPosition().y(),preStepPoint->GetPosition().z());
       printf("(%10f,%10f,%10f)",postStepPoint->GetPosition().x(),postStepPoint->GetPosition().y(),postStepPoint->GetPosition().z());
  LogDebug("TotemRP") << " positions " << "(" <<postStepPoint->GetPosition().x()<<","<<postStepPoint->GetPosition().y()<<","<<postStepPoint->GetPosition().z()<<")"
       << " For Track  " << theTrack->GetTrackID()
       << " which is a " << theTrack->GetDefinition()->GetParticleName();

  if(theTrack->GetTrackID()==1)
  {
    LogDebug("TotemRP") << " primary particle ";
  }
  else
  {
    LogDebug("TotemRP") << " daughter of part. " << theTrack->GetParentID();
  }

  LogDebug("TotemRP")  << " and created by " ;

  if(theTrack->GetCreatorProcess()!=NULL)
    LogDebug("TotemRP") << theTrack->GetCreatorProcess()->GetProcessName() ;
  else
    LogDebug("TotemRP") << "NO process";

  LogDebug("TotemRP") << std::endl;
}

uint32_t TotemSD::setDetUnitId(G4Step * aStep) { 

  return (numberingScheme == 0 ? 0 : numberingScheme->GetUnitID(aStep));
}

void TotemSD::Initialize(G4HCofThisEvent * HCE) {
  LogDebug("TotemRP") << "TotemSD : Initialize called for " << name;

  theHC = new TotemG4HitCollection(name, collectionName[0]);
  G4SDManager::GetSDMpointer()->AddNewCollection(name, collectionName[0]);

  if (hcID<0)
    hcID = G4SDManager::GetSDMpointer()->GetCollectionID(theHC);
  HCE->AddHitsCollection(hcID, theHC);
}


void TotemSD::EndOfEvent(G4HCofThisEvent* )
{
  // here we loop over transient hits and make them persistent
  for (int j=0; j<theHC->entries() && j<15000; j++) {
    TotemG4Hit* aHit = (*theHC)[j];

    Local3DPoint Entrata(aHit->getLocalEntry().x(),
       aHit->getLocalEntry().y(),
       aHit->getLocalEntry().z());
    Local3DPoint Uscita(aHit->getLocalExit().x(),
       aHit->getLocalExit().y(),
       aHit->getLocalExit().z());
    slave->processHits(PSimHit(Entrata,Uscita,
             aHit->getPabs(), aHit->getTof(),
             aHit->getEnergyLoss(), aHit->getParticleType(),
             aHit->getUnitID(), aHit->getTrackID(),
             aHit->getThetaAtEntry(),aHit->getPhiAtEntry()));
  }
  Summarize();
}

void TotemSD::clear() {
} 

void TotemSD::DrawAll() {
} 

void TotemSD::PrintAll() {
  LogDebug("ForwardSim") << "TotemSD: Collection " << theHC->GetName();
  theHC->PrintAllHits();
} 

void TotemSD::fillHits(edm::PSimHitContainer& c, std::string n) {
  if (slave->name() == n) c=slave->hits();
}

void TotemSD::update (const BeginOfEvent * i) {
  LogDebug("ForwardSim") << " Dispatched BeginOfEvent for " << GetName()
                       << " !" ;
   clearHits();
   eventno = (*i)()->GetEventID();
}

void TotemSD::update (const ::EndOfEvent*) {
}

void TotemSD::clearHits(){
  slave->Initialize();
}

void TotemSD::SetNumberingScheme(TotemVDetectorOrganization* scheme)
{
  if (numberingScheme)
    delete numberingScheme;
  numberingScheme = scheme;
}


G4ThreeVector TotemSD::SetToLocal(G4ThreeVector global)
{
  G4ThreeVector localPoint;
  const G4VTouchable* touch= preStepPoint->GetTouchable();
  localPoint = touch->GetHistory()->GetTopTransform().TransformPoint(global);

  return localPoint;
}


void TotemSD::GetStepInfo(G4Step* aStep)
{
  preStepPoint = aStep->GetPreStepPoint();
  postStepPoint = aStep->GetPostStepPoint();
  theTrack = aStep->GetTrack();
  hitPoint = preStepPoint->GetPosition();
  exitPoint = postStepPoint->GetPosition();
  currentPV = preStepPoint->GetPhysicalVolume();
  theLocalEntryPoint = SetToLocal(hitPoint);
  theLocalExitPoint = SetToLocal(exitPoint);

  G4String name = currentPV->GetName();
  name.assign(name,0,4);

  G4String particleType = theTrack->GetDefinition()->GetParticleName();

  tSlice = (postStepPoint->GetGlobalTime() )/nanosecond;
  tSliceID = (int) tSlice;
  unitID = setDetUnitId(aStep);

  if(verbosity_)
    LogDebug("TotemRP") << "UNITa " << unitID <<std::endl;

  primaryID = theTrack->GetTrackID();

  Pabs = (aStep->GetPreStepPoint()->GetMomentum().mag())/GeV;
  p_x = (aStep->GetPreStepPoint()->GetMomentum().x())/GeV;
  p_y = (aStep->GetPreStepPoint()->GetMomentum().y())/GeV;
  p_z = (aStep->GetPreStepPoint()->GetMomentum().z())/GeV;

  Tof = aStep->GetPostStepPoint()->GetGlobalTime()/nanosecond;
  Eloss = aStep->GetTotalEnergyDeposit()/GeV;
  ParticleType = theTrack->GetDefinition()->GetPDGEncoding();

  //corrected phi and theta treatment
  G4ThreeVector gmd  = aStep->GetPreStepPoint()->GetMomentumDirection();
  // convert it to local frame
  G4ThreeVector lmd = ((G4TouchableHistory *)(aStep->GetPreStepPoint()->GetTouchable()))->GetHistory()
    ->GetTopTransform().TransformAxis(gmd);
  Local3DPoint lnmd = ConvertToLocal3DPoint(lmd);
  ThetaAtEntry = lnmd.theta();
  PhiAtEntry = lnmd.phi();


//  LogDebug("TotemRP") << "UUUUUUUNNNNNNNNNNIIIIIIIIIITTTTTTTTTTTTTIIIIDDDD " <<
//    numberingScheme->GetUnitID(aStep) << std::endl ;
  if(IsPrimary(theTrack))
    ParentId = 0;
  else ParentId = theTrack->GetParentID();

  Vx = theTrack->GetVertexPosition().x()/mm;
  Vy = theTrack->GetVertexPosition().y()/mm;
  Vz = theTrack->GetVertexPosition().z()/mm;
}

void TotemSD::CreateNewHit()
{
  currentHit = new TotemG4Hit;
  currentHit->setTrackID(primaryID);
  currentHit->setTimeSlice(tSlice);
  currentHit->setUnitID(unitID);
  currentHit->setIncidentEnergy(incidentEnergy);

  currentHit->setPabs(Pabs);
  currentHit->setTof(Tof);
  currentHit->setEnergyLoss(Eloss);
  currentHit->setParticleType(ParticleType);
  currentHit->setThetaAtEntry(ThetaAtEntry);
  currentHit->setPhiAtEntry(PhiAtEntry);

  currentHit->setEntry(hitPoint);
  currentHit->setExit(exitPoint);
  currentHit->setLocalEntry(theLocalEntryPoint);
  currentHit->setLocalExit(theLocalExitPoint);

  currentHit->setParentId(ParentId);
  currentHit->setVx(Vx);
  currentHit->setVy(Vy);
  currentHit->setVz(Vz);

  currentHit->set_p_x(p_x);
  currentHit->set_p_y(p_y);
  currentHit->set_p_z(p_z);

  StoreHit(currentHit);
// LogDebug("TotemRP") << "STORED HIT IN: " << unitID << std::endl;
}
 
G4ThreeVector TotemSD::PosizioEvo(const G4ThreeVector& Pos, double vx, double vy,
				  double vz, double pabs, int& accettanza) {
  accettanza=0;
  G4ThreeVector PosEvo;
  double ThetaX=atan((Pos.x()-vx)/(Pos.z()-vz));                 
  double ThetaY=atan((Pos.y()-vy)/(Pos.z()-vz));                
  double X_at_0 =(vx-((Pos.x()-vx)/(Pos.z()-vz))*vz)/1000.;   
  double Y_at_0 =(vy-((Pos.y()-vy)/(Pos.z()-vz))*vz)/1000.;
 
  //csi=-dp/d
  double csi = fabs((7000.-pabs)/7000.);

  // all in m 
  const int no_rp=4;
  double x_par[no_rp+1];
  double y_par[no_rp+1];
  //rp z position
  double rp[no_rp]={141.,149.,198.,220.};
  //{lx0,mlx} for each rp; Lx=lx0+mlx*csi
  double leffx[][2]={{122.5429,-46.9312},{125.4194,-49.1849},{152.6,-81.157},{98.8914,-131.8390}};
  //{ly0,mly} for each rp; Ly=ly0+mly*csi
  double leffy[][2]={{124.2314,-55.4852},{127.7825,-57.4503},{179.455,-76.274},{273.0931,-40.4626}};
  //{vx0,mvx0} for each rp; vx=vx0+mvx*csi
  double avx[][2]={{0.515483,-1.0123},{0.494122,-1.0534},{0.2217,-1.483},{0.004633,-1.0719}};
  //{vy0,mvy0} for each rp; vy=vy0+mvy*csi
  double avy[][2]={{0.371418,-1.6327},{0.349035,-1.6955},{0.0815,-2.59},{0.007592,-4.0841}};                
  //{D0,md,a,b} for each rp; D=D0+(md+a*thetax)*csi+b*thetax
  double ddx[][4]= {{-0.082336,-0.092513,112.3436,-82.5029},{-0.086927,-0.097670,114.9513,-82.9835},
		    {-0.092117,-0.0915,180.6236,-82.443},{-0.050470,0.058837,208.1106,20.8198}};
  // {10sigma_x+0.5mm,10sigma_y+0.5mm}
  double detlim[][2]={{0,0},{0.0028,0.0021},{0,0},{0.0008,0.0013}};   
  //{rmax,dmax}
  double pipelim[][2]={{0.026,0.026},{0.04,0.04},{0.0226,0.0177},{0.04,0.04}};
  
  
  for(int j=0; j<no_rp ; j++)  { 
    //y=Ly*thetay+vy*y0
    //x=Lx*thetax+vx*x0-csi*D   
    y_par[j]=ThetaY*(leffy[j][0]+leffy[j][1]*csi)+(avy[j][0]+avy[j][1]*csi)*Y_at_0;
    x_par[j]=ThetaX*(leffx[j][0]+leffx[j][1]*csi)+(avx[j][0]+avx[j][1]*csi)*X_at_0-
      csi*(ddx[j][0]+(ddx[j][1]+ddx[j][2]*ThetaX)*csi+ddx[j][3]*ThetaX);
  }
   
   
  //pass TAN@141
  if (fabs(y_par[0])<pipelim[0][1] && sqrt((y_par[0]*y_par[0])+(x_par[0]*x_par[0]))<pipelim[0][0])  {
    //pass 149
    if ((sqrt((y_par[1]*y_par[1])+(x_par[1]*x_par[1]))<pipelim[1][0]) &&
	(fabs(y_par[1])>detlim[1][1] || x_par[1]>detlim[1][0]))  {
      accettanza = 1;
    }
  }

      
  //pass TAN@141
  if (fabs(y_par[0])<pipelim[0][1] && sqrt((y_par[0])*(y_par[0])+(x_par[0])*(x_par[0]))<pipelim[0][0]) {
    //pass Q5@198
    if (fabs(y_par[2])<pipelim[2][1] && sqrt((y_par[2]*y_par[2])+(x_par[2]*x_par[2]))<pipelim[2][0]) {
      //pass 220
      if ((sqrt((y_par[3]*y_par[3])+(x_par[3]*x_par[3]))<pipelim[3][0]) &&
	  (fabs(y_par[3])>detlim[3][1] || x_par[3]>detlim[3][0])) {
	accettanza = 1;
	
	PosEvo.setX(1000*x_par[3]);
	PosEvo.setY(1000*y_par[3]);
	PosEvo.setZ(1000*rp[3]);	  
	if(Pos.z()<vz)PosEvo.setZ(-1000*rp[3]);
      }
    }
    
  }
/*
  LogDebug("ForwardSim") << "\n"
			 << "ACCETTANZA: "<<accettanza << "\n" 
			 << "CSI: "<< csi << "\n"
			 << "Theta_X: " << ThetaX << "\n"
			 << "Theta_Y: " << ThetaY << "\n"
			 << "X_at_0: "<< X_at_0 << "\n"
			 << "Y_at_0: "<< Y_at_0 << "\n" 
			 << "x_par[3]: "<< x_par[3] << "\n"
			 << "y_par[3]: "<< y_par[3] << "\n"
			 << "pos " << Pos.x() << " " << Pos.y() << " " 
			 << Pos.z() << "\n" << "V "<< vx << " " << vy << " "
			 << vz << "\n"
*/
// --------------
  return PosEvo;
}

void TotemSD::StoreHit(TotemG4Hit* hit)
{
  if (hit == 0 )
  {
    if(verbosity_)
      LogDebug("TotemRP") << "TotemSD: hit to be stored is NULL !!" <<std::endl;
    return;
  }

  theHC->insert( hit );
}

void TotemSD::ResetForNewPrimary() {
  
  entrancePoint  = SetToLocal(hitPoint);
  incidentEnergy = preStepPoint->GetKineticEnergy();
}

void TotemSD::Summarize() {
}

bool TotemSD::IsPrimary(const G4Track * track)
{
  TrackInformation* info
    = dynamic_cast<TrackInformation*>( track->GetUserInformation() );
  return info && info->isPrimary();
}
