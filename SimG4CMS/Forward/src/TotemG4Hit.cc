// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemG4Hit
//
// Implementation:
//     <Notes on implementation>
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
//

// system include files

// user include files
#include "SimG4CMS/Forward/interface/TotemG4Hit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

TotemG4Hit::TotemG4Hit():entry(0)
{
  theIncidentEnergy = 0.0;
  theTrackID = -1;
  theUnitID  =  0;
  theTimeSlice = 0.0;


  theX              = 0.;
  theY              = 0.;
  theZ              = 0.;
  thePabs           = 0.;
  theTof            = 0.;
  theEnergyLoss     = 0.;
  theParticleType   = 0;
  theParentId=0;
  theVx = 0.0;
  theVy = 0.0;
  theVz = 0.0;
  p_x = p_y = p_z = 0.0;
}

TotemG4Hit::~TotemG4Hit() {}

TotemG4Hit::TotemG4Hit(const TotemG4Hit &right)
{
  theIncidentEnergy  = right.theIncidentEnergy;
  theTrackID = right.theTrackID;
  theUnitID = right.theUnitID;
  theTimeSlice = right.theTimeSlice;
  entry    = right.entry;

  thePabs =right.thePabs;
  theTof=right.theTof ;
  theEnergyLoss=right.theEnergyLoss   ;
  theParticleType=right.theParticleType ;
  theX = right.theX;
  theY = right.theY;
  theZ = right.theZ;

  theVx = right.theVx;
  theVy = right.theVy;
  theVz = right.theVz;

  theParentId = right.theParentId;
}


const TotemG4Hit& TotemG4Hit::operator=(const TotemG4Hit &right)
{

  entry             = right.entry;
  theIncidentEnergy = right.theIncidentEnergy;
  theTrackID        = right.theTrackID;
  theUnitID         = right.theUnitID;
  theTimeSlice      = right.theTimeSlice;

  theX              = right.theX;
  theY              = right.theY;
  theZ              = right.theZ;
  thePabs           = right.thePabs;
  theTof            = right.theTof ;
  theEnergyLoss     = right.theEnergyLoss   ;
  theParticleType   = right.theParticleType ;


  theVx = right.theVx;
  theVy = right.theVy;
  theVz = right.theVz;

  theParentId = right.theParentId;

  return *this;
}

void TotemG4Hit::Print() {
  edm::LogInfo("TotemRP") << (*this);
}


Hep3Vector TotemG4Hit::getEntry() const {return entry;}
void TotemG4Hit::setEntry(Hep3Vector xyz) { entry = xyz; }

Hep3Vector TotemG4Hit::getExit() const {return exit;}
void TotemG4Hit::setExit(Hep3Vector xyz) { exit = xyz; }

Hep3Vector TotemG4Hit::getLocalEntry() const {return local_entry;}
void TotemG4Hit::setLocalEntry(const Hep3Vector &xyz) { local_entry = xyz;}
Hep3Vector TotemG4Hit::getLocalExit() const {return local_exit;}
void TotemG4Hit::setLocalExit(const Hep3Vector &xyz) { local_exit = xyz;}

double     TotemG4Hit::getIncidentEnergy() const  {return theIncidentEnergy; }
void       TotemG4Hit::setIncidentEnergy(double e) {theIncidentEnergy  = e; }

unsigned int TotemG4Hit::getTrackID() const {return theTrackID; }
void       TotemG4Hit::setTrackID (int i)         { theTrackID = i; }

int TotemG4Hit::getUnitID() const {return theUnitID; }
void TotemG4Hit::setUnitID (unsigned int i) { theUnitID = i; }

double     TotemG4Hit::getTimeSlice() const       {return theTimeSlice; }
void       TotemG4Hit::setTimeSlice (double d)    { theTimeSlice = d; }
int        TotemG4Hit::getTimeSliceID() const     {return (int)theTimeSlice;}

double TotemG4Hit::getPabs() const {return thePabs;}
double TotemG4Hit::getTof() const {return theTof;}
double TotemG4Hit::getEnergyLoss() const {return theEnergyLoss;}
int TotemG4Hit::getParticleType() const {return theParticleType;}

void TotemG4Hit::setPabs(double e) {thePabs = e;}
void TotemG4Hit::setTof(double e) {theTof = e;}
void TotemG4Hit::setEnergyLoss(double e) {theEnergyLoss = e;}
void TotemG4Hit::addEnergyLoss(double e) {theEnergyLoss += e;}
void TotemG4Hit::setParticleType(short i) {theParticleType = i;}

double TotemG4Hit::getThetaAtEntry() const {return theThetaAtEntry;}
double TotemG4Hit::getPhiAtEntry() const{ return thePhiAtEntry;}

void TotemG4Hit::setThetaAtEntry(double t) {theThetaAtEntry = t;}
void TotemG4Hit::setPhiAtEntry(double f) {thePhiAtEntry = f ;}

double TotemG4Hit::getX() const{ return theX;}
void TotemG4Hit::setX(double t){theX = t;}

double TotemG4Hit::getY() const{ return theY;}
void TotemG4Hit::setY(double t){theY = t;}

double TotemG4Hit::getZ() const{ return theZ;}
void TotemG4Hit::setZ(double t){theZ = t;}

int TotemG4Hit::getParentId() const {return theParentId;}
void TotemG4Hit::setParentId(int p){theParentId = p;}

double TotemG4Hit::getVx() const{ return theVx;}
void TotemG4Hit::setVx(double t){theVx = t;}

double TotemG4Hit::getVy() const{ return theVy;}
void TotemG4Hit::setVy(double t){theVy = t;}

double TotemG4Hit::getVz() const{ return theVz;}
void TotemG4Hit::setVz(double t){theVz = t;}

void TotemG4Hit::set_p_x(double p) {p_x = p;}
void TotemG4Hit::set_p_y(double p) {p_y = p;}
void TotemG4Hit::set_p_z(double p) {p_z = p;}

double TotemG4Hit::get_p_x() const {return p_x;}
double TotemG4Hit::get_p_y() const {return p_y;}
double TotemG4Hit::get_p_z() const {return p_z;}


std::ostream& operator<<(std::ostream& os, const TotemG4Hit& hit) {
  os << " Data of this TotemG4Hit are:" << std::endl
     << " Time slice ID: " << hit.getTimeSliceID() << std::endl
     << " EnergyDeposit = " << hit.getEnergyLoss() << std::endl
     << " Energy of primary particle (ID = " << hit.getTrackID()
     << ") = " << hit.getIncidentEnergy() << " (MeV)"<<std::endl
     << " Entry point in Totem unit number " << hit.getUnitID()
     << " is: " << hit.getEntry() << " (mm)" << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
     << std::endl;
  return os;
}


