#ifndef IOMC_ParticleGuns_FlatProtonLogKsiLogTGun_h
#define IOMC_ParticleGuns_FlatProtonLogKsiLogTGun_h

#include <string>
#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include "boost/shared_ptr.hpp"

using namespace edm;
namespace HepMC {
	class GenEvent;
}


class FlatProtonLogKsiLogTGun : public edm::EDProducer
{
  public:
    FlatProtonLogKsiLogTGun(const edm::ParameterSet &);
    virtual ~FlatProtonLogKsiLogTGun();
  
  private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void beginRun(Run const&, EventSetup const&);
    bool ComputeTheta(double ksi, double t, double &theta);  //[-1..0], [GeV], [rad] 
    bool GenerateProton(double &px, double &py, double &pz, double &E);  //[GeV x 4]
  
  protected:
    // data members
    bool right_arm_;
    bool left_arm_;
    double nominalEnergy_;
    double min_t_;
    double max_t_;
    double min_ksi_;
    double max_ksi_;
    double min_phi_;
    double max_phi_;
    int max_iter_number_;
    
    bool log_ksi_distribution_;
    bool log_t_distribution_;
    bool z_symetric_mode_;
    
    double log10_t_min_;
    double log10_t_max_;
    double log10_ksi_min_;
    double log10_ksi_max_;
    
    HepMC::GenEvent* fEvt_;
    edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable;
    int verbosity_;
    CLHEP::HepRandomEngine* fRandomEngine;
    CLHEP::RandFlat* fRandomGenerator;
    
    double m0_;
    double p1_;
    int PartID_;

    std::string instanceLabel;
};


#endif
