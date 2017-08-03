#ifndef IOMC_GaussEvtVtxEnergyGenerator_H
#define IOMC_GaussEvtVtxEnergyGenerator_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "IOMC/EventVertexGenerators/interface/GaussEvtVtxGenerator.h"
#include "HepMC/GenEvent.h"



class GaussEvtVtxEnergyGenerator : public edm::EDProducer{
public:
    GaussEvtVtxEnergyGenerator(const edm::ParameterSet &p);
    virtual ~GaussEvtVtxEnergyGenerator();
    virtual void produce( edm::Event&, const edm::EventSetup&);
    void applyBeamSmearing(const HepMC::GenEvent *evt, CLHEP::HepRandomEngine* engine);
    void beginRun(edm::Run const&iR, edm::EventSetup const&es);
    HepMC::FourVector* newVertex(CLHEP::HepRandomEngine* engine);

private:

    double MeanXi, SiXi;                            ///< Xi distribution
    double Al;                                      ///< Alpha (half crossing angle) d.
    double SiThX, SiThY;                            ///< Theta (beam divergence) d.

    double E_CMS;

    double MeanX, MeanY, MeanZ;                     ///< [mm] means
    double SiX, SiY, SiZ;                           ///< [mm] variance

    edm::InputTag sourceLabel;
    unsigned int verbosity;
};

#endif
