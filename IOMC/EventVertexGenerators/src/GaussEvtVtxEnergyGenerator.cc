#include "IOMC/EventVertexGenerators/interface/GaussEvtVtxEnergyGenerator.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "HepMC/GenEvent.h"

#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "HepMC/SimpleVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "TLorentzVector.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"

#include <iostream>

using namespace std;

GaussEvtVtxEnergyGenerator::GaussEvtVtxEnergyGenerator(const edm::ParameterSet & p ):
        sourceLabel(p.getParameter<edm::InputTag>("src")),
        verbosity(p.getUntrackedParameter<unsigned int>("verbosity", 1))
{
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
        throw cms::Exception("Configuration")
            << "The BaseEvtVtxGenerator requires the RandomNumberGeneratorService\n"
            "which is not present in the configuration file. \n"
            "You must add the service\n"
            "in the configuration file or remove the modules that require it.";
    }

    consumes<edm::HepMCProduct>(sourceLabel);
    produces<edm::HepMCProduct>();
}

GaussEvtVtxEnergyGenerator::~GaussEvtVtxEnergyGenerator()
{
}

void GaussEvtVtxEnergyGenerator::produce(edm::Event& event, const edm::EventSetup& es)
{
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(event.streamID());

    edm::Handle<edm::HepMCProduct> HepUnsmearedMCEvt;

    event.getByLabel( sourceLabel, HepUnsmearedMCEvt );

    // Copy the HepMC::GenEvent
    HepMC::GenEvent* genevt = new HepMC::GenEvent(*HepUnsmearedMCEvt->GetEvent());
    std::unique_ptr<edm::HepMCProduct> HepMCEvt(new edm::HepMCProduct(genevt));

    //apply beam smearing
    applyBeamSmearing(HepMCEvt->GetEvent(), engine);

    //apply vertex smearing
    HepMCEvt->applyVtxGen( newVertex(engine) ) ;

    event.put(std::move(HepMCEvt)) ;

    return;
}

void GaussEvtVtxEnergyGenerator::beginRun(edm::Run const&iR, edm::EventSetup const&es) {
    edm::ESHandle<BeamOpticsParams> par;
    es.get<BeamOpticsParamsRcd>().get(par);

    if(!par.isValid()){
        cout<<"<<SmearingGenerator:: Invalid BeamOpticsParams"<<endl;
        throw cms::Exception("SmearingGenerator") << " edm::ESHandle<BeamOpticsParams> is invalid";
    }
    // angles in rad and positions in m
    Al = par->GetCrossingAngleX();
    SiThX = par->GetBeamDivergenceX();
    SiThY = par->GetBeamDivergenceY();
    MeanXi = par->GetMeanXi();
    SiXi = par->GetSigmaXi();
    E_CMS = par->GetBeamEnergy();
    MeanX = par->GetBeamDisplacementX() * 1000.;
    SiX = par->GetPrimVertSizeX() * 1000.;
    MeanY = par->GetBeamDisplacementY() * 1000.;
    SiY = par->GetPrimVertSizeY() * 1000.;
    MeanZ = par->GetBeamDisplacementZ() * 1000.;
    SiZ = par->GetPrimVertSizeZ() * 1000.;

    if (verbosity) {
        printf(">> SmearingGenerator::beginJob\n");
        printf("Al = %E\n", Al);
        printf("SiThX = %E\n", SiThX);
        printf("SiThY = %E\n", SiThY);
        printf("MeanXi = %E\n", MeanXi);
        printf("SiXi = %E\n", SiXi);
        printf("E_CMS = %E GeV\n", E_CMS);
        printf("MeanX = %E mm\n", MeanX);
        printf("SiX = %E mm\n", SiX);
        printf("MeanY = %E mm\n", MeanY);
        printf("SiY = %E mm\n", SiY);
        printf("MeanZ = %E mm\n", MeanZ);
        printf("SiZ = %E mm\n", SiZ);
    }
}

HepMC::FourVector* GaussEvtVtxEnergyGenerator::newVertex(CLHEP::HepRandomEngine* engine) {
    double X,Y,Z,T;
    X = CLHEP::RandGaussQ::shoot(engine, MeanX, SiX);
    Y = CLHEP::RandGaussQ::shoot(engine, MeanY, SiY);
    Z = CLHEP::RandGaussQ::shoot(engine, MeanZ, SiZ);
    T = 0.;

    HepMC::FourVector* fVertex = new HepMC::FourVector() ;
    fVertex->set( X, Y, Z, T);

    return fVertex;
}

void GaussEvtVtxEnergyGenerator::applyBeamSmearing(const HepMC::GenEvent *evt, CLHEP::HepRandomEngine* engine)
{
    // generate energy/angular smearing
    double al1 = Al;
    double al2 = -Al;
    double thx1 = CLHEP::RandGaussQ::shoot(engine, 0., SiThX);
    double thy1 = CLHEP::RandGaussQ::shoot(engine, 0., SiThY);
    double thx2 = CLHEP::RandGaussQ::shoot(engine, 0., SiThX);
    double thy2 = CLHEP::RandGaussQ::shoot(engine, 0., SiThY);
    double xi1 = CLHEP::RandGaussQ::shoot(engine, MeanXi, SiXi);
    double xi2 = CLHEP::RandGaussQ::shoot(engine, MeanXi, SiXi);

    // compute transform parameters
    double m = 0.938271; //proton mass
    double p_nom = sqrt(E_CMS * E_CMS - m*m);

    double thz1 = sqrt(1. - thx1*thx1 - thy1*thy1);
    double thz2 = sqrt(1. - thx2*thx2 - thy2*thy2);

    TVector3 p1(cos(al1)*thx1 + sin(al1)*thz1, thy1, -sin(al1)*thx1 + cos(al1)*thz1);
    p1 *= (1. + xi1) * p_nom;
    double E1 = sqrt(p1.Mag()*p1.Mag() + m*m);
    TLorentzVector P1(p1.x(), p1.y(), p1.z(), E1);

    TVector3 p2(cos(al2)*thx2 + sin(al2)*thz2, thy2, -sin(al2)*thx2 + cos(al2)*thz2);
    p2 *= -(1. + xi2) * p_nom;
    double E2 = sqrt(p2.Mag()*p2.Mag() + m*m);
    TLorentzVector P2(p2.X(), p2.Y(), p2.Z(), E2);

    double factor = (P1 + P2).Mag() / 2. / E_CMS;             // energy correction factor

    TVector3 boost = (p1 + p2) * (1. / (E1 + E2));            // boost vector (direction * beta)
    double beta = (p1 + p2).Mag() / (E1 + E2);                // beta of the boost
    P1.Boost(-boost);
    TVector3 axis(P1.Y(), -P1.X(), 0.);                       // rotation axis
    double angle = -acos( P1.Z() / P1.Vect().Mag() );         // angle of rotation

    if (verbosity > 5) {
        cout << ">> SmearingGenerator::ApplyMomentumSmearing" << endl;
    }

    // apply transform to all particles
    for (HepMC::GenEvent::particle_const_iterator particle = evt->particles_begin(); particle != evt->particles_end(); ++particle) {
        HepMC::FourVector P = (*particle)->momentum();

        if (verbosity > 5)
            printf("\tbefore: (%+.0f| %+4.2f; %+4.2f; %+.0f)\n", P.t(), P.x(), P.y(), P.z());

        // energy scaling
        TVector3 p(P.x(), P.y(), P.z());
        double E = P.e() * factor, m = P.m();
        p = sqrt(E*E - m*m) / p.Mag() * p;
        TLorentzVector PP(p, E);

        // rotation
        if (fabs(angle) > 1E-8)
            PP.Rotate(angle, axis);

        // boost
        if (fabs(beta) > 1E-8)
            PP.Boost(boost);

        if (verbosity > 5)
            printf("\t after: (%+.0f| %+4.2f; %+4.2f; %+.0f), th_x = %.2E, th_y = %.2E\n\n",
                   PP.T(), PP.X(), PP.Y(), PP.Z(), PP.X()/PP.Z(), PP.Y()/PP.Z());

        (*particle)->set_momentum(HepMC::FourVector(PP.X(), PP.Y(), PP.Z(), PP.T()));
    }
}

