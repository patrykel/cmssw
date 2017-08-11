#include "IOMC/ParticleGuns/interface/FlatProtonLogKsiLogTGun.h"
#include <ostream>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <iostream>
#include <TMath.h>



FlatProtonLogKsiLogTGun::FlatProtonLogKsiLogTGun(const edm::ParameterSet& pset)
 :  fEvt_(0)
{
  right_arm_ = pset.getUntrackedParameter<bool>("RightArm");
  left_arm_ = pset.getUntrackedParameter<bool>("LeftArm");
  nominalEnergy_ = pset.getUntrackedParameter<double>("NominalEnergy");
  min_t_ = pset.getUntrackedParameter<double>("MinT");
  max_t_ = pset.getUntrackedParameter<double>("MaxT");
  min_ksi_ = pset.getUntrackedParameter<double>("MinKsi");
  max_ksi_ = pset.getUntrackedParameter<double>("MaxKsi");
  min_phi_ = pset.getUntrackedParameter<double>("MinPhi");
  max_phi_ = pset.getUntrackedParameter<double>("MaxPhi");
  verbosity_ = pset.getUntrackedParameter<int>("Verbosity");
  max_iter_number_ = pset.getUntrackedParameter<int>("MaximumIterationNumber");
  
  log_ksi_distribution_ = pset.getUntrackedParameter<bool>("LogKsiDistShape");
  log_t_distribution_ = pset.getUntrackedParameter<bool>("LogTDistShape");
  z_symetric_mode_ = pset.getUntrackedParameter<bool>("GenerateMirroredProtons");
  instanceLabel = pset.getUntrackedParameter<std::string>("instanceLabel");

  produces<edm::HepMCProduct>();

  log10_t_min_ = TMath::Log10( TMath::Abs(min_t_) );
  log10_t_max_ = TMath::Log10( TMath::Abs(max_t_) );
  log10_ksi_min_ = TMath::Log10( TMath::Abs(min_ksi_) );
  log10_ksi_max_ = TMath::Log10( TMath::Abs(max_ksi_) );
  if(verbosity_){
    std::cout.precision(5);
    std::cout<<" log10_t_min_="<<log10_t_min_<<" log10_t_max_="<<log10_t_max_
    <<" log10_ksi_min_="<<log10_ksi_min_<<" log10_ksi_max_="<<log10_ksi_max_
    <<std::endl;
  }
}

FlatProtonLogKsiLogTGun::~FlatProtonLogKsiLogTGun()
{
  if(fRandomGenerator)
    delete fRandomGenerator;
}


bool FlatProtonLogKsiLogTGun::ComputeTheta(double ksi, double t, double &theta)
{
  long double term1 = 1+ksi;
  long double term2 = term1*term1;
  long double term3 = m0_*m0_;
  long double term4 = p1_*p1_;
  long double term5 = term3*term3;
  long double term6 = ksi*ksi;
  long double sqrt1 = TMath::Sqrt(term3+term4);
  long double sqrt2 = TMath::Sqrt(term3+term2*term4);
  long double denom1 = term2*term4*term4;
  long double denom2 = term2*term4;
  
  long double bracket1 = -2*term5/denom1 - 2*term3/denom2 - 2*ksi*term3/denom2
      - term6*term3/denom2 + 2*term3*sqrt1*sqrt2/denom1;
      
  long double bracket2 = t*(term3/denom1 - t/(4*denom1) - sqrt1*sqrt2/denom1);
  
  long double theta_squared = bracket1 + bracket2;
    
  bool res = false;
  theta = 0.0;
      
  if(theta_squared>=0.0 && theta_squared<0.01)
  {
    theta=TMath::Sqrt(theta_squared);
    res = true;
  }
  
  if(verbosity_)
    std::cout<<" ksi="<<ksi<<" t="<<t<<" theta_squared="<<theta_squared<<std::endl;
  
  res = res || (ksi==0.0 && t==0.0);
  return res;
}

void FlatProtonLogKsiLogTGun::beginRun(Run const & r, EventSetup const& es){
  es.getData(fPDGTable);
  PartID_ = 2212;
  const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(PartID_));
  m0_ = PData->mass().value();
  p1_ = TMath::Sqrt(nominalEnergy_*nominalEnergy_ - m0_*m0_);
}

bool FlatProtonLogKsiLogTGun::GenerateProton(double &px, double &py, double &pz, double &E)
{
  if(verbosity_)
    std::cout<<"FlatProtonLogKsiLogTGun::GenerateProton"<<std::endl;
    
  long double log10t, log10ksi, phi;
  double ksi, t, theta;
  long double p2, E2;
  
  bool result;
  
  int counter = 0;
  do
  { 
    phi = fRandomGenerator->fire(min_phi_, max_phi_);
    if(log_ksi_distribution_)
    {
      log10ksi = fRandomGenerator->fire(log10_ksi_min_, log10_ksi_max_);
      ksi = -TMath::Power(10.0, (double) log10ksi);
    }
    else
      ksi = fRandomGenerator->fire(min_ksi_, max_ksi_);
    
    if(log_t_distribution_)
    {
      log10t = fRandomGenerator->fire(log10_t_min_, log10_t_max_);
      t = -TMath::Power(10.0, (double)log10t);
    }
    else
      t = fRandomGenerator->fire(min_t_, max_t_);
    
    p2 = (1+ksi)*p1_;
    E2 = TMath::Sqrt(p2*p2 + m0_*m0_);
    result = ComputeTheta(ksi, t, theta);
    ++counter;
    if(verbosity_)
    {
      std::cout.precision(5);
      std::cout<<"counter="<<counter<<"  result="<<result<<std::endl;
      std::cout<<" log10t="<<log10t<<" log10ksi="<<log10ksi<<" phi="<<phi
          <<" ksi="<<ksi<<" t="<<t<<" p2="<<p2<<" E2="<<E2<<" beam_momentum="<<p1_<<std::endl;
    }
  }
  while(!result && counter<max_iter_number_);

  px = p2*sin(theta)*cos(phi);
  py = p2*sin(theta)*sin(phi);
  pz = p2*cos(theta);
  E = E2;
  
  if(verbosity_ && result)
  {
    std::cout<<"Particle added: "<<"px:"<<px<<" py:"<<py<<" pz:"<<pz<<" E:"<<E<<std::endl;
    std::cout<<"mad_thx:"<<px/p1_<<" mad_thy:"<<py/p1_<<std::endl;
  }
  
  return result;
}


void FlatProtonLogKsiLogTGun::produce(edm::Event& e, const edm::EventSetup& es)
{
  if(verbosity_)
    std::cout<<" FlatProtonLogKsiLogTGun : Begin New Event Generation"<<std::endl;


  edm::Service<edm::RandomNumberGenerator> rng;
  fRandomEngine = &rng->getEngine(e.streamID());
  fRandomGenerator = new CLHEP::RandFlat(fRandomEngine);

  fEvt_ = new HepMC::GenEvent() ;
  HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));

  int barcode = 1;
  double px1, py1, pz1, E1;
  double px2, py2, pz2, E2;
  int part_count=0;
  bool right_generated = false; 
  bool left_generated = false;
  
  if(z_symetric_mode_ && (right_arm_||left_arm_))
  {
    right_generated = left_generated = 
          GenerateProton(px1, py1, pz1, E1);
    px2=-px1; py2=-py1; pz2=-pz1; E2=E1;
    right_generated = right_generated && right_arm_;
    left_generated = left_generated && left_arm_;
  }
  else
  {
    if(right_arm_)
      right_generated = GenerateProton(px1, py1, pz1, E1);
    if(left_arm_)
    {
      left_generated = GenerateProton(px2, py2, pz2, E2);
      pz2 = -pz2;
    }
  }
  
  if(right_generated)
  {
    HepMC::ThreeVector p(px1,py1,pz1);
    HepMC::GenParticle* Part =
          new HepMC::GenParticle(HepMC::FourVector(px1,py1,pz1,E1),PartID_,1);
    Part->suggest_barcode(barcode);
    barcode++ ;
    Vtx->add_particle_out(Part);
    ++part_count;
  }
  
  if(left_generated)
  {
    HepMC::ThreeVector p(px2,py2,pz2);
    HepMC::GenParticle* Part =
          new HepMC::GenParticle(HepMC::FourVector(px2,py2,pz2,E2),PartID_,1);
    Part->suggest_barcode(barcode);
    barcode++ ;
    Vtx->add_particle_out(Part);
    ++part_count;
  }

  fEvt_->add_vertex(Vtx);
  fEvt_->set_event_number(e.id().event());
  fEvt_->set_signal_process_id(20);

  if(verbosity_)
  {
    fEvt_->print();
  }

  std::unique_ptr<edm::HepMCProduct> BProduct(new edm::HepMCProduct()) ;
  BProduct->addHepMCData(fEvt_);
  e.put(std::move(BProduct));

  if(verbosity_)
  {
    std::cout << "FlatProtonLogKsiLogTGun : Event Generation Done, "<<part_count<< std::endl;
  }
}

DEFINE_FWK_MODULE(FlatProtonLogKsiLogTGun);
