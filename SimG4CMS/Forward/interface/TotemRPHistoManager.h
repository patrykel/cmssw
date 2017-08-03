#ifndef SimG4CMS_Forward_TotemRPHistoManager_H
#define SimG4CMS_Forward_TotemRPHistoManager_H 1
// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemRPHistoManager
//
/**
 Description: Manages Root file creation for Totem RP Tests
 
 Usage:
    Used in testing Totem RP simulation
 
*/
 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"

#include <string>

class TotemRPHistoClass;

class TotemRPHistoManager {

public: 

  TotemRPHistoManager(const std::string &);
  virtual ~TotemRPHistoManager();

  void fillTree(TotemRPHistoClass *  histos);

private:
  edm::Service<TFileService> fs;
  TTree                     *tree;
  TotemRPHistoClass         *h;
  int                       kount;
};

#endif  //TotemRP_TotemRPHistoManager_h
