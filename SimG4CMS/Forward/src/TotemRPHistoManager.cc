// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemRPHistoManager
//
// Implementation:
//     <Notes on implementation>
//

// system include files
#include <iostream>
#include <cmath>

// user include files
#include "SimG4CMS/Forward/interface/TotemRPHistoManager.h"
#include "SimG4CMS/Forward/interface/TotemRPHistoClass.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/PluginManager/interface/PluginManager.h"

//
// constructors and destructor
//

TotemRPHistoManager::TotemRPHistoManager(const std::string & file):
	tree(0), h(0), kount(0) {
	if (fs.isAvailable()) {
		h    = new TotemRPHistoClass();

		tree = fs->make<TTree>("TotemRP", "TotemRP");
		tree->SetAutoSave(10000);
		tree->Branch("TotemRPHisto", "TotemRPHistoClass", &h);
		edm::LogInfo("TotemRP") << "TotemRPHistoManager:===>>>  Book the Tree";
	} else {
		edm::LogInfo("TotemRP") << "TotemRPHistoManager:===>>> No file provided";
	}
}

TotemRPHistoManager::~TotemRPHistoManager() {

	edm::LogInfo("TotemRP") << "============================================="
	                        << "========================================\n"
	                        << "=== TotemRPHistoManager: Start writing user "
	                        << "histograms after " << kount << " events ";
	if (h) delete h;
}

void TotemRPHistoManager::fillTree(TotemRPHistoClass *  histos) {

	kount++;
	LogDebug("TotemRP") << "TotemRPHistoManager: tree pointer for " << kount
	                    << " = " << histos;
	if (tree) {
	    h = histos;
	    tree->Fill();
	}
}
