{
	gROOT->Reset();

	#include <string>
	#include <iostream>
	using namespace std;

	#define PLANES_IN_RP 1
	#define RPS_IN_STATION 2

	/////////////////////////////////////////////////////////////////////////
	// DESCRIPTION
	/////////////////////////////////////////////////////////////////////////
	// Input file	- the ROOT file to open
	// Option 
	//		- PLANES_IN_RP 		- on ten different planes of single RP
	//		- RPS_IN_STATION	- on six different RPs from selected station
	// In order to select arm, station, RP, plane you specify:
	// 		- armId			[0,1] 
	//		- stationId 	[0 | 2]
	//		- RomanPotId 	[0,6]
	//		- planeId		[0,9]
	//	
	// usage: 
	// > root
	// > .x [path to the script]
	/////////////////////////////////////////////////////////////////////////

	char in_root_file[]="eliza_2m.root";
	int option = RPS_IN_STATION;

	int c_arId = 0x0;
	int c_stId = 0x0;
	int c_rpId = 0x0;
	int c_plId = 0x0;

	TFile file(in_root_file);
	TTreeReader reader("Events", &file);
	TTreeReaderArray<PSimHit> ps(reader, "PSimHits_g4SimHits_TotemHitsRP_TestFlatGun.obj");


	// script for plotting all planes in given RP
	if(option == PLANES_IN_RP){
		string plNames[10] = {"0","1","2","3","4","5","6","7","8","9"};
		TH2D** h = new TH2D*[10];
		for(int i=0; i<10; ++i) {
			h[i] = new TH2D(plNames[i].c_str(), plNames[i].c_str(), 100, -20, 20, 100, -20, 20);
		}

		while (reader.Next()) {
			for (int i=0; i<ps.GetSize(); i++) {

				int arId = (ps[i].detUnitId()>> 24) & 0x1;
				int stId = (ps[i].detUnitId()>> 22) & 0x3;
				int rpId = (ps[i].detUnitId()>> 19) & 0x7;
				int plId = (ps[i].detUnitId()>> 15) & 0xF;

				if(arId!=c_arId || stId!=c_stId || rpId!=c_rpId) continue;

				double xVal = ps[i].entryPoint().x();
				double yVal = ps[i].entryPoint().y();
				h[plId]->Fill(xVal, yVal);
		  	}
		}

		TCanvas *cc = new TCanvas("Planes","Planes",1500,1000); 
		cc->Divide(5, 2);

		for(int i=0; i<10; ++i) {
			cc->cd(i+1);
			h[i]->Draw("COLZ");
		}
	} else if (option == RPS_IN_STATION) {
		// script for all RPs in given arm, station at given plane
		string rpNames[6] = {"near top", "near bottom", "near horizontal", "far horizontal", "far top", "far bottom"};
		TH2D** h = new TH2D*[6]; 
		for(int i=0; i<6; ++i) {
			h[i] = new TH2D(rpNames[i].c_str(), rpNames[i].c_str(), 100, -20, 20, 100, -20, 20);
		}

		while (reader.Next()) {
			for (int i=0; i<ps.GetSize(); i++) {

				int arId = (ps[i].detUnitId()>> 24) & 0x1;
				int stId = (ps[i].detUnitId()>> 22) & 0x3;
				int rpId = (ps[i].detUnitId()>> 19) & 0x7;
				int plId = (ps[i].detUnitId()>> 15) & 0xF;

				if(arId!=c_arId || stId!=c_stId || plId!=c_plId) continue;

				double xVal = ps[i].entryPoint().x();
				double yVal = ps[i].entryPoint().y();
				h[rpId]->Fill(xVal, yVal);
		  	}
		}

		TCanvas *cc = new TCanvas("RPs","RPs",1500,1000); 
		cc->Divide(3, 2);

		for(int i=0; i<6; ++i) {
			cc->cd(i+1);
			h[i]->Draw("COLZ");
		}
	}

}


