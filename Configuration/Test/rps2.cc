{
gROOT->Reset();

#include <string>
#include <iostream>
#include <math.h>
using namespace std;

// draw 2D histograms showing the events on six different RPs
//
// > root
// > .x [path to the script]

// ========= PARAMS =============

// mode:
// 	0 - draw all 6 RPs separately
// 	1 - draw plots for near/far RPs separately
// 	2 - draw plots for near/far horizontal RPs separately
const int mode = 0;

// filename:
TFile file("eliza_2m.root");

// select the RPs to draw - leave (-1) to draw all of them
const int arm = 1;
const int station = 2;
// in this script we only use station 0 - the others were not used during our run
const int romanPot = -1;
const int plane = 1; 

// ========= ENDPARAMS =========

// =========== CONST ===========
// first 6 RPs for right arm, then for the left
// distances taken from RP_Dist_Beam_Cent.xml plus 11mm of th distance from the detector cut edge to its center
const double xt[12] = { 0.0, 0.0, 17.605, 16.319, 0.0, 0.0, 0.0, 0.0, 17.480, 17.181, 0.0, 0.0 };
const double yt[12] = { 20.052, -19.017, 0.0, 0.0, 18.553, -18.363, 20.129, -19.885, 0.0, 0.0, 19.043, -19.320 };
double degree[12] = { 135.0, -45.0, 45.0, 45.0, 135.0, -45.0, -45.0, 45.0, -45.0, -45.0, -135.0, 45.0 };	// todo the tilted ones
for(int i=0; i<12; ++i) degree[i] = degree[i] * 3.14159 / 180;
string rpNames0[6] = {"near top", "near bottom", "near horizontal", "far horizontal", "far top", "far bottom"};
if(mode == 1){
	rpNames0[0] = "near";
	rpNames0[1] = "far";
}
string rpNames1[6] = {"", "", "", "", "", ""};	//todo
const int divX = mode==0 ? 3 : 2;
const int divY = mode==0 ? 2 : 1;
const int RPCount = divX*divY;
// =========== ENDCONST ========






TTreeReader reader("Events", &file);
TTreeReaderArray<PSimHit> ps(reader, "PSimHits_g4SimHits_TotemHitsRP_TestFlatGun.obj");

TH2D** h = new TH2D*[RPCount]; 

for(int i=0; i<RPCount; ++i) {
	h[i] = new TH2D(rpNames0[i].c_str(), rpNames0[i].c_str(), 1000, -50, 50, 1000, -50, 50);
//	h[i]->GetXaxis()->SetTitle("x [mm]");
//	h[i]->GetYaxis()->SetTitle("y [mm]");
}

while (reader.Next()) {
	for (int i=0; i<ps.GetSize(); i++) {

		int arId = (ps[i].detUnitId()>> 24) & 0x1;
		int stId = (ps[i].detUnitId()>> 22) & 0x3;
		int rpId = (ps[i].detUnitId()>> 19) & 0x7;
		int plId = (ps[i].detUnitId()>> 15) & 0xF;

		if(mode==2 && rpId!=2 && rpId!=3) continue;
		if((arm>=0 && arId!=arm) || (stId!=station) || (romanPot>=0 && rpId!=romanPot) || (plane>=0 && plId!=plane)) continue;

		double xVal = ps[i].entryPoint().x();
		double yVal = ps[i].entryPoint().y();
	
		if(mode==1) {
			int rpPosId = rpId;//6*arId+rpId;
			double xVal2 = xVal;
			double yVal2 = yVal;
			xVal = xVal2*cos(degree[rpPosId]) + yVal2*sin(degree[rpPosId]) + xt[rpPosId];
			yVal = xVal2*sin(degree[rpPosId]) - yVal2*cos(degree[rpPosId]) + yt[rpPosId];
		}

		int graphId = mode==0 ? rpId : (mode==1 ? rpId/3 : rpId-2);
		h[graphId]->Fill(xVal, yVal);
  	}
}

TCanvas *cc = new TCanvas("RPs", "RPs", 1800, 1800*divY/divX); 
cc->Divide(divX, divY);

for(int i=0; i<RPCount; ++i) {
	cc->cd(i+1);
	h[i]->SetStats(false);
	h[i]->Draw("COLZ");
}



}

