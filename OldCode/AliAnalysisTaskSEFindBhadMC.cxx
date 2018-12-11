/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskSEFindBhadMC.cxx 53668 2011-12-16 14:14:00Z prino $ */

/////////////////////////////////////////////////////////////
// Authors: J.Stiller
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEFindBhadMC.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliCFVertexingHF.h"
#include "AliCFVertexingHF2Prong.h"

#include "AliGenHijingEventHeader.h"


ClassImp(AliAnalysisTaskSEFindBhadMC)


//________________________________________________________________________
AliAnalysisTaskSEFindBhadMC::AliAnalysisTaskSEFindBhadMC():
AliAnalysisTaskSE(),
fCuts(0),
fhijingGenHeader(0),
fOutputList(0),
fComparList(0),
fMClimAcc(0),
fFiducialCut(0),
fNentries(0),
fMCTruth_Bplus(0),
fMCTruth_Bminus(0),
fMCTruth2(0),
hEtaD0(0),
hEtaPi(0),
hEtaBoth(0),
hEtaBothiR(0),
hPtD0(0),
hPtPi(0),
hPtVsEtaD0(0),
hPtVsEtaPi(0),
hSgn_CosThetaS_0(0),
hSgn_CosThetaS_1(0),
hSgn_CosThetaS_2(0),
hSgn_CosThetaS_3(0),
hSgn_CosThetaS_4(0),
hSgn_CosThetaS_5(0),
hSgn_CosThetaS_6(0),
hSgn_CosThetaS_7(0),
hSgn_CosThetaS_8(0),
hSgn_CosThetaS_9(0),
hSgn_CosThetaS_10(0),
hSgn_CosThetaS_11(0),
hSgn_CosThetaS_12(0),
hSgn_CosThetaS_D_0(0),
hSgn_CosThetaS_D_1(0),
hSgn_CosThetaS_D_2(0),
hSgn_CosThetaS_D_3(0),
hSgn_CosThetaS_D_4(0),
hSgn_CosThetaS_D_5(0),
hSgn_CosThetaS_D_6(0),
hSgn_CosThetaS_D_7(0),
hSgn_CosThetaS_D_8(0),
hSgn_CosThetaS_D_9(0),
hSgn_CosThetaS_D_10(0),
hSgn_CosThetaS_D_11(0),
hSgn_CosThetaS_D_12(0),
hThetaD0(0),
hThetaPi(0),
hThetaD0Ka(0),
hThetaD0Pi(0)
{
	// Default constructor	
}
//________________________________________________________________________
AliAnalysisTaskSEFindBhadMC::AliAnalysisTaskSEFindBhadMC(const char *name,AliRDHFCutsD0toKpi* cuts):
AliAnalysisTaskSE(name),
fCuts(0),
fhijingGenHeader(0),
fOutputList(0),
fComparList(0),
fMClimAcc(0),
fFiducialCut(0),
fNentries(0),
fMCTruth_Bplus(0),
fMCTruth_Bminus(0),
fMCTruth2(0),
hEtaD0(0),
hEtaPi(0),
hEtaBoth(0),
hEtaBothiR(0),
hPtD0(0),
hPtPi(0),
hPtVsEtaD0(0),
hPtVsEtaPi(0),
hSgn_CosThetaS_0(0),
hSgn_CosThetaS_1(0),
hSgn_CosThetaS_2(0),
hSgn_CosThetaS_3(0),
hSgn_CosThetaS_4(0),
hSgn_CosThetaS_5(0),
hSgn_CosThetaS_6(0),
hSgn_CosThetaS_7(0),
hSgn_CosThetaS_8(0),
hSgn_CosThetaS_9(0),
hSgn_CosThetaS_10(0),
hSgn_CosThetaS_11(0),
hSgn_CosThetaS_12(0),
hSgn_CosThetaS_D_0(0),
hSgn_CosThetaS_D_1(0),
hSgn_CosThetaS_D_2(0),
hSgn_CosThetaS_D_3(0),
hSgn_CosThetaS_D_4(0),
hSgn_CosThetaS_D_5(0),
hSgn_CosThetaS_D_6(0),
hSgn_CosThetaS_D_7(0),
hSgn_CosThetaS_D_8(0),
hSgn_CosThetaS_D_9(0),
hSgn_CosThetaS_D_10(0),
hSgn_CosThetaS_D_11(0),
hSgn_CosThetaS_D_12(0),
hThetaD0(0),
hThetaPi(0),
hThetaD0Ka(0),
hThetaD0Pi(0)
{
	// Constructor	
	
	fCuts=cuts;
	
	DefineOutput(1,TList::Class());  //My private output
	DefineOutput(2,TList::Class());  //My private output
	
}

//________________________________________________________________________
AliAnalysisTaskSEFindBhadMC::~AliAnalysisTaskSEFindBhadMC()
{
	
	if (fOutputList) {
		delete fOutputList;
		fOutputList = 0;
	}
	if (fComparList) {
		delete fComparList;
		fComparList = 0;
	}
}  

//________________________________________________________________________
void AliAnalysisTaskSEFindBhadMC::Init()
{
	// Initialization
}

//________________________________________________________________________
void AliAnalysisTaskSEFindBhadMC::UserCreateOutputObjects()
{
	
	// Create the output container
	//
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass::UserCreateOutputObjects() \n");
	
	fOutputList = new TList();
	fOutputList->SetOwner();
	fComparList = new TList();
	fComparList->SetOwner();
	
	fNentries=new TH1I("MCHistoForChecks", "MCHistoForChecks",40,0,40);
	fNentries->GetXaxis()->SetBinLabel(1,"B injected");
	fNentries->GetXaxis()->SetBinLabel(2,"B from b");
	fNentries->GetXaxis()->SetBinLabel(3,"B(inj):B+");
	fNentries->GetXaxis()->SetBinLabel(4,"B(inj):B-");
	fNentries->GetXaxis()->SetBinLabel(5,"B(inj):!2dghts");
	fNentries->GetXaxis()->SetBinLabel(6,"B(inj):!D");
	fNentries->GetXaxis()->SetBinLabel(7,"B(inj):211 cut");
	fNentries->GetXaxis()->SetBinLabel(8,"B(inj):!211");	
	fNentries->GetXaxis()->SetBinLabel(9,"B(inj):D!->2dghts");
	fNentries->GetXaxis()->SetBinLabel(10,"B(inj):broken dghts");
	fNentries->GetXaxis()->SetBinLabel(11,"B(inj):!Ddghts cut");
	fNentries->GetXaxis()->SetBinLabel(12,"B(inj):Remaining B cand.");
	fNentries->GetXaxis()->SetBinLabel(13,"B(b):B+");
	fNentries->GetXaxis()->SetBinLabel(14,"B(b):B-");
	fNentries->GetXaxis()->SetBinLabel(15,"B(b):!2dghts");
	fNentries->GetXaxis()->SetBinLabel(16,"B(b):!D");
	fNentries->GetXaxis()->SetBinLabel(17,"B(b):211 cut");
	fNentries->GetXaxis()->SetBinLabel(18,"B(b):!211");	
	fNentries->GetXaxis()->SetBinLabel(19,"B(b):D!->2dghts");
	fNentries->GetXaxis()->SetBinLabel(20,"B(b):broken dghts");
	fNentries->GetXaxis()->SetBinLabel(21,"B(b):!Ddghts cut");
	fNentries->GetXaxis()->SetBinLabel(22,"B(b):Remaining B cand.");
	
	fNentries->GetXaxis()->SetBinLabel(30,"B daughter 0 not in theta");
	fNentries->GetXaxis()->SetBinLabel(31,"B daughter 1 not in theta");
	
	fMCTruth_Bplus = new TH1F("fMCTruth_Bplus","Bplus",32,0,32);
	
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(1,"Events analyzed");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(2,"B found");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(3,"B -> !2 daughters");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(4,"B -> 2 daughters");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(5,"1st daughter is D");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(6,"2nd daughter is D");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(7,"no daughter is D");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(8,"2nd daughter is Pi");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(9,"1st daughter is Pi");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(10,"no daughter is Pi");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(11,"D -> !2 daughters");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(12,"D -> 2 daughters");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(13,"B+ -> D+Pi -> Ka+Pi");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(14,"B- -> D+Pi -> Ka+Pi");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(15,"id B+ daughters in eta");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(16,"id B- daughters in eta");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(17,"0.0<pt<0.5 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(18,"0.5<pt<1.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(19,"1.0<pt<2.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(20,"2.0<pt<3.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(21,"3.0<pt<4.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(22,"4.0<pt<5.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(23,"5.0<pt<6.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(24,"6.0<pt<8.0 ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(25,"8.0<pt<12. ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(26,"12.<pt<16. ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(27,"16.<pt<20. ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(28,"20.<pt<24. ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(29,"24.<pt<9999. ");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(30,"id B+ from b");
	fMCTruth_Bplus->GetXaxis()->SetBinLabel(31,"id B+ by hand");
	
	fMCTruth_Bminus = new TH1F("fMCTruth_Bminus","Bminus",31,0,31);
	
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(1,"Events analyzed");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(2,"B found");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(3,"B -> !2 daughters");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(4,"B -> 2 daughters");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(5,"1st daughter is D");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(6,"2nd daughter is D");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(7,"no daughter is D");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(8,"2nd daughter is Pi");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(9,"1st daughter is Pi");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(10,"no daughter is Pi");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(11,"D -> !2 daughters");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(12,"D -> 2 daughters");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(13,"B+ -> D+Pi -> Ka+Pi");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(14,"B- -> D+Pi -> Ka+Pi");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(15,"id B+ daughters in eta");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(16,"id B- daughters in eta");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(17,"0.0<pt<0.5 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(18,"0.5<pt<1.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(19,"1.0<pt<2.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(20,"2.0<pt<3.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(21,"3.0<pt<4.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(22,"4.0<pt<5.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(23,"5.0<pt<6.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(24,"6.0<pt<8.0 ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(25,"8.0<pt<12. ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(26,"12.<pt<16. ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(27,"16.<pt<20. ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(28,"20.<pt<24. ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(29,"24.<pt<9999. ");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(30,"id B- from b");
	fMCTruth_Bminus->GetXaxis()->SetBinLabel(31,"id B- by hand");
	
	
	fMCTruth2 = new TH1F("fMCTruth2","daughters",20000,-10000,10000);
	hEtaD0 = new TH1F("hEtaD0","hEta of D0",400,-2,2);
	hEtaPi = new TH1F("hEtaPi","hEta of Pi",400,-2,2);
	hEtaBoth = new TH1F("hEtaBoth","Rapidity y of B",400,-2,2);
	hEtaBothiR = new TH2F("hEtaBothiR","hEta of D0 vs Pi",160,-0.8,0.8,160,-0.8,0.8);
	hPtD0 = new TH1F("hPtD0","pt of D0",20000,0,200);
	hPtPi = new TH1F("hPtPi","pt of Pi",20000,0,200);
	hPtVsEtaD0 = new TH2F("hPtVsEtaD0","D0: Pt vs Eta",2000,0,200,400,-2,2);
	hPtVsEtaPi = new TH2F("hPtVsEtaPi","Pion: Pt vs Eta",2000,0,200,400,-2,2);
	hSgn_CosThetaS_0  = new TH1F("hSgn_CosThetaS_0","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_1  = new TH1F("hSgn_CosThetaS_1","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_2  = new TH1F("hSgn_CosThetaS_2","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_3  = new TH1F("hSgn_CosThetaS_3","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_4  = new TH1F("hSgn_CosThetaS_4","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_5  = new TH1F("hSgn_CosThetaS_5","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_6  = new TH1F("hSgn_CosThetaS_6","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_7  = new TH1F("hSgn_CosThetaS_7","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_8  = new TH1F("hSgn_CosThetaS_8","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_9  = new TH1F("hSgn_CosThetaS_9","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_10 = new TH1F("hSgn_CosThetaS_10","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_11 = new TH1F("hSgn_CosThetaS_11","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_12 = new TH1F("hSgn_CosThetaS_12","Cos(#theta^{*})",2000,-1,1);

	hSgn_CosThetaS_D_0  = new TH1F("hSgn_CosThetaS_D_0","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_1  = new TH1F("hSgn_CosThetaS_D_1","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_2  = new TH1F("hSgn_CosThetaS_D_2","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_3  = new TH1F("hSgn_CosThetaS_D_3","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_4  = new TH1F("hSgn_CosThetaS_D_4","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_5  = new TH1F("hSgn_CosThetaS_D_5","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_6  = new TH1F("hSgn_CosThetaS_D_6","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_7  = new TH1F("hSgn_CosThetaS_D_7","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_8  = new TH1F("hSgn_CosThetaS_D_8","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_9  = new TH1F("hSgn_CosThetaS_D_9","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_10 = new TH1F("hSgn_CosThetaS_D_10","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_11 = new TH1F("hSgn_CosThetaS_D_11","Cos(#theta^{*})",2000,-1,1);
	hSgn_CosThetaS_D_12 = new TH1F("hSgn_CosThetaS_D_12","Cos(#theta^{*})",2000,-1,1);

	hThetaD0 = new TH1F("hThetaD0","D0 from B #theta",3200,0,3.2);
	hThetaPi = new TH1F("hThetaPi","#pi from B #theta",3200,0,3.2);
	hThetaD0Ka = new TH1F("hThetaD0Ka","Ka from D #theta",3200,0,3.2);
	hThetaD0Pi = new TH1F("hThetaD0Pi","#pi from D #theta",3200,0,3.2);

	TH1F *hPdgInvMassB = new TH1F("hPdgInvMassB","B;PDG inv. mass;counts",1000,0,10);
	TH1F *hCalcInvMassB = new TH1F("hCalcInvMassB","B;MC calc. inv. mass;counts",1000,0,10);
	TH1F *hPdgInvMassD = new TH1F("hPdgInvMassD","D;PDG inv. mass;counts",1000,0,10);
	TH1F *hCalcInvMassD = new TH1F("hCalcInvMassD","D;MC calc. inv. mass;counts",1000,0,10);
	
	Int_t NBinsPt = 13;   // number of bins
	Double_t BinsPt[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 24.0, 40.0};  // bin limits
	TH1F* hSgn_CountsY_0[16];
	TH1F* hSgn_CountsY_PtBins_0[16];
		
	for(Int_t i=1; i<16;i++){
		hSgn_CountsY_0[i] = new TH1F(Form("hSgn_CountsY_PtBins_0[%i]",i),Form("B counts within |y|<%i.%i",i/10,i%10),100,0,40);
		hSgn_CountsY_PtBins_0[i] = (TH1F*)hSgn_CountsY_0[i]->Rebin(NBinsPt,"",BinsPt);
	}
	
	TH1F* MCLimAcc = new TH1F("MCLimAcc_PtBins","B counts within |y|<0.5",100,0,40);
	TH1F* MCLimAcc_PtBins = (TH1F*)MCLimAcc->Rebin(NBinsPt,"",BinsPt);	
	
	/*
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(1,"MC stack called");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(2,"MC: B- -> D0(bar)");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(3,"MC: B- -> D0(bar)+X");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(4,"MC: B- -> D0 & Pi-");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(5,"MC: B+ -> D0bar & Pi+");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(6,"MC: B- -> D0+Pi-");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(7,"MC: B+ -> D0bar+Pi+");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(8,"MC: B- in |y|<0.9");
	 fMCTruth_Bplus->GetXaxis()->SetBinLabel(9,"MC: B+ in |y|<0.9");*/
	
	fOutputList->Add(fNentries);
	fOutputList->Add(fMCTruth_Bplus);
	fOutputList->Add(fMCTruth_Bminus);
	fOutputList->Add(fMCTruth2);
	fOutputList->Add(hPdgInvMassB);
	fOutputList->Add(hCalcInvMassB);
	fOutputList->Add(hPdgInvMassD);
	fOutputList->Add(hCalcInvMassD);
	fOutputList->Add(hEtaD0);
	fOutputList->Add(hEtaPi);
	fOutputList->Add(hEtaBoth);
	fOutputList->Add(hEtaBothiR);
	fOutputList->Add(hPtD0);
	fOutputList->Add(hPtPi);
	fOutputList->Add(hPtVsEtaD0);
	fOutputList->Add(hPtVsEtaPi);
	fOutputList->Add(hSgn_CosThetaS_0);
	fOutputList->Add(hSgn_CosThetaS_1);
	fOutputList->Add(hSgn_CosThetaS_2);
	fOutputList->Add(hSgn_CosThetaS_3);
	fOutputList->Add(hSgn_CosThetaS_4);
	fOutputList->Add(hSgn_CosThetaS_5);
	fOutputList->Add(hSgn_CosThetaS_6);
	fOutputList->Add(hSgn_CosThetaS_7);
	fOutputList->Add(hSgn_CosThetaS_8);
	fOutputList->Add(hSgn_CosThetaS_9);
	fOutputList->Add(hSgn_CosThetaS_10);
	fOutputList->Add(hSgn_CosThetaS_11);
	fOutputList->Add(hSgn_CosThetaS_12);
	fOutputList->Add(hSgn_CosThetaS_D_0);
	fOutputList->Add(hSgn_CosThetaS_D_1);
	fOutputList->Add(hSgn_CosThetaS_D_2);
	fOutputList->Add(hSgn_CosThetaS_D_3);
	fOutputList->Add(hSgn_CosThetaS_D_4);
	fOutputList->Add(hSgn_CosThetaS_D_5);
	fOutputList->Add(hSgn_CosThetaS_D_6);
	fOutputList->Add(hSgn_CosThetaS_D_7);
	fOutputList->Add(hSgn_CosThetaS_D_8);
	fOutputList->Add(hSgn_CosThetaS_D_9);
	fOutputList->Add(hSgn_CosThetaS_D_10);
	fOutputList->Add(hSgn_CosThetaS_D_11);
	fOutputList->Add(hSgn_CosThetaS_D_12);
	fOutputList->Add(hThetaD0);
	fOutputList->Add(hThetaPi);
	fOutputList->Add(hThetaD0Ka);
	fOutputList->Add(hThetaD0Pi);
	
	for(Int_t i=1; i<16;i++){
		fOutputList->Add(hSgn_CountsY_PtBins_0[i]);
	}
	fOutputList->Add(MCLimAcc_PtBins);
	
	
	TH1F* hSgn_CountsY_byHand = new TH1F(Form("hSgn_CountsY_byHand_PtBins"),Form("B(by Hand) counts within |y|<0.8"),100,0,40);
	TH1F* hSgn_CountsY_byHand_PtBins = (TH1F*)hSgn_CountsY_byHand->Rebin(NBinsPt,"",BinsPt);
	TH1F* hSgn_CountsY_bQuark = new TH1F(Form("hSgn_CountsY_bQuark_PtBins"),Form("B(b quark) counts within |y|<0.8"),100,0,40);
	TH1F* hSgn_CountsY_bQuark_PtBins = (TH1F*)hSgn_CountsY_bQuark->Rebin(NBinsPt,"",BinsPt);;
	
	fComparList->Add(hSgn_CountsY_byHand_PtBins);
	fComparList->Add(hSgn_CountsY_bQuark_PtBins);
	
	TH1F* hComp_byHand = new TH1F(Form("hComp_byHand"),Form("B(by Hand) rapidity in pt;y;counts"),400,-2.,2.);
	TH1F* hComp_bQuark = new TH1F(Form("hComp_bQuark"),Form("B(b quark) rapidity in pt;y;counts"),400,-2.,2.);
	
	fComparList->Add(hComp_byHand);
	fComparList->Add(hComp_bQuark);
	
	TH1F* hComp_byHand_PtBins[13];
	TH1F* hComp_bQuark_PtBins[13];
	TH2F* declBvsD[13];
	
	TH1F* hSgn_Ctp_PtBin[13];
	TH1F* hSgn_EtaPi_PtBin[13];
	TH1F* hSgn_EtaD0_PtBin[13];

	
	for(Int_t ipt=0; ipt<13; ipt++){
		hComp_byHand_PtBins[ipt] = new TH1F(Form("hComp_byHand_PtBins[%i]",ipt),Form("B(by Hand) rapidity in pt bin %i;y;counts",ipt),400,-2.,2.);
		hComp_bQuark_PtBins[ipt] = new TH1F(Form("hComp_bQuark_PtBins[%i]",ipt),Form("B(b quark) rapidity in pt bin %i;y;counts",ipt),400,-2.,2.);
		declBvsD[ipt] = new TH2F(Form("DecayLenghtBvsD_[%i]",ipt),"Decay Length B vs D;D Decay Length [cm]; B Decay Length [cm]",1000,0.,0.1,1000,0.,0.1);
		fComparList->Add(declBvsD[ipt]);
		hSgn_Ctp_PtBin[ipt] = new TH1F(Form("hSgn_Ctp_PtBin[%i]",ipt),Form("D cos(#theta_{p}) in pt bin %i;;counts",ipt),20000,-1.,1.);
		hSgn_EtaPi_PtBin[ipt] = new TH1F(Form("hSgn_EtaPi_PtBin[%i]",ipt),Form("Pi from B #eta in pt bin %i;;counts",ipt),400,-2.,2.);
		hSgn_EtaD0_PtBin[ipt] = new TH1F(Form("hSgn_EtaD0_PtBin[%i]",ipt),Form("D0 from B #eta in pt bin %i;;counts",ipt),400,-2.,2.);
		fOutputList->Add(hSgn_Ctp_PtBin[ipt]);		
	}
	
	for(Int_t i=0; i<13;i++){
		fComparList->Add(hComp_byHand_PtBins[i]);
		fComparList->Add(hComp_bQuark_PtBins[i]);
		fOutputList->Add(hSgn_EtaPi_PtBin[i]);		
		fOutputList->Add(hSgn_EtaD0_PtBin[i]);
	}
	
	// Post the data
	PostData(1,fOutputList);
	PostData(2,fComparList);
	
	return;
}

//________________________________________________________________________
void AliAnalysisTaskSEFindBhadMC::UserExec(Option_t */*option*/)
{
	// Execute analysis for current event:
	// heavy flavor candidates association to MC truth
	//cout<<"I'm in UserExec"<<endl;
	//AliLog::SetGlobalDebugLevel(1);
	
	// Analysis Flags
	Bool_t MCLimAcc			= fMClimAcc;
	Bool_t withFiducialCut	= fFiducialCut;
	Bool_t TestwD   = kFALSE;
	Bool_t oldMC	= kFALSE;
	
	Float_t thmin = (180./TMath::Pi())*2.*atan(exp(-1.))*TMath::Pi()/180.;
	Float_t thmax = (180./TMath::Pi())*2.*atan(exp( 1.))*TMath::Pi()/180.;
	
	Bool_t fromBquark = kTRUE;
	Int_t candidateCounter = 0;
	AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
	fMCTruth_Bplus->Fill(0.);
	Bool_t allDaughtersInRange = kTRUE;
	
	TClonesArray *arrayMC=static_cast<TClonesArray*>(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
	if(!arrayMC)return;
	
	// load MC header
	AliAODMCHeader *mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
	if(!mcHeader) {
		printf("AliAnalysisTaskSEBpmMass::UserExec: MC header branch not found!\n");
		return;
	}
	/*
	AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(mcHeader);

	for (Int_t i=0; i<(Int_t)mcHeader->GetNCocktailHeaders(); i++) {
		hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(mcHeader->GetCocktailHeader(i));
		if (hijingGenHeader){
			fhijingGenHeader=hijingGenHeader;
			break;
		}
	}
	if(!hijingGenHeader){
		cout << "Hijing event header not found" << endl;
		return;
	}
	*/
	AliAODVertex *aodVtx = (AliAODVertex*)aod->GetPrimaryVertex();
	if (!aodVtx) return;
	
	for(Int_t ic=0; ic<arrayMC->GetEntriesFast(); ic++){
		
		AliAODMCParticle *BmesonCandidate=static_cast<AliAODMCParticle*>(arrayMC->At(ic));
		if(!BmesonCandidate) continue;

		// For B Meson
		if(MCLimAcc==kFALSE && withFiducialCut==kTRUE){
			// Fiducial acceptance cut of B
			if(TMath::Abs(BmesonCandidate->Y())>=0.8){
				//cout << "FIDUCIAL CUT ACTIVE" << endl;
				BmesonCandidate=0;
				continue;
			}
		}
		else if(MCLimAcc==kFALSE && withFiducialCut==kFALSE){
			//cout << "FIDUCIAL CUT INACTIVE" << endl;
			// nothing to be done
		}
		else if(MCLimAcc==kTRUE){
			// MC limited acceptance cut
			if(TMath::Abs(BmesonCandidate->Y())>=0.5){
				BmesonCandidate=0;
				continue;
			}
		}
		//---
		// Reject if comming from B quark
		if(!TestwD){
			if(IsTrackInjected(BmesonCandidate->GetLabel(),mcHeader,arrayMC)==kTRUE)fromBquark=kFALSE;
			else continue;
		}
		//---------------------------------------------------------------------------------------------------------------------------
		// Selection of B meson, if in channel
		//
		Int_t pdgBcand=BmesonCandidate->GetPdgCode();

		// Check if candidate is a B Meson
		if(pdgBcand==521){
			if(!fromBquark){
				fNentries->Fill(2);	
			}
			if(fromBquark){	
				fNentries->Fill(12);
			}
			fMCTruth_Bplus->Fill(1);
		}
		else if(pdgBcand==-521){
			if(!fromBquark){
				fNentries->Fill(3);	
			}
			if(fromBquark){	
				fNentries->Fill(13);
			}
			fMCTruth_Bminus->Fill(1);
		}
		else{
			BmesonCandidate=0;
			//cout << "Rejected: Label is: " << ic << endl;
			//cout << "This is not a B meson" << endl;
			continue;
		}
		
		// Check if B Meson has 2 daughters & load both daughters
		if(BmesonCandidate->GetNDaughters()!=2){
			if(!fromBquark){
				fNentries->Fill(4);	
			}
			if(fromBquark){	
				fNentries->Fill(14);
			}
			fMCTruth_Bplus->Fill(2);
			BmesonCandidate=0;
			// cout << "Rejected: Label is: " << ic << endl;
			// cout << "B does not have 2 daughters" << endl;
			continue;
		}
		if(pdgBcand==521) fMCTruth_Bplus->Fill(3);
		if(pdgBcand==-521)fMCTruth_Bminus->Fill(3);
		fMCTruth2->Fill(pdgBcand);
		Int_t d1 = BmesonCandidate->GetDaughter(0);
		Int_t d2 = BmesonCandidate->GetDaughter(1);
		AliAODMCParticle *Bdaughter1=(AliAODMCParticle*)arrayMC->At(d1);
		AliAODMCParticle *Bdaughter2=(AliAODMCParticle*)arrayMC->At(d2);
		Int_t pdg1=Bdaughter1->GetPdgCode();
		Int_t pdg2=Bdaughter2->GetPdgCode();
		fMCTruth2->Fill(pdg1);
		fMCTruth2->Fill(pdg2);
		//---

		// Check which daughter is D Meson and if D Meson has two daughters
		if(TMath::Abs(pdg1)==421){
			if(TestwD) {
				if(IsTrackInjected(Bdaughter1->GetLabel(),mcHeader,arrayMC)==kFALSE)fromBquark=kTRUE;
			}
			if(pdgBcand==521) fMCTruth_Bplus->Fill(4);
			if(pdgBcand==-521)fMCTruth_Bminus->Fill(4);
			if(Bdaughter1->GetNDaughters()!=2){
				if(pdgBcand==521) fMCTruth_Bplus->Fill(10);
				if(pdgBcand==-521)fMCTruth_Bminus->Fill(10);
				Bdaughter1=0;
				Bdaughter2=0;
				BmesonCandidate=0;
				// cout << "Rejected: Label is: " << ic << endl;
				// cout << "D meson does not have 2 daughters" << endl;
				if(!fromBquark){
					fNentries->Fill(8);	
				}
				if(fromBquark){	
					fNentries->Fill(18);
				}
				continue;
			}			
		}
		else if(TMath::Abs(pdg2)==421){
			if(TestwD) {
				if(IsTrackInjected(Bdaughter2->GetLabel(),mcHeader,arrayMC)==kFALSE)fromBquark=kTRUE;
			}
			if(pdgBcand==521) fMCTruth_Bplus->Fill(5);
			if(pdgBcand==-521)fMCTruth_Bminus->Fill(5);	
			if(Bdaughter1->GetNDaughters()!=2){
				if(pdgBcand==521) fMCTruth_Bplus->Fill(10);
				if(pdgBcand==-521)fMCTruth_Bminus->Fill(10);
				Bdaughter1=0;
				Bdaughter2=0;
				BmesonCandidate=0;
				// cout << "Rejected: Label is: " << ic << endl;
				// cout << "D meson does not have 2 daughters" << endl;
				if(!fromBquark){
					fNentries->Fill(8);	
				}
				if(fromBquark){	
					fNentries->Fill(18);
				}
				continue;
			}
		}
		else{
			if(!fromBquark){
				fNentries->Fill(5);	
			}
			if(fromBquark){	
				fNentries->Fill(15);
			}
			if(pdgBcand==521) fMCTruth_Bplus->Fill(6);
			if(pdgBcand==-521)fMCTruth_Bminus->Fill(6);
			Bdaughter1=0;
			Bdaughter2=0;
			BmesonCandidate=0;
			// cout << "Rejected: Label is: " << ic << endl;
			// cout << "No D meson at all" << endl;
			continue;
		}
		//---
		// Theta & Y acceptance cut
		if((Bdaughter1->Theta()<thmin || Bdaughter1->Theta()>thmax || TMath::Abs(Bdaughter1->Y())>1.) &&
		   (TMath::Abs(pdg1)!=421)){
			Bdaughter1=0;
			Bdaughter2=0;
			BmesonCandidate=0;
			continue;
		}
		if((Bdaughter2->Theta()<thmin || Bdaughter2->Theta()>thmax || TMath::Abs(Bdaughter2->Y())>1.) &&
		   (TMath::Abs(pdg2)!=421)){
			Bdaughter1=0;
			Bdaughter2=0;
			BmesonCandidate=0;
			continue;
		}
		//---

		// Check if Pion passes pt and eta cut, and if there is a Pion at all?
		if(TMath::Abs(pdg2)==211){
			if(MCLimAcc==kFALSE){
				if(TMath::Abs(Bdaughter2->Eta())>0.9 || Bdaughter2->Pt()<0.1) {
					Bdaughter1=0;
					Bdaughter2=0;
					BmesonCandidate=0;
					// cout << "Rejected: Label is: " << ic << endl;
					// cout << "Pion does not pass pt and eta cut" << endl;
					if(!fromBquark){
						fNentries->Fill(6);	
					}
					if(fromBquark){	
						fNentries->Fill(16);
					}
					continue;
				}
			}
			if(pdgBcand==521) fMCTruth_Bplus->Fill(7);
			if(pdgBcand==-521)fMCTruth_Bminus->Fill(7);
		}
		else if(TMath::Abs(pdg1)==211){
			if(MCLimAcc==kFALSE){
				if(TMath::Abs(Bdaughter1->Eta())>0.9 || Bdaughter1->Pt()<0.1){
					Bdaughter1=0;
					Bdaughter2=0;
					BmesonCandidate=0;
					// cout << "Rejected: Label is: " << ic << endl;
					// cout << "Pion does not pass pt and eta cut" << endl;
					if(!fromBquark){
						fNentries->Fill(6);	
					}
					if(fromBquark){	
						fNentries->Fill(16);
					}
					continue;
				}
			}
			if(pdgBcand==521) fMCTruth_Bplus->Fill(8);
			if(pdgBcand==-521)fMCTruth_Bminus->Fill(8);
		}
		else{
			if(!fromBquark){
				fNentries->Fill(7);	
			}
			if(fromBquark){	
				fNentries->Fill(17);
			}
			if(pdgBcand==521) fMCTruth_Bplus->Fill(9);
			if(pdgBcand==-521)fMCTruth_Bminus->Fill(9);
			Bdaughter1=0;
			Bdaughter2=0;
			BmesonCandidate=0;
			// cout << "Rejected: Label is: " << ic << endl;
			// cout << "No Pion in B decay" << endl;
			continue;
		}
		//---
		if(pdgBcand==521) fMCTruth_Bplus->Fill(11);
		if(pdgBcand==-521)fMCTruth_Bminus->Fill(11);
		
		// Examine D daughters
		Int_t d3(-99),d4(-99);
		if(TMath::Abs(pdg1)==421){
			d3 = Bdaughter1->GetDaughter(0);
			d4 = Bdaughter1->GetDaughter(1);
		}
		if(TMath::Abs(pdg2)==421){
			d3 = Bdaughter2->GetDaughter(0);
			d4 = Bdaughter2->GetDaughter(1);
		}
		// Check if D daughters exist
		AliAODMCParticle *part3=(AliAODMCParticle*)arrayMC->At(d3);
		if(!part3){
			Bdaughter1=0;
			Bdaughter2=0;
			BmesonCandidate=0;
			// cout << "Rejected: Label is: " << ic << endl;
			// cout << "1) D daughter does not exist" << endl;
			if(!fromBquark){
				fNentries->Fill(9);	
			}
			if(fromBquark){	
				fNentries->Fill(19);
			}
			continue;
		}
		AliAODMCParticle *part4=(AliAODMCParticle*)arrayMC->At(d4);
		if(!part4){
			Bdaughter1=0;
			Bdaughter2=0;
			part3=0;
			BmesonCandidate=0;
			// cout << "Rejected: Label is: " << ic << endl;
			// cout << "2) D daughter does not exist" << endl;
			if(!fromBquark){
				fNentries->Fill(9);	
			}
			if(fromBquark){	
				fNentries->Fill(19);
			}
			continue;
		}
		Int_t pdg3=part3->GetPdgCode();
		Int_t pdg4=part4->GetPdgCode();
		fMCTruth2->Fill(pdg3);
		fMCTruth2->Fill(pdg4);
		
		// Check if D daughters pass pt and eta cut
		// Theta acceptance cut
		if(part3->Theta()<thmin || part3->Theta()>thmax){
			Bdaughter1=0;
			Bdaughter2=0;
			part3=0;
			part4=0;
			BmesonCandidate=0;
			allDaughtersInRange=kFALSE;
			continue;
		}
		if(part4->Theta()<thmin || part4->Theta()>thmax){
			Bdaughter1=0;
			Bdaughter2=0;
			part3=0;
			part4=0;
			BmesonCandidate=0;
			allDaughtersInRange=kFALSE;
			continue;
		}
		//---
		// Y acceptance cut
		if(TMath::Abs(part3->Y())>1.){
			Bdaughter1=0;
			Bdaughter2=0;
			part3=0;
			part4=0;
			BmesonCandidate=0;
			allDaughtersInRange=kFALSE;
			continue;
		}
		if(TMath::Abs(part4->Y())>1.){
			Bdaughter1=0;
			Bdaughter2=0;
			part3=0;
			part4=0;
			BmesonCandidate=0;
			allDaughtersInRange=kFALSE;
			continue;
		}

		if(MCLimAcc==kFALSE){
			if(TMath::Abs(part3->Eta())>0.9 || part3->Pt()<0.1 || TMath::Abs(part4->Eta())>0.9 || part4->Pt()<0.1){
				Bdaughter1=0;
				Bdaughter2=0;
				part3=0;
				part4=0;
				BmesonCandidate=0;
				// cout << "Rejected: Label is: " << ic << endl;
				// cout << "D daughters did not pass pt and eta cut " << endl;
				if(!fromBquark){
					fNentries->Fill(10);	
				}
				if(fromBquark){	
					fNentries->Fill(20);
				}
				continue;
			}
		}
		if(!fromBquark){
			fNentries->Fill(11);	
		}
		if(fromBquark){	
			fNentries->Fill(21);
		}
		
		if((pdgBcand== 521 && pdg1==-421 && pdg2== 211 && pdg3== 321 && pdg4==-211) ||
		   (pdgBcand== 521 && pdg1==-421 && pdg2== 211 && pdg3==-211 && pdg4== 321) ||
		   (pdgBcand== 521 && pdg1== 211 && pdg2==-421 && pdg3== 321 && pdg4==-211) ||
		   (pdgBcand== 521 && pdg1== 211 && pdg2==-421 && pdg3==-211 && pdg4== 321)){
			
			if(!fromBquark){
				fNentries->Fill(0);	
			}
			if(fromBquark){	
				fNentries->Fill(1);
			}
			
			fMCTruth_Bplus->Fill(12);
			hEtaBoth->Fill(BmesonCandidate->Y());
			
			fMCTruth_Bplus->Fill(14);
			
			Double_t BdecayLength = 0;
			Double_t DdecayLength = 0;
			if(pdg1==-421){
				hEtaD0->Fill(Bdaughter1->Eta());
				hEtaPi->Fill(Bdaughter2->Eta());
				hPtD0->Fill(Bdaughter1->Pt());
				hPtPi->Fill(Bdaughter2->Pt());
				hPtVsEtaD0->Fill(Bdaughter1->Pt(),Bdaughter1->Eta());
				hPtVsEtaPi->Fill(Bdaughter2->Pt(),Bdaughter2->Eta());
				hEtaBothiR->Fill(Bdaughter1->Eta(),Bdaughter2->Eta());
				BdecayLength = CalculateDecayLength(Bdaughter1,aodVtx);
				hThetaD0->Fill(Bdaughter1->Theta());
				hThetaPi->Fill(Bdaughter2->Theta());
				if(TMath::Abs(pdg3)==321){
				hThetaD0Ka->Fill(part3->Theta());
				hThetaD0Pi->Fill(part4->Theta());
				}
				else{
					hThetaD0Ka->Fill(part4->Theta());
					hThetaD0Pi->Fill(part3->Theta());
				}
			}
			else{
				hEtaPi->Fill(Bdaughter1->Eta());
				hEtaD0->Fill(Bdaughter2->Eta());	
				hPtPi->Fill(Bdaughter1->Pt());
				hPtD0->Fill(Bdaughter2->Pt());
				hPtVsEtaPi->Fill(Bdaughter1->Pt(),Bdaughter1->Eta());
				hPtVsEtaD0->Fill(Bdaughter2->Pt(),Bdaughter2->Eta());
				hEtaBothiR->Fill(Bdaughter2->Eta(),Bdaughter1->Eta());
				BdecayLength = CalculateDecayLength(Bdaughter2,aodVtx);
				
				hThetaD0->Fill(Bdaughter2->Theta());
				hThetaPi->Fill(Bdaughter1->Theta());
				if(TMath::Abs(pdg3)==321){
					hThetaD0Ka->Fill(part3->Theta());
					hThetaD0Pi->Fill(part4->Theta());
				}
				else{
					hThetaD0Ka->Fill(part4->Theta());
					hThetaD0Pi->Fill(part3->Theta());
				}
			}
			if(TMath::Abs(pdg3)==321) DdecayLength = CalculateDecayLength(part3,aodVtx);
			if(TMath::Abs(pdg4)==321) DdecayLength = CalculateDecayLength(part4,aodVtx);
			
			TString fillthis="";
			fillthis="hPdgInvMassB";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(BmesonCandidate->M());
			fillthis="hCalcInvMassB";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(BmesonCandidate->GetCalcMass());

			Double_t bMesonPt = BmesonCandidate->Pt();
			if(fromBquark==kTRUE){
				fillthis="hSgn_CountsY_bQuark_PtBins";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			else{
				fillthis="hSgn_CountsY_byHand_PtBins";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(oldMC){
				if(!fromBquark) continue;
			}
			if(TMath::Abs(pdg1)==421) {
				fillthis="hPdgInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->M());
				fillthis="hCalcInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->GetCalcMass());
			}
			if(TMath::Abs(pdg2)==421) {
				fillthis="hPdgInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->M());
				fillthis="hCalcInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->GetCalcMass());
			}

			// for tests with CF framework
			if(TestwD==kTRUE){
				if(TMath::Abs(pdg1)==421) {
					bMesonPt = Bdaughter1->Pt();
				}
				if(TMath::Abs(pdg2)==421) {
					bMesonPt = Bdaughter2->Pt();
				}
			}
			// Checking number of B versus Rapidity
			if(MCLimAcc==kTRUE){
				fillthis="MCLimAcc_PtBins";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			
			Double_t dMesonAbsRapidity = TMath::Abs(BmesonCandidate->Y());
			
			if(TMath::Abs(dMesonAbsRapidity)<=0.1 && TMath::Abs(dMesonAbsRapidity)>0.0){
				fillthis="hSgn_CountsY_PtBins_0[1]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.2 && TMath::Abs(dMesonAbsRapidity)>0.1){
				fillthis="hSgn_CountsY_PtBins_0[2]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);	
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.3 && TMath::Abs(dMesonAbsRapidity)>0.2){
				fillthis="hSgn_CountsY_PtBins_0[3]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.4 && TMath::Abs(dMesonAbsRapidity)>0.3){
				fillthis="hSgn_CountsY_PtBins_0[4]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.5 && TMath::Abs(dMesonAbsRapidity)>0.4){
				fillthis="hSgn_CountsY_PtBins_0[5]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.6 && TMath::Abs(dMesonAbsRapidity)>0.5){
				fillthis="hSgn_CountsY_PtBins_0[6]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.7 && TMath::Abs(dMesonAbsRapidity)>0.6){
				fillthis="hSgn_CountsY_PtBins_0[7]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.8 && TMath::Abs(dMesonAbsRapidity)>0.7){
				fillthis="hSgn_CountsY_PtBins_0[8]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.9 && TMath::Abs(dMesonAbsRapidity)>0.8){
				fillthis="hSgn_CountsY_PtBins_0[9]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.0 && TMath::Abs(dMesonAbsRapidity)>0.9){
				fillthis="hSgn_CountsY_PtBins_0[10]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.1 && TMath::Abs(dMesonAbsRapidity)>1.0){
				fillthis="hSgn_CountsY_PtBins_0[11]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.2 && TMath::Abs(dMesonAbsRapidity)>1.1){
				fillthis="hSgn_CountsY_PtBins_0[12]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.3 && TMath::Abs(dMesonAbsRapidity)>1.2){
				fillthis="hSgn_CountsY_PtBins_0[13]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.4 && TMath::Abs(dMesonAbsRapidity)>1.3){
				fillthis="hSgn_CountsY_PtBins_0[14]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)>1.4){
				fillthis="hSgn_CountsY_PtBins_0[15]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			
			// Calculating Cos(ThetaStar) B
			Double_t massvtx = 5.279;
			Double_t massp[2];
			massp[0] = 1.864;
			massp[1] = 0.139;
			Double_t pStar = TMath::Sqrt(TMath::Power(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1],2.)-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);
			Double_t beta = BmesonCandidate->P()/BmesonCandidate->E();
			Double_t gamma = BmesonCandidate->E()/massvtx;
			Double_t cts =0;
			
			TVector3 momTot(BmesonCandidate->Px(),BmesonCandidate->Py(),BmesonCandidate->Pz());

			if(pdg1==211){
				TVector3 mom(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
				Double_t Q1 = mom.Dot(momTot)/momTot.Mag();
				
				cts = (Q1/gamma-beta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;	 				
			}
			if(pdg2==211){
				TVector3 mom(Bdaughter2->Px(),Bdaughter2->Py(),Bdaughter2->Pz());
				Double_t Q1 = mom.Dot(momTot)/momTot.Mag();
				
				cts = (Q1/gamma-beta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;
			}

			// Calculating Cos(ThetaStar) D
			Double_t massvtxD = 1.864;
			Double_t masspD[2];
			masspD[0] = 0.493;
			masspD[1] = 0.139;
			Double_t pStarD = TMath::Sqrt(TMath::Power(massvtxD*massvtxD-masspD[0]*masspD[0]-masspD[1]*masspD[1],2.)-4.*masspD[0]*masspD[0]*masspD[1]*masspD[1])/(2.*massvtxD);
			Double_t betaD = Bdaughter1->P()/Bdaughter1->E();
			Double_t gammaD = Bdaughter1->E()/massvtxD;
			Double_t ctsD =0;
			
			TVector3 momTotD(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
			if(pdg2==-421){
				betaD = Bdaughter2->P()/Bdaughter2->E();
				gammaD = Bdaughter2->E()/massvtxD;
				TVector3 momTotD(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
			}
			if(pdg3==321){
				TVector3 momD(part3->Px(),part3->Py(),part3->Pz());
				Double_t Q1D = momD.Dot(momTotD)/momTotD.Mag();
				
				ctsD = (Q1D/gammaD-betaD*TMath::Sqrt(pStarD*pStarD+masspD[0]*masspD[0]))/pStarD;
			}
			if(pdg4==321){
				TVector3 momD(part4->Px(),part4->Py(),part4->Pz());
				Double_t Q1D = momD.Dot(momTotD)/momTotD.Mag();
				
				ctsD = (Q1D/gammaD-betaD*TMath::Sqrt(pStarD*pStarD+masspD[0]*masspD[0]))/pStarD;
			}
			// ---
			// Cosine of pointing angle  of D in space assuming it is produced at primary vtx
			Double_t cos =0.;
			if(pdg1==-421){
				TVector3 momD(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
				Double_t daughterVtx[3]={0,0,0};
				part3->XvYvZv(daughterVtx); // cm
				
				TVector3 fline(daughterVtx[0]-aodVtx->GetX(),
							   daughterVtx[1]-aodVtx->GetY(),
							   daughterVtx[2]-aodVtx->GetZ());
				
				Double_t ptot2 = momD.Mag2()*fline.Mag2();
				if(ptot2 <= 0) {
					break;
				} 
				else {
					cos = momD.Dot(fline)/TMath::Sqrt(ptot2);
					if(cos >  1.0) cos =  1.0;
					if(cos < -1.0) cos = -1.0;
				}
			}
			else if(pdg2==-421){
				TVector3 momD(Bdaughter2->Px(),Bdaughter2->Py(),Bdaughter2->Pz());
				Double_t daughterVtx[3]={0,0,0};
				part3->XvYvZv(daughterVtx); // cm
				
				TVector3 fline(daughterVtx[0]-aodVtx->GetX(),
							   daughterVtx[1]-aodVtx->GetY(),
							   daughterVtx[2]-aodVtx->GetZ());
				
				Double_t ptot2 = momD.Mag2()*fline.Mag2();
				if(ptot2 <= 0) {
					break;
				} 
				else {
					cos = momD.Dot(fline)/TMath::Sqrt(ptot2);
					if(cos >  1.0) cos =  1.0;
					if(cos < -1.0) cos = -1.0;
				}
			}
			// ---
			// Cut on CTP as in AOD filtering
			if(cos<0.5) continue;
			
			
			if(fromBquark==kTRUE){
				fillthis="hComp_bQuark";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
			}
			else{
				fillthis="hComp_byHand";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
			}
			
			if(bMesonPt<0.5){
				fMCTruth_Bplus->Fill(16);
				hSgn_CosThetaS_0->Fill(cts); hSgn_CosThetaS_D_0->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[0]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[0]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[0]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[0]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=0.5 && bMesonPt<1.0){
				fMCTruth_Bplus->Fill(17);
				hSgn_CosThetaS_1->Fill(cts); hSgn_CosThetaS_D_1->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[1]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[1]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[1]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[1]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=1.0 && bMesonPt<2.0){
				fMCTruth_Bplus->Fill(18);
				hSgn_CosThetaS_2->Fill(cts); hSgn_CosThetaS_D_2->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[2]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[2]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[2]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[2]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=2.0 && bMesonPt<3.0){
				fMCTruth_Bplus->Fill(19);
				hSgn_CosThetaS_3->Fill(cts); hSgn_CosThetaS_D_3->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[3]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[3]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[3]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[3]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=3.0 && bMesonPt<4.0){
				fMCTruth_Bplus->Fill(20);
				hSgn_CosThetaS_4->Fill(cts); hSgn_CosThetaS_D_4->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[4]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[4]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[4]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[4]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=4.0 && bMesonPt<5.0){
				fMCTruth_Bplus->Fill(21);
				hSgn_CosThetaS_5->Fill(cts); hSgn_CosThetaS_D_5->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[5]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[5]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[5]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[5]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=5.0 && bMesonPt<6.0){
				fMCTruth_Bplus->Fill(22);
				hSgn_CosThetaS_6->Fill(cts); hSgn_CosThetaS_D_6->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[6]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[6]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[6]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[6]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=6.0 && bMesonPt<8.0){
				fMCTruth_Bplus->Fill(23);
				hSgn_CosThetaS_7->Fill(cts); hSgn_CosThetaS_D_7->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[7]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[7]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[7]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[7]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=8.0 && bMesonPt<12.){
				fMCTruth_Bplus->Fill(24);
				hSgn_CosThetaS_8->Fill(cts); hSgn_CosThetaS_D_8->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[8]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[8]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[8]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[8]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=12. && bMesonPt<16.){
				fMCTruth_Bplus->Fill(25);
				hSgn_CosThetaS_9->Fill(cts); hSgn_CosThetaS_D_9->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[9]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[9]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[9]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[9]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=16. && bMesonPt<20.){
				fMCTruth_Bplus->Fill(26);
				hSgn_CosThetaS_10->Fill(cts); hSgn_CosThetaS_D_10->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[10]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[10]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[10]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[10]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=20. && bMesonPt<24.){
				fMCTruth_Bplus->Fill(27);
				hSgn_CosThetaS_11->Fill(cts); hSgn_CosThetaS_D_11->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[11]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[11]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[11]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[11]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=24.){
				fMCTruth_Bplus->Fill(28);	
				hSgn_CosThetaS_12->Fill(cts); hSgn_CosThetaS_D_12->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[12]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[12]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[12]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[12]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
		}
		else if((pdgBcand==-521 && pdg1== 421 && pdg2==-211 && pdg3==-321 && pdg4== 211) ||
				(pdgBcand==-521 && pdg1== 421 && pdg2==-211 && pdg3== 211 && pdg4==-321) ||
				(pdgBcand==-521 && pdg1==-211 && pdg2== 421 && pdg3==-321 && pdg4== 211) ||
				(pdgBcand==-521 && pdg1==-211 && pdg2== 421 && pdg3== 211 && pdg4==-321)){
			
			if(!fromBquark){
				fNentries->Fill(0);	
			}
			if(fromBquark){	
				fNentries->Fill(1);
			}
			
			fMCTruth_Bminus->Fill(13);
			hEtaBoth->Fill(BmesonCandidate->Y());
			
			fMCTruth_Bminus->Fill(15);
			
			Double_t BdecayLength = 0;
			Double_t DdecayLength = 0;
			if(pdg1==421){
				hEtaD0->Fill(Bdaughter1->Eta());
				hEtaPi->Fill(Bdaughter2->Eta());
				hPtD0->Fill(Bdaughter1->Pt());
				hPtPi->Fill(Bdaughter2->Pt());
				hPtVsEtaD0->Fill(Bdaughter1->Pt(),Bdaughter1->Eta());
				hPtVsEtaPi->Fill(Bdaughter2->Pt(),Bdaughter2->Eta());
				hEtaBothiR->Fill(Bdaughter1->Eta(),Bdaughter2->Eta());
				BdecayLength = CalculateDecayLength(Bdaughter1,aodVtx);
				hThetaD0->Fill(Bdaughter1->Theta());
				hThetaPi->Fill(Bdaughter2->Theta());
				if(TMath::Abs(pdg3)==321){
					hThetaD0Ka->Fill(part3->Theta());
					hThetaD0Pi->Fill(part4->Theta());
				}
				else{
					hThetaD0Ka->Fill(part4->Theta());
					hThetaD0Pi->Fill(part3->Theta());
				}

			}
			else{
				hEtaPi->Fill(Bdaughter1->Eta());
				hEtaD0->Fill(Bdaughter2->Eta());	
				hPtPi->Fill(Bdaughter1->Pt());
				hPtD0->Fill(Bdaughter2->Pt());
				hPtVsEtaPi->Fill(Bdaughter1->Pt(),Bdaughter1->Eta());
				hPtVsEtaD0->Fill(Bdaughter2->Pt(),Bdaughter2->Eta());
				hEtaBothiR->Fill(Bdaughter2->Eta(),Bdaughter1->Eta());
				BdecayLength = CalculateDecayLength(Bdaughter2,aodVtx);
				hThetaD0->Fill(Bdaughter2->Theta());
				hThetaPi->Fill(Bdaughter1->Theta());
				if(TMath::Abs(pdg3)==321){
					hThetaD0Ka->Fill(part3->Theta());
					hThetaD0Pi->Fill(part4->Theta());
				}
				else{
					hThetaD0Ka->Fill(part4->Theta());
					hThetaD0Pi->Fill(part3->Theta());
				}

			}
			if(TMath::Abs(pdg3)==321) DdecayLength = CalculateDecayLength(part3,aodVtx);
			if(TMath::Abs(pdg4)==321) DdecayLength = CalculateDecayLength(part4,aodVtx);
			
			TString fillthis="";
			fillthis="hPdgInvMassB";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(BmesonCandidate->M());
			fillthis="hCalcInvMassB";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(BmesonCandidate->GetCalcMass());

			Double_t bMesonPt = BmesonCandidate->Pt();
			if(fromBquark==kTRUE){
				fillthis="hSgn_CountsY_bQuark_PtBins";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			else{
				fillthis="hSgn_CountsY_byHand_PtBins";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(oldMC){
				if(!fromBquark) continue;
			}
			
			if(TMath::Abs(pdg1)==421) {
				fillthis="hPdgInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->M());
				fillthis="hCalcInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->GetCalcMass());
			}
			if(TMath::Abs(pdg2)==421) {
				fillthis="hPdgInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->M());
				fillthis="hCalcInvMassD";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->GetCalcMass());
			}
			// for tests with CF framework
			if(TestwD==kTRUE){
				if(TMath::Abs(pdg1)==421) {
					bMesonPt = Bdaughter1->Pt();
				}
				if(TMath::Abs(pdg2)==421) {
					bMesonPt = Bdaughter2->Pt();
				}
			}
			// Checking number of B versus Rapidity
			if(MCLimAcc){
				fillthis="MCLimAcc_PtBins";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}

			Double_t dMesonAbsRapidity = TMath::Abs(BmesonCandidate->Y());

			if(TMath::Abs(dMesonAbsRapidity)<=0.1 && TMath::Abs(dMesonAbsRapidity)>0.0){
				fillthis="hSgn_CountsY_PtBins_0[1]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.2 && TMath::Abs(dMesonAbsRapidity)>0.1){
				fillthis="hSgn_CountsY_PtBins_0[2]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);	
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.3 && TMath::Abs(dMesonAbsRapidity)>0.2){
				fillthis="hSgn_CountsY_PtBins_0[3]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.4 && TMath::Abs(dMesonAbsRapidity)>0.3){
				fillthis="hSgn_CountsY_PtBins_0[4]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.5 && TMath::Abs(dMesonAbsRapidity)>0.4){
				fillthis="hSgn_CountsY_PtBins_0[5]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.6 && TMath::Abs(dMesonAbsRapidity)>0.5){
				fillthis="hSgn_CountsY_PtBins_0[6]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.7 && TMath::Abs(dMesonAbsRapidity)>0.6){
				fillthis="hSgn_CountsY_PtBins_0[7]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.8 && TMath::Abs(dMesonAbsRapidity)>0.7){
				fillthis="hSgn_CountsY_PtBins_0[8]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=0.9 && TMath::Abs(dMesonAbsRapidity)>0.8){
				fillthis="hSgn_CountsY_PtBins_0[9]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.0 && TMath::Abs(dMesonAbsRapidity)>0.9){
				fillthis="hSgn_CountsY_PtBins_0[10]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.1 && TMath::Abs(dMesonAbsRapidity)>1.0){
				fillthis="hSgn_CountsY_PtBins_0[11]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.2 && TMath::Abs(dMesonAbsRapidity)>1.1){
				fillthis="hSgn_CountsY_PtBins_0[12]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.3 && TMath::Abs(dMesonAbsRapidity)>1.2){
				fillthis="hSgn_CountsY_PtBins_0[13]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)<=1.4 && TMath::Abs(dMesonAbsRapidity)>1.3){
				fillthis="hSgn_CountsY_PtBins_0[14]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			if(TMath::Abs(dMesonAbsRapidity)>1.4){
				fillthis="hSgn_CountsY_PtBins_0[15]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(bMesonPt);
			}
			
			// Calculating Cos(ThetaStar) B
			Double_t massvtx = 5.279;
			Double_t massp[2];
			massp[0] = 1.864;
			massp[1] = 0.139;
			Double_t pStar = TMath::Sqrt(TMath::Power(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1],2.)-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);
			Double_t beta = BmesonCandidate->P()/BmesonCandidate->E();
			Double_t gamma = BmesonCandidate->E()/massvtx;
			Double_t cts =0;
			
			TVector3 momTot(BmesonCandidate->Px(),BmesonCandidate->Py(),BmesonCandidate->Pz());
			
			if(pdg1==-211){
				TVector3 mom(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
				Double_t Q1 = mom.Dot(momTot)/momTot.Mag();
				
				cts = (Q1/gamma-beta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;	 				
			}
			if(pdg2==-211){
				TVector3 mom(Bdaughter2->Px(),Bdaughter2->Py(),Bdaughter2->Pz());
				Double_t Q1 = mom.Dot(momTot)/momTot.Mag();
				
				cts = (Q1/gamma-beta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;
			}
			// Calculating Cos(ThetaStar) D
			Double_t massvtxD = 1.864;
			Double_t masspD[2];
			masspD[0] = 0.493;
			masspD[1] = 0.139;
			Double_t pStarD = TMath::Sqrt(TMath::Power(massvtxD*massvtxD-masspD[0]*masspD[0]-masspD[1]*masspD[1],2.)-4.*masspD[0]*masspD[0]*masspD[1]*masspD[1])/(2.*massvtxD);
			Double_t betaD = Bdaughter1->P()/Bdaughter1->E();
			Double_t gammaD = Bdaughter1->E()/massvtxD;
			Double_t ctsD =0;
			
				TVector3 momTotD(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
			if(pdg2==421){
				TVector3 momTotD(Bdaughter2->Px(),Bdaughter2->Py(),Bdaughter2->Pz());
				betaD = Bdaughter2->P()/Bdaughter2->E();
				gammaD = Bdaughter2->E()/massvtxD;
			}
			if(pdg3==-321){
				TVector3 momD(part3->Px(),part3->Py(),part3->Pz());
				Double_t Q1D = momD.Dot(momTotD)/momTotD.Mag();
				
				ctsD = (Q1D/gammaD-betaD*TMath::Sqrt(pStarD*pStarD+masspD[0]*masspD[0]))/pStarD;
			}
			if(pdg4==-321){
				TVector3 momD(part4->Px(),part4->Py(),part4->Pz());
				Double_t Q1D = momD.Dot(momTotD)/momTotD.Mag();
				
				ctsD = (Q1D/gammaD-betaD*TMath::Sqrt(pStarD*pStarD+masspD[0]*masspD[0]))/pStarD;
			}

			// ---
			// Cosine of pointing angle  of D in space assuming it is produced at primary vtx
			Double_t cos =0.;
			if(pdg1==421){
				TVector3 momD(Bdaughter1->Px(),Bdaughter1->Py(),Bdaughter1->Pz());
				Double_t daughterVtx[3]={0,0,0};
				part3->XvYvZv(daughterVtx); // cm
				
				TVector3 fline(daughterVtx[0]-aodVtx->GetX(),
							   daughterVtx[1]-aodVtx->GetY(),
							   daughterVtx[2]-aodVtx->GetZ());
				
				Double_t ptot2 = momD.Mag2()*fline.Mag2();
				if(ptot2 <= 0) {
					break;
				} 
				else {
					cos = momD.Dot(fline)/TMath::Sqrt(ptot2);
					if(cos >  1.0) cos =  1.0;
					if(cos < -1.0) cos = -1.0;
				}
			}
			else if(pdg2==421){
				TVector3 momD(Bdaughter2->Px(),Bdaughter2->Py(),Bdaughter2->Pz());
				Double_t daughterVtx[3]={0,0,0};
				part3->XvYvZv(daughterVtx); // cm
				
				TVector3 fline(daughterVtx[0]-aodVtx->GetX(),
							   daughterVtx[1]-aodVtx->GetY(),
							   daughterVtx[2]-aodVtx->GetZ());
				
				Double_t ptot2 = momD.Mag2()*fline.Mag2();
				if(ptot2 <= 0) {
					break;
				} 
				else {
					cos = momD.Dot(fline)/TMath::Sqrt(ptot2);
					if(cos >  1.0) cos =  1.0;
					if(cos < -1.0) cos = -1.0;
				}
			}
			// ---
			// Cut on CTP as in AOD filtering
			if(cos<0.5) continue;
			
			if(fromBquark==kTRUE){
				fillthis="hComp_bQuark";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
			}
			else{
				fillthis="hComp_byHand";
				((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
			}
			if(bMesonPt<0.5){
				fMCTruth_Bminus->Fill(16);
				hSgn_CosThetaS_0->Fill(cts); hSgn_CosThetaS_D_0->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[0]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[0]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[0]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[0]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[0]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=0.5 && bMesonPt<1.0){
				fMCTruth_Bminus->Fill(17);
				hSgn_CosThetaS_1->Fill(cts); hSgn_CosThetaS_D_1->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[1]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[1]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[1]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());			
				}
				else{
					fillthis="hComp_byHand_PtBins[1]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());				
					fillthis="DecayLenghtBvsD_[1]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}	
			if(bMesonPt>=1.0 && bMesonPt<2.0){
				fMCTruth_Bminus->Fill(18);
				hSgn_CosThetaS_2->Fill(cts); hSgn_CosThetaS_D_2->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[2]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[2]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[2]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[2]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[2]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=2.0 && bMesonPt<3.0){
				fMCTruth_Bminus->Fill(19);
				hSgn_CosThetaS_3->Fill(cts); hSgn_CosThetaS_D_3->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[3]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[3]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[3]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[3]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[3]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=3.0 && bMesonPt<4.0){
				fMCTruth_Bminus->Fill(20);
				hSgn_CosThetaS_4->Fill(cts); hSgn_CosThetaS_D_4->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[4]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[4]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[4]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[4]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[4]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=4.0 && bMesonPt<5.0){
				fMCTruth_Bminus->Fill(21);
				hSgn_CosThetaS_5->Fill(cts); hSgn_CosThetaS_D_5->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[5]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[5]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[5]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[5]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[5]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=5.0 && bMesonPt<6.0){
				fMCTruth_Bminus->Fill(22);
				hSgn_CosThetaS_6->Fill(cts); hSgn_CosThetaS_D_6->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[6]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[6]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[6]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[6]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[6]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=6.0 && bMesonPt<8.0){
				fMCTruth_Bminus->Fill(23);
				hSgn_CosThetaS_7->Fill(cts); hSgn_CosThetaS_D_7->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[7]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[7]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[7]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[7]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[7]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=8.0 && bMesonPt<12.){
				fMCTruth_Bminus->Fill(24);
				hSgn_CosThetaS_8->Fill(cts); hSgn_CosThetaS_D_8->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[8]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[8]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[8]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[8]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[8]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=12. && bMesonPt<16.){
				fMCTruth_Bminus->Fill(25);
				hSgn_CosThetaS_9->Fill(cts); hSgn_CosThetaS_D_9->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[9]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[9]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[9]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[9]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[9]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=16. && bMesonPt<20.){
				fMCTruth_Bminus->Fill(26);
				hSgn_CosThetaS_10->Fill(cts); hSgn_CosThetaS_D_10->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[10]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[10]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[10]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[10]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[10]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=20. && bMesonPt<24.){
				fMCTruth_Bminus->Fill(27);
				hSgn_CosThetaS_11->Fill(cts); hSgn_CosThetaS_D_11->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[11]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[11]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[11]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
				}
				else{
					fillthis="hComp_byHand_PtBins[11]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[11]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
			if(bMesonPt>=24.){
				fMCTruth_Bminus->Fill(28);	
				hSgn_CosThetaS_12->Fill(cts); hSgn_CosThetaS_D_12->Fill(ctsD);
				fillthis="hSgn_Ctp_PtBin[12]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(cos);
				if(TMath::Abs(pdg1)==421){
					fillthis="hSgn_EtaD0_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());
					fillthis="hSgn_EtaPi_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
				}
				else{
					fillthis="hSgn_EtaD0_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter2->Eta());
					fillthis="hSgn_EtaPi_PtBin[12]";
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(Bdaughter1->Eta());	
				}
				if(fromBquark==kTRUE){
					fillthis="hComp_bQuark_PtBins[12]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());				
				}
				else{
					fillthis="hComp_byHand_PtBins[12]";
					((TH1F*)(fComparList->FindObject(fillthis)))->Fill(BmesonCandidate->Y());
					fillthis="DecayLenghtBvsD_[12]";
					((TH2F*)(fComparList->FindObject(fillthis)))->Fill(DdecayLength,BdecayLength);
				}
			}
		}
		else{
			Bdaughter1=0;
			Bdaughter2=0;
			part3=0;
			part4=0;
			BmesonCandidate=0;
			// cout << "Rejected: Label is: " << ic << endl;
			continue;
		}
		
		// cout << "Mother is: " << pdgBcand << endl;
		// cout << "Label is: " << ic << endl;
		// cout << "Decays into : " << pdg1 << " + " << pdg2 << endl;
		// cout << "Further into: " << pdg3 << " + " << pdg4 << endl;
		// cout << d1 << " Bdaughter1: " << pdg1 << " Pt: " << Bdaughter1->Pt() << " Eta: " << Bdaughter1->Eta() << endl; 
		// cout << d2 << " Bdaughter2: " << pdg2 << " Pt: " << Bdaughter2->Pt() << " Eta: " << Bdaughter2->Eta() << endl; 
		// cout << d3 << " Part3: " << pdg3 << " Pt: " << part3->Pt() << " Eta: " << part3->Eta() << endl; 
		// cout << d4 << " Part4: " << pdg4 << " Pt: " << part4->Pt() << " Eta: " << part4->Eta() << endl; 
		// cout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
		candidateCounter++;
		// Post the data
	}
	if (!allDaughtersInRange)fMCTruth_Bplus->Fill(31,1.);
	PostData(1,fOutputList);
	PostData(2,fComparList);

	// cout << ">>>>>>>>> B selected in this event: " << candidateCounter << endl;
	return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEFindBhadMC::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
	//
	// Checking if particle is in fiducial acceptance region 
	//
	
	if(pt > 5.) {
		// applying cut for pt > 5 GeV
		AliDebug(2,Form("pt of particle = %f (> 5), cutting at |y| < 0.8\n",pt)); 
		if (TMath::Abs(y) > 0.8){
			return kFALSE;
		}
	} else {
		// appliying smooth cut for pt < 5 GeV
		Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
		Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
		AliDebug(2,Form("pt of particle = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
		if (y < minFiducialY || y > maxFiducialY){
			return kFALSE;
		}
	}
	
	return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskSEFindBhadMC::Terminate(Option_t */*option*/)
{
//	TCanvas *c1 = new TCanvas();
//	fMCTruth_Bplus->Draw("htext0");
//	TCanvas *c2 = new TCanvas();
//	fMCTruth_Bminus->Draw("htext0");
	// Terminate analysis
	//
}
//____________________________________________________________________________
Int_t AliAnalysisTaskSEFindBhadMC::IsTrackInjected(Int_t lab,AliAODMCHeader *header,TClonesArray *arrayMC){
	
	if(lab<0) return 0;
	AliVertexingHFUtils* ggg=new  AliVertexingHFUtils();

	TString bbb=ggg->GetGenerator(lab,header);
	TString empty="";
	//cout << " FIRST CALL " << bbb << endl;
	
	while(bbb.IsWhitespace()){
		if(lab>=arrayMC->GetSize()){delete ggg; return 1;}
		AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(TMath::Abs(lab));
		if(!mcpart){delete ggg; return 1;}
		Int_t mother = mcpart->GetMother();
		lab=mother;
		bbb=ggg->GetGenerator(mother,header);
	//	cout << "Testing " << bbb << " " << endl;
	}
	//cout << " FINAL CALL " << bbb << endl;
	
	if(bbb.Contains("ijing")){delete ggg; return 0;}
	
	delete ggg;
	return 1;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEFindBhadMC::CalculateDecayLength(AliAODMCParticle *MCParticle, AliAODVertex *aodVtx) const{		
	//
	// Decay length of B meson, evaluated via distance between primary vertex and production vertex of daughter
	// or 
	// Decay length of D meson, evaluated via distance between primary vertex and production vertex of daughter
	//
	
	Double_t daughterVtx[3]={0,0,0};
	
	MCParticle->XvYvZv(daughterVtx); // cm
	
	Double_t xVtx = aodVtx->GetX();
	Double_t yVtx = aodVtx->GetY();
	Double_t zVtx = aodVtx->GetZ();
	
	//printf(" --- Primary Vertex at \n x: %f \n y: %f \n z: %f \n",xVtx,yVtx,zVtx);
	//printf(" --- Decay Vertex at \n x: %f \n y: %f \n z: %f \n",daughterVtx[0],daughterVtx[1],daughterVtx[2]);
	
	return TMath::Sqrt((daughterVtx[0]-xVtx)
					   *(daughterVtx[0]-xVtx)
					   +(daughterVtx[1]-yVtx)
					   *(daughterVtx[1]-yVtx)
					   +(daughterVtx[2]-zVtx)
					   *(daughterVtx[2]-zVtx));
	
	
}
