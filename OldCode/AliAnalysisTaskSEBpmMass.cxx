/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *::
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskSEBpmMass.cxx 53668 2011-12-16 14:14:00Z prino $ */

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for B candidates invariant mass histogram
// and comparison with the MC truth and cut variables distributions.
//
// Authors: Johannes Stiller, jstiller@cern.ch
// Track Rotations: Alessandro Grelli
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2I.h>
#include "THnSparse.h"
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
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
#include "AliAnalysisTaskSEBpmMass.h"
#include "AliNormalizationCounter.h"
#include "AliVertexerTracks.h"
#include "AliNeutralTrackParam.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliGenHijingEventHeader.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliVertexingHFUtils.h"

ClassImp(AliAnalysisTaskSEBpmMass)

//________________________________________________________________________
AliAnalysisTaskSEBpmMass::AliAnalysisTaskSEBpmMass():
AliAnalysisTaskSE(),
fPIDResponse(0),
fbHIJINGavailable(kFALSE),
fbSelectHIJING(kFALSE),
fbtrackRotation(kFALSE),
fbCentrality(kFALSE),
fhijingGenHeader(0),
fNentries(0),
massWindowB(0),
massWindowD(0),
fCuts(0),
fArray(0),
fReadMC(0),
fKillCandidates(0),
fSys(0),
fBzkG(0.),
fvtx1(0x0),
fOutputList(0),
fOutputListVar(0),
fOutputListMCVar(0),
fOutputListDSandB(0),
fOutputListDSonly(0),
fInvMassBpm(0),
fInvMassDiff(0),
fCosThetaStar(0),
fCosPoinAngle(0),
fCosPoinAngleXY(0),
fTransvMome(0),
fDCA(0),
fTransMomD(0),
fTransMomPi(0),
fd0D(0),
fd0Pi(0),
fProdd0d0(0),
fNormDecLthXY(0),
fImpParXY(0),
fctau(0),
fMCTruth(0),
fMCInvMassBpm(0),
fMCInvMassDiff(0),
fMCCosThetaStar(0),
fMCCosPoinAngle(0),
fMCCosPoinAngleXY(0),
fMCTransvMome(0),
fMCDCA(0),
fMCTransMomD(0),
fMCTransMomPi(0),
fMCd0D(0),
fMCd0Pi(0),
fMCProdd0d0(0),
fMCNormDecLthXY(0),
fMCImpParXY(0),
fMCctau(0),
fRunNumber(0),
hThetaD0(0),
hThetaPi(0),
hThetaD0Ka(0),
hThetaD0Pi(0),
fSelectionVariables(0),
fRotatedCandidates(0),
fmcHeader(0)
,fMcArray(0),
fRot(13),
fAngleFirst(3.14),
fAngle(0.087),
fNDebug(0),
fAODMapSize(0),
fAODMap(0)
{
	
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass::DefaultConstructor \n");
	
	// Default constructor
	massWindowB = 0.5;//0.1;
	massWindowD = 0.10;//0.1;
	
	// Default cut values loaded here
	for(Int_t ptbin=0;ptbin<13;ptbin++){
		fBCtpCut[ptbin] = 0.99;
		fD0PtCut[ptbin] = 0.;
		fDNdlCut[ptbin] = 5.;
		fBIxyCut[ptbin] = 0.008;
		fPiPtCut[ptbin] = 0.5;
	}
	
	fBCtpCut[2] = 0.99650;	fD0PtCut[2] = 1.40;		fDNdlCut[2] = 13.;		fBIxyCut[2] = 0.00635;	fPiPtCut[2] = 1.8;
	fBCtpCut[3] = 0.99800;	fD0PtCut[3] = 0.70;		fDNdlCut[3] = 29.;		fBIxyCut[3] = 0.00455;	fPiPtCut[3] = 1.4;
	fBCtpCut[4] = 0.99875;	fD0PtCut[4] = 1.00;		fDNdlCut[4] = 20.;		fBIxyCut[4] = 0.00625;	fPiPtCut[4] = 1.4;
	fBCtpCut[5] = 0.99975;	fD0PtCut[5] = 0.10;		fDNdlCut[5] = 13.;		fBIxyCut[5] = 0.00080;	fPiPtCut[5] = 1.1;
	fBCtpCut[6] = 0.99975;	fD0PtCut[6] = 0.40;		fDNdlCut[6] = 12.;		fBIxyCut[6] = 0.00470;	fPiPtCut[6] = 1.1;
	fBCtpCut[7] = 0.99890;	fD0PtCut[7] = 1.80;		fDNdlCut[7] = 19.;		fBIxyCut[7] = 0.00200;	fPiPtCut[7] = 1.2;
	fBCtpCut[8] = 0.999875;	fD0PtCut[8] = 0.00;		fDNdlCut[8] = 20.;		fBIxyCut[8] = 0.00620;	fPiPtCut[8] = 0.9;
	fBCtpCut[9] = 0.99980;	fD0PtCut[9] = 1.30;		fDNdlCut[9] = 9.;			fBIxyCut[9] = 0.00650;	fPiPtCut[9] = 1.3;
	fBCtpCut[10] = 0.99850;	fD0PtCut[10] = 2.25;	fDNdlCut[10] = 13.;		fBIxyCut[10] = 0.00350;	fPiPtCut[10] = 1.2;
	fBCtpCut[11] = 0.99800;	fD0PtCut[11] = 3.40;	fDNdlCut[11] = 3.75;	fBIxyCut[11] = 0.00425;	fPiPtCut[11] = 1.1;
	fBCtpCut[12] = 0.99800;	fD0PtCut[12] = 4.40;	fDNdlCut[12] = 5.0;		fBIxyCut[12] = 0.00600;	fPiPtCut[12] = 0.8;
}
//________________________________________________________________________
AliAnalysisTaskSEBpmMass::AliAnalysisTaskSEBpmMass(const char *name,AliRDHFCutsD0toKpi* cuts):
AliAnalysisTaskSE(name),
fPIDResponse(0),
fbHIJINGavailable(kFALSE),
fbSelectHIJING(kFALSE),
fbtrackRotation(kFALSE),
fbCentrality(kFALSE),
fhijingGenHeader(0),
fNentries(0),
massWindowB(0),
massWindowD(0),
fCuts(0),
fArray(0),
fReadMC(0),
fKillCandidates(0),
fSys(0),
fBzkG(0.),
fvtx1(0x0),
fOutputList(0),
fOutputListVar(0),
fOutputListMCVar(0),
fOutputListDSandB(0),
fOutputListDSonly(0),
fInvMassBpm(0),
fInvMassDiff(0),
fCosThetaStar(0),
fCosPoinAngle(0),
fCosPoinAngleXY(0),
fTransvMome(0),
fDCA(0),
fTransMomD(0),
fTransMomPi(0),
fd0D(0),
fd0Pi(0),
fProdd0d0(0),
fNormDecLthXY(0),
fImpParXY(0),
fctau(0),
fMCTruth(0),
fMCInvMassBpm(0),
fMCInvMassDiff(0),
fMCCosThetaStar(0),
fMCCosPoinAngle(0),
fMCCosPoinAngleXY(0),
fMCTransvMome(0),
fMCDCA(0),
fMCTransMomD(0),
fMCTransMomPi(0),
fMCd0D(0),
fMCd0Pi(0),
fMCProdd0d0(0),
fMCNormDecLthXY(0),
fMCImpParXY(0),
fMCctau(0),
fRunNumber(0),
hThetaD0(0),
hThetaPi(0),
hThetaD0Ka(0),
hThetaD0Pi(0),
fSelectionVariables(0)
,fRotatedCandidates(0)
,fmcHeader(0)
,fMcArray(0),
fRot(13),
fAngleFirst(3.14),
fAngle(0.087),
fNDebug(0),
fAODMapSize(0),
fAODMap(0)
{
	
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass::AdvancedConstructor \n");
	
	massWindowB = 0.5;//0.1;
	massWindowD = 0.10;//0.1;
	
	// Default cut values loaded here
	for(Int_t ptbin=0;ptbin<13;ptbin++){
		fBCtpCut[ptbin] = 0.99;
		fD0PtCut[ptbin] = 0.;
		fDNdlCut[ptbin] = 5.;
		fBIxyCut[ptbin] = 0.008;
		fPiPtCut[ptbin] = 0.5;
	}
	
	fBCtpCut[2] = 0.99650;	fD0PtCut[2] = 1.40;		fDNdlCut[2] = 13.;		fBIxyCut[2] = 0.00635;	fPiPtCut[2] = 1.8;
	fBCtpCut[3] = 0.99800;	fD0PtCut[3] = 0.70;		fDNdlCut[3] = 29.;		fBIxyCut[3] = 0.00455;	fPiPtCut[3] = 1.4;
	fBCtpCut[4] = 0.99875;	fD0PtCut[4] = 1.00;		fDNdlCut[4] = 20.;		fBIxyCut[4] = 0.00625;	fPiPtCut[4] = 1.4;
	fBCtpCut[5] = 0.99975;	fD0PtCut[5] = 0.10;		fDNdlCut[5] = 13.;		fBIxyCut[5] = 0.00080;	fPiPtCut[5] = 1.1;
	fBCtpCut[6] = 0.99975;	fD0PtCut[6] = 0.40;		fDNdlCut[6] = 12.;		fBIxyCut[6] = 0.00470;	fPiPtCut[6] = 1.1;
	fBCtpCut[7] = 0.99890;	fD0PtCut[7] = 1.80;		fDNdlCut[7] = 19.;		fBIxyCut[7] = 0.00200;	fPiPtCut[7] = 1.2;
	fBCtpCut[8] = 0.999875;	fD0PtCut[8] = 0.00;		fDNdlCut[8] = 20.;		fBIxyCut[8] = 0.00620;	fPiPtCut[8] = 0.9;
	fBCtpCut[9] = 0.99980;	fD0PtCut[9] = 1.30;		fDNdlCut[9] = 9.;		fBIxyCut[9] = 0.00650;	fPiPtCut[9] = 1.3;
	fBCtpCut[10] = 0.99850;	fD0PtCut[10] = 2.25;	fDNdlCut[10] = 13.;		fBIxyCut[10] = 0.00350;	fPiPtCut[10] = 1.2;
	fBCtpCut[11] = 0.99800;	fD0PtCut[11] = 3.40;	fDNdlCut[11] = 3.75;	fBIxyCut[11] = 0.00425;	fPiPtCut[11] = 1.1;
	fBCtpCut[12] = 0.99800;	fD0PtCut[12] = 4.40;	fDNdlCut[12] = 5.0;		fBIxyCut[12] = 0.00600;	fPiPtCut[12] = 0.8;
	
	fCuts=cuts;
	
	DefineOutput(1,TList::Class());  //My private output
	DefineOutput(2,TList::Class());  //My private output
	DefineOutput(3,TList::Class());  //My private output
	DefineOutput(4,TList::Class());  //My private output
	DefineOutput(5,TList::Class());  //My private output
	DefineOutput(6,TTree::Class());  //My private output
	DefineOutput(7,TTree::Class());  //My private output
	DefineOutput(8,AliRDHFCutsD0toKpi::Class());  //My private output
}
//________________________________________________________________________
AliAnalysisTaskSEBpmMass::~AliAnalysisTaskSEBpmMass()
{
	
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass::~AliAnalysisTaskSEBpmMass() \n");
	
	if(fAODMap) { delete [] fAODMap; fAODMap=0; }
	
	
	if (fOutputList) {
		delete fOutputList;
		fOutputList = 0;
	}
	if (fOutputListVar) {
		delete fOutputListVar;
		fOutputListVar = 0;
	}
	if (fOutputListMCVar) {
		delete fOutputListMCVar;
		fOutputListMCVar = 0;
	}
	if (fOutputListDSandB) {
		delete fOutputListDSandB;
		fOutputListDSandB = 0;
	}
	if (fOutputListDSonly) {
		delete fOutputListDSonly;
		fOutputListDSonly = 0;
	}
	if(fvtx1){
		delete fvtx1;
		fvtx1=0;
	}
}
//________________________________________________________________________
void AliAnalysisTaskSEBpmMass::Init()
{
	// Initialization
	
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass::Init() \n");
	
	AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
	const char* nameoutput=GetOutputSlot(8)->GetContainer()->GetName();
	copyfCuts->SetName(nameoutput);
	// Post the data
	PostData(8,copyfCuts);
	
	return;
}
//________________________________________________________________________
void AliAnalysisTaskSEBpmMass::UserCreateOutputObjects()
{
	
	// Create the output container
	//
	
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass::UserCreateOutputObjects() \n");
	
	// Load PID Response
	AliAnalysisManager *man							= AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler	= (AliInputEventHandler*) (man->GetInputEventHandler());
	fPIDResponse												= inputHandler->GetPIDResponse();
	
	// Tree with selection variables
	// Open Output File, since TTree is created
	OpenFile(6);
	
	// Called once
	if(fSelectionVariables != NULL){
		delete fSelectionVariables;
		fSelectionVariables = NULL;
	}
	if(fSelectionVariables == NULL){
		fSelectionVariables = new TTree("fSelectionVariables","fSelectionVariables");
		fSelectionVariables->Branch("bBInvMass",&bBInvMass);
		fSelectionVariables->Branch("bBInvMassDiff",&bBInvMassDiff);
		fSelectionVariables->Branch("bPiPhysPrimary",&bPiPhysPrimary);
		fSelectionVariables->Branch("bDDgh0PhysPrimary",&bDDgh0PhysPrimary);
		fSelectionVariables->Branch("bDDgh1PhysPrimary",&bDDgh1PhysPrimary);
		fSelectionVariables->Branch("bBPiTOFnSigma",&bBPiTOFnSigma);
		fSelectionVariables->Branch("bBPiTPCnSigma",&bBPiTPCnSigma);
		fSelectionVariables->Branch("bPiWeak",&bPiWeak);
		fSelectionVariables->Branch("bPiPdgCode",&bPiPdgCode);
		fSelectionVariables->Branch("bDDgh0Weak",&bDDgh0Weak);
		fSelectionVariables->Branch("bDDgh0PdgCode",&bDDgh0PdgCode);
		fSelectionVariables->Branch("bDDgh1Weak",&bDDgh1Weak);
		fSelectionVariables->Branch("bDDgh1PdgCode",&bDDgh1PdgCode);
		//		fSelectionVariables->Branch("bPiMaterial",&bPiMaterial);
		//		fSelectionVariables->Branch("bDDgh0Material",&bDDgh0Material);
		//		fSelectionVariables->Branch("bDDgh1Material",&bDDgh1Material);
		fSelectionVariables->Branch("bBpT",&bBpT);
		fSelectionVariables->Branch("bBctp",&bBctp);
		fSelectionVariables->Branch("bBctpXY",&bBctpXY);
		fSelectionVariables->Branch("bBd0XY",&bBd0XY);
		fSelectionVariables->Branch("bBPipt",&bBPipt);
		fSelectionVariables->Branch("bd0d0",&bd0d0);
		fSelectionVariables->Branch("bBndlXY",&bBndlXY);
		fSelectionVariables->Branch("bVtxChi2",&bVtxChi2);
		fSelectionVariables->Branch("bBdl",&bBdl);
		fSelectionVariables->Branch("bBNormIpPi",&bBNormIpPi);
		fSelectionVariables->Branch("bBNormIpD0",&bBNormIpD0);
		fSelectionVariables->Branch("bDpT",&bDpT);
		fSelectionVariables->Branch("bDndl",&bDndl);
		fSelectionVariables->Branch("bDdl",&bDdl);
		fSelectionVariables->Branch("bDInvMass",&bDInvMass);
		fSelectionVariables->Branch("bDdca",&bDdca);
		fSelectionVariables->Branch("bDcts",&bDcts);
		fSelectionVariables->Branch("bDKapt",&bDKapt);
		fSelectionVariables->Branch("bDPipt",&bDPipt);
		fSelectionVariables->Branch("bDd0Ka",&bDd0Ka);
		fSelectionVariables->Branch("bDd0Pi",&bDd0Pi);
		fSelectionVariables->Branch("bDd0d0",&bDd0d0);
		fSelectionVariables->Branch("bDctp",&bDctp);
		fSelectionVariables->Branch("bDctpXY",&bDctpXY);
		fSelectionVariables->Branch("bDndlXY",&bDndlXY);
		 fSelectionVariables->Branch("dcaZBPi",&dcaZBPi);
		 fSelectionVariables->Branch("dcaZD0",&dcaZD0);
		 fSelectionVariables->Branch("dcaZD0Ka",&dcaZD0Ka);
		 fSelectionVariables->Branch("dcaZD0Pi",&dcaZD0Pi);
		
	}
	
	// Tree with background information
	// Open Output File, since TTree is created
	OpenFile(7);
	
	if(fRotatedCandidates != NULL){
		delete fRotatedCandidates;
		fRotatedCandidates = NULL;
	}
	if(fbtrackRotation){
		if(fRotatedCandidates == NULL){
			fRotatedCandidates = new TTree("fRotatedCandidates","fRotatedCandidates");
			fRotatedCandidates->Branch("bRotBInvMass",&bRotBInvMass);
			fRotatedCandidates->Branch("bRotBInvMassDiff",&bRotBInvMassDiff);
			fRotatedCandidates->Branch("bRotBpT",&bRotBpT);
			fRotatedCandidates->Branch("bRotBctp",&bRotBctp);
			fRotatedCandidates->Branch("bRotBctpXY",&bRotBctpXY);
			fRotatedCandidates->Branch("bRotBd0XY",&bRotBd0XY);
			fRotatedCandidates->Branch("bRotBPipt",&bRotBPipt);
			fRotatedCandidates->Branch("bRotd0d0",&bRotd0d0);
			fRotatedCandidates->Branch("bRotBndlXY",&bRotBndlXY);
			fRotatedCandidates->Branch("bRotVtxChi2",&bRotVtxChi2);
			fRotatedCandidates->Branch("bRotBdl",&bRotBdl);
			fRotatedCandidates->Branch("bRotDpT",&bRotDpT);
			fRotatedCandidates->Branch("bRotDndl",&bRotDndl);
			fRotatedCandidates->Branch("bRotDdl",&bRotDdl);
			fRotatedCandidates->Branch("bRotDInvMass",&bRotDInvMass);
			fRotatedCandidates->Branch("bRotDdca",&bRotDdca);
			fRotatedCandidates->Branch("bRotDcts",&bRotDcts);
			fRotatedCandidates->Branch("bRotDKapt",&bRotDKapt);
			fRotatedCandidates->Branch("bRotDPipt",&bRotDPipt);
			fRotatedCandidates->Branch("bRotDd0Ka",&bRotDd0Ka);
			fRotatedCandidates->Branch("bRotDd0Pi",&bRotDd0Pi);
			fRotatedCandidates->Branch("bRotDd0d0",&bRotDd0d0);
			fRotatedCandidates->Branch("bRotDctp",&bRotDctp);
			fRotatedCandidates->Branch("bRotDctpXY",&bRotDctpXY);
			fRotatedCandidates->Branch("bRotDndlXY",&bRotDndlXY);
			
			 fRotatedCandidates->Branch("dcaZBPi",&dcaZBPi);
			 fRotatedCandidates->Branch("dcaZD0",&dcaZD0);
			 fRotatedCandidates->Branch("dcaZD0Ka",&dcaZD0Ka);
			 fRotatedCandidates->Branch("dcaZD0Pi",&dcaZD0Pi);
			
		}
	}
	else {fRotatedCandidates=0x0;}
	
	// Johannes' Histos
	fOutputList = new TList();
	fOutputList->SetOwner();
	
	fOutputListVar = new TList();
	fOutputListVar->SetOwner();
	fOutputListVar->SetName("listVarPtBins");
	
	fOutputListMCVar = new TList();
	fOutputListMCVar->SetOwner();
	fOutputListMCVar->SetName("listVarMCPtBins");
	
	fOutputListDSandB = new TList();
	fOutputListDSandB->SetOwner();
	fOutputListDSandB->SetName("listDSandB");
	
	fOutputListDSonly = new TList();
	fOutputListDSonly->SetOwner();
	fOutputListDSonly->SetName("listDSonly");
	
	// D variable histograms - data
	TH1F* InvMassD = new TH1F("InvMassD","m(KaPi); m(KaPi) (GeV/c^{2}); entries",120,1.864-massWindowD,1.864+massWindowD);
	TH1F* InvMassDbar = new TH1F("InvMassDbar","m(KaPi); m(KaPi) (GeV/c^{2}); entries",120,1.864-massWindowD,1.864+massWindowD);
	TH1F* DCosThetaStar = new TH1F("DCosThetaStar","D cos(#theta^{*}) w.r.t. Pion; ; entries",200,-1.0,1.0);
	TH1F* DCosPoinAngle = new TH1F("DCosPoinAngle","D cos(#theta_{point}); ; entries",200,-1.0,1.0);
	TH1F* DCosPoinAngleXY = new TH1F("DCosPoinAngleXY","D cos(#theta_{point}) XY; ; entries",200,-1.0,1.0);
	TH1F* DTransvMome = new TH1F("DTransvMome","D transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	TH1F* DDCA = new TH1F("DDCA","D dca; dca (cm) ; entries",100,0.0,10.0);
	TH1F* DTransMomKa  = new TH1F("DTransMomKa","Ka^{+/-} from D transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	TH1F* DTransMomPi = new TH1F("DTransMomPi","Pi^{+/-} from D transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	TH1F* Dd0Ka  = new TH1F("Dd0Ka","d0 of Ka^{+/-} from D; dca (cm) ; entries",200,-1.0,1.0);
	TH1F* Dd0Pi = new TH1F("Dd0Pi","d0 of Pi^{+/-} from D ; dca (cm) ; entries",200,-1.0,1.0);
	TH1F* DProdd0d0 = new TH1F("DProdd0d0","d0xd0 of Pi^{+/-} and Ka^{+/-}  ; (cm^{2}) ; entries",200,-0.001,0.001);
	TH1F* Ddeclgt = new TH1F("Ddeclgt", "Decay Length distribution; Decay Length [cm]",100,0,10);
	TH1F* Dnormdeclgt = new TH1F("Dnormdeclgt", "Normalized Decay Length distribution;(Decay Length/Err) ",100,0,100.);
	TH1F* Ddeclgtxy = new TH1F("Ddeclgtxy","Decay Length XY distribution;Decay Length XY (cm)",100,0,10);;
	TH1F* DnormdeclgtXY = new TH1F("DnormdeclgtXY","normalized decay length of D ; l (cm) ; entries",100,0,100);
	TH1F* DImpParXY = new TH1F("DImpParXY","impact parameter XY of D ; l (cm) ; entries",200,-1000.0,1000.0);
	TH1F* Dctau = new TH1F("Dctau","c#tau of D ; c#tau (cm) ; entries",200,-0.001,1.001);
	
	fOutputListDSandB->Add(InvMassD);
	fOutputListDSandB->Add(InvMassDbar);
	fOutputListDSandB->Add(DCosThetaStar);
	fOutputListDSandB->Add(DCosPoinAngle);
	fOutputListDSandB->Add(DCosPoinAngleXY);
	fOutputListDSandB->Add(DTransvMome);
	fOutputListDSandB->Add(DDCA);
	fOutputListDSandB->Add(DTransMomKa);
	fOutputListDSandB->Add(DTransMomPi);
	fOutputListDSandB->Add(Dd0Ka);
	fOutputListDSandB->Add(Dd0Pi);
	fOutputListDSandB->Add(DProdd0d0);
	fOutputListDSandB->Add(Ddeclgt);
	fOutputListDSandB->Add(Dnormdeclgt);
	fOutputListDSandB->Add(Ddeclgtxy);
	fOutputListDSandB->Add(DnormdeclgtXY);
	fOutputListDSandB->Add(DImpParXY);
	fOutputListDSandB->Add(Dctau);
	
	if(fReadMC){
		
		// D variable histograms - data
		TH1F* MCInvMassD = new TH1F("MCInvMassD","m(KaPi); m(KaPi) (GeV/c^{2}); entries",120,1.864-massWindowD,1.864+massWindowD);
		TH1F* MCInvMassDbar = new TH1F("MCInvMassDbar","m(KaPi); m(KaPi) (GeV/c^{2}); entries",120,1.864-massWindowD,1.864+massWindowD);
		TH1F* MCDCosThetaStar = new TH1F("MCDCosThetaStar","D cos(#theta^{*}) w.r.t. Pion; ; entries",200,-1.0,1.0);
		TH1F* MCDCosPoinAngle = new TH1F("MCDCosPoinAngle","D cos(#theta_{point}); ; entries",200,-1.0,1.0);
		TH1F* MCDCosPoinAngleXY = new TH1F("MCDCosPoinAngleXY","D cos(#theta_{point}) XY; ; entries",200,-1.0,1.0);
		TH1F* MCDTransvMome = new TH1F("MCDTransvMome","D transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* MCDDCA = new TH1F("MCDDCA","D dca; dca (cm) ; entries",100,0.0,10.0);
		TH1F* MCDTransMomKa  = new TH1F("MCDTransMomKa","Ka^{+/-} from D transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* MCDTransMomPi = new TH1F("MCDTransMomPi","Pi^{+/-} from D transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* MCDd0Ka  = new TH1F("MCDd0Ka","d0 of Ka^{+/-} from D; dca (cm) ; entries",200,-1.0,1.0);
		TH1F* MCDd0Pi = new TH1F("MCDd0Pi","d0 of Pi^{+/-} from D ; dca (cm) ; entries",200,-1.0,1.0);
		TH1F* MCDProdd0d0 = new TH1F("MCDProdd0d0","d0xd0 of Pi^{+/-} and Ka^{+/-}  ; (cm^{2}) ; entries",200,-0.001,0.001);
		TH1F* MCDdeclgt = new TH1F("MCDdeclgt", "Decay Length distribution; Decay Length (cm)",100,0,10);
		TH1F* MCDnormdeclgt = new TH1F("MCDnormdeclgt", "Normalized Decay Length distribution;(Decay Length/Err) ",100,0,100.);
		TH1F* MCDdeclgtxy = new TH1F("MCDdeclgtxy","Decay Length XY distribution;Decay Length XY [cm]",100,0,1.);
		TH1F* MCDnormdeclgtXY = new TH1F("MCDnormdeclgtXY","normalized decay length of D ; l (cm) ; entries",100,0,100);
		TH1F* MCDImpParXY = new TH1F("MCDImpParXY","impact parameter XY of D ; l (cm) ; entries",200,-1000.0,1000.0);
		TH1F* MCDctau = new TH1F("MCDctau","c#tau of D ; c#tau (cm) ; entries",200,-0.001,1.001);
		
		fOutputListDSonly->Add(MCInvMassD);
		fOutputListDSonly->Add(MCInvMassDbar);
		fOutputListDSonly->Add(MCDCosThetaStar);
		fOutputListDSonly->Add(MCDCosPoinAngle);
		fOutputListDSonly->Add(MCDCosPoinAngleXY);
		fOutputListDSonly->Add(MCDTransvMome);
		fOutputListDSonly->Add(MCDDCA);
		fOutputListDSonly->Add(MCDTransMomKa);
		fOutputListDSonly->Add(MCDTransMomPi);
		fOutputListDSonly->Add(MCDd0Ka);
		fOutputListDSonly->Add(MCDd0Pi);
		fOutputListDSonly->Add(MCDProdd0d0);
		fOutputListDSonly->Add(MCDdeclgt);
		fOutputListDSonly->Add(MCDnormdeclgt);
		fOutputListDSonly->Add(MCDdeclgtxy);
		fOutputListDSonly->Add(MCDnormdeclgtXY);
		fOutputListDSonly->Add(MCDImpParXY);
		fOutputListDSonly->Add(MCDctau);
	}
	
	// B variable histograms - data
	fInvMassBpm = new TH1F("fInvMassBpm","m(PiKPi); m(PiKPi) (GeV/c^{2}); entries",(Int_t)(2*100*massWindowB),5.279-massWindowB,5.279+massWindowB);
	fInvMassDiff = new TH1F("fInvMassDiff","m(D^{0}#pi_{B})-m(K#pi_{D}); m(D^{0}#pi_{B})-m(K#pi_{D}) (GeV/#it{c}^{2}); entries",(Int_t)(2*100*massWindowB),3.415-massWindowB,3.415+massWindowB);
	fCosThetaStar = new TH1F("fCosThetaStar","B^{+/-} cos(#theta^{*}) w.r.t. Pion; ; entries",200,-1.0,1.0);
	fCosPoinAngle = new TH1F("fCosPoinAngle","B^{+/-} cos(#theta_{point}); ; entries",200,-1.0,1.0);
	fCosPoinAngleXY = new TH1F("fCosPoinAngleXY","B^{+/-} cos(#theta_{point}) XY; ; entries",200,-1.0,1.0);
	fTransvMome = new TH1F("fTransvMome","B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	fDCA = new TH1F("fDCA","B^{+/-} dca; dca (cm) ; entries",100,0.0,10.0);
	fTransMomD  = new TH1F("fTransMomD","D^{0}(bar) from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	fTransMomPi = new TH1F("fTransMomPi","Pi^{+/-} from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	fd0D  = new TH1F("fd0D","d0 of D^{0}(bar) from B^{+/-}; dca (cm) ; entries",200,-1.0,1.0);
	fd0Pi = new TH1F("fd0Pi","d0 of Pi^{+/-} from B^{+/-} ; dca (cm) ; entries",200,-1.0,1.0);
	fProdd0d0 = new TH1F("fProdd0d0","d0xd0 of Pi^{+/-} and D^{0}(bar) ; (cm^{2}) ; entries",200,-0.001,0.001);
	// new***
	fNormDecLthXY = new TH1F("fNormDecLthXY","normalized decay length of B^{+/-} ; l (cm) ; entries",2000,-0.0,20.0);
	
	TH1F *hdeclgt = new TH1F("hdeclgt", "Decay Length distribution; Decay Length [cm]",100,0,0.1);
	TH1F *hnormdeclgt = new TH1F("hnormdeclgt", "Normalized Decay Length distribution;(Decay Length/Err) ",100,0,100.);
	TH1F *hdeclgtxy = new TH1F("hdeclgtxy","Decay Length XY distribution;Decay Length XY [cm]",100,0,1.);
	TH1F *hnormdeclgtXY = new TH1F("hnormdeclgtXY","normalized decay length of B^{+/-} ; l (cm) ; entries",100,0,10);
	TH2F *hdeclxyd0d0 = new TH2F("hdeclxyd0d0","Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",200,0,0.15,200,-0.001,0.001);
	TH2F *hnormdeclxyd0d0 = new TH2F("hnormdeclxyd0d0","Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",200,0,6,200,-0.001,0.001);
	
	fImpParXY = new TH1F("fImpParXY","impact parameter XY of B^{+/-} ; l (cm) ; entries",200,-1000.0,1000.0);
	fctau = new TH1F("fctau","c#tau of B^{+/-} ; c#tau (cm) ; entries",1000,-0.001,1.001);
	fRunNumber = new TH1F("fRunNumber","run number with MC B ; run X ; entries",300000,1000,400000);
	TH1F *hcentrality = new TH1F("hcentrality","centrality",500,0,5000.);
	TH1F *hSDpt  = new TH1F("hSDpt","D^{0}(bar) from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
	TH2F *hBgkPionPID = new TH2F("hBgkPionPID","TPC Vs TOF PID when available; TPC # of |#sigma|; TOF # of |#sigma|",100,0,100,100,0,100);
	TH1F *Bfiducial_before = new TH1F("Bfiducial_before","Bfiducial_before",200,-1.,+1.);
	TH1F *Bfiducial_after = new TH1F("Bfiducial_after","Bfiducial_after",200,-1.,+1.);
	TH2F *hSgnPionPID = new TH2F("hSgnPionPID","TPC Vs TOF PID when available; TPC # of |#sigma|; TOF # of |#sigma|",100,0,100,100,0,100);
	
	fMCTruth = new TH1F("fMCTruth","candidates; ; Entries",30,0,30);
	fMCTruth->GetXaxis()->SetBinLabel(1,"Called");
	fMCTruth->GetXaxis()->SetBinLabel(2,"D0");
	fMCTruth->GetXaxis()->SetBinLabel(3,"D0bar");
	fMCTruth->GetXaxis()->SetBinLabel(4,"either");
	fMCTruth->GetXaxis()->SetBinLabel(5,"D0 in MC");
	fMCTruth->GetXaxis()->SetBinLabel(6,"D0bar in MC");
	fMCTruth->GetXaxis()->SetBinLabel(7,"either->D0");
	fMCTruth->GetXaxis()->SetBinLabel(8,"either->D0bar");
	fMCTruth->GetXaxis()->SetBinLabel(9,"MCtrue D0fB");
	fMCTruth->GetXaxis()->SetBinLabel(10,"MCtrue D0barfB");
	fMCTruth->GetXaxis()->SetBinLabel(12,"FillMChistosCalls");
	fMCTruth->GetXaxis()->SetBinLabel(13,"MC D0");
	fMCTruth->GetXaxis()->SetBinLabel(14,"MC D0bar");
	fMCTruth->GetXaxis()->SetBinLabel(15,"Filled D0 MC");
	fMCTruth->GetXaxis()->SetBinLabel(16,"Filled DObar MC");
	fMCTruth->GetXaxis()->SetBinLabel(17,"either");
	fMCTruth->GetXaxis()->SetBinLabel(18,"Filled e D0 MC");
	fMCTruth->GetXaxis()->SetBinLabel(19,"Filled e D0bar MC");
	fMCTruth->GetXaxis()->SetBinLabel(20,"Finished Fill MC");
	
	fOutputList->Add(fInvMassBpm);
	fOutputList->Add(fInvMassDiff);
	fOutputList->Add(fCosThetaStar);
	fOutputList->Add(fCosPoinAngle);
	fOutputList->Add(fCosPoinAngleXY);
	fOutputList->Add(fTransvMome);
	fOutputList->Add(fDCA);
	fOutputList->Add(fTransMomD);
	fOutputList->Add(fTransMomPi);
	fOutputList->Add(fd0D);
	fOutputList->Add(fd0Pi);
	fOutputList->Add(fProdd0d0);
	fOutputList->Add(fNormDecLthXY);
	fOutputList->Add(hdeclgt);
	fOutputList->Add(hnormdeclgt);
	fOutputList->Add(hdeclgtxy);
	fOutputList->Add(hnormdeclgtXY);
	fOutputList->Add(hdeclxyd0d0);
	fOutputList->Add(hnormdeclxyd0d0);
	fOutputList->Add(fImpParXY);
	fOutputList->Add(fctau);
	fOutputList->Add(hSDpt);
	
	fOutputList->Add(hcentrality);
	fOutputList->Add(hBgkPionPID);
	fOutputList->Add(hSgnPionPID);
	fOutputList->Add(fRunNumber);
	fOutputList->Add(Bfiducial_before);
	fOutputList->Add(Bfiducial_after);
	
	// B variable histograms - MC
	if(fReadMC){
		fMCInvMassBpm = new TH1F("fMCInvMassBpm","m(PiKPi); m(PiKPi) (GeV/c^{2}); entries",(Int_t)(2*100*massWindowB),5.279-massWindowB,5.279+massWindowB);
		fMCInvMassDiff = new TH1F("fMCInvMassDiff","m(D^{0}#pi_{B})-m(K#pi_{D}) - MC; m(D^{0}#pi_{B})-m(K#pi_{D}) (GeV/#it{c}^{2}); entries",(Int_t)(2*100*massWindowB),3.415-massWindowB,3.415+massWindowB);
		fMCCosThetaStar = new TH1F("fMCCosThetaStar","B^{+/-} cos(#theta^{*}) w.r.t. Pion; ; entries",200,-1.0,1.0);
		fMCCosPoinAngle = new TH1F("fCMCosPoinAngle","B^{+/-} cos(#theta_{point}); ; entries",200,-1.0,1.0);
		fMCCosPoinAngleXY = new TH1F("fMCCosPoinAngleXY","B^{+/-} cos(#theta_{point}) XY; ; entries",200,-1.0,1.0);;
		fMCTransvMome = new TH1F("fMCTransvMome","B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		fMCDCA = new TH1F("fMCDCA","B^{+/-} dca; dca (cm) ; entries",100,0.0,10.0);
		fMCTransMomD  = new TH1F("fMCTransMomD","D^{0}(bar) from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		fMCTransMomPi = new TH1F("fMCTransMomPi","Pi^{+/-} from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		fMCd0D  = new TH1F("fMCd0D","d0 of D^{0}(bar) from B^{+/-}; dca (cm) ; entries",200,-1.0,1.0);
		fMCd0Pi = new TH1F("fMCd0Pi","d0 of Pi^{+/-} from B^{+/-} ; dca (cm) ; entries",200,-1.0,1.0);
		fMCProdd0d0 = new TH1F("fMCProdd0d0","d0xd0 of Pi^{+/-} and D^{0}(bar) ; dca (cm^{2}) ; entries",200,-0.001,0.001);
		fMCNormDecLthXY = new TH1F("fMCNormDecLthXY","normalized decay length of B^{+/-} ; (cm) ; entries",2000,-0.0,20.0);
		TH1F *hMCdeclgt = new TH1F("hMCdeclgt", "Decay Length distribution; Decay Length [cm]",100,0,0.1);
		TH1F *hMCnormdeclgt = new TH1F("hMCnormdeclgt", "Normalized Decay Length distribution;(Decay Length/Err) ",100,0,100.);
		TH1F *hMCdeclgtxy = new TH1F("hMCdeclgtxy","Decay Length XY distribution;Decay Length XY [cm]",100,0,1);
		TH1F *hMCnormdeclgtXY = new TH1F("hMCnormdeclgtXY","normalized decay length of B^{+/-} ; l (cm) ; entries",100,0,10);
		TH2F *hMCdeclxyd0d0 = new TH2F("hMCdeclxyd0d0","Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",500,0,1.,200,-0.001,0.001);
		TH2F *hMCnormdeclxyd0d0 = new TH2F("hMCnormdeclxyd0d0","Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",500,0,10.,200,-0.001,0.001);
		fMCImpParXY = new TH1F("fMCImpParXY","impact parameter XY of B^{+/-} ; l (cm) ; entries",200,-1000.0,1000.0);
		fMCctau = new TH1F("fMCctau","c#tau of B^{+/-} ; c#tau (cm) ; entries",1000,-0.001,1.001);
		hThetaD0 = new TH1F("hThetaD0","D0 from B #theta",3200,0,3.2);
		hThetaPi = new TH1F("hThetaPi","#pi from B #theta",3200,0,3.2);
		hThetaD0Ka = new TH1F("hThetaD0Ka","Ka from D #theta",3200,0,3.2);
		hThetaD0Pi = new TH1F("hThetaD0Pi","#pi from D #theta",3200,0,3.2);
		
		fOutputList->Add(fMCTruth);
		fOutputList->Add(fMCInvMassBpm);
		fOutputList->Add(fMCInvMassDiff);
		fOutputList->Add(fMCCosThetaStar);
		fOutputList->Add(fMCCosPoinAngle);
		fOutputList->Add(fMCCosPoinAngleXY);
		fOutputList->Add(fMCTransvMome);
		fOutputList->Add(fMCDCA);
		fOutputList->Add(fMCTransMomD);
		fOutputList->Add(fMCTransMomPi);
		fOutputList->Add(fMCd0D);
		fOutputList->Add(fMCd0Pi);
		fOutputList->Add(fMCProdd0d0);
		fOutputList->Add(fMCNormDecLthXY);
		fOutputList->Add(hMCdeclgt);
		fOutputList->Add(hMCnormdeclgt);
		fOutputList->Add(hMCdeclgtxy);
		fOutputList->Add(hMCnormdeclgtXY);
		fOutputList->Add(hMCdeclxyd0d0);
		fOutputList->Add(hMCnormdeclxyd0d0);
		fOutputList->Add(fMCImpParXY);
		fOutputList->Add(fMCctau);
		fOutputList->Add(hThetaD0);
		fOutputList->Add(hThetaPi);
		fOutputList->Add(hThetaD0Ka);
		fOutputList->Add(hThetaD0Pi);
	}
	
	
	
	TString nameBplusMass=" ", nameBMassDiff=" ", nameCTS=" ", nameCPA=" ", nameCPAxy=" ", nameBpt=" ", nameBdca=" ", nameDpt=" ", namePiPt=" ", nameDd0=" ", namePid0=" ", named0d0=" ", nameDl2=" ", namenormDl2=" ", nameDlXY=" ", namenormDlXY=" ", nameDlXYd0d0=" ", namenormDlXYd0d0=" ",nameBbXY=" ", nameCtau=" ", nameSDPt=" ", nameLengthBvsD="", nameNormLengthBvsD="", nameMomentumAngle="", nameadvancedCutAll="",nameDnormDl="";
	
	for(Int_t i=0;i<fCuts->GetNPtBins();i++){
		
		// For D
		TString nameDbarMass=" ", nameDmass=" ", nameDCTS=" ", nameDCPA=" ", nameDCPAxy=" ", nameDBpt=" ", nameDBdca=" ", nameDDpt=" ", nameDPiPt=" ", nameDKaP=" ", nameDPiP=" ", nameDDd0=" ", nameDPid0=" ", nameDd0d0=" ", nameDDl2=" ", nameDnormDl2=" ", nameDDlXY=" ", nameDnormDlXY=" ", nameDBbXY=" ", nameDCtau=" ";
		
		nameDbarMass="InvMassDbar_";
		nameDbarMass+=i;
		nameDmass="InvMassD_";
		nameDmass+=i;
		nameDCTS="DCosThetaStar_";
		nameDCTS+=i;
		nameDCPA="DCosPointAngle_";
		nameDCPA+=i;
		nameDCPAxy="DCosPointAngleXY_";
		nameDCPAxy+=i;
		nameDBpt="D0pt_";
		nameDBpt+=i;
		nameDBdca="D0dca_";
		nameDBdca+=i;
		nameDDpt="DKapt_";
		nameDDpt+=i;
		nameDPiPt="DPipt_";
		nameDPiPt+=i;
		nameDKaP="DKap_";
		nameDKaP+=i;
		nameDPiP="DPip_";
		nameDPiP+=i;
		nameDDd0="DKad0_";
		nameDDd0+=i;
		nameDPid0="DPid0_";
		nameDPid0+=i;
		nameDd0d0="DProdd0d0_";
		nameDd0d0+=i;
		nameDDl2="DDecayLength_";
		nameDDl2+=i;
		nameDnormDl2="DnormDecayLength_";
		nameDnormDl2+=i;
		nameDDlXY="DDecayLengthXY_";
		nameDDlXY+=i;
		nameDnormDlXY="DnormDecayLengthXY_";
		nameDnormDlXY+=i;
		nameDBbXY="DBimpParXY_";
		nameDBbXY+=i;
		nameDCtau="Dctau_";
		nameDCtau+=i;
		
		TH1F* tmpD = new TH1F(nameDbarMass.Data(), "m(KPi) - MC; m(PiKPi) (GeV/c^{2}); Entries",(Int_t)(massWindowD*200),1.864-massWindowD,1.864+massWindowD);
		TH1F* tmpDbar = new TH1F(nameDmass.Data(),"m(KPi) - MC; m(PiKPi) (GeV/c^{2}); entries",(Int_t)(massWindowD*200),1.864-massWindowD,1.864+massWindowD);
		
		fOutputListDSandB->Add(tmpD);
		fOutputListDSandB->Add(tmpDbar);
		
		TH1F* tmpDCTS   = new TH1F(nameDCTS.Data(),"D^{+/-} cos(#theta^{*}) w.r.t. Kaon - MC; ; entries",200,-1.0,1.0);
		TH1F* tmpDCPA   = new TH1F(nameDCPA.Data(),"D^{+/-} cos(#theta_{point}) - MC; ; entries",200,-1.0,1.0);
		TH1F* tmpDCPAxy = new TH1F(nameDCPAxy.Data(),"D^{+/-} cos(#theta_{point}) XY - MC; ; entries",200,-1.0,1.0);
		TH1F* tmpDDCA   = new TH1F(nameDBdca.Data(),"D^{+/-} dca - MC; dca (cm) ; entries",100,0.0,10.0);
		TH1F* tmpDDpt   = new TH1F(nameDBpt.Data(),"D^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* tmpKaDpt   = new TH1F(nameDDpt.Data(),"Ka^{+/-} from D^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* tmpDPipt  = new TH1F(nameDPiPt.Data(),"Pi^{+/-} from D^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* tmpKaDp   = new TH1F(nameDKaP.Data(),"Ka^{+/-} from D^{+/-} momentum - MC; p (GeV/c); entries",100,-0.0,100);
		TH1F* tmpDPip  = new TH1F(nameDPiP.Data(),"Pi^{+/-} from D^{+/-} momentum - MC; p (GeV/c); entries",100,-0.0,100);
		TH1F* tmpDDd0   = new TH1F(nameDDd0.Data(),"d0 of Ka^{0}(bar) from D^{+/-} - MC; dca (cm) ; entries",200,-1.001,1.001);
		TH1F* tmpDPid0  = new TH1F(nameDPid0.Data(), "d0 of Pi^{+/-} from D^{+/-} - MC; dca (cm) ; entries",200,-1.001,1.001);
		TH1F* tmpDd0d0  = new TH1F(nameDd0d0.Data(),"d0xd0 of Pi^{+/-} and Ka^{+/-} - MC; (cm^{2}) ; entries",200,-0.001,0.001);
		TH1F* tmpDDbXY  = new TH1F(nameDBbXY.Data(),"impact parameter XY of D - MC; l (cm) ; entries",2000,-100.0,100.0);
		TH1F* tmpDCtau  = new TH1F(nameDCtau.Data(),"c#tau of D - MC; c#tau (cm) ; entries",100,-0.0,1.00);
		TH1F *tmpDdeclgt = new TH1F(nameDDl2.Data(), "Decay Length distribution; Decay Length [cm]",100,0,10);
		TH1F *tmpDnormdeclgt = new TH1F(nameDnormDl2.Data(), "Normalized Decay Length distribution;(Decay Length/Err) ",100,0,100.);
		TH1F *tmpDDlXY = new TH1F(nameDDlXY.Data(),"Decay Length XY distribution;Decay Length XY [cm]",100,0,1.);
		TH1F *tmpDnormDlXY = new TH1F(nameDnormDlXY.Data(),"normalized decay length of D ; l (cm) ; entries",100,0,100);
		
		fOutputListDSandB->Add(tmpDCTS);
		fOutputListDSandB->Add(tmpDCPA);
		fOutputListDSandB->Add(tmpDCPAxy);
		fOutputListDSandB->Add(tmpDDCA);
		fOutputListDSandB->Add(tmpDDpt);
		fOutputListDSandB->Add(tmpKaDpt);
		fOutputListDSandB->Add(tmpDPipt);
		fOutputListDSandB->Add(tmpKaDp);
		fOutputListDSandB->Add(tmpDPip);
		fOutputListDSandB->Add(tmpDDd0);
		fOutputListDSandB->Add(tmpDPid0);
		fOutputListDSandB->Add(tmpDd0d0);
		fOutputListDSandB->Add(tmpDdeclgt);
		fOutputListDSandB->Add(tmpDnormdeclgt);
		fOutputListDSandB->Add(tmpDDlXY);
		fOutputListDSandB->Add(tmpDnormDlXY);
		fOutputListDSandB->Add(tmpDDbXY);
		fOutputListDSandB->Add(tmpDCtau);
		
		//For B
		nameBplusMass="InvMassBplus_";
		nameBplusMass+=i;
		nameBMassDiff="InvBMassDiff_";
		nameBMassDiff+=i;
		nameCTS="CosThetaStar_";
		nameCTS+=i;
		nameCPA="CosPointAngle_";
		nameCPA+=i;
		nameCPAxy="CosPointAngleXY_";
		nameCPAxy+=i;
		nameBpt="Bpt_";
		nameBpt+=i;
		nameBdca="Bdca_";
		nameBdca+=i;
		nameDpt="Dpt_";
		nameDpt+=i;
		namePiPt="Pipt_";
		namePiPt+=i;
		nameDd0="Dd0_";
		nameDd0+=i;
		namePid0="Pid0_";
		namePid0+=i;
		named0d0="Prodd0d0_";
		named0d0+=i;
		nameDl2="DecayLength_";
		nameDl2+=i;
		namenormDl2="normDecayLength_";
		namenormDl2+=i;
		nameDlXY="DecayLengthXY_";
		nameDlXY+=i;
		namenormDlXY="normDecayLengthXY_";
		namenormDlXY+=i;
		nameDlXYd0d0="DLXYd0d0_";
		nameDlXYd0d0+=i;
		namenormDlXYd0d0="normDLXYd0d0_";
		namenormDlXYd0d0+=i;
		nameBbXY="BimpParXY_";
		nameBbXY+=i;
		nameCtau="ctau_";
		nameCtau+=i;
		nameSDPt="countDPt_";
		nameSDPt+=i;
		nameLengthBvsD="LengthBvsD_";
		nameLengthBvsD+=i;
		nameNormLengthBvsD="NormLengthBvsD_";
		nameNormLengthBvsD+=i;
		nameMomentumAngle="MometumAngleDwrtB_";
		nameMomentumAngle+=i;
		nameadvancedCutAll="advancedCutAll_";
		nameadvancedCutAll+=i;
		nameDnormDl="normDDecayLength_";
		nameDnormDl+=i;
		
		TH1F* tmpCTS   = new TH1F(nameCTS.Data(),"B^{+/-} cos(#theta^{*}) w.r.t. Pion; ; entries",200,-1.0,1.0);
		TH1F* tmpCPA   = new TH1F(nameCPA.Data(),"B^{+/-} cos(#theta_{point}); ; entries",200,-1.0,1.0);
		TH1F* tmpCPAxy = new TH1F(nameCPAxy.Data(),"B^{+/-} cos(#theta_{point}) XY; ; entries",200,-1.0,1.0);
		TH1F* tmpDCA	 = new TH1F(nameBdca.Data(),"B^{+/-} dca; dca (cm) ; entries",100,0.0,10.0);
		TH1F* tmpBpt	 = new TH1F(nameBpt.Data(),"B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* tmpDpt   = new TH1F(nameDpt.Data(),"D^{0}(bar) from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* tmpPipt  = new TH1F(namePiPt.Data(),"Pi^{+/-} from B^{+/-} transverse momentum; p_{t} (GeV/c); entries",100,-0.0,100);
		TH1F* tmpDd0   = new TH1F(nameDd0.Data(),"d0 of D^{0}(bar) from B^{+/-}; dca (cm) ; entries",200,-1.001,1.001);
		TH1F* tmpPid0  = new TH1F(namePid0.Data(),"d0 of Pi^{+/-} from B^{+/-} ; dca (cm) ; entries",200,-1.001,1.001);
		TH1F* tmpd0d0  = new TH1F(named0d0.Data(),"d0xd0 of Pi^{+/-} and D^{0}(bar) ; (cm^{2}) ; entries",200,-0.001,.001);
		TH1F *tmpdeclgt = new TH1F(nameDl2.Data(), "Decay Length^{2} distribution; Decay Length^{2} [cm]",100,0,0.1);
		TH1F *tmpnormdeclgt = new TH1F(namenormDl2.Data(), "Normalized Decay Length^{2} distribution;(Decay Length/Err)^{2} ",100,0,100.);
		TH1F *tmpDlXY = new TH1F(nameDlXY.Data(),"Decay Length XY distribution;Decay Length XY [cm]",100,0,1);
		TH1F *tmpnormDlXY = new TH1F(namenormDlXY.Data(),"normalized decay length of B^{+/-} ; l (cm) ; entries",100,0,10);
		TH2F *tmpdeclxyd0d0 = new TH2F(nameDlXYd0d0.Data(),"Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",500,0,1.,200,-0.001,0.001);
		TH2F *tmpnormdeclxyd0d0 = new TH2F(namenormDlXYd0d0.Data(),"Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",500,0,10.,200,-0.001,0.001);
		TH1F* tmpBbXY  = new TH1F(nameBbXY.Data(),"impact parameter XY of B^{+/-} ; l (cm) ; entries",200,-100.0,100.0);
		TH1F* tmpCtau  = new TH1F(nameCtau.Data(),"c#tau of B^{+/-} ; c#tau (cm) ; entries",100,-0.0,1.0);
		TH1F* tmpSDpt  = new TH1F(nameSDPt.Data(),"D^{0}(bar) from B^{+/-} p_{t} in bins from Dp_t; p_{t} (GeV/c); entries",100,-0.0,100);
		TH2F *tmpdeclBvsD = new TH2F(nameLengthBvsD.Data(),"Decay Length B vs D;D Decay Length [cm]; B Decay Length [cm]",100,0.,0.1,100,0.,0.1);
		TH2F *tmpNormdeclBvsD = new TH2F(nameNormLengthBvsD.Data(),"Norm. Decay Length B vs D;Norm. D Decay Length [cm]; Norm. B Decay Length [cm]",1000,0.,10.,1000,0.,10.);
		TH1F* tmpMomAngle   = new TH1F(nameMomentumAngle.Data(),"cos(Angle between p^{D} #cdot p^{B}); ; entries",200,-1.0,1.0);
		TH1F *tmpAdvancedCutAll = new TH1F(nameadvancedCutAll.Data(),"Advanced Candidate Selection - D decay point wrt B decay point and momentum;;entries",2,0,2);
		tmpAdvancedCutAll->GetXaxis()->SetBinLabel(1,"rejected");
		tmpAdvancedCutAll->GetXaxis()->SetBinLabel(2,"passed");
		TH1F *tmpDnormDl = new TH1F(nameDnormDl.Data(),"normalized decay length of #bar{D}^{0} ; l (cm) ; entries",10000,0,100);
		
		TH1F* tmpBm = new TH1F(nameBplusMass.Data(), "m(PiKPi); m(PiKPi) (GeV/c^{2}); Entries",(Int_t)(2*100*massWindowB),5.279-massWindowB,5.279+massWindowB);
		TH1F* tmpBmdiff = new TH1F(nameBMassDiff.Data(),"m(D^{0}#pi_{B})-m(K#pi_{D}); m(D^{0}#pi_{B})-m(K#pi_{D}) (GeV/#it{c}^{2}); entries",(Int_t)(2*100*massWindowB),3.415-massWindowB,3.415+massWindowB);
		
		fOutputListVar->Add(tmpBm);
		fOutputListVar->Add(tmpBmdiff);
		fOutputListVar->Add(tmpCTS);
		fOutputListVar->Add(tmpCPA);
		fOutputListVar->Add(tmpCPAxy);
		fOutputListVar->Add(tmpDCA);
		fOutputListVar->Add(tmpBpt);
		fOutputListVar->Add(tmpDpt);
		fOutputListVar->Add(tmpPipt);
		fOutputListVar->Add(tmpDd0);
		fOutputListVar->Add(tmpPid0);
		fOutputListVar->Add(tmpd0d0);
		fOutputListVar->Add(tmpdeclgt);
		fOutputListVar->Add(tmpnormdeclgt);
		fOutputListVar->Add(tmpDlXY);
		fOutputListVar->Add(tmpnormDlXY);
		fOutputListVar->Add(tmpdeclxyd0d0);
		fOutputListVar->Add(tmpnormdeclxyd0d0);
		fOutputListVar->Add(tmpBbXY);
		fOutputListVar->Add(tmpCtau);
		fOutputListVar->Add(tmpSDpt);
		fOutputListVar->Add(tmpdeclBvsD);
		fOutputListVar->Add(tmpNormdeclBvsD);
		fOutputListVar->Add(tmpMomAngle);
		fOutputListVar->Add(tmpAdvancedCutAll);
		fOutputListVar->Add(tmpDnormDl);
		
		if(fReadMC){
			
			// For B
			TString nameMCBplusMass=" ", nameMCBMassDiff=" ", nameMCCTS=" ", nameMCCPA=" ", nameMCCPAxy=" ", nameMCBpt=" ", nameMCBdca=" ", nameMCDpt=" ", nameMCPiPt=" ", nameMCDd0=" ", nameMCPid0=" ", nameMCd0d0=" ", nameMCDl2=" ", nameMCnormDl2=" ", nameMCDlXY=" ", nameMCnormDlXY=" ", nameMCDlXYd0d0=" ", nameMCnormDlXYd0d0=" ", nameMCBbXY=" ", nameMCCtau=" ",  nameMCLengthBvsD="", nameNormMCLengthBvsD="", nameMCMomentumAngle="", nameadvancedCutMC="", nameDiffPosX="", nameDiffPosY="", nameDiffPosZ="", nameDiffPosDX="", nameDiffPosDY="", nameDiffPosDZ="";
			
			nameMCBplusMass="InvMassBplus_";
			nameMCBplusMass+=i;
			nameMCBMassDiff="InvBMassDiff_";
			nameMCBMassDiff+=i;
			nameMCCTS="CosThetaStar_";
			nameMCCTS+=i;
			nameMCCPA="CosPointAngle_";
			nameMCCPA+=i;
			nameMCCPAxy="CosPointAngleXY_";
			nameMCCPAxy+=i;
			nameMCBpt="Bpt_";
			nameMCBpt+=i;
			nameMCBdca="Bdca_";
			nameMCBdca+=i;
			nameMCDpt="Dpt_";
			nameMCDpt+=i;
			nameMCPiPt="Pipt_";
			nameMCPiPt+=i;
			nameMCDd0="Dd0_";
			nameMCDd0+=i;
			nameMCPid0="Pid0_";
			nameMCPid0+=i;
			nameMCd0d0="Prodd0d0_";
			nameMCd0d0+=i;
			nameMCDl2="DecayLength_";
			nameMCDl2+=i;
			nameMCnormDl2="normDecayLength_";
			nameMCnormDl2+=i;
			nameMCDlXY="DecayLengthXY_";
			nameMCDlXY+=i;
			nameMCnormDlXY="normDecayLengthXY_";
			nameMCnormDlXY+=i;
			nameMCDlXYd0d0="DLXYd0d0_";
			nameMCDlXYd0d0+=i;
			nameMCnormDlXYd0d0="normDLXYd0d0_";
			nameMCnormDlXYd0d0+=i;
			nameMCBbXY="BimpParXY_";
			nameMCBbXY+=i;
			nameMCCtau="ctau_";
			nameMCCtau+=i;
			nameMCLengthBvsD="LengthBvsD_";
			nameMCLengthBvsD+=i;
			nameNormMCLengthBvsD="NormLengthBvsD_";
			nameNormMCLengthBvsD+=i;
			nameMCMomentumAngle="MometumAngleDwrtB_";
			nameMCMomentumAngle+=i;
			nameadvancedCutMC="advancedCutMC_";
			nameadvancedCutMC+=i;
			nameDiffPosX="DiffPosX_";
			nameDiffPosX+=i;
			nameDiffPosY="DiffPosY_";
			nameDiffPosY+=i;
			nameDiffPosZ="DiffPosZ_";
			nameDiffPosZ+=i;
			nameDiffPosDX="DiffPosDX_";
			nameDiffPosDX+=i;
			nameDiffPosDY="DiffPosDY_";
			nameDiffPosDY+=i;
			nameDiffPosDZ="DiffPosDZ_";
			nameDiffPosDZ+=i;
			
			TH1F* tmpMCCTS   = new TH1F(nameMCCTS.Data(),"B^{+/-} cos(#theta^{*}) w.r.t. Pion - MC; ; entries",200,-1.0,1.0);
			TH1F* tmpMCCPA   = new TH1F(nameMCCPA.Data(),"B^{+/-} cos(#theta_{point}) - MC; ; entries",200,-1.0,1.0);
			TH1F* tmpMCCPAxy = new TH1F(nameMCCPAxy.Data(),"B^{+/-} cos(#theta_{point}) XY - MC; ; entries",200,-1.0,1.0);
			TH1F* tmpMCDCA   = new TH1F(nameMCBdca.Data(),"B^{+/-} dca - MC; dca (cm) ; entries",100,0.0,10.0);
			TH1F* tmpMCBpt   = new TH1F(nameMCBpt.Data(),"B^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCDpt   = new TH1F(nameMCDpt.Data(),"D^{0}(bar) from B^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCPipt  = new TH1F(nameMCPiPt.Data(),"Pi^{+/-} from B^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCDd0   = new TH1F(nameMCDd0.Data(),"d0 of D^{0}(bar) from B^{+/-} - MC; dca (cm) ; entries",200,-1.001,1.001);
			TH1F* tmpMCPid0  = new TH1F(nameMCPid0.Data(), "d0 of Pi^{+/-} from B^{+/-} - MC; dca (cm) ; entries",200,-1.001,1.001);
			TH1F* tmpMCd0d0  = new TH1F(nameMCd0d0.Data(),"d0xd0 of Pi^{+/-} and D^{0}(bar) - MC; (cm^{2}) ; entries",200,-0.001,0.001);
			TH1F* tmpMCBbXY  = new TH1F(nameMCBbXY.Data(),"impact parameter XY of B^{+/-} - MC; l (cm) ; entries",2000,-100.0,100.0);
			TH1F* tmpMCCtau  = new TH1F(nameMCCtau.Data(),"c#tau of B^{+/-} - MC; c#tau (cm) ; entries",100,-0.0,1.00);
			TH1F *tmpMCdeclgt = new TH1F(nameMCDl2.Data(), "Decay Length^{2} distribution; Decay Length^{2} [cm]",100,0,0.1);
			TH1F *tmpMCnormdeclgt = new TH1F(nameMCnormDl2.Data(), "Normalized Decay Length^{2} distribution;(Decay Length/Err)^{2} ",100,0,100.);
			TH1F *tmpMCDlXY = new TH1F(nameMCDlXY.Data(),"Decay Length XY distribution;Decay Length XY [cm]",100,0,1.);
			TH1F *tmpMCnormDlXY = new TH1F(nameMCnormDlXY.Data(),"normalized decay length of B^{+/-} ; l (cm) ; entries",100,0,10);
			TH2F *tmpMCdeclxyd0d0 = new TH2F(nameMCDlXYd0d0.Data(),"Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",500,0,1.,200,-0.001,0.001);
			TH2F *tmpMCnormdeclxyd0d0 = new TH2F(nameMCnormDlXYd0d0.Data(),"Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",500,0,10.,200,-0.001,0.001);
			TH2F *tmpMCdeclBvsD = new TH2F(nameMCLengthBvsD.Data(),"Decay Length B vs D;D Decay Length [cm];B Decay Length [cm]",100,0.,0.1,100,0.,0.1);
			TH2F *tmpMCNormdeclBvsD = new TH2F(nameNormMCLengthBvsD.Data(),"Norm. Decay Length B vs D;Norm. D Decay Length [cm];Norm. B Decay Length [cm]",1000,0.,10.,1000,0.,10.);
			TH1F* tmpMCMomAngle   = new TH1F(nameMCMomentumAngle.Data(),"cos(Angle between p^{D} #cdot p^{B}) - MC; ; entries",200,-1.0,1.0);
			TH1F *tmpAdvancedCutMC = new TH1F(nameadvancedCutMC.Data(),"MC Advanced Candidate Selection - D decay point wrt B decay point and momentum;;entries",2,0,2);
			tmpAdvancedCutMC->GetXaxis()->SetBinLabel(1,"rejected");
			tmpAdvancedCutMC->GetXaxis()->SetBinLabel(2,"passed");
			TH1F *tmpDiffPosX = new TH1F(nameDiffPosX.Data(),"X (Reconstructed-True) Vertex Position;cm;entries",100,-0.05,0.05);
			TH1F *tmpDiffPosY = new TH1F(nameDiffPosY.Data(),"Y (Reconstructed-True) Vertex Position;cm;entries",100,-0.05,0.05);
			TH1F *tmpDiffPosZ = new TH1F(nameDiffPosZ.Data(),"Z (Reconstructed-True) Vertex Position;cm;entries",100,-0.05,0.05);
			TH1F *tmpDiffPosDX = new TH1F(nameDiffPosDX.Data(),"D X (Reconstructed-True) Vertex Position;cm;entries",100,-0.05,0.05);
			TH1F *tmpDiffPosDY = new TH1F(nameDiffPosDY.Data(),"D Y (Reconstructed-True) Vertex Position;cm;entries",100,-0.05,0.05);
			TH1F *tmpDiffPosDZ = new TH1F(nameDiffPosDZ.Data(),"D Z (Reconstructed-True) Vertex Position;cm;entries",100,-0.05,0.05);
			
			TH1F* tmpMCBm = new TH1F(nameMCBplusMass.Data(), "m(PiKPi) - MC; m(PiKPi) (GeV/c^{2}); Entries",(Int_t)(2*100*massWindowB),5.279-massWindowB,5.279+massWindowB);
			TH1F* tmpMCBmdiff = new TH1F(nameMCBMassDiff.Data(),"m(D^{0}#pi_{B})-m(K#pi_{D}) - MC; m(D^{0}#pi_{B})-m(K#pi_{D}) (GeV/#it{c}^{2}); entries",(Int_t)(2*100*massWindowB),3.415-massWindowB,3.415+massWindowB);
			
			fOutputListMCVar->Add(tmpMCBm);
			fOutputListMCVar->Add(tmpMCBmdiff);
			
			fOutputListMCVar->Add(tmpMCCTS);
			fOutputListMCVar->Add(tmpMCCPA);
			fOutputListMCVar->Add(tmpMCCPAxy);
			fOutputListMCVar->Add(tmpMCDCA);
			fOutputListMCVar->Add(tmpMCBpt);
			fOutputListMCVar->Add(tmpMCDpt);
			fOutputListMCVar->Add(tmpMCPipt);
			fOutputListMCVar->Add(tmpMCDd0);
			fOutputListMCVar->Add(tmpMCPid0);
			fOutputListMCVar->Add(tmpMCd0d0);
			fOutputListMCVar->Add(tmpMCdeclgt);
			fOutputListMCVar->Add(tmpMCnormdeclgt);
			fOutputListMCVar->Add(tmpMCDlXY);
			fOutputListMCVar->Add(tmpMCnormDlXY);
			fOutputListMCVar->Add(tmpMCdeclxyd0d0);
			fOutputListMCVar->Add(tmpMCnormdeclxyd0d0);
			fOutputListMCVar->Add(tmpMCBbXY);
			fOutputListMCVar->Add(tmpMCCtau);
			fOutputListMCVar->Add(tmpMCdeclBvsD);
			fOutputListMCVar->Add(tmpMCNormdeclBvsD);
			fOutputListMCVar->Add(tmpMCMomAngle);
			fOutputListMCVar->Add(tmpAdvancedCutMC);
			fOutputListMCVar->Add(tmpDiffPosX);
			fOutputListMCVar->Add(tmpDiffPosY);
			fOutputListMCVar->Add(tmpDiffPosZ);
			fOutputListMCVar->Add(tmpDiffPosDX);
			fOutputListMCVar->Add(tmpDiffPosDY);
			fOutputListMCVar->Add(tmpDiffPosDZ);
			
			// For D
			TString nameMCDbarMass=" ", nameMCDmass=" ", nameMCDCTS=" ", nameMCDCPA=" ", nameMCDCPAxy=" ", nameMCDBpt=" ", nameMCDBdca=" ", nameMCDDpt=" ", nameMCDPiPt=" ", nameMCDKaP=" ", nameMCDPiP=" ", nameMCDDd0=" ", nameMCDPid0=" ", nameMCDd0d0=" ", nameMCDDl2=" ", nameMCDnormDl2=" ", nameMCDDlXY=" ", nameMCDnormDlXY=" ", nameMCDBbXY=" ", nameMCDCtau=" ", nameDRecoPosX="",nameDRecoPosY="",nameDRecoPosZ="";//, nameDTruePosX="", nameDTruePosY="", nameDTruePosZ="", nameDDiffPosX="", nameDDiffPosY="", nameDDiffPosZ="";
			
			nameMCDbarMass="InvMassDbar_";
			nameMCDbarMass+=i;
			nameMCDmass="InvMassD_";
			nameMCDmass+=i;
			nameMCDCTS="DCosThetaStar_";
			nameMCDCTS+=i;
			nameMCDCPA="DCosPointAngle_";
			nameMCDCPA+=i;
			nameMCDCPAxy="DCosPointAngleXY_";
			nameMCDCPAxy+=i;
			nameMCDBpt="D0pt_";
			nameMCDBpt+=i;
			nameMCDBdca="D0dca_";
			nameMCDBdca+=i;
			nameMCDDpt="DKapt_";
			nameMCDDpt+=i;
			nameMCDPiPt="DPipt_";
			nameMCDPiPt+=i;
			nameMCDKaP="DKap_";
			nameMCDKaP+=i;
			nameMCDPiP="DPip_";
			nameMCDPiP+=i;
			nameMCDDd0="DKad0_";
			nameMCDDd0+=i;
			nameMCDPid0="DPid0_";
			nameMCDPid0+=i;
			nameMCDd0d0="DProdd0d0_";
			nameMCDd0d0+=i;
			nameMCDDl2="DDecayLength_";
			nameMCDDl2+=i;
			nameMCDnormDl2="DnormDecayLength_";
			nameMCDnormDl2+=i;
			nameMCDDlXY="DDecayLengthXY_";
			nameMCDDlXY+=i;
			nameMCDnormDlXY="DnormDecayLengthXY_";
			nameMCDnormDlXY+=i;
			nameMCDBbXY="DBimpParXY_";
			nameMCDBbXY+=i;
			nameMCDCtau="Dctau_";
			nameMCDCtau+=i;
			
			TH1F* tmpMCD = new TH1F(nameMCDbarMass.Data(), "m(KPi) - MC; m(PiKPi) (GeV/c^{2}); Entries",120,1.864-massWindowD,1.864+massWindowD);
			TH1F* tmpMCDbar = new TH1F(nameMCDmass.Data(),"m(KPi) - MC; m(PiKPi) (GeV/c^{2}); entries",120,1.864-massWindowD,1.864+massWindowD);
			
			fOutputListDSonly->Add(tmpMCD);
			fOutputListDSonly->Add(tmpMCDbar);
			
			TH1F* tmpMCDCTS   = new TH1F(nameMCDCTS.Data(),"D^{+/-} cos(#theta^{*}) w.r.t. Kaon - MC; ; entries",200,-1.0,1.0);
			TH1F* tmpMCDCPA   = new TH1F(nameMCDCPA.Data(),"D^{+/-} cos(#theta_{point}) - MC; ; entries",200,-1.0,1.0);
			TH1F* tmpMCDCPAxy = new TH1F(nameMCDCPAxy.Data(),"D^{+/-} cos(#theta_{point}) XY - MC; ; entries",200,-1.0,1.0);
			TH1F* tmpMCDDCA   = new TH1F(nameMCDBdca.Data(),"D^{+/-} dca - MC; dca (cm) ; entries",100,0.0,10.0);
			TH1F* tmpMCDDpt   = new TH1F(nameMCDBpt.Data(),"D^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCKaDpt   = new TH1F(nameMCDDpt.Data(),"Ka^{+/-} from D^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCDPipt  = new TH1F(nameMCDPiPt.Data(),"Pi^{+/-} from D^{+/-} transverse momentum - MC; p_{t} (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCKaDp   = new TH1F(nameMCDKaP.Data(),"Ka^{+/-} from D^{+/-} momentum - MC; p (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCDPip  = new TH1F(nameMCDPiP.Data(),"Pi^{+/-} from D^{+/-} momentum - MC; p (GeV/c); entries",100,-0.0,100);
			TH1F* tmpMCDDd0   = new TH1F(nameMCDDd0.Data(),"d0 of Ka^{0}(bar) from D^{+/-} - MC; dca (cm) ; entries",200,-1.001,1.001);
			TH1F* tmpMCDPid0  = new TH1F(nameMCDPid0.Data(), "d0 of Pi^{+/-} from D^{+/-} - MC; dca (cm) ; entries",200,-1.001,1.001);
			TH1F* tmpMCDd0d0  = new TH1F(nameMCDd0d0.Data(),"d0xd0 of Pi^{+/-} and Ka^{+/-} - MC; (cm^{2}) ; entries",200,-0.001,0.001);
			TH1F* tmpMCDDbXY  = new TH1F(nameMCDBbXY.Data(),"impact parameter XY of D - MC; l (cm) ; entries",2000,-100.0,100.0);
			TH1F* tmpMCDCtau  = new TH1F(nameMCDCtau.Data(),"c#tau of D - MC; c#tau (cm) ; entries",100,-0.0,1.00);
			TH1F *tmpMCDdeclgt = new TH1F(nameMCDDl2.Data(), "Decay Length distribution; Decay Length [cm]",100,0,10);
			TH1F *tmpMCDnormdeclgt = new TH1F(nameMCDnormDl2.Data(), "Normalized Decay Length distribution;(Decay Length/Err) ",100,0,100.);
			TH1F *tmpMCDDlXY = new TH1F(nameMCDDlXY.Data(),"Decay Length XY distribution;Decay Length XY [cm]",100,0,1.);
			TH1F *tmpMCDnormDlXY = new TH1F(nameMCDnormDlXY.Data(),"normalized decay length of D ; l (cm) ; entries",100,0,100);
			
			fOutputListDSonly->Add(tmpMCDCTS);
			fOutputListDSonly->Add(tmpMCDCPA);
			fOutputListDSonly->Add(tmpMCDCPAxy);
			fOutputListDSonly->Add(tmpMCDDCA);
			fOutputListDSonly->Add(tmpMCDDpt);
			fOutputListDSonly->Add(tmpMCKaDpt);
			fOutputListDSonly->Add(tmpMCDPipt);
			fOutputListDSonly->Add(tmpMCKaDp);
			fOutputListDSonly->Add(tmpMCDPip);
			fOutputListDSonly->Add(tmpMCDDd0);
			fOutputListDSonly->Add(tmpMCDPid0);
			fOutputListDSonly->Add(tmpMCDd0d0);
			fOutputListDSonly->Add(tmpMCDdeclgt);
			fOutputListDSonly->Add(tmpMCDnormdeclgt);
			fOutputListDSonly->Add(tmpMCDDlXY);
			fOutputListDSonly->Add(tmpMCDnormDlXY);
			fOutputListDSonly->Add(tmpMCDDbXY);
			fOutputListDSonly->Add(tmpMCDCtau);
		}
	}
	if(fbtrackRotation==kTRUE){
		// for rotated tracks
		TH1F *hRotationMonitor = new TH1F("hRotationMonitor", "hRotationMonitor",10,0,10);
		hRotationMonitor->GetXaxis()->SetBinLabel(1,"Rotations Performed");
		hRotationMonitor->GetXaxis()->SetBinLabel(2,"Rotations Inside DeltaMass");
		hRotationMonitor->GetXaxis()->SetBinLabel(3,"Tries to reconstruct vertex");
		hRotationMonitor->GetXaxis()->SetBinLabel(4,"Fails  to reconstruct vertex");
		hRotationMonitor->GetXaxis()->SetBinLabel(5,"No primary Vertex");
		hRotationMonitor->GetXaxis()->SetBinLabel(6,"Sucessfull rotations");
		hRotationMonitor->GetXaxis()->SetBinLabel(7,"Rejected by 'Kill Candidates'");
		
		fOutputList->Add(hRotationMonitor);
		
	}
	fNentries = new TH1F("HistoforChecks", "HistoforChecks",140,0,140);
	
	fNentries->GetXaxis()->SetBinLabel(1,"Events Analyzed");
	if(fReadMC){
		fNentries->GetXaxis()->SetBinLabel(2,"D(B+) -> 0. - 0.5");
		fNentries->GetXaxis()->SetBinLabel(3,"D(B+) -> 0.5 - 1.");
		fNentries->GetXaxis()->SetBinLabel(4,"D(B+)-> 1. - 2.");
		fNentries->GetXaxis()->SetBinLabel(5,"D(B+)-> 2. - 3.");
		fNentries->GetXaxis()->SetBinLabel(6,"D(B+)-> 3. - 4.");
		fNentries->GetXaxis()->SetBinLabel(7,"D(B+)-> 4. - 5.");
		fNentries->GetXaxis()->SetBinLabel(8,"D(B+)-> 5. - 6.");
		fNentries->GetXaxis()->SetBinLabel(9,"D(B+)-> 6. - 8.");
		fNentries->GetXaxis()->SetBinLabel(10,"D(B+)-> 8. - 12.");
		fNentries->GetXaxis()->SetBinLabel(11,"D(B+)-> 12. - 16.");
		fNentries->GetXaxis()->SetBinLabel(12,"D(B+)-> 16. - 20.");
		fNentries->GetXaxis()->SetBinLabel(13,"D(B+)-> 20. - 24.");
		fNentries->GetXaxis()->SetBinLabel(14,"D(B+)-> 24. - inf");
	}
	else{
		fNentries->GetXaxis()->SetBinLabel(2,"-");
		fNentries->GetXaxis()->SetBinLabel(3,"-");
		fNentries->GetXaxis()->SetBinLabel(4,"-");
		fNentries->GetXaxis()->SetBinLabel(5,"-");
		fNentries->GetXaxis()->SetBinLabel(6,"-");
		fNentries->GetXaxis()->SetBinLabel(7,"-");
		fNentries->GetXaxis()->SetBinLabel(8,"-");
		fNentries->GetXaxis()->SetBinLabel(9,"-");
		fNentries->GetXaxis()->SetBinLabel(10,"-");
		fNentries->GetXaxis()->SetBinLabel(11,"-");
		fNentries->GetXaxis()->SetBinLabel(12,"-");
		fNentries->GetXaxis()->SetBinLabel(13,"-");
		fNentries->GetXaxis()->SetBinLabel(14,"-");
	}
	fNentries->GetXaxis()->SetBinLabel(15,"-");
	fNentries->GetXaxis()->SetBinLabel(16,"D0=1");
	fNentries->GetXaxis()->SetBinLabel(17,"D0bar=2");
	fNentries->GetXaxis()->SetBinLabel(18,"either=3");
	fNentries->GetXaxis()->SetBinLabel(19,"-");
	fNentries->GetXaxis()->SetBinLabel(20,"FillBpmHistos");
	fNentries->GetXaxis()->SetBinLabel(21,"D0 in #Deltam_{B}");
	fNentries->GetXaxis()->SetBinLabel(22,"D0bar in #Deltam_{B}");
	fNentries->GetXaxis()->SetBinLabel(23,"e D0 in #Deltam_{B}");
	fNentries->GetXaxis()->SetBinLabel(24,"e D0bar in #Deltam_{B}");
	fNentries->GetXaxis()->SetBinLabel(25,"Finished Bpm hists");
	fNentries->GetXaxis()->SetBinLabel(26,"-");
	if(fReadMC){
		fNentries->GetXaxis()->SetBinLabel(27,"D(B-)-> 0. - 0.5");
		fNentries->GetXaxis()->SetBinLabel(28,"D(B-)-> 0.5 - 1.");
		fNentries->GetXaxis()->SetBinLabel(29,"D(B-)-> 1. - 2.");
		fNentries->GetXaxis()->SetBinLabel(30,"D(B-)-> 2. - 3.");
		fNentries->GetXaxis()->SetBinLabel(31,"D(B-)-> 3. - 4.");
		fNentries->GetXaxis()->SetBinLabel(32,"D(B-)-> 4. - 5.");
		fNentries->GetXaxis()->SetBinLabel(33,"D(B-)-> 5. - 6.");
		fNentries->GetXaxis()->SetBinLabel(34,"D(B-)-> 6. - 8.");
		fNentries->GetXaxis()->SetBinLabel(35,"D(B-)-> 8. - 12.");
		fNentries->GetXaxis()->SetBinLabel(36,"D(B-)-> 12. - 16.");
		fNentries->GetXaxis()->SetBinLabel(37,"D(B-)-> 16. - 20.");
		fNentries->GetXaxis()->SetBinLabel(38,"D(B-)-> 20. - 24.");
		fNentries->GetXaxis()->SetBinLabel(39,"D(B-)-> 24. - inf");
	}
	else{
		fNentries->GetXaxis()->SetBinLabel(27,"-");
		fNentries->GetXaxis()->SetBinLabel(28,"-");
		fNentries->GetXaxis()->SetBinLabel(29,"-");
		fNentries->GetXaxis()->SetBinLabel(30,"-");
		fNentries->GetXaxis()->SetBinLabel(31,"-");
		fNentries->GetXaxis()->SetBinLabel(32,"-");
		fNentries->GetXaxis()->SetBinLabel(33,"-");
		fNentries->GetXaxis()->SetBinLabel(34,"-");
		fNentries->GetXaxis()->SetBinLabel(35,"-");
		fNentries->GetXaxis()->SetBinLabel(36,"-");
		fNentries->GetXaxis()->SetBinLabel(37,"-");
		fNentries->GetXaxis()->SetBinLabel(38,"-");
		fNentries->GetXaxis()->SetBinLabel(39,"-");
	}
	if(fReadMC){
		fNentries->GetXaxis()->SetBinLabel(40,"B from b");
		fNentries->GetXaxis()->SetBinLabel(41,"B injected");
		
		fNentries->GetXaxis()->SetBinLabel(42,"B(b): D0");
		fNentries->GetXaxis()->SetBinLabel(43,"B(b,D0):mum !521");
		fNentries->GetXaxis()->SetBinLabel(44,"B(b,D0):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(45,"B(b,D0):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(46,"B(b,D0):!211");
		fNentries->GetXaxis()->SetBinLabel(47,"B(b,D0):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(48,"B(b,D0):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(49,"B(b): D0bar");
		fNentries->GetXaxis()->SetBinLabel(50,"B(b,D0-):mum !521");
		fNentries->GetXaxis()->SetBinLabel(51,"B(b,D0-):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(52,"B(b,D0-):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(53,"B(b,D0-):!211");
		fNentries->GetXaxis()->SetBinLabel(54,"B(b,D0-):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(55,"B(b,D0-):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(56,"B(b): both (D0)");
		fNentries->GetXaxis()->SetBinLabel(57,"B(b,D0b):mum !521");
		fNentries->GetXaxis()->SetBinLabel(58,"B(b,D0b):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(59,"B(b,D0b):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(60,"B(b,D0b):!211");
		fNentries->GetXaxis()->SetBinLabel(61,"B(b,D0b):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(62,"B(b,D0b):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(63,"B(b): both (D0bar)");
		fNentries->GetXaxis()->SetBinLabel(64,"B(b,D0-b):mum !521");
		fNentries->GetXaxis()->SetBinLabel(65,"B(b,D0-b):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(66,"B(b,D0-b):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(67,"B(b,D0-b):!211");
		fNentries->GetXaxis()->SetBinLabel(68,"B(b,D0-b):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(69,"B(b,D0-b):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(70,"B(b,D0): no mum match");
		fNentries->GetXaxis()->SetBinLabel(71,"B(b,D0-): no mum match");
		fNentries->GetXaxis()->SetBinLabel(72,"B(b,D0b): no mum match");
		fNentries->GetXaxis()->SetBinLabel(73,"B(b,D0-b): no mum match");
		
		fNentries->GetXaxis()->SetBinLabel(74,"B(inj): D0");
		fNentries->GetXaxis()->SetBinLabel(75,"B(inj,D0):mum !521");
		fNentries->GetXaxis()->SetBinLabel(76,"B(inj,D0):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(77,"B(inj,D0):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(78,"B(inj,D0):!211");
		fNentries->GetXaxis()->SetBinLabel(79,"B(inj,D0):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(80,"B(inj,D0):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(81,"B(inj): D0bar");
		fNentries->GetXaxis()->SetBinLabel(82,"B(inj,D0-):mum !521");
		fNentries->GetXaxis()->SetBinLabel(83,"B(inj,D0-):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(84,"B(inj,D0-):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(85,"B(inj,D0-):!211");
		fNentries->GetXaxis()->SetBinLabel(86,"B(inj,D0-):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(87,"B(inj,D0-):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(88,"B(inj): both (D0)");
		fNentries->GetXaxis()->SetBinLabel(89,"B(inj,D0b):mum !521");
		fNentries->GetXaxis()->SetBinLabel(90,"B(inj,D0b):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(91,"B(inj,D0b):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(92,"B(inj,D0b):!211");
		fNentries->GetXaxis()->SetBinLabel(93,"B(inj,D0b):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(94,"B(inj,D0b):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(95,"B(inj): both (D0bar)");
		fNentries->GetXaxis()->SetBinLabel(96,"B(inj,D0-b):mum !521");
		fNentries->GetXaxis()->SetBinLabel(97,"B(inj,D0-b):!2dghts");
		fNentries->GetXaxis()->SetBinLabel(98,"B(inj,D0-b):noMCpart");
		fNentries->GetXaxis()->SetBinLabel(99,"B(inj,D0-b):!211");
		fNentries->GetXaxis()->SetBinLabel(100,"B(inj,D0-b):no 211-mum");
		fNentries->GetXaxis()->SetBinLabel(101,"B(inj,D0-b):mum !521");
		
		fNentries->GetXaxis()->SetBinLabel(102,"B(inj,D0): no mum match");
		fNentries->GetXaxis()->SetBinLabel(103,"B(inj,D0-): no mum match");
		fNentries->GetXaxis()->SetBinLabel(104,"B(inj,D0b): no mum match");
		fNentries->GetXaxis()->SetBinLabel(105,"B(inj,D0-b): no mum match");
		
	}
	else{
		fNentries->GetXaxis()->SetBinLabel(40,"-");
		fNentries->GetXaxis()->SetBinLabel(41,"-");
	}
	fNentries->GetXaxis()->SetBinLabel(106,"Kill Candidate Cuts?");
	fNentries->GetXaxis()->SetBinLabel(107,"B: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(108,"B: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(109,"B: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(110,"MC_B: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(111,"MC_B: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(112,"MC_B: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(113,"D: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(114,"D: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(115,"D: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(116,"MC_D: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(117,"MC_D: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(118,"MC_D: D not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(119,"D: B not in #Delta mB");
	fNentries->GetXaxis()->SetBinLabel(120,"MC_D: B not in #Delta mB");
	fNentries->GetXaxis()->SetBinLabel(121,"D is Hijing");
	fNentries->GetXaxis()->SetBinLabel(122,"D is not Hijing");
	fNentries->GetXaxis()->SetBinLabel(123,"Pi is Hijing");
	fNentries->GetXaxis()->SetBinLabel(124,"Pi is not Hijing");
	fNentries->GetXaxis()->SetBinLabel(125,"rejected by advanced cut");
	fNentries->GetXaxis()->SetBinLabel(126,"passed advanced cut");
	fNentries->GetXaxis()->SetBinLabel(127,"MC rejected by advanced cut");
	fNentries->GetXaxis()->SetBinLabel(128,"MC passed advanced cut");
	fNentries->GetXaxis()->SetBinLabel(129,"Tries to reconstruct B vertex");
	fNentries->GetXaxis()->SetBinLabel(130,"Fails to reconstruct B vertex");
	fNentries->GetXaxis()->SetBinLabel(131,"D0 Not Selected");
	fNentries->GetXaxis()->SetBinLabel(132,"D0 not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(133,"D0bar not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(134,"Either not in #Delta mD");
	fNentries->GetXaxis()->SetBinLabel(135,"Tries to select D0");
	fNentries->GetXaxis()->SetBinLabel(136,"Success to select D0");
	
	fNentries->GetXaxis()->SetNdivisions(1,kFALSE);
	
	fOutputList->Add(fNentries);
	
	// Post the data
	PostData(1,fOutputList);
	PostData(2,fOutputListVar);
	PostData(3,fOutputListMCVar);
	PostData(4,fOutputListDSandB);
	PostData(5,fOutputListDSonly);
	PostData(6,fSelectionVariables);
	PostData(7,fRotatedCandidates);
	return;
}
//________________________________________________________________________
void AliAnalysisTaskSEBpmMass::UserExec(Option_t */*option*/){
	//
	// Execute analysis for current event:
	// Full reconstruction of B Mesons with and without MC verification
	//
	if(fDebug>1) cout<<"I'm in UserExec"<<endl;
	
	// array for optional use of HIJING flag in FillBpmHists()
	
	AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
	
	if(fReadMC){
		TClonesArray *mc_check=static_cast<TClonesArray*>(aod->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
		if(!mc_check)return;
	}
	
	TString bname="D0toKpi";
	
	fBzkG = aod->GetMagneticField();
	
	TClonesArray *inputArray=0;
	
	if(!aod && AODEvent() && IsStandardAOD()) {
		// In case there is an AOD handler writing a standard AOD, use the AOD
		// event in memory rather than the input (ESD) event.
		aod = dynamic_cast<AliAODEvent*> (AODEvent());
		// in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
		// have to taken from the AOD event hold by the AliAODExtension
		AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
		
		if(aodHandler->GetExtensions()) {
			AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
			AliAODEvent* aodFromExt = ext->GetAOD();
			inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
		}
	} else if(aod) {
		inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
	}
	
	if(!inputArray || !aod) {
		printf("AliAnalysisTaskSEBpmMass::UserExec: input branch not found!\n");
		return;
	}
	
	if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;
	
	// Centrality Selection
	if(fbCentrality==kTRUE){
		// Do not use for LHC12c4 because V0 not simulated
		if(aod->GetCentrality()->GetCentralityPercentile("V0M")<0 || aod->GetCentrality()->GetCentralityPercentile("V0M")>10) return;
	}
	
	if(fReadMC) {
		// load MC particles
		fMcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
		if(!fMcArray) {
			printf("AliAnalysisTaskSEBpmMass::UserExec: MC particles branch not found!\n");
			return;
		}
		// load MC header
		fmcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		if(!fmcHeader) {
			printf("AliAnalysisTaskSEBpmMass::UserExec: MC header branch not found!\n");
			return;
		}
		AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(fmcHeader);
		if(fbSelectHIJING && !fbHIJINGavailable){
			if(fDebug>1)printf("#################### You want HIJING but it is not available #################### \n");
			return;
		}
		if(fbHIJINGavailable && !fReadMC){
			if(fDebug>1)printf("#################### You want HIJING but it is not event MC ##################### \n");
			return;
		}
		if(fbHIJINGavailable){
			cout << "####################" << endl;
			cout << "# Hijing Available #" << endl;
			if(fbSelectHIJING) cout << "# Hijing  Selected #" << endl;
			cout << "####################" << endl;
			
			for (Int_t i=0; i<(Int_t)fmcHeader->GetNCocktailHeaders(); i++) {
				hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(fmcHeader->GetCocktailHeader(i));
				if (hijingGenHeader){
					fhijingGenHeader=hijingGenHeader;
					break;
				}
			}
			if(!hijingGenHeader){
				cout << "Hijing event header not found" << endl;
				return;
			}
		}
	}
	
	if(fbtrackRotation){
		cout << "#########################" << endl;
		cout << "# Track Rotation Active #" << endl;
		cout << "#########################" << endl;
	}
	if(fKillCandidates){
		cout << "########################" << endl;
		cout << "# Cut Selection Active #" << endl;
		cout << "########################" << endl;
	}
	//histogram filled with 1 for every AOD
	fNentries->Fill(0);
	
	// AOD primary vertex
	fvtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
	
	
	fCuts->IsEventSelected(aod);
	
	// loop over D0(bar)candidates
	Int_t nInD0toKpi = inputArray->GetEntriesFast();
//	if(fDebug>2)
		printf("Number of D0->Kpi: %d\n",nInD0toKpi);
	
	for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
		AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);
		
		//______________________________________________________________________________________
		// Select suitable D0 candidates
		// use this cut to select mass hypothesis within mass range [and cos(ts) hypothesis] before vertex reconstruction
		if(fDebug>2) cout << "Applying D0 Mass cuts" << endl;
		
		// 0 = rejected		--> Get next candidate
		// 1 = D0			--> Check if D0 is in mass range
		// 2 = D0bar		--> Check if D0bar is in mass range
		// 3 = both			--> Check if either is in mass range
		fNentries->Fill(134);
		Int_t selectedD0cand = fCuts->IsSelected(d,AliRDHFCuts::kAll,aod); // alternatively use AliRDHFCuts::kCandidate  AliRDHFCuts::kAll kPID  kTracks
		fNentries->Fill(135);
		// Set primary vertex if not available by default
		if(!d->GetOwnPrimaryVtx()) {
			d->SetOwnPrimaryVtx(fvtx1); // needed to compute all variables
		}
		
		Double_t massTrueD0 = 1.8648;//TDatabasePDG::Instance()->GetParticle(421)->Mass();
		
		if(selectedD0cand==0){
			fNentries->Fill(130);
			continue;
		}
		else if(selectedD0cand==1){
			if(TMath::Abs(d->InvMassD0() - massTrueD0) > massWindowD){
				fNentries->Fill(131);
				continue;
			}
			fNentries->Fill(15);
		}
		else if(selectedD0cand==2){
			if(TMath::Abs(d->InvMassD0bar() - massTrueD0) > massWindowD){
				fNentries->Fill(132);
				continue;
			}
			fNentries->Fill(16);
		}
		else if(selectedD0cand==3){
			if((TMath::Abs(d->InvMassD0() - massTrueD0) > massWindowD) && (TMath::Abs(d->InvMassD0bar() - massTrueD0) > massWindowD)){
				fNentries->Fill(133);
				continue;
			}
			fNentries->Fill(17);
		}
		
		// If requested and available, select only Hijing particles (or not)
		if(fbHIJINGavailable){
			Bool_t bInjectedD = IsCandidateInjected(d,fmcHeader);
			if(fbSelectHIJING && bInjectedD==kTRUE) {
				//if(fDebug>1) cout << "D candidate rejected: Not from Hijing (and we want HIJING)" << endl;
				fNentries->Fill(121);
				continue;
			}
			if(fbSelectHIJING && bInjectedD==kFALSE){
				//if(fDebug>1) cout << "D candidate from HIJING (and we want HIJING)" << endl;
				fNentries->Fill(120);
			}
			if(!fbSelectHIJING && bInjectedD==kFALSE) {
				//if(fDebug>1) cout << "D candidate rejected: It is from Hijing (and we don't want HIJING)" << endl;
				fNentries->Fill(120);
				continue;
			}
			
			if(!fbSelectHIJING && bInjectedD==kTRUE){
				//if(fDebug>1) cout << "D candidate not from HIJING (and we don't want HIJING)" << endl;
				fNentries->Fill(121);
			}
		}
		
		
		// check that there are only two daughters
		if(fDebug>2) cout << "D candidates has how many daughters? " << d->GetNDaughters() << endl;
		
		if(d->GetNDaughters()!=2 || d->GetNProngs()!=2) {
			continue;
		}
		// check the track label, so it cannot be attached later to the D
		Int_t Ddaughter1 = d->GetProngID(0);
		Int_t Ddaughter2 = d->GetProngID(1);
		if(fDebug>2){
			cout << "D daughter labels are: " << Ddaughter1 << " and " << Ddaughter2 << endl;
			cout << "selectedD0cand " << selectedD0cand << endl;
			cout << "After D0 Mass cuts" << endl;
		}
		// Centrality Selection
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(d->Pt());
		fillthis="hcentrality";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(aod/*->GetHeader()*/->GetCentrality()->GetCentralityPercentile("V0M"));
		fillthis="hSDpt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(d->Pt());
		fillthis="countDPt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->Pt());
		if(fDebug>2) cout << "After centrality selection" << endl;
		
		///if(!fReadMC)  fNentries->Fill(1);
		//if(fReadMC) FillMCTruthD(d, mcArray, selectedD0cand);
		//if(!fReadMC)  fNentries->Fill(2);
		
		//Loop over Pion tracks for B-reconstruction
		for(Int_t k = 0 ; k < aod->GetNumberOfTracks() ; k++) {
			if(fDebug>2) cout << "Beginning Pion Loop" << endl;
			Double_t xdummy,ydummy;
			Double_t dzdummy[2];
			Double_t covardummy[3];
			
			AliAODTrack *HPiAODtrk = (AliAODTrack*)aod->GetTrack(k);
			
			if(fbHIJINGavailable){
				Bool_t bInjectedPi = IsTrackInjected(HPiAODtrk,fmcHeader);
				if(fbSelectHIJING && bInjectedPi==kTRUE) {
					//if(fDebug>1) cout << "Pion candidate rejected: Not from Hijing (and we want HIJING)" << endl;
					fNentries->Fill(123);
					continue;
				}
				if(fbSelectHIJING && bInjectedPi==kFALSE) {
					//if(fDebug>1) cout << "Pion candidate from HIJING (and we want HIJING)" << endl;
					fNentries->Fill(122);
				}
				if(!fbSelectHIJING && bInjectedPi==kFALSE) {
					//if(fDebug>1) cout << "Pion candidate rejected: It is from Hijing (and we don't want HIJING)" << endl;
					fNentries->Fill(122);
					continue;
				}
				if(!fbSelectHIJING && bInjectedPi==kTRUE) {
					//if(fDebug>1) cout << "Pion candidate not from HIJING (and we don't want HIJING)" << endl;
					fNentries->Fill(123);
				}
			}
			
			
			// Checking if this track was already daughter of D0
			if(HPiAODtrk->GetID()==Ddaughter1 || HPiAODtrk->GetID()==Ddaughter2){
				if(fDebug>2)cout << "Track Rejected: This tracks is a daughter of the D candidate: " << HPiAODtrk->GetLabel() << endl;
				if(fDebug>2)cout << "D daughter labels are: " << Ddaughter1 << " and " << Ddaughter2 << endl;
				continue;
			}
			
			
			if(!fReadMC)  fNentries->Fill(3);
			Double_t BpionCharge = HPiAODtrk->Charge();
			// D0 + Pi+ --> Reject
			if(selectedD0cand==1){
				if(BpionCharge>=0) continue;
			}
			// D0bar + Pi- --> Reject
			if(selectedD0cand==2){
				if(BpionCharge<=0) continue;
			}
			// Pi0 --> Reject
			else if(BpionCharge==0) continue;
			
			// Testing Pion PID
			if ((HPiAODtrk->GetStatus()&AliESDtrack::kTPCpid ) && (HPiAODtrk->GetStatus()&AliESDtrack::kTPCin) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFpid) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFout) && (HPiAODtrk->GetStatus()&AliESDtrack::kTIME) && !(HPiAODtrk->GetStatus()&AliESDtrack::kTOFmismatch)) {
				fillthis="hBgkPionPID";
				((TH2F*)(fOutputList->FindObject(fillthis)))->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(HPiAODtrk, AliPID::kPion)),TMath::Abs(fPIDResponse->NumberOfSigmasTOF(HPiAODtrk, AliPID::kPion)));
				//  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(HPiAODtrk, AliPID::kPion)) > 2 || TMath::Abs(fPIDResponse->NumberOfSigmasTOF(HPiAODtrk, AliPID::kPion)) > 3) continue;
			} // 2.; 3.
			
			if(!fReadMC)  fNentries->Fill(4);
			/*
			 //___
			 // D0 KF vertexing selection here, as Pion(B) charge clear case for D0/D0bar
			 Int_t pdgsD[2]={0,0};
			 if(BpionCharge>0){
			 pdgsD[0] = +321; // positive Kaon candidate
			 pdgsD[1] = -211; // negative Pion candidate
			 }
			 else{
			 pdgsD[0] = +211; // positive Pion candidate
			 pdgsD[1] = -321; // negative Kaon candidate
			 }
			 if(!VertexingKFD0(d,pdgsD)) continue;
			 */
			
			
			//______________________________________________________________________________________
			// Track cuts
			
			if(fDebug>2) cout << "Applying Track Cuts" << endl;
			
			AliESDtrackCuts *trackCutsHPi = fCuts->GetTrackCuts();
			if(!fCuts->IsDaughterSelected(HPiAODtrk,(AliESDVertex*)fvtx1,trackCutsHPi)){
				if(fDebug>2) cout << " *** Track Rejected *** " << endl;
				continue;
			}
			if(fDebug>2) cout << " *** Track Accepted *** " << endl;
			
			if(!fReadMC)  fNentries->Fill(5);
			
			//______________________________________________________________________________________
			// Transforming candidates to TrackParameters
            AliExternalTrackParam *trackPi;
            trackPi->CopyFromVTrack(HPiAODtrk);
			AliNeutralTrackParam  *trackD0 = new AliNeutralTrackParam(d);
			
			//______________________________________________________________________________________
			// Construction of Secondary vertex
			if(fDebug>2) cout << "Building Secondary Vertex" << endl;
			
			TObjArray *recoArray = new TObjArray(2);
			recoArray->SetOwner(kTRUE);
			Double_t dispersion;
			
			recoArray->AddAt(trackPi,0);
			recoArray->AddAt(trackD0,1);
			
			fNentries->Fill(128); // Tries to reconstruc B vertex
			AliAODVertex *vtxAODNew = ReconstructSecondaryVertex(recoArray,dispersion,kFALSE);
			if(fDebug>2) cout << "Secondary Vertex Build" << endl;
			
			if(!vtxAODNew) {
				fNentries->Fill(129); // Fails to reconstruc B vertex
				recoArray->Clear();
				HPiAODtrk=0;
				continue;
			}
			if(!fReadMC)  fNentries->Fill(6);
			
			//______________________________________________________________________________________
			// construction of B (with secondary vertex)
			if(fDebug>2) cout << "Constructing B-Meson Candidate" << endl;
			
			const Double_t maxd = 1;
			
			//___
			// Propagate candidates to secondary vertex
			if(fDebug>2) cout << "Propagating Daughter Tracks to Secondary Vertex" << endl;
			
			Double_t px[2],py[2],pz[2],d0[2],d0err[2],dcaCand;
			
			trackPi->PropagateToDCA(vtxAODNew,fBzkG,maxd,dzdummy,covardummy);
			trackD0->PropagateToDCA(vtxAODNew,fBzkG,maxd,dzdummy,covardummy);
			
			// Calculate momenta
			if(fDebug>2) cout << "Calculating Momenta" << endl;
			Double_t momentum[3];
			trackPi->GetPxPyPz(momentum);
			px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
			trackD0->GetPxPyPz(momentum);
			px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];
			
			// Calculate impact parameters
			if(fDebug>2) cout << "Calculating Impact Parameters" << endl;
			Double_t d0z0[2],covd0z0[3];
			trackPi->PropagateToDCA(fvtx1,fBzkG,maxd,d0z0,covd0z0);
			d0[0] = d0z0[0];
			d0err[0] = TMath::Sqrt(covd0z0[0]);
			trackD0->PropagateToDCA(fvtx1,fBzkG,maxd,d0z0,covd0z0);
			d0[1] = d0z0[0];
			d0err[1] = TMath::Sqrt(covd0z0[0]);
			
			// Calculate DCA between both tracks
			if(fDebug>2) cout << "Calculating DCA between the Tracks" << endl;
			
			dcaCand = trackPi->GetDCA(trackD0,fBzkG,xdummy,ydummy);
			
			// Create Bcandidate as AliAODRecoDecayHF2Prong
			if(fDebug>2) cout << "Assigning B-Meson as HF2Prong" << endl;
			
			AliAODRecoDecayHF2Prong *BcandProng = new AliAODRecoDecayHF2Prong(vtxAODNew,px,py,pz,d0,d0err,dcaCand);
			BcandProng->SetOwnPrimaryVtx(fvtx1);
			UShort_t id[2]={(UShort_t)HPiAODtrk->GetID(),(UShort_t)trackD0->GetID()};
			BcandProng->SetProngIDs(2,id);
			
			BcandProng->SetSecondaryVtx(vtxAODNew);
			vtxAODNew->SetParent(BcandProng);
			AddDaughterRefs(vtxAODNew,(AliAODEvent*)aod,recoArray);
			// Charge for B candidate from Pion track
			BcandProng->SetCharge(HPiAODtrk->Charge());
			
			// adding daughter refs
			Int_t idHardPi=(Int_t)HPiAODtrk->GetID();
			if (idHardPi > -1) {
				vtxAODNew->AddDaughter(HPiAODtrk);
			}
			vtxAODNew->AddDaughter(trackD0); // never ever use "d" here; totally changes background
			
			if(!fReadMC)  fNentries->Fill(7);
			
			// Selection of candidates via charge combination
			if(selectedD0cand==0) {
				//delete trackPi;
				//delete trackD0;
				delete recoArray;
				if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
				delete BcandProng;
				continue;
			} // should not be possible at this stage
			if(HPiAODtrk->Charge()==0) {
				//delete trackPi;
				//delete trackD0;
				delete recoArray;
				if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
				delete BcandProng;
				continue;
			}// should not be possible at this stage
			if((selectedD0cand==1 && HPiAODtrk->Charge() > 0) || (selectedD0cand==2 && HPiAODtrk->Charge() < 0)) {
				//delete trackPi;
				//delete trackD0;
				delete recoArray;
				if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
				delete BcandProng;
				continue;
			}// (D0 & Pi+) or (D0bar & Pi-) (<-- cannot be B+/- decay)
			
			
			if(fbtrackRotation==kTRUE){
				fillthis="hRotationMonitor";
				if(fDebug>1) printf("Rotating %i times with angle %f \n",fRot,fAngle);
				TObjArray *rotatedHF2Prongs = GetArrayCandRotated(BcandProng,HPiAODtrk);
				rotatedHF2Prongs->SetOwner(kTRUE);
				for(Int_t iRot = 0; iRot<rotatedHF2Prongs->GetEntriesFast();iRot++){
					if(fDebug>1) printf("Rotation number %i has angle %f \n",iRot,fAngleFirst+ iRot*(TMath::Pi()-fAngle*(fRot - 1.)/2.));
					AliAODRecoDecayHF2Prong *BRotProng = (AliAODRecoDecayHF2Prong*)rotatedHF2Prongs->At(iRot);
					if(!BRotProng){
						if(fDebug>1)printf("Rotation Candidate Rejected \n");
						delete BRotProng; BRotProng=0x0;
						continue;
					}
					BRotProng->SetOwnPrimaryVtx(fvtx1);
					if(fDebug>1) printf("Info on candidate \n");
					if(fDebug>1) BRotProng->Dump();
					if(fDebug>1) printf("Done: Info on candidate \n");
					if(fDebug>1) BRotProng->GetSecondaryVtx()->Print();
					((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(0);
					
					if(fKillCandidates){
						if(KillCandidates(BRotProng,selectedD0cand,d)==kFALSE){//,HPiAODtrk,d)==kFALSE){
							((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(6);
							BRotProng->UnsetOwnPrimaryVtx();
							AliAODVertex* vtxt = (AliAODVertex*)BRotProng->GetSecondaryVtx();
							delete vtxt; vtxt=0x0;
							continue;
						}
					}
					FillRotateTrackHists(BRotProng,d,selectedD0cand);
					BRotProng->UnsetOwnPrimaryVtx();
					AliAODVertex* vtxt = (AliAODVertex*)BRotProng->GetSecondaryVtx();
					delete vtxt; vtxt=0x0;
				}
				if(fDebug>1)printf("Clearing rotation array \n");
				rotatedHF2Prongs->Delete();
				delete rotatedHF2Prongs;
				rotatedHF2Prongs=0x0;
			}
			
			// B Meson selection cuts
			if(fKillCandidates){
				fNentries->Fill(105);
				if(KillCandidates(BcandProng,selectedD0cand,d)==kFALSE){//,HPiAODtrk,d)==kFALSE){
					//delete trackPi;
					//delete trackD0;
					delete recoArray;
					if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
					delete BcandProng;
					continue;
				}
			}
			// ---
			// Acceptance Cut on B
			fillthis="Bfiducial_before";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(BcandProng->Y(521));
			
			// After analysis of 'B vs Pt vs Y', we accept all B |y|<0.8
			if(TMath::Abs(BcandProng->Y(521))>0.8) continue;
			
			fillthis="Bfiducial_after";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(BcandProng->Y(521));
			// ---
			
			// Filling candidate histos
			if(!fReadMC){
				if(FillBpmHists(BcandProng,selectedD0cand,HPiAODtrk,d)==kFALSE){
					//delete trackPi;
					//delete trackD0;
					delete recoArray;
					if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
					delete BcandProng;
					continue;
				}
				FillDHists(BcandProng,selectedD0cand,HPiAODtrk,d);
			}
			if(fReadMC){
				if(FillMCTruthHistos(d,selectedD0cand,HPiAODtrk,BcandProng)==kTRUE){
					if((HPiAODtrk->GetStatus()&AliESDtrack::kTPCpid ) && (HPiAODtrk->GetStatus()&AliESDtrack::kTPCin) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFpid) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFout) && (HPiAODtrk->GetStatus()&AliESDtrack::kTIME) && !(HPiAODtrk->GetStatus()&AliESDtrack::kTOFmismatch)) {
						fillthis="hSgnPionPID";
						((TH2F*)(fOutputList->FindObject(fillthis)))->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(HPiAODtrk, AliPID::kPion)),TMath::Abs(fPIDResponse->NumberOfSigmasTOF(HPiAODtrk, AliPID::kPion)));
					}
					if(fbSelectHIJING==kFALSE){
						fRunNumber->Fill(aod->GetRunNumber());
						FillBpmMCHists(BcandProng,selectedD0cand,HPiAODtrk,d);
						FillMCTruthD(d,selectedD0cand);
						FillMCtruthDHists(BcandProng,selectedD0cand,HPiAODtrk,d);
					}
					else{
						delete recoArray;
						if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
						delete BcandProng;
						continue;
					}
				}
				else{
					if(fbSelectHIJING==kTRUE){
						FillDHists(BcandProng,selectedD0cand,HPiAODtrk,d);
						if(FillBpmHists(BcandProng,selectedD0cand,HPiAODtrk,d)==kTRUE){
							//CheckMCHistory(d,fMcArray,HPiAODtrk);
						}
					}
					else{
						//delete trackPi;
						//delete trackD0;
						delete recoArray;
						if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
						delete BcandProng;
						continue;
					}
				}
			}
			
			BcandProng->UnsetOwnPrimaryVtx();
			if(fDebug>2) cout << "__________________________________________*Done*__________________________________________" << endl;
			
			//
			//delete trackPi;
			//delete trackD0;
			delete recoArray;
			if(vtxAODNew){delete vtxAODNew; vtxAODNew=NULL;}
			delete BcandProng;
		}
	} //end for prongs
	
	// Post the data
	PostData(1,fOutputList);
	PostData(2,fOutputListVar);
	PostData(3,fOutputListMCVar);
	PostData(4,fOutputListDSandB);
	PostData(5,fOutputListDSonly);
	PostData(6,fSelectionVariables);
	PostData(7,fRotatedCandidates);
	return;
}
//____________________________________________________________________________
void AliAnalysisTaskSEBpmMass::FillMCTruthD(AliAODRecoDecayHF2Prong* d, Double_t sD0c){//,AliAODRecoDecayHF2Prong* bMeson){
	
	// MC Truth for D Mesons
	TString fillthis="";
	Int_t nbinsD=fCuts->PtBin(d->Pt());
	
	Int_t pdgDaughtersD0[2]={321,211};
	Int_t labDpD0 = d->MatchToMC(421,fMcArray,2,pdgDaughtersD0);
	if(labDpD0>0){
		AliAODMCParticle *partD0= (AliAODMCParticle*)fMcArray->At(labDpD0);
		if(!partD0){
			return;
		}
		Int_t D0mum = partD0->GetMother();
		if(D0mum==-1)return;
		AliAODMCParticle *partMD0= (AliAODMCParticle*)fMcArray->At(D0mum);
		if(!partMD0){
			return;
		}
		Int_t pdgsBm = TMath::Abs(partMD0->GetPdgCode());
		if(pdgsBm!=521){
			return;
		}
		if(partMD0->GetNDaughters()!=2){
			return;
		}
		Int_t dID = 0;
		AliAODMCParticle *dcheck = (AliAODMCParticle*)fMcArray->At(partMD0->GetDaughter(0));
		if(dcheck->GetLabel()==partD0->GetLabel()){
			dID=1;
		}
		AliAODMCParticle *partPi = (AliAODMCParticle*)fMcArray->At(partMD0->GetDaughter(dID));
		if(TMath::Abs(partPi->GetPdgCode())!=211){
			return;
		}
		else {
			Double_t posRecoD[3],posTrueD[3];
			
			d->GetSecondaryVtx(posRecoD);
			AliAODTrack *trk = (AliAODTrack*)d->GetDaughter(0);
			AliAODMCParticle *partDghD= (AliAODMCParticle*)fMcArray->At(TMath::Abs(trk->GetLabel()));
			if(!partDghD){ return;}
			partDghD->XvYvZv(posTrueD);
			
			
			if(sD0c==0) return;
			else if(sD0c==1 || (sD0c==3 && partPi->GetPdgCode()==-211)){
				fNentries->Fill(26+nbinsD);
				
				fillthis="DiffPosDX_";
				fillthis+=nbinsD;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posRecoD[0]-posTrueD[0]);
				fillthis="DiffPosDY_";
				fillthis+=nbinsD;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posRecoD[1]-posTrueD[1]);
				fillthis="DiffPosDZ_";
				fillthis+=nbinsD;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posRecoD[2]-posTrueD[2]);
			}
			else if(sD0c==2 || (sD0c==3 && partPi->GetPdgCode()== 211)){
				fNentries->Fill(1+nbinsD);
				
				fillthis="DiffPosDX_";
				fillthis+=nbinsD;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posRecoD[0]-posTrueD[0]);
				fillthis="DiffPosDY_";
				fillthis+=nbinsD;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posRecoD[1]-posTrueD[1]);
				fillthis="DiffPosDZ_";
				fillthis+=nbinsD;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posRecoD[2]-posTrueD[2]);
			}
		}
	}
	else {
		return;
	}
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::KillDCandidates(AliAODRecoDecayHF2Prong* Bmeson, AliAODRecoDecayHF2Prong* Dmeson, Double_t sD0c){
	
	// Initialize cut variables for candidates
	
	//Double_t cosThStaD,PtPiD,PtKaD,d0PiD,d0KaD;
	
	//Double_t cosPtAngD		= Dmeson->CosPointingAngleXY();
	//Double_t cosPtAngDXY	= Dmeson->CosPointingAngleXY();
	//Double_t DCAD			= Dmeson->GetDCA();
	//Double_t Prodd0d0D		= Dmeson->Prodd0d0();
	//Double_t LthD			= Dmeson->DecayLength();
	//Double_t normDecLthD	= Dmeson->NormalizedDecayLength();
	//Double_t LthXYD			= Dmeson->DecayLengthXY();
	//Double_t normDecLthXYD	= Dmeson->NormalizedDecayLengthXY();
	Double_t Dpt			= Dmeson->Pt();
	
	if(sD0c==0) return 0;
	/*
	 else if(sD0c==1 || sD0c==3){
	 cosThStaD	= Dmeson->CosThetaStarD0();
	 PtPiD		= Dmeson->PtProng(0);
	 PtKaD		= Dmeson->PtProng(1);
	 d0PiD		= Dmeson->Getd0Prong(0);
	 d0KaD		= Dmeson->Getd0Prong(1);
	 }
	 else if(sD0c==2){
	 cosThStaD   = Dmeson->CosThetaStarD0bar();
	 PtPiD		= Dmeson->PtProng(1);
	 PtKaD		= Dmeson->PtProng(0);
	 d0PiD		= Dmeson->Getd0Prong(1);
	 d0KaD		= Dmeson->Getd0Prong(0);
	 }*/
	
	// first cuts tuned individually such, that <10% of Signal reduction is seen
	
	Int_t nbins = fCuts->PtBin(Bmeson->Pt());
	if((nbins==0 || nbins==1) /*&& Dpt<2.0*/) return 1;			// In first two pt bins, only D0 pt cut applied
	if(nbins==2 &&
		 (//cosPtAngD<0.85		||
			//cosPtAngDXY<0.95		||
			//DCAD>0.018			||
			//LthD<0.02				||
			//LthXYD<0.02			||
			//normDecLthD<5.0		||
			//normDecLthXYD<3.3	||
			//Prodd0d0D>-0.000016	||
			Dpt<2.8)) return 0;
	if(nbins==3 &&
		 (//cosPtAngD<0.85		||
			//cosPtAngDXY<0.95	||
			//DCAD>0.018			||
			//LthD<0.02			||
			//LthXYD<0.02			||
			//normDecLthD<6.0		||
			//normDecLthXYD<3.3	||
			//Prodd0d0D>-0.0002	||
			Dpt<3.6)) return 0;
	if(nbins==4 &&
		 (//cosPtAngD<0.90		||
			//cosPtAngDXY<0.95	||
			//DCAD>0.018			||
			//LthD<0.02			||
			//LthXYD<0.02			||
			//normDecLthD<5.0		||
			//normDecLthXYD<3.3	||
			//Prodd0d0D>-0.00004	||
			Dpt<2.8)) return 0;
	if(nbins==5 &&
		 (//cosPtAngD<0.90		||
			//cosPtAngDXY<0.95	||
			//DCAD>0.018			||
			//LthD<0.03			||
			//LthXYD<0.03			||
			//normDecLthD<7.0		||
			//normDecLthXYD<3.5	||
			//Prodd0d0D>-0.00004	||
			Dpt<3.6)) return 0;
	if(nbins==6 &&
		 (//cosPtAngD<0.90		||
			//cosPtAngDXY<0.95	||
			//DCAD>0.015			||
			//LthD<0.03			||
			//LthXYD<0.03			||
			//normDecLthD<6.0		||
			//normDecLthXYD<3.5	||
			//Prodd0d0D>-0.00012	||
			Dpt<2.8)) return 0;
	if(nbins==7 &&
		 (//cosPtAngD<0.92		||
			//cosPtAngDXY<0.95	||
			//DCAD>0.015			||
			//LthD<0.03			||
			//LthXYD<0.03			||
			//normDecLthD<6.0		||
			//normDecLthXYD<3.5	||
			//Prodd0d0D>-0.00012	||
			Dpt<2.8)) return 0;
	if(nbins==8 &&
		 (//cosPtAngD<0.92		||
			//cosPtAngDXY<0.96	||
			//DCAD>0.015			||
			//LthD<0.03			||
			//LthXYD<0.03			||
			//normDecLthD<5.0		||
			//normDecLthXYD<3.5	||
			//Prodd0d0D>-0.00016	||
			Dpt<2.8)) return 0;
	if(nbins==9 &&
		 (//cosPtAngD<0.95		||
			//cosPtAngDXY<0.96	||
			//DCAD>0.015			||
			//LthD<0.05			||
			//LthXYD<0.05			||
			//normDecLthD<5.0		||
			//normDecLthXYD<5.0	||
			//Prodd0d0D>-0.00016	||
			Dpt<3.6)) return 0;
	if(nbins==10 &&
		 (//cosPtAngD<0.95		||
			//cosPtAngDXY<0.98	||
			//DCAD>0.013			||
			//LthD<0.05			||
			//LthXYD<0.05			||
			//normDecLthD<5.0		||
			//normDecLthXYD<5.0	||
			//Prodd0d0D>-0.00024	||
			Dpt<2.8)) return 0;
	if(nbins==11 &&
		 (//cosPtAngD<0.95		||
			//cosPtAngDXY<0.98	||
			//DCAD>0.013			||
			//LthD<0.06			||
			//LthXYD<0.06			||
			//normDecLthD<5.0		||
			//normDecLthXYD<5.0	||
			//Prodd0d0D>-0.00012	||
			Dpt<3.6)) return 0;
	if(nbins==12 &&
		 (//cosPtAngD<0.95		||
			//cosPtAngDXY<0.98	||
			//DCAD>0.013			||
			//LthD<0.1			||
			//LthXYD<0.1			||
			//normDecLthD<8.0		||
			//normDecLthXYD<6.5	||
			//Prodd0d0D>-0.0002	||
			Dpt<2.8)) return 0;
	
	return 1;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::KillCandidates(AliAODRecoDecayHF2Prong* Bmeson, Double_t sD0c, AliAODRecoDecayHF2Prong* /*Dmeson*/){
	
	// first cuts estimated from MC vs Bgk comparison
	// Initialize cut variables for candidates
	
	if(sD0c==0) return 0;
	
	UInt_t pdgB[2]={0,0};
	pdgB[0] = 211;
	pdgB[1] = 421;
	
	Double_t threeSigmaB[13]={0.0678618,0.0726664,0.072042,0.0716596,0.0713123,0.0774529,0.0765073,0.0835032,0.0854299,0.093386,0.0963455,0.103999,0.124467};
	
	Double_t cut_B_ctp[13]={0.7,0.9,0.98,0.99,0.995,0.995,0.995,0.995,0.995,0.99,0.99,0.99,0.99};
	
	Double_t B_invMass	= Bmeson->InvMass(2,pdgB); // three sigma invariant mass selection
	Int_t nbins = fCuts->PtBin(Bmeson->Pt());
	if(TMath::Abs(B_invMass-5.279)>threeSigmaB[nbins])return 0;
	
	/*
	Double_t B_ctp = Bmeson->CosPointingAngle();
	if(B_ctp<cut_B_ctp[nbins])return 0;
	
	Double_t D_dl = Dmeson->DecayLength();
	if(D_dl<0.00)return 0;
	
	Double_t B_ptD0 = Bmeson->PtProng(1);
	if(B_ptD0<0.0)return 0;
	
	Double_t B_ptPi = Bmeson->PtProng(0);
	if(B_ptPi<0.50)return 0;
	
	Double_t B_d0xy = TMath::Abs(Bmeson->ImpParXY());
	if(B_d0xy>0.1)return 0;
	
	Double_t B_d0d0 = Bmeson->Prodd0d0();
	if(B_d0d0>1.E-3)return 0;
	*/
	return 1;
	
	/*
	 if(cosPtAngB<0.90)return 0;
	 if(LthNormD<0.0)return 0;
	 if(PtDB<0.5)return 0;
	 if(PtPiB<0.5)return 0;
	 //	if(ImpParXY>0.01)return 0;
	 
	 return 1;
	 */
	/*
	 Int_t nbins = fCuts->PtBin(Bmeson->Pt());
	 if(nbins==0) return 0;			// In first pt bin, no cuts applied, because no FONLL prediction available
	 if(nbins==1) return 0;			// In first pt bin, no cuts applied, because no FONLL prediction available
	 if(nbins==2) return 0;			// In first pt bin, no cuts applied, because no FONLL prediction available
	 if(nbins==1 && (cosPtAngB<fBCtpCut[1]	|| LthNormD<fDNdlCut[1]	|| PtDB<fD0PtCut[1]	|| ImpParXY>fBIxyCut[1] || PtPiB<fPiPtCut[1])) return 0;
	 if(nbins==2 && (cosPtAngB<fBCtpCut[2]	|| LthNormD<fDNdlCut[2]	|| PtDB<fD0PtCut[2]	|| ImpParXY>fBIxyCut[2] || PtPiB<fPiPtCut[2])) return 0;
	 if(nbins==3 && (cosPtAngB<fBCtpCut[3]	|| LthNormD<fDNdlCut[3]	|| PtDB<fD0PtCut[3]	|| ImpParXY>fBIxyCut[3] || PtPiB<fPiPtCut[3])) return 0;
	 if(nbins==4 && (cosPtAngB<fBCtpCut[4]	|| LthNormD<fDNdlCut[4]	|| PtDB<fD0PtCut[4]	|| ImpParXY>fBIxyCut[4] || PtPiB<fPiPtCut[4])) return 0;
	 if(nbins==5 && (cosPtAngB<fBCtpCut[5]	|| LthNormD<fDNdlCut[5]	|| PtDB<fD0PtCut[5]	|| ImpParXY>fBIxyCut[5] || PtPiB<fPiPtCut[5])) return 0;
	 if(nbins==6 && (cosPtAngB<fBCtpCut[6]	|| LthNormD<fDNdlCut[6]	|| PtDB<fD0PtCut[6]	|| ImpParXY>fBIxyCut[6] || PtPiB<fPiPtCut[6])) return 0;
	 if(nbins==7 && (cosPtAngB<fBCtpCut[7]	|| LthNormD<fDNdlCut[7]	|| PtDB<fD0PtCut[7]	|| ImpParXY>fBIxyCut[7] || PtPiB<fPiPtCut[7])) return 0;
	 if(nbins==8 && (cosPtAngB<fBCtpCut[8]	|| LthNormD<fDNdlCut[8]	|| PtDB<fD0PtCut[8]	|| ImpParXY>fBIxyCut[8] || PtPiB<fPiPtCut[8])) return 0;
	 if(nbins==9 && (cosPtAngB<fBCtpCut[9]	|| LthNormD<fDNdlCut[9]	|| PtDB<fD0PtCut[9]	|| ImpParXY>fBIxyCut[9] || PtPiB<fPiPtCut[9])) return 0;
	 if(nbins==10 && (cosPtAngB<fBCtpCut[10]	|| LthNormD<fDNdlCut[10] || PtDB<fD0PtCut[10] || ImpParXY>fBIxyCut[10] || PtPiB<fPiPtCut[10])) return 0;
	 if(nbins==11 && (cosPtAngB<fBCtpCut[11]	|| LthNormD<fDNdlCut[11] || PtDB<fD0PtCut[11] || ImpParXY>fBIxyCut[11] || PtPiB<fPiPtCut[11])) return 0;
	 if(nbins==12 && (cosPtAngB<fBCtpCut[12]	|| LthNormD<fDNdlCut[12] || PtDB<fD0PtCut[12] || ImpParXY>fBIxyCut[12] || PtPiB<fPiPtCut[12])) return 0;
	 */
	// intermediate cuts
	/*
	 if(nbins==1   && (cosPtAngB<0.50000 || LthNormD<0. || PtDB<0.00 || PtPiB<0.5)) return 0;
	 if(nbins==2   && (cosPtAngB<0.99000 || LthNormD<0.  || PtDB<0.00 || ImpParXY>0.00875 || PtPiB<0.5)) return 0;
	 if(nbins==3   && (cosPtAngB<0.80000 || LthNormD<10. || PtDB<0.00 || ImpParXY>0.00450 || PtPiB<0.5)) return 0;
	 if(nbins==4   && (cosPtAngB<0.99500 || LthNormD<5.0 || PtDB<0.00 || ImpParXY>0.02000 || PtPiB<0.5)) return 0;
	 if(nbins==5   && (cosPtAngB<0.99900 || LthNormD<0.  || PtDB<0.00 || ImpParXY>0.00800 || PtPiB<0.5)) return 0;
	 if(nbins==6   && (cosPtAngB<0.99950 || LthNormD<5.  || PtDB<0.00 || ImpParXY>0.0055  || PtPiB<0.5)) return 0;
	 if(nbins==7   && (cosPtAngB<0.99950 || LthNormD<10. || PtDB<0.00 || ImpParXY>0.005   || PtPiB<0.5)) return 0;
	 if(nbins==8   && (cosPtAngB<0.99950 || LthNormD<15. || PtDB<1.00 || ImpParXY>0.00500 || PtPiB<0.5)) return 0;
	 if(nbins==9   && (cosPtAngB<0.99900 || LthNormD<5.  || PtDB<0.00 || ImpParXY>0.00500 || PtPiB<0.5)) return 0;
	 if(nbins==10  && (cosPtAngB<0.99950 || LthNormD<10. || PtDB<0.00 || ImpParXY>0.00625 || PtPiB<0.5)) return 0;
	 if(nbins==11  && (cosPtAngB<0.99500 || LthNormD<4.0 || PtDB<0.00 || ImpParXY>0.01 || PtPiB<0.5)) return 0;
	 if(nbins==12  && (cosPtAngB<0.99500 || LthNormD<0.0 || PtDB<0.00 || PtPiB<0.5)) return 0;
	 */
	// analysis cuts
	/*	if(nbins==1   && (cosPtAngB<0.99300 || LthNormD<23. || PtDB<0.00 || ImpParXY>0.00840 || PtPiB<0.5)) return 0;
	 if(nbins==2   && (cosPtAngB<0.99650 || LthNormD<13. || PtDB<1.40 || ImpParXY>0.00635 || PtPiB<1.8)) return 0;
	 if(nbins==3   && (cosPtAngB<0.99800 || LthNormD<29. || PtDB<0.70 || ImpParXY>0.00455 || PtPiB<1.4)) return 0;
	 if(nbins==4   && (cosPtAngB<0.99875 || LthNormD<20. || PtDB<1.00 || ImpParXY>0.00625 || PtPiB<1.4)) return 0;
	 if(nbins==5   && (cosPtAngB<0.99975 || LthNormD<13. || PtDB<0.10 || ImpParXY>0.00080 || PtPiB<1.1)) return 0;
	 if(nbins==6   && (cosPtAngB<0.99975 || LthNormD<12. || PtDB<0.40 || ImpParXY>0.00470 || PtPiB<1.1)) return 0;
	 if(nbins==7   && (cosPtAngB<0.99890 || LthNormD<19. || PtDB<1.80 || ImpParXY>0.00200 || PtPiB<1.2)) return 0;
	 if(nbins==8   && (cosPtAngB<0.999875|| LthNormD<20. || PtDB<0.00 || ImpParXY>0.00620 || PtPiB<0.9)) return 0;
	 if(nbins==9   && (cosPtAngB<0.99980 || LthNormD<9.  || PtDB<1.30 || ImpParXY>0.00650 || PtPiB<1.3)) return 0;
	 if(nbins==10  && (cosPtAngB<0.99850 || LthNormD<13. || PtDB<2.25 || ImpParXY>0.00350 || PtPiB<1.2)) return 0;
	 if(nbins==11  && (cosPtAngB<0.99800 || LthNormD<3.75|| PtDB<3.40 || ImpParXY>0.00425 || PtPiB<1.1)) return 0;
	 if(nbins==12  && (cosPtAngB<0.99800 || LthNormD<5.0 || PtDB<4.40 || ImpParXY>0.00600 || PtPiB<0.8)) return 0;
	 */
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::KillCandidatesAdvanced(AliAODRecoDecayHF2Prong* Bmeson, AliAODRecoDecayHF2Prong* Dmeson, Bool_t MCtruth){
	//
	// Advanced cut selection based on topology of D decay point with respect ot B decay point and momentum direction
	//
	
	TString fillthis="";
	Int_t nbins=fCuts->PtBin(Bmeson->Pt());
	
	TVector3 decVecB(Bmeson->GetSecVtxX(),Bmeson->GetSecVtxY(),Bmeson->GetSecVtxZ());
	TVector3 decVecD(Dmeson->GetSecVtxX(),Dmeson->GetSecVtxY(),Dmeson->GetSecVtxZ());
	
	Double_t tempSP = (decVecD.Dot(decVecB));
	
	TVector3 Dprojection(tempSP/(decVecB.Dot(decVecB))*Bmeson->Px(),tempSP/(decVecB.Dot(decVecB))*Bmeson->Py(),tempSP/(decVecB.Dot(decVecB))*Bmeson->Pz());
	
	// If sign changed, vector projection is in opposite direction, namely the D decay point is in wrong direction w.r.t. B decay point and B momentum (momentum defines axis)
	if(tempSP <= 0. || Dprojection.Mag() <= decVecB.Mag()){
		
		// Fill Histogram, and special hist if MC availa
		if(MCtruth){
			fNentries->Fill(126);
			fillthis="advancedCutMC_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(0);
		}
		else{
			fNentries->Fill(124);
			fillthis="advancedCutAll_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(0);
		}
		return kTRUE;
	}
	// Fill Histogram
	// for now always return kTRUE to monitor distribution
	if(MCtruth){
		fNentries->Fill(127);
		fillthis="advancedCutMC_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(1);
	}
	else{
		fNentries->Fill(125);
		fillthis="advancedCutAll_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(1);
	}
	
	return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::FillMCTruthHistos(AliAODRecoDecayHF2Prong* d,Double_t sD0c, AliAODTrack* HPiAODtrk,AliAODRecoDecayHF2Prong* bMeson){
	//
	// function used for MC truth of B candidate
	// check MC truth of D0 candidate and GetMother()
	//
	// Match to MC
	
	if(fDebug>2) cout << "Filling MCTruthHistos" << endl;
	if(fDebug>2) cout << "____________________________________________________________________________________" << endl;
	
	if(HPiAODtrk->Charge()==0) {
		if(fDebug>2) cout << "Pion charge is 0, which should not be possible!" << endl;
		return 0;
	}
	fMCTruth->Fill(0);
	
	Double_t posReco[3],posTrue[3];
	TString fillthis="";
	Int_t nbins  = fCuts->PtBin(bMeson->Pt());
	
	if(fReadMC) {
		AliAODMCParticle *partPi= (AliAODMCParticle*)fMcArray->At(TMath::Abs(HPiAODtrk->GetLabel()));
		if(partPi){
			bPiPhysPrimary	= partPi->IsPhysicalPrimary();
			bPiWeak					= partPi->IsSecondaryFromWeakDecay();
			bPiPdgCode			= partPi->GetPdgCode();

//			bPiFastWeak			= IsFastMcFromWeakDecay(partPi);
//			if(partPi->IsSecondaryFromMaterial()) bPiFastWeak =0;
			//			bPiMaterial			= partPi->IsSecondaryFromMaterial();
		}
		AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
		AliAODMCParticle *DaughterPart0 = (AliAODMCParticle*)fMcArray->At(TMath::Abs(trk0->GetLabel()));
		if(DaughterPart0){
			bDDgh0PhysPrimary		= DaughterPart0->IsPhysicalPrimary();
			bDDgh0Weak					= DaughterPart0->IsSecondaryFromWeakDecay();
			bDDgh0PdgCode				= DaughterPart0->GetPdgCode();
			//			bDDgh0FastWeak			= IsFastMcFromWeakDecay(DaughterPart0);
//			if(DaughterPart0->IsSecondaryFromMaterial()) bDDgh0FastWeak =0;
	//			bDDgh0Material			= DaughterPart0->IsSecondaryFromMaterial();
		}
		
		AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
		AliAODMCParticle *DaughterPart1 = (AliAODMCParticle*)fMcArray->At(TMath::Abs(trk1->GetLabel()));
		if(DaughterPart1){
			bDDgh1PhysPrimary		= DaughterPart1->IsPhysicalPrimary();
			bDDgh1Weak					= DaughterPart1->IsSecondaryFromWeakDecay();
			bDDgh1PdgCode				= DaughterPart1->GetPdgCode();
//			bDDgh1FastWeak			= IsFastMcFromWeakDecay(DaughterPart1);
//			if(DaughterPart1->IsSecondaryFromMaterial()) bDDgh1FastWeak =0;
			//			bDDgh1Material			= DaughterPart1->IsSecondaryFromMaterial();
		}
	}
	
	if((HPiAODtrk->GetStatus()&AliESDtrack::kTPCpid ) && (HPiAODtrk->GetStatus()&AliESDtrack::kTPCin) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFpid) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFout) && (HPiAODtrk->GetStatus()&AliESDtrack::kTIME) && !(HPiAODtrk->GetStatus()&AliESDtrack::kTOFmismatch)) {
		bBPiTOFnSigma = fPIDResponse->NumberOfSigmasTOF(HPiAODtrk, AliPID::kPion);
		bBPiTPCnSigma = fPIDResponse->NumberOfSigmasTPC(HPiAODtrk, AliPID::kPion);
	}
	
	if(sD0c==0){
		return 0;
	}
	else if(sD0c==1){
		fMCTruth->Fill(1);
		Int_t pdgDaughtersD0[2]={321,211};
		Int_t pdgsBm=0;
		Int_t pdgsMPi=0;
		Int_t pdgsPi=0;
		Int_t labDpD0 = d->MatchToMC(421,fMcArray,2,pdgDaughtersD0);
		if(labDpD0>0){
			if(fDebug>2) printf("Identified D0 Meson \n");
			fMCTruth->Fill(4);
			AliAODMCParticle *partD0= (AliAODMCParticle*)fMcArray->At(labDpD0);
			if(!partD0){
				return 0;
			}
			if(fDebug>2) printf("D0 candidate found at label: %i \n",partD0->GetLabel());
			// Get mother of the D0	from MC stack
			Int_t labMD0 = partD0->GetMother();
			AliAODMCParticle *partMD0= (AliAODMCParticle*)fMcArray->At(labMD0);
			if(!partMD0){
				return 0;
			}
			
			Bool_t bfromB = kFALSE;
			if(fbHIJINGavailable){bfromB = IsCandidateInjected(d,fmcHeader);}
			
			if(fDebug>2) printf("D0 Mother candidate found at label: %i \n",partMD0->GetLabel());
			// Check if mother is B-; otherwise --> reject
			pdgsBm=partMD0->GetPdgCode();
			if(pdgsBm!=-521){
				if(bfromB==kTRUE) fNentries->Fill(42);
				if(bfromB==kFALSE) fNentries->Fill(74);
				return 0;
			}
			// Check if B- has 2 daughters; otherwise --> reject
			if(partMD0->GetNDaughters()!=2){
				if(bfromB==kTRUE) fNentries->Fill(43);
				if(bfromB==kFALSE) fNentries->Fill(75);
				return 0;
			}
			if(fDebug>2) printf("D0 Mother candidate is %i and decays into %i daughters \n",pdgsBm,partMD0->GetNDaughters());
			
			// Get charged track from MC stack
			AliAODMCParticle *partPi= (AliAODMCParticle*)fMcArray->At(TMath::Abs(HPiAODtrk->GetLabel()));
			if(!partPi){
				if(bfromB==kTRUE) fNentries->Fill(44);
				if(bfromB==kFALSE) fNentries->Fill(76);
				return 0;
			}
			// Check if Pi- otherwise --> reject
			pdgsPi=partPi->GetPdgCode();
			if(pdgsPi!=-211){
				if(bfromB==kTRUE) fNentries->Fill(45);
				if(bfromB==kFALSE) fNentries->Fill(77);
				return 0;
			}
			if(fDebug>2) printf("%i track found at label: %i \n",pdgsPi,partPi->GetLabel());
			// Get mother of the Pi- from MC stack
			Int_t labMPi = partPi->GetMother();
			AliAODMCParticle *partMPi= (AliAODMCParticle*)fMcArray->At(labMPi);
			if(!partMPi){
				if(bfromB==kTRUE) fNentries->Fill(46);
				if(bfromB==kFALSE) fNentries->Fill(78);
				return 0;
			}
			// Check if mother is B-; otherwise --> reject
			pdgsMPi=partMPi->GetPdgCode();
			if(pdgsMPi!=-521){
				if(bfromB==kTRUE) fNentries->Fill(47);
				if(bfromB==kFALSE) fNentries->Fill(79);
				return 0;
			}
			if(fDebug>2) printf("Pion mother candidate found at label: %i \n",partMPi->GetLabel());
			if(fDebug>2) printf("Pion mother candidate is %i and decays into %i daughters \n",pdgsMPi,partMPi->GetNDaughters());
			if(partMD0->GetLabel()==partMPi->GetLabel()) {
				fMCTruth->Fill(8);
				printf("Complete B- --> D0 + Pi- --> Pi+ + K- detected!!! \n");
				if(bfromB==kFALSE){
					cout << "### Particle injected by hand!" << endl;
					fNentries->Fill(40);
				}
				else if(bfromB==kTRUE){
					cout << "### Particle from quark!" << endl;
					fNentries->Fill(39);
				}
				// Fill Plots for Secondary Vertex Resolution
				// for B meson
				partPi->XvYvZv(posTrue);
				bMeson->GetSecondaryVtx(posReco);
				
				fillthis="DiffPosX_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[0]-posTrue[0]);
				fillthis="DiffPosY_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[1]-posTrue[1]);
				fillthis="DiffPosZ_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[2]-posTrue[2]);
				
				return 1;
			}
			else{
				if(bfromB==kTRUE) fNentries->Fill(69);
				if(bfromB==kFALSE) fNentries->Fill(101);
				return 0;
			}
		}
	}
	else if(sD0c==2){
		fMCTruth->Fill(2);
		Int_t pdgDaughtersD0bar[2]={321,211};
		Int_t pdgsBp=0;
		Int_t pdgsMPi=0;
		Int_t pdgsPi=0;
		Int_t labDpD0bar = d->MatchToMC(421,fMcArray,2,pdgDaughtersD0bar);
		if(labDpD0bar>0){
			if(fDebug>2)   printf("Identified D0bar Meson \n");
			fMCTruth->Fill(5);
			AliAODMCParticle *partD0bar= (AliAODMCParticle*)fMcArray->At(labDpD0bar);
			if(!partD0bar) return 0;
			if(fDebug>2) printf("D0bar candidate found at label: %i \n",partD0bar->GetLabel());
			// Get mother of the D0	from MC stack
			Int_t labMD0bar = partD0bar->GetMother();
			AliAODMCParticle *partMD0bar= (AliAODMCParticle*)fMcArray->At(labMD0bar);
			if(!partMD0bar) return 0;
			
			Bool_t bfromB = kFALSE;
			if(fbHIJINGavailable){bfromB = IsCandidateInjected(d,fmcHeader);}
			
			if(fDebug>2)   printf("D0bar Mother candidate found at label: %i \n",partMD0bar->GetLabel());
			// Check if mother is B-; otherwise --> reject
			pdgsBp=partMD0bar->GetPdgCode();
			if(pdgsBp!=521){
				if(bfromB==kTRUE) fNentries->Fill(49);
				if(bfromB==kFALSE) fNentries->Fill(81);
				return 0;
			}
			// Check if B- has 2 daughters; otherwise --> reject
			if(partMD0bar->GetNDaughters()!=2){
				if(bfromB==kTRUE) fNentries->Fill(50);
				if(bfromB==kFALSE) fNentries->Fill(82);
				return 0;
			}
			if(fDebug>2)   printf("D0bar Mother candidate is %i and decays into %i daughters \n",pdgsBp,partMD0bar->GetNDaughters());
			
			// Get charged track from MC stack
			AliAODMCParticle *partPi= (AliAODMCParticle*)fMcArray->At(TMath::Abs(HPiAODtrk->GetLabel()));
			if(!partPi){
				if(bfromB==kTRUE) fNentries->Fill(51);
				if(bfromB==kFALSE) fNentries->Fill(83);
				return 0;
			}
			// Check if Pi+ otherwise --> reject
			pdgsPi=partPi->GetPdgCode();
			if(pdgsPi!=211){
				if(bfromB==kTRUE) fNentries->Fill(52);
				if(bfromB==kFALSE) fNentries->Fill(84);
				return 0;
			}
			if(fDebug>2)   printf("%i track found at label: %i \n",pdgsPi,partPi->GetLabel());
			// Get mother of the Pi+ from MC stack
			Int_t labMPi = partPi->GetMother();
			AliAODMCParticle *partMPi= (AliAODMCParticle*)fMcArray->At(labMPi);
			if(!partMPi){
				if(bfromB==kTRUE) fNentries->Fill(53);
				if(bfromB==kFALSE) fNentries->Fill(85);
				return 0;
			}
			// Check if mother is B+; otherwise --> reject
			pdgsMPi=partMPi->GetPdgCode();
			if(pdgsMPi!=521){
				if(bfromB==kTRUE) fNentries->Fill(54);
				if(bfromB==kFALSE) fNentries->Fill(86);
				return 0;
			}
			if(fDebug>2)   printf("Pion mother candidate found at label: %i \n",partMPi->GetLabel());
			if(fDebug>2)   printf("Pion mother candidate is %i and decays into %i daughters \n",pdgsMPi,partMPi->GetNDaughters());
			if(partMD0bar->GetLabel()==partMPi->GetLabel()) {
				fMCTruth->Fill(9);
				printf("Complete B+ --> D0bar + Pi+ --> K+ + Pi- detected!!! \n");
				if(bfromB==kFALSE){
					cout << "### Particle injected by hand!" << endl;
					fNentries->Fill(40);
				}
				else if(bfromB==kTRUE){
					cout << "### Particle from quark!" << endl;
					fNentries->Fill(39);
				}
				// Fill Plots for Secondary Vertex Resolution
				// for B meson
				partPi->XvYvZv(posTrue);
				bMeson->GetSecondaryVtx(posReco);
				
				fillthis="DiffPosX_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[0]-posTrue[0]);
				fillthis="DiffPosY_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[1]-posTrue[1]);
				fillthis="DiffPosZ_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[2]-posTrue[2]);
				
				return 1;
				
			}
			else{
				if(bfromB==kTRUE) fNentries->Fill(70);
				if(bfromB==kFALSE) fNentries->Fill(102);
				return 0;
			}
		}
	}
	else if(sD0c==3 && (HPiAODtrk->Charge()<0)){
		fMCTruth->Fill(3);
		Int_t pdgDaughtersD0[2]={321,211};
		Int_t pdgsBm=0;
		Int_t pdgsMPi=0;
		Int_t pdgsPi=0;
		Int_t labDpD0 = d->MatchToMC(421,fMcArray,2,pdgDaughtersD0);
		if(labDpD0>0){
			if(fDebug>2)   printf("Identified D0 Meson; could have been either \n");
			fMCTruth->Fill(6);
			AliAODMCParticle *partD0= (AliAODMCParticle*)fMcArray->At(labDpD0);
			if(!partD0) return 0;
			if(fDebug>2) printf("D0 candidate found at label: %i \n",partD0->GetLabel());
			// Get mother of the D0	from MC stack
			Int_t labMD0 = partD0->GetMother();
			AliAODMCParticle *partMD0= (AliAODMCParticle*)fMcArray->At(labMD0);
			if(!partMD0) return 0;
			
			Bool_t bfromB = kFALSE;
			if(fbHIJINGavailable){bfromB = IsCandidateInjected(d,fmcHeader);}
			
			if(fDebug>2)   printf("D0 Mother candidate found at label: %i \n",partMD0->GetLabel());
			// Check if mother is B-; otherwise --> reject
			pdgsBm=partMD0->GetPdgCode();
			if(pdgsBm!=-521){
				if(bfromB==kTRUE) fNentries->Fill(56);
				if(bfromB==kFALSE) fNentries->Fill(88);
				return 0;
			}
			// Check if B- has 2 daughters; otherwise --> reject
			if(partMD0->GetNDaughters()!=2){
				if(bfromB==kTRUE) fNentries->Fill(57);
				if(bfromB==kFALSE) fNentries->Fill(89);
				return 0;
			}
			if(fDebug>2)   printf("D0 Mother candidate is %i and decays into %i daughters \n",pdgsBm,partMD0->GetNDaughters());
			// Get charged track from MC stack
			AliAODMCParticle *partPi= (AliAODMCParticle*)fMcArray->At(TMath::Abs(HPiAODtrk->GetLabel()));
			if(!partPi){
				if(bfromB==kTRUE) fNentries->Fill(58);
				if(bfromB==kFALSE) fNentries->Fill(90);
				return 0;
			}
			// Check if Pi- otherwise --> reject
			pdgsPi=partPi->GetPdgCode();
			if(pdgsPi!=-211){
				if(bfromB==kTRUE) fNentries->Fill(59);
				if(bfromB==kFALSE) fNentries->Fill(91);
				return 0;
			}
			if(fDebug>2)   printf("%i track found at label: %i \n",pdgsPi,partPi->GetLabel());
			// Get mother of the Pi- from MC stack
			Int_t labMPi = partPi->GetMother();
			//if(labMPi<=0) return 0;
			AliAODMCParticle *partMPi= (AliAODMCParticle*)fMcArray->At(labMPi);
			if(!partMPi){
				if(bfromB==kTRUE) fNentries->Fill(60);
				if(bfromB==kFALSE) fNentries->Fill(92);
				return 0;
			}
			// Check if mother is B-; otherwise --> reject
			pdgsMPi=partMPi->GetPdgCode();
			if(pdgsMPi!=-521){
				if(bfromB==kTRUE) fNentries->Fill(61);
				if(bfromB==kFALSE) fNentries->Fill(93);
				return 0;
			}
			if(fDebug>2)   printf("Pion mother candidate found at label: %i \n",partMPi->GetLabel());
			if(fDebug>2)   printf("Pion mother candidate is %i and decays into %i daughters \n",pdgsMPi,partMPi->GetNDaughters());
			if(partMD0->GetLabel()==partMPi->GetLabel()) {
				fMCTruth->Fill(8);
				printf("Complete B- --> D0 + Pi- --> Pi+ + K- detected!!! \n");
				if(bfromB==kFALSE){
					cout << "### Particle injected by hand!" << endl;
					fNentries->Fill(40);
				}
				else if(bfromB==kTRUE){
					cout << "### Particle from quark!" << endl;
					fNentries->Fill(39);
				}
				// Fill Plots for Secondary Vertex Resolution
				// for B meson
				partPi->XvYvZv(posTrue);
				bMeson->GetSecondaryVtx(posReco);
				
				fillthis="DiffPosX_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[0]-posTrue[0]);
				fillthis="DiffPosY_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[1]-posTrue[1]);
				fillthis="DiffPosZ_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[2]-posTrue[2]);
				
				return 1;
				
			}
			else{
				if(bfromB==kTRUE) fNentries->Fill(71);
				if(bfromB==kFALSE) fNentries->Fill(103);
				return 0;
			}
		}
	}
	else if(sD0c==3 && (HPiAODtrk->Charge()>0)){
		fMCTruth->Fill(3);
		Int_t pdgDaughtersD0bar[2]={321,211};
		Int_t pdgsBp=0;
		Int_t pdgsMPi=0;
		Int_t pdgsPi=0;
		Int_t labDpD0bar = d->MatchToMC(421,fMcArray,2,pdgDaughtersD0bar);
		if(labDpD0bar>0){
			if(fDebug>2)   printf("Identified D0bar Meson; could have been either \n");
			fMCTruth->Fill(7);
			AliAODMCParticle *partD0bar= (AliAODMCParticle*)fMcArray->At(labDpD0bar);
			if(!partD0bar) return 0;
			if(fDebug>2)printf("D0bar candidate found at label: %i \n",partD0bar->GetLabel());
			// Get mother of the D0	from MC stack
			Int_t labMD0bar = partD0bar->GetMother();
			//if(labMD0bar<=0) return 0;
			AliAODMCParticle *partMD0bar= (AliAODMCParticle*)fMcArray->At(labMD0bar);
			if(!partMD0bar) return 0;
			
			Bool_t bfromB = kFALSE;
			if(fbHIJINGavailable){bfromB = IsCandidateInjected(d,fmcHeader);}
			
			if(fDebug>2)   printf("D0bar Mother candidate found at label: %i \n",partMD0bar->GetLabel());
			// Check if mother is B-; otherwise --> reject
			pdgsBp=partMD0bar->GetPdgCode();
			if(pdgsBp!=521){
				if(bfromB==kTRUE) fNentries->Fill(63);
				if(bfromB==kFALSE) fNentries->Fill(95);
				return 0;
			}
			// Check if B- has 2 daughters; otherwise --> reject
			if(partMD0bar->GetNDaughters()!=2){
				if(bfromB==kTRUE) fNentries->Fill(64);
				if(bfromB==kFALSE) fNentries->Fill(96);
				return 0;
			}
			if(fDebug>2)   printf("D0bar Mother candidate is %i and decays into %i daughters \n",pdgsBp,partMD0bar->GetNDaughters());
			// Get charged track from MC stack
			AliAODMCParticle *partPi= (AliAODMCParticle*)fMcArray->At(TMath::Abs(HPiAODtrk->GetLabel()));
			if(!partPi){
				if(bfromB==kTRUE) fNentries->Fill(65);
				if(bfromB==kFALSE) fNentries->Fill(97);
				return 0;
			}
			// Check if Pi+ otherwise --> reject
			pdgsPi=partPi->GetPdgCode();
			if(pdgsPi!=211){
				if(bfromB==kTRUE) fNentries->Fill(66);
				if(bfromB==kFALSE) fNentries->Fill(98);
				return 0;
			}
			if(fDebug>2)   printf("%i track found at label: %i \n",pdgsPi,partPi->GetLabel());
			// Get mother of the Pi+ from MC stack
			Int_t labMPi = partPi->GetMother();
			//if(labMPi<=0) return 0;
			AliAODMCParticle *partMPi= (AliAODMCParticle*)fMcArray->At(labMPi);
			if(!partMPi){
				if(bfromB==kTRUE) fNentries->Fill(67);
				if(bfromB==kFALSE) fNentries->Fill(99);
				return 0;
			}
			// Check if mother is B+; otherwise --> reject
			pdgsMPi=partMPi->GetPdgCode();
			if(pdgsMPi!=521){
				if(bfromB==kTRUE) fNentries->Fill(68);
				if(bfromB==kFALSE) fNentries->Fill(100);
				return 0;
			}
			if(fDebug>2)   printf("Pion mother candidate found at label: %i \n",partMPi->GetLabel());
			if(fDebug>2)   printf("Pion mother candidate is %i and decays into %i daughters \n",pdgsMPi,partMPi->GetNDaughters());
			if(partMD0bar->GetLabel()==partMPi->GetLabel()) {
				fMCTruth->Fill(9);
				printf("Complete B+ --> D0bar + Pi+ --> K+ + Pi- detected!!! \n");
				if(bfromB==kFALSE){
					cout << "### Particle injected by hand!" << endl;
					fNentries->Fill(40);
				}
				else if(bfromB==kTRUE){
					cout << "### Particle from quark!" << endl;
					fNentries->Fill(39);
				}
				// Fill Plots for Secondary Vertex Resolution
				// for B meson
				partPi->XvYvZv(posTrue);
				bMeson->GetSecondaryVtx(posReco);
				
				fillthis="DiffPosX_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[0]-posTrue[0]);
				fillthis="DiffPosY_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[1]-posTrue[1]);
				fillthis="DiffPosZ_";
				fillthis+=nbins;
				((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(posReco[2]-posTrue[2]);
				
				return 1;
				
			}
			else{
				if(bfromB==kTRUE) fNentries->Fill(72);
				if(bfromB==kFALSE) fNentries->Fill(104);
				return 0;
			}
		}
	}
	return 0;
}
//____________________________________________________________________________
void AliAnalysisTaskSEBpmMass::FillMCtruthDHists(AliAODRecoDecayHF2Prong *Bmeson, Double_t sD0c, AliAODTrack* HPiAODtrk, AliAODRecoDecayHF2Prong *Dmeson){
	
	if(fDebug>2) cout << "Filling MC Truth D Histos" << endl;
	
	// For safety, should not be possible at this stage
	Double_t massTrueD0 = 1.8648;
	if(sD0c==0) return;
	else if(sD0c==1){
		if(TMath::Abs(Dmeson->InvMassD0() - massTrueD0) > massWindowD){
			fNentries->Fill(115);
			return;
		}
	}
	else if(sD0c==2){
		if(TMath::Abs(Dmeson->InvMassD0bar() - massTrueD0) > massWindowD){
			fNentries->Fill(116);
			return;
		}
	}
	else if(sD0c==3){
		if((TMath::Abs(Dmeson->InvMassD0() - massTrueD0) > massWindowD) && (TMath::Abs(Dmeson->InvMassD0bar() - massTrueD0) > massWindowD)){
			fNentries->Fill(117);
			return;
		}
	};
	
	Double_t cosThStaD		= Dmeson->CosThetaStarD0();
	Double_t cosThStaDbar	= Dmeson->CosThetaStarD0bar();
	Double_t cosPtAngD		= Dmeson->CosPointingAngle();
	Double_t cosPtAngDXY	= Dmeson->CosPointingAngleXY();
	Double_t transMomD		= Dmeson->Pt();
	Double_t DCAD			= Dmeson->GetDCA();
	Double_t PtPiD			= Dmeson->PtProng(0);
	Double_t PtKaD			= Dmeson->PtProng(1);
	Double_t PPiD			= Dmeson->PProng(0);
	Double_t PKaD			= Dmeson->PProng(1);
	Double_t d0PiD			= Dmeson->Getd0Prong(0);
	Double_t d0KaD			= Dmeson->Getd0Prong(1);
	Double_t PtPiDbar		= Dmeson->PtProng(1);
	Double_t PtKaDbar		= Dmeson->PtProng(0);
	Double_t PPiDbar		= Dmeson->PProng(1);
	Double_t PKaDbar		= Dmeson->PProng(0);
	Double_t d0PiDbar		= Dmeson->Getd0Prong(1);
	Double_t d0KaDbar		= Dmeson->Getd0Prong(0);
	Double_t Prodd0d0D		= Dmeson->Prodd0d0();
	Double_t LthD			= Dmeson->DecayLength();
	Double_t normDecLthD	= Dmeson->NormalizedDecayLength();
	Double_t LthXYD			= Dmeson->DecayLengthXY();
	Double_t normDecLthXYD	= Dmeson->NormalizedDecayLengthXY();
	Double_t impparXYD		= Dmeson->ImpParXY()*1000.;
	Double_t ctauD			= Dmeson->Ct(421);
	
	Double_t massTrueB = 5.279;
	UInt_t pdgB[2]={0,0};
	pdgB[0] = 211;
	pdgB[1] = 421;
	if(TMath::Abs(Bmeson->InvMass(2,pdgB)-massTrueB)>massWindowB){
		fNentries->Fill(119);
		return;
	}
	
	// D0
	if(sD0c==1 || (sD0c==3 && HPiAODtrk->Charge()<0)){
		
		UInt_t pdgD[2]={0,0};
		pdgD[0] = 211;
		pdgD[1] = 321;
		Double_t massCandD = Dmeson->InvMass(2,pdgD);
		
		TString fillthis="";
		//Int_t nbins=fCuts->PtBin(Bmeson->Pt());
		// changed to pt binning of D for analysis of RMS
		Int_t nbinsD=fCuts->PtBin(Dmeson->Pt());
		
		fillthis="MCInvMassD";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(massCandD);
		fillthis="MCDCosThetaStar";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosThStaD);
		fillthis="MCDCosPoinAngle";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="MCDCosPoinAngleXY";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="MCDTransvMome";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(transMomD);
		fillthis="MCDDCA";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(DCAD);
		fillthis="MCDTransMomKa";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtKaD);
		fillthis="MCDTransMomPi";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtPiD);
		fillthis="MCDd0Ka";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0KaD);
		fillthis="MCDd0Pi";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0PiD);
		fillthis="MCDProdd0d0";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="MCDdeclgt";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthD);
		fillthis="MCDnormdeclgt";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="MCDdeclgtxy";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="MCDnormdeclgtXY";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthXYD);
		fillthis="MCDImpParXY";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="MCDctau";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(ctauD);
		
		fillthis="InvMassD_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(massCandD);
		fillthis="DCosThetaStar_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosThStaD);
		fillthis="DCosPointAngle_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="DCosPointAngleXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="D0pt_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(transMomD);
		fillthis="D0dca_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(DCAD);
		fillthis="DKapt_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtKaD);
		fillthis="DPipt_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtPiD);
		fillthis="DKap_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PKaD);
		fillthis="DPip_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PPiD);
		fillthis="DKad0_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0KaD);
		fillthis="DPid0_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0PiD);
		fillthis="DProdd0d0_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="DBimpParXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="Dctau_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(ctauD);
		fillthis="DDecayLength_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthD);
		fillthis="DnormDecayLength_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="DDecayLengthXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="DnormDecayLengthXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthXYD);
		
	}
	// D0bar
	else if(sD0c==2 || (sD0c==3 && HPiAODtrk->Charge()>0)){
		
		UInt_t pdgD[2]={0,0};
		pdgD[0] = 321;
		pdgD[1] = 211;
		Double_t massCandD = Dmeson->InvMass(2,pdgD);
		
		TString fillthis="";
		//Int_t nbins=fCuts->PtBin(Bmeson->Pt());
		// changed to pt binning of D for analysis of RMS
		Int_t nbinsD=fCuts->PtBin(Dmeson->Pt());
		
		fillthis="MCInvMassDbar";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(massCandD);
		fillthis="MCDCosThetaStar";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosThStaDbar);
		fillthis="MCDCosPoinAngle";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="MCDCosPoinAngleXY";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="MCDTransvMome";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(transMomD);
		fillthis="MCDDCA";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(DCAD);
		fillthis="MCDTransMomKa";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtKaDbar);
		fillthis="MCDTransMomPi";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtPiDbar);
		fillthis="MCDd0Ka";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0KaDbar);
		fillthis="MCDd0Pi";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0PiDbar);
		fillthis="MCDProdd0d0";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="MCDdeclgt";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthD);
		fillthis="MCDnormdeclgt";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="MCDdeclgtxy";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="MCDnormdeclgtXY";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthXYD);
		fillthis="MCDImpParXY";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="MCDctau";
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(ctauD);
		
		fillthis="InvMassDbar_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(massCandD);
		fillthis="DCosThetaStar_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosThStaDbar);
		fillthis="DCosPointAngle_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="DCosPointAngleXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="D0pt_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(transMomD);
		fillthis="D0dca_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(DCAD);
		fillthis="DKapt_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtKaDbar);
		fillthis="DPipt_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PtPiDbar);
		fillthis="DKap_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PKaDbar);
		fillthis="DPip_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(PPiDbar);
		fillthis="DKad0_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0KaDbar);
		fillthis="DPid0_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(d0PiDbar);
		fillthis="DProdd0d0_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="DBimpParXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="Dctau_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(ctauD);
		fillthis="DDecayLength_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthD);
		fillthis="DnormDecayLength_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="DDecayLengthXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="DnormDecayLengthXY_";
		fillthis+=nbinsD;
		((TH1F*)(fOutputListDSonly->FindObject(fillthis)))->Fill(normDecLthXYD);
	}
	return;
}
//____________________________________________________________________________
void AliAnalysisTaskSEBpmMass::FillDHists(AliAODRecoDecayHF2Prong *Bmeson, Double_t sD0c, AliAODTrack* HPiAODtrk, AliAODRecoDecayHF2Prong *Dmeson){
	
	if(fDebug>2) cout << "Filling MC Truth D Histos" << endl;
	
	// For safety, should not be possible at this stage
	Double_t massTrueD0 = 1.8648;
	if(sD0c==0) return;
	else if(sD0c==1){
		if(TMath::Abs(Dmeson->InvMassD0() - massTrueD0) > massWindowD){
			fNentries->Fill(112);
			return;
		}
	}
	else if(sD0c==2){
		if(TMath::Abs(Dmeson->InvMassD0bar() - massTrueD0) > massWindowD){
			fNentries->Fill(113);
			return;
		}
	}
	else if(sD0c==3){
		if((TMath::Abs(Dmeson->InvMassD0() - massTrueD0) > massWindowD) && (TMath::Abs(Dmeson->InvMassD0bar() - massTrueD0) > massWindowD)){
			fNentries->Fill(114);
			return;
		}
	}
	
	Double_t cosThStaD		= Dmeson->CosThetaStarD0();
	Double_t cosThStaDbar	= Dmeson->CosThetaStarD0bar();
	Double_t cosPtAngD		= Dmeson->CosPointingAngle();
	Double_t cosPtAngDXY	= Dmeson->CosPointingAngleXY();
	Double_t transMomD		= Dmeson->Pt();
	Double_t DCAD			= Dmeson->GetDCA();
	Double_t PtPiD			= Dmeson->PtProng(0);
	Double_t PtKaD			= Dmeson->PtProng(1);
	Double_t PPiD			= Dmeson->PProng(0);
	Double_t PKaD			= Dmeson->PProng(1);
	Double_t d0PiD			= Dmeson->Getd0Prong(0);
	Double_t d0KaD			= Dmeson->Getd0Prong(1);
	Double_t PtPiDbar		= Dmeson->PtProng(1);
	Double_t PtKaDbar		= Dmeson->PtProng(0);
	Double_t PPiDbar		= Dmeson->PProng(1);
	Double_t PKaDbar		= Dmeson->PProng(0);
	Double_t d0PiDbar		= Dmeson->Getd0Prong(1);
	Double_t d0KaDbar		= Dmeson->Getd0Prong(0);
	Double_t Prodd0d0D		= Dmeson->Prodd0d0();
	Double_t LthD			= Dmeson->DecayLength();
	Double_t normDecLthD	= Dmeson->NormalizedDecayLength();
	Double_t LthXYD			= Dmeson->DecayLengthXY();
	Double_t normDecLthXYD	= Dmeson->NormalizedDecayLengthXY();
	Double_t impparXYD		= Dmeson->ImpParXY()*1000.;
	Double_t ctauD			= Dmeson->Ct(421);
	
	Double_t massTrueB = 5.279;
	UInt_t pdgB[2]={0,0};
	pdgB[0] = 211;
	pdgB[1] = 421;
	if(TMath::Abs(Bmeson->InvMass(2,pdgB)-massTrueB)>massWindowB){
		fNentries->Fill(118);
		return;
	}
	// D0
	if(sD0c==1 || (sD0c==3 && HPiAODtrk->Charge()<0)){
		
		UInt_t pdgD[2]={0,0};
		pdgD[0] = 211;
		pdgD[1] = 321;
		Double_t massCandD = Dmeson->InvMass(2,pdgD);
		
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(Bmeson->Pt());
		// changed to pt binning of D for analysis of RMS
		//		Int_t nbins=fCuts->PtBin(Dmeson->Pt());
		
		fillthis="InvMassD";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(massCandD);
		fillthis="DCosThetaStar";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosThStaD);
		fillthis="DCosPoinAngle";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="DCosPoinAngleXY";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="DTransvMome";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(transMomD);
		fillthis="DDCA";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(DCAD);
		fillthis="DTransMomKa";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtKaD);
		fillthis="DTransMomPi";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtPiD);
		fillthis="Dd0Ka";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0KaD);
		fillthis="Dd0Pi";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0PiD);
		fillthis="DProdd0d0";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="Ddeclgt";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthD);
		fillthis="Dnormdeclgt";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="Ddeclgtxy";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="DnormdeclgtXY";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthXYD);
		fillthis="DImpParXY";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="Dctau";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(ctauD);
		
		
		fillthis="InvMassD_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(massCandD);
		fillthis="DCosThetaStar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosThStaD);
		fillthis="DCosPointAngle_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="DCosPointAngleXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="D0pt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(transMomD);
		fillthis="D0dca_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(DCAD);
		fillthis="DKapt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtKaD);
		fillthis="DPipt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtPiD);
		fillthis="DKap_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PKaD);
		fillthis="DPip_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PPiD);
		fillthis="DKad0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0KaD);
		fillthis="DPid0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0PiD);
		fillthis="DProdd0d0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="DBimpParXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="Dctau_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(ctauD);
		fillthis="DDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthD);
		fillthis="DnormDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="DDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="DnormDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthXYD);
		
	}
	// D0bar
	else if(sD0c==2 || (sD0c==3 && HPiAODtrk->Charge()>0)){
		
		UInt_t pdgD[2]={0,0};
		pdgD[0] = 321;
		pdgD[1] = 211;
		Double_t massCandD = Dmeson->InvMass(2,pdgD);
		
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(Bmeson->Pt());
		// changed to pt binning of D for analysis of RMS
		//Int_t nbins=fCuts->PtBin(Dmeson->Pt());
		
		fillthis="InvMassDbar";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(massCandD);
		fillthis="DCosThetaStar";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosThStaDbar);
		fillthis="DCosPoinAngle";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="DCosPoinAngleXY";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="DTransvMome";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(transMomD);
		fillthis="DDCA";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(DCAD);
		fillthis="DTransMomKa";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtKaDbar);
		fillthis="DTransMomPi";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtPiDbar);
		fillthis="Dd0Ka";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0KaDbar);
		fillthis="Dd0Pi";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0PiDbar);
		fillthis="DProdd0d0";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="Ddeclgt";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthD);
		fillthis="Dnormdeclgt";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="Ddeclgtxy";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="DnormdeclgtXY";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthXYD);
		fillthis="DImpParXY";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="Dctau";
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(ctauD);
		
		fillthis="InvMassDbar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(massCandD);
		fillthis="DCosThetaStar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosThStaDbar);
		fillthis="DCosPointAngle_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngD);
		fillthis="DCosPointAngleXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(cosPtAngDXY);
		fillthis="D0pt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(transMomD);
		fillthis="D0dca_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(DCAD);
		fillthis="DKapt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtKaDbar);
		fillthis="DPipt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PtPiDbar);
		fillthis="DKap_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PKaDbar);
		fillthis="DPip_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(PPiDbar);
		fillthis="DKad0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0KaDbar);
		fillthis="DPid0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(d0PiDbar);
		fillthis="DProdd0d0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(Prodd0d0D);
		fillthis="DBimpParXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(impparXYD);
		fillthis="Dctau_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(ctauD);
		fillthis="DDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthD);
		fillthis="DnormDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthD);
		fillthis="DDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(LthXYD);
		fillthis="DnormDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListDSandB->FindObject(fillthis)))->Fill(normDecLthXYD);
		
	}
	return;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::FillBpmHists(AliAODRecoDecayHF2Prong *part, Double_t sD0c, AliAODTrack* HPiAODtrk,AliAODRecoDecayHF2Prong* d){
	//
	// function to fill B Meson histograms:
	//
	if(fDebug>2) cout << "Filling Bpm Histos" << endl;
	
	fNentries->Fill(19);
	
	
	TVector3 pVectorB(part->Px(),part->Py(),part->Pz());
	TVector3 pVectorD(d->Px(),d->Py(),d->Pz());
	Double_t pAngle = pVectorD.Angle(pVectorB); // angle of D momentum wrt B momentum
	Double_t pCosine = TMath::Cos(pAngle);
	
	Double_t massTrueB = 5.279;//TDatabasePDG::Instance()->GetParticle(521)->Mass();
	Double_t massCandBp = 0;
	Double_t massCandBm = 0;
	Double_t massCandDiff = 0;
	Double_t cosThStaBp = 0;
	Double_t cosThStaBm = 0;
	Double_t cosPtAngBp = 0;
	Double_t cosPtAngBm = 0;
	Double_t cosPtAngBpXY = 0;
	Double_t cosPtAngBmXY = 0;
	Double_t transMomBp = 0;
	Double_t transMomBm = 0;
	Double_t DCABp = 0;
	Double_t DCABm = 0;
	Double_t PtDBp = 0;
	Double_t PtPiBp = 0;
	Double_t PtDBm = 0;
	Double_t PtPiBm = 0;
	Double_t d0DBp = 0;
	Double_t d0PiBp = 0;
	Double_t d0DBm = 0;
	Double_t d0PiBm = 0;
	Double_t Prodd0d0Bp = 0;
	Double_t Prodd0d0Bm = 0;
	Double_t normDecLthXYBp = 0;
	Double_t normDecLthXYBm = 0;
	Double_t impparXYBp = 0;
	Double_t impparXYBm = 0;
	Double_t ctauBp = 0;
	Double_t ctauBm = 0;
	
	UInt_t pdgBp[2]={0,0};
	UInt_t pdgBm[2]={0,0};
	
	Double_t massTrueD0 = 1.8648;
	if(sD0c==0) return kFALSE;            // should not be possible at this stage
	else if(sD0c==1){
		if(TMath::Abs(d->InvMassD0() - massTrueD0) > massWindowD){
			fNentries->Fill(106);
			return kFALSE;
		}
	}
	else if(sD0c==2){
		if(TMath::Abs(d->InvMassD0bar() - massTrueD0) > massWindowD){
			fNentries->Fill(107);
			return kFALSE;
		}
	}
	else if(sD0c==3){
		if((TMath::Abs(d->InvMassD0() - massTrueD0) > massWindowD) && (TMath::Abs(d->InvMassD0bar() - massTrueD0) > massWindowD)){
			fNentries->Fill(108);
			return kFALSE;
		}
	}
	
	if(fReadMC) {
		AliAODMCParticle *partPi= (AliAODMCParticle*)fMcArray->At(TMath::Abs(HPiAODtrk->GetLabel()));
		if(partPi){
			bPiPhysPrimary	= partPi->IsPhysicalPrimary();
			bPiWeak					= partPi->IsSecondaryFromWeakDecay();
			bPiPdgCode			= partPi->GetPdgCode();
//			bPiFastWeak			= IsFastMcFromWeakDecay(partPi);
//			if(partPi->IsSecondaryFromMaterial()) bPiFastWeak =0;
			//			bPiMaterial			= partPi->IsSecondaryFromMaterial();
		}
		AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
		AliAODMCParticle *DaughterPart0 = (AliAODMCParticle*)fMcArray->At(TMath::Abs(trk0->GetLabel()));
		if(DaughterPart0){
			bDDgh0PhysPrimary		= DaughterPart0->IsPhysicalPrimary();
			bDDgh0Weak					= DaughterPart0->IsSecondaryFromWeakDecay();
			bDDgh0PdgCode				= DaughterPart0->GetPdgCode();
			//			bDDgh0FastWeak			= IsFastMcFromWeakDecay(DaughterPart0);
//			if(DaughterPart0->IsSecondaryFromMaterial()) bDDgh0FastWeak =0;
	//			bDDgh0Material			= DaughterPart0->IsSecondaryFromMaterial();
		}
		
		AliAODTrack *trk1 = (AliAODTrack*)d->GetDaughter(1);
		AliAODMCParticle *DaughterPart1 = (AliAODMCParticle*)fMcArray->At(TMath::Abs(trk1->GetLabel()));
		if(DaughterPart1){
			bDDgh1PhysPrimary		= DaughterPart1->IsPhysicalPrimary();
			bDDgh1Weak					= DaughterPart1->IsSecondaryFromWeakDecay();
			bDDgh1PdgCode				= DaughterPart1->GetPdgCode();
//			bDDgh1FastWeak			= IsFastMcFromWeakDecay(DaughterPart1);
//			if(DaughterPart1->IsSecondaryFromMaterial()) bDDgh1FastWeak =0;
			//			bDDgh1Material			= DaughterPart1->IsSecondaryFromMaterial();
		}
	}
	
	if((HPiAODtrk->GetStatus()&AliESDtrack::kTPCpid ) && (HPiAODtrk->GetStatus()&AliESDtrack::kTPCin) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFpid) && (HPiAODtrk->GetStatus()&AliESDtrack::kTOFout) && (HPiAODtrk->GetStatus()&AliESDtrack::kTIME) && !(HPiAODtrk->GetStatus()&AliESDtrack::kTOFmismatch)) {
		bBPiTOFnSigma = fPIDResponse->NumberOfSigmasTOF(HPiAODtrk, AliPID::kPion);
		bBPiTPCnSigma = fPIDResponse->NumberOfSigmasTPC(HPiAODtrk, AliPID::kPion);
	}
	//
	
	//
	// new topomatic cut by Andrea Rossi
	//
	Double_t diffIP, errdiffIP;
	part->Getd0MeasMinusExpProng(0,fBzkG,diffIP,errdiffIP);
	bBNormIpPi = diffIP/errdiffIP;
	part->Getd0MeasMinusExpProng(1,fBzkG,diffIP,errdiffIP);
	bBNormIpD0 = diffIP/errdiffIP;
	//
	// ---
	//
	
	if(sD0c==1){              // B- --> D0 + Pi-
		pdgBm[0] = 211;
		pdgBm[1] = 421;
		massCandBm = part->InvMass(2,pdgBm);
		if((TMath::Abs(massCandBm - massTrueB)>massWindowB)) return kFALSE;
		if(fDebug>2) cout << "MASS TRUE B " << massTrueB << endl;
		if(fDebug>2) cout << "Invariant Mass B-: " << massCandBm << endl;
		
		if(fKillCandidates){
			// Advanced cut tuning applied here
			//KillCandidatesAdvanced(part,d,kFALSE);
		}
		
		fNentries->Fill(20);
		
		
		massCandDiff = massCandBm-d->InvMassD0();
		//		massCandDiff = CalculateInvMass(part,d->InvMassD0())-d->InvMassD0();
		cosThStaBm = part->CosThetaStar(0,521,211,421);
		cosPtAngBm = part->CosPointingAngleXY();
		cosPtAngBmXY = part->CosPointingAngleXY();
		transMomBm = part->Pt();
		DCABm = part->GetDCA();
		PtDBm = part->PtProng(1);
		PtPiBm = part->PtProng(0);
		d0DBm = part->Getd0Prong(1);
		d0PiBm = part->Getd0Prong(0);
		Prodd0d0Bm = part->Prodd0d0();
		normDecLthXYBm = part->NormalizedDecayLengthXY();
		impparXYBm = part->ImpParXY()*1000.;
		ctauBm = part->Ct(521);
		
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(part->Pt());
		
		fillthis="hdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="hnormdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="hdeclgtxy";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="hnormdeclgtXY";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="hdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
		fillthis="hnormdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
		
		fillthis="InvMassBplus_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandBm);
		fillthis="InvBMassDiff_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandDiff);
		fillthis="CosThetaStar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosThStaBm);
		fillthis="CosPointAngle_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBm);
		fillthis="CosPointAngleXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBmXY);
		fillthis="Bpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(transMomBm);
		fillthis="Bdca_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(DCABm);
		fillthis="Dpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtDBm);
		fillthis="Pipt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtPiBm);
		fillthis="Dd0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0DBm);
		fillthis="Pid0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0PiBm);
		fillthis="Prodd0d0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(Prodd0d0Bm);
		fillthis="BimpParXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(impparXYBm);
		fillthis="ctau_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(ctauBm);
		fillthis="DecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="normDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="DecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="normDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="DLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
		fillthis="normDLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
		fillthis="LengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
		fillthis="NormLengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
		fillthis="MometumAngleDwrtB_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(pCosine);
		fillthis="normDDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
		
		bBInvMass			= part->InvMass(2,pdgBm);
		bBInvMassDiff	= massCandDiff;
		bBpT					= part->Pt();
		bBctp					= part->CosPointingAngle();
		bBctpXY				= part->CosPointingAngleXY();
		bBd0XY				= TMath::Abs(part->ImpParXY());
		bBPipt				= part->PtProng(0);
		bd0d0					= part->Prodd0d0();
		bBndlXY				= part->NormalizedDecayLengthXY();
		bVtxChi2			= part->GetReducedChi2();
		bBdl					= part->DecayLength();
		bDpT					= part->PtProng(1);
		bDndl					= d->NormalizedDecayLength();
		bDdl					= d->DecayLength();
		bDInvMass			= d->InvMassD0();
		bDdca					= d->GetDCA();
		bDcts					= d->CosThetaStarD0();;
		bDKapt				= d->PtProng(1);
		bDPipt				= d->PtProng(0);
		bDd0Ka				= d->Getd0Prong(1);
		bDd0Pi				= d->Getd0Prong(0);
		bDd0d0				= d->Prodd0d0();
		bDctp					= d->CosPointingAngle();
		bDctpXY				= d->CosPointingAngleXY();
		bDndlXY				= d->NormalizedDecayLengthXY();
		
		AliAODTrack *daughterBPion	= (AliAODTrack*)part->GetDaughter(0);	// Pi
		AliAODTrack *daughterD0			= (AliAODTrack*)d->GetDaughter(1);		// Ka
		AliAODTrack *daughterD1			= (AliAODTrack*)d->GetDaughter(0);		// Pi
		Double_t d0BPion[2],covd0BPion[3];
		Double_t d0BD0[2],covd0BD0[3];
		Double_t d0Ddgh0[2],covd0Ddgh0[3];
		Double_t d0Ddgh1[2],covd0Ddgh1[3];

		daughterBPion->PropagateToDCA(fvtx1,fBzkG,100.,d0BPion,covd0BPion);
		d->PropagateToDCA(fvtx1,fBzkG,100.,d0BD0,covd0BD0);
		daughterD0->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh0,covd0Ddgh0);
		daughterD1->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh1,covd0Ddgh1);
		

		 dcaZBPi		= d0BPion[1];
		 dcaZD0			= d0BD0[1];
		 dcaZD0Ka		= d0Ddgh0[1];
		 dcaZD0Pi		= d0Ddgh1[1];
		
		
		fInvMassBpm->Fill(massCandBm);
		fInvMassDiff->Fill(massCandBm-(d->InvMassD0()));
		fCosThetaStar->Fill(cosThStaBm);
		fCosPoinAngle->Fill(cosPtAngBm);
		fCosPoinAngleXY->Fill(cosPtAngBmXY);
		fTransvMome->Fill(transMomBm);
		fDCA->Fill(DCABm);
		fTransMomD->Fill(PtDBm);
		fTransMomPi->Fill(PtPiBm);
		fd0D->Fill(d0DBm);
		fd0Pi->Fill(d0PiBm);
		fProdd0d0->Fill(Prodd0d0Bm);
		fNormDecLthXY->Fill(normDecLthXYBm);
		fImpParXY->Fill(impparXYBm);
		fctau->Fill(ctauBm);
	}
	else if(sD0c==2) {                  // B+
		pdgBp[0] = 211;
		pdgBp[1] = 421;
		massCandBp = part->InvMass(2,pdgBp);
		if((TMath::Abs(massCandBp - massTrueB)>massWindowB)) return kFALSE;
		if(fDebug>2) cout << "MASS TRUE B " << massTrueB << endl;
		if(fDebug>2) cout << "Invariant Mass B+: " << massCandBp << endl;
		
		if(fKillCandidates){
			// Advanced cut tuning applied here
			//KillCandidatesAdvanced(part,d,kFALSE);
		}
		
		fNentries->Fill(21);
		
		massCandDiff = massCandBp-d->InvMassD0bar();
		//		massCandDiff = CalculateInvMass(part,d->InvMassD0bar())-d->InvMassD0bar();
		cosThStaBp = part->CosThetaStar(0,521,211,421);
		cosPtAngBp = part->CosPointingAngleXY();
		cosPtAngBpXY = part->CosPointingAngleXY();
		transMomBp = part->Pt();
		DCABp = part->GetDCA();
		PtDBp = part->PtProng(1);
		PtPiBp = part->PtProng(0);
		d0DBp = part->Getd0Prong(1);
		d0PiBp = part->Getd0Prong(0);
		Prodd0d0Bp = part->Prodd0d0();
		normDecLthXYBp = part->NormalizedDecayLengthXY();
		impparXYBp = part->ImpParXY()*1000.;
		ctauBp = part->Ct(521);
		
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(part->Pt());
		
		fillthis="hdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="hnormdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="hdeclgtxy";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="hnormdeclgtXY";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="hdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
		fillthis="hnormdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
		
		fillthis="InvMassBplus_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandBp);
		fillthis="InvBMassDiff_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandDiff);
		fillthis="CosThetaStar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosThStaBp);
		fillthis="CosPointAngle_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBp);
		fillthis="CosPointAngleXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBpXY);
		fillthis="Bpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(transMomBp);
		fillthis="Bdca_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(DCABp);
		fillthis="Dpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtDBp);
		fillthis="Pipt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtPiBp);
		fillthis="Dd0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0DBp);
		fillthis="Pid0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0PiBp);
		fillthis="Prodd0d0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(Prodd0d0Bp);
		fillthis="BimpParXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(impparXYBp);
		fillthis="ctau_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(ctauBp);
		fillthis="DecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="normDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="DecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="normDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="DLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
		fillthis="normDLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
		fillthis="LengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
		fillthis="NormLengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
		fillthis="MometumAngleDwrtB_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(pCosine);
		fillthis="normDDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
		
		bBInvMass			= part->InvMass(2,pdgBp);
		bBInvMassDiff	= massCandDiff;
		bBpT					= part->Pt();
		bBctp					= part->CosPointingAngle();
		bBctpXY				= part->CosPointingAngleXY();
		bBd0XY				= TMath::Abs(part->ImpParXY());
		bBPipt				= part->PtProng(0);
		bd0d0					= part->Prodd0d0();
		bBndlXY				= part->NormalizedDecayLengthXY();
		bVtxChi2			= part->GetReducedChi2();
		bBdl					= part->DecayLength();
		bDpT					= part->PtProng(1);
		bDndl					= d->NormalizedDecayLength();
		bDdl					= d->DecayLength();
		bDInvMass			= d->InvMassD0bar();
		bDdca					= d->GetDCA();
		bDcts					= d->CosThetaStarD0bar();;
		bDKapt				= d->PtProng(0);
		bDPipt				= d->PtProng(1);
		bDd0Ka				= d->Getd0Prong(0);
		bDd0Pi				= d->Getd0Prong(1);
		bDd0d0				= d->Prodd0d0();
		bDctp					= d->CosPointingAngle();
		bDctpXY				= d->CosPointingAngleXY();
		bDndlXY				= d->NormalizedDecayLengthXY();
		
		AliAODTrack *daughterBPion	= (AliAODTrack*)part->GetDaughter(0);	// Pi
		AliAODTrack *daughterD0		= (AliAODTrack*)d->GetDaughter(0);		// Ka
		AliAODTrack *daughterD1		= (AliAODTrack*)d->GetDaughter(1);		// Pi
		Double_t d0BPion[2],covd0BPion[3];
		Double_t d0BD0[2],covd0BD0[3];
		Double_t d0Ddgh0[2],covd0Ddgh0[3];
		Double_t d0Ddgh1[2],covd0Ddgh1[3];

		daughterBPion->PropagateToDCA(fvtx1,fBzkG,100.,d0BPion,covd0BPion);
		d->PropagateToDCA(fvtx1,fBzkG,100.,d0BD0,covd0BD0);
		daughterD0->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh0,covd0Ddgh0);
		daughterD1->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh1,covd0Ddgh1);
		
		
		dcaZBPi			= d0BPion[1];
		dcaZD0			= d0BD0[1];
		dcaZD0Ka		= d0Ddgh0[1];
		dcaZD0Pi		= d0Ddgh1[1];
		
		
		fInvMassBpm->Fill(massCandBp);
		fInvMassDiff->Fill(massCandBp-(d->InvMassD0bar()));
		fCosThetaStar->Fill(cosThStaBp);
		fCosPoinAngle->Fill(cosPtAngBp);
		fCosPoinAngleXY->Fill(cosPtAngBpXY);
		fTransvMome->Fill(transMomBp);
		fDCA->Fill(DCABp);
		fTransMomD->Fill(PtDBp);
		fTransMomPi->Fill(PtPiBp);
		fd0D->Fill(d0DBp);
		fd0Pi->Fill(d0PiBp);
		fProdd0d0->Fill(Prodd0d0Bp);
		fNormDecLthXY->Fill(normDecLthXYBp);
		fImpParXY->Fill(impparXYBp);
		fctau->Fill(ctauBp);
	}
	else if(sD0c==3) {
		pdgBm[0] = 211; // Bm
		pdgBm[1] = 421;
		pdgBp[0] = 211; // Bp
		pdgBp[1] = 421;
		massCandBm = part->InvMass(2,pdgBm);
		massCandBp = part->InvMass(2,pdgBp);
		if((TMath::Abs(massCandBp - massTrueB)>massWindowB) && (TMath::Abs(massCandBm - massTrueB)>massWindowB)) return kFALSE;
		if((HPiAODtrk->Charge() < 0) && (TMath::Abs(massCandBm - massTrueB)<massWindowB)){
			if(fDebug>2) cout << "Invariant Mass B-: " << massCandBm << endl;
			
			if(fKillCandidates){
				// Advanced cut tuning applied here
				//KillCandidatesAdvanced(part,d,kFALSE);
			}
			
			massCandDiff = massCandBm-d->InvMassD0();
			//			massCandDiff = CalculateInvMass(part,d->InvMassD0())-d->InvMassD0();
			cosThStaBm = part->CosThetaStar(0,521,211,421);
			cosPtAngBm = part->CosPointingAngleXY();
			cosPtAngBmXY = part->CosPointingAngleXY();
			transMomBm = part->Pt();
			DCABm = part->GetDCA();
			PtDBm = part->PtProng(1);
			PtPiBm = part->PtProng(0);
			d0DBm = part->Getd0Prong(1);
			d0PiBm = part->Getd0Prong(0);
			Prodd0d0Bm = part->Prodd0d0();
			normDecLthXYBm = part->NormalizedDecayLengthXY();
			impparXYBm = part->ImpParXY()*1000.;
			ctauBm = part->Ct(521);
			
			TString fillthis="";
			Int_t nbins=fCuts->PtBin(part->Pt());
			
			fillthis="hdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="hnormdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="hdeclgtxy";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="hnormdeclgtXY";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="hdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
			fillthis="hnormdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
			
			fillthis="InvMassBplus_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandDiff);
			fillthis="InvBMassDiff_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandBm);
			fillthis="CosThetaStar_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosThStaBm);
			fillthis="CosPointAngle_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBm);
			fillthis="CosPointAngleXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBmXY);
			fillthis="Bpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(transMomBm);
			fillthis="Bdca_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(DCABm);
			fillthis="Dpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtDBm);
			fillthis="Pipt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtPiBm);
			fillthis="Dd0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0DBm);
			fillthis="Pid0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0PiBm);
			fillthis="Prodd0d0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(Prodd0d0Bm);
			fillthis="BimpParXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(impparXYBm);
			fillthis="ctau_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(ctauBm);
			fillthis="DecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="normDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="DecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="normDecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="DLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
			fillthis="normDLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
			fillthis="LengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
			fillthis="NormLengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
			fillthis="MometumAngleDwrtB_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(pCosine);
			fillthis="normDDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
			
			bBInvMass			= part->InvMass(2,pdgBm);
			bBInvMassDiff	= massCandDiff;
			bBpT					= part->Pt();
			bBctp					= part->CosPointingAngle();
			bBctpXY				= part->CosPointingAngleXY();
			bBd0XY				= TMath::Abs(part->ImpParXY());
			bBPipt				= part->PtProng(0);
			bd0d0					= part->Prodd0d0();
			bBndlXY				= part->NormalizedDecayLengthXY();
			bVtxChi2			= part->GetReducedChi2();
			bBdl					= part->DecayLength();
			bDpT					= part->PtProng(1);
			bDndl					= d->NormalizedDecayLength();
			bDdl					= d->DecayLength();
			bDInvMass			= d->InvMassD0();
			bDdca					= d->GetDCA();
			bDcts					= d->CosThetaStarD0();;
			bDKapt				= d->PtProng(1);
			bDPipt				= d->PtProng(0);
			bDd0Ka				= d->Getd0Prong(1);
			bDd0Pi				= d->Getd0Prong(0);
			bDd0d0					= d->Prodd0d0();
			bDctp					= d->CosPointingAngle();
			bDctpXY				= d->CosPointingAngleXY();
			bDndlXY				= d->NormalizedDecayLengthXY();
			
			AliAODTrack *daughterBPion	= (AliAODTrack*)part->GetDaughter(0);	// Pi
			AliAODTrack *daughterD0		= (AliAODTrack*)d->GetDaughter(1);		// Ka
			AliAODTrack *daughterD1		= (AliAODTrack*)d->GetDaughter(0);		// Pi
			Double_t d0BPion[2],covd0BPion[3];
			Double_t d0BD0[2],covd0BD0[3];
			Double_t d0Ddgh0[2],covd0Ddgh0[3];
			Double_t d0Ddgh1[2],covd0Ddgh1[3];

			daughterBPion->PropagateToDCA(fvtx1,fBzkG,100.,d0BPion,covd0BPion);
			d->PropagateToDCA(fvtx1,fBzkG,100.,d0BD0,covd0BD0);
			daughterD0->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh0,covd0Ddgh0);
			daughterD1->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh1,covd0Ddgh1);
			
			
		 dcaZBPi		= d0BPion[1];
		 dcaZD0			= d0BD0[1];
		 dcaZD0Ka		= d0Ddgh0[1];
		 dcaZD0Pi		= d0Ddgh1[1];
			
			
			fInvMassBpm->Fill(massCandBm);
			fInvMassDiff->Fill(massCandBm-(d->InvMassD0()));
			fCosThetaStar->Fill(cosThStaBm);
			fCosPoinAngle->Fill(cosPtAngBm);
			fCosPoinAngleXY->Fill(cosPtAngBmXY);
			fTransvMome->Fill(transMomBm);
			fDCA->Fill(DCABm);
			fTransMomD->Fill(PtDBm);
			fTransMomPi->Fill(PtPiBm);
			fd0D->Fill(d0DBm);
			fd0Pi->Fill(d0PiBm);
			fProdd0d0->Fill(Prodd0d0Bm);
			fNormDecLthXY->Fill(normDecLthXYBm);
			fImpParXY->Fill(impparXYBm);
			fctau->Fill(ctauBm);
			fNentries->Fill(22);
		}
		else if((HPiAODtrk->Charge() > 0) && (TMath::Abs(massCandBp - massTrueB)<massWindowB)){
			if(fDebug>2) cout << "Invariant Mass B+: " << massCandBp << endl;
			
			if(fKillCandidates){
				// Advanced cut tuning applied here
				//KillCandidatesAdvanced(part,d,kFALSE);
			}
			
			massCandDiff = massCandBp-d->InvMassD0bar();
			//		massCandDiff = CalculateInvMass(part,d->InvMassD0bar())-d->InvMassD0bar();
			cosThStaBp = part->CosThetaStar(0,521,211,421);
			cosPtAngBp = part->CosPointingAngleXY();
			cosPtAngBpXY = part->CosPointingAngleXY();
			transMomBp = part->Pt();
			DCABp = part->GetDCA();
			PtDBp = part->PtProng(1);
			PtPiBp = part->PtProng(0);
			d0DBp = part->Getd0Prong(1);
			d0PiBp = part->Getd0Prong(0);
			Prodd0d0Bp = part->Prodd0d0();
			normDecLthXYBp = part->NormalizedDecayLengthXY();
			impparXYBp = part->ImpParXY()*1000.;
			ctauBp = part->Ct(521);
			
			TString fillthis="";
			Int_t nbins=fCuts->PtBin(part->Pt());
			
			fillthis="hdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="hnormdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="hdeclgtxy";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="hnormdeclgtXY";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="hdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
			fillthis="hnormdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
			
			fillthis="InvMassBplus_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandBp);
			fillthis="InvBMassDiff_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(massCandDiff);
			fillthis="CosThetaStar_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosThStaBp);
			fillthis="CosPointAngle_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBp);
			fillthis="CosPointAngleXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(cosPtAngBpXY);
			fillthis="Bpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(transMomBp);
			fillthis="Bdca_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(DCABp);
			fillthis="Dpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtDBp);
			fillthis="Pipt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(PtPiBp);
			fillthis="Dd0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0DBp);
			fillthis="Pid0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d0PiBp);
			fillthis="Prodd0d0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(Prodd0d0Bp);
			fillthis="BimpParXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(impparXYBp);
			fillthis="ctau_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(ctauBp);
			fillthis="DecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="normDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="DecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="normDecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="DLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
			fillthis="normDLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
			fillthis="LengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
			fillthis="NormLengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
			fillthis="MometumAngleDwrtB_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(pCosine);
			fillthis="normDDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
			
			bBInvMass			= part->InvMass(2,pdgBp);
			bBInvMassDiff	= massCandDiff;
			bBpT					= part->Pt();
			bBctp					= part->CosPointingAngle();
			bBctpXY				= part->CosPointingAngleXY();
			bBd0XY				= TMath::Abs(part->ImpParXY());
			bBPipt				= part->PtProng(0);
			bd0d0					= part->Prodd0d0();
			bBndlXY				= part->NormalizedDecayLengthXY();
			bVtxChi2			= part->GetReducedChi2();
			bBdl					= part->DecayLength();
			bDpT					= part->PtProng(1);
			bDndl					= d->NormalizedDecayLength();
			bDdl					= d->DecayLength();
			bDInvMass			= d->InvMassD0bar();
			bDdca					= d->GetDCA();
			bDcts					= d->CosThetaStarD0bar();;
			bDKapt				= d->PtProng(0);
			bDPipt				= d->PtProng(1);
			bDd0Ka				= d->Getd0Prong(0);
			bDd0Pi				= d->Getd0Prong(1);
			bDd0d0					= d->Prodd0d0();
			bDctp					= d->CosPointingAngle();
			bDctpXY				= d->CosPointingAngleXY();
			bDndlXY				= d->NormalizedDecayLengthXY();
			
			AliAODTrack *daughterBPion	= (AliAODTrack*)part->GetDaughter(0);	// Pi
			AliAODTrack *daughterD0		= (AliAODTrack*)d->GetDaughter(0);		// Ka
			AliAODTrack *daughterD1		= (AliAODTrack*)d->GetDaughter(1);		// Pi
			Double_t d0BPion[2],covd0BPion[3];
			Double_t d0BD0[2],covd0BD0[3];
			Double_t d0Ddgh0[2],covd0Ddgh0[3];
			Double_t d0Ddgh1[2],covd0Ddgh1[3];

			daughterBPion->PropagateToDCA(fvtx1,fBzkG,100.,d0BPion,covd0BPion);
			d->PropagateToDCA(fvtx1,fBzkG,100.,d0BD0,covd0BD0);
			daughterD0->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh0,covd0Ddgh0);
			daughterD1->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh1,covd0Ddgh1);
			
			
		 dcaZBPi			= d0BPion[1];
		 dcaZD0			= d0BD0[1];
		 dcaZD0Ka		= d0Ddgh0[1];
		 dcaZD0Pi		= d0Ddgh1[1];
			
			
			fInvMassBpm->Fill(massCandBp);
			fInvMassDiff->Fill(massCandBp-(d->InvMassD0bar()));
			fCosThetaStar->Fill(cosThStaBp);
			fCosPoinAngle->Fill(cosPtAngBp);
			fCosPoinAngleXY->Fill(cosPtAngBpXY);
			fTransvMome->Fill(transMomBp);
			fDCA->Fill(DCABp);
			fTransMomD->Fill(PtDBp);
			fTransMomPi->Fill(PtPiBp);
			fd0D->Fill(d0DBp);
			fd0Pi->Fill(d0PiBp);
			fProdd0d0->Fill(Prodd0d0Bp);
			fNormDecLthXY->Fill(normDecLthXYBp);
			fImpParXY->Fill(impparXYBp);
			fNentries->Fill(23);
		}
		else return kFALSE;
	}
	fNentries->Fill(24);
	fSelectionVariables->Fill();
	return kTRUE;
}
//____________________________________________________________________________
void AliAnalysisTaskSEBpmMass::FillBpmMCHists(AliAODRecoDecayHF2Prong *part, Double_t sD0c, AliAODTrack* HPiAODtrk,AliAODRecoDecayHF2Prong* d){
	//
	// function used to fill B Meson histograms:
	//
	if(fDebug>2) cout << "Filling BpmMC Histos" << endl;
	
	fMCTruth->Fill(11);
	
	
	TVector3 pVectorB(part->Px(),part->Py(),part->Pz());
	TVector3 pVectorD(d->Px(),d->Py(),d->Pz());
	Double_t pAngle = pVectorD.Angle(pVectorB); // angle of D momentum wrt B momentum
	Double_t pCosine = TMath::Cos(pAngle);
	
	Double_t massTrueB = 5.279;//TDatabasePDG::Instance()->GetParticle(521)->Mass();
	Double_t massCandBp = 0;
	Double_t massCandBm = 0;
	Double_t massCandDiff = 0;
	Double_t cosThStaBp = 0;
	Double_t cosThStaBm = 0;
	Double_t cosPtAngBp = 0;
	Double_t cosPtAngBm = 0;
	Double_t cosPtAngBpXY = 0;
	Double_t cosPtAngBmXY = 0;
	Double_t transMomBp = 0;
	Double_t transMomBm = 0;
	Double_t DCABp = 0;
	Double_t DCABm = 0;
	Double_t PtDBp = 0;
	Double_t PtPiBp = 0;
	Double_t PtDBm = 0;
	Double_t PtPiBm = 0;
	Double_t d0DBp = 0;
	Double_t d0PiBp = 0;
	Double_t d0DBm = 0;
	Double_t d0PiBm = 0;
	Double_t Prodd0d0Bp = 0;
	Double_t Prodd0d0Bm = 0;
	Double_t normDecLthXYBp = 0;
	Double_t normDecLthXYBm = 0;
	Double_t impparXYBp = 0;
	Double_t impparXYBm = 0;
	Double_t ctauBp = 0;
	Double_t ctauBm = 0;
	
	UInt_t pdgBp[2]={0,0};
	UInt_t pdgBm[2]={0,0};
	
	Double_t massTrueD0 = 1.8648;
	if(sD0c==0) return;            // should not be possible at this stage
	else if(sD0c==1){
		if(TMath::Abs(d->InvMassD0() - massTrueD0) > massWindowD){
			fNentries->Fill(109);
			return;
		}
	}
	else if(sD0c==2){
		if(TMath::Abs(d->InvMassD0bar() - massTrueD0) > massWindowD){
			fNentries->Fill(110);
			return;
		}
	}
	else if(sD0c==3){
		if((TMath::Abs(d->InvMassD0() - massTrueD0) > massWindowD) && (TMath::Abs(d->InvMassD0bar() - massTrueD0) > massWindowD)){
			fNentries->Fill(111);
			return;
		}
	}
	
	//
	// new topomatic cut by Andrea Rossi
	//
	Double_t diffIP, errdiffIP;
	part->Getd0MeasMinusExpProng(0,fBzkG,diffIP,errdiffIP);
	bBNormIpPi = diffIP/errdiffIP;
	part->Getd0MeasMinusExpProng(1,fBzkG,diffIP,errdiffIP);
	bBNormIpD0 = diffIP/errdiffIP;
	//
	// ---
	//
	if(sD0c==1){                   // B-
		pdgBm[0] = 211;
		pdgBm[1] = 421;
		massCandBm = part->InvMass(2,pdgBm);
		fMCTruth->Fill(12);
		if((TMath::Abs(massCandBm - massTrueB)>massWindowB)) return;
		if(fDebug>2) cout << "MASS TRUE B " << massTrueB << endl;
		if(fDebug>2) cout << "Invariant Mass B-: " << massCandBm << endl;
		
		if(fKillCandidates){
			// Advanced cut tuning applied here
			//KillCandidatesAdvanced(part,d,kTRUE);
		}
		
		massCandDiff = massCandBm-d->InvMassD0();
		//	massCandDiff = CalculateInvMass(part,d->InvMassD0())-d->InvMassD0();
		cosThStaBm = part->CosThetaStar(0,521,211,421);
		cosPtAngBm = part->CosPointingAngleXY();
		cosPtAngBmXY = part->CosPointingAngleXY();
		transMomBm = part->Pt();
		DCABm = part->GetDCA();
		PtDBm = part->PtProng(1);
		PtPiBm = part->PtProng(0);
		d0DBm = part->Getd0Prong(1);
		d0PiBm = part->Getd0Prong(0);
		Prodd0d0Bm = part->Prodd0d0();
		normDecLthXYBm = part->NormalizedDecayLengthXY();
		impparXYBm = part->ImpParXY()*1000.;
		ctauBm = part->Ct(521);
		
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(part->Pt());
		
		fillthis="hMCdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="hMCnormdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="hMCdeclgtxy";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="hMCnormdeclgtXY";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="hMCdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
		fillthis="hMCnormdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
		
		fillthis="InvMassBplus_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandBm);
		fillthis="InvBMassDiff_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandDiff);
		fillthis="CosThetaStar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosThStaBm);
		fillthis="CosPointAngle_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBm);
		fillthis="CosPointAngleXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBmXY);
		fillthis="Bpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(transMomBm);
		fillthis="Bdca_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(DCABm);
		fillthis="Dpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtDBm);
		fillthis="Pipt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtPiBm);
		fillthis="Dd0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0DBm);
		fillthis="Pid0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0PiBm);
		fillthis="Prodd0d0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(Prodd0d0Bm);
		fillthis="BimpParXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(impparXYBm);
		fillthis="ctau_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(ctauBm);
		fillthis="DecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="normDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="DecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="normDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="DLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
		fillthis="normDLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
		fillthis="LengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
		fillthis="NormLengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
		fillthis="MometumAngleDwrtB_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(pCosine);
		fillthis="normDDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
		
		bBInvMass			= part->InvMass(2,pdgBm);
		bBInvMassDiff	= massCandDiff;
		bBpT					= part->Pt();
		bBctp					= part->CosPointingAngle();
		bBctpXY				= part->CosPointingAngleXY();
		bBd0XY				= TMath::Abs(part->ImpParXY());
		bBPipt				= part->PtProng(0);
		bd0d0					= part->Prodd0d0();
		bBndlXY				= part->NormalizedDecayLengthXY();
		bVtxChi2			= part->GetReducedChi2();
		bBdl					= part->DecayLength();
		bDpT					= part->PtProng(1);
		bDndl					= d->NormalizedDecayLength();
		bDdl					= d->DecayLength();
		bDInvMass			= d->InvMassD0();
		bDdca					= d->GetDCA();
		bDcts					= d->CosThetaStarD0();;
		bDKapt				= d->PtProng(1);
		bDPipt				= d->PtProng(0);
		bDd0Ka				= d->Getd0Prong(1);
		bDd0Pi				= d->Getd0Prong(0);
		bDd0d0					= d->Prodd0d0();
		bDctp					= d->CosPointingAngle();
		bDctpXY				= d->CosPointingAngleXY();
		bDndlXY				= d->NormalizedDecayLengthXY();
		
		fMCTruth->Fill(14);
		fMCInvMassBpm->Fill(massCandBm);
		fMCInvMassDiff->Fill(massCandBm-(d->InvMassD0()));
		fMCCosThetaStar->Fill(cosThStaBm);
		fMCCosPoinAngle->Fill(cosPtAngBm);
		fMCCosPoinAngleXY->Fill(cosPtAngBmXY);
		fMCTransvMome->Fill(transMomBm);
		fMCDCA->Fill(DCABm);
		fMCTransMomD->Fill(PtDBm);
		fMCTransMomPi->Fill(PtPiBm);
		fMCd0D->Fill(d0DBm);
		fMCd0Pi->Fill(d0PiBm);
		fMCProdd0d0->Fill(Prodd0d0Bm);
		fMCNormDecLthXY->Fill(normDecLthXYBm);
		fMCImpParXY->Fill(impparXYBm);
		fMCctau->Fill(ctauBm);
		hThetaD0->Fill(part->ThetaProng(1));
		hThetaPi->Fill(part->ThetaProng(0));
		hThetaD0Ka->Fill(d->ThetaProng(1));
		hThetaD0Pi->Fill(d->ThetaProng(0));
	}
	else if(sD0c==2) {                  // B+
		pdgBp[0] = 211;
		pdgBp[1] = 421;
		massCandBp = part->InvMass(2,pdgBp);
		fMCTruth->Fill(13);
		if((TMath::Abs(massCandBp - massTrueB)>massWindowB)) return;
		if(fDebug>2) cout << "MASS TRUE B " << massTrueB << endl;
		if(fDebug>2) cout << "Invariant Mass B+: " << massCandBp << endl;
		
		if(fKillCandidates){
			// Advanced cut tuning applied here
			//KillCandidatesAdvanced(part,d,kTRUE);
		}
		
		massCandDiff = massCandBp-d->InvMassD0bar();
		//	massCandDiff = CalculateInvMass(part,d->InvMassD0bar())-d->InvMassD0bar();
		cosThStaBp = part->CosThetaStar(0,521,211,421);
		cosPtAngBp = part->CosPointingAngleXY();
		cosPtAngBpXY = part->CosPointingAngleXY();
		transMomBp = part->Pt();
		DCABp = part->GetDCA();
		PtDBp = part->PtProng(1);
		PtPiBp = part->PtProng(0);
		d0DBp = part->Getd0Prong(1);
		d0PiBp = part->Getd0Prong(0);
		Prodd0d0Bp = part->Prodd0d0();
		normDecLthXYBp = part->NormalizedDecayLengthXY();
		impparXYBp = part->ImpParXY()*1000.;
		ctauBp = part->Ct(521);
		
		TString fillthis="";
		Int_t nbins=fCuts->PtBin(part->Pt());
		
		fillthis="hMCdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="hMCnormdeclgt";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="hMCdeclgtxy";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="hMCnormdeclgtXY";
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="hMCdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
		fillthis="hMCnormdeclxyd0d0";
		((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
		
		fillthis="InvMassBplus_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandBp);
		fillthis="InvBMassDiff_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandDiff);
		fillthis="CosThetaStar_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosThStaBp);
		fillthis="CosPointAngle_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBp);
		fillthis="CosPointAngleXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBpXY);
		fillthis="Bpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(transMomBp);
		fillthis="Bdca_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(DCABp);
		fillthis="Dpt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtDBp);
		fillthis="Pipt_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtPiBp);
		fillthis="Dd0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0DBp);
		fillthis="Pid0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0PiBp);
		fillthis="Prodd0d0_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(Prodd0d0Bp);
		fillthis="BimpParXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(impparXYBp);
		fillthis="ctau_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(ctauBp);
		fillthis="DecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLength());
		fillthis="normDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
		fillthis="DecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
		fillthis="normDecayLengthXY_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
		fillthis="DLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
		fillthis="normDLXYd0d0_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
		fillthis="LengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
		fillthis="NormLengthBvsD_";
		fillthis+=nbins;
		((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
		fillthis="MometumAngleDwrtB_";
		fillthis+=nbins;
		((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(pCosine);
		fillthis="normDDecayLength_";
		fillthis+=nbins;
		((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
		
		bBInvMass			= part->InvMass(2,pdgBp);
		bBInvMassDiff	= massCandDiff;
		bBpT					= part->Pt();
		bBctp					= part->CosPointingAngle();
		bBctpXY				= part->CosPointingAngleXY();
		bBd0XY				= TMath::Abs(part->ImpParXY());
		bBPipt				= part->PtProng(0);
		bd0d0					= part->Prodd0d0();
		bBndlXY				= part->NormalizedDecayLengthXY();
		bVtxChi2			= part->GetReducedChi2();
		bBdl					= part->DecayLength();
		bDpT					= part->PtProng(1);
		bDndl					= d->NormalizedDecayLength();
		bDdl					= d->DecayLength();
		bDInvMass			= d->InvMassD0bar();
		bDdca					= d->GetDCA();
		bDcts					= d->CosThetaStarD0bar();;
		bDKapt				= d->PtProng(0);
		bDPipt				= d->PtProng(1);
		bDd0Ka				= d->Getd0Prong(0);
		bDd0Pi				= d->Getd0Prong(1);
		bDd0d0					= d->Prodd0d0();
		bDctp					= d->CosPointingAngle();
		bDctpXY				= d->CosPointingAngleXY();
		bDndlXY				= d->NormalizedDecayLengthXY();
		
		fMCTruth->Fill(15);
		fMCInvMassBpm->Fill(massCandBp);
		fMCInvMassDiff->Fill(massCandBp-(d->InvMassD0bar()));
		fMCCosThetaStar->Fill(cosThStaBp);
		fMCCosPoinAngle->Fill(cosPtAngBp);
		fMCCosPoinAngleXY->Fill(cosPtAngBpXY);
		fMCTransvMome->Fill(transMomBp);
		fMCDCA->Fill(DCABp);
		fMCTransMomD->Fill(PtDBp);
		fMCTransMomPi->Fill(PtPiBp);
		fMCd0D->Fill(d0DBp);
		fMCd0Pi->Fill(d0PiBp);
		fMCProdd0d0->Fill(Prodd0d0Bp);
		fMCNormDecLthXY->Fill(normDecLthXYBp);
		fMCImpParXY->Fill(impparXYBp);
		fMCctau->Fill(ctauBp);
		hThetaD0->Fill(part->ThetaProng(1));
		hThetaPi->Fill(part->ThetaProng(0));
		hThetaD0Ka->Fill(d->ThetaProng(0));
		hThetaD0Pi->Fill(d->ThetaProng(1));
		
	}
	else if(sD0c==3) {
		fMCTruth->Fill(16);
		pdgBm[0] = 211; // Bm
		pdgBm[1] = 421;
		pdgBp[0] = 211; // Bp
		pdgBp[1] = 421;
		massCandBm = part->InvMass(2,pdgBm);
		massCandBp = part->InvMass(2,pdgBp);
		if((TMath::Abs(massCandBp - massTrueB)>massWindowB) && (TMath::Abs(massCandBm - massTrueB)>massWindowB)) return;
		if((HPiAODtrk->Charge() < 0) && (TMath::Abs(massCandBm - massTrueB)<massWindowB)){
			if(fDebug>2) cout << "Invariant Mass B-: " << massCandBm << endl;
			
			if(fKillCandidates){
				// Advanced cut tuning applied here
				//KillCandidatesAdvanced(part,d,kTRUE);
			}
			
			massCandDiff = massCandBm-d->InvMassD0();
			//		massCandDiff = CalculateInvMass(part,d->InvMassD0())-d->InvMassD0();
			cosThStaBm = part->CosThetaStar(0,521,211,421);
			cosPtAngBm = part->CosPointingAngleXY();
			cosPtAngBmXY = part->CosPointingAngleXY();
			transMomBm = part->Pt();
			DCABm = part->GetDCA();
			PtDBm = part->PtProng(1);
			PtPiBm = part->PtProng(0);
			d0DBm = part->Getd0Prong(1);
			d0PiBm = part->Getd0Prong(0);
			Prodd0d0Bm = part->Prodd0d0();
			normDecLthXYBm = part->NormalizedDecayLengthXY();
			impparXYBm = part->ImpParXY()*1000.;
			ctauBm = part->Ct(521);
			
			TString fillthis="";
			Int_t nbins=fCuts->PtBin(part->Pt());
			
			fillthis="hMCdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="hMCnormdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="hMCdeclgtxy";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="hMCnormdeclgtXY";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="hMCdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
			fillthis="hMCnormdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
			
			fillthis="InvMassBplus_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandBm);
			fillthis="InvBMassDiff_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandDiff);
			fillthis="CosThetaStar_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosThStaBm);
			fillthis="CosPointAngle_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBm);
			fillthis="CosPointAngleXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBmXY);
			fillthis="Bpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(transMomBm);
			fillthis="Bdca_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(DCABm);
			fillthis="Dpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtDBm);
			fillthis="Pipt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtPiBm);
			fillthis="Dd0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0DBm);
			fillthis="Pid0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0PiBm);
			fillthis="Prodd0d0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(Prodd0d0Bm);
			fillthis="BimpParXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(impparXYBm);
			fillthis="ctau_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(ctauBm);
			fillthis="DecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="normDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="DecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="normDecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="DLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bm);
			fillthis="normDLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(normDecLthXYBm,Prodd0d0Bm);
			fillthis="LengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
			fillthis="NormLengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
			fillthis="MometumAngleDwrtB_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(pCosine);
			fillthis="normDDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
			
			bBInvMass			= part->InvMass(2,pdgBm);
			bBInvMassDiff	= massCandDiff;
			bBpT					= part->Pt();
			bBctp					= part->CosPointingAngle();
			bBctpXY				= part->CosPointingAngleXY();
			bBd0XY				= TMath::Abs(part->ImpParXY());
			bBPipt				= part->PtProng(0);
			bd0d0					= part->Prodd0d0();
			bBndlXY				= part->NormalizedDecayLengthXY();
			bVtxChi2			= part->GetReducedChi2();
			bBdl					= part->DecayLength();
			bDpT					= part->PtProng(1);
			bDndl					= d->NormalizedDecayLength();
			bDdl					= d->DecayLength();
			bDInvMass			= d->InvMassD0();
			bDdca					= d->GetDCA();
			bDcts					= d->CosThetaStarD0();;
			bDKapt				= d->PtProng(1);
			bDPipt				= d->PtProng(0);
			bDd0Ka				= d->Getd0Prong(1);
			bDd0Pi				= d->Getd0Prong(0);
			bDd0d0					= d->Prodd0d0();
			bDctp					= d->CosPointingAngle();
			bDctpXY				= d->CosPointingAngleXY();
			bDndlXY				= d->NormalizedDecayLengthXY();
			
			AliAODTrack *daughterBPion	= (AliAODTrack*)part->GetDaughter(0);	// Pi
			AliAODTrack *daughterD0		= (AliAODTrack*)d->GetDaughter(1);		// Ka
			AliAODTrack *daughterD1		= (AliAODTrack*)d->GetDaughter(0);		// Pi
			Double_t d0BPion[2],covd0BPion[3];
			Double_t d0BD0[2],covd0BD0[3];
			Double_t d0Ddgh0[2],covd0Ddgh0[3];
			Double_t d0Ddgh1[2],covd0Ddgh1[3];

			daughterBPion->PropagateToDCA(fvtx1,fBzkG,100.,d0BPion,covd0BPion);
			d->PropagateToDCA(fvtx1,fBzkG,100.,d0BD0,covd0BD0);
			daughterD0->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh0,covd0Ddgh0);
			daughterD1->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh1,covd0Ddgh1);
			
			
		 dcaZBPi			= d0BPion[1];
		 dcaZD0			= d0BD0[1];
		 dcaZD0Ka		= d0Ddgh0[1];
		 dcaZD0Pi		= d0Ddgh1[1];
			
			
			/*
			 Double_t pv[3] = {0.,0.,0.};
			 fvtx1->GetXYZ(pv);
			 TVector3 fline(part->GetSecVtxX()-pv[0],
			 part->GetSecVtxY()-pv[1],
			 part->GetSecVtxZ()-pv[2]);
			 TVector3 mom(HPiAODtrk->Px(),HPiAODtrk->Py(),HPiAODtrk->Pz());
			 Double_t pta = mom.Angle(fline);
			 
			 TVector3 flineXY(part->GetSecVtxX()-pv[0],
			 part->GetSecVtxY()-pv[1],
			 0.);
			 TVector3 momXY(HPiAODtrk->Px(),HPiAODtrk->Py(),0.);
			 Double_t ptaXY = momXY.Angle(flineXY);
			 
			 bPctp = TMath::Cos(pta);
			 bPctpXY = TMath::Cos(ptaXY);
			 */
			fMCTruth->Fill(17);
			fMCInvMassBpm->Fill(massCandBm);
			fMCInvMassDiff->Fill(massCandBm-(d->InvMassD0()));
			fMCCosThetaStar->Fill(cosThStaBm);
			fMCCosPoinAngle->Fill(cosPtAngBm);
			fMCCosPoinAngleXY->Fill(cosPtAngBmXY);
			fMCTransvMome->Fill(transMomBm);
			fMCDCA->Fill(DCABm);
			fMCTransMomD->Fill(PtDBm);
			fMCTransMomPi->Fill(PtPiBm);
			fMCd0D->Fill(d0DBm);
			fMCd0Pi->Fill(d0PiBm);
			fMCProdd0d0->Fill(Prodd0d0Bm);
			fMCNormDecLthXY->Fill(normDecLthXYBm);
			fMCImpParXY->Fill(impparXYBm);
			fMCctau->Fill(ctauBm);
			hThetaD0->Fill(part->ThetaProng(1));
			hThetaPi->Fill(part->ThetaProng(0));
			hThetaD0Ka->Fill(d->ThetaProng(1));
			hThetaD0Pi->Fill(d->ThetaProng(0));
		}
		if((HPiAODtrk->Charge() > 0) && (TMath::Abs(massCandBp - massTrueB)<massWindowB)){
			if(fDebug>2) cout << "Invariant Mass B+: " << massCandBp << endl;
			
			if(fKillCandidates){
				// Advanced cut tuning applied here
				//KillCandidatesAdvanced(part,d,kTRUE);
			}
			
			massCandDiff = massCandBp-d->InvMassD0bar();
			//	massCandDiff = CalculateInvMass(part,d->InvMassD0bar())-d->InvMassD0bar();
			cosThStaBp = part->CosThetaStar(0,521,211,421);
			cosPtAngBp = part->CosPointingAngleXY();
			cosPtAngBpXY = part->CosPointingAngleXY();
			transMomBp = part->Pt();
			DCABp = part->GetDCA();
			PtDBp = part->PtProng(1);
			PtPiBp = part->PtProng(0);
			d0DBp = part->Getd0Prong(1);
			d0PiBp = part->Getd0Prong(0);
			Prodd0d0Bp = part->Prodd0d0();
			normDecLthXYBp = part->NormalizedDecayLengthXY();
			impparXYBp = part->ImpParXY()*1000.;
			ctauBp = part->Ct(521);
			
			TString fillthis="";
			Int_t nbins=fCuts->PtBin(part->Pt());
			
			fillthis="hMCdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="hMCnormdeclgt";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="hMCdeclgtxy";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="hMCnormdeclgtXY";
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="hMCdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
			fillthis="hMCnormdeclxyd0d0";
			((TH2F*)fOutputList->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
			
			fillthis="InvMassBplus_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandBp);
			fillthis="InvBMassDiff_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(massCandDiff);
			fillthis="CosThetaStar_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosThStaBp);
			fillthis="CosPointAngle_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBp);
			fillthis="CosPointAngleXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(cosPtAngBpXY);
			fillthis="Bpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(transMomBp);
			fillthis="Bdca_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(DCABp);
			fillthis="Dpt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtDBp);
			fillthis="Pipt_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(PtPiBp);
			fillthis="Dd0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0DBp);
			fillthis="Pid0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(d0PiBp);
			fillthis="Prodd0d0_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(Prodd0d0Bp);
			fillthis="BimpParXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(impparXYBp);
			fillthis="ctau_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(ctauBp);
			fillthis="DecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLength());
			fillthis="normDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLength());
			fillthis="DecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->DecayLengthXY());
			fillthis="normDecayLengthXY_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(part->NormalizedDecayLengthXY());
			fillthis="DLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(part->DecayLengthXY(),Prodd0d0Bp);
			fillthis="normDLXYd0d0_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(normDecLthXYBp,Prodd0d0Bp);
			fillthis="LengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->DecayLength(),part->DecayLength());
			fillthis="NormLengthBvsD_";
			fillthis+=nbins;
			((TH2F*)fOutputListMCVar->FindObject(fillthis))->Fill(d->NormalizedDecayLength(),part->NormalizedDecayLength());
			fillthis="MometumAngleDwrtB_";
			fillthis+=nbins;
			((TH1F*)(fOutputListMCVar->FindObject(fillthis)))->Fill(pCosine);
			fillthis="normDDecayLength_";
			fillthis+=nbins;
			((TH1F*)(fOutputListVar->FindObject(fillthis)))->Fill(d->NormalizedDecayLength());
			
			// Checking number of B versus Rapidity
			/*
			 if(TMath::Abs(part->Y(521))>0.0 && TMath::Abs(part->Y(521))<=0.1){
				fillthis="hSgn_CountsY_PtBins_0[1]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.1 && TMath::Abs(part->Y(521))<=0.2){
				fillthis="hSgn_CountsY_PtBins_0[2]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.2 && TMath::Abs(part->Y(521))<=0.3){
				fillthis="hSgn_CountsY_PtBins_0[3]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.3 && TMath::Abs(part->Y(521))<=0.4){
				fillthis="hSgn_CountsY_PtBins_0[4]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.4 && TMath::Abs(part->Y(521))<=0.5){
				fillthis="hSgn_CountsY_PtBins_0[5]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.5 && TMath::Abs(part->Y(521))<=0.6){
				fillthis="hSgn_CountsY_PtBins_0[6]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.6 && TMath::Abs(part->Y(521))<=0.7){
				fillthis="hSgn_CountsY_PtBins_0[7]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.7 && TMath::Abs(part->Y(521))<=0.8){
				fillthis="hSgn_CountsY_PtBins_0[8]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.8 && TMath::Abs(part->Y(521))<=0.9){
				fillthis="hSgn_CountsY_PtBins_0[9]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>0.9 && TMath::Abs(part->Y(521))<=1.0){
				fillthis="hSgn_CountsY_PtBins_0[10]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>1.0 && TMath::Abs(part->Y(521))<=1.1){
				fillthis="hSgn_CountsY_PtBins_0[11]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>1.1 && TMath::Abs(part->Y(521))<=1.2){
				fillthis="hSgn_CountsY_PtBins_0[12]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>1.2 && TMath::Abs(part->Y(521))<=1.3){
				fillthis="hSgn_CountsY_PtBins_0[13]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>1.3 && TMath::Abs(part->Y(521))<=1.4){
				fillthis="hSgn_CountsY_PtBins_0[14]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 if(TMath::Abs(part->Y(521))>1.4){
				fillthis="hSgn_CountsY_PtBins_0[15]";
				((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(part->Pt());
			 }
			 */
			bBInvMass			= part->InvMass(2,pdgBp);
			bBInvMassDiff	= massCandDiff;
			bBpT					= part->Pt();
			bBctp					= part->CosPointingAngle();
			bBctpXY				= part->CosPointingAngleXY();
			bBd0XY				= TMath::Abs(part->ImpParXY());
			bBPipt				= part->PtProng(0);
			bd0d0					= part->Prodd0d0();
			bBndlXY				= part->NormalizedDecayLengthXY();
			bVtxChi2			= part->GetReducedChi2();
			bBdl					= part->DecayLength();
			bDpT					= part->PtProng(1);
			bDndl					= d->NormalizedDecayLength();
			bDdl					= d->DecayLength();
			bDInvMass			= d->InvMassD0bar();
			bDdca					= d->GetDCA();
			bDcts					= d->CosThetaStarD0bar();;
			bDKapt				= d->PtProng(0);
			bDPipt				= d->PtProng(1);
			bDd0Ka				= d->Getd0Prong(0);
			bDd0Pi				= d->Getd0Prong(1);
			bDd0d0					= d->Prodd0d0();
			bDctp					= d->CosPointingAngle();
			bDctpXY				= d->CosPointingAngleXY();
			bDndlXY				= d->NormalizedDecayLengthXY();
			
			AliAODTrack *daughterBPion	= (AliAODTrack*)part->GetDaughter(0);	// Pi
			AliAODTrack *daughterD0		= (AliAODTrack*)d->GetDaughter(0);		// Ka
			AliAODTrack *daughterD1		= (AliAODTrack*)d->GetDaughter(1);		// Pi
			Double_t d0BPion[2],covd0BPion[3];
			Double_t d0BD0[2],covd0BD0[3];
			Double_t d0Ddgh0[2],covd0Ddgh0[3];
			Double_t d0Ddgh1[2],covd0Ddgh1[3];

			daughterBPion->PropagateToDCA(fvtx1,fBzkG,100.,d0BPion,covd0BPion);
			d->PropagateToDCA(fvtx1,fBzkG,100.,d0BD0,covd0BD0);
			daughterD0->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh0,covd0Ddgh0);
			daughterD1->PropagateToDCA(fvtx1,fBzkG,100.,d0Ddgh1,covd0Ddgh1);
			
			
		 dcaZBPi			= d0BPion[1];
		 dcaZD0			= d0BD0[1];
		 dcaZD0Ka		= d0Ddgh0[1];
		 dcaZD0Pi		= d0Ddgh1[1];
			
			
			
			/*
			 Double_t pv[3] = {0.,0.,0.};
			 fvtx1->GetXYZ(pv);
			 TVector3 fline(part->GetSecVtxX()-pv[0],
			 part->GetSecVtxY()-pv[1],
			 part->GetSecVtxZ()-pv[2]);
			 TVector3 mom(HPiAODtrk->Px(),HPiAODtrk->Py(),HPiAODtrk->Pz());
			 Double_t pta = mom.Angle(fline);
			 
			 TVector3 flineXY(part->GetSecVtxX()-pv[0],
			 part->GetSecVtxY()-pv[1],
			 0.);
			 TVector3 momXY(HPiAODtrk->Px(),HPiAODtrk->Py(),0.);
			 Double_t ptaXY = momXY.Angle(flineXY);
			 
			 bPctp = TMath::Cos(pta);
			 bPctpXY = TMath::Cos(ptaXY);
			 */
			fMCTruth->Fill(18);
			fMCInvMassBpm->Fill(massCandBp);
			fMCInvMassDiff->Fill(massCandBp-(d->InvMassD0bar()));
			fMCCosThetaStar->Fill(cosThStaBp);
			fMCCosPoinAngle->Fill(cosPtAngBp);
			fMCCosPoinAngleXY->Fill(cosPtAngBpXY);
			fMCTransvMome->Fill(transMomBp);
			fMCDCA->Fill(DCABp);
			fMCTransMomD->Fill(PtDBp);
			fMCTransMomPi->Fill(PtPiBp);
			fMCd0D->Fill(d0DBp);
			fMCd0Pi->Fill(d0PiBp);
			fMCProdd0d0->Fill(Prodd0d0Bp);
			fMCNormDecLthXY->Fill(normDecLthXYBp);
			fMCImpParXY->Fill(impparXYBp);
			fMCctau->Fill(ctauBp);
			hThetaD0->Fill(part->ThetaProng(1));
			hThetaPi->Fill(part->ThetaProng(0));
			hThetaD0Ka->Fill(d->ThetaProng(0));
			hThetaD0Pi->Fill(d->ThetaProng(1));
		}
	}
	fMCTruth->Fill(19);
	fSelectionVariables->Fill();
	return;
}
//____________________________________________________________________________
void AliAnalysisTaskSEBpmMass::Terminate(Option_t */*option*/){
	// Terminate analysis
	//
	if(fDebug > 1) printf("AnalysisTaskSEBpmMass: Terminate() \n");
	
	//	TCanvas *c1 = new TCanvas();
	//fNentries->Draw("htext0");
	//	TCanvas *c2 = new TCanvas();
	//	fMCTruth->Draw("htext0");
	
	return;
}
//____________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEBpmMass::ReconstructSecondaryVertex(TObjArray *trkArray, Double_t &dispersion,Bool_t useTRefArray) const{
	// Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle
	
	AliESDVertex *vertexESD = 0;
	AliAODVertex *vertexAOD = 0;
	
	// AliVertexerTracks
	
	AliVertexerTracks *vertexer = new AliVertexerTracks(fBzkG);
	vertexer->SetVtxStart((AliESDVertex*)fvtx1);
	vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
	delete vertexer; vertexer=NULL;
	
	if(!vertexESD) return vertexAOD;
	if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) {
		//AliDebug(2,"vertexing failed");
		delete vertexESD; vertexESD=NULL;
		return vertexAOD;
	}
	
	Double_t vertRadius2=vertexESD->GetX()*vertexESD->GetX()+vertexESD->GetY()*vertexESD->GetY();
	if(vertRadius2>8.){
		// vertex outside beam pipe, reject candidate to avoid propagation through material
		delete vertexESD; vertexESD=NULL;
		return vertexAOD;
	}
	// convert to AliAODVertex
	Double_t pos[3],cov[6],chi2perNDF;
	vertexESD->GetXYZ(pos); // position
	vertexESD->GetCovMatrix(cov); //covariance matrix
	chi2perNDF = vertexESD->GetChi2toNDF();
	dispersion = vertexESD->GetDispersion();
	delete vertexESD; vertexESD=NULL;
	
	Int_t nprongs= (useTRefArray ? 0 : trkArray->GetEntriesFast());
	vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
	
	return vertexAOD;
}
//____________________________________________________________________________
void AliAnalysisTaskSEBpmMass::AddDaughterRefs(AliAODVertex *v,const AliVEvent *event,const TObjArray *trkArray) const{
	// Add the AOD tracks as daughters of the vertex (TRef)
	
	Int_t nDg = v->GetNDaughters();
	TObject *dg = 0;
	if(nDg) dg = v->GetDaughter(0);
	
	if(dg) return; // daughters already added
	
	Int_t nTrks = trkArray->GetEntriesFast();
	
	AliExternalTrackParam *track = 0;
	AliAODTrack *aodTrack = 0;
	Int_t id;
	
	for(Int_t i=0; i<nTrks; i++) {
		track = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
		id = (Int_t)track->GetID();
		//printf("---> %d\n",id);
		if(id<0) continue; // this track is a AliAODRecoDecay
		aodTrack = (AliAODTrack*)event->GetTrack(id);
		v->AddDaughter(aodTrack);
	}
	
	return;
}
//____________________________________________________________________________
Int_t AliAnalysisTaskSEBpmMass::IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header){
	
	AliVertexingHFUtils* ggg=new  AliVertexingHFUtils();
	
	Int_t lab=part->GetLabel();
	if(lab<0) {delete ggg;return 1;} // discuss this with rossella
	TString nameGen=ggg->GetGenerator(lab,header);
	TString empty="";
	//cout << " FIRST CALL " << nameGen << endl;
	Int_t countControl =0;
	
	while(nameGen.IsWhitespace()){
		AliAODMCParticle *mcpart= (AliAODMCParticle*)fMcArray->At(lab);
		if(!mcpart){
			printf("AliVertexingHFUtils::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
			break;
		}
		Int_t mother = mcpart->GetMother();
		if(mother<0){
			printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
			break;
		}
		lab=mother;
		nameGen=ggg->GetGenerator(mother,header);
		countControl++;
		if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
			printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
			break;
		}
		//cout << "Testing " << nameGen << " with PDG code " << mcpart->GetPdgCode() << " and is primary? " << mcpart->IsPrimary() << " count control: " << countControl << " ! " << endl;
	}
	//cout << " FINAL CALL " << nameGen << endl;
	
	if(nameGen.IsWhitespace() || nameGen.Contains("ijing")){delete ggg; return 0;}
	
	
	// old
	/*
	 while(bbb.IsWhitespace()){
	 AliAODMCParticle *mcpart= (AliAODMCParticle*)fMcArray->At(lab);
	 if(!mcpart){delete ggg; return 1;}
	 Int_t mother = mcpart->GetMother();
	 lab=mother;
	 bbb=ggg->GetGenerator(mother,header);
	 //cout << "Testing " << bbb << " with PDG code " << mcpart->GetPdgCode() << " and is primary? " << mcpart->IsPrimary() << " ! " << endl;
	 }
	 //cout << " FINAL CALL " << bbb << endl;
	 
	 if(bbb.Contains("ijing")){delete ggg; return 0;}
	 */
	delete ggg;
	return 1;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::IsCandidateInjected(AliAODRecoDecayHF *part, AliAODMCHeader *header){
	
	
	Int_t nprongs=part->GetNProngs();
	for(Int_t i=0;i<nprongs;i++){
		AliAODTrack *daugh=(AliAODTrack*)part->GetDaughter(i);
		Int_t lab=daugh->GetLabel();
		if(lab<0) return 0;
		if(IsTrackInjected(daugh,header)) return kTRUE;
	}
	return kFALSE;
}

//____________________________________________________________________________
/*Bool_t AliAnalysisTaskSEBpmMass::VertexingKFD0(AliAODRecoDecayHF2Prong *D0cand,Int_t *pdgs) const{
	// apply vertexing KF for D0 Meson
	Int_t iprongs[2]={0,1};
	Double_t mass[2]={1.864,3.};
	Bool_t constraint=kTRUE;
	
	AliKFParticle *D0meson = D0cand->ApplyVertexingKF(iprongs,2,pdgs,constraint,fBzkG,mass);
	if(!D0meson) return kFALSE;
	//if(lambdac->GetChi2()/lambdac->GetNDF()>fCutsKF[1]) return kFALSE;
	
	return kTRUE;
 }*/

// Classes from Alex Grelli
//________________________________________________________________________________________-
AliAnalysisTaskSEBpmMass& AliAnalysisTaskSEBpmMass::operator=(const AliAnalysisTaskSEBpmMass &mbkg){
	
	//assignment operator
	
	if(&mbkg == this) return *this;
	fRot=mbkg.fRot;
	fAngleFirst=mbkg.fAngleFirst;
	fAngle=mbkg.fAngle;
	fAODMapSize = mbkg.fAODMapSize;
	
	return *this;
}
//_________________________________________________________________________________________
void AliAnalysisTaskSEBpmMass::SetRotationAngle(Double_t Rangle) {
	
	fAngle = Rangle;
	cout<< "you setted manually the rotation angle to = "<< fAngle<<endl;
	
}
//_________________________________________________________________________________________
void AliAnalysisTaskSEBpmMass::SetNRotations(Int_t Nrotations) {
	
	fRot = Nrotations;
	//cout<< "you setted manually the number of rotations to Nrot = "<< fRot<<endl;
}
//_________________________________________________________________________________________
TObjArray* AliAnalysisTaskSEBpmMass::GetArrayCandRotated(/*AliAODEvent* ev,*/AliAODRecoDecayHF2Prong *decay,AliAODTrack *pionTrack) {
	//
	// fill a TObjArray with the rotated candidates
	//
	
	if(fDebug>2) printf("I am in rotation function \n");
	TString fillthis="hRotationMonitor";
	
	Bool_t rotateFirst = kTRUE;
	Bool_t rotateSecond = kFALSE;//kFALSE;
	
	// magnetic field
	Double_t bz=fBzkG;//  ev->GetMagneticField();
	//primary vertex
	//	AliVVertex *primaryVertex=fvtx1;// ev->GetPrimaryVertex();
	//	if(!primaryVertex){
	TObjArray *charmArray = new TObjArray(fRot);
	charmArray->SetOwner(kTRUE);
	
	if(!fvtx1){
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(4); // No primary vertex
		//delete primaryVertex;
		return charmArray;
	}
	
	Double_t pseudoX2[3], pseudoP2[3];
	Double_t CovPseudo2[21],CovPseudo1[21];
	Double_t pseudoX1[3],pseudoP1[3],Xrot[3];
	
	// for rotations
	AliExternalTrackParam * et1;
	AliNeutralTrackParam * et2;
	
	if(fDebug>2) printf("Prepared Parameters \n");
	
	// positive track to be rotated
	AliAODTrack* positive = (AliAODTrack*)decay->GetDaughter(0);
	if(fDebug>2) {
		printf("Got daughter 0 \n");
		positive->Print();
	}
	
	positive->GetCovarianceXYZPxPyPz(CovPseudo2);
	positive->GetXYZ(pseudoX2);
	positive->GetPxPyPz(pseudoP2);
	Short_t sign = positive->Charge();
	// if you like to rotate the first daughter
	/*	if(rotateSecond)  {
		//et1 = new AliExternalTrackParam(pionTrack);
		et1 = new AliExternalTrackParam(pseudoX2,pseudoP2,CovPseudo2,sign);
		if(fDebug>2) printf("Prepared daughter 1 with charge %i \n",positive->Charge());
	 }*/
	// negative track // modified to AliNeutralTrackParam by Johannes
	AliAODTrack* negative = (AliAODTrack*)decay->GetDaughter(1);
	if(fDebug>2){
		printf("Got daughter 1 \n");
		negative->Print();
	}
	negative->GetCovarianceXYZPxPyPz(CovPseudo1);
	if(fDebug>2) printf("Got daughter 1 \n");
	negative->GetXYZ(pseudoX1);
	if(fDebug>2) printf("Got daughter 1 \n");
	negative->GetPxPyPz(pseudoP1);
	if(fDebug>2) printf("Got daughter 1 \n");
	Short_t sign1 = negative->Charge();
	/*	if(rotateFirst) {
		if(fDebug>2) printf("Getting D0 \n");
		et2 = new AliNeutralTrackParam(pseudoX1,pseudoP1,CovPseudo1,sign1);
		if(fDebug>2) printf("Prepared daughter 2 with charge %i \n",negative->Charge());
	 }*/
	Double_t Prot[3];
	
	// vector with rotated candidates
	//	TObjArray *charmArray = new TObjArray(fRot);
	//	charmArray->SetOwner(kTRUE);
	
	Double_t d0z0[2],covd0z0[3],d0[2],d0err[2];
	Double_t xdummy=0.,ydummy=0.,dca;
	
	Double_t Angle = (TMath::Pi() - fAngle*(fRot - 1.)/2.);
	
	//Double_t px[2],py[2],pz[2];
	UShort_t id[2];
	
	// parameters of the actual 2Prong
	//	for (Int_t i=0;i<2;++i) {
	//		const AliAODTrack *t=static_cast<AliAODTrack*>(decay->GetDaughter(i));
	//		px[i]=t->Px();
	//		py[i]=t->Py();
	//		pz[i]=t->Pz();
	//		d0[i]= decay->Getd0Prong(i);           // momenta of the pion change. It need to be recalculated
	//		d0err[i]= decay->Getd0errProng(i);     // same as before
	//		id[i]= decay->GetProngID(i);
	//	}	dca = decay->GetDCA(); // will be recalculated
	
	id[0]=pionTrack->GetID();
	id[1]=negative->GetID();
	
	if(fDebug>2) printf("Rotating tracks \n");
	
	// rotate the first track in the xy plane
	for (Int_t r = 0; r < fRot; r++)
	{
		if(fDebug>2) printf("Rotating %i \n",r);
		
		if(rotateFirst){
			if(fDebug>2) printf("Rotating First \n");
			
			// rotate momenta in XY
			Prot[0] = TMath::Cos(fAngleFirst+ r*Angle)*pseudoP2[0] - TMath::Sin(fAngleFirst+ r*Angle)*pseudoP2[1];
			Prot[1] = TMath::Sin(fAngleFirst+ r*Angle)*pseudoP2[0] + TMath::Cos(fAngleFirst+ r*Angle)*pseudoP2[1];
			Prot[2] = pseudoP2[2];
			// position, not needed at the moment
			Xrot[0] = TMath::Cos(fAngleFirst+ r*Angle)*pseudoX2[0] - TMath::Sin(fAngleFirst+ r*Angle)*pseudoX2[1];
			Xrot[1] = TMath::Sin(fAngleFirst+ r*Angle)*pseudoX2[0] + TMath::Cos(fAngleFirst+ r*Angle)*pseudoX2[1];
			Xrot[2] = pseudoX2[2];
			et2 = new AliNeutralTrackParam(pseudoX1,pseudoP1,CovPseudo1,sign1);
			
			et1 = new AliExternalTrackParam(pseudoX2,Prot,CovPseudo2,sign);
			if(fDebug>2) et1->Print();
		}
		else{
			if(fDebug>2) printf("Rotating Second \n");
			
			// rotate momenta in XY
			Prot[0] = TMath::Cos(fAngleFirst+ r*Angle)*pseudoP1[0] - TMath::Sin(fAngleFirst+ r*Angle)*pseudoP1[1];
			Prot[1] = TMath::Sin(fAngleFirst+ r*Angle)*pseudoP1[0] + TMath::Cos(fAngleFirst+ r*Angle)*pseudoP1[1];
			Prot[2] = pseudoP1[2];
			// position, not needed at the moment
			Xrot[0] = TMath::Cos(fAngleFirst+ r*Angle)*pseudoX2[0] - TMath::Sin(fAngleFirst+ r*Angle)*pseudoX2[1];
			Xrot[1] = TMath::Sin(fAngleFirst+ r*Angle)*pseudoX2[0] + TMath::Cos(fAngleFirst+ r*Angle)*pseudoX2[1];
			Xrot[2] = pseudoX2[2];
			
			et1 = new AliExternalTrackParam(pseudoX2,pseudoP2,CovPseudo2,sign);
			
			et2 = new AliNeutralTrackParam(pseudoX1,Prot,CovPseudo1,sign1);
			if(fDebug>2) et2->Print();
		}
		
		if(fDebug>2) printf("Recalculating \n");
		TObjArray ta12;
		//	ta12.SetOwner();
		ta12.Add(et1); ta12.Add(et2);
		
		//recalculate the secondary vertex
		// not crucial now but needed if you rotate the momenta
		
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(2); // Tries to reconstruct B vertex
		
		AliAODVertex *vtxt =RecalculateVertex((AliVVertex*)fvtx1,&ta12,bz);
		if(!vtxt) {
			((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(3); // Fails to reconstruct B vertex
			if(fDebug>1) printf("Vertex Recalculation Failed at slot %i \n",r);
			
			charmArray->AddAt(0x0,r);
			if(rotateFirst) et1->Reset();
			if(rotateSecond) et2->Reset();
			ta12.Delete();ta12=0x0;
			continue;
		}
		
		// this are the new impact prameters
		// with relative errors
		//			et1->PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
		et1->PropagateToDCA((AliVVertex*)fvtx1,bz,100.,d0z0,covd0z0);
		d0[0]=-d0z0[0];
		d0err[0] = TMath::Sqrt(covd0z0[0]);
		//			et2->PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
		et2->PropagateToDCA((AliVVertex*)fvtx1,bz,100.,d0z0,covd0z0);
		d0[1]=d0z0[0];
		d0err[1] = TMath::Sqrt(covd0z0[0]);
		
		//this is the new DCA
		dca=et1->GetDCA(et2,bz,xdummy,ydummy);
		
		if(fDebug>2) printf("Daughters Momenta \n");
		
		//daughters momenta
		
		Double_t px1[2],py1[2],pz1[2];
		
		if(rotateFirst){
			px1[1]= pseudoP1[0];
			py1[1]= pseudoP1[1];
			pz1[1]= pseudoP1[2];
			px1[0]= Prot[0];
			py1[0]= Prot[1];
			pz1[0]= Prot[2];
		}else{
			px1[0]= pseudoP2[0];
			py1[0]= pseudoP2[1];
			pz1[0]= pseudoP2[2];
			px1[1]= Prot[0];
			py1[1]= Prot[1];
			pz1[1]= Prot[2];
		}
		// make the rotated D0 candidate and add to the candidate vector
		
		if(fDebug>2) printf("make the rotated candidate and add to the candidate vector \n");
		
		// this is the new rotated candidate
		AliAODRecoDecayHF2Prong *the2Prong = new AliAODRecoDecayHF2Prong(vtxt,px1,py1,pz1,d0,d0err,dca);
		if(fDebug>2) {
			printf("**************************************************** \n");
			vtxt->Dump();
			the2Prong->Dump();
			printf("**************************************************** \n");
			printf("this is the new rotated candidate \n");
		}
		vtxt->SetParent(the2Prong);
		the2Prong->SetCharge(sign);
		if(fDebug>2) printf("Charge assigned\n");
		the2Prong->GetSecondaryVtx()->AddDaughter(positive);
		if(fDebug>2) printf("1st Daughter assigned\n");
		the2Prong->GetSecondaryVtx()->AddDaughter(negative);
		if(fDebug>2) printf("2nd Daughter assigned\n");
		the2Prong->SetProngIDs(2,id);
		if(fDebug>2) printf("Adding candidate to array \n");
		
		charmArray->AddAt(the2Prong,r);
		if(fDebug>2) printf("Clearing ta12 after slot %i \n",r);
		//the2Prong->Print();
		//		ta12.Clear();
		ta12.Delete();ta12=0x0;
		//	vtxt->Delete();vtxt=0x0;
		//	if(rotateFirst) et1->Reset();
		//	if(rotateSecond) et2->Reset();
		((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(5);
	}
	//	if(rotateSecond)et1->Delete(); et1=0x0;
	//	if(rotateFirst)et2->Delete(); et2=0x0;
	
	if(fDebug>2) printf("Done with rotation \n");
	
	return charmArray;
}
//___________________________________________________________________________-
AliAODVertex* AliAnalysisTaskSEBpmMass::RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField) {
	//
	// Helper function to recalculate a vertex.
	//
	
	AliESDVertex *vertexESD = 0;
	AliAODVertex *vertexAOD = 0;
	//Double_t covmatrix[6];
	// AliVertexerTracks
	AliVertexerTracks *vertexer = new AliVertexerTracks(bField);
	vertexer->SetVtxStart((AliESDVertex*)primary);//primary vertex
	vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(tracks);
	delete vertexer; vertexer=NULL;
	
	if(!vertexESD){ if(fDebug>2) printf("No ESD vertex"); return vertexAOD;}
	
	if(vertexESD->GetNContributors()!=tracks->GetEntriesFast()) {
		//tracks->Print();
		if(fDebug>2){
			printf("No vertex - number of contributors (%d != %d) wrong \n",vertexESD->GetNContributors(),tracks->GetEntriesFast());
		}
		delete vertexESD; vertexESD=NULL;
		return vertexAOD;
	}
	
	Double_t vertRadius2=vertexESD->GetX()*vertexESD->GetX()+vertexESD->GetY()*vertexESD->GetY();
	if(vertRadius2>8.){//(2.82)^2 radius beam pipe
		if(fDebug>2){
			printf("vertex outside beam pipe, reject candidate to avoid propagation through material \n ");
		}
		delete vertexESD; vertexESD=NULL;
		return vertexAOD;
	}
	// convert to AliAODVertex
	//
	Double_t dispersion;
	Double_t pos[3],cov[6],chi2perNDF;
	for(Int_t a=0;a<3;a++)pos[a]=0.;
	for(Int_t b=0;b<6;b++)cov[b]=0.;
	chi2perNDF=0;
	//
	vertexESD->GetXYZ(pos); // position
	vertexESD->GetCovMatrix(cov); //covariance matrix
	chi2perNDF = vertexESD->GetChi2toNDF();
	dispersion = vertexESD->GetDispersion();
	delete vertexESD; vertexESD=NULL;
	Int_t nprongs= tracks->GetEntriesFast();
	vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
	return vertexAOD;
	
}
//_______________________________________________________________
void AliAnalysisTaskSEBpmMass::FillRotateTrackHists(AliAODRecoDecayHF2Prong* BRotProng, AliAODRecoDecayHF2Prong* D0cand, Double_t sD0c){
	
	if(fDebug>1) printf("Filling Rotation Histograms \n");
	
	TString fillthis="";
	Double_t massTrueB = 5.279;
	Double_t massTrueD0 = 1.8648;
	
	UInt_t pdgB[2]={211,421};
	
	Double_t massCandB = BRotProng->InvMass(2,pdgB);
	Double_t massCandD = 0.;
	
	if(BRotProng->Charge()==0) return;	// should not be possible at this stage
	if(sD0c==0) return;					// should not be possible at this stage
	else if(sD0c==1 || (sD0c==3 && BRotProng->Charge()<0)){
		if(TMath::Abs(D0cand->InvMassD0() - massTrueD0) > massWindowD){
			return;
		}
		if((TMath::Abs(massCandB - massTrueB)>massWindowB)){
			return;
		}
		if(fDebug>2) cout << "Invariant Mass B-: " << massCandB << endl;
		massCandD=D0cand->InvMassD0();
	}
	else if(sD0c==2 ||(sD0c==3 && BRotProng->Charge()>0)){
		if(TMath::Abs(D0cand->InvMassD0bar() - massTrueD0) > massWindowD){
			return;
		}
		if((TMath::Abs(massCandB - massTrueB)>massWindowB)){
			return;
		}
		if(fDebug>2) cout << "Invariant Mass B+: " << massCandB << endl;
		massCandD=D0cand->InvMassD0bar();
	}
	else return;
	
	Double_t Bpt	= BRotProng->Pt();
	Double_t cts	= BRotProng->CosThetaStar(0,521,211,421);
	Double_t ctp	= BRotProng->CosPointingAngleXY();
	Double_t ctpXY	= BRotProng->CosPointingAngleXY();
	Double_t dca	= BRotProng->GetDCA();
	Double_t ptD	= BRotProng->PtProng(1);
	Double_t ptPi	= BRotProng->PtProng(0);
	Double_t d0D	= BRotProng->Getd0Prong(1);
	Double_t d0Pi	= BRotProng->Getd0Prong(0);
	Double_t d0d0	= BRotProng->Prodd0d0();
	Double_t decL	= BRotProng->DecayLength();
	Double_t nDecL	= BRotProng->NormalizedDecayLength();
	Double_t decLXY	= BRotProng->DecayLengthXY();
	Double_t nDecLXY= BRotProng->NormalizedDecayLengthXY();
	Double_t iParXY	= BRotProng->ImpParXY()*1000.;
	Double_t ctau	= BRotProng->Ct(521);
	
	bRotBInvMass			= massCandB;
	bRotBpT						= BRotProng->Pt();
	bRotBctp					= BRotProng->CosPointingAngle();
	bRotBctpXY				= BRotProng->CosPointingAngleXY();
	bRotBd0XY					= TMath::Abs(BRotProng->ImpParXY());
	bRotBPipt					= BRotProng->PtProng(0);
	bRotd0d0					= BRotProng->Prodd0d0();
	bRotBndlXY				= BRotProng->NormalizedDecayLengthXY();
	bRotVtxChi2				= BRotProng->GetReducedChi2();
	bRotBdl						= BRotProng->DecayLength();
	bRotDpT						= BRotProng->PtProng(1);
	bRotDndl					= D0cand->NormalizedDecayLength();
	bRotDdl						= D0cand->DecayLength();
	bRotDdca					= D0cand->GetDCA();
	bRotDd0d0					= D0cand->Prodd0d0();
	bRotDctp					= D0cand->CosPointingAngle();
	bRotDctpXY				= D0cand->CosPointingAngleXY();
	bRotDndlXY				= D0cand->NormalizedDecayLengthXY();
	
	
	Int_t nbins=fCuts->PtBin(BRotProng->Pt());
	fillthis="hRotationMonitor";
	((TH1F*)(fOutputList->FindObject(fillthis)))->Fill(1);
	
	if(BRotProng->Charge()<0){
		bRotBInvMassDiff  = bRotBInvMass-D0cand->InvMassD0();;
		bRotDInvMass			= D0cand->InvMassD0();
		bRotDcts					= D0cand->CosThetaStarD0();;
		bRotDKapt					= D0cand->PtProng(1);
		bRotDPipt					= D0cand->PtProng(0);
		bRotDd0Ka					= D0cand->Getd0Prong(1);
		bRotDd0Pi					= D0cand->Getd0Prong(0);
	}
	else if(BRotProng->Charge()>0){
		bRotBInvMassDiff  = bRotBInvMass-D0cand->InvMassD0bar();;
		bRotDInvMass			= D0cand->InvMassD0bar();
		bRotDcts					= D0cand->CosThetaStarD0bar();;
		bRotDKapt					= D0cand->PtProng(0);
		bRotDPipt					= D0cand->PtProng(1);
		bRotDd0Ka					= D0cand->Getd0Prong(0);
		bRotDd0Pi					= D0cand->Getd0Prong(1);
	}
	fRotatedCandidates->Fill();
	return;
}
//_______________________________________________________________
Double_t AliAnalysisTaskSEBpmMass::CalculateInvMass(AliAODRecoDecayHF2Prong *BMeson, Double_t mD0){
	
	Double_t ePi = BMeson->EProng(0,211);
	Double_t pD0 = BMeson->PProng(1);
	Double_t eD0 = TMath::Sqrt(mD0*mD0+pD0*pD0);
	Double_t eSum = ePi+eD0;
	
	Double_t pB = BMeson->P();
	Double_t mB = TMath::Sqrt(eSum*eSum-pB*pB);
	//	cout << "This is the D mass " << mD0 << endl;
	//	cout << "This is the Pi energy " << ePi << endl;
	//	cout << "This is the B mass " << mB << endl;
	return mB;
}
//_______________________________________________________________
Bool_t AliAnalysisTaskSEBpmMass::IsFastMcFromWeakDecay(AliAODMCParticle *mcParticle){

	if(fMcArray) {
		if (!mcParticle) return kFALSE;
		Float_t codepart = (Float_t)TMath::Abs(mcParticle->GetPdgCode());
		Int_t motherLabel = mcParticle->GetMother();
		if(motherLabel==-1) return kFALSE;
		
		AliAODMCParticle *mcParticleMum= (AliAODMCParticle*)fMcArray->At(motherLabel);
		if(mcParticleMum){
			Int_t nDgh = mcParticleMum->GetNDaughters();
			Float_t codemoth = (Float_t)TMath::Abs(mcParticleMum->GetPdgCode());
			if(codemoth>10000000) return kFALSE;
			Int_t mfl = Int_t (codemoth / TMath::Power(10, Int_t(TMath::Log10(codemoth))));
			Double_t massMum = TDatabasePDG::Instance()->GetParticle(codemoth)->Mass();
			if(mfl>3) {return kFALSE;}
			//	if(mfl!=3) return kFALSE;
			if(codemoth>=331 && codemoth<=337){return kFALSE;}
			for(Int_t i=0;i<nDgh;i++){
				AliAODMCParticle *mcParticleDgh= (AliAODMCParticle*)fMcArray->At(mcParticleMum->GetLastDaughter()-i);
				if(mcParticleDgh){
					Int_t dfl = 0;
					Float_t codedgh = (Float_t)TMath::Abs(mcParticleDgh->GetPdgCode());
					if(codedgh>10000000 && codepart==13){return kTRUE;} // muon decays --> Weak decays ala AliStack
					if(codedgh>10000000 && codepart!=13){return kFALSE;}
					dfl = Int_t (codedgh / TMath::Power(10, Int_t(TMath::Log10(codedgh))));
					Double_t massDgh = TDatabasePDG::Instance()->GetParticle(codedgh)->Mass();
					if(codepart==11 && codemoth==13){return kTRUE;}
					if(codepart==211 && mfl==dfl && ((codedgh<=3334 && codedgh>=3300) || (codedgh<=4332 && codedgh>=4334) || (codemoth<=3334 && codemoth>=3300) || (codemoth<=4334 && codemoth>=4332))){
						// double/triple strange content
						if(((codedgh<3334 && codedgh>=3300) || (codedgh<=4332 && codedgh>=4334)) && ((codemoth<3334 && codemoth>=3300) || (codemoth<=4334 && codemoth>=4332))){return kFALSE;} // double strange particles decay into double strange particles
						if((codedgh<3334 && codedgh>=3300) && (codemoth<=3334)){return kTRUE; }// triple strange into double strange + weak decay
						return kTRUE; // double strange particles decay into single strange + weak decay
					}
					if ((dfl>=mfl))/* && codepart!=11) || (codepart==11 && mfl!=1 && dfl>=mfl && nDgh==1))*/{return	kFALSE;}
					if(massDgh>massMum){return kFALSE;}
				}
				else return kFALSE;
			}
			return kTRUE;
		}
	}
	return kFALSE;
}
