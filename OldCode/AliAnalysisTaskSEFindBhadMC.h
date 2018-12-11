#ifndef ALIANALYSISTASKSEFINDBHADMC_H
#define ALIANALYSISTASKSEFINDBHADMC_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskSEFindBhadMC.h 53532 2011-12-12 11:11:37Z prino $ */ 

//*************************************************************************
// Class AliAnalysisTaskSEFindBhadMC
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison to MC truth (kinematics stored in the AOD) and cut variables
// distributions
// Authors: J.Stiller
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>

#include "AliRDHFCutsD0toKpi.h"
#include "AliAnalysisTaskSE.h"
#include "AliNormalizationCounter.h"

#include "AliGenHijingEventHeader.h"
class AliAODEvent;
class AliRDHFCuts;

class AliAnalysisTaskSEFindBhadMC : public AliAnalysisTaskSE
{
public:
	
	AliAnalysisTaskSEFindBhadMC();
	AliAnalysisTaskSEFindBhadMC(const char *name,AliRDHFCutsD0toKpi* cuts);

	virtual ~AliAnalysisTaskSEFindBhadMC();
	
	
	// Implementation of interface methods
	void SetMCLimAcc(Bool_t limAcc=kFALSE) {fMClimAcc=limAcc;}
	void SetFiducialCut(Bool_t fidCut=kTRUE) {fFiducialCut=fidCut;}
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void LocalInit() {Init();}
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *option);

	
private:
	
	AliAnalysisTaskSEFindBhadMC(const AliAnalysisTaskSEFindBhadMC &source);
	AliAnalysisTaskSEFindBhadMC& operator=(const AliAnalysisTaskSEFindBhadMC& source); 

	AliRDHFCutsD0toKpi *fCuts;      //  Cuts - sent to output slot 4
	AliGenHijingEventHeader *fhijingGenHeader;
	TList	 *fOutputList;		  //! output slot 8
	TList	 *fComparList;
	
	Bool_t fMClimAcc;
	Bool_t fFiducialCut;
	// MC histograms
	TH1I     *fNentries;	
	TH1F     *fMCTruth_Bplus;	
	TH1F     *fMCTruth_Bminus;			  
	TH1F     *fMCTruth2;			  
	TH1F     *hEtaD0;			  
	TH1F     *hEtaPi;			  
	TH1F     *hEtaBoth;
	TH2F     *hEtaBothiR;
	TH1F	 *hPtD0;
	TH1F	 *hPtPi;
	TH2F	 *hPtVsEtaD0;
	TH2F	 *hPtVsEtaPi;
	TH1F	 *hSgn_CosThetaS_0;
	TH1F	 *hSgn_CosThetaS_1;
	TH1F	 *hSgn_CosThetaS_2;
	TH1F	 *hSgn_CosThetaS_3;
	TH1F	 *hSgn_CosThetaS_4;
	TH1F	 *hSgn_CosThetaS_5;
	TH1F	 *hSgn_CosThetaS_6;
	TH1F	 *hSgn_CosThetaS_7;
	TH1F	 *hSgn_CosThetaS_8;
	TH1F	 *hSgn_CosThetaS_9;
	TH1F	 *hSgn_CosThetaS_10;
	TH1F	 *hSgn_CosThetaS_11;
	TH1F	 *hSgn_CosThetaS_12;
	TH1F	 *hSgn_CosThetaS_D_0;
	TH1F	 *hSgn_CosThetaS_D_1;
	TH1F	 *hSgn_CosThetaS_D_2;
	TH1F	 *hSgn_CosThetaS_D_3;
	TH1F	 *hSgn_CosThetaS_D_4;
	TH1F	 *hSgn_CosThetaS_D_5;
	TH1F	 *hSgn_CosThetaS_D_6;
	TH1F	 *hSgn_CosThetaS_D_7;
	TH1F	 *hSgn_CosThetaS_D_8;
	TH1F	 *hSgn_CosThetaS_D_9;
	TH1F	 *hSgn_CosThetaS_D_10;
	TH1F	 *hSgn_CosThetaS_D_11;
	TH1F	 *hSgn_CosThetaS_D_12;
	TH1F     *hThetaD0;
	TH1F     *hThetaPi;
	TH1F     *hThetaD0Ka;
	TH1F     *hThetaD0Pi;

	
	
	Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y) const;
	Int_t IsTrackInjected(Int_t lab,AliAODMCHeader *header,TClonesArray *arrayMC);
	Double_t CalculateDecayLength(AliAODMCParticle *MCParticle, AliAODVertex *aodVtx) const;
	//AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const; //Reconstruct vertex of B
	
	ClassDef(AliAnalysisTaskSEFindBhadMC,16);
};

#endif

