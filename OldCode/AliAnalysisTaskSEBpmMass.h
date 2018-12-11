#ifndef ALIANALYSISTASKSEBPMMASS_H
#define ALIANALYSISTASKSEBPMMASS_H

//*************************************************************************
// Class AliAnalysisTaskSEBpmMass
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison to MC truth (kinematics stored in the AOD) and cut variables
// distributions
// Authors: A.Dainese, andrea.dainese@ln.infn.it
// and C.Bianchin, chiara.bianchin@pd.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>
#include "TClonesArray.h"
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliNormalizationCounter.h"
#include "AliGenHijingEventHeader.h"

class AliAODEvent;

class AliAnalysisTaskSEBpmMass : public AliAnalysisTaskSE
{
public:
	
	AliAnalysisTaskSEBpmMass();
	AliAnalysisTaskSEBpmMass(const char *name,AliRDHFCutsD0toKpi* cuts);
	virtual ~AliAnalysisTaskSEBpmMass();
	
	// Implementation of interface methods
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void LocalInit() {Init();}
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *option);
	void SetReadMC(Bool_t readMC=kTRUE)		  {fReadMC=readMC;}
	void SetKillCandidates(Bool_t kill=kTRUE) {fKillCandidates=kill;}
	void SelectHIJINGOnly(Bool_t b)		   	  {fbSelectHIJING = b;}
	void IsHIJINGAvailable(Bool_t isit)		  {fbHIJINGavailable = isit;}
	void PerformRotation(Bool_t rotate)		  {fbtrackRotation = rotate;}
	void SetUseCentrality(Bool_t bCent)		  {fbCentrality = bCent;}
	
	void SetRotationAngle(Double_t Rangle);
	void SetNRotations(Int_t Nrotations); 
	TObjArray* GetArrayCandRotated(/*AliAODEvent* ev,*/AliAODRecoDecayHF2Prong *decay, AliAODTrack* HPiAODtrk);
	void FillRotateTrackHists(AliAODRecoDecayHF2Prong* BRotProng, AliAODRecoDecayHF2Prong* D0cand, Double_t sD0c);
	
private:	
	AliAnalysisTaskSEBpmMass(const AliAnalysisTaskSEBpmMass &source);
	AliAnalysisTaskSEBpmMass& operator=(const AliAnalysisTaskSEBpmMass& source); 
	void	FillMCTruthD(AliAODRecoDecayHF2Prong* d, Double_t sD0c);//,AliAODRecoDecayHF2Prong* bMeson);
	Bool_t  KillDCandidates(AliAODRecoDecayHF2Prong *Bmeson,AliAODRecoDecayHF2Prong* d, Double_t sD0c);
	Bool_t	KillCandidates(AliAODRecoDecayHF2Prong *part, Double_t sD0c,AliAODRecoDecayHF2Prong* d);//, AliAODTrack* HPiAODtrk, AliAODRecoDecayHF2Prong *d);
	Bool_t	KillCandidatesAdvanced(AliAODRecoDecayHF2Prong* Bmeson, AliAODRecoDecayHF2Prong* Dmeson, Bool_t MCtruth);
	void    FillMCtruthDHists(AliAODRecoDecayHF2Prong *part, Double_t sD0c, AliAODTrack* HPiAODtrk, AliAODRecoDecayHF2Prong *Dmeson);
	void    FillDHists(AliAODRecoDecayHF2Prong *Bmeson, Double_t sD0c, AliAODTrack *AODtrack, AliAODRecoDecayHF2Prong* d);
	Bool_t  FillBpmHists(AliAODRecoDecayHF2Prong *Bmeson, Double_t sD0c, AliAODTrack *AODtrack, AliAODRecoDecayHF2Prong* d);
	void    FillBpmMCHists(AliAODRecoDecayHF2Prong *part, Double_t sD0c, AliAODTrack *AODtrack, AliAODRecoDecayHF2Prong* d);
	Bool_t	FillMCTruthHistos(AliAODRecoDecayHF2Prong* d, Double_t sD0c, AliAODTrack* HPiAODtrk,AliAODRecoDecayHF2Prong* bMeson);
    void	AddDaughterRefs(AliAODVertex *v,const AliVEvent *event, const TObjArray *trkArray) const;
	Bool_t  IsCandidateInjected(AliAODRecoDecayHF *part,AliAODMCHeader *header);
	Int_t	IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header);
//	Bool_t	VertexingKFD0(AliAODRecoDecayHF2Prong *D0cand,Int_t *pdgs) const;
	Double_t CalculateInvMass(AliAODRecoDecayHF2Prong *BMeson, Double_t mD0);
	Bool_t IsFastMcFromWeakDecay(AliAODMCParticle *mcParticle);
	
	AliPIDResponse					*fPIDResponse;
	Bool_t									fbHIJINGavailable;		// Is HIJING available?
	Bool_t									fbSelectHIJING;			// Select only particles from HIJING event
	Bool_t									fbtrackRotation;
	Bool_t									fbCentrality;			// Flag==kTRUE if only V0M centrality 0-10% used
	AliGenHijingEventHeader *fhijingGenHeader;
	TH1F										*fNentries;            //! histogram with number of events on output slot 3
	Double_t								massWindowB;
	Double_t								massWindowD;
	
	AliRDHFCutsD0toKpi		*fCuts;      //  Cuts - sent to output slot 4
	Int_t					fArray;               //  can be D0 or Like Sign candidates
	Bool_t					fReadMC;              //  flag for MC array: kTRUE = read it, kFALSE = do not read it
	Bool_t					fKillCandidates;              //  flag for MC array: kTRUE = read it, kFALSE = do not read it
	Int_t					fSys;                 // fSys=0 -> p-p; fSys=1 ->PbPb (in this case fFillVarHists=kFALSE by default: set it to kTRUE *after* if needed)
	
	Double_t				fBzkG;                     // z component of magnetic field
	AliAODVertex			*fvtx1;                // primary vertex
	TList					*fOutputList;		      //! list send on output slot 8
	TList					*fOutputListPt;          //! list send on output slot 9
	TList    *fOutputListVar;         //! list send on output slot 10
	TList    *fOutputListMCVar;         //! list send on output slot 12
	TList	 *fOutputListDSandB;		//! output off D mesons (Sgn+Bgk)
	TList	 *fOutputListDSonly;		//! output off D mesons (Sgn only)
	// Data histograms
	
	TH1F     *fInvMassBpm;        //! histogram with invariant mass of Bpm on output slot 8
	TH1F     *fInvMassDiff;       //! histogram with invariant mass difference of Bpm-D on output slot 8
	TH1F	 *fCosThetaStar;	  //! histogram with cosine theta star of Bpm on output slot 8
	TH1F	 *fCosPoinAngle;	  //! histogram with cosine theta pointing of Bpm on output slot 8
	TH1F	 *fCosPoinAngleXY;	  //! histogram with cosine theta pointing XY of Bpm on output slot 8
	TH1F	 *fTransvMome;        //! histogram with pt of Bpm on output slot 8
	TH1F	 *fDCA;				  //! histogram with dca of Bpm on output slot 8
	TH1F	 *fTransMomD;   	  //! histogram with pt of D from Bpm on output slot 8
	TH1F	 *fTransMomPi;		  //! histogram with pt of Pi from Bpm on output slot 8
	TH1F	 *fd0D;				  //! histogram with d0 of D from Bpm on output slot 8
	TH1F	 *fd0Pi;			  //! histogram with d0 of Pi from Bpm on output slot 8
	TH1F	 *fProdd0d0;		  //! histogram with d0xd0 of Pi and D on output slot 8
	TH1F	 *fNormDecLthXY;	  //! histogram with normalized decay length of Bpm on output slot 8
	TH1F	 *fImpParXY;		  //! histogram with impact paramter XY of Bpm on output slot 8
	TH1F	 *fctau;     		  //! histogram with ctau of Bpm on output slot 8
	
	// MC histograms
	
	TH1F     *fMCTruth;			  
	TH1F     *fMCInvMassBpm;      //! MC histogram with invariant mass of Bpm on output slot 8
	TH1F     *fMCInvMassDiff;        //! histogram with invariant mass difference of Bpm-D on output slot 8
	TH1F	 *fMCCosThetaStar;	  //! MC histogram with cosine theta star of Bpm on output slot 8
	TH1F	 *fMCCosPoinAngle;	  //! MC histogram with cosine theta pointing of Bpm on output slot 8
	TH1F	 *fMCCosPoinAngleXY;  //! MC histogram with cosine theta pointing XY of Bpm on output slot 8
	TH1F	 *fMCTransvMome;      //! MC histogram with pt of Bpm on output slot 8
	TH1F	 *fMCDCA;			  //! MC histogram with dca of Bpm on output slot 8
	TH1F	 *fMCTransMomD;   	  //! MC histogram with pt of D from Bpm on output slot 8
	TH1F	 *fMCTransMomPi;	  //! MC histogram with pt of Pi from Bpm on output slot 8
	TH1F	 *fMCd0D;			  //! MC histogram with d0 of D from Bpm on output slot 8
	TH1F	 *fMCd0Pi;			  //! MC histogram with d0 of Pi from Bpm on output slot 8
	TH1F	 *fMCProdd0d0;		  //! MC histogram with d0xd0 of Pi and D on output slot 8
	TH1F	 *fMCNormDecLthXY;	  //! MC histogram with normalized decay length of Bpm on output slot 8
	TH1F	 *fMCImpParXY;		  //! MC histogram with impact paramter XY of Bpm on output slot 8
	TH1F	 *fMCctau;     		  //! MC histogram with ctau of Bpm on output slot 8
	TH1F	 *fRunNumber;		  //! run number of MC identified Bpm
	TH1F   *hThetaD0;
	TH1F   *hThetaPi;
	TH1F   *hThetaD0Ka;
	TH1F   *hThetaD0Pi;
	TTree	 *fSelectionVariables;
	TTree  *fRotatedCandidates;
	AliAODMCHeader *fmcHeader;
	TClonesArray *fMcArray;

	// by A. Grelli
	// Helper functions
	AliAODVertex* RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField);
	
	TList   *fDebugOutput; //! collection of debug output
	Int_t    fRot;            //!
	Double_t fAngleFirst;  //!
	Double_t fAngle;       //!
	Int_t    fNDebug;       // Max number of debug entries into Ntuple
	Int_t    fAODMapSize; // size of fAODMap 
	Int_t    *fAODMap; //[fAODMapSize] map between index and ID for AOD tracks
	//-----	
	
	
	Double_t dcaZBPi;
	Double_t dcaZD0;
	Double_t dcaZD0Ka;
	Double_t dcaZD0Pi;

	
	Double_t bBInvMass;
	Double_t bBInvMassDiff;
	Double_t bBpT;
	Double_t bBctp;
	Double_t bBctpXY;
	Double_t bBd0XY;
	Double_t bBPipt;
	Double_t bd0d0;
	Double_t bBndlXY;
	Double_t bVtxChi2;
	Double_t bBdl;
	Double_t bDpT;
	Double_t bDndl;
	Double_t bDdl;

	Double_t bDInvMass;
	Double_t bDdca;
	Double_t bDcts;
	Double_t bDKapt;
	Double_t bDPipt;
	Double_t bDd0Ka;
	Double_t bDd0Pi;
	Double_t bDd0d0;
	Double_t bDctp;
	Double_t bDctpXY;
	Double_t bDndlXY;
	
	Double_t bBPiTOFnSigma;
	Double_t bBPiTPCnSigma;

	Double_t bBNormIpD0;
	Double_t bBNormIpPi;
	Double_t bpdgPiB;
	Double_t bPiPhysPrimary;
	Double_t bPiWeak;
	Double_t bPiPdgCode;
	Double_t bDDgh0PhysPrimary;
	Double_t bDDgh0Weak;
	Double_t bDDgh0PdgCode;
	Double_t bDDgh1PhysPrimary;
	Double_t bDDgh1Weak;
	Double_t bDDgh1PdgCode;

	Double_t bRotBInvMass;
	Double_t bRotBInvMassDiff;
	Double_t bRotBpT;
	Double_t bRotBctp;
	Double_t bRotBctpXY;
	Double_t bRotBd0XY;
	Double_t bRotBPipt;
	Double_t bRotd0d0;
	Double_t bRotBndlXY;
	Double_t bRotVtxChi2;
	Double_t bRotBdl;
	Double_t bRotDpT;
	Double_t bRotDndl;
	Double_t bRotDdl;
	
	Double_t bRotDInvMass;
	Double_t bRotDdca;
	Double_t bRotDcts;
	Double_t bRotDKapt;
	Double_t bRotDPipt;
	Double_t bRotDd0Ka;
	Double_t bRotDd0Pi;
	Double_t bRotDd0d0;
	Double_t bRotDctp;
	Double_t bRotDctpXY;
	Double_t bRotDndlXY;

	Double_t bpdgPiBmum;
	Double_t bpdgGPiBmum;
	Double_t bpdgGGPiBmum;
	
	Double_t bpdgDdgh0;
	Double_t bpdgDdgh1;
	Double_t bpdgDcandSame;
	Double_t bpdgDmumSame;
	Double_t bpdgGDmumSame;
	Double_t bpdgDcand0;
	Double_t bpdgDmum0;
	Double_t bpdgGDmum0;
	Double_t bpdgDcand1;
	Double_t bpdgDmum1;
	Double_t bpdgGDmum1;
	Double_t bBinChannel;
	
	Double_t fBCtpCut[13];
	Double_t fD0PtCut[13];
	Double_t fDNdlCut[13];
	Double_t fBIxyCut[13];
	Double_t fPiPtCut[13];
	
	AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const; //Reconstruct vertex of B
	
	ClassDef(AliAnalysisTaskSEBpmMass,16); // AliAnalysisTaskSE for D0->Kpi
};

#endif

