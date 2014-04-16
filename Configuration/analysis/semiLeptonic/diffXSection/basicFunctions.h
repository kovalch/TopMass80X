#ifndef basicFunctions_h
#define basicFunctions_h

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <map>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <typeinfo>

#include <TROOT.h>
#include <TSystem.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TKey.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TTree.h>
#include <TLine.h>
#include <TBox.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h> 
#include <TBox.h>
#include <TGaxis.h>
#include <TError.h>
#include <TMath.h>
#include <TVectorD.h>

#include "HHStyle.h"

#ifdef DILEPTON_MACRO
namespace semileptonic {
#endif
  
  // group space folder with all rootfiles and other analysis relevant stuff
  TString groupSpace="/nfs/dust/cms/group/tophh/";
  //TString groupSpace="/afs/naf.desy.de/group/cms/scratch/tophh/";
  TString AnalysisFolder="RecentAnalysisRun8TeV_doubleKinFit"; // used mostly as default argument
  // large madspin or smaller madgraph samples?
  bool MadSpin=true; 
  // adapt several things like removing official labels for my PHD
  bool PHD=false;//false;
  // exclude ttV contribution from all plots and in calculation
  bool excludettV=true;
  // Use global switch to rm/add Preliminary to plots 
  bool globalPrelim   =false; // value used in DrawCMSLabels
  bool useGlobalPrelim=true;  // for false the value above is ignored
  // label extension for center of mass in topPt in ttbar center of mass system
  TString cmsExt="com";

  // ==========================================
  //  Cross-section variables and labels
  // ==========================================

  // basic variables
  TString xSecVariablesKinFit[]       = {"topPtLead"    , "topPtSubLead"    , "topPtTtbarSys"    , "topPt"    , "topY"    , "ttbarPt"    , "ttbarY"    , "ttbarMass"    , "ttbarDelPhi"    , "ttbarPhiStar"    };
  TString xSecVariablesKinFitNorm[]   = {"topPtLeadNorm", "topPtSubLeadNorm", "topPtTtbarSysNorm", "topPtNorm", "topYNorm", "ttbarPtNorm", "ttbarYNorm", "ttbarMassNorm", "ttbarDelPhiNorm", "ttbarPhiStarNorm"};
  TString xSecVariablesFinalState[]     = {"lepPt"    , "lepEta"    , "bqPt"    , "bqEta"    , "bbbarMass"    , "bbbarPt"    , "lbMass"    , "Njets"    , "rhos"    };

  TString xSecVariablesFinalStateNorm[] = {"lepPtNorm", "lepEtaNorm", "bqPtNorm", "bqEtaNorm", "bbbarMassNorm", "bbbarPtNorm", "lbMassNorm", "NjetsNorm", "rhosNorm"};

  TString xSecVariablesIncl[] = {"inclusive"};

  TString xSecLabelKinFit[]     = {"p_{T}^{lead t}/[GeV]", "p_{T}^{sublead t}/[GeV]", "p_{T}^{t} #scale[0.8]{(t#bar{t} "+cmsExt+")}/[GeV]", "p_{T}^{t}/[GeV]", "y^{t}/ ", "p_{T}^{t#bar{t}}/[GeV]", "y^{t#bar{t}}/ ", "m^{t#bar{t}}/[GeV]", "#Delta#phi(t,#bar{t})/ ", "#Phi^{#lower[-0.9]{* }}(t,#bar{t})/ "};
  TString xSecLabelFinalState[] = {"p_{T}^{l}/[GeV]", "#eta^{l}/ ", "p_{T}^{b}/[GeV]", "#eta^{b}/ ", "m^{b#bar{b}}/[GeV]", "p_{T}^{b#bar{b}}/[GeV]", "m^{lb}/[GeV]", "N_{jets}/ ", "#rho_{S}/ "};

  // cross-check variables
 
  TString xSecVariablesCCVar[]     = {"topPtPlus",    "topPtMinus",    "topYPlus",    "topYMinus",    "lepEtaPlus",    "lepEtaMinus"    };
  TString xSecVariablesCCVarNorm[] = {"topPtPlusNorm","topPtMinusNorm","topYPlusNorm","topYMinusNorm","lepEtaPlusNorm","lepEtaMinusNorm"};

  TString xSecLabelCCVar[] = {"p_{T}^{t}/[GeV]", "p_{T}^{#bar{t}}/[GeV]", "y^{t}/ ", "y^{#bar{t}}/ ", "#eta^{l^{+}}/ ", "#eta^{l^{-}}/ "};

  // =================================================
  //  Parameters for samples and style options 
  // =================================================
  
  // Enumerator to assign an unique integer value to each sample

                   /*0:*/    /*1:*/    /*2:*/    /*3:*/    
  enum samples    {kSig    , kBkg    , kZjets  , kWjets  ,
                   /*4b:*/
                   kTTVjets, //######
		   /*4:*/    /*5:*/     
		   kQCD    , kSTop   ,   
                    /*6:*/    /*7:*/                 
                   kDiBos  , kData,   
                   /*8b*/     /*9b*/     /*10b*/
                   kTTGjets , kTTZjets, kTTWjets,  //###########
		   /*8*/     /*9*/     /*10*/    /*11*/    /*12*/    /*13*/
		   kQCDEM1 , kQCDEM2 , kQCDEM3 , kQCDBCE1, kQCDBCE2, kQCDBCE3,  
		   /*14*/    /*15*/    /*16*/
		   kWW     , kWZ     , kZZ     , 
		   /*17*/    /*18*/    /*19*/    /*20*/    /*21*/    /*22*/
		   kSTops  , kSATops , kSTopt  , kSATopt , kSToptW , kSAToptW,
		   /*23*/    /*24*/    /*25*/    /*26*/    /*27*/    /*28*/
		   kSigPow , kBkgPow,  kSigMca , kBkgMca , kSigPowHer , kBkgPowHer,
                   /*29*/
		   ENDOFSAMPLEENUM,
		   // stop scale sample is splitted in 3 subsamples:
		   /*30*/    /*31*/     /*32*/    /*33*/     /*34*/    /*35*/
		   kSToptW1, kSAToptW1, kSToptW2, kSAToptW2, kSToptW3, kSAToptW3,
		   /*35*/
		   ENDOFSAMPLEENUM2
  };

  // Colors for event samples (<=kSAToptW)

  int color_[] ={ kRed+1  , kRed-7  , kAzure-2, kGreen-3,
                  kOrange-2,//########## 
		  kYellow , kMagenta, 
                  10      , kBlack  , 
                  kOrange-2 , kOrange-2, kOrange-2, //####### 
		  kYellow , kYellow , kYellow , kYellow , kYellow , kYellow ,
		  10      , 10      , 10      , 
		  kMagenta, kMagenta, kMagenta, kMagenta, kMagenta, kMagenta };

  // Colors for data and theory curves

  const unsigned int constDataColor     = kBlack;
  const unsigned int constMadgraphColor = kRed+1;
  const unsigned int constMcatnloColor  = kBlue;
  const unsigned int constPowhegColor   = kGreen+1;
  const unsigned int constPowhegColor2  = kGreen+3;
  const unsigned int constNnloColor     = kViolet-6;//kOrange+4
  const unsigned int constNnloColor2    = kMagenta+2;
  const unsigned int constMadgraphPerugiaColor = kRed-7;
  const unsigned int constMadgraphPerugiaMpiHiColor = kCyan;
  const unsigned int constMadgraphPerugiaNoCRColor  = kRed-5;

  // Line style for theory curves

  const unsigned int constPowhegStyle  = 7;
  const unsigned int constPowhegStyle2 = 4;
  const unsigned int constNnloStyle    = 2;
  const unsigned int constMcatnloStyle = 5;
  const unsigned int constNnloStyle2   = 6;//10
  const unsigned int constMadgraphPerugiaStyle = 3;
  const unsigned int constMadgraphPerugiaNoCRStyle  = 3;
  const unsigned int constMadgraphPerugiaMpiHiStyle = 3;

  // legend entries for theory curves
  const TString constMadGraphPythiaLabel   = "MadGraph+Pythia";
  const TString constMadGraphPythiaLabelSC = "MadGraph+Pythia SC";
  const TString constPowhegPythiaLabel  = "Powheg+Pythia";
  const TString constPowhegHerwigLabel  = "Powheg+Herwig";
  const TString constNnloLabelKidonakis = "Approx. NNLO";
  const TString constMcatnloHerwigLabel = "MC@NLO+Herwig";
  const TString constNloNNLLLabelAhrens = "NLO+NNLL";
  const TString constMadGraphPythiaPerugiaLabel      = "MadGraph+Pythia(P11)";
  const TString constMadGraphPythiaPerugiaNoCRLabel  = "MadGraph+Pythia(P11,noCR)";
  const TString constMadGraphPythiaPerugiaMpiHiLabel = "MadGraph+Pythia(P11,MpiHi)";

  // theory arxiv numbers
  const TString constApproxNNLOLabel="arXiv:1210.7813";
  const TString constNLONNLLLabel   ="arXiv:1003.5827";

  // Marker style (<=kSAToptW)

  int marker_[] = {20, 22, 29, 23, 
                   23,//############ 
		   21, 27,
                   28, 20, 
                   23, 23, 23,//######  
		   21, 21, 21, 21, 21, 21,
		   28, 28, 28,
		   27, 27, 27, 27, 27, 27};

  // Enumerator to assign an unique integer value to each systematic variation

  enum systematicVariation { /* 0:*/ sysNo,
			     /* 1:*/ sysLumiUp,                  /* 2:*/ sysLumiDown,                
			     /* 3:*/ sysPUUp,                    /* 4:*/ sysPUDown,                  
			     /* 5:*/ sysJESUp,                   /* 6:*/ sysJESDown,                 
			     /* 7:*/ sysJERUp,                   /* 8:*/ sysJERDown,                 
			     /* 9:*/ sysLepEffSFNormUp,          /*10:*/ sysLepEffSFNormDown,        
			     /*11:*/ sysLepEffSFShapeUpEta,      /*12:*/ sysLepEffSFShapeDownEta,
			     /*13:*/ sysLepEffSFShapeUpPt,       /*14:*/ sysLepEffSFShapeDownPt,  
			     //      sysLepEffSFShapeUpPt40,             sysLepEffSFShapeDownPt40,
		             //      sysLepEffSFNormUpStat,              sysLepEffSFNormUpStat,
			     // old: sysTriggerEffSFJetNormUp,           sysTriggerEffSFJetNormDown,
			     //      sysTriggerEffSFJetShapeUp,          sysTriggerEffSFJetShapeDown,    
			     /*15:*/ sysBtagSFUp,                /*16:*/ sysBtagSFDown, 
			     // old: sysBtagSFShapeUpPt100,              sysBtagSFShapeDownPt100,
			     //      sysBtagSFShapeUpEta1p2,             sysBtagSFShapeDownEta1p2,
			     /*17:*/ sysBtagSFShapeUpPt65,       /*18:*/ sysBtagSFShapeDownPt65,
			     /*19:*/ sysBtagSFShapeUpEta0p7,     /*20:*/ sysBtagSFShapeDownEta0p7,
			     /*21:*/ sysMisTagSFUp,              /*22:*/ sysMisTagSFDown,
			     /*23:*/ sysTopScaleUp,              /*24:*/ sysTopScaleDown,            
			     /*25:*/ sysVBosonScaleUp,           /*26:*/ sysVBosonScaleDown,          
			     /*27:*/ sysSingleTopScaleUp,        /*28:*/ sysSingleTopScaleDown,     
			     /*29:*/ sysTopMatchUp,              /*30:*/ sysTopMatchDown,            
			     /*31:*/ sysVBosonMatchUp,           /*32:*/ sysVBosonMatchDown,         
			     /*33:*/ sysTopMassUp,               /*34:*/ sysTopMassDown,    
			     /*35:*/ sysTopMassUp2,              /*36:*/ sysTopMassDown2,
			     /*37:*/ sysTopMassUp3,              /*38:*/ sysTopMassDown3,
			     /*39:*/ sysTopMassUp4,              /*40:*/ sysTopMassDown4,
			     /*41:*/ sysQCDUp,                   /*42:*/ sysQCDDown,                 
			     /*43:*/ sysSTopUp,                  /*44:*/ sysSTopDown,               
			     /*45:*/ sysDiBosUp,                 /*46:*/ sysDiBosDown,  
			     /*47:*/ sysVjetsUp,                 /*48:*/ sysVjetsDown, //##### ttV is added to these
			     /*49:*/ sysBRUp   ,                 /*50:*/ sysBRDown,
			     /*51:*/ sysPDFUp,                   /*52:*/ sysPDFDown,
			     /*53:*/ sysUnf,                     /*54:*/ sysMad,
			     /*55:*/ sysHadUp,                   /*56:*/ sysHadDown,			     
			     /*57:*/ sysGenMCatNLO,              /*58:*/ sysGenPowheg,
			     /*59:*/ sysGenPowhegHerwig,         /*60:*/ ENDOFSYSENUM, 
			     /*61:*/ sysTest, 			 /*62:*/ sysTestMCatNLO,             
			     /*63:*/ sysTestPowheg,     	 /*64:*/ sysTestPowhegHerwig,        
			     /*65:*/ ENDOFSYSENUM2,

  };

  // ============================
  //  Numerical Constants
  // ============================
  // remove official CMS in labels (needed for thesis)
  bool rmOff=false;
  // draw preliminary in labels for all plots?
  bool prelim=true;
  
  //const double ttbarCrossSection=234;                      // approx.NNNLO Kidonakis, recalculated for mtop=172.5 GeV (cf. TOP-11-008)
  //const double ttbarCrossSectionError=sqrt(12.*12.+((10.+7.)/2.)*((10.+7.)/2.)); // Scale and PDF uncertainties on NNLO value
                                                             // --> the scale contributions has been symetrized

  const double ttbarCrossSection=245.8;// NNNLO+NNLL Mitov, calculated for mtop=173.3 GeV (arXiv:1303.6254)
  const double ttbarCrossSectionError=sqrt(((6.2+8.4)/2.)*((6.2+8.4)/2.)+((6.2+6.4)/2.)*((6.2+6.4)/2.)); // Scale and PDF uncertainties 
                                                             // --> up/down contributions have been symetrized

  const double SF_TopMassDownUncertainty=MadSpin ? 1.0 : 0.9/3.0; // scale factors for top mass uncertainty
  const double SF_TopMassUpUncertainty  =MadSpin ? 1.0 : 0.9/3.0; // --> world average is presently known at +/-0.9 GeV (arXiv:1107.5255v3 [hep-ex])
                                                  // --> systematic samples are varied by +/-3.0 GeV 
                                                  // --> linearily rescale uncertainty on top mass in combineTopDiffXSecUncertainties.C

  const double constHadUncertainty   = 0.050; // relative uncertainty // outdated and only used as placeholder for bquark quantities
  const double globalLumiUncertainty = sqrt(0.025*0.025+0.005*0.005); // relative uncertainty // from CMS-PAS-LUM-13-001
	
  const double constLumiElec = 19712.0; // luminosity of Jan22ReReco ABCD dataset, see https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM, LUM-13-001
  const double constLumiMuon = 19712.0; // see above
  
  const double BRwlnu   =0.1080;
  const double BRwlnuUnc=0.0009;
  const double BRwqq    =0.6760;
  const double BRwqqUnc =0.0027;

  // ============== 
  //  Functions
  // ============== 

  double BRPDG(int sys=sysNo)
  {
    // This function returns the PDG BR for ttbar->e/mu+jets=2*[BR(W->lnu)*BR(W->qq)]
    // modified quantities: none
    // used quantities: BRwlnu, BRwlnuUnc, BRwqq, BRwqqUnc
    // used enumerators: systematicVariation
    double out = 2.;
    if(     sys==sysBRUp  ) out*=(BRwlnu+BRwlnuUnc);
    else if(sys==sysBRDown) out*=(BRwlnu-BRwlnuUnc);
    else out*=BRwlnu;
    if(     sys==sysBRUp  ) out*=(BRwqq+BRwqqUnc  );
    else if(sys==sysBRDown) out*=(BRwqq-BRwqqUnc  );
    else out*=BRwqq;
    return out;
  }

  double BRSF(int sys=sysNo, TString finalState="ljets")
  {
    // This function returns event SFs to get to the PDG 
    // BR for ttbar decays [BR(W->lnu)*BR(W->qq)]
    // modified quantities: none
    // used quantities: BRwlnu, BRwlnuUnc, BRwqq, BRwqqUnc
    // used enumerators: systematicVariation
    // finalState="ljets","alljets","dilepton"
    double out = 1.; 

    // BR W-> lnu
    double wlnuSF=3./(1./3.); // 3 leptons / (W->lnu BR used in simulation)
    if(     sys==sysBRUp  ) wlnuSF*=(BRwlnu+BRwlnuUnc);
    else if(sys==sysBRDown) wlnuSF*=(BRwlnu-BRwlnuUnc);
    else wlnuSF*=BRwlnu;
    // BR W->qq
    double wqqSF=1./(2./3.); // 1. / (W->all qqbar BR used in simulation)
    if(     sys==sysBRUp  ) wqqSF*=(BRwqq+BRwqqUnc  );
    else if(sys==sysBRDown) wqqSF*=(BRwqq-BRwqqUnc  );
    else wqqSF*=BRwqq;
    // SF for final state 
    if(     finalState=="ljets"   ) out*=(wlnuSF*wqqSF );
    else if(finalState=="alljets" ) out*=(wqqSF *wqqSF );
    else if(finalState=="dilepton") out*=(wlnuSF*wlnuSF);
    else{
      std::cout << "ERROR in BRSF: unknow input finalState=" << finalState << std::endl;
      exit(0);
    }
    return out;
  }

  double BRSFBG(int sys=sysNo, double fracLjets=0.6, double fracDilep=0.4, double fracAllHad=0.){
    // This function returns event SFs to get to the PDG 
    // BR for ttbar a mixed sample of ttbar events specified
    // via fracLjets, fracDilep & fracAllHad
    // modified quantities: none
    // used functions: BRSF
    // used enumerators: systematicVariation
    // NOTE: the 60%l+jets and 40%dilepton composition 
    // corresponds to the 8TeV doubleKinFit+prob analysis
    return fracLjets*BRSF(sys,"ljets")+fracDilep*BRSF(sys,"dilepton")+fracAllHad*BRSF(sys,"alljets");
  }

  TString sysLabel(unsigned int sys)
  {
    // this function returns a TString that corresponds
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // modified quantities: none
    // used functions: none
    // used enumerators: none (label correspond to systematicVariation)

    switch (sys)
    {
      case sysNo                       : return "sysNo";
      case sysLumiUp                   : return "sysLumiUp";
      case sysLumiDown                 : return "sysLumiDown";                
      case sysPUUp                     : return "sysPUUp";                  
      case sysPUDown                   : return "sysPUDown";                  
      case sysJESUp                    : return "sysJESUp";
      case sysJESDown                  : return "sysJESDown";
      case sysJERUp                    : return "sysJERUp";                 
      case sysJERDown                  : return "sysJERDown";                
      case sysLepEffSFNormUp           : return "sysLepEffSFNormUp";
      case sysLepEffSFNormDown         : return "sysLepEffSFNormDown";
      case sysLepEffSFShapeUpEta       : return "sysLepEffSFShapeUpEta";
      case sysLepEffSFShapeDownEta     : return "sysLepEffSFShapeDownEta";
      case sysLepEffSFShapeUpPt        : return "sysLepEffSFShapeUpPt";
      case sysLepEffSFShapeDownPt      : return "sysLepEffSFShapeDownPt"; 
   // case sysLepEffSFShapeUpPt40      : return "sysLepEffSFShapeUpPt40";
   // case sysLepEffSFShapeDownPt40    : return "sysLepEffSFShapeDownPt40";
   // case sysLepEffSFNormUpStat       : return "sysLepEffSFNormUpStat";
   // case sysLepEffSFNormDownStat     : return "sysLepEffSFNormDownStat";
   // case sysTriggerEffSFJetNormUp    : return "sysTriggerEffSFJetNormUp";
   // case sysTriggerEffSFJetNormDown  : return "sysTriggerEffSFJetNormDown";
   // case sysTriggerEffSFJetShapeUp   : return "sysTriggerEffSFJetShapeUp";
   // case sysTriggerEffSFJetShapeDown : return "sysTriggerEffSFJetShapeDown"; 
      case sysBtagSFUp                 : return "sysBtagSFUp";
      case sysBtagSFDown               : return "sysBtagSFDown";
   // case sysBtagSFShapeUpPt100       : return "sysBtagSFShapeUpPt100";
   // case sysBtagSFShapeDownPt100     : return "sysBtagSFShapeDownPt100";
   // case sysBtagSFShapeUpEta1p2      : return "sysBtagSFShapeUpEta1p2";
   // case sysBtagSFShapeDownEta1p2    : return "sysBtagSFShapeDownEta1p2";
      case sysBtagSFShapeUpPt65        : return "sysBtagSFShapeUpPt65";
      case sysBtagSFShapeDownPt65      : return "sysBtagSFShapeDownPt65";
      case sysBtagSFShapeUpEta0p7      : return "sysBtagSFShapeUpEta0p7";
      case sysBtagSFShapeDownEta0p7    : return "sysBtagSFShapeDownEta0p7"; 
   // case sysBTagSFFullShapeUpPt65    : return "sysBTagSFFullShapeUpPt65";
   // case sysBTagSFFullShapeDownPt65  : return "sysBTagSFFullShapeDownPt65";
   // case sysBTagSFFullShapeUpEta0p7  : return "sysBTagSFFullShapeUpEta0p7";
   // case sysBTagSFFullShapeDownEta0p7: return "sysBTagSFFullShapeDownEta0p7"; 
      case sysMisTagSFUp               : return "sysMisTagSFUp";
      case sysMisTagSFDown             : return "sysMisTagSFDown";
      case sysTopScaleUp               : return "sysTopScaleUp";
      case sysTopScaleDown             : return "sysTopScaleDown";
      case sysVBosonScaleUp            : return "sysVBosonScaleUp";
      case sysVBosonScaleDown          : return "sysVBosonScaleDown";
      case sysSingleTopScaleUp         : return "sysSingleTopScaleUp";
      case sysSingleTopScaleDown       : return "sysSingleTopScaleDown";
      case sysTopMatchUp               : return "sysTopMatchUp";
      case sysTopMatchDown             : return "sysTopMatchDown";            
      case sysVBosonMatchUp            : return "sysVBosonMatchUp";          
      case sysVBosonMatchDown          : return "sysVBosonMatchDown";
      case sysTopMassUp                : return "sysTopMassUp";  
      case sysTopMassDown              : return "sysTopMassDown";   
      case sysTopMassUp2               : return "sysTopMassUp2";  
      case sysTopMassDown2             : return "sysTopMassDown2";   
      case sysTopMassUp3               : return "sysTopMassUp3";  
      case sysTopMassDown3             : return "sysTopMassDown3";   
      case sysTopMassUp4               : return "sysTopMassUp4";  
      case sysTopMassDown4             : return "sysTopMassDown4";   
      case sysQCDUp                    : return "sysQCDUp";           
      case sysQCDDown                  : return "sysQCDDown";
      case sysSTopUp                   : return "sysSTopUp";               
      case sysSTopDown                 : return "sysSTopDown";
      case sysDiBosUp                  : return "sysDiBosUp";              
      case sysDiBosDown                : return "sysDiBosDown";
      case sysVjetsUp                  : return "sysVjetsUp";      
      case sysVjetsDown                : return "sysVjetsDown";      
      case sysBRUp                     : return "sysBRUp";              
      case sysBRDown                   : return "sysBRDown";    
      case sysPDFUp                    : return "sysPDFUp";              
      case sysPDFDown                  : return "sysPDFDown"; 
      case sysUnf                      : return "sysUnf";   
      case sysMad                      : return "sysMad";    
      case sysHadUp                    : return "sysHadronizationUp";
      case sysHadDown                  : return "sysHadronizationDown";
      case sysGenMCatNLO               : return "sysGenMCatNLO";
      case sysGenPowheg                : return "sysGenPowheg";  
      case sysGenPowhegHerwig          : return "sysGenPowhegHerwig";  
      case sysTest                     : return "sysTest";  
      default                          : std::cout << "ERROR: the chosen input for function sysLabel is not valid" << std::endl;
                                         std::cout << "chosen systematic variation:  " << sys            << std::endl;
                                         std::cout << "maximum systematic variation: " << ENDOFSYSENUM-1 << std::endl;
                                         return "";
    }
  }

  TString sampleLabel(unsigned int sample, const std::string decayChannel)
  {
    // this function returns the name of the entered MC process
    // corresponding to the enumerator entry "sample" as defined in samples
    // modified quantities: NONE
    // used functions:      NONE
    // used enumerators:    sample
    
    // list all MC process/sample names
    if(sample==kSig      ) return "t#bar{t} Signal";
    if(sample==kSigPow   ) return "t#bar{t} Signal Powheg+Pythia";
    if(sample==kSigPowHer) return "t#bar{t} Signal Powheg+Herwig";
    if(sample==kSigMca   ) return "t#bar{t} Signal MC@NLO+Herwig";
    if(sample==kBkg       ) return "t#bar{t} Other";
    if(sample==kBkgPow    ) return "t#bar{t} Other Powheg+Pythia";
    if(sample==kBkgPowHer ) return "t#bar{t} Other Powheg+Herwig";
    if(sample==kBkgMca    ) return "t#bar{t} Other MC@NLO+Herwig";
    if(sample==kSTop   ) return "Single Top";
    if(sample==kSToptW ) return "Single Top tW";
    if(sample==kSTops  ) return "Single Top s";
    if(sample==kSTopt  ) return "Single Top t";
    if(sample==kSAToptW) return "Single Antitop tW";
    if(sample==kSATops ) return "Single Antitop s";
    if(sample==kSATopt ) return "Single Antitop t";
    if(sample==kSToptW1 ) return "Single Top tW dilepton";
    if(sample==kSAToptW1) return "Single Antitop tW dilepton";
    if(sample==kSToptW1 ) return "Single Top tW t->had W->lep";
    if(sample==kSAToptW1) return "Single Antitop tW t->had W->lep";
    if(sample==kSToptW1 ) return "Single Top tW t->lep W->had";
    if(sample==kSAToptW1) return "Single Antitop tW t->lep W->had";
    if(sample==kWjets  ) return "W+Jets"; 
    if(sample==kTTGjets  ) return "t#bar{t}+#gamma"; //###############
    if(sample==kTTZjets  ) return "t#bar{t}+Z"; //###############
    if(sample==kTTWjets  ) return "t#bar{t}+W"; //###############
    if(sample==kTTVjets  ) return "t#bar{t}+V"; //###############
    if(sample==kZjets  ) return "Z / #gamma_{#lower[-0.7]{*}}+Jets";
    if(sample==kDiBos  ) return "Diboson";
    if(sample==kWW     ) return "WW";
    if(sample==kWZ     ) return "WZ";
    if(sample==kZZ     ) return "ZZ";
    if(sample==kQCD    ) return "QCD Multijet";
    if(sample==kQCDEM1 ) return "QCDEM1";
    if(sample==kQCDEM2 ) return "QCDEM2";
    if(sample==kQCDEM3 ) return "QCDEM3";
    if(sample==kQCDBCE1) return "QCDBCE1"; 
    if(sample==kQCDBCE2) return "QCDBCE2";
    if(sample==kQCDBCE3) return "QCDBCE3"; 
    if(sample==kData   ) return "Data";
    else return "Default "+decayChannel;
  }

  TString extractxSecName(TString in)
  {
    // this function returns a TString that corresponds
    // to the xSec quantity contained in "in", none is 
    // returned if "in" does not contain any quantities
    // modified quantities: none
    // used functions: none
    // used enumerators: none 

    // define output
    TString out="none";
    // parton level full PS quantities
    if(     in.Contains("topPtSubLead" )) out="topPtSubLead" ;
    else if(in.Contains("topPtLead"    )) out="topPtLead"    ;
    else if(in.Contains("topPtTtbarSys")) out="topPtTtbarSys";
    else if(in.Contains("topPt"        )) out="topPt"        ;
    else if(in.Contains("topY"         )) out="topY"         ;
    else if(in.Contains("ttbarPt"      )) out="ttbarPt"      ;
    else if(in.Contains("ttbarY"       )) out="ttbarY"       ;
    else if(in.Contains("ttbarMass"    )) out="ttbarMass"    ;
    else if(in.Contains("ttbarDelPhi"  )) out="ttbarDelPhi"  ;
    else if(in.Contains("ttbarPhiStar" )) out="ttbarPhiStar" ;
    // visible hadron/particle level quantities
    else if(in.Contains("lepPt"    )) out="lepPt";
    else if(in.Contains("lepEta"   )) out="lepEta";
    else if(in.Contains("bqPt"     )) out="bqPt";
    else if(in.Contains("bqEta"    )) out="bqEta";
    else if(in.Contains("bbbarMass")) out="bbbarMass";
    else if(in.Contains("bbbarPt"  )) out="bbbarPt";
    else if(in.Contains("lbMass"   )) out="lbMass";
    else if(in.Contains("rho"      )) out="rhos";
    else if(in.Contains("Njets")||in.Contains("Ngenjets")) out="Njets";
    // return output
    return out;
  }



  double effSFAB(int sys=sysNo, std::string decayChannel="unset")
  {
    // this function returns the muon eff SF
    // as derived from Z->mumu tag and probe method
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    // "sys": if sysMuEffSFup/sysMuEffSFdown is
    // chosen, the corresponding systematic shifted
    // SF is returned

    if(sys>=ENDOFSYSENUM2){
      std::cout << "ERROR in effSFAB:" << std::endl;
      std::cout << "sys must be smaller than " << ENDOFSYSENUM << std::endl;
    }

    // SF is now applied as event weight in analyzer
    //return 1.0;

    // combined single muon SF Run A+B from tag and probe
    double result = -1.;
    if (decayChannel.compare("muon")==0) {
      result = 1.0; // SF is now applied as event weight in analyzer
    }
    else if (decayChannel.compare("electron")==0) {
      // result = 0.9658; 
      result = 1.0; // SF is now applied as event weight in analyzer
    }
    else if (decayChannel.compare("combined")==0) result = 1.0;
    // errors for the derived SF
    //double errorUp   = 0.03*result;
    //double errorDown = 0.03*result;
    // return requestet SF
    // SF uncertainty is now applied as event weight in analyzer
    //if(sys==sysMuEffSFUp  ) result+=errorUp;
    //if(sys==sysMuEffSFDown) result-=errorDown;
    return result;
  }
  
  void histogramStyle(TH1& hist, int sampleType, bool filled=true, double markersize=1.2, unsigned int color=0)
  {
    // this function configures the style of a TH1 histogram "hist"
    // using "sampleTyp" to identify the corresponding sample from
    // the enumerator samples definition and assign the corresponding
    // marker and color
    // modified quantities: hist
    // used functions: marker_, color_
    // used enumerators: samples
    // "filled": 0 - line only
    //           1 - area under plot filled

    // identify subsamples with main sample
    if(sampleType>=kSTops&&sampleType<=kSAToptW ) sampleType=kSTop;
    if(sampleType>=kQCDEM1&&sampleType<=kQCDBCE3) sampleType=kQCD;
    if(sampleType>=kTTGjets&&sampleType<=kTTWjets) sampleType=kTTVjets; //######

    hist.SetStats(kFALSE);
    if(sampleType==kData || !filled){
      if(!filled)hist.SetLineWidth(2);
      hist.SetLineColor(color_[sampleType]);
      hist.SetMarkerColor(color_[sampleType]);
      if(sampleType==kQCD){
	hist.SetLineColor(kYellow+1);
	hist.SetMarkerColor(kYellow+1);
      }
      if(sampleType==kDiBos){
	hist.SetLineColor(kGray);
	hist.SetMarkerColor(kGray);
      }
      hist.SetMarkerStyle(marker_[sampleType]);
      hist.SetMarkerSize(markersize);
      hist.SetFillStyle(0);
    }
    else{
      hist.SetFillColor(color_[sampleType]);
    }
    if(color!=0){
      hist.SetLineColor(color);
      hist.SetFillColor(color);
      hist.SetMarkerColor(color);
    }
  }

  void histStyle2D(TH2& hist, const TString titleHisto, const TString titleX, const TString titleY)
  {
    // this function configures the style of a TH2 histogram "hist"
    // modified quantities: hist
    // used functions: NONE
    // used enumerators: NONE
    // "titleHisto": set tite of the histo (drawn in upper left)
    // "titleX": title of the x-axis
    // "titleY": title of the y-axis

    hist.SetTitle(titleHisto);
    if(titleX!="") hist.GetXaxis()->SetTitle(titleX);
    if(titleY!="") hist.GetYaxis()->SetTitle(titleY);
    hist.SetStats(kFALSE);
  }

  double readValFromLineFromFile(int line, TString file="crossSectionCalculation.txt")
  {
    // this function to reads and returns a double value
    // from a specific line "line" of the .txt-like file "file"
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE

    // define variables
    std::ifstream finDouble (file);
    std::string readIn;
    // check if file exists
    if (!finDouble){
      std::cout << "can not open file " << file << std::endl;
      return -1;
    }
    // loop lines of the file
    for(int l=1; !finDouble.eof(); ++l){
      // save line content in readIn
      getline(finDouble, readIn);
      // convert your chosen line into double and return it
      if(l==line) return atof(readIn.c_str());
    }
    // if line is not found
    std::cout << "can not find line" << line << std::endl;
    return -1.;
  }

  TString readLineFromFile(int line, TString file="crossSectionCalculation.txt")
  {
    // this function to reads and returns a TString expression
    // from a specific line "line" of the .txt-like file "file"
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE

    // define variables
    std::ifstream finDouble (file);
    std::string readIn;
    // check if file exists
    if (!finDouble){
      std::cout << "can not open file " << file << std::endl;
      return -1;
    }
    // loop lines of the file
    for(int l=1; !finDouble.eof(); ++l){
      // save line content in readIn
      getline(finDouble, readIn);
      // convert your chosen line into double and return it
      if(l==line) return TString(readIn.c_str());
    }
    // if line is not found
    std::cout << "can not find line" << line << std::endl;
    return "";
  }

  //void canvasStyle(TCanvas& canv) re-definition of header to avoid warnings when compiling, must be changed once function body is not empty anymore
  void canvasStyle()
  {
    // function is not called anymore - just kept if indiviudal configurations might be required
    // could be considered to be deleted
  }

  void axesStyle(TH1& hist, const char* titleX, const char* titleY, float yMin=-123, float yMax=-123)
  {
    // this function configures the axis style of a TH1 histogram "hist"
    // modified quantities: hist
    // used functions: NONE
    // used enumerators: NONE
    // "titleX": title of the x-axis
    // "titleY": title of the y-axis
    // "yMin": minimum of the y axis
    // "yMax": maximum of the y axis

    hist.SetTitle("");
    hist.GetXaxis()->SetTitle(titleX);
    hist.GetYaxis()->SetTitle(titleY);
    if(yMin!=-123) hist.SetMinimum(yMin);
    if(yMax!=-123) hist.SetMaximum(yMax);
  }

  void legendStyle(TLegend& leg, TString header="", double textsize=0.03, unsigned int aligncode=12)
  {
    // this function configures the style of legends
    // modified quantities: leg
    // used functions: NONE
    // used enumerators: NONE
    // "header": header of the legend

    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(textsize);
    leg.SetTextAlign(aligncode);
    if (header!="") leg.SetHeader(header); // setting the header adds additional row even if contents is empty (""), important when using TLegend::GetNRows()
  }

  template <class T>
    void writeToFile(T output, TString file, bool append)
    {
      // this function writes the content of "output"
      // into "file" which is of a type like .txt
      // "output" might be of type TString or double
      // modified quantities: file
      // used functions: NONE
      // used enumerators: NONE
      // "append": if false, "file" will be recreated and
      // the existing "file" will be deleted

      // a) write into file
      if(!append){
	std::ofstream fout(file);
	fout << output << std::endl;
	fout.close();
      }
      // b) add output to the end of the file
      if(append){
	std::ofstream fapp(file, std::ios::app);
	fapp << output << std::endl;
	fapp.close();
      }
    }

  void drawLine(const double xmin, const double ymin, const double xmax, const double ymax, const unsigned int color=kBlack, const double lineWidth=2, const unsigned int lineStyle=1)
  {
    // this function draws a line withe the chosen coordinates,
    // color and width into the active canvas
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    TLine *line = new TLine();
    line->SetLineWidth(lineWidth);
    line->SetLineStyle(lineStyle);
    line->SetLineColor(color);
    line->DrawLine(xmin, ymin, xmax, ymax);
  }

  void DrawSteps(TH1F* histo, TString drawOptions=""){
    // this function draws "histo" in the
    // active canvas like hist but without 
    // any vertical lines
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    
    // adjust markers to be invisible wrt line
    histo->SetMarkerSize(0);
    histo->SetMarkerColor(histo->GetLineColor());
    // remove errors from markers
    for(int bin=1; bin<=histo->GetNbinsX(); ++bin){
      histo->SetBinError(bin, 1E-9);
    }
    // draw not as histo
    drawOptions.ReplaceAll("hist","");
    // draw as marker!
    drawOptions+=" E";
    //std::cout << drawOptions << std::endl;
    // ensure that x range is drawn around marker
    float_t originalValue=gStyle->GetErrorX();
    if(originalValue!=0.5){
      std::cout << "WARNING in DrawSteps: will change";
      std::cout << " gStyle->SetErrorX from " << originalValue;
      std::cout << " to 0.5 to make DrawSteps working!!!" << std::endl; 
      gStyle->SetErrorX(0.5);
    }
    //drawing
    histo->Draw(drawOptions);
    // restore original configuration
    //gStyle->SetErrorX(originalValue);
  }
  
  int roundToInt(double value, bool roundDown=false)
  {
    // function to round a double "value"
    // to an int and return this one
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    // "roundDown": choose if you always want to round down

    // create output
    int outputInt=0;
    // take care of negative numbers
    if(value<0){
      std::cout << "no negative numbers allowed in roundToInt" << std::endl;
      return 0;
    }
    // get number before the digit
    for(int x=0; value>x; ++x){
      outputInt=x;
    }
    // see if the number behind the digit is > 0.5
    if((roundDown==false)&&(value-outputInt) >=0.5) return (outputInt+1);
    // return result
    return outputInt;
  }

  Double_t round(Double_t value, Int_t k)
  {
    // function to round a double "value"
    // to a double with "k" digits, k>0

    // no need for rounding if the value is 0
    if(value==0.)return value;
    if(k==0) return value;
    Double_t result = value;
    // search window: 10^20 to 10^-20
    for(Int_t i = 20; i > -20; i--){
      Double_t exp = pow(10.,i);
      // searching for the 'k'th non-zero digit
      if(round(value/(exp*pow(10.,k-1))) != 0){
	// round to the found digit
	result = round(value/(exp))*exp;
	break;
      }
    }
    return result;
  }
  
  TString getTStringFromInt(int i)
  {
    // function to convert an int "i" to
    // a TString and return this one
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    char result[20];
    sprintf(result, "%i", i);
    return (TString)result;
  }

  int getRelevantDigits(double value){
    // this function calculates the number 
    // of digits of "value"+1 if following digit is nonzero
    // ATTENTION: only the first 10 digits are checked!
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    if(value==0.) return 0.;
    double modvalue=std::abs(value);
    for(int i=0; i<10; ++i){
      if(modvalue>=1){
	double returnvalue=i;
	// check if next digit is nonzero
	if(10*(modvalue-int(modvalue))>=1) returnvalue+=1;
	return returnvalue;
      }
      else modvalue*=10;
    }
    return -1;
  }

  TString getTStringFromDouble(double d, int precision=2, bool output=false)
  {
    // function to convert an double "d" to
    // a TString and return this one
    // precision: number of digits to be considered
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    TString conv="%.";
    conv+=getTStringFromInt(precision);
    conv+="f";
    char result[30];
    if(output) std::cout << "sprintf(result, conv, d)=sprintf(" << result << ", " << conv << ", " << d << ")" << std::endl;
    sprintf(result, conv, d);
    return (TString)result;
  }

  TH1F* divideByBinwidth(TH1F* histo, bool calculateError=true)
    {
      // function divides the #entries in every bin of the input plot "histo"
      // by its binwidth and returns the result
      // the errors are recalculated if "calculateError"=1 is chosen
      // careful: not done for underflow/overflow
      // modified quantities: NONE
      // used functions: NONE
      // used enumerators: NONE

      // create output histo
      TH1F* outputHisto = (TH1F*)histo->Clone();
      // loop bins, excluding underflow and overflow
      for(int i=1; i<= histo->GetNbinsX(); i++){
	double binwidth=(double)(histo->GetBinWidth(i));
	// take care of bins width width 0
	// (these could not be filled and
	// exist because of technical reasons only)
	if(binwidth==0)binwidth=1;
	// recalculate errors
	double error = (double)(histo->GetBinError(i)) / binwidth;
	// recalculate bin content
	outputHisto->SetBinContent(i, ( (double)(histo->GetBinContent(i)) / binwidth) );
	if(calculateError) outputHisto->SetBinError(i, error);
      }
      // return result
      return outputHisto;
    }

  int getBigitFromDouble(double d)
  {
    // function to determine the number of 
    // digits from d before the comma, 
    // for negative numbers the minus sign 
    // is counted too
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    
    int out=1;
    while(std::abs(d)>=pow(10, out)) out++;
    if(d<0.) out++;
    return out;

  }

  TString fillspace(double number, int lengthref=3){
    // function to determine the number of 
    // spaces of number wrt. lengthref
    // another space is added for positive values
    // -> used to make uniform indented output, e.g.
    // -0.5 to  0.0
    //  0.0 to 10.5
    // modified quantities: NONE
    // used functions: getBigitFromDouble
    // used enumerators: NONE
    TString out="";   
    int lenght=getBigitFromDouble(number);
    if(number<0) lenght-=1;
    //std::cout << "fillspace(" << number << "->" << lenght << "," << lengthref << ")" << std::endl;
    if(lenght<lengthref){
      for(int space=(lengthref-lenght); space>0; space--){ out+=" "; }
    }
    if(number>=0){ out+=" "; } 
    return out;
  }

  TString fillspaceT(TString label, int lengthref=3){
    // as fillspace but for TStrings instead of doubles
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    TString out="";   
    int lenght=label.Length();
    //std::cout << "fillspaceT(" << label << "->" << lenght << "," << lengthref << ")" << std::endl;
    if(lenght<lengthref){
      for(int space=(lengthref-lenght); space>0; space--){ out+=" "; }
    }
    return out;
  }

  void DrawDecayChLabel(TString decaychannel="", double textSize=0.04)
  {
    // Draw label for Decay Channel in upper left corner of plot

    TPaveText *decch = new TPaveText();

    decch -> AddText(decaychannel);

    decch -> SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch -> SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch -> SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch -> SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch -> SetFillStyle(0);
    decch -> SetBorderSize(0);
    decch -> SetTextSize(textSize);
    decch -> SetTextAlign(12);
    decch -> Draw("same");
  }

  void DrawCMSLabels(bool cmsprelim=true, double luminosity=0.0, double textSize=0.04, bool cmssimulation=false, bool privateWork=false, bool sevenTeV=false)
  {
    // Draw official labels (CMS Preliminary/Simulation, luminosity and CM energy) above plot
    bool rmOff2= PHD ? true : rmOff;

    TPaveText *label = new TPaveText();

    bool cmsprelim2= (useGlobalPrelim) ? globalPrelim : cmsprelim;

    label -> SetX1NDC(gStyle->GetPadLeftMargin());
    label -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label -> SetY2NDC(1.0);
    label -> SetTextFont(42);

    TString extension="CMS";
    TString comE="8";
    if(sevenTeV) comE="7";

    if(cmssimulation) extension+=" Simulation";
    if(cmsprelim2&&!PHD) extension+=" Preliminary";
    if(privateWork) extension+=" (Private Work)";
    
    TString finalLabel="";
    if (cmssimulation||luminosity==0.){
      if(!rmOff2) finalLabel=extension+", #sqrt{s} = "+comE+" TeV";
      else  finalLabel="Simulation at  #sqrt{s} = "+comE+" TeV";
    }
    else{
      if(!rmOff2) finalLabel=TString(Form(extension+", %2.1f fb^{-1} at #sqrt{s} = "+comE+" TeV",luminosity/1000));
      else finalLabel=TString(Form("%2.1f fb^{-1} at #sqrt{s} = "+comE+" TeV",luminosity/1000));
    }

    
    label-> AddText(finalLabel);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
  }

  void DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering=12, double textSize=0.04)
  {
    // function to directly draw a label into the active canvas
    // the coordinates of the window are: "x1", "y1", "x2", "y2"
    // the label will contain the text "text" and the textsize is "textSize"
    // modified quantities: the active canvas
    // used functions: NONE
    // used enumerators: NONE

    TPaveText *label = new TPaveText(x1, y1, x2, y2, "br NDC");
    label->AddText(text);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if(textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(centering);
    label->Draw("same");
  }

  void scaleByFactor(TH1F*& histo, const double scaleValue)
  {
    // function to scale the TH1F plot "histo" and its errors
    // by an entered "scaleValue" correctly for each bin
    // N'=N/L, sN'=sN/L
    // modified quantities: histo
    // used functions: NONE
    // used enumerators: NONE

    // scaling factor = 0?
    if(scaleValue==0) std::cout << "WARNING: scaling histogram by 0" << std::endl;
    // loop bins (including underflow and overflow)
    for(int bin=0; bin<=histo->GetNbinsX()+1; ++bin){
      double oldError = histo->GetBinError(bin);
      // scale bin value and error
      histo->SetBinContent(bin, histo->GetBinContent(bin)/scaleValue);
      histo->SetBinError  (bin, oldError                 /scaleValue);
    }
  }

  double lumiweight(unsigned int sample, double luminosity, unsigned int kSys, const std::string decayChannel, bool MS=MadSpin)
  {
    // this function derives the lumiweight for every standard MC
    // sample "sample" based on the theoretical cross section, the
    // number of generated events and the chosen "luminosity"
    // Furthermore, the BR correction is considered for ttbarSG/BG
    // NOTE: enter luminosity IN / pb!!!!
    // modified quantities: NONE
    // used functions: BRSFBG, BRSF
    // used enumerators: samples, systematicVariation
    TString calledfunction="lumiweight(unsigned int sample=";
    calledfunction+=sample;
    calledfunction+=", double luminosity=";
    calledfunction+=luminosity;
    calledfunction+=", unsigned int kSys=";
    calledfunction+=kSys;
    calledfunction+=", const std::string decayChannel=";
    calledfunction+=decayChannel;
    calledfunction+=")";

    // large MadSpin sample?
    bool MS2= (((MS&&kSys!=sysMad)||(!MS&&kSys==sysMad))&&((sample==kSig)||(sample==kBkg))&&kSys<=sysHadDown)  ? true : false;

    // function internal detail level of text output
    int verbose=0;
    // a) check if input is valid
    // sample existing?
    if(sample>=ENDOFSAMPLEENUM2||sample==ENDOFSAMPLEENUM){
      std::cout << "ERROR: invalid sample label for lumiweight calculation, no scaling will be done" << std::endl;
      std::cout << "       ( " << calledfunction << " )" << std::endl;
      return 1.;
    }
    // valid luminosity entered?
    if(luminosity<=0){
      std::cout << "ERROR: chosen luminosity for lumiweight calculation is <= 0, no scaling will be done" << std::endl;
      std::cout << "       ( " << calledfunction << " )" << std::endl;
      return 1.;
    }
    // valid decayChannel entered?
    if(!(decayChannel.compare("muon")==0 || decayChannel.compare("electron")==0)){
      // for combined channel: no weighting but also no text output
      if(!(decayChannel.compare("combined")==0)){ 
	std::cout << "ERROR: chosen decaychannel does not correspond to electron or muon, no scaling will be done" << std::endl;
	std::cout << "       ( " << calledfunction << " )" << std::endl;
      }
      return 1.;
    }
    // b) define in/output for weight calculation
    double crossSection=1;
    double Nevents=0;
    double weight =0;

    // c) fill #events in sample and cross sections (in / pb)
    //
    // The second number always indicates number without kinematics filter for non-Pythia samples
    //
    //  MadGraph: ttbar signal + background
    if((sample==kSig)||(sample==kBkg)){
      crossSection=ttbarCrossSection; 
      // Summer12
      Nevents = MS2 ? 62131965 : 6923750;
      // Systematic samples
      if(kSys!=sysNo){
	if     (kSys==sysTopScaleUp  ) Nevents= MS2 ?  41908271 : 5009488;
	else if(kSys==sysTopScaleDown) Nevents= MS2 ?  39286663 : 5387181;
	else if(kSys==sysTopMatchUp  ) Nevents= MS2 ?  37083003 : 5415010;
	else if(kSys==sysTopMatchDown) Nevents= MS2 ?  34053113 : 5476728;
	else if(kSys==sysTopMassUp   ) Nevents= MS2 ?  26489020 : 5249525;
	else if(kSys==sysTopMassDown ) Nevents= MS2 ?  24439341 : 5369214;
	else if(kSys==sysTopMassUp2  ) Nevents= MS2 ?  40244328 : 4733483;
	else if(kSys==sysTopMassDown2) Nevents= MS2 ?  39518234 : 4469095;
	else if(kSys==sysTopMassUp3  ) Nevents= MS2 ?  24359161 : 5145140;
	else if(kSys==sysTopMassDown3) Nevents= MS2 ?  27078777 : 5365348;
	else if(kSys==sysTopMassUp4  ) Nevents= MS2 ? -1 : 5249525;
	else if(kSys==sysTopMassDown4) Nevents= MS2 ? -1 : 5369214;
      }
    }
    else if((sample==kSigPow)||(sample==kBkgPow)){
      crossSection=ttbarCrossSection; 
      // Summer12
      Nevents = 21675970*0.999886;// 8 TEV PU SF from gen ttbarSG Integral/entries
    }
    else if((sample==kSigPowHer)||(sample==kBkgPowHer)){
      crossSection=ttbarCrossSection; 
      // Summer12
      Nevents = 27684235*0.999985;// 8 TEV PU SF from gen ttbarSG Integral/entries
    }
    else if((sample==kSigMca)||(sample==kBkgMca)){      crossSection=ttbarCrossSection; 
      // Summer12
      // MC@NLO: total events is sum of weights
      //    21745196 total     weights
      //    19277872 positive  weights (88.65%)
      //     2467324 negative  weights (11.35%)
      //    16810548 effective weights
      // each |weight| = 190.41256 -> number directly from LHE
      Nevents = 32852589.*210.9532; // 8 TEV convoluted PU*weight SF from gen ttbarSG Integral/entries
    }
    // MadGraph: W->lnu+jets
    else if(sample==kWjets){
      crossSection=36257;
      // Summer12
      Nevents = 57709905;
      // Systematic samples:
      if(kSys==sysVBosonScaleUp  ) Nevents=20784770; 
      if(kSys==sysVBosonScaleDown) Nevents=20760884;  
      if(kSys==sysVBosonMatchUp  ) Nevents=20976082;
      if(kSys==sysVBosonMatchDown) Nevents=21364637;  
    }
    //###########################
    else if(sample==kTTGjets){ 
      crossSection=1.8;
      // Summer12
      Nevents = 1719954;
    }
    else if(sample==kTTZjets){
      crossSection=0.2057;
      // Summer12
      Nevents = 210160;
    }
    else if(sample==kTTWjets){
      crossSection=0.232;
      Nevents = 196046;
    } 
    //###########################
    // MadGraph: DY->ll+jets
    else if(sample==kZjets){
      crossSection=3503;
      // Summer12
      Nevents =30459503;
      // Summer12 systematic samples:
      if(kSys==sysVBosonScaleUp  ) Nevents=2170270;
      if(kSys==sysVBosonScaleDown) Nevents=1934901;
      if(kSys==sysVBosonMatchUp  ) Nevents=1985529;
      if(kSys==sysVBosonMatchDown) Nevents=2112387;
    }
    // Pythia6: QCD Mu enriched
    else if(sample==kQCD&&decayChannel.compare("muon")==0){
      // Summer12
      crossSection=3.64E8*3.7E-4; // generator crossSection * prefilter efficiency
      Nevents     =21484602;
    }
    // QCD e+jets channel 
    // Summer12
    // a) subsamples
    else if(sample==kQCDEM1){
      crossSection= 2.886E8*1.01E-2;  // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 35040695;
    }
    else if(sample==kQCDEM2){
      crossSection= 7.433E7*6.21E-2;    // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 33088888;
    }
    else if(sample==kQCDEM3){
      crossSection= 1.191E6*0.1539;      // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 34542763;
    } //EM4: 31697066 3.099E4*0.148 //EM5: 34611322 4.25E3*0.131 //EM6: 34080562 8.1E2*0.11
    else if(sample==kQCDBCE1){
      crossSection= 2.886E8*5.8E-4; // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 1740229;
    }
    else if(sample==kQCDBCE2){
      crossSection= 7.424E7*2.25E-3;  // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 2048152; 
    }
    else if(sample==kQCDBCE3){
      crossSection= 1.191E6*1.09E-2;     // generator prefilter efficiency * LO PYTHIA crossSection
      Nevents     = 1945525; 
    } //BCE4: 1948112 3.098E4*2.04E-2 //BCE5: 2026521 4.25E3*2.43E-2  //BCE6: 1948532 8.11E2*2.95E-2 
    // b)combined
    else if(sample==kQCD&&decayChannel.compare("electron")==0){
      // Summer12: already added in combineMCsamples.C
      // with cross section as weight, only
      // lumi normalization is done here
      Nevents     =1.;
      crossSection=1;
    }
    // Powheg: singleTop/Antitop Samples
    // a) subsamples
    else if(sample==kSATops){
      crossSection= 1.8;
      Nevents     = 139974;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =999812;
      if(kSys==sysTopScaleDown) Nevents =999849;
    }
    else if(sample==kSTops){
      crossSection= 3.8;
      Nevents     = 259961;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =1998845;
      if(kSys==sysTopScaleDown) Nevents =1945749;
    }
    else if(sample==kSATopt){
      crossSection= 30.7;
      Nevents     = 1935072;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =979898;
      if(kSys==sysTopScaleDown) Nevents =979359;
    }
    else if(sample==kSTopt){
      crossSection= 56.4;
      Nevents     = 3758227;
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =1945116;
      if(kSys==sysTopScaleDown) Nevents =1951907;
    }
    else if(sample==kSAToptW){
      crossSection= 11.1;
      Nevents     = 493460;
      // scale variation
      // if(kSys==sysTopScaleUp  ) Nevents =; // not existing for 8 TEV
      // if(kSys==sysTopScaleDown) Nevents =; // not existing for 8 TEV
    }
    else if(sample==kSToptW){
      crossSection= 11.1; 
      Nevents     = 497658;
      // scale variation
      // if(kSys==sysTopScaleUp  ) Nevents =; // not existing 8 TEV
      // if(kSys==sysTopScaleDown) Nevents =; // not existing 8 TEV
    }
    // subsamples for stop scale tW channel
    else if(sample==kSAToptW1){
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =1492534;
      if(kSys==sysTopScaleDown) Nevents =1493101;
    }
    else if(sample==kSToptW1){
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =1492816;
      if(kSys==sysTopScaleDown) Nevents =1493130;
    }
    else if(sample==kSAToptW2){
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =497676;
      if(kSys==sysTopScaleDown) Nevents =497682;
    }
    else if(sample==kSToptW2){
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =455270;
      if(kSys==sysTopScaleDown) Nevents =496818;
    }
    else if(sample==kSAToptW3){
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =497376;
      if(kSys==sysTopScaleDown) Nevents =497674;
    }
    else if(sample==kSToptW3){
      // scale variation
      if(kSys==sysTopScaleUp  ) Nevents =442237;
      if(kSys==sysTopScaleDown) Nevents =453233;
    }
    // b) combined single top sample
    else if(sample==kSTop){
      // already added in combineMCsamples.C
      // with cross section as weight,
      // lumi normalization is done here
      Nevents=1;
      crossSection=1.;
    }
    //################
    else if(sample==kTTVjets){
      // already added in combineMCsamples.C
      // with cross section as weight,
      // lumi normalization is done here
      Nevents=1;
      crossSection=1.;
      if(excludettV) crossSection=0.;
   }
    //#################
    // Pythia6: DiBoson Samples
    // Summer12
    // a) subsamples
    else if(sample==kWW){
      crossSection= 54.8;
      Nevents     = 10000431;
    }
    else if(sample==kWZ){
      crossSection= 33.2;
      Nevents     = 10000283;
    }
    else if(sample==kZZ){
      crossSection= 8.1;
      Nevents     = 9799908;
    }
    // b) combined diboson: WW,WZ,ZZ
    else if(sample==kDiBos){
      // Summer12: already added in combineMCsamples.C
      // with cross section as weight,
      // lumi normalization is done here
      crossSection=1.;
      Nevents=1;
    }
    // Data
    else if(sample==kData){
      return 1.;
      // this leads to a weight
      // of 1.0 as data needs no weight
    }
    // unknown input
    else{
      std::cout << "ERROR: unknown configuration sample=" << sample;
      std::cout << ", decayChannel=" << decayChannel << std::endl;
      exit(0);
    }
    // d) calculate weight
    weight = effSFAB(kSys,decayChannel) * luminosity / ( Nevents / crossSection );
    if(verbose>0){
      std::cout << "sample: " << sampleLabel(sample,decayChannel) << std::endl;
      std::cout << "systematic var.: " << sysLabel(kSys) << std::endl;
    }
    // e) systematic shifts
    double weight2=weight;
    // e0) ttbar data basd reweighting
    if(kSys==sysTest){ // FIXME UPDATE FOR 8 TEV (update derived from ttbar signal)
      if(sample==kSig&&decayChannel.compare("electron")==0) weight*=1.02674242/1.02511720;//549.162/557.161;
      if(sample==kSig&&decayChannel.compare("muon"    )==0) weight*=1.02521244/1.02378495;//546.5/567.8;
      if(sample==kBkg&&decayChannel.compare("electron")==0) weight*=1.02674242/1.02511720;//3148.5/3193.74;
      if(sample==kBkg&&decayChannel.compare("muon"    )==0) weight*=1.02521244/1.02378495;//3151.23/3196.51;
    }
    if(verbose>1) std::cout << "weight before scaling: " << weight2 << std::endl;
    // e1) for ttbar/W BR correction
    if((sample==kSig)||(sample==kSigPow)||(sample==kSigMca)||(sample==kSigPowHer)) weight *= BRSF  (kSys);
    if((sample==kBkg)||(sample==kBkgPow)||(sample==kBkgMca)||(sample==kBkgPowHer)) weight *= BRSFBG(kSys);
    // e2) systematic higher/lower BG
    double scale=0;
    // (i) more/less DiBoson
    if(sample==kWW||sample==kWZ||sample==kZZ||sample==kDiBos){
      scale=0.3;
      if(kSys==sysDiBosUp  ) weight*=(1.0+scale);
      if(kSys==sysDiBosDown) weight*=(1.0-scale);
    }
    // (ii) more/less QCD
    if(sample==kQCD){
      scale=0.5;
      if(kSys==sysQCDUp  ) weight*=(1.0+scale);
      if(kSys==sysQCDDown) weight*=(1.0-scale);
    }
    // (iii) more/less single top
    if(sample==kSATops||sample==kSATopt||sample==kSAToptW||sample==kSTops||sample==kSTopt||sample==kSToptW||sample==kSTop){
      scale=0.3;
      if(kSys!=sysNo){
	if(kSys==sysSTopUp  ) weight*=(1.0+scale);
	if(kSys==sysSTopDown) weight*=(1.0-scale);
      }
    }
    // (iv) more/less v+jets
    if(sample==kWjets||sample==kZjets||sample==kTTVjets||sample==kTTGjets||sample==kTTZjets||sample==kTTWjets){ //#### added ttV samples
      scale=1.0;
      if(kSys!=sysNo){
	if(kSys==sysVjetsUp  ) weight*=(1.0+scale);
	if(kSys==sysVjetsDown) weight*=(1.0-scale);
	if(scale!=0&&verbose>0) std::cout << "possible scale factor: " << scale << std::endl;
      }
    }
    // printout for systematic variation
    if(verbose==1){
      std::cout << "weight";
      if(verbose>1){
	if(weight!=weight2) std::cout << "(scaled)";
	if(weight==weight2) std::cout << "(not scaled)";
      }
      std::cout << ": " << weight << std::endl;
      if(verbose>1) std::cout << "ratio: " << weight/weight2 << std::endl;
      if(weight!=weight2&&((sample==kSig)||(sample==kSigPow)||(sample==kSigPowHer)||(sample==kSigMca))) std::cout << "(BR correction applied)" << std::endl;
    }

    // return result
    if(sample!=kData&&weight==1){
      std::cout << "WARNING: function lumiweight";
      std::cout << " gives result 1 for sample: " << std::endl;
      std::cout << sampleLabel(sample,decayChannel) << ", " << decayChannel << " channel" << std::endl;
      std::cout << "used xSec:       " << crossSection << std::endl;
      std::cout << "used Nevents:    " << Nevents      << std::endl;
      std::cout << "used luminosity: " << luminosity   << std::endl;
    }
    return weight;
  }

  TString getStringEntry(const TString inputTString, unsigned int entry=42, const TString seperateBy="/")
  {
    // this function devides "inputTString" using "seperateBy"
    // as seperator and returns the "entry"-th element
    // 42: default, means last element after seperator
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    char *path, *element;
    // enable output for debugging
    bool verbose=false;
    // save inputstring in char* path
    std::string inputString=(std::string)inputTString;
    path = new char [inputString.size()+1];
    strcpy (path, inputString.c_str());
    //  path = (char*)inputTString.Data();
    // get pointer element to elements of path using seperateBy as seperator
    element=strtok (path,(char*)(seperateBy.Data()));
    // get vector with splittet elements
    std::vector<TString> result_;
    if(verbose){
      std::cout << "input: " << inputTString << std::endl;
      std::cout << "seperator: " << seperateBy << std::endl;
      std::cout << "splitted in the following elements:" << std::endl;
    }
    while(element!=NULL){
      if(verbose) std::cout << element << std::endl;
      std::string output="";
      output=output+element;
      result_.push_back((TString)output);
      element=strtok(NULL,(char*)(seperateBy.Data()));
    }
    delete[] path;
    delete[] element;
    // default configuration: get last entry
    if(entry==42) entry = (result_.size());
    // get requestet element
    if(verbose) std::cout << "requested " << "element #" << entry << " of " << result_.size() << std::endl;
    TString result="";
    // check if element exists
    if(entry<=result_.size()) result=result_[entry-1];
    else{
      std::cout << "warning - the requested element (" << entry <<") does not exist" << std::endl;
      std::cout << "the string "  << inputTString << " with seperator " << seperateBy << " gives only " <<  result_.size() << " elements!" << std::endl;
    }
    if(verbose) std::cout << "will return: " <<  result_[entry-1] << std::endl;
    // return output
    return result;
  }

  TString TopFilename(unsigned int sample, unsigned int sys, const std::string decayChannel, bool MS=MadSpin)
  {
    // this function contains the basic convention for the MC
    // .root files and returns the correct names for chosen samplesample"
    // and systematic variation "sys"
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: samples, systematicVariation

    TString fileName = "decayChannel_undefined" ;
    if      (decayChannel.compare("electron")==0) fileName ="elecDiffXSec";
    else if (decayChannel.compare("muon")==0)     fileName ="muonDiffXSec";
    else if (decayChannel.compare("combined")==0) fileName ="combinedDiffXSec"; 
    else fileName ="muonDiffXSec"; // default
    // name of data file is given directly in the .C file
    if(sample==kData) return "";
    // standard MC filenames
    if(sample==kSig    )fileName += "Sig";
    if(sample==kBkg    )fileName += "Bkg";
    // MadSpin: ttbar MadGraph sample, MS==true OR sysMad
    if(((MS&&sys!=sysMad)||(!MS&&sys==sysMad))&&((sample==kSig)||(sample==kBkg))&&sys<=sysHadDown) fileName+="MadSpin";
    if(sample==kSigPow )fileName += "SigPowheg";
    if(sample==kBkgPow )fileName += "BkgPowheg";
    if(sample==kSigPowHer )fileName += "SigPowhegHerwig";
    if(sample==kBkgPowHer )fileName += "BkgPowhegHerwig";
    if(sample==kSigMca )fileName += "SigMcatnlo";
    if(sample==kBkgMca )fileName += "BkgMcatnlo";
    if(sample==kWjets  )fileName += "Wjets";
    if(sample==kTTVjets  )fileName += "TTVjets"; //########## 
    if(sample==kZjets  )fileName += "Zjets";
    if(sample==kDiBos  )fileName += "VV";
    if(sample==kQCD    )fileName += "QCD";
    if(sample==kSTop   )fileName += "SingleTop";
    // subsamples are all located in subfolder called "MergedFiles"
    if(sample==kWW     )fileName = "MergedFiles/"+fileName+"WW";
    if(sample==kWZ     )fileName = "MergedFiles/"+fileName+"WZ";
    if(sample==kZZ     )fileName = "MergedFiles/"+fileName+"ZZ";
    if(sample==kQCDEM1 )fileName = "MergedFiles/"+fileName+"QCDEM1";
    if(sample==kQCDEM2 )fileName = "MergedFiles/"+fileName+"QCDEM2";
    if(sample==kQCDEM3 )fileName = "MergedFiles/"+fileName+"QCDEM3";
    if(sample==kQCDBCE1)fileName = "MergedFiles/"+fileName+"QCDBCE1"; 
    if(sample==kQCDBCE2)fileName = "MergedFiles/"+fileName+"QCDBCE2";
    if(sample==kQCDBCE3)fileName = "MergedFiles/"+fileName+"QCDBCE3"; 
    //####################
    if(sample==kTTGjets)fileName = "MergedFiles/"+fileName+"TTGjets";
    if(sample==kTTZjets)fileName = "MergedFiles/"+fileName+"TTZjets";
    if(sample==kTTWjets)fileName = "MergedFiles/"+fileName+"TTWjets";
    //################### 
    if(sample==kSToptW ||sample==kSToptW1 ||sample==kSToptW2 ||sample==kSToptW3 )fileName = "MergedFiles/"+fileName+"SingleTopTW";
    if(sample==kSTops  )fileName = "MergedFiles/"+fileName+"SingleTopS";
    if(sample==kSTopt  )fileName = "MergedFiles/"+fileName+"SingleTopT";
    if(sample==kSAToptW||sample==kSAToptW1||sample==kSAToptW2||sample==kSAToptW3)fileName = "MergedFiles/"+fileName+"SingleAntiTopTW";
    if(sample==kSATops )fileName = "MergedFiles/"+fileName+"SingleAntiTopS";
    if(sample==kSATopt )fileName = "MergedFiles/"+fileName+"SingleAntiTopT";
    // take care of systematic variations
    // they are located in dedicated subfolders
    // JES
    if(sys==sysJESUp  ) fileName = "JESUp/"+fileName+"JESUp"; 
    if(sys==sysJESDown) fileName = "JESDown/"+fileName+"JESDown";
    // JER
    if(sys==sysJERUp  ) fileName = "JERUp/"+fileName+"JERUp"; 
    if(sys==sysJERDown) fileName = "JERDown/"+fileName+"JERDown";
    // PDF Shape variation
    // only for new MC and ttbar signal
    if(sample==kSig){
      if(sys==sysPDFUp  ) fileName = "PDFUp/"+fileName+"PdfVarUp";
      if(sys==sysPDFDown) fileName = "PDFDown/"+fileName+"PdfVarDown";
    }
    // data based pt ttbar reweighting
    if( (sys==sysTest&&(sample==kSig||sample==kBkg))||(sys==sysTestPowheg&&(sample==kSigPow||sample==kBkgPow))||(sys==sysTestPowhegHerwig&&(sample==kSigPowHer||sample==kBkgPowHer))||(sys==sysTestMCatNLO&&(sample==kSigMca||sample==kBkgMca)) ) fileName = "ttbarReweight/"+fileName+"SysDistortdata";    
    // Q2-Scale
    // (a) top (ttbar+single top)
    if((sample==kSig)||(sample==kBkg)||(sample==kSTop)||(sample==kSToptW)||(sample==kSTops)||(sample==kSTopt)||(sample==kSAToptW)||(sample==kSATops)||(sample==kSATopt)||(sample==kSToptW1)||(sample==kSAToptW1)||(sample==kSToptW2)||(sample==kSAToptW2)||(sample==kSToptW3)||(sample==kSAToptW3)){
      if     (sys==sysTopScaleUp  ) fileName = "ScaleUp/"+fileName+"ScaleUp";    
      else if(sys==sysTopScaleDown) fileName = "ScaleDown/"+fileName+"ScaleDown";
      // additional prefix for tW scale subsamples
      if     (sample==kSToptW1||sample==kSAToptW1)fileName += "1";
      else if(sample==kSToptW2||sample==kSAToptW2)fileName += "2";
      else if(sample==kSToptW3||sample==kSAToptW3)fileName += "3";
    }
    // (b) V+jets
    if((sys==sysVBosonScaleUp  )&&(sample==kWjets)) fileName = "ScaleUp/"+fileName+"ScaleUp";    
    if((sys==sysVBosonScaleDown)&&(sample==kWjets)) fileName = "ScaleDown/"+fileName+"ScaleDown";
    if((sys==sysVBosonScaleUp  )&&(sample==kZjets)) {fileName = "ScaleUp/"+fileName+"ScaleUp";     }
    if((sys==sysVBosonScaleDown)&&(sample==kZjets)) {fileName = "ScaleDown/"+fileName+"ScaleDown"; }
    // (c) joint SingleTop
    if((sys==sysSingleTopScaleUp)  &&(sample==kSTop)) fileName = "ScaleUp/"+fileName+"ScaleUp";   
    if((sys==sysSingleTopScaleDown)&&(sample==kSTop)) fileName = "ScaleDown/"+fileName+"ScaleDown";
    // Matching Scale
    // (a) top   
    if((sys==sysTopMatchUp  )&&((sample==kSig)||(sample==kBkg))) fileName = "MatchUp/"+fileName+"MatchUp";   
    if((sys==sysTopMatchDown)&&((sample==kSig)||(sample==kBkg))) fileName = "MatchDown/"+fileName+"MatchDown";
    // (b) V+jets
    if((sys==sysVBosonMatchUp  )&&(sample==kWjets)) fileName = "MatchUp/"+fileName+"MatchUp";
    if((sys==sysVBosonMatchDown)&&(sample==kWjets)) fileName = "MatchDown/"+fileName+"MatchDown";
    if((sys==sysVBosonMatchUp  )&&(sample==kZjets)) {fileName = "MatchUp/"+fileName+"MatchUp";     }
    if((sys==sysVBosonMatchDown)&&(sample==kZjets)) {fileName = "MatchDown/"+fileName+"MatchDown"; }
    // Top Mass
    if((sysLabel(sys).Contains("sysTopMassUp"  ))&&((sample==kSig)||(sample==kBkg))){ 
      fileName = "TopMassUp/"+fileName+"TopMassUp";   
      if(     sys==sysTopMassUp2) fileName+="2";
      else if(sys==sysTopMassUp3) fileName+="3";
      else if(sys==sysTopMassUp4) fileName+="4"; 
    }
    else if((sysLabel(sys).Contains("sysTopMassDown"))&&((sample==kSig)||(sample==kBkg))){
      fileName = "TopMassDown/"+fileName+"TopMassDown";
      if(     sys==sysTopMassDown2) fileName+="2";
      else if(sys==sysTopMassDown3) fileName+="3";
      else if(sys==sysTopMassDown4) fileName+="4";      
    }
    // label for MC production cycle
    fileName += "Summer12";
    fileName += "PF.root";
    // debug
    if(sys==sysTestMCatNLO) std::cout << fileName << std::endl;
    // return output
    return fileName;
  }

  void saveCanvas(const std::vector<TCanvas*> MyCanvas, const TString outputFolder, const TString pdfName, const bool savePdf=true, const bool saveEps=true, const bool savePng=false )
  {
    // introduce function that saves every single canvas in
    // MyCanvas as ./outputFolder/CanvasTitle.eps and in addition
    // all in one rootFile: ./outputFolder/pdfName.pdf
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE

    // a) save all plots in one pdf
    if(savePdf){
      MyCanvas[0]->Print(outputFolder+pdfName+".pdf(", "pdf");
      for(unsigned int idx=1; idx<MyCanvas.size()-1; idx++){
	MyCanvas[idx]->Print(outputFolder+pdfName+".pdf", "pdf");
      }
      MyCanvas[MyCanvas.size()-1]->Print(outputFolder+pdfName+".pdf)", "pdf");
    }
    // b) save every plot as pdf with title as name
    if(saveEps){
      for(unsigned int idx=0; idx<MyCanvas.size(); idx++){
	MyCanvas[idx]->Print(outputFolder+(TString)(MyCanvas[idx]->GetTitle())+".eps");
      }
    }
    // c) save every plot as png with title as name
    if(savePng){
      for(unsigned int idx=0; idx<MyCanvas.size(); idx++){
	MyCanvas[idx]->Print(outputFolder+(TString)(MyCanvas[idx]->GetTitle())+".png");
      }
    }
  }

  std::map<unsigned int, TFile*> getStdTopAnalysisFiles( const TString inputFolder, const unsigned int systematicVariation, const TString dataFile, const std::string decayChannel, TString ttbarMC="Madgraph", bool MS=MadSpin)
    {
      // this function returns a map containing all existing .root in "inputFolder"
      // corresponding to the chosen systematic variation "systematicVariation"
      // and the data file "dataFile" for which the direct path/name.root is needed
      // The MC .root file names correspond to the standard names as defined in the
      // function TopFilename
      // modified quantities: NONE
      // used functions: TopFilename
      // used enumerators: samples
      // ttbarMC: use "Madgraph", "Mcatnlo", "Powheg" or "PowhegHerwig" for the respective samples

      bool debug=false;
      // open our standard analysis files and save them in a map
      std::map<unsigned int, TFile*> files_;
      // loop samples
      for(int sample = kSig; sample<=kSAToptW; sample++){
	// there are no QCD subsamples for muon channel
	if(!(sample>=kQCDEM1&&sample<=kQCDBCE3&&decayChannel.compare("muon")==0)){
	  TString fileName =""; 
	  TString fileName2=""; 
	  TString fileName3="";
	  int sampleflag=sample;
	  int sampleflag2=-1;
	  int sampleflag3=-1;
	  if(sample!=kData){
	    if(((sample==kSig)||(sample==kBkg))&&ttbarMC=="Powheg") {
	      if(sample==kSig) fileName = inputFolder+"/"+TopFilename(kSigPow, systematicVariation, decayChannel);
	      if(sample==kBkg) fileName = inputFolder+"/"+TopFilename(kBkgPow, systematicVariation, decayChannel);
	    }
	    else if(((sample==kSig)||(sample==kBkg))&&ttbarMC=="PowhegHerwig") {
	      if(sample==kSig) fileName = inputFolder+"/"+TopFilename(kSigPowHer, systematicVariation, decayChannel);
	      if(sample==kBkg) fileName = inputFolder+"/"+TopFilename(kBkgPowHer, systematicVariation, decayChannel);
	    }
	    else if(((sample==kSig)||(sample==kBkg))&&ttbarMC=="Mcatnlo") {
	      if(sample==kSig) fileName = inputFolder+"/"+TopFilename(kSigMca, systematicVariation, decayChannel);
	      if(sample==kBkg) fileName = inputFolder+"/"+TopFilename(kBkgMca, systematicVariation, decayChannel);
	    }
	    else{
	      fileName = inputFolder+"/"+TopFilename(sample, systematicVariation, decayChannel, MS);
	    }
	    // take care of splitted single top tW sample for scale systematics
	    if((systematicVariation==sysSingleTopScaleUp||systematicVariation==sysSingleTopScaleDown)&&(sample==kSToptW||sample==kSAToptW)){
	      if(sample==kSToptW){
		sampleflag =kSToptW1;
		sampleflag2=kSToptW2;
		sampleflag3=kSToptW3;
	      }
	      else if(sample==kSAToptW){
		sampleflag =kSAToptW1;
		sampleflag2=kSAToptW2;
		sampleflag3=kSAToptW3;
	      }
	      fileName =inputFolder+"/"+TopFilename(sampleflag , systematicVariation, decayChannel);
	      fileName2=inputFolder+"/"+TopFilename(sampleflag2, systematicVariation, decayChannel);
	      fileName3=inputFolder+"/"+TopFilename(sampleflag3, systematicVariation, decayChannel);	      
	    }
	  }
	  else fileName = dataFile;
	  // if file exists - save in map
	  if((fileName!="no")&&(fileName!="")){
	    TFile* file = TFile::Open(fileName);
	    if(file&&!(file->IsZombie())){
	      files_[sampleflag]= file;
	      if(debug) std::cout << "found " << fileName << std::endl;
	    }
	  }
	  if(sampleflag2>-1&&(fileName2!="no")&&(fileName2!="")){
	    TFile* file = TFile::Open(fileName);
	    if(file&&!(file->IsZombie())) files_[sampleflag2]= file;
	  }
	  if(sampleflag3>-1&&(fileName3!="no")&&(fileName3!="")){
	    TFile* file = TFile::Open(fileName);
	    if(file&&!(file->IsZombie())) files_[sampleflag3]= file;
	  }
	}
      }
      return files_;
    }

  void getAllPlots( std::map<unsigned int, TFile*> files_, const std::vector<TString> plotList_,  std::map< TString, std::map <unsigned int, TH1F*> >& histo_, std::map< TString, std::map <unsigned int, TH2F*> >& histo2_, const unsigned int N1Dplots, int& Nplots, const int verbose=0, const std::string decayChannel = "unset", std::vector<TString> *vecRedundantPartOfNameInData = 0, bool SSV=false, TString ignorePartNameInMC="", std::vector<TString> *vecRedundantPartOfNameInNonTTbarSG = 0)
  {
    // this function searches for every plot listed in "plotList_" in all files listed in "files_",
    // saves all 1D histos into "histo_" and all 2D histos into "histo2_"
    // empty plots will be neglected
    // modified quantities: "histo_", "histo2_", "Nplots"
    // used functions: sampleLabel
    // used enumerators: samples
    // "N1Dplots": the #1D plots is needed as input to distinguish between 1D and 2D plots
    // "Nplots": the total # of plots is calculated
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "redundantPartOfNameInData": delete this part in the (folder)name when loading the plots from data
    //                              (needed to handle systematic variations where foldername in data and MC is different)
    // "SSV": for all btagging plots SSV control plots are used instead of the default (CSV) plots
    // "ignorePartNameInMC": like redundantPartOfNameInData, but for MC instead of data
    // "ignorePartNameInNonTTbarSG": like ignorePartNameInMC, but for non ttbar MC only, expects the format plotname/ignoreString

    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      TString testPlotNameForTrackingDownErrors="NotUsedAtTheMoment";//"analyzeTopRecoKinematicsKinFitProbSel/topPt";//"NotUsedAtTheMoment";
      //TString testPlotNameForTrackingDownErrors="analyzeTopRecoKinematicsKinFitPUupProbSel/topPt";//"NotUsedAtTheMoment";
      bool hugo=plotList_[plot].Contains(testPlotNameForTrackingDownErrors) ? true : false;
      // check if plot exists in any sample
      bool existsInAnySample=false;
      if(hugo){
	std::cout << "name of target plot: " << plotList_[plot] << std::endl;
	std::cout << "ignore in name for MC: " << ignorePartNameInMC << std::endl;
	if (vecRedundantPartOfNameInData != 0 && vecRedundantPartOfNameInData->size() != 0){
	  std::cout << "ignore in name for data: ";
	  std::vector<TString>::iterator iter;
	  for (iter = (vecRedundantPartOfNameInData->begin()); iter != (vecRedundantPartOfNameInData->end()); iter++){
	     std::cout << *iter;
	     if(iter != (vecRedundantPartOfNameInData->end()-1))  std::cout << ", ";
	     else  std::cout << std::endl;
	  }
	}
      }
      // loop samples
      for(unsigned int sample=kSig; sample<=kSAToptW; ++sample){
	// check existence of sample
	if(hugo) std::cout << "-> file " << files_[sample]->GetName() << std::endl;
	if(files_.count(sample)!=0){
	  if(hugo) std::cout << "-> found!" << std::endl;
	  // create plot container
	  TH1* targetPlot;
	  // delete additional part of MC foldername
	  // that does not exist in data 
	  TString plotname=plotList_[plot];	  
	  if(hugo) std::cout << "plot: " << plotList_[plot] << std::endl;
	  if(sample==kData){
	    if (vecRedundantPartOfNameInData != 0 && vecRedundantPartOfNameInData->size() != 0){
	      std::vector<TString>::iterator iter;
	      for (iter = (vecRedundantPartOfNameInData->begin()); iter != (vecRedundantPartOfNameInData->end()); iter++){
		plotname.ReplaceAll((*iter), ""); 
	      }
	    }
	    // by hand fix for probability label in mixed object analyzer data folder
	    if(ignorePartNameInMC=="ProbSel"&&plotname.Contains("compositedKinematicsKinFit")){
	      plotname.ReplaceAll("compositedKinematicsKinFit","compositedKinematicsProbSel");
	    }
	  }
	  else{
	    plotname.ReplaceAll(ignorePartNameInMC, ""); 
	    if(sample!=kSig){
	      if (vecRedundantPartOfNameInNonTTbarSG != 0 && vecRedundantPartOfNameInNonTTbarSG->size() != 0){
		std::vector<TString>::iterator iter;
		for (iter = (vecRedundantPartOfNameInNonTTbarSG->begin()); iter != (vecRedundantPartOfNameInNonTTbarSG->end()); iter++){
		  TString targetPlotString = getStringEntry((*iter), 1);
		  TString ignorePartNameInNonttbar= getStringEntry((*iter), 2);
		  std::cout << ignorePartNameInNonttbar << std::endl;
		  if(plotname.Contains(targetPlotString)){
		    plotname.ReplaceAll(ignorePartNameInNonttbar, "");
		  } // end if plot is target plot
		} // end loop vecRedundantPartOfNameInNonTTbarSG elements
	      } // end if vecRedundantPartOfNameInNonTTbarSG non empty
	    } // end if !SG
	  } // end else kData
	  if(hugo) std::cout << "-> searching plot " << plotname << std::endl;
	  files_[sample]->GetObject(plotname, targetPlot);
	  // Check existence of plot
	  if(targetPlot){ 
	    existsInAnySample=true;
	    // go to end of loop
	    sample=kSAToptW;
	    if(hugo) std::cout << "-> found! " << std::endl;
	  }
	  else if(hugo){std::cout << "-> NOT found! " << std::endl;}
	}
      }
      // end program and draw error if plot does not exist at all
      if(!existsInAnySample){
	std::cout << "ERROR: plot "+plotList_[plot]+" does not exist in any file!" << std::endl;
	exit(1);
      }
      else{
	// otherwise: get plots from sample
	// loop samples
	for(unsigned int sample=kSig; sample<=kSAToptW; ++sample){
	  // delete additional part of MC foldername
	  // that does not exist in data 
	  TString plotname=plotList_[plot];
	  // replace CSV with SSV plots
	  if(SSV&&!plotname.Contains("bottomJet")) plotname.ReplaceAll("Tagged", "SSV");
	  if(sample==kData){
	    if (vecRedundantPartOfNameInData != 0 && vecRedundantPartOfNameInData->size() != 0){
	      std::vector<TString>::iterator iter;
	      for (iter = (vecRedundantPartOfNameInData->begin()); iter != (vecRedundantPartOfNameInData->end()); iter++){
		if (plotname.Contains((*iter))) plotname.ReplaceAll((*iter), "");
	      }
	    }
	    // by hand fix for probability label in mixed object analyzer data folder
	    if(ignorePartNameInMC=="ProbSel"&&plotname.Contains("compositedKinematicsKinFit")){
	      plotname.ReplaceAll("compositedKinematicsKinFit","compositedKinematicsProbSel");
	    }
	  }
	  else {
	    plotname.ReplaceAll(ignorePartNameInMC, ""); 
	    if(sample!=kSig){
	      if (vecRedundantPartOfNameInNonTTbarSG != 0 && vecRedundantPartOfNameInNonTTbarSG->size() != 0){
		std::vector<TString>::iterator iter;
		for (iter = (vecRedundantPartOfNameInNonTTbarSG->begin()); iter != (vecRedundantPartOfNameInNonTTbarSG->end()); iter++){
		  TString targetPlotString = getStringEntry((*iter), 1);
		  TString ignorePartNameInNonttbar= getStringEntry((*iter), 2);
		  if(plotname.Contains(targetPlotString)){
		    // debug
		    if(verbose>1) std::cout << "INFO: plot " << plotname << " in sample " << sampleLabel(sample, decayChannel) << " contains " << targetPlotString << " - will ignore " << ignorePartNameInNonttbar << " in name" << std::endl;
		    plotname.ReplaceAll(ignorePartNameInNonttbar, "");
		  } // end if plot is target plot
		} // end loop vecRedundantPartOfNameInNonTTbarSG elements
	      } // end if vecRedundantPartOfNameInNonTTbarSG non empty
	    } // end if !SG
	  } // end else kData
	  // check if file exists
	  // give warning if file does not exist
	  if((files_.count(sample)==0)&&(plot==0)&&(verbose>0)) std::cout << "file for " << sampleLabel(sample,decayChannel) << " does not exist- continue and neglect this sample" << std::endl;
	  if(files_.count(sample)>0){
	    // create plot container
	    TH1* targetPlot;
	    if(verbose>1){
	      std::cout << std::endl << "sample: " << sample << ", " << files_[sample]->GetName() << std::endl;
	      std::cout << "plot: " << plot << ", " << plotname << std::endl;
	    }
	    files_[sample]->GetObject(plotname, targetPlot);
	    if(hugo) std::cout << files_[sample]->GetName() << " -> " << plotname << std::endl;
	    // Check if plot exits
	    // give warning if plot does not exist
	    if((!targetPlot)&&(verbose>0)) std::cout << "can not find plot "+plotname << " in file "+(TString)(files_[sample]->GetName()) << " - continue and neglect this plot" << std::endl;
	    if(targetPlot){
	      // check if plot is empty
	      bool emptyPlot=false;
	      if((plot<N1Dplots )&&((((TH1*)(files_[sample]->Get(plotname)))->GetEntries())==0.)) emptyPlot=true;
	      if((plot>=N1Dplots)&&((((TH2*)(files_[sample]->Get(plotname)))->GetEntries())==0.)) emptyPlot=true;
	      if(emptyPlot && verbose>1) std::cout << "plot "+plotList_[plot] << " in file "+(TString)(files_[sample]->GetName()) << " is empty- continue and neglect this plot" << std::endl;
	      // to avoid problems with samples where no event is passing the selection we will drop this requirement by now
	      emptyPlot=false;
	      if(!emptyPlot){
		// save plot in corresponding map
		if(plot<N1Dplots ) histo_ [plotList_[plot]][sample] = (TH1F*)(files_[sample]->Get(plotname)->Clone());
		if(plot>=N1Dplots) histo2_[plotList_[plot]][sample] = (TH2F*)(files_[sample]->Get(plotname)->Clone());
		// count every existing 2D plot (every sample is counted separetly as it will be drawn into an own canvas)
		if(plot>=N1Dplots) Nplots++;
	      }
	    }
	  }
	}
	// count every existing type of 1D plots (neglect #samples here as they well be drawn into the same canvas)
	if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)) Nplots++;
      }
    }
  }

  void scaleByLuminosity(const std::vector<TString> plotList_,  std::map< TString, std::map <unsigned int, TH1F*> >& histo_, std::map< TString, std::map <unsigned int, TH2F*> >& histo2_, const unsigned int N1Dplots, const double luminosity, const int verbose=1, const int systematicVariation=sysNo, const std::string decayChannel = "unset", TString ttbarMC="Madgraph")
  {
    // this function scales every histo in histo_ and histo2_ to the corresponding luminosity
    // Additionally the mu eff SF is applied
    // modified quantities: "histo_", "histo2_"
    // used functions: sampleLabel, effSFAB
    // used enumerators: samples
    // "N1Dplots": the #1D plots is needed as input to destinguish between 1D and 2D plots
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "luminosity": [/pb]
    // "systematicVariation": specify systematic variation corresponding to enum systematicVariation
    // ttbarMC: use "Madgraph", "Mcatnlo", "Powheg" or "PowhegHerwig" for the respective samples

    // loop samples
    for(unsigned int sample=kSig; sample<=kSAToptW; ++sample) {
      unsigned int sample2=sample;
      if(ttbarMC=="Powheg"){
	if(sample==kSig) sample2=kSigPow;
	if(sample==kBkg) sample2=kBkgPow;
      }
      else if(ttbarMC=="PowhegHerwig"){
	if(sample==kSig) sample2=kSigPowHer;
	if(sample==kBkg) sample2=kBkgPowHer;
      }
      else if(ttbarMC=="Mcatnlo"){
	if(sample==kSig) sample2=kSigMca;
	if(sample==kBkg) sample2=kBkgMca;
      }
      // loop plots
      for(unsigned int plot=0; plot<plotList_.size(); ++plot){
	// a) 1D
	// check if 1D plot exists
	if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	  if(verbose>1) std::cout << std::endl << "plot: "+plotList_[plot] << " for sample " << sampleLabel(sample2,decayChannel) << ":" << std::endl;
	  if(verbose>1) std::cout << "#events before weighting: " << histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()+1) << std::endl;
	  // scale MC samples to same luminosity
	  double weight = lumiweight(sample2, luminosity, systematicVariation, decayChannel);
	  histo_[plotList_[plot]][sample]->Scale(weight);
	  if(verbose>1) std::cout << "weight: " << weight << std::endl;
	  if(verbose>1) std::cout << "#events after weighting: " << histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()+1) << std::endl;
	  // apply eff. SF for MC
	  //if(sample!=kData) histo_[plotList_[plot]][sample]->Scale(effSFAB(systematicVariation,decayChannel));
	}
	// b) 2D
	// check if 2D plot exists
	if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	  if(verbose>1) std::cout << std::endl << "plot: "+plotList_[plot] << " for sample " << sampleLabel(sample2,decayChannel) << ":" << std::endl;
	  if(verbose>1) std::cout << "#events before weighting: " << histo2_[plotList_[plot]][sample]->Integral(0, histo2_[plotList_[plot]][sample]->GetNbinsX()+1, 0, histo2_[plotList_[plot]][sample]->GetNbinsY()+1) << std::endl;
	  // scale MC samples to same luminosity
	  double weight = lumiweight(sample2, luminosity, systematicVariation, decayChannel);
	  histo2_[plotList_[plot]][sample]->Scale(weight);
	  if(verbose>1) std::cout << "weight: " << weight << std::endl;
	  if(verbose>1) std::cout << "#events after weighting: " << histo2_[plotList_[plot]][sample]->Integral(0, histo2_[plotList_[plot]][sample]->GetNbinsX()+1, 0, histo2_[plotList_[plot]][sample]->GetNbinsY()+1) << std::endl;
	  // apply eff. SF for MC
	  //if(sample!=kData) histo2_[plotList_[plot]][sample]->Scale(effSFAB(systematicVariation,decayChannel));
	}
      }
    }
  }

  void AddSingleTopAndDiBoson(const std::vector<TString> plotList_,  std::map< TString, std::map <unsigned int, TH1F*> >& histo_, std::map< TString, std::map <unsigned int, TH2F*> >& histo2_, const unsigned int N1Dplots, const int verbose=1, bool reCreate=false, const std::string decayChannel = "unset")
  {
    // this function creates plots for all diboson and all single
    // top samples combined if the combined SingleTop and DiBoson
    // samples do not exist as .root file.
    // Every plot in "plotList_" existing in each of the three
    // sTop/ diBoson files will be combined and saved in the histo_
    // and histo2_ maps
    // NOTE: all plots from the samples are considered to be weighted
    // to the same luminosity before
    // modified quantities: "histo_", "histo2_"
    // used functions: sampleLabel
    // used enumerators: samples
    // "N1Dplots": the #1D plots is needed as input to distinguish between 1D and 2D plots
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "reCreate": choose if you want to create the combined plot from the
    // single plots if it alredy exists. Careful: the old one will be deleted

    if(verbose>0) std::cout << std::endl;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // loop STop and DiBoson
      for(unsigned int sample=kSTop; sample<=kDiBos; ++sample){
	if(verbose>1) std::cout << "plot " << plotList_[plot] << ": " << std::endl;
	// mark first plot found
	bool first=true;
	// check that combined plot does not already exist
	// if existing and reCreate==1, it will be overwritten
	bool combinedplotExists=false;
	if((plot<N1Dplots )&&(histo_ [plotList_[plot]].count(sample)>0)) combinedplotExists=true;
	if((plot>=N1Dplots)&&(histo2_[plotList_[plot]].count(sample)>0)) combinedplotExists=true;
	if(reCreate==1) combinedplotExists=false;
	if((combinedplotExists==true)&&(verbose>0)){
	  std::cout << "combined plot " << plotList_[plot];
	  std::cout << " for sample " << sampleLabel(sample,decayChannel);
	  std::cout << " already exists, will keep this one" << std::endl;
	}
	if(!combinedplotExists){
	  // loop all three subchannels
	  unsigned int firstSubChannel=kSTops;
	  if(sample==kDiBos) firstSubChannel=kWW;
	  for(unsigned int subSample=firstSubChannel; subSample<=firstSubChannel+5; ++subSample){
	    // single top has 6, diboson only 2 subsamples
	    if((sample==kDiBos&&subSample<firstSubChannel+3)||sample==kSTop){
	      // check if plot exists for the subchannel
	      bool subPlotExists=false;
	      if((plot<N1Dplots )&&(histo_ [plotList_[plot]].count(subSample)>0)) subPlotExists=true;
	      if((plot>=N1Dplots)&&(histo2_[plotList_[plot]].count(subSample)>0)) subPlotExists=true;
	      if((subPlotExists==false)&&(verbose>0)){
		std::cout << "plot " << plotList_[plot];
		std::cout << " does not exist for subSample ";
		std::cout << sampleLabel(subSample,decayChannel) << std::endl;
	      }
	      if(subPlotExists){
		// add histo
		// a) 1D
		if(plot<N1Dplots){
		  if(first ) histo_[plotList_[plot]][sample]   =  (TH1F*)histo_[plotList_[plot]][subSample]->Clone();
		  if(!first) histo_[plotList_[plot]][sample]->Add((TH1F*)histo_[plotList_[plot]][subSample]->Clone());
		}
		// b) 2D
		if(plot>=N1Dplots){
		  if(first ) histo2_[plotList_[plot]][sample]   =  (TH2F*)histo2_[plotList_[plot]][subSample]->Clone();
		  if(!first) histo2_[plotList_[plot]][sample]->Add((TH2F*)histo2_[plotList_[plot]][subSample]->Clone());
		}
		// indicate that already one plot is found
		first=0;
		// print out information
		if(verbose>0){
		  std::cout << "will add " << plotList_[plot];
		  std::cout << " from " << sampleLabel(subSample,decayChannel) << std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void createStackPlot(const std::vector<TString> plotList_, std::map< TString, std::map <unsigned int, TH1F*> >& histo_, const unsigned int plot, const unsigned int N1Dplots, const int verbose=1, const std::string decayChannel="muon")
  {
    // this function will transform the histogram "plotList_"["plot"]
    // within "histo_" into stacked plots.
    // modified quantities: "histo_"
    // used functions: sampleLabel
    // used enumerators: samples
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    // "N1Dplots": the #1D plots is needed as input to destinguish between 1D and 2D plots

    // loop samples backwards
    for(int sample=kDiBos-1; sample>=kSig; --sample){
      // check if plot exists
      if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	// check if previous plot also exists
	if(histo_[plotList_[plot]].count(sample+1)>0){
	  histo_[plotList_[plot]][sample]->Add((TH1F*)histo_[plotList_[plot]][sample+1]->Clone());
	  if(verbose>1) std::cout << "add "+plotList_[plot] << " plot for " << sampleLabel(sample,decayChannel) << " and " << sampleLabel(sample+1,decayChannel) << std::endl;
	  if(verbose>1) std::cout << "this new stacked plot contains now " <<  histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()+1) << " events" << std::endl;
	}
      }
    }
  }

  float modulo(const float a, const float b)
  {
    // this function calculates a modulo b
    // where a and b are float
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    // "a": dividend
    // "b": divisor

    // prevent errors due to rounding effects
    // a) get the last relevant digit "dig" 
    int dig=getRelevantDigits(a);
    // b) convert double to int value by multiplying 10^"dig"
    if(dig<getRelevantDigits(b)) dig=getRelevantDigits(b);
    float atemp= a;
    float btemp= b;
    double shift=1.0;
    for(int n=0; n < dig; ++n) {
      atemp*=10.;
      btemp*=10.;
      shift*=10.;
    }
    // c) use int rounding to get rid of inprecise values
    int astar= roundToInt(atemp);
    int bstar= roundToInt(btemp);
    // control printout
    //std::cout << std::endl << a << " modulo ";
    //std::cout << b << std::endl;
    //std::cout << "-> " << astar << " modulo " << bstar << " / " << shift << std::endl;
    float rest=a-b;
    if(rest<0){
      std::cout << "can not compute " << a << " modulo " << b << std::endl;
      std::cout << "= " << astar << " modulo " << bstar << " / " << shift << std::endl;
      std::cout << "because " << a << " < " << b << std::endl;
      exit(1);
    }
    // use int modulo function
    rest=float(astar%bstar)/shift;
    if(astar%bstar==0) rest=0.;
    return rest;
  }

  TH1F* cutOffTH1(TH1F* histoOri, const double cutOff, const int verbose=0){
    // this function removes all bins with x>cutOff and moves their content to the overflow
    // modified quantities: "histoOri"
    // ATTENTION: histoUnbinned needs to have an initial equidistant binning
    // used functions: none
    // used enumerators: none
    // "histoUnbinned": plot to be binned
    // "cutOff": new maximum of the x Axis
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
    
    // new #bins
    double Nbins=-1;
    // find last remaining bin
    for(int bin=1; bin<=histoOri->GetNbinsX(); ++bin){
      if(histoOri->GetBinLowEdge(bin+1)>cutOff){
	Nbins=bin-1;
	break;
      }
    }
    if(Nbins>0&&Nbins<histoOri->GetNbinsX()){
      // create new histo
      TH1F* result=new TH1F(TString(histoOri->GetName())+"cutoff", TString(histoOri->GetTitle())+"cutoff", Nbins, histoOri->GetBinLowEdge(1), histoOri->GetBinLowEdge(Nbins+1));
      // fill new histo
      for(int bin=1; bin<=histoOri->GetNbinsX(); ++bin){
	if(bin<=Nbins) result->SetBinContent(bin, histoOri->GetBinContent(bin));
	else  result->SetBinContent(result->GetNbinsX()+1, result->GetBinContent(result->GetNbinsX()+1)+histoOri->GetBinContent(bin));
      }
      // histo with new range
      return result;
    }
    else if(verbose>0){
      std::cout << "WARNING in cutOffTH1: histo " << histoOri->GetName() << " not changed (target cutOff was " << cutOff << ")" << std::endl;
    }
    // histo stays untouched
    return histoOri;
  }


  void reBinTH1F(TH1F& histoUnbinned, const std::vector<double> &binlowerEdges_, const int verbose=0)
  {
    // this function rebins an histogram using a variable binning
    // modified quantities: "histoUnbinned"
    // ATTENTION: histoUnbinned needs to have an initial equidistant binning
    // used functions: modulo
    // used enumerators: none
    // "histoUnbinned": plot to be binned
    // "binlowerEdges": vector containing all lower bin edges starting at xaxis.min, ending with xaxis.max
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )

    unsigned int NinitialBins=histoUnbinned.GetNbinsX();
    double xMin=histoUnbinned.GetBinLowEdge(1);
    double xMax=histoUnbinned.GetBinLowEdge(NinitialBins+1);
    if(verbose>1){
      std::cout << "plot:            " << histoUnbinned.GetName() << std::endl;
      std::cout << "min:             " << xMin << std::endl;
      std::cout << "max:             " << xMax << std::endl;
      std::cout << "N(initial bins): " << NinitialBins << std::endl;
    }
    // check if chosen binning is valid
    // 1) N(bins)
    if(NinitialBins<=binlowerEdges_.size()){
      std::cout << "histo " << histoUnbinned.GetName() << " cannot be rebinned" << std::endl;
      std::cout << "N(chosen bins) < N(initial bins)!" << std::endl;
      exit(1);
    }
    // 2) coarseness of chosen binning
    //  double initialBinWidth=(xMax-xMin)/(double)NinitialBins;
    double initialBinWidth = round(histoUnbinned.GetBinWidth(1), getRelevantDigits(histoUnbinned.GetBinWidth(1)));
    //double initialBinWidth = (double)(roundToInt(10000.*histoUnbinned.GetBinWidth(1)))/10000.;
    
    if(verbose>1) std::cout << "initial binwidth: " << initialBinWidth << std::endl;
    if(binlowerEdges_.size()>1){
      for(unsigned int finalBin=1; finalBin<binlowerEdges_.size()-1; ++finalBin){
	double finalBinWidth=binlowerEdges_[finalBin]-binlowerEdges_[finalBin-1];
	//finalBinWidth=round(finalBinWidth, getRelevantDigits(finalBinWidth));
	if(verbose>1){
	  std::cout << "bin #" << finalBin << ": ";
	  std::cout << std::setprecision(20) << std::fixed << finalBinWidth << " modulo ";
	  std::cout << std::setprecision(20) << std::fixed << initialBinWidth;
	  std::cout << " = " << modulo(finalBinWidth, initialBinWidth) << std::endl;
	}
	if(modulo(finalBinWidth, initialBinWidth)!=0&&verbose>-1){
	  std::cout << "WARNING: histo " << histoUnbinned.GetName() << " may not be rebinned correctly!" << std::endl;
	  std::cout << "the ininital binning is to coarse for the chosen binning!" << std::endl;
	  //std::cout << "attention: probably error in modulo function," << std::endl;
	  std::cout << "bin #" << finalBin << " of the chosen binning has ";
	  std::cout << "finalBinWidth modulo initialBinWidth of:" <<std::endl;
	  std::cout << std::setprecision(20) << std::fixed << finalBinWidth;
	  std::cout << " modulo " << std::setprecision(20) << std::fixed << initialBinWidth << " = ";
	  std::cout << std::setprecision(20) << std::fixed << modulo(finalBinWidth,initialBinWidth) << std::endl;
	  //exit(1);
	}
      }
    }
    // fill vector entries in array
    // because TH1F constructor needs an array
    const unsigned int arraySize=binlowerEdges_.size();
    double *binLowerEdgeArray = new double[arraySize];
 
    for(unsigned int arrayEntry=0; arrayEntry<arraySize; ++arrayEntry){
      binLowerEdgeArray[arrayEntry]=binlowerEdges_[arrayEntry];
      if(verbose>1) std::cout << "array entry #" << arrayEntry << ": " << binLowerEdgeArray[arrayEntry] << " ";
    }
    if(verbose>1) std::cout << std::endl;
   
    if(verbose>1) std::cout << "N(bins) before rebinning: " << histoUnbinned.GetNbinsX() << std::endl;
    histoUnbinned = *(TH1F*)histoUnbinned.Rebin(arraySize-1, histoUnbinned.GetName(), binLowerEdgeArray);
    if(verbose>1){
      std::cout << "N(bins) after rebinning: " << histoUnbinned.GetNbinsX() << " (should be " << arraySize-1 << ")" << std::endl;
      std::cout << "from which only " << binlowerEdges_.size() << " bins are used" << std::endl;
      std::cout << "the rest will be ignored and exists only for technical reasons" << std::endl;
    }
    delete binLowerEdgeArray;
  }


  TH2F* reBinTH2F(TH2F& histoUnbinned, const std::vector<double> &binlowerEdges_, const int verbose=0)
  {
    // this function rebins an 2D histogram using a variable binning
    // modified quantities: "histoUnbinned"
    // used functions: none
    // used enumerators: none
    // "histoUnbinned": plot to be binned
    // "binlowerEdges": vector containing all lower bin edges starting at xaxis.min, ending with xaxis.max
    // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )

    // fill vector entries in array
    // because TH1F constructor needs an array
    const unsigned int arraySize=binlowerEdges_.size();
    double *binLowerEdgeArray = new double[arraySize];
    for(unsigned int arrayEntry=0; arrayEntry<arraySize; ++arrayEntry){
      binLowerEdgeArray[arrayEntry]=binlowerEdges_[arrayEntry];
      if(verbose>1) std::cout << "array entry #" << arrayEntry << ": " << binLowerEdgeArray[arrayEntry] << " ";
    }
    if(verbose>1) std::cout << std::endl;
    
    // create 2D histo with new binning
    if(verbose>1) std::cout << "create binned histo" << std::endl;
    TH2F* binnedHisto= new TH2F(TString(histoUnbinned.GetName())+"Binned", histoUnbinned.GetTitle(), arraySize-1, binLowerEdgeArray, arraySize-1, binLowerEdgeArray);
    if(verbose>1) std::cout << "get old axis" << std::endl;
    TAxis *xaxis = histoUnbinned.GetXaxis();
    TAxis *yaxis = histoUnbinned.GetYaxis();
    if(verbose>1) std::cout << "refill values" << std::endl;
    for(int j=1;j<=yaxis->GetNbins();j++) {
      for (int i=1;i<=xaxis->GetNbins();i++) {
	binnedHisto->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),histoUnbinned.GetBinContent(i,j));
      }
    }
    delete binLowerEdgeArray;
    return binnedHisto;
  }


    // ===============================================================
    // ===============================================================

    template <class histoType, class varType> 
      histoType* reBinTH1FIrregularNewBinning(histoType *histoOldBinning, const std::vector<varType> &vecBinning, TString plotname, bool rescale=1, int verbose=-1)
    {
	//  This function rebins a histogram using a variable binning
	// 
	//  (1) It is not required to have an equidistant binning.
	//  (2) Any type of ROOT-histgramme can be used, so the template 
	//      arguments should be 
	//      (a) histoT = TH1D,   TH1F,  ....
	//      (b) varT   = double, float, ....
	//  
	//  modified quantities: none
	//  used functions:      none
	//  used enumerators:    none
	//  
	//  "histoOldBinning":   plot to be re-binned
	//  "vecBinning":        vector containing all bin edges 
	//                       from xaxis.min to xaxis.max
	//  "rescale":           rescale the rebinned histogramme
	//                       (applicable for cross-section e.g.) 
        //  "verbose":           level of output 
      
        unsigned int vecIndex=0;
	
	// fill vector into array to use appropriate constructor of TH1-classes
	const varType *binArray = vecBinning.data();
	
	// create histo with new binning
	TString name = (TString)histoOldBinning->GetName();

	// suppress error massage
	// FIXME: source of error should be revised at some point
        int initialIgnoreLevel=gErrorIgnoreLevel;
	if(verbose==-1) gErrorIgnoreLevel=kError;
	histoType histoNewBinning = histoType ("histoNewBinning"+plotname+name,"histoNewBinning"+plotname+name,vecBinning.size()-1,binArray);
	gErrorIgnoreLevel=initialIgnoreLevel;

	// fill contents of histoOldBinning into histoNewBinning and rescale if applicable
	for (vecIndex = 0; vecIndex < vecBinning.size()-1; vecIndex++){
	    
	    varType lowEdge      = vecBinning[vecIndex]; 
	    //if (plotname=="topPt"&&vecIndex==0&&lowEdge==0.0) lowEdge+=10;  // adhoc fix to compensate for minimum top-Pt cut in NNLO curve
	    varType highEdge     = vecBinning[vecIndex+1];
	    varType newBinWidth  = highEdge - lowEdge;
	    varType newBinCenter = 0.5*(highEdge+lowEdge);
	    varType binSum       = 0.0;	    	  
	    
	    for (int j=0; j<histoOldBinning->GetNbinsX(); j++){
		
		varType oldBin = histoOldBinning->GetBinLowEdge(j+1);
		
		if( (oldBin>lowEdge) && (oldBin<=highEdge) ){		   
		    if (rescale) binSum+=histoOldBinning->GetBinContent(j+1) * histoOldBinning->GetBinWidth(j+1);
		    else         binSum+=histoOldBinning->GetBinContent(j+1);
		}
	    }

	    if (rescale) histoNewBinning.Fill(newBinCenter,binSum/newBinWidth);
	    else histoNewBinning.Fill(newBinCenter,binSum);
	}
	histoType * result=(histoType *)histoNewBinning.Clone();
	return result;
    }

    // ===============================================================
    // ===============================================================


  bool plotExists(std::map< TString, std::map <unsigned int, TH1F*> > histo_, const TString plotName, const unsigned int sample)
  {
    // this function checks the existence of an specific
    // entry "histo_[plotName][sample]" in the map "histo_"
    // that contains all 1D plots
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    bool result=false;
    if((histo_.count(plotName)>0)&&(histo_[plotName].count(sample)>0)){
      result=true;
    }
    return result;
  }

  bool plotExists(std::map< TString, std::map <unsigned int, TH2F*> > histo_, const TString plotName, const unsigned int sample)
  {
    // this function checks the existence of an specific
    // entry "histo_[plotName][sample]" in the map "histo_"
    // that contains all 1D plots
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    bool result=false;
    if((histo_.count(plotName)>0)&&(histo_[plotName].count(sample)>0)){
      result=true;
    }
    return result;
  }

  void equalReBinTH1(const int reBinFactor, std::map< TString, std::map <unsigned int, TH1F*> > histo_, const TString plotName, const unsigned int sample)
  {
    // this uses an equal rebinning (factor "reBinFactor")
    // for an specific entry "histo_[plotName][sample]"
    // in the map "histo_" that contains all 1D plots
    // modified quantities: histo_[plotName][sample]
    // used functions: plotExists
    // used enumerators: none

    if(plotExists(histo_, plotName, sample)){
      (histo_[plotName][sample])->Rebin(reBinFactor);
    }

  }
  
  std::map<TString, std::vector<double> > makeVariableBinning(bool addCrossCheckVariables=false)
    {
      // this function creates a map with the hard coded
      // bin values for variable binning
      // NOTE: it is important to quote the overflow bin
      // of the initial binning as last bin here!!!
      // otherwise dividing per binwidth function
      // might later give nan and histo wouldn't be drawn
      // modified quantities: none
      // used functions: none
      // used enumerators: none

      std::map<TString, std::vector<double> > result;
      std::vector<double> bins_;

      // pt(top)
      double topPtBins[]={0.0, 60.0, 100.0, 150.0, 200.0, 260.0, 320.0, 400.0, 500.0};
      //double topPtBins[]={0.0, 50.0, 85.0, 120.0, 160.0, 205.0, 255.0, 320.0, 400.0, 500.0};//alternative: smaller bins
      bins_.insert( bins_.begin(), topPtBins, topPtBins + sizeof(topPtBins)/sizeof(double) );
      result["topPt"        ] = bins_;
      result["topPtLead"    ] = bins_;
      result["topPtSubLead" ] = bins_;
      result["topPtTtbarSys"] = bins_;
      if (addCrossCheckVariables){
	result["topPtPlus" ] = bins_;
	result["topPtMinus"] = bins_;
      }
      bins_.clear();

      // y(top)
      double topYBins[]={-2.5, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.5}; 
      // PAS binning: double topYBins[]={-5., -2.5, -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.5, 5.};
      bins_.insert( bins_.begin(), topYBins, topYBins + sizeof(topYBins)/sizeof(double) );
      result["topY"   ]   = bins_;
      if (addCrossCheckVariables){
	result["topYPlus"]  = bins_;
	result["topYMinus"] = bins_;
      }
      bins_.clear();

      // pt(ttbar)
      double ttbarPtBins[]={0.0, 20.0, 45.0, 75.0, 120.0, 190.0, 300.0};
      // PAS binning: double ttbarPtBins[]={0., 20., 60., 110., 200., 300.}; // PAS Binning
      bins_.insert( bins_.begin(), ttbarPtBins, ttbarPtBins + sizeof(ttbarPtBins)/sizeof(double) );
      result["ttbarPt"]=bins_;
      bins_.clear();

      // y(ttbar)
      // old: double ttbarYBins[]={-5., -1.3, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.3, 5.};
      double ttbarYBins[]={-2.5, -1.3, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.3, 2.5};
      // PAS binning: double ttbarYBins[]={-5., -1.3, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.3, 5.};
      bins_.insert( bins_.begin(), ttbarYBins, ttbarYBins + sizeof(ttbarYBins)/sizeof(double) );
      result["ttbarY"]=bins_;
      bins_.clear();

      // m(ttbar)
      double ttbarMassBins[]={345.0, 400.0, 470.0, 550.0, 650.0, 800.0, 1100.0, 1600.0};
      //double ttbarMassBins[]={250.0, 450.0, 550.0, 700.0, 950.0, 2700.0};// ATLAS binning 2500 -> 2700
      // First option: double ttbarMassBins[]={0.0, 345.0, 400.0, 470.0, 550.0, 650.0, 800.0, 1200.0};  
      // Korea:        double ttbarMassBins[]={0.0, 345.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 1200.0}; 
      // PAS binning:  double ttbarMassBins[]={0., 345., 410., 480., 580., 750., 1200.};
      bins_.insert( bins_.begin(), ttbarMassBins, ttbarMassBins + sizeof(ttbarMassBins)/sizeof(double) );
      result["ttbarMass"]=bins_;
      bins_.clear();

      // phistar(t,tbar)
      double ttbarPhiStarBins[]={0., 0.04, 0.13, 0.38, 2.0};
      bins_.insert( bins_.begin(), ttbarPhiStarBins, ttbarPhiStarBins + sizeof(ttbarPhiStarBins)/sizeof(double) );
      result["ttbarPhiStar"]=bins_;
      bins_.clear();

      // pt(lepton)
      // double lepPtBins[]={30., 35., 40., 45., 50., 60., 70., 80., 100., 120., 150., 200.};
      double lepPtBins[]={30., 37., 45., 55., 68., 80., 100., 135., 200.};
      bins_.insert( bins_.begin(), lepPtBins, lepPtBins + sizeof(lepPtBins)/sizeof(double) );
      result["lepPt"]=bins_;
      bins_.clear();

      // eta(lepton)
      double lepEtaBins[]={-2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
      bins_.insert( bins_.begin(), lepEtaBins, lepEtaBins + sizeof(lepEtaBins)/sizeof(double) );
      result["lepEta"]      = bins_;
      if (addCrossCheckVariables){
	result["lepEtaPlus"]  = bins_;
	result["lepEtaMinus"] = bins_;
      }
      bins_.clear();

      // pt(bjet)
      // double bqPtBins[]={30., 60., 95., 140., 200., 400.};
      double bqPtBins[]={30., 48., 75., 180., 400.};
      bins_.insert( bins_.begin(), bqPtBins, bqPtBins + sizeof(bqPtBins)/sizeof(double) );
      result["bqPt"]=bins_;
      bins_.clear();

      // eta(bquark)
      double bqEtaBins[]={-2.4, -1.5, -1, -0.5, 0., 0.5, 1., 1.5, 2.4};
      bins_.insert( bins_.begin(), bqEtaBins, bqEtaBins + sizeof(bqEtaBins)/sizeof(double) );
      result["bqEta"]=bins_;
      bins_.clear();


      // pt(bbbar)
      //double bbbarPtBins[]={0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 125.0, 160.0, 200.0, 275.0, 500.0};
      double bbbarPtBins[]={0.0, 35.0, 75.0, 115.0, 155.0, 280.0, 500.0};
      bins_.insert( bins_.begin(), bbbarPtBins, bbbarPtBins + sizeof(bbbarPtBins)/sizeof(double) );
      result["bbbarPt"]=bins_;
      bins_.clear();

      // m(bbbar)
      //double bbbarMassBins[]={0.0, 45.0, 60.0, 80.0, 105.0, 135.0, 165.0, 205.0, 255.0, 315.0, 390.0, 485.0, 620.0, 800.0};
      double bbbarMassBins[]={0.0, 85.0, 135.0, 190.0, 255.0, 325.0, 415.0, 505.0, 630.0, 800.0};
      bins_.insert( bins_.begin(), bbbarMassBins, bbbarMassBins + sizeof(bbbarMassBins)/sizeof(double) );
      result["bbbarMass"]=bins_;
      bins_.clear();

      // m(lb)
      // double lbMassBins[]={0.0, 70.0, 95.0, 120.0, 140.0, 500.0};
      double lbMassBins[]={0.0, 90.0, 122., 500.0};
      bins_.insert( bins_.begin(), lbMassBins, lbMassBins + sizeof(lbMassBins)/sizeof(double) );
      result["lbMass"]=bins_;
      bins_.clear();

      // Delta phi(tops)
      double topDelPhiBins[]={0., 2., 2.75, 3., 3.15};
      bins_.insert( bins_.begin(), topDelPhiBins, topDelPhiBins + sizeof(topDelPhiBins)/sizeof(double) );
      result["ttbarDelPhi"]   = bins_;
      bins_.clear();

      // N(jets)
      double NjetsBins[]={3.5, 4.5, 5.5, 6.5, 7.5, 9.5};
      bins_.insert( bins_.begin(), NjetsBins, NjetsBins + sizeof(NjetsBins)/sizeof(double) );
      result["Njets"]=bins_;
      bins_.clear();

      // rhos=tt+1jet/170.GeV
      double rhoSBins[]={0., 0.25, 0.55, 0.7, 1.1};
      bins_.insert( bins_.begin(), rhoSBins, rhoSBins + sizeof(rhoSBins)/sizeof(double) );
      result["rhos"]=bins_;
      bins_.clear();

      return result;
    }

  void setXAxisRange(TH1* his, TString variable)
  {
    // this function restricts the drawn range of the x axis in the differential variable
    // modified quantities: xAxis of his
    // used functions: none
    // used enumerators: none
      
    // restrict axis
    if(     variable.Contains("topPt"    )) his->GetXaxis()->SetRangeUser(0.   , 499. ); 
    else if(variable.Contains("topY"     )) his->GetXaxis()->SetRangeUser(-2.5 , 2.49 );
    else if(variable.Contains("ttbarY"   )) his->GetXaxis()->SetRangeUser(-2.5 , 2.49 );
    else if(variable.Contains("ttbarMass")) his->GetXaxis()->SetRangeUser(346. , 1599.);
    else if(variable.Contains("ttbarDelPhi" )) his->GetXaxis()->SetRangeUser(0., 3.149);
    else if(variable.Contains("ttbarPhiStar")) his->GetXaxis()->SetRangeUser(0., 1.99 );
    else if(variable.Contains("ttbarPt"  )) his->GetXaxis()->SetRangeUser(0.   , 299. );
    else if(variable.Contains("lepPt"    )) his->GetXaxis()->SetRangeUser(30   , 199. );
    else if(variable.Contains("lepEta"   )) his->GetXaxis()->SetRangeUser(-2.1 , 2.09 );
    else if(variable.Contains("bqPt"     )) his->GetXaxis()->SetRangeUser(30.  , 399. );
    else if(variable.Contains("bqEta"    )) his->GetXaxis()->SetRangeUser(-2.4 , 2.39 );
    else if(variable.Contains("lbMass"   )) his->GetXaxis()->SetRangeUser(0.   , 500  );
    else if(variable.Contains("Njets"    )) his->GetXaxis()->SetRangeUser(3.5  , 14.5 );
    else if(variable.Contains("rho"      )) his->GetXaxis()->SetRangeUser(0.   , 1.1  );
  }
    
  template <class T>
    unsigned int positionInVector(std::vector<T> vec_, T object)
    {
      // this function returns the position
      // of element object in vector vec_
      // modified quantities: none
      // used functions: none
      // used enumerators: none

      // loop vector elements
      for(unsigned int element=0; element<vec_.size(); ++element)
	{
	  // ask for requested element
	  if(vec_[element]==object){
	    return element;
	  }
	}
      // return -1 if requested element is not found
      std::cout << "requested element " << object << " is not found in  vector" << std::endl;
      return -1;
    }

  void DivideYieldByEfficiencyAndLumi(TH1F* yield, TH1F* efficiency, double luminosity, bool includeEffError)
  {
    // this function divides yield by efficiency
    // and calculates the correct error in each bin based
    // on gaussian error propagation
    // modified quantities: yield
    // used functions: none
    // used enumerators: none

    // check if #bins are the same
    if(yield->GetNbinsX()!=efficiency->GetNbinsX()){
      std::cout << "#bins in yield and efficiency plots are not the same!" << std::endl;
      exit(1);
    }
    // loop all bins
    for(int bin=1; bin<=yield->GetNbinsX(); ++bin){
      double content  = yield     ->GetBinContent(bin);
      double eff      = efficiency->GetBinContent(bin);
      double binwidth = yield     ->GetBinWidth  (bin);
      if(binwidth==0) binwidth=1;
      double events   = binwidth*content;
      // value (events/binwidth/eff/lumi)
      double xSec = content/(eff*luminosity);
      // error
      double effError=efficiency->GetBinError(bin);
      if(!includeEffError) effError=0;
      double xSecError = 1/(binwidth*luminosity)*sqrt( ((events)/(eff*eff)) + ((events*effError)/(eff*eff))* ((events*effError)/(eff*eff)) );
      // check if efficiency is 0
      if(eff==0){
	xSec=0;
	xSecError=0;
      }
      // set value and error
      yield->SetBinContent(bin, xSec     );
      yield->SetBinError  (bin, xSecError);
    }
  }

 void DivideTwoYields(TH1F* nom, TH1F* denom)
  {
    // this function divides nom by denom
    // and calculates the correct error in each bin based
    // on gaussian error propagation
    // modified quantities: nom
    // used functions: none
    // used enumerators: none

    // check if #bins are the same
    if(nom->GetNbinsX()!=denom->GetNbinsX()){
      std::cout << "#bins in nom and denom plots are not the same!" << std::endl;
      exit(1);
    }
    // loop all bins
    for(int bin=1; bin<=nom->GetNbinsX(); ++bin){
      double Nnom   = nom  ->GetBinContent(bin);
      double Ndenom = denom->GetBinContent(bin);
      // value
      double ratio = Nnom/Ndenom;
      // error
      double ratioError = sqrt((ratio*(1-ratio))/Ndenom);
      // check if denom is 0
      if(Ndenom==0){
	ratio=0;
	ratioError=0;
      }
      // set value and error
      nom->SetBinContent(bin, ratio     );
      nom->SetBinError  (bin, ratioError);
    }
  }

  template <class T>
  void saveToRootFile(const TString& outputFile, const T& object, const bool& overwrite=false, const int& verbose=1, const TString& folder="")
    {
      // this function saves objects of class T
      // (such as TH1) in an output rootfile with name outputFile
      // modified quantities: outputFile
      // used functions: none
      // used enumerators: none


      // check if file is already open
      TFile* openfile = 0;
      TSeqCollection *openfiles= gROOT->GetListOfFiles();
      for(int i=0; i<openfiles->GetSize(); ++i){
	 if(verbose>2) std::cout << openfiles->At(i)->GetName() << std::endl;
	if(TString(openfiles->At(i)->GetName())==outputFile){
	  openfile=(TFile*)(openfiles->At(i));
	  if(verbose>1) std::cout << "file " << outputFile << " already open, use reopen" << std::endl;
	}
      }

      bool saveObject=true;
      // check if file exist
      TFile* file = openfile ? openfile : TFile::Open(outputFile, "UPDATE");
      if(openfile) file->ReOpen("UPDATE");
      // if not exist: create
      if(!file){
	if(verbose>1) std::cout << "file " << outputFile << " does not exist, will be created" << std::endl;
	file = new TFile(outputFile, "RECREATE");
      }
      // check if file is broken
      if(file->IsZombie()){
	std:: cout << "file " << outputFile << " is broken" << std::endl;
	// if broken: don't save object
	saveObject=false;
      }
      else{
	// if not broken: open file
	file->cd();
	// create folder if it does not exist
	if(folder!=""){
	  if(verbose>1) std::cout << "check existence of directory " << folder << std::endl;
	  // see how many subfolders are within folder
	  unsigned int Nfolders=folder.CountChar(*"/")+1;
	  if(verbose>1) std::cout << "these are " << Nfolders << " subdirectories" << std::endl;
	  // loop subdirectories
	  for(unsigned int subfolderNumber=1; subfolderNumber<=Nfolders; ++subfolderNumber){
	    TString subfolder= getStringEntry(folder, subfolderNumber, "/");
	    if(gDirectory->GetDirectory(subfolder)==0){
	      if(verbose>1) std::cout << "subfolder " << subfolder  << " not existing - will create it" << std::endl;
	      gDirectory->mkdir(subfolder);
	    }
	    else if(verbose>1) std::cout << "subfolder " << subfolder  << " is already existing" << std::endl;
	    // go to directory
	    gDirectory->cd(subfolder);
	  }
	}
	// if you don't want to overwrite, check if object exists
	int count=-1;
	if(!overwrite){
	  TString saveObjectName=(TString)object->GetTitle();
	  if(verbose>1) std::cout << "searching for object " << saveObjectName << std::endl;
	  // loop all objects in file
	  TList * list = gDirectory->GetListOfKeys();
	  if(list->At(0)){
	    do{
	      ++count;
	      TObject *folderObject = list->At(count);
	      TString folderObjectName = (TString)folderObject->GetName();
	      if(verbose>1) std::cout << "candidate #" << count+1 << ": " << folderObjectName << std::endl;
	      // check if object you want to save is already existing
	      // by comparing the names
	      if(folderObjectName==saveObjectName){
		if(verbose>1){
		  std::cout << "already exists in file " << outputFile << std::endl;
		  std::cout << "will keep the old one!" << std::endl << std::endl;
		}
		saveObject=false;
		break;
	      }
	    }while( list->At( count ) != list->Last());
	  }
	}
      }
      // save object
      if(saveObject){
	if(verbose>0){
	  std::cout << "saving object " << (TString)object->GetTitle();
	  std::cout << " to file " << outputFile;
	  if(folder!="") std::cout << " ( folder " << folder << " )";
	  std::cout << std::endl << std::endl;
	}
	// overwrite existing
	object->Write(object->GetTitle(), TObject::kOverwrite);
      }
      // close file
      if(!openfile) file->Close();
    }

  template <class T>
    void saveToRootFileAll(const TString& outputFile, const std::vector< std::pair< T*,TString > > &objects_, const bool& overwrite=false, const int& verbose=1)
    {
      // this function saves all objects of class T in saveCanvas_
      // (such as TH1) in an output rootfile with name outputFile
      // modified quantities: outputFile
      // used functions: none
      // used enumerators: none

      bool saveObject=true;
      // check if file exist
      TFile* file = TFile::Open(outputFile, "UPDATE");
      // if not exist: create
      if(!file){
	if(verbose>1) std::cout << "file " << outputFile << " does not exist, will be created" << std::endl;
	file = new TFile(outputFile, "RECREATE");
      }
      // check if file is broken
      if(file->IsZombie()){
	std:: cout << "file " << outputFile << " is broken" << std::endl;
	// if broken: don't save object
	saveObject=false;
      }
      else{
	// if not broken: open file
	file->cd();
	// loop list of plots to be saved
	for(unsigned int idx=0; idx<objects_.size(); ++idx){  
	  //get object and folder
	  TString folder=(TString)(objects_[idx].second);
	  T* object=(T*)(objects_[idx].first);
	  // create folder if it does not exist
	  if(folder!=""){
	    if(verbose>1) std::cout << "check existence of directory " << folder << std::endl;
	    // see how many subfolders are within folder
	    unsigned int Nfolders=folder.CountChar(*"/")+1;
	    if(verbose>1) std::cout << "these are " << Nfolders << " subdirectories" << std::endl;
	    // loop subdirectories
	    for(unsigned int subfolderNumber=1; subfolderNumber<=Nfolders; ++subfolderNumber){
	      if(verbose>1) std::cout << "#" << subfolderNumber << std::endl;
	      TString subfolder= getStringEntry(folder, subfolderNumber, "/");
	      if(gDirectory->GetDirectory(subfolder)==0){
		if(verbose>1) std::cout << "subfolder " << subfolder  << " not existing - will create it" << std::endl;
		gDirectory->mkdir(subfolder);
	      }
	      else if(verbose>1) std::cout << "subfolder " << subfolder  << " is already existing" << std::endl;
	      // go to directory
	      gDirectory->cd(subfolder);
	    }
	  }
	  // if you don't want to overwrite, check if object exists
	  int count=-1;
	  if(!overwrite){
	    TString saveObjectName=(TString)object->GetTitle();
	    if(verbose>1) std::cout << "searching for object " << saveObjectName << std::endl;
	    // loop all objects in file
	    TList * list = gDirectory->GetListOfKeys();
	    while( list->At( count+1 ) != list->Last()){
	      ++count;
	      TObject *folderObject = list->At(count);
	      TString folderObjectName = (TString)folderObject->GetName();
	      if(verbose>1) std::cout << "candidate #" << count+1 << ": " << folderObjectName << std::endl;
	      // check if object you want to save is already existing
	      // by comparing the names
	      if(folderObjectName==saveObjectName){
		if(verbose>1){
		  std::cout << "already exists in file " << outputFile << std::endl;
		  std::cout << "will keep the old one!" << std::endl << std::endl;
		}
		saveObject=false;
		break;
	      }
	    }
	  }
	  // save object
	  if(saveObject){
	    if(verbose>0){
	      std::cout << "saving object " << (TString)object->GetTitle();
	      std::cout << " to file " << outputFile;
	      if(folder!="") std::cout << " ( folder " << folder << " )";
	      std::cout << std::endl << std::endl;
	    }
	    // overwrite existing
	    object->Write(object->GetTitle(), TObject::kOverwrite);
	  }
	  // return to main folder
	  for(int nfolder=1; nfolder<=folder.CountChar(*"/")+1; ++nfolder){
	    gDirectory->cd("..");
	  }
	}
      }
      // close file
      file->Close();
    }

  void CloseOpenFiles(int verbose=0){   
    // this function closes all open rootfiles
    // use it at the begin of all macros to prevent problems with saving
    // modified quantities: none
    // used functions: none
    // used enumerators: none 

    TSeqCollection *openfiles= gROOT->GetListOfFiles();
    for(int i=0; i<openfiles->GetSize(); ++i){
      if(verbose>2) std::cout << openfiles->At(i)->GetName() << std::endl;
      ((TFile*)(openfiles->At(i)))->Close();
    }
  }

  int drawRatio(const TH1* histNumerator, TH1* histDenominator, const Double_t& ratioMin, const Double_t& ratioMax, TStyle myStyle, int verbose=0, const std::vector<double> err_=std::vector<double>(0), TString ratioLabelNominator="N_{data}", TString ratioLabelDenominator="N_{MC}", TString ratioDrawOption="p e X0", int ratioDrawColor=kBlack, bool error=true, double ratioMarkersize=1.2, int ndiv=505, bool invert=false, bool drawOne=true, bool linFit=false)
  {
    // this function draws a pad with the ratio of 'histNumerator' and 'histDenominator'
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // per default only the gaussian error of the 'histNumerator' is considered:
    // (error(bin i) = sqrt(histNumerator->GetBinContent(i))/histDenominator->GetBinContent(i))
    // if 'err_' is present and its size equals the number of bins in the histos,
    // its values are considered as error for the ratio (if error=true) 
    // (first entry in err_ - index 0 - corresponds to bin 1)
    // 'ndiv' defines the Ndivisions for the ratio y axis
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // invert: draw histDenominator/histNumerator instead of the other way round 
    //         (this allows to switch fast without additional fixing of axes labels etc)
    // NOTE: ratioLabelDenominator and ratioLabelNominator are assumed to correspond to 
    //       histNumerator and histDenominator, if invert is used, they are automatically 
    //       switched, wrong treatment leads to wrong errors for empty err_ and error=true
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    
    // invert ==false: build histNumerator/histDenominator
    // invert ==true : build histDenominator/histNumerator
    // -> build in any case tempNumerator/tempDenominator
    TH1F* tempNumerator   = invert ? (TH1F*)histDenominator->Clone() : (TH1F*)histNumerator  ->Clone();
    TH1F* tempDenominator = invert ? (TH1F*)histNumerator  ->Clone() : (TH1F*)histDenominator->Clone();
    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
      std::cout << " ERROR when calling drawRatio - histos have different number of bins" << std::endl;
      return -1;
    }
    if(verbose>1){
      std::cout << "building ratio plot of " << tempNumerator->GetName();
      std::cout << " and " << tempDenominator->GetName() << std::endl;
    }
    // create ratio
    TH1F* ratio = (TH1F*)histNumerator->Clone();
    // separate calculation of each bin
    for(int bin=1; bin<=tempNumerator->GetNbinsX(); bin++){
      double num  =tempNumerator  ->GetBinContent(bin);
      double denom=tempDenominator->GetBinContent(bin);
      if(denom==0.) denom=0.00001;
      ratio->SetBinContent(bin, num/denom);
    }
    // calculate error for ratio
    // a) from err_
    if(err_.size()==(unsigned int)histNumerator->GetNbinsX()){
      if(verbose>0) std::cout << "ratio error from vector" << std::endl;
      for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
	ratio->SetBinError(bin, err_[bin-1]);
      }
    }
    // b) from data statistics
    else if (error==true){
      // get data (condition is true if data is finally used as denominator)
      TH1F* data=tempNumerator;
      if((!invert&& (ratioLabelDenominator.Contains("data")||ratioLabelDenominator.Contains("Data"))&&!(ratioLabelNominator.Contains("data")||ratioLabelNominator.Contains("Data")))||( invert&&!(ratioLabelDenominator.Contains("data")||ratioLabelDenominator.Contains("Data"))&& (ratioLabelNominator.Contains("data")||ratioLabelNominator.Contains("Data")))){
	data=tempDenominator;
      }
      // uncertainty from gaussian propagated error of data
      //     = ratio*relative poissonean error of Ndata -> ratio(data,MC)*sqrt(Ndata)/Ndata
      for(int bin=1; bin<=tempDenominator->GetNbinsX(); bin++){
	ratio->SetBinError(bin, ratio->GetBinContent(bin)*(sqrt(data->GetBinContent(bin)))/data->GetBinContent(bin));
      }
    }
    // get some values from old pad
    //Int_t    logx = gPad->GetLogx();
    //Double_t left = gPad->GetLeftMargin();
    //Double_t right = gPad->GetRightMargin();
      
    Int_t    logx  = myStyle.GetOptLogx();
    //Double_t left  = myStyle.GetPadLeftMargin();
    //Double_t right = myStyle.GetPadRightMargin();
    Double_t left  = gPad->GetLeftMargin();
    Double_t right = gPad->GetRightMargin();
    if(!histDenominator->GetXaxis()->GetNoExponent()&&TString(histDenominator->GetName()).Contains("shift")&&(TString(histDenominator->GetName()).Contains("Eta")||TString(histDenominator->GetName()).Contains("Phi"))&&!(TString(histDenominator->GetName()).Contains("Nu"))) right=0.11;

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.36;
    double scaleFactor = 1./(canvAsym*ratioSize);
    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);

    if(ratioDrawOption.Contains("same")&&verbose>1) std::cout << "drawing into same ratio" << std::endl;
    if(!ratioDrawOption.Contains("same")){
      TPad *whitebox3 = new TPad("whitebox3","",0.2,0.,1.,1.);
      whitebox3->SetFillColor(10);
      whitebox3->SetFillStyle(1001);
      whitebox3->SetBorderSize(0);
      whitebox3->SetBorderMode(0);
      // create new pad for ratio plot
      TPad *rPad;
      rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
#ifdef DILEPTON_MACRO
      rPad->SetFillColor(10);
#else
      rPad->SetFillStyle(0);
      rPad->SetFillColor(0);
#endif
      rPad->SetBorderSize(0);
      rPad->SetBorderMode(0);
      rPad->Draw();
      // configure ratio plot
      rPad->cd();
      //whitebox3->Draw();
      rPad->SetLogy(0);
      rPad->SetLogx(logx);
      rPad->SetTicky(1);
      rPad->SetTopMargin(0.0);
      rPad->SetBottomMargin(0.15*scaleFactor);
      rPad->SetRightMargin(right);
      // draw grid
      rPad->SetGrid(0,1);

      ratio->SetStats(kFALSE);
      ratio->SetTitle("");
      ratio->SetMaximum(ratioMax);
      ratio->SetMinimum(ratioMin);
      // configure axis of ratio plot
      ratio->GetXaxis()->SetTitleSize(histDenominator->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
      ratio->GetXaxis()->SetTitleOffset(histDenominator->GetXaxis()->GetTitleOffset()*0.9);
      ratio->GetXaxis()->SetLabelSize(histDenominator->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
      ratio->GetXaxis()->SetTitle(histDenominator->GetXaxis()->GetTitle());
      ratio->GetXaxis()->SetNdivisions(histDenominator->GetNdivisions());
      ratio->GetYaxis()->CenterTitle();
      TString ratioAxisTitle="";
      //if(ratioLabelDenominator!=""&&ratioLabelNominator!="") ratioAxisTitle="#frac{"+ratioLabelNominator+"}{"+ratioLabelDenominator+"}";
      if(ratioLabelDenominator!=""&&ratioLabelNominator!="") ratioAxisTitle="#frac{"+(!invert ? ratioLabelDenominator : ratioLabelNominator)+"}{"+(!invert ? ratioLabelNominator : ratioLabelDenominator)+"}";
      else if(ratioLabelDenominator!="") ratioAxisTitle=ratioLabelDenominator;
      else if(ratioLabelNominator  !="") ratioAxisTitle=ratioLabelNominator  ;
      ratio->GetYaxis()->SetTitle(ratioAxisTitle);
      ratio->GetYaxis()->SetTitleSize(histDenominator->GetYaxis()->GetTitleSize()*scaleFactor);
      ratio->GetYaxis()->SetTitleOffset(histDenominator->GetYaxis()->GetTitleOffset()/scaleFactor);
      ratio->GetYaxis()->SetLabelSize(histDenominator->GetYaxis()->GetLabelSize()*scaleFactor);
      ratio->GetYaxis()->SetLabelOffset(histDenominator->GetYaxis()->GetLabelOffset()*3.3);
      ratio->GetYaxis()->SetTickLength(0.03);
      ratio->GetYaxis()->SetNdivisions(ndiv);
      ratio->GetXaxis()->SetRange(histDenominator->GetXaxis()->GetFirst(), histDenominator->GetXaxis()->GetLast());
      ratio->GetXaxis()->SetNoExponent(histDenominator->GetXaxis()->GetNoExponent());
      // delete axis of initial plot
      histDenominator->GetXaxis()->SetLabelSize(0);
      histDenominator->GetXaxis()->SetTitleSize(0);
      rPad->RedrawAxis("g");
    }
    // draw ratio plot
    ratio->SetLineWidth(2);
    ratio->SetLineColor(ratioDrawColor);
    ratio->SetMarkerSize(ratioMarkersize);
    ratio->SetMarkerColor(ratioDrawColor);
    //gPad->Update();
    ratio->Draw(ratioDrawOption);
    if(linFit)  ratio->Fit("pol1", "q", "", histDenominator->GetBinLowEdge(histDenominator->GetXaxis()->GetFirst()), histDenominator->GetBinLowEdge(histDenominator->GetXaxis()->GetLast()+1));
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis("g");

    if(!ratioDrawOption.Contains("same")){
      // draw a horizontal lines on a given histogram
      // a) at 1
      Double_t xmin = ratio->GetXaxis()->GetXmin();
      Double_t xmax = ratio->GetXaxis()->GetXmax();
      TString height = ""; height += 1;
      TF1 *f = new TF1("f", height, xmin, xmax);
      f->SetLineStyle(1);
      f->SetLineWidth(1);
      f->SetLineColor(kBlack);
      if(drawOne) f->Draw("L same");
      // b) at upper end of ratio pad
      TString height2 = ""; height2 += ratioMax;
      TF1 *f2 = new TF1("f2", height2, xmin, xmax);
      f2->SetLineStyle(1);
      f2->SetLineWidth(1);
      f2->SetLineColor(kBlack);
      f2->Draw("L same");
      // c) at 0.5
      TString height05 = ""; height05 += 0.5;
      TF1 *f05 = new TF1("f05", height05, xmin, xmax);
      f05->SetLineStyle(3);
      f05->SetLineWidth(1);
      f05->SetLineColor(kBlack);
      //f05->Draw("L same");
      // d) at 1.5
      TString height15 = ""; height15 += 1.5;
      TF1 *f15 = new TF1("f15", height15, xmin, xmax);
      f15->SetLineStyle(3);
      f15->SetLineWidth(1);
      f15->SetLineColor(kBlack);
      //f15->Draw("L same");
      // e) at 0.8
      TString height08 = ""; height08 += 0.8;
      TF1 *f08 = new TF1("f08", height08, xmin, xmax);
      f08->SetLineStyle(3);
      f08->SetLineWidth(1);
      f08->SetLineColor(kBlack);
      //f08->Draw("L same");
      // f) at 1.2
      TString height12 = ""; height12 += 1.2;
      TF1 *f12 = new TF1("f12", height12, xmin, xmax);
      f12->SetLineStyle(3);
      f12->SetLineWidth(1);
      f12->SetLineColor(kBlack);
      //f12->Draw("L same");
    }
    return 0;    
  }

  unsigned int theoryColor(TString theo="madgraph"){
    // this function returns the default color for a given theory specified by "theo"
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    if(theo.Contains("mcatnlo")||theo.Contains("MC@NLO")||theo.Contains("mc@nlo")||theo.Contains("McAtNlo")||theo.Contains("Mc@Nlo")) return constMcatnloColor; 
    if(theo.Contains("powheg" )||theo.Contains("Powheg")||theo.Contains("POWHEG")||theo.Contains("PowHeg")){
      if(theo.Contains("herwig" )||theo.Contains("Herwig")||theo.Contains("HERWIG")) return constPowhegColor2; 
      else return constPowhegColor; 
    }
    if(theo.Contains("nnlo")){
      if(theo.Contains("ahrens")||theo.Contains("ttbarMass")||theo.Contains("ttbarPt")) return constNnloColor2;
      return constNnloColor;
    }
    if(theo.Contains("data")) return kBlack;
    return constMadgraphColor;
  }

  unsigned int theoryStyle(TString theo="madgraph"){
    // this function returns the default line style for a given theory specified by "theo"
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    if(theo.Contains("mcatnlo")||theo.Contains("MC@NLO")||theo.Contains("mc@nlo")||theo.Contains("McAtNlo")||theo.Contains("Mc@Nlo")) return constMcatnloStyle;
    if(theo.Contains("powheg" )||theo.Contains("Powheg")||theo.Contains("POWHEG")||theo.Contains("PowHeg")){
      if(theo.Contains("herwig" )||theo.Contains("Herwig")||theo.Contains("HERWIG")) return constPowhegStyle2; 
      else  return constPowhegStyle; 
    }
    if(theo.Contains("nnlo")){
      if(theo.Contains("ahrens")||theo.Contains("ttbarMass")||theo.Contains("ttbarPt")) return constNnloStyle2;
      return constNnloStyle;
    }
    if(theo.Contains("nnlo"    )||theo.Contains("kidonakis")) return constNnloStyle;
    if(theo.Contains("madgraph")||theo.Contains("Madgraph")||theo.Contains("MadGraph")) return 1;
    return 42; 
  }
  
  TString xSecLabelName(TString variable="", bool noUnit=false){
    // this function returns the x axis label for a given quantity specified by "variable"
    // if noUnit==true -> skip unit
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    
    TString strUnitGeV = " #left[GeV#right]";
    if(noUnit) strUnitGeV="";

    if(variable == "topPt"        ) return "p_{T}^{t}"+strUnitGeV;
    if(variable == "topPtTtbarSys") return "p_{T}^{t} #scale[0.8]{(t#bar{t} "+cmsExt+")}"+strUnitGeV;
    if(variable == "topPtLead"    ) return "p_{T}^{lead t}"+strUnitGeV;
    if(variable == "topPtSubLead" ) return "p_{T}^{sublead t}"+strUnitGeV;
    if(variable == "topPtPlus"    ) return "p_{T}^{t}"+strUnitGeV;
    if(variable == "topPtMinus"   ) return "p_{T}^{#bar{t}}"+strUnitGeV;
    if(variable == "topY"         ) return "y^{t}";
    if(variable == "topYPlus"     ) return "y^{t}";
    if(variable == "topYMinus"    ) return "y^{#bar{t}}";
    if(variable == "ttbarPt"      ) return "p_{T}^{t#bar{t}}"+strUnitGeV;
    if(variable == "ttbarY"       ) return "y^{t#bar{t}}";
    if(variable == "ttbarMass"    ) return "m^{t#bar{t}}"+strUnitGeV;
    if(variable == "ttbarPhiStar" ) return "#Phi^{#lower[-0.9]{* }}(t,#bar{t})";
    if(variable == "ttbarDelPhi"  ) return "#Delta#phi(t,#bar{t})";
    if(variable == "lepPt"        ) return "p_{T}^{l}"+strUnitGeV;
    if(variable == "lepEta"       ) return "#eta^{l}"; 
    if(variable == "lepEtaPlus"   ) return "#eta^{l^{+}}";
    if(variable == "lepEtaMinus"  ) return "#eta^{l^{-}}";
    if(variable == "bqPt"         ) return "p_{T}^{b}"+strUnitGeV;
    if(variable == "bqEta"        ) return "#eta^{b}";
    if(variable == "bbbarPt"      ) return "p_{T}^{b#bar{b}}"+strUnitGeV;
    if(variable == "bbbarMass"    ) return "m^{b#bar{b}}"+strUnitGeV;
    if(variable == "lbMass"       ) return "m^{lb}"+strUnitGeV;
    if(variable == "Njets"        ) return "N_{jets} (p_{T}>30 GeV, |#eta|<2.4)";
    if(variable == "rhos"         ) return "#rho_{S}= #frac{#scale[0.55]{2#upoint170 GeV}}{#scale[0.55]{m(t#bar{t}+1jet)}}";
    return "Default Label for variable "+variable;
  }

  double xSecMaximum(TString variable=""){
    // this function returns the maximum of the y axis 
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    if(variable.Contains("topPtTtbarSys")) return 9.0E-03;
    if(variable.Contains("topPtLead"    )) return 8.5E-03;
    if(variable.Contains("topPtSubLead" )) return 9.0E-03;
    if(variable.Contains("topPt"        )) return 8.5E-03;
    if(variable.Contains("topPtPlus"    )) return -1.0;
    if(variable.Contains("topPtMinus"   )) return -1.0;
    if(variable.Contains("topY"         )) return 9.0E-01;
    if(variable.Contains("topYPlus"     )) return -1.0;
    if(variable.Contains("topYMinus"    )) return -1.0;
    if(variable.Contains("ttbarPt"      )) return 2.0E-02;
    if(variable.Contains("ttbarY"       )) return 9.0E-01;
    if(variable.Contains("ttbarMass"    )) return 1.0;
    if(variable.Contains("ttbarPhiStar" )) return 10.0;
    if(variable.Contains("ttbarDelPhi"  )) return 3.0;
    if(variable.Contains("lepPt"        )) return 3.0E-02;
    if(variable.Contains("lepEta"       )) return 0.55;
    if(variable.Contains("lepEtaPlus"   )) return -1.0;
    if(variable.Contains("lepEtaMinus"  )) return -1.0;
    if(variable.Contains("bqPt"         )) return 1.8E-02;
    if(variable.Contains("bqEta"        )) return 0.5;
    if(variable.Contains("bbbarPt"      )) return 1.2E-02;
    if(variable.Contains("bbbarMass"    )) return 8.0E-03;
    if(variable.Contains("lbMass"       )) return 1.2E-02;
    if(variable.Contains("Njets"        )) return 0.75;
    if(variable.Contains("rho"          )) return 5.5;
    return -1.0;
  }

  TCanvas* drawFinalResultRatio(TH1F* histNumeratorData, const Double_t& ratioMin, const Double_t& ratioMax, TStyle myStyle, int verbose=0, std::vector<TH1F*> histDenominatorTheory_=std::vector<TH1F*>(0), TCanvas* canv=0, double rangeMin=-1., double rangeMax=-1., TGraphAsymmErrors* histStatData=0, bool addXBinGrid=true, bool addYGrid=false, bool bands=true, unsigned int ndiv=505)
  {
    // this function draws a pad with the ratio "histNumeratorData" over "histDenominatorTheoryX" 
    // for up to five specified theory curves, using "histNumeratorDataDown" and "histNumeratorDataUp"
    // the range of the ratio is 'ratioMin' to 'ratioMax'
    // to the systematic variation "sys" of the enumerator "systematicVariation"
    // the uncertainty of "histNumeratorData" only is considered for the uncertainty of the ratio:
    // ( error(bin i) = sqrt(1/histNumerator->GetBinContent(i)) / histDenominator->GetBinContent(i) )
    // 'addXBinGrid': a vertical grid is drawn at the bin boundaries if true
    // 'histStatData': statistical errors are drawn in the ratio and the BCC x
    //                 position from this object is used also for the ratio plots
    // 'bands'       : draw data uncertainty as errorbands if set to true
    // NOTE: x Axis is transferred from histDenominator to the bottom of the canvas
    // modified quantities: none
    // used functions: setXAxisRange, theoryColor, drawLine
    // used enumerators: none
    
    // some basic printout
    if(histDenominatorTheory_.size()>0&&histNumeratorData){
      if(verbose>0) std::cout << "calling drawFinalResultRatio for " << histNumeratorData->GetName() << " and " << histDenominatorTheory_.size() << " theory curves" << std::endl;
      // check that histos have the same binning
      for(unsigned int nTheory=0; nTheory<histDenominatorTheory_.size(); ++nTheory){
	if(!histDenominatorTheory_[nTheory]||histNumeratorData->GetNbinsX()!=histDenominatorTheory_[nTheory]->GetNbinsX()){
	  std::cout << "error when calling drawRatio - data and theory" << nTheory+1 << " histo have different number of bins!" << std::endl;
	  std::cout << histNumeratorData->GetNbinsX() << " (data) vs. " << histDenominatorTheory_[nTheory]->GetNbinsX() << " (theory)" << std::endl;
	  std::cout << "set this ratio to 1!" << std::endl;
	  // if not: theory is set to data -> ratio=1
	  histDenominatorTheory_[nTheory]=(TH1F*)histNumeratorData->Clone("WARNINGerrorInRatio");
	}
      }
      if(verbose>1) std::cout << "range: " << histNumeratorData->GetNbinsX() << " bins from " <<  histNumeratorData->GetBinLowEdge(histNumeratorData->GetXaxis()->GetFirst()) << " to " << histNumeratorData->GetBinLowEdge(histNumeratorData->GetXaxis()->GetLast()+1) << std::endl;

      // enter canvas
      if(!canv){
	std::cout << "canvas handed over in drawFinalResultRatio as canv is empty" << std::endl;
	return canv;
      }	       
      canv->cd();
      canv->Draw();
      //std::cout << "width =" << canv->GetWindowWidth()  << std::endl;
      //std::cout << "height=" << canv->GetWindowHeight() << std::endl;

      //canv->SetGrayscale(); // FIXME: for black-white compatibility check

      // create ratios
      std::vector<TH1*> ratio_;
      // sort for drawing order: MADGRAPH, MC@NLO, POWHEG, POWHEG+HERWIG, higher order theory calculation
      std::vector<TH1F*> histDenominatorTheoryOrdered_=histDenominatorTheory_;
/*       unsigned int count =0; */
/*       while(histDenominatorTheoryOrdered_.size()!=histDenominatorTheory_.size()){ */
/* 	++count; */
/* 	for(unsigned int nTheory=0; nTheory<histDenominatorTheory_.size(); ++nTheory){  */
/* 	  TString theo=(TString)histDenominatorTheory_[nTheory]->GetName(); */
/* 	  if((count==1)&&!(theo.Contains("errorBand"))) histDenominatorTheoryOrdered_.push_back((TH1*)(histDenominatorTheory_[nTheory]->Clone())); */
/* 	  if((count==2)&&!(theo.Contains("mcatnlo")||theo.Contains("MC@NLO")||theo.Contains("mc@nlo")||theo.Contains("McAtNlo")||theo.Contains("Mc@Nlo")||theo.Contains("powheg" )||theo.Contains("Powheg")||theo.Contains("POWHEG")||theo.Contains("PowHeg")||theo.Contains("nnlo"))) histDenominatorTheoryOrdered_.push_back((TH1*)(histDenominatorTheory_[nTheory]->Clone())); */
/* 	  else if((count==3)&&!(theo.Contains("errorBand"))&&(theo.Contains("mcatnlo")||theo.Contains("MC@NLO")||theo.Contains("mc@nlo")||theo.Contains("McAtNlo")||theo.Contains("Mc@Nlo"))) histDenominatorTheoryOrdered_.push_back((TH1*)(histDenominatorTheory_[nTheory]->Clone())); */
/* 	  else if((count==4)&&(theo.Contains("powheg" )||theo.Contains("Powheg")||theo.Contains("POWHEG")||theo.Contains("PowHeg"))) histDenominatorTheoryOrdered_.push_back((TH1*)(histDenominatorTheory_[nTheory]->Clone())); */
/* 	  else if((count>=5)&&(theo.Contains("nnlo"))) histDenominatorTheoryOrdered_.push_back((TH1*)(histDenominatorTheory_[nTheory]->Clone())); */
/* 	  if(count>histDenominatorTheory_.size()){ std::cout << "ERROR in drawFinalResultRatio when trying to order theory curves- new or unknown curve???" << std::endl; exit(0);} */
/* 	} */
/*       } */
      
      for(unsigned int nTheory=0; nTheory<histDenominatorTheoryOrdered_.size(); ++nTheory){
	if(verbose>1){
	  std::cout << "ratio for theory #" << nTheory;
	  std::cout << "  (" << histDenominatorTheoryOrdered_[nTheory]->GetName();
	  std::cout << " / " << histNumeratorData->GetName() << ")" << std::endl;
	}
	// clone theory curve
	TH1F* tempRatio = (TH1F*)histDenominatorTheoryOrdered_[nTheory]->Clone();
	tempRatio->Reset("icms");
	// loop bins
	for(int bin=1; bin<=histDenominatorTheoryOrdered_[nTheory]->GetNbinsX(); ++bin){
	  //std::cout << "bin: " << bin << std::endl;
	  if(histNumeratorData->GetBinWidth(bin)!=0.&& histNumeratorData->GetBinContent(bin)!=0){
	    // set content
	    //std::cout << "histDenominatorTheoryOrdered_[nTheory]->GetBinContent(bin): "<< histDenominatorTheoryOrdered_[nTheory]->GetBinContent(bin) << std::endl;
	    //std::cout << "histNumeratorData->GetBinContent(bin): " << histNumeratorData->GetBinContent(bin) << std::endl;
	    tempRatio->SetBinContent(bin, histDenominatorTheoryOrdered_[nTheory]->GetBinContent(bin) / histNumeratorData->GetBinContent(bin));
	    // set error: simple gaussian error propagation neglecting the theory error
	    tempRatio->SetBinError  (bin, ( histDenominatorTheoryOrdered_[nTheory]->GetBinContent(bin) / ( histNumeratorData->GetBinContent(bin)* histNumeratorData->GetBinContent(bin)) ) * histNumeratorData->GetBinError(bin));
	  }
	  // set technical bins with width==0 to zero
	  else{
	    tempRatio->SetBinError  (bin, 0.5);
	    tempRatio->SetBinContent(bin, 1  );
	  }
	  if(verbose>1) {
	    std::cout << histDenominatorTheoryOrdered_[nTheory]->GetName() << " bin " << bin << ":" << tempRatio->GetBinContent(bin) << "+-" << tempRatio->GetBinError(bin) << std::endl;
	    std::cout << "rel err ratio: " << tempRatio->GetBinError(bin)/tempRatio->GetBinContent(bin) << std::endl;
	    std::cout << "rel err data:  " << histNumeratorData->GetBinError(bin)/histNumeratorData->GetBinContent(bin) << std::endl;
	  }
	}
	ratio_.push_back((TH1F*)tempRatio->Clone());
      }
      if(verbose>3) std::cout << "debug A" << std::endl;
      // get some values from style
      Int_t    logx  = myStyle.GetOptLogx();
      Double_t left  = myStyle.GetPadLeftMargin();
      Double_t right = myStyle.GetPadRightMargin();
      setXAxisRange(histNumeratorData, (TString)histNumeratorData->GetName());
      if(rangeMin!=-1.&&rangeMax!=-1.){
	histNumeratorData->GetXaxis()->SetRange(histNumeratorData->FindBin(rangeMin),histNumeratorData->FindBin(rangeMax));
      }
      if(verbose>3) std::cout << "debug B" << std::endl;
      //Double_t xmin = histNumeratorData->GetXaxis()->GetXmin();
      //Double_t xmax = histNumeratorData->GetXaxis()->GetXmax();
      // y:x size ratio for canvas
      double canvAsym = 4./3.;
      // ratio size of pad with plot and pad with ratio
      double ratioSize = 0.36;
      double scaleFactor = 1./(canvAsym*ratioSize);
      // change old pad
      gPad->SetBottomMargin(ratioSize);
      gPad->SetRightMargin(right);
      gPad->SetLeftMargin(left);
      gPad->SetBorderMode(0);
      gPad->SetBorderSize(0);
      gPad->SetFillColor(10);
      if(verbose>3) std::cout << "debug C" << std::endl;
      // create pad to hide old axis
      TPad *whitebox = new TPad("whitebox","",0.97*left,ratioSize-0.2,1,ratioSize+0.001);
      whitebox->SetFillColor(10);
      whitebox->SetFillStyle(1001);
      whitebox->SetBorderSize(0);
      whitebox->SetBorderMode(0);
      //whitebox->Draw("");
      // create pad to hide the outer left part of the
      // rapidity x-axis without masking the 0 of the y-axis
      TPad *whitebox2 = new TPad("whitebox2","",0.176332,0.144974,0.982758,0.346332);
      whitebox2->SetFillColor(10);
      whitebox2->SetFillStyle(1001);
      whitebox2->SetBorderSize(0);
      whitebox2->SetBorderMode(0);
      //if(((TString)histNumeratorData->GetName()).Contains("Y"))	whitebox2->Draw("");
      // create pad to hide the outer left part of the
      // rapidity x-axis without masking the 0 of the y-axis
      TPad *whitebox3 = new TPad("whitebox3","",0.2,0.,1.,1.);
      whitebox3->SetFillColor(10);
      whitebox3->SetFillStyle(1001);
      whitebox3->SetBorderSize(0);
      whitebox3->SetBorderMode(0);
      // create new pad for ratio plot
      TPad *rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
#ifdef DILEPTON_MACRO
      rPad->SetFillColor(10);
#else
      rPad->SetFillStyle(0);
      rPad->SetFillColor(0);
#endif
      rPad->SetBorderSize(0);
      rPad->SetBorderMode(0);
      rPad->SetLogy(0);
      rPad->SetLogx(logx);
      rPad->SetTicky(1);

      if(verbose>3) std::cout << "debug D" << std::endl;
      // draw grid
      TPad *grid =(TPad*)rPad->Clone("grid");
      grid->SetGrid(0,1);
      if(addYGrid){
	grid->Draw("");
	rPad->SetGrid(0,1);
      }
      rPad->Draw("");
      rPad->cd();
      //whitebox3->Draw();
      TH1F* one2=0;
      TGraphAsymmErrors* errorband2=0;
      TGraphAsymmErrors* errorbandStat2=0;
      TLegend *leg  = new TLegend(); 
      TLegend *legBand = bands ? new TLegend() : 0;
      if(verbose>3) std::cout << "debug E" << std::endl;
      for(unsigned int nTheory=0; nTheory<ratio_.size(); ++nTheory){
	if(verbose>3) std::cout << "  -" << nTheory << std::endl;
	ratio_[nTheory]->SetStats(kFALSE);
	ratio_[nTheory]->SetTitle("");
	ratio_[nTheory]->SetMaximum(ratioMax);
	ratio_[nTheory]->SetMinimum(ratioMin);
	unsigned int color = histDenominatorTheoryOrdered_[nTheory]->GetMarkerColor();
	if(color==kBlack) color=theoryColor((TString)histDenominatorTheoryOrdered_[nTheory]->GetName());
	ratio_[nTheory]->SetLineColor(color);
	double linewidth=histDenominatorTheoryOrdered_[nTheory]->GetLineWidth();
	linewidth=2;
	ratio_[nTheory]->SetLineWidth(linewidth);
	ratio_[nTheory]->SetMarkerColor(color);
	ratio_[nTheory]->SetFillColor(color-4);
	ratio_[nTheory]->SetFillStyle(0);
	ratio_[nTheory]->SetMarkerSize(0.2);
	unsigned int style =theoryStyle((TString)histDenominatorTheoryOrdered_[nTheory]->GetName());
	if(style==42) style=histDenominatorTheoryOrdered_[nTheory]->GetLineStyle();
	//std::cout << histDenominatorTheoryOrdered_[nTheory]->GetName() << ": " << style << std::endl;
	ratio_[nTheory]->SetLineStyle(style);
	// configure axis of ratio_[nTheory] plot
	ratio_[nTheory]->GetXaxis()->SetTitleSize(histDenominatorTheoryOrdered_[nTheory]->GetXaxis()->GetTitleSize()*scaleFactor*1.3);
	ratio_[nTheory]->GetXaxis()->SetTitleOffset(histDenominatorTheoryOrdered_[nTheory]->GetXaxis()->GetTitleOffset()*0.9);
	ratio_[nTheory]->GetXaxis()->SetLabelSize(histDenominatorTheoryOrdered_[nTheory]->GetXaxis()->GetLabelSize()*scaleFactor*1.4);
	ratio_[nTheory]->GetXaxis()->SetTitle(histDenominatorTheoryOrdered_[nTheory]->GetXaxis()->GetTitle());
	ratio_[nTheory]->GetXaxis()->SetNdivisions(histDenominatorTheoryOrdered_[nTheory]->GetNdivisions());
	ratio_[nTheory]->GetYaxis()->CenterTitle();
	ratio_[nTheory]->GetYaxis()->SetTitle("#frac{Theory}{Data}");
	ratio_[nTheory]->GetYaxis()->SetTitleSize(histDenominatorTheoryOrdered_[nTheory]->GetYaxis()->GetTitleSize()*scaleFactor);
	ratio_[nTheory]->GetYaxis()->SetTitleOffset(histDenominatorTheoryOrdered_[nTheory]->GetYaxis()->GetTitleOffset()/scaleFactor);
	ratio_[nTheory]->GetYaxis()->SetLabelSize(histDenominatorTheoryOrdered_[nTheory]->GetYaxis()->GetLabelSize()*scaleFactor);
	ratio_[nTheory]->GetYaxis()->SetLabelOffset(histDenominatorTheoryOrdered_[nTheory]->GetYaxis()->GetLabelOffset()*3.3);
	ratio_[nTheory]->GetYaxis()->SetTickLength(0.03);
	ratio_[nTheory]->GetYaxis()->SetNdivisions(ndiv);
	ratio_[nTheory]->GetXaxis()->SetRange(histNumeratorData->GetXaxis()->GetFirst(), histNumeratorData->GetXaxis()->GetLast());
	// first plot
	if(nTheory==0){
	  if(verbose>3) std::cout << "debug F" << std::endl;
	  setXAxisRange(ratio_[nTheory], (TString)ratio_[nTheory]->GetName());
	  if(rangeMin!=-1.&&rangeMax!=-1.){
	    if(verbose>0) std::cout << "manual range chosen: " << ratio_[nTheory]->FindBin(rangeMin) << " .. " << ratio_[nTheory]->FindBin(rangeMax) << std::endl;
	    ratio_[nTheory]->GetXaxis()->SetRange(ratio_[nTheory]->FindBin(rangeMin),ratio_[nTheory]->FindBin(rangeMax));
	  }
	}
	if(verbose>3) std::cout << "debug G" << std::endl;
	// MC@NLO errorbands
	if(((TString)ratio_[nTheory]->GetName()).Contains("errorBand")){
	  int color2=kGray;
	  if(((TString)ratio_[nTheory]->GetName()).Contains("Dn"))color2=10;
	  ratio_[nTheory]->SetLineColor(color2);
	  ratio_[nTheory]->SetMarkerColor(color2);
	  ratio_[nTheory]->SetFillColor(color2);
	  ratio_[nTheory]->SetFillStyle(1001);
	  ratio_[nTheory]->SetLineWidth(1);
	  ratio_[nTheory]->SetLineStyle(1);
	}
	TString titleX=histNumeratorData->GetXaxis()->GetTitle();
	if(titleX!="") ratio_[nTheory]->GetXaxis()->SetTitle(titleX);
	titleX=(TString)(histNumeratorData->GetXaxis()->GetTitle());
	if(titleX!="") ratio_[nTheory]->GetXaxis()->SetTitle(titleX);
	ratio_[nTheory]->GetXaxis()->SetNoExponent(true);
	// delete axis of initial plot
	histDenominatorTheoryOrdered_[nTheory]->GetXaxis()->SetLabelSize(0);
	histDenominatorTheoryOrdered_[nTheory]->GetXaxis()->SetTitleSize(0);
	// draw ratio_[nTheory] plot
	if(nTheory==0){
          if(verbose>3) std::cout << "debug G2" << std::endl;
	  // line at 1
	  TH1F* one=(TH1F*)ratio_[nTheory]->Clone();
	  one->Divide(one);
	  one->SetLineWidth(1);
	  one->SetLineStyle(1);
	  one->SetLineColor(kBlack);
	  one->SetFillStyle(0);
	  // FIXME: Printout for synchronization
	  //std::cout << "ratio->GetXaxis()->GetTitleSize  : " << one->GetXaxis()->GetTitleSize()   << std::endl;
	  //std::cout << "ratio->GetXaxis()->GetTitleOffset: " << one->GetXaxis()->GetTitleOffset() << std::endl;
	  //std::cout << "ratio->GetXaxis()->GetLabelSize  : " << one->GetXaxis()->GetLabelSize()   << std::endl;
	  //std::cout << "ratio->GetXaxis()->GetLabelOffset: " << one->GetXaxis()->GetLabelOffset() << std::endl;
	  //std::cout << "ratio->GetXaxis()->GetNdivisions : " << one->GetXaxis()->GetNdivisions()  << std::endl;
	  //std::cout << "ratio->GetXaxis()->GetTickLength : " << one->GetXaxis()->GetTickLength()  << std::endl;
	  //std::cout << "ratio->GetYaxis()->GetTitleSize  : " << one->GetYaxis()->GetTitleSize()   << std::endl;
	  //std::cout << "ratio->GetYaxis()->GetTitleOffset: " << one->GetYaxis()->GetTitleOffset() << std::endl;
	  //std::cout << "ratio->GetYaxis()->GetLabelSize  : " << one->GetYaxis()->GetLabelSize()   << std::endl;
	  //std::cout << "ratio->GetYaxis()->GetLabelOffset: " << one->GetYaxis()->GetLabelOffset() << std::endl;
	  //std::cout << "ratio->GetYaxis()->GetTickLength : " << one->GetYaxis()->GetTickLength()  << std::endl;
	  //std::cout << "ratio->GetYaxis()->GetNdivisions : " << one->GetYaxis()->GetNdivisions()  << std::endl;
	  TString add= nTheory==0 ? "" : " same";
	  one2=(TH1F*)one->Clone();
	  one->DrawClone("hist"+add);
	  // line at upper border
	  TH1F* up=(TH1F*)one->Clone("up");
	  up->SetLineWidth(1);
	  up->SetLineColor(kBlack);
	  up->Scale(ratioMax);
	  up->DrawClone("hist same");
	  // errorband for data total error around 1
          if(verbose>3) std::cout << "debug G3" << std::endl;
	  TGraphAsymmErrors* errorband=histStatData ? (TGraphAsymmErrors*)(histStatData->Clone("errorband")) : new TGraphAsymmErrors(one->GetNbinsX()+1);
	  for(int bin=1; bin<=one->GetNbinsX(); ++bin){
	    //int bin2=histStatData->GetX()[bin]==0. ? bin : bin-1; 
	    double totY=histNumeratorData->GetBinContent(bin);
	    double totX=histStatData ? histStatData->GetX()[bin] : histNumeratorData->GetBinCenter(bin);
	    double totErrY=histNumeratorData->GetBinError(bin);
	    errorband->SetPoint(bin, totX, 1);
	    double errUp=bands ? 0.5*histNumeratorData->GetBinWidth(bin) : 0;
	    double errDn=bands ? 0.5*histNumeratorData->GetBinWidth(bin) : 0;
	    errorband->SetPointError(bin, errDn, errUp,totErrY/totY, totErrY/totY);
	    //errorband->SetBinContent(bin, 1);
	    //errorband->SetBinError  (bin, histNumeratorData->GetBinError(bin)/histNumeratorData->GetBinContent(bin));	  
	    //std::cout << bin << ": (x,dy/y)=(" << totX << "," << totErrY/totY << "), (x,y)=(" << totX << "," << totY << "+/-" << totErrY << ")" << std::endl;
	  }
          if(verbose>3) std::cout << "debug G4" << std::endl;
	  int errcolor=kBlack;
	  errorband->SetLineWidth(3.5);
	  errorband->SetLineColor(errcolor);
	  errorband2=(TGraphAsymmErrors*)errorband->Clone();
	  // errorband for data stat error around 1
	  if(histStatData){
	    if(verbose>3) std::cout << "debug G5" << std::endl;
	    TGraphAsymmErrors* errorbandStat=(TGraphAsymmErrors*)(histStatData->Clone("errorbandStat"));
	    for(int bin=1; bin<=one->GetNbinsX(); ++bin){
	      double statY=histStatData->GetY()[bin];
	      double statX=histStatData->GetX()[bin];
	      double statErrY=histStatData->GetEYlow()[bin];
	      errorbandStat->SetPoint(bin, statX, 1);
	      double errUp=bands ? 0.5*histNumeratorData->GetBinWidth(bin) : 0;
	      double errDn=bands ? 0.5*histNumeratorData->GetBinWidth(bin) : 0;
	      errorbandStat->SetPointError(bin, errUp, errDn, statErrY/statY, statErrY/statY);
	      //std::cout << bin << ": (x,dy/y)=(" << statX << "," << statErrY/statY << "), (x,y)=(" << statX << "," << statY << "+/-" << statErrY << ")" << std::endl;
	    }
	    errorbandStat->SetLineWidth(2.);
	    errorbandStat->SetLineColor(errcolor);
	    errorbandStat2=(TGraphAsymmErrors*)errorbandStat->Clone();
	  }
	  // adapt style for errorbands if they are drawn as bands
	  if(bands){
	    // style
	    int totColor =kOrange-4;//kOrange-4;
	    int statColor=kGray+1;//18;
	    int Errstyle=1001;
	    // total error
	    errorband2->SetLineStyle(1);
	    errorband2->SetLineColor(totColor);
	    errorband2->SetLineWidth(1.);
	    errorband2->SetFillStyle(Errstyle);
	    errorband2->SetFillColor(totColor);
	    errorband2->SetMarkerColor(totColor);
	    errorband2->SetMarkerSize(0.01);
	    // stat error
	    if(histStatData){
	      errorbandStat2->SetLineStyle(1);
	      errorbandStat2->SetLineColor(statColor);
	      errorbandStat2->SetLineWidth(1.);
	      errorbandStat2->SetFillStyle(Errstyle);
	      errorbandStat2->SetFillColor(statColor);
	      errorbandStat2->SetMarkerColor(statColor);
	      errorbandStat2->SetMarkerSize(0.01);
	    }
	    // legend
	    legBand->SetX1NDC(0.22);
	    legBand->SetY1NDC(0.97);
	    legBand->SetX2NDC(0.46);
	    legBand->SetY2NDC(0.77);
	    legBand->SetFillStyle(1001);
	    //legBand->SetFillStyle(4050);
	    legBand->SetFillColor(10);
	    legBand->SetBorderSize(0);
	    legBand->SetTextSize(0.1);
	    legBand->SetTextAlign(12);
	    legBand->AddEntry(errorband2, "Stat. #oplus Syst.", "f"); 	    
	    if(errorbandStat2) legBand->AddEntry(errorbandStat2, "Stat.", "f");
	  }
          if(verbose>3) std::cout << "debug G6" << std::endl;
	  leg->SetX1NDC(0.22);
	  leg->SetY1NDC(0.79);
	  leg->SetX2NDC(0.5 );
	  leg->SetY2NDC(0.99);
	  leg ->SetFillStyle(1001);
	  leg ->SetFillColor(10);
	  leg ->SetBorderSize(0);
	  leg ->SetTextSize(0.1);
	  leg ->SetTextAlign(12);
	  leg ->AddEntry(errorband, "data stat+sys error", "F");
	  if(verbose>3) std::cout << "debug G7" << std::endl;
	  if(bands){
	    errorband2    ->Draw("2 same");
	    errorbandStat2->Draw("2 same");
	    one->DrawClone("hist same");
	  }
	}
	ratio_[nTheory]->DrawClone("hist same");
	if(((TString)ratio_[nTheory]->GetName()).Contains("errorBandDn")) one2->Draw("hist same");
	if(nTheory==ratio_.size()-1){
	  gPad->RedrawAxis("g");
	  gPad->RedrawAxis();
	  if(!bands) errorband2->Draw("p z same");
	  gStyle->SetEndErrorSize(8);
	  if(!bands&&errorbandStat2) errorbandStat2->Draw("p e same");
	  if(bands&&legBand) legBand->Draw("same");
	}
      }
      if(verbose>3) std::cout << "debug H" << std::endl;
      if(addXBinGrid){
	for(int bin=2; bin<=histNumeratorData->GetNbinsX(); ++bin){
	  double xBorder=histNumeratorData->GetBinLowEdge(bin);
	  drawLine(xBorder, ratioMin, xBorder, ratioMax, kBlack, 1, 3);
	} 
      }
      rPad->SetTopMargin(0.0);
      rPad->SetBottomMargin(ratioSize);
      rPad->SetRightMargin(right);
      //rPad->SetLeftMargin(left);
      //leg->DrawClone("same");
      //rPad->Print("./"+(TString)(histNumeratorData->GetName())+".png");
      gPad->cd();
      rPad->Draw();
      if(verbose>3) std::cout << "debug H" << std::endl;
      return canv;
    }
    else{
      if(verbose>0){ 
	std::cout << "calling drawFinalResultRatio for ";
	if(histNumeratorData) std::cout << histNumeratorData->GetName();
	else  std::cout << " !WARNING! missing data curve";
	std::cout << " and " << histDenominatorTheory_.size() << " theory curves" << std::endl;
      }
    }
    return (TCanvas*)0;
  }

  TH1F* killEmptyBins(TH1F* histo, int verbose=0){
    // this function removes all bins with width==0 from "histo"
    // verbose: level of text output
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    
    // collect nonemty bins
    std::vector<double> binEdges_;
    std::vector<double> binValues_;
    std::vector<double> blockList_;
    // loop original histogram
    if(verbose>0) std::cout << "scanning " << histo->GetName() << " for technical bins with width 0" << std::endl;
    for(int bin =1; bin<=histo->GetNbinsX(); ++bin){
       if(verbose>1) std::cout << "bin #" << bin << " ( lowedge " << histo->GetBinLowEdge(bin) << ")" << std::endl;
      if(histo->GetBinWidth(bin)!=0.){
	if(verbose>1) std::cout << "nonzero width!"<< std::endl;
	binValues_.push_back(histo->GetBinContent(bin));
	if(binEdges_.size()==0||binEdges_[binEdges_.size()-1]!=histo->GetBinLowEdge(bin)){
	  if(verbose>1) std::cout << "add "<< histo->GetBinLowEdge(bin) << std::endl;
	  binEdges_.push_back(histo->GetBinLowEdge(bin));
	}
	if(binEdges_[binEdges_.size()-1]!=histo->GetBinLowEdge(bin+1)){
	  binEdges_.push_back(histo->GetBinLowEdge(bin+1));
	  if(verbose>1) std::cout << "add "<< histo->GetBinLowEdge(bin+1) << std::endl;
	}
      }
      else{
	blockList_.push_back(bin);
	if(verbose>1) std::cout << "blocked!" << std::endl;
      }
    }
    //if(binEdges_[binEdges_.size()-1]=!histo->GetBinLowEdge(histo->GetNbinsX()+1)) binEdges_.push_back(histo->GetBinLowEdge(histo->GetNbinsX()+1));
    if(verbose>0) std::cout << "scanning 2" << std::endl;
    // create histo without empty bins
    int Nbins=binValues_.size();
    double bins[binEdges_.size()];
    if(verbose>0) std::cout << Nbins << " bins: ";
    for(unsigned int entry=0; entry<binEdges_.size(); ++entry){
      bins[entry]=binEdges_[entry];
      if(verbose>0) std::cout << bins[entry] << ", ";
    }
    if(verbose>0) std::cout << "end" << std::endl;
    if(verbose>0) std::cout << "scanning 3" << std::endl;
    TString lab="mod"+(TString)histo->GetName();
    TH1F* result = new TH1F(lab, lab, Nbins, bins);
    // fill new histo
    int resultbin=1;
    for(int bin =1; bin<=histo->GetNbinsX(); ++bin){
      // check if bin has to be refilled
      bool fill=true;
      for(unsigned int blocked=0; blocked<blockList_.size(); ++blocked){
	if(blockList_[blocked]==bin) fill=false;
      }
      if(fill){
	result->SetBinContent(resultbin,histo->GetBinContent(bin));
	result->SetBinError  (resultbin,histo->GetBinError  (bin));
	resultbin++;
      }
    }
    // first bin
    result->SetBinContent(0,histo->GetBinContent(0));
    result->SetBinError  (0,histo->GetBinError  (0));
    // last bin
    result->SetBinContent(result->GetNbinsX()+1,histo->GetBinContent(histo->GetNbinsX()+1));
    result->SetBinError  (result->GetNbinsX()+1,histo->GetBinError  (histo->GetNbinsX()+1));
    // print out final plot
    if(verbose>1){
      for(int bin =0; bin<=result->GetNbinsX()+1; ++bin){
	std::cout << "bin #" << bin << ": " <<  result->GetBinContent(bin) << " +- " << result->GetBinError(bin);
	std::cout << "(" << result->GetBinLowEdge(bin) << ")" << std::endl;
      }
    }
    // return cleaned plot
    return result;
  }

  double getInclusiveXSec(TH1* hist, int verbose=0)
  {
    // this function integrates a given
    // histogramm 'hist' by summing up the product
    // of binwidth and bincontent for all bins.
    // NOTE: binwidth==0 is treated as binwidth==1
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    double value=0;
    if(verbose>0) std::cout << "histo: " << hist->GetTitle() << std::endl;
    for(int bin=0; bin<=hist->GetNbinsX()+1; ++bin){

      double binWidth=hist->GetBinWidth(bin);
      if(binWidth==0.||bin==0||bin==hist->GetNbinsX()+1) binWidth=1.;
      double thisBin=binWidth*hist->GetBinContent(bin);
      value+=thisBin;
      if(verbose>1){ 
	std::cout << "xSec bin" << bin << ": " << thisBin << std::endl;
	std::cout << "( with lower edge" << hist->GetBinLowEdge(bin) << ")" << std::endl;
      }
    }
    if(verbose>0) std::cout << "integrated xSec: " << value << std::endl;
    return value;
  }


  void whipEmptyBinsAway(TGraphAsymmErrors* hist, int verbose=0)
  {
    // this function shifts all unset points
    // (at x,y=0,0 with error 0,0)
    // out of sight to 0, -1000
    // modified quantities: none
    // used functions: none
    // used enumerators: none
    if(verbose>1) std::cout << "checking " << hist->GetName() << " for empty bins" << std::endl;
    // loop points
    for(int bin=0; bin<=hist->GetN()-1; ++bin){
      // search for unset points
      if( hist->GetX()[bin]==0 &&
	  hist->GetY()[bin]==0 &&
	  hist->GetErrorXhigh(bin)==0 &&
	  hist->GetErrorXlow (bin)==0 &&
	  hist->GetErrorYhigh(bin)==0 &&
	  hist->GetErrorYlow (bin)==0 ){
	hist->SetPoint(bin, 0, -1000);
	if(verbose>1) std::cout << "bin #" << bin << " set to 0,-1000" << std::endl;
      }
    }
  }

  TH1F* getTheoryPrediction(TString plotName="", TString fileName="")
  {

    // check if input is valid
    if(plotName==""||fileName==""){
      std::cout << "ERROR: no plotName or filename chosen in function getTheoryPrediction" << std::endl;
      exit(0);
    }
    // load file
    TFile* file = new (TFile)(fileName);
    // check if file is valid
    if(file->IsZombie()){
      std::cout << "ERROR: input file " << fileName << " in function getTheoryPrediction is zombie" << std::endl;
      exit(0);
    }
    // get plot
    TH1* targetPlot;
    file->GetObject(plotName, targetPlot);
    // check if plot exists
    if(!targetPlot){
      std::cout << "ERROR in function getTheoryPrediction: plot " << plotName;
      std::cout << " does not exist in file " << fileName << std::endl;
      exit(0);
    }
    gROOT->cd();
    TH1F* result = (TH1F*)(targetPlot->Clone());
    result->Sumw2();
    file->Close();
    delete file;
    return result;
  }

  double signum(double value){
    // this function is a
    // simple signum function
    // modified quantities: none
    // used functions: none
    // used enumerators: none

    if(value>0) return  1.0;
    if(value<0) return -1.0;
    return 0;
  }

  void closeStdTopAnalysisFiles(std::map<unsigned int, TFile*> files_)
    {
      // this function closes and deletes the files saved in "files_"
      // modified quantities: files_
      // used functions: NONE
      // used enumerators: NONE

      // close and delete input files
      for(std::map<unsigned int, TFile*>::const_iterator file=files_.begin(); file!=files_.end(); ++file){
	if(file->second){
	  file->second->Close();	  
	  delete file->second;
	}
      }
    }

  // ===========================================================================================
  // ===========================================================================================

  void makeUncertaintyBands(std::map< TString, std::map <unsigned int, TH1F*> >& histo_,
			    std::map< TString, TH1F* >& histoErrorBand_,
			    std::vector<TString>& plotList_,
			    unsigned int& Nplots)
  {
    // this function generates set of errorbands ("histoErrorBand_") for all samples 
    // and histogrammes in "histo_", using the "Nplots" plots listed in "plotList_" 
    // taking into account luminosity and cross section error
    // 
    // modified quantities: histoErrorBand_
    // used functions:      none
    // used enumerators:    samples

    for(unsigned int plot=0; plot<Nplots; plot++){
      
      TString plotName = plotList_[plot];
      
      // Initialize and reset histograms
      
      TH1F* histoSumRef       = (TH1F*)histo_[plotName][kSig]->Clone();
      TH1F* histoSumTTbarOnly = (TH1F*)histo_[plotName][kSig]->Clone();
      
      histoSumRef       -> Reset("ICESM");
      histoSumTTbarOnly -> Reset("ICESM");
      
      // Integral over all samples before accessing the differences
      
      for(unsigned int sample=kSig; sample<kData; ++sample){

	if((plot<Nplots)&&(plotExists(histo_, plotName, sample))){

	  histoSumRef -> Add(histo_[plotName][sample]);

	  if (sample == kSig || sample == kBkg) 
	    histoSumTTbarOnly -> Add(histo_[plotName][sample]);
	}
      }
      
      // Compare summed histograms, symmetrize deviations and store relative deviation to error histogram after adding constant contributions
      
      double relXSecError = ttbarCrossSectionError/ttbarCrossSection;
      
      for (int bin = 0; bin < histoSumRef->GetNbinsX(); bin++){

	double xSecError       = (histoSumTTbarOnly->GetBinContent(bin+1))*relXSecError;
	double luminosityError = (histoSumRef->GetBinContent(bin+1))*globalLumiUncertainty;
	
	double totalError      = sqrt(xSecError*xSecError+luminosityError*luminosityError);
	
	histoSumRef->SetBinError(bin+1,totalError);
      }
      
      histoErrorBand_[plotName] = (TH1F*)histoSumRef->Clone();
    }
  }


  void makeTheoryUncertaintyBands(std::map< TString, std::map <unsigned int, TH1F*> >& histo_,
				  std::map< TString, TH1F* >& histoErrorBand_,
				  std::vector<TString>& plotList_, std::vector<TString>& plotListEl_, std::vector<TString>& plotListMu_,
				  int& Nplots, unsigned int& N1Dplots,
				  TString& inputFolder, TString& dataFileMu, TString& dataFileEl,
				  std::vector<TString> vecRedundantPartOfNameInData,
				  double& luminosityMu, double& luminosityEl, TString type ="model", TString addSel=""
				  )
  {
    // this function generates set of errorbands ("histoErrorBand_") for all samples 
    // and histogrammes in "histo_", using the "Nplots" plots listed in "plotList_" 
    // taking into account ttbar modeling uncertainties specified with "type"
    // type=="model" -> Q2 & matching scale
    // type=="JE" -> JES & JER    
    

    // modified quantities: histoErrorBand_
    // used functions:      getStdTopAnalysisFiles
    // used enumerators:    samples, systematicVariation, 
    
    // internal configurations
    int verbose =0; // printout within makeTheoryUncertaintyBands
    int verbose2=0; // argument passed to called functions
    unsigned int sysNoBG =42;
    unsigned int sysNoAll=42*42;
    TString testvar="";
    //TString testvar="compositedKinematicsKinFit/mHbb";
    //TString testvar="analyzeTopRecoKinematicsKinFitProbSel/topPt";
    
    // debug output
    if(verbose>0) std::cout << "executing function makeTheoryUncertaintyBands" << std::endl;

    // A collect relevant systematics
    // WARNING: systematics are expected to be listed as up/down pairs in this order!!!
    if(verbose>0) std::cout << std::endl << "A collect relevant systematics" << std::endl;
    std::vector<int> RelevantSys_;
    int sysListJE   [ ] = { sysJESUp, sysJESDown, sysJERUp, sysJERDown}; // JE related uncertainties
    int sysListModel[ ] = { sysTopMatchUp, sysTopMatchDown, sysTopScaleUp, sysTopScaleDown, sysTopMassUp, sysTopMassDown, sysHadUp, sysHadDown}; // ttbar modeling
    int sysListMain [ ] = { sysTopMatchUp, sysTopMatchDown, sysTopScaleUp, sysTopScaleDown, sysHadUp, sysHadDown, sysJESUp, sysJESDown, sysJERUp, sysJERDown, sysBtagSFUp, sysBtagSFDown, sysBtagSFShapeUpPt65, sysBtagSFShapeDownPt65, sysBtagSFShapeUpEta0p7, sysBtagSFShapeDownEta0p7 }; // ttbar modeling
    //int sysListPaper[ ] = {sysBtagSFUp, sysBtagSFDown, sysBtagSFShapeDownPt65, sysBtagSFShapeUpEta0p7};
    int sysListPaper[ ] = { sysTopMatchUp, sysTopMatchDown, sysTopScaleUp, sysTopScaleDown, sysHadUp, sysHadDown, sysTopMassUp, sysTopMassDown, sysPDFUp, sysPDFDown}; // ttbar main modeling unc.
    int sysListPHD[ ] = {  sysTopMatchUp, sysTopMatchDown, sysTopScaleUp, sysTopScaleDown, sysHadUp, sysHadDown, sysJESUp, sysJESDown, sysJERUp, sysJERDown, sysBtagSFUp, sysBtagSFDown, sysBtagSFShapeUpPt65, sysBtagSFShapeDownPt65, sysBtagSFShapeUpEta0p7, sysBtagSFShapeDownEta0p7, sysMisTagSFUp, sysMisTagSFDown}; // separate list used for martins thesis
    int sysListAll  [ ] = { sysTopMatchUp, sysTopMatchDown, sysTopScaleUp, sysTopScaleDown, sysTopMassUp, sysTopMassDown, sysJESUp, sysJESDown ,sysJERUp, sysJERDown, sysPUUp, sysPUDown, sysLepEffSFNormUp, sysLepEffSFNormDown, sysLepEffSFShapeUpEta, sysLepEffSFShapeDownEta, sysLepEffSFShapeUpPt, sysLepEffSFShapeDownPt, sysBtagSFUp, sysBtagSFDown, sysBtagSFShapeUpPt65, sysBtagSFShapeDownPt65, sysBtagSFShapeUpEta0p7, sysBtagSFShapeDownEta0p7, sysMisTagSFUp, sysMisTagSFDown, sysBRUp, sysBRDown, sysPDFUp, sysPDFDown, sysHadUp, sysHadDown};
    
    if(     type =="model") RelevantSys_.insert(RelevantSys_.begin(), sysListModel, sysListModel+ sizeof(sysListModel)/sizeof(int));
    else if(type =="JE"   ) RelevantSys_.insert(RelevantSys_.begin(), sysListJE   , sysListJE   + sizeof(sysListJE   )/sizeof(int));
    else if(type =="all"  ) RelevantSys_.insert(RelevantSys_.begin(), sysListAll  , sysListAll  + sizeof(sysListAll  )/sizeof(int));
    else if(type =="main" ) RelevantSys_.insert(RelevantSys_.begin(), sysListMain , sysListMain + sizeof(sysListMain )/sizeof(int));
    else if(type =="paper") RelevantSys_.insert(RelevantSys_.begin(), sysListPaper, sysListPaper+ sizeof(sysListPaper)/sizeof(int));
    else if(type =="phd"  ) RelevantSys_.insert(RelevantSys_.begin(), sysListPHD  , sysListPHD  + sizeof(sysListPHD)  /sizeof(int));
    else{ 
      std::cout << "ERROR in function  " << std::endl;
      exit(0);
    }
    if(verbose>1){
      std::cout << "considered systematics: makeTheoryUncertaintyBands, provided input argument type=" << type << " is invalid" << std::endl;
      for(unsigned int i=0; i<RelevantSys_.size(); ++i){
	std::cout << sysLabel(RelevantSys_[i]) << std::endl;
      }
    }
    // define plot container: one histo for each sample and systematic variation
    std::map< TString, std::map <unsigned int, TH1F*> > histoSys_;

    // B Get ttbar sysNo Reference plot that is drawn as central MC prediction
    if(verbose>0) std::cout << std::endl << "B Get ttbar sysNo Reference plot" << std::endl;
    // loop plots
    for(unsigned int plot=0; plot<N1Dplots; plot++){
      // name of the plot
      TString plotName = plotList_[plot];  
      // debug output
      if(verbose>1) std::cout << "plot #" << plot+1 << "/" << N1Dplots << ": " << plotName;
      // consider only plots available for SG+BG (+FIXME: also Z as there is a bug with the binning in one case...)
      if(histo_[plotName].count(kSig)>0&&histo_[plotName].count(kBkg)>0&&histo_[plotName].count(kZjets)>0){
	if(verbose>1)  std::cout << " (ok) ";
	// get ttbar (SG+BG) sysNo
	histoSys_[plotName][sysNo]=     (TH1F*)histo_[plotName][kSig]->Clone();
	histoSys_[plotName][sysNo]->Add((TH1F*)histo_[plotName][kBkg]->Clone());
	// get non-ttbar sysNo
	histoSys_[plotName][sysNoBG ]=(TH1F*)histoSys_[plotName][sysNo]->Clone();
	histoSys_[plotName][sysNoBG ]->Reset("ICESM");
	for(unsigned int sample=kZjets; sample<kData; ++sample){
	  histoSys_[plotName][sysNoBG ]->Add((TH1F*)histo_[plotName][sample]->Clone());
	}
	// get combined sysNo
	histoSys_[plotName][sysNoAll]=     (TH1F*)histoSys_[plotName][sysNo  ]->Clone();
	histoSys_[plotName][sysNoAll]->Add((TH1F*)histoSys_[plotName][sysNoBG]->Clone());
      }
      else if(verbose>1) std::cout << " (skipped) ";	
      if(verbose>1) std::cout << std::endl;
    }
    
    // C Get ttbar (SG+BG) systematically shifted plots
    if(verbose>0) std::cout << std::endl << "C Get ttbar systematic shifted plots (note: takes centuries)" << std::endl;
    // loop systematics
    for(unsigned int sys=0; sys<RelevantSys_.size(); ++sys){
      if(verbose==1) std::cout << ".";
      int sysNow=RelevantSys_[sys];
      TString inputFolder2=inputFolder;
      TString sysInputFolderExtension="";
      TString addSelData="";
      std::vector<TString> plotListElMod_=plotListEl_;
      std::vector<TString> plotListMuMod_=plotListMu_;
      // special treatment for hadronization uncertainty
      TString ttbarMC="Madgraph";
      int sysNow2=sysNow;
      if(     sysNow==sysHadUp  ){ sysNow2=sysNo; ttbarMC="Mcatnlo";}
      else if(sysNow==sysHadDown){ sysNow2=sysNo; ttbarMC="Powheg" ;}
      // for systematic variations done in central signal sample (e.g. PU, effSF, btagSF):
      // use MC rootfiles from subfolder /Prob and ignore addSel plotname extension
      else if(addSel=="ProbSel"&&((sysNow2>=sysPUUp&&sysNow2<=sysPUDown)||(sysNow2>=sysLepEffSFNormUp&&sysNow2<=sysMisTagSFDown))){
	inputFolder2+="/Prob";
	addSelData=addSel;
	switch (sysNow2)
	  {
	  case sysPUUp                     : sysInputFolderExtension="PUup";                  break;
	  case sysPUDown                   : sysInputFolderExtension="PUdown";                break;
	  case sysLepEffSFNormUp           : sysInputFolderExtension="EffSFNormUpSys";        break;
	  case sysLepEffSFNormDown         : sysInputFolderExtension="EffSFNormDownSys";      break;
	  case sysLepEffSFShapeUpEta       : sysInputFolderExtension="EffSFShapeUpEta";       break;
	  case sysLepEffSFShapeDownEta     : sysInputFolderExtension="EffSFShapeDownEta";     break;
	  case sysLepEffSFShapeUpPt        : sysInputFolderExtension="EffSFShapeUpPt";        break;
	  case sysLepEffSFShapeDownPt      : sysInputFolderExtension="EffSFShapeDownPt";      break;
	  case sysBtagSFUp                 : sysInputFolderExtension="BtagSFup";              break;
	  case sysBtagSFDown               : sysInputFolderExtension="BtagSFdown";            break;
	  case sysBtagSFShapeUpPt65        : sysInputFolderExtension="BTagSFShapeUpPt65";     break;
	  case sysBtagSFShapeDownPt65      : sysInputFolderExtension="BTagSFShapeDownPt65";   break;
	  case sysBtagSFShapeUpEta0p7      : sysInputFolderExtension="BTagSFShapeUpEta0p7";   break;
	  case sysBtagSFShapeDownEta0p7    : sysInputFolderExtension="BTagSFShapeDownEta0p7"; break;
	  case sysMisTagSFUp               : sysInputFolderExtension="MisTagSFup";            break;
	  case sysMisTagSFDown             : sysInputFolderExtension="MisTagSFdown";          break;
	  default: break;
	  }
	for(unsigned int i=0; i<plotListMuMod_.size(); ++i){
	  TString plotNameMu=plotListMu_[i];
	  if(plotNameMu.Contains(addSel)){
	    plotNameMu.ReplaceAll(addSel, sysInputFolderExtension);
	    if(plotNameMu.Contains("composited")) plotNameMu.ReplaceAll(sysInputFolderExtension, "KinFit"+sysInputFolderExtension);
	    plotListMuMod_[i]=plotNameMu;
	  }
	}
	for(unsigned int i=0; i<plotListElMod_.size(); ++i){ 
	  TString plotNameEl=plotListEl_[i];
	  if(plotNameEl.Contains(addSel)){
	    plotNameEl.ReplaceAll(addSel, sysInputFolderExtension);
	    if(plotNameEl.Contains("composited")) plotNameEl.ReplaceAll(sysInputFolderExtension, "KinFit"+sysInputFolderExtension);
	    plotListElMod_[i]=plotNameEl;
	  }
	}
      }
      if(verbose>1) std::cout << "- " << sysNow << ":" << std::endl;
      // c1 get files
      if(verbose>1) std::cout << "c1 get rootfiles" << std::endl;
      std::map<unsigned int, TFile*> filesMu_, filesEl_;
      filesMu_ = getStdTopAnalysisFiles(inputFolder2, sysNow2, dataFileMu, "muon"    , ttbarMC);
      filesEl_ = getStdTopAnalysisFiles(inputFolder2, sysNow2, dataFileEl, "electron", ttbarMC);
      // c2 get plots 
      if(verbose>1) std::cout << "c2 get all plots" << std::endl;
      std::map< TString, std::map <unsigned int, TH1F*> > histoEl_ , histoMu_ , histoComb_;
      std::map< TString, std::map <unsigned int, TH2F*> > histo2El_, histo2Mu_, histo2Comb_;
      // store plots NOT into the files
      gROOT->cd();
      // get plots
      getAllPlots(filesEl_, plotListElMod_, histoEl_, histo2El_, N1Dplots, Nplots, verbose2, "electron", &vecRedundantPartOfNameInData, false, addSelData);
      getAllPlots(filesMu_, plotListMuMod_, histoMu_, histo2Mu_, N1Dplots, Nplots, verbose2, "muon"    , &vecRedundantPartOfNameInData, false, addSelData);
      // take care of adapted naming
      for(unsigned int i=0; i<plotListElMod_.size(); ++i){ 
	TString plotNameEl=plotListEl_[i];
	if(plotNameEl.Contains(addSel)&&(plotListElMod_[i]!=plotListEl_[i])){
	  for(unsigned int sample=kSig; sample<kData; ++sample)
	    if(histoEl_.count(plotListElMod_[i])>0&&histoEl_[plotListElMod_[i]].count(sample)>0) histoEl_[plotListEl_[i]][sample]=histoEl_[plotListElMod_[i]][sample];
	}
      }
      for(unsigned int i=0; i<plotListMuMod_.size(); ++i){ 
	TString plotNameMu=plotListMu_[i];
	if(plotNameMu.Contains(addSel)&&(plotListMuMod_[i]!=plotListMu_[i])){
	  for(unsigned int sample=kSig; sample<kData; ++sample)
	    if(histoMu_.count(plotListMuMod_[i])>0&&histoMu_[plotListMuMod_[i]].count(sample)>0) histoMu_[plotListMu_[i]][sample]=histoMu_[plotListMuMod_[i]][sample];
	}
      }
      // close files
      closeStdTopAnalysisFiles(filesMu_);
      closeStdTopAnalysisFiles(filesEl_);
      // c3 lumi scale
      if(verbose>1) std::cout << "c3 lumiscaling" << std::endl;
      scaleByLuminosity(plotListMu_, histoMu_, histo2Mu_, N1Dplots, luminosityMu, verbose2, sysNow2, "muon"    , ttbarMC);
      scaleByLuminosity(plotListEl_, histoEl_, histo2El_, N1Dplots, luminosityEl, verbose2, sysNow2, "electron", ttbarMC);
      // c4 add channels, collect plots in temporary existing map histoComb_
     if(verbose>1) std::cout << "c4 add mu and el channel plots" << std::endl;
      // loop samples
      for(unsigned int sample=kSig; sample<kData; ++sample){
	if(verbose>2) std::cout << sampleLabel(sample, "combined") << std::endl;
	// loop plots
	for(unsigned int plot=0; plot<plotList_.size(); ++plot){
	  if(verbose>2) std::cout << "comb: " << plotList_[plot] << ", mu: " << plotListMu_[plot] << ", el: " << plotListEl_[plot] << std::endl;
	  // c4.1) 1D, existing in both
	  if((plot<N1Dplots)&&(histoMu_.count(plotListMu_[plot])>0)&&(histoMu_[plotListMu_[plot]].count(sample)>0)&&(histoEl_.count(plotListEl_[plot])>0)&&(histoEl_[plotListEl_[plot]].count(sample)>0)){ 
	    if(verbose>2) std::cout << "-> 1D" << std::endl;
	    // add channels
	    histoComb_[plotList_[plot]][sample]=     (TH1F*)(histoMu_[plotListMu_[plot]][sample]->Clone());
	    // take care of lepton pt from different folders
	    if (plotListMu_[plot]=="tightMuonKinematics/pt"){
	      histoComb_[plotList_[plot]][sample]->Add((TH1F*)(histoEl_["tightElectronKinematics/et"][sample]->Clone()));
	    }
	    else histoComb_[plotList_[plot]][sample]->Add((TH1F*)(histoEl_[plotListEl_[plot]][sample]->Clone()));
	  }
	  // c4.2) 2D, existing in both
	  else if((plot>=N1Dplots)&&(histo2Comb_.count(plotList_[plot])>0)&&(histo2Comb_[plotList_[plot]].count(sample)>0)){
	    if(verbose>2) std::cout << "-> 2D" << std::endl;
	    // add them
	    histo2Comb_[plotList_[plot]][sample]=     (TH2F*)(histo2Mu_[plotListMu_[plot]][sample]->Clone());
	    histo2Comb_[plotList_[plot]][sample]->Add((TH2F*)(histo2El_[plotListEl_[plot]][sample]->Clone()));
	  }
	  else{
	    if(verbose>2){
	      std::cout << "WARNING: " << plotList_[plot] << "NOT FOUND" << std::endl;
	      if(!histoMu_.count(plotListMu_[plot])>0) std::cout << "in histoMu" << std::endl;
	      if( histoMu_.count(plotListMu_[plot])>0&&!histoMu_[plotListMu_[plot]].count(sample)>0) std::cout << "in mu sample " << sampleLabel(sample, "muon"    ) << std::endl;
	      if(!histoEl_.count(plotListEl_[plot])>0) std::cout << "in histoEl" << std::endl;
	      if( histoEl_.count(plotListEl_[plot])>0&&!histoEl_[plotListEl_[plot]].count(sample)>0) std::cout << "in el sample " << sampleLabel(sample, "electron") << std::endl;
	    } // end if verbose
	  } // end else (plot not found)
	} // end loop plots
      } // end loop sample
      // c5 rebinning and storing of systematically shifted ttbar SG&BG plots in histoSys_
      if(verbose>1) std::cout << "c5 rebinning and storing of systematically shifted plots" << std::endl;
      // loop plots
      for(unsigned int plot=0; plot<plotList_.size(); ++plot){
	// plot name
	TString plotName = plotList_[plot];
	// debug output
	if(verbose>2) std::cout << "- plot #" << plot+1 << "/" << plotList_.size() << ": " << plotName;
	// existing 1D plots
	if(  histoComb_.count(plotList_[plot])>0 && histoComb_[plotName].count(kSig)>0 && histo_.count(plotList_[plot])>0 && plot<N1Dplots ){
	  // do rebinning (based on binning for sysNo that is drawn as central MC value)
	  int reBinFactor=roundToInt((histoComb_[plotName][kSig]->GetNbinsX())/(histo_[plotName][sysNo]->GetNbinsX()));
	  if(reBinFactor>1){
	    equalReBinTH1(reBinFactor, histoComb_, plotName, kSig);
	    if(histoComb_[plotName].count(kBkg)>0) equalReBinTH1(reBinFactor, histoComb_, plotName, kBkg);
	    else std::cout << "WARNING: ttbar BG NOT found for plot " << plotName << " in makeTheoryUncertaintyBands for variation " << sysLabel(sysNow) << std::endl;
	  }
	  // store Plot
	  histoSys_[plotName][sysNow]=                                            (TH1F*)histoComb_[plotName][kSig]->Clone();
	  if(histoComb_[plotName].count(kBkg)>0) histoSys_[plotName][sysNow]->Add((TH1F*)histoComb_[plotName][kBkg]->Clone());
	  // ensure same normalization of systematically shifted plot as for sysNo case
	  // (N.B.: the idea is to adress only shape uncertainties, as the ttbar component is
	  //        anyhow scaled wrt. the inclusive cross section derived from data)	  
	  histoSys_[plotName][sysNow]->Scale((histoSys_[plotName][sysNo]->Integral(0,histoSys_[plotName][sysNo]->GetNbinsX()+1))/(histoSys_[plotName][sysNow]->Integral(0,histoSys_[plotName][sysNow]->GetNbinsX()+1)));
	}
      } // for plots
    } // end loop systematics
    if(verbose==1) std::cout << std::endl;

    // D Symmetrize and Combine Variations, create error band  
    if(verbose>0) std::cout << std::endl << "D Symmetrize and Combine Variations, create error band" << std::endl;    
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // name of plot
      TString plotName = plotList_[plot];
      bool test= testvar!=""&&plotName.Contains(testvar) ? true : false;
      // debug output
      if(verbose>1) std::cout << "- plot #" << plot+1 << "/" << plotList_.size() << ": " << plotName;
      // process only 1D plots existing for sysNo and for systematically shifted plots
      if( histoSys_.count(plotName)>0 && histo_.count(plotList_[plot])>0 && plot<N1Dplots){
	// clone plot for total MC (SG+BG+nottbarBG) sysNo for errorband
	histoErrorBand_[plotName]=(TH1F*)histoSys_[plotName][sysNoAll]->Clone();
	if(test) std::cout << plotName << std::endl;
	// loop bins
	for (int bin=1; bin<=histoErrorBand_[plotName]->GetNbinsX(); ++bin){
	  // debug output
	  if(verbose>2||test){
	    std::cout << "  - bin #" << bin << "/" << histoErrorBand_[plotName]->GetNbinsX();
	    std::cout << ", <x>=" << histoErrorBand_[plotName]->GetBinCenter(bin) << ":" << std::endl;
	  }
	  double NttbarSysNo=histoSys_[plotName][sysNo  ]->GetBinContent(bin);
	  double NBGSysNo   =histoSys_[plotName][sysNoBG]->GetBinContent(bin);
	  double uncUp=0;
	  double uncDn=0;
	  // loop systematics
	  for(unsigned int sys=0; sys<RelevantSys_.size(); ++sys){
	    int sysNow=RelevantSys_[sys];
	    // debug output
	    if(verbose>2||test) std::cout << "    - " << sysLabel(sysNow);
	    // check if considered systematics exists for this plot 
	    if(histoSys_[plotName].count(sysNow)>0){
	      if(verbose>2||test) std::cout << " (considered)" << std::endl;
	      // process up and down at the same time (during dn is looped)
	      if(verbose>2) std::cout << sys+1 << "%2=" << (sys+1)%2 << std::endl;
	      if(((sys+1)%2)==0){
		if(verbose>2) std::cout << "down systematics" << std::endl;
		int sysUp=RelevantSys_[sys-1];
		int sysDn=RelevantSys_[sys  ];
		double NttbarSysUp=histoSys_[plotName][sysUp]->GetBinContent(bin);
		double NttbarSysDn=histoSys_[plotName][sysDn]->GetBinContent(bin);
		double relDiffUp=std::abs((NttbarSysUp-NttbarSysNo)/(NttbarSysNo+NBGSysNo));
		double relDiffDn=std::abs((NttbarSysDn-NttbarSysNo)/(NttbarSysNo+NBGSysNo));
		// apply SF for mass uncertainty
		if(sysUp==sysTopMassUp&&sysDn==sysTopMassDown){
		  relDiffDn*=SF_TopMassDownUncertainty;
		  relDiffUp*=SF_TopMassUpUncertainty;
		}
		// special treatment for hadronization uncertainty
		if(sysUp==sysHadDown||sysUp==sysHadUp){
		  relDiffUp=std::abs((NttbarSysUp-NttbarSysDn)/(NttbarSysNo+NBGSysNo));
		  relDiffDn=relDiffUp;
		}
		// 100% unc. for empty bins
		if((NttbarSysNo+NBGSysNo)==0.){
		  relDiffDn=1.0;
		  relDiffUp=1.0;
		}
		if(verbose>2||test){
		  if(!(sysUp==sysHadDown||sysUp==sysHadUp)) std::cout << "      (NttbarSysNo=" << NttbarSysNo << ")" << std::endl; 
		  std::cout << "      (NBGSysNo   =" << NBGSysNo    << ")" << std::endl;
		  std::cout << "       NttbarSysUp=" << NttbarSysUp << std::endl; 
		  std::cout << "       NttbarSysDn=" << NttbarSysDn << std::endl; 
		  std::cout << "       relDiffUp  =" << relDiffUp*100 << "%" << std::endl; 
		  std::cout << "       relDiffDn  =" << relDiffDn*100 << "%" << std::endl; 
		}
		// collect ^2 of all uncertainties
		uncUp+=relDiffUp*relDiffUp;
		uncDn+=relDiffDn*relDiffDn;
	      } // end if modulo 2 - processing Dn systematics
	    } // end if sys is successfully loaded
	    else if(verbose>2) std::cout << " (skipped)" << std::endl;
	  } // end loop systematics
	  // symmetrize total up and total down uncertainty
	  double uncSymm=(sqrt(uncUp)+sqrt(uncDn))/2;
	  // Set uncertainty value for errorband
	  if(testvar!=""&&plotList_[plot].Contains(testvar)){
	    std::cout << "  -> tot. symm. unc.=" << uncSymm*100 << "%" << std::endl;
	  }
	  histoErrorBand_[plotName]->SetBinError(bin, uncSymm*histoErrorBand_[plotName]->GetBinContent(bin));
	} // end loop bins
      } // end if plot exists
      else if(verbose>1) std::cout << " (skipped)" << std::endl;
    } // end loop plots    
  }

  // ===========================================================================================
  // ===========================================================================================
  
  void addCanvas(std::vector<TCanvas*> &plotCanvas_){
    // this function adds a new canvas to the vector of canvas "plotCanvas_"
    // modified quantities: plotCanvas_
    // used functions: NONE
    // used enumerators: NONE
    char canvname[10];
    int number =plotCanvas_.size()+1;
    sprintf(canvname,"canv%i",number);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    //canvasStyle(*plotCanvas_[number-1]);
  }

  TH1F* useFittedFunctions(TH1F* input, TString model="", TString plotname="", int verbose=0){
    // this function modifys "input" in a way that wiggly parts 
    // with large uncertainties are replaced by fitted functions
    // modified quantities: input
    // used functions: NONE
    // used enumerators: NONE
    // input: fine binned theory curve
    // model: indicates theory (madgraph, powheg or mcatnalo)
    // verbose: level of output
    int initialIgnoreLevel=gErrorIgnoreLevel;
    // check chosen PS
    TString PS="";
    if(!plotname.Contains("PhaseSpace")) PS+="full";
    else if(plotname.Contains("Hadron")) PS+="had";
    else if(plotname.Contains("Parton")) PS+="part";
    else std::cout << "error in useFittedFunctions, can not get PS from plotname " << plotname << std::endl;
    // if large sample:
    double largeSampleSF=1;
    bool largeSample=false;
    if(model=="madgraph"&&((PS=="full"&&input->GetEntries()>10000000)||(PS!="full"&&input->GetEntries()>3000000))){
      largeSampleSF=25424818/6923750;
      largeSample=true;
    }
    if(verbose==0) gErrorIgnoreLevel=kFatal;
    if(verbose>1) std::cout << plotname << "(" << model << ")" << std::endl;
    // clone histo for output
    TH1F* result=(TH1F*)input->Clone();
    // define fit function A
    double fitLowEdge=0;
    double fitHighEdge=0;
    TString def="";
    double a=0.;
    double b=0.;
    double c=0.;
    double d=0.;
    double e=0.;
    TString addOpt="";
    // define fit function B
    double fitLowEdgeB=0;
    double fitHighEdgeB=0;
    TString defB="";
    double aB=0.;
    double bB=0.;
    double cB=0.;
    double dB=0.;    double eB=0.;
    TString addOptB="";
    // define fit function C
    double fitLowEdgeC=0;
    double fitHighEdgeC=0;
    TString defC="";
    double aC=0.;
    double bC=0.;
    double cC=0.;
    double dC=0.;
    double eC=0.;
    TString addOptC="";
    // configure parameter for all specific plots
    if(model=="mcatnlo"){
      if(plotname.Contains("TTbarM")){
	//tail:
	fitLowEdge=430.;
	fitHighEdge=1600.;
	def="[0]*exp([1]*x)+[2]";
	a=0.11;
	b=-0.0074;
	c=0.00001;
	addOpt="LL";
	// start:
	fitLowEdgeB=345.;
	fitHighEdgeB=400.;
	defB="TMath::GammaDist(x,[0],[1],[2])";
	aB=1.53;
	bB=345.;
	cB=96.2;
	addOptB="LL";
	// very early start:
	fitLowEdgeC=317.;
	fitHighEdgeC=336.;
	defC="[0]*exp([1]*x)+[2]";
	aC=0.00000373757;
	bC=0.0125775;
	cC=-0.000194989;
	addOptC="LL";
      }
      else if(plotname.Contains("TTbarPt")){
	//tail:
	fitLowEdge=50.;
	fitHighEdge=300.;
	def="[0]*exp([1]*x)+[2]*x+[3]*x*x+[4]";
	a= 0.0214186;
	b=-0.0274519;
	c=-0.0000135818;
	d= 0.0000000213453;
	e= 0.00227380;
	addOpt="LL";
      }
     else if(plotname.Contains("TTbarY")){
	fitLowEdge=0.5;
	fitHighEdge=-0.75;
	def="TMath::Exp(x*x*[0]+x*x*x*x*[1]+x*x*x*x*x*x*[2])*[3]";
	a=-0.52;
	b=-0.047;
	c=-0.094;
	d=812149;
	addOpt="LL";
     }
     else if(plotname.Contains("TopY")){
       //tail:
       fitLowEdge=-1.;
       fitHighEdge=1.;
       def="[0]*TMath::Gaus(x,0,[1])";
       a=0.00402075;
       b=1.14871;
       addOpt="LL";
     }
     //else if(plotname.Contains("LepPt")){
     //  fitLowEdge=45.;
     //  fitHighEdge=210.;
     //  def="TMath::Exp(x*[0]+x*x*[1])*[2]";
     //  a=-0.0414389;
     //  b=0.0000334359;
     //  c=0.133893;
     //}
    }
    else if(model=="madgraph"){
      if(plotname.Contains("ttbarMass")){
	//tail:
	fitLowEdge=470.0;
	fitHighEdge=720.0;
	def="[0]*exp([1]*x)+[2]";
	a= (PS=="full" ? 164748.0 : 28165.)*largeSampleSF;
	b= (PS=="full" ? -0.00810181 : -0.00756);
	c= (PS=="full" ? 3.59495 : 2.77)*largeSampleSF;
	addOpt="LL";
        fitLowEdgeB=680.0;
        fitHighEdgeB=920.0;
        defB="[0]*exp([1]*x)+[2]";
        aB= (PS=="full" ? 164748.0 : 28165.)*largeSampleSF;
        bB= (PS=="full" ? -0.00810181 : -0.00756);
        cB= (PS=="full" ? 3.59495 : 2.77)*largeSampleSF;
        addOptB="LL";  
	//end of tail:
	fitLowEdgeC=880.0;
	fitHighEdgeC=1600.0;
	defC="[0]*exp([1]*x)+[2]";
	aC= 54699.0*largeSampleSF;  
	bC=-0.0068;
	cC= 0.979;
	addOptC="LL";
	// start:
	//if(!largeSample){
	//  fitLowEdgeC=345.;
	//  fitHighEdgeC=400.;
	//  defC="[3]*TMath::GammaDist(x,[0],[1],[2])";
	//  aC=(PS=="full" ? 1.50 : 1.56 );
	//  bC=(PS=="full" ? 345. : 344.8);
	//  cC=(PS=="full" ? 79.3 : 80.8 );
	//  dC=(PS=="full" ? 1010990 : 1000000.);
	//  addOptC="LL";
	//}
      }
/*       else if(plotname.Contains("ttbarPt")){ */
/* 	//tail: */
/* 	fitLowEdge=50.; */
/* 	fitHighEdge=300.; */
/* 	def="[0]*exp([1]*x)+[2]*x+[3]*x*x+[4]"; */
/* 	a= 5180.*largeSampleSF; */
/* 	b=-0.03; */
/* 	c=-4.48; */
/* 	d= 0.007; */
/* 	e= 745.5; */
/* 	addOpt="LL"; */
/*       } */
/*       else if(plotname.Contains("topY")){ */
/* 	//tail: */
/* 	if(!largeSample){ */
/* 	  fitLowEdge=-1.; */
/* 	  fitHighEdge=1.; */
/* 	  def="[0]*TMath::Gaus(x,0,[1])"; */
/* 	  a=19845.2; */
/* 	  b=1.09186; */
/* 	  addOpt="LL"; */
/* 	} */
/*       } */
      else if(plotname.Contains("TTbarY")){
	fitLowEdge=-0.5;
	fitHighEdge=0.75;
	def="TMath::Exp(x*x*[0]+x*x*x*x*[1]+x*x*x*x*x*x*[2])*[3]";
	a=-0.52;
	b=-0.047;
	c=-0.094;
	d=812149;
	addOpt="LL";
      }
      else if(plotname.Contains("lepPt")){ 
/* 	// head */
/* 	fitLowEdge = 30.; */
/* 	fitHighEdge= 58.; */
/* 	def="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]"; */
/* 	a=(largeSample ?  5.53215e-02 :  4.26320e-02); */
/* 	b=(largeSample ? -1.36559e-03 : -1.11187e-03); */
/* 	c=(largeSample ?  7.65236e-06 :  6.03412e-06); */
/* 	d=(largeSample ?  6.78192e+04 :  4.84434e+03); */
 	// tail 
 	fitLowEdgeB =50.; 
 	fitHighEdgeB=200.; 
 	defB="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]"; 
 	aB=(largeSample ? -2.03019e-02 : -1.39410e-02); 
 	bB=(largeSample ? -9.19906e-05 : -1.49428e-04); 
 	cB=(largeSample ?  2.63030e-07 :  4.11490e-07); 
 	dB=(largeSample ?  3.18800e+05 :  1.51488e+04); 
      } 
      else if(plotname.Contains("bqPt")){
 	// head
 	fitLowEdge =30.0;
 	fitHighEdge=93.0;
 	def="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]";
 	a= 6.96386e-02;
 	b=-1.01480e-03;
 	c= 3.86642e-06;
 	d= 1.63953e+03;
 	// tail
 	fitLowEdgeB = 92.0;
 	fitHighEdgeB=410.0;
 	defB="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]";
 	aB=-1.87272e-02;
 	bB=-3.43912e-05;
 	cB= 7.16568e-08;
 	dB= 2.67318e+04;
      }
/*       else if(plotname.Contains("bqPt")&&plotname.Contains("Hadron")){ */
/* 	// tail */
/* 	fitLowEdge =100.0; */
/* 	fitHighEdge=400.0; */
/* 	def="TMath::Exp(x*[0]+x*x*[1]+x*x*x*[2])*[3]"; */
/* 	a=-2.08896e-02; */
/* 	b=-2.63966e-05; */
/* 	c= 6.52234e-08; */
/* 	d= 5.32653e+05; */
/*       } */
    }
    else if(model=="powheg"){
      if(plotname.Contains("ttbarMass")){
      fitLowEdge=480.;
      fitHighEdge=1600.;
      def="[0]*exp([1]*x)+[2]";
      a=40565.;
      b=-0.007;
      c=3.1;
      addOpt="LL";
      // start:
      fitLowEdgeB=345.;
      fitHighEdgeB=500.;
      defB="[3]*TMath::GammaDist(x,[0],[1],[2])";
      aB=1.58;
      bB=344.;
      cB=86.6;
      dB=234473.;
      addOptB="LL";
      // very early start:
      fitLowEdgeC=295.;
      fitHighEdgeC=337.;
      defC="[0]*exp([1]*x)+[2]";
      aC=0.000000000530501;
      bC=0.075;
      cC=3.81;
      addOptC="LL";
      }
      else if(plotname.Contains("ttbarPt")){
	//tail:
	fitLowEdge=50.;
	fitHighEdge=300.;
	def="[0]*exp([1]*x)+[2]*x+[3]*x*x+[4]";
	a= 5180.;
	b=-0.03;
	c=-4.48;
	d= 0.007;
	e= 745.5;
	addOpt="LL";
      }
     else if(plotname.Contains("ttbarY")){
	fitLowEdge=-1.5;
	fitHighEdge=1.5;
	def="TMath::Exp(x*x*[0]+x*x*x*x*[1]+x*x*x*x*x*x*[2])*[3]";
	a=-0.598570;
	b=-0.120246;
	c=-0.00770785;
	d=12987.6;
	addOpt="LL";
      }
     //else if(plotname.Contains("lepPt")){
     //  fitLowEdge=45.;
     //  fitHighEdge=210.;
     //  def="TMath::Exp(x*[0]+x*x*[1])*[2]";
     //  a=-0.0414389;
     //  b=0.0000334359;
     //  c=0.133893;
     //}
    }
    // use fit function A
    if(def!=""&&fitLowEdge!=fitHighEdge){
      if(verbose>0) std::cout << "will perform fit" << std::endl;
      // empty bins within fit range for output histo
      bool lastBin=false;
      for(int bin=1; bin<=result->GetNbinsX()+1; ++bin){
	bool thisBin=false;
	double binX=result->GetBinCenter(bin);
	if(binX>=fitLowEdge&&binX<=fitHighEdge){ 
	  result->SetBinContent(bin, 0.);
	  thisBin=true;
	}
	// adapt fit range to binning
	// low
	if(lastBin==false&&thisBin==true) fitLowEdge=result->GetBinLowEdge(bin);
	// high
	if(lastBin==true&&thisBin==false) fitHighEdge=result->GetBinLowEdge(bin);	
	lastBin=thisBin;
      }
      if(verbose>0){
	std::cout << "fitting " << def << " in: " << fitLowEdge << ".." << fitHighEdge << std::endl;
	if(a!=0) std::cout << "start value [0]:" << a << std::endl;
	if(b!=0) std::cout << "start value [1]:" << b << std::endl;
	if(c!=0) std::cout << "start value [2]:" << c << std::endl;
	if(d!=0) std::cout << "start value [3]:" << d << std::endl;
	if(e!=0) std::cout << "start value [4]:" << e << std::endl;
      }
      // perform fit
      TF1* function=new TF1("function",def,fitLowEdge, fitHighEdge);
      if(a!=0.){
	function->SetParameter(0,a);
	//if(def.Contains("GammaDist")) function->SetParLimits(0,0.6*a,1.4*a);
      }
      if(b!=0.){
	function->SetParameter(1,b);
	if(def.Contains("GammaDist")) function->SetParLimits(1,0.85*b,1.15*b);
      }
      if(c!=0.){
	function->SetParameter(2,c);
	//if(def.Contains("GammaDist")) function->SetParLimits(2,0.6*c,1.4*c);
      }
      if(d!=0.){
	function->SetParameter(3,d);
	//if(def.Contains("GammaDist")) function->SetParLimits(3,0.6*d,1.4*d);
      }
      if(e!=0.){
	function->SetParameter(4,e);
	//if(def.Contains("GammaDist")) function->SetParLimits(4,0.6*e,1.4*e);
      }
      TString option="Q";
      if(verbose>0)option="";
      if(verbose>1)option="V";
      option+="N";
      option+=addOpt;
      option+="R";
      //option+="M";
      input->Fit(function,option,"same",fitLowEdge, fitHighEdge);
      // add part from fit
      result->Add(function);
    }
    // use fit function B
    if(defB!=""&&fitLowEdgeB!=fitHighEdgeB){
      if(verbose>0) std::cout << "will perform fit" << std::endl;
      // empty bins within fit range for output histo
      bool lastBin=false;
      for(int bin=0; bin<=result->GetNbinsX()+1; ++bin){
	bool thisBin=false;
	double binX=result->GetBinCenter(bin);
	if(binX>=fitLowEdgeB&&binX<=fitHighEdgeB){ 
	  result->SetBinContent(bin, 0.);
	  thisBin=true;
	}
	// adapt fit range to binning
	// low
	if(lastBin==false&&thisBin==true) fitLowEdgeB=result->GetBinLowEdge(bin);
	// high
	if(lastBin==true&&thisBin==false) fitHighEdgeB=result->GetBinLowEdge(bin);	
	lastBin=thisBin;
      }
      if(verbose>0){
	std::cout << "fitting " << defB << " in: " << fitLowEdgeB << ".." << fitHighEdgeB << std::endl;
	if(aB!=0) std::cout << "start value [0]:" << aB << std::endl;
	if(bB!=0) std::cout << "start value [1]:" << bB << std::endl;
	if(cB!=0) std::cout << "start value [2]:" << cB << std::endl;
	if(dB!=0) std::cout << "start value [3]:" << dB << std::endl;
	if(eB!=0) std::cout << "start value [4]:" << eB << std::endl;
      }
      // perform fit
      TF1* functionB=new TF1("functionB",defB,fitLowEdgeB, fitHighEdgeB);
      if(aB!=0.){
	functionB->SetParameter(0,aB);
	if(defB.Contains("GammaDist")){
	  if(aB<0.)functionB->SetParLimits(0,1.4*aB,0.6*aB);
	  if(aB>0.)functionB->SetParLimits(0,0.6*aB,1.4*aB);	  
	}
      }
      if(bB!=0.){
	functionB->SetParameter(1,bB);
	if(defB.Contains("GammaDist")){
	  if(bB<0.) functionB->SetParLimits(1,1.15*bB,0.85*bB);
	  if(bB>0.) functionB->SetParLimits(1,0.85*bB,1.15*bB);
	}
      }
      if(cB!=0.){
	functionB->SetParameter(2,cB);
	if(defB.Contains("GammaDist")){
	  if(cB<0.) functionB->SetParLimits(2,1.4*cB,0.6*cB);
	  if(cB>0.) functionB->SetParLimits(2,0.6*cB,1.4*cB);
	}
      }
      if(dB!=0.){
	functionB->SetParameter(3,dB);
	if(defB.Contains("GammaDist")){
	  if(dB<0.) functionB->SetParLimits(3,1.4*dB,0.6*dB);
	  if(dB>0.) functionB->SetParLimits(3,0.6*dB,1.4*dB);
	}
      }
      if(eB!=0.){
	functionB->SetParameter(4,eB);
	if(defB.Contains("GammaDist")){
	  if(eB<0.) functionB->SetParLimits(4,1.4*eB,0.6*eB);
	  if(eB>0.) functionB->SetParLimits(4,0.6*eB,1.4*eB);
	}
      }
      TString option="Q";
      if(verbose>0)option="";
      if(verbose>1)option="V";
      option+="N";
      option+=addOptB;
      option+="R";
      //option+="M";
      input->Fit(functionB,option,"same",fitLowEdgeB, fitHighEdgeB);
      // add part from fit
      result->Add(functionB);
    }
   // use fit function C
    if(defC!=""&&fitLowEdgeC!=fitHighEdgeC){
      if(verbose>0) std::cout << "will perform fit" << std::endl;
      // empty bins within fit range for output histo
      bool lastBin=false;
      for(int bin=1; bin<=result->GetNbinsX()+1; ++bin){
	bool thisBin=false;
	double binX=result->GetBinCenter(bin);
	if(binX>=fitLowEdgeC&&binX<=fitHighEdgeC){ 
	  result->SetBinContent(bin, 0.);
	  thisBin=true;
	}
	// adapt fit range to binning
	// low
	if(lastBin==false&&thisBin==true) fitLowEdgeC=result->GetBinLowEdge(bin);
	// high
	if(lastBin==true&&thisBin==false) fitHighEdgeC=result->GetBinLowEdge(bin);	
	lastBin=thisBin;
      }
      if(verbose>0){
	std::cout << "fitting " << defC << " in: " << fitLowEdgeC << ".." << fitHighEdgeC << std::endl;
	if(aC!=0) std::cout << "start value [0]:" << aC << std::endl;
	if(bC!=0) std::cout << "start value [1]:" << bC << std::endl;
	if(cC!=0) std::cout << "start value [2]:" << cC << std::endl;
	if(dC!=0) std::cout << "start value [3]:" << dC << std::endl;
	if(eC!=0) std::cout << "start value [4]:" << eC << std::endl;
      }
      // perform fit
      TF1* functionC=new TF1("functionC",defC,fitLowEdgeC, fitHighEdgeC);
      if(aC!=0.){
	functionC->SetParameter(0,aC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(0,0.6*aC,1.4*aC);
      }
      if(bC!=0.){
	functionC->SetParameter(1,bC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(1,0.85*bC,1.15*bC);
      }
      if(cC!=0.){
	functionC->SetParameter(2,cC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(2,0.6*cC,1.4*cC);
      }
      if(dC!=0.){
	functionC->SetParameter(3,dC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(3,0.6*dC,1.4*dC);
      }
      if(eC!=0.){
	functionC->SetParameter(4,eC);
	if(defC.Contains("GammaDist")) functionC->SetParLimits(4,0.6*eC,1.4*eC);
      }
      TString option="Q";
      if(verbose>0)option="";
      if(verbose>1)option="V";
      option+="N";
      option+=addOptC;
      option+="R";
      //option+="M";
      input->Fit(functionC,option,"same",fitLowEdgeC, fitHighEdgeC);
      // add part from fit
      result->Add(functionC);
    }
    gErrorIgnoreLevel=initialIgnoreLevel;
    // return result
    return result;
  }

  void DrawTheoryCurve(TString filename="", TString plotname="", bool normalize=true, int smoothFactor=0, int rebinFactor=0, int color=kBlack, int linestyle=1, double rangeLow=-1., double rangeHigh=-1., bool errorbands=false, int errorRebinFactor=0, int errorSmoothFactor=0, int verbose=0, bool drawOnlyErrors=false, bool drawRawPlot=false, TString model="", bool smoothcurves=true, TString PS="Parton")
  {
    // this function draws "plot" from "file" into the active canvas
    // modified quantities: NONE
    // used functions: getTheoryPrediction, histogramStyle, useFittedFunctions, makeVariableBinning
    // used enumerators: NONE
    // filename: name of input file
    // plotname: name of histogram
    // normalize: indicates whether absolute or normalized plot is shown
    // smoothFactor: argument of ->Smooth()
    // rebinFactor:  argument of ->Rebin()
    // color: color of histogram
    // rangeLow:  lower plot range of histo
    // rangeHigh: upper plot range of histo
    // errorbands: additional errorbands are calculated and drawn
    // INFORMATION: if errorbands==true is chosen, the varied up/down plots are expected
    //              to have the same name with an additional "_Up"/"_Down" in the end
    // errorRebinFactor: argument of ->Rebin()
    // errorSmoothFactor: argument of ->Smooth()
    // verbose: level of output
    // drawOnlyErrors: draw only error bands but no central curve
    // drawRawPlot: draw in addition the unbinned and unsmoothed plot 
    //              to see whether rebinning and smoothing has changed the shape 
    // model: indicates theory (madgraph, powhegherwig, powheg or mcatnalo)
    // smoothcurves: indicates wheter smooth or binned curve is drawn
    // PS: "Parton", "Hadron"
    
    // get unified variable name from specific file dependend plotnames
    TString plotname2=plotname;
    if(plotname2.BeginsWith("h")) plotname2.Replace(0,1,"");
    plotname2.ReplaceAll("Vis"    ,"");
    plotname2.ReplaceAll("MC@NLO" ,"");
    plotname2.ReplaceAll("MC@NLO2","");
    plotname2.ReplaceAll("POWHEG","" );
    plotname2.ReplaceAll("Herwig","" );
    plotname2.ReplaceAll("HERWIG","" );
    plotname2.ReplaceAll("herwig","" );
    plotname2.ReplaceAll("Top","top" );
    plotname2.ReplaceAll("TTbarM","ttbarMass" );
    plotname2.ReplaceAll("TTbar","ttbar");
    plotname2.ReplaceAll("Lep","lep");
    plotname2.ReplaceAll("Bottom","bq");
    plotname2.ReplaceAll("Gen","");
    if(plotname2.Contains("/")){
      plotname2.ReplaceAll(getStringEntry(plotname2,1)+"/","");
    }
    //if(plotname2.Contains("bq")&&model=="madgraph"&&smoothcurves) verbose=1;
    // create variable bin edges for non smooth curves
    std::map<TString, std::vector<double> > binning_ = makeVariableBinning();
    // output
    if(verbose>1){
      std::cout << std::endl;
      std::cout << "variable:           " << plotname     << std::endl;
      std::cout << "variable (modifed): " << plotname2    << std::endl;
      if(rangeLow!=-1.||rangeHigh!=-1.) std::cout << "range:              " << rangeLow << " .... " << rangeHigh << std::endl;
      std::cout << "file:               " << filename     << std::endl;
      std::cout << "model:              " << model        << std::endl;
      std::cout << "normalize:          " << normalize    << std::endl;
      std::cout << "smoothFactor:       " << smoothFactor << std::endl;
      std::cout << "rebinFactor:        " << rebinFactor  << std::endl;
      std::cout << "color:              " << color        << std::endl;
      std::cout << "draw errorbands:    " << errorbands   << std::endl;
      std::cout << "draw smooth curves: " << smoothcurves << std::endl;
      if(errorbands){
	std::cout << "errorRebinFactor:   " << errorRebinFactor << std::endl;
	std::cout << "errorSmoothFactor:  " << errorSmoothFactor << std::endl;
      }
    }
    // --- 
    // absolute normalization
    // --- 
    // scale to inclusive xSec
    double weight=ttbarCrossSection;
    // consider efficiency in phase space
    //if(model=="powheg" ) weight*=255743./1227541.;
    //if(model=="mcatnlo") weight*=119134./546459; // MaGraph value for the moment
    // BR correction = <BR W->e/mu> * <BR W->qq> (from pdg)
    weight*=0.5*(0.1075+0.1057)*0.676;
    // N-> <e+mu>
    weight*=0.5;
    // output
    if(!normalize&&(model=="powheg"||model=="powhegherwig"||model=="mcatnlo")&&verbose>0){
      std::cout << "incl xSec(PS): " << 2*weight << " (" << model << ")"<< std::endl;
    }
    // --- 
    // draw raw plot to estimate effect of rebinning and smoothing
    // --- 
    if(drawRawPlot){
      TH1F* raw=getTheoryPrediction(plotname, filename);
      raw->SetMarkerColor(kBlack);
      raw->SetLineColor(kBlack);
      //raw->SetMarkerSize(2);
      //raw->SetMarkerStyle(29);
      if(normalize) raw->Scale(1./(raw->Integral(0,raw->GetNbinsX()+1)));
      // absolute normalization for POWHEG and MCatNlo curve
      if(!normalize&&(model=="powheg"||model=="powhegherwig"||model=="mcatnlo")) raw->Scale(weight/getInclusiveXSec(raw));
      raw->Scale(1./(raw->GetBinWidth(1)));
      raw->Draw("c same"); 
    }
    // --- 
    // central prediction
    // --- 
    // get plot
    TH1F* load=getTheoryPrediction(plotname, filename);
    TH1F* result=0;
    // replace low statistic parts with fitted curve
    if(smoothcurves) result=useFittedFunctions(load, model, plotname, verbose);
    // variable rebinning for binned curve
    TString plotname3 = plotname2;
    plotname3.ReplaceAll("Gen","");
    plotname3.ReplaceAll("gen","");
    if(!smoothcurves){
      // check if maxima of x-axes of loaded plot and new binning is the same
      double newMax=binning_[plotname3][binning_[plotname3].size()-1];
      double oldMax=load->GetBinLowEdge(load->GetNbinsX()+1);
      if(newMax<=oldMax){
	result=(TH1F*)(load->Clone());
	reBinTH1F(*result, binning_[plotname3], verbose);
      }
      // if not: use different rebinning function and ignore bin content in the missing range
      else{
	if(verbose>0){
	  std::cout << "use alternative rebinning function" << std::endl;
	  std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	  std::cout << " theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	}
	result= reBinTH1FIrregularNewBinning(load, binning_[plotname3] , plotname3, false);
      }
      // force under/ and overflow bins to be empty
      result->SetBinContent(0, 0.);
      result->SetBinContent(result->GetNbinsX()+1, 0.);
    }
    // rename
    TString name=plotname2;
    if(name.Contains("Ngenjets")) name.ReplaceAll("Ngenjets", "Njets");
    if(name.Contains("rhosGen" )) name.ReplaceAll("rhosGen" , "rhos" );
    if(model=="powheg") name+="POWHEG";
    else if(model=="powhegherwig") name+="POWHEGHERWIG";
    else if(model=="mcatnlo"){
      if(filename.Contains("catnlo" )) name+="MC@NLO";
      else if(filename.Contains("NtupleCteq6m" )) name+="MC@NLO2";
    }
    else if(model=="madgraph") name+="MADGRAPH";
    if(verbose>0) std::cout << "name of theory curve: " << name << std::endl;
    result->SetName(name);
    result->SetTitle(name);
     // configure style
    histogramStyle(*result, kSig, false, 1.2, color);
    result->SetLineStyle(linestyle);
    // smoothing 1
    if(smoothcurves&&smoothFactor){ 
      if(!plotname.Contains("bqEta")) result->Smooth(smoothFactor);
      else{
	TH1F* duplicate=(TH1F*)result->Clone();
	result->Smooth(smoothFactor); 
	// no smoothing for |eta|>1.6
	for(int bin=0; bin<=result->GetNbinsX(); ++bin){
	  if(std::abs(result->GetBinCenter(bin))>1.6) result->SetBinContent(bin, duplicate->GetBinContent(bin));	  
	}
      }
    }  
    // equal rebinning before
    if(smoothcurves&&rebinFactor) result->Rebin(rebinFactor);
    // normalize to area
    if(normalize) result->Scale(1./(result->Integral(0,result->GetNbinsX()+1)));
    // absolute normalization for POWHEG and MCatNlo curve
    if(!normalize&&(model=="powheg"||model=="powhegherwig"||model=="mcatnlo")) result->Scale(weight/getInclusiveXSec(result));
    // divide by binwidth
    if(smoothcurves) result->Scale(1.0/result->GetBinWidth(1));
    else result=divideByBinwidth(result, false);
    // smoothing 2
    if(smoothcurves&&smoothFactor){
      if(!plotname.Contains("bqEta")) result->Smooth(smoothFactor);
      else{
        TH1F* duplicate=(TH1F*)result->Clone(); 
        result->Smooth(smoothFactor);
        // no smoothing for |eta|>1.6
        for(int bin=0; bin<=result->GetNbinsX(); ++bin){
          if(std::abs(result->GetBinCenter(bin))>1.6) result->SetBinContent(bin, duplicate->GetBinContent(bin));
        }
      }                    
    }   
    // set range
    if(rangeLow!=-1.&&rangeHigh!=-1.) result->GetXaxis()->SetRangeUser(rangeLow, rangeHigh);
    // --- 
    // error bands
    // --- 
    if(errorbands){
      // get values   
      TH1F* loadcentral   =getTheoryPrediction(plotname        , filename);
      TH1F* loadErrorUp   =getTheoryPrediction(plotname+"_Up"  , filename);
      TH1F* loadErrorDown =getTheoryPrediction(plotname+"_Down", filename);
      // for bq and lep error bands: use mcatnlo central curve for correct PS from default running
      // and transfer relative error from parton level studies
      TString plotname4="analyzeTop"+PS+"LevelKinematics";     
      if(PS=="Hadron"){
	if     (plotname.Contains("Lep"   )) plotname4+="Lepton";
	else if(plotname.Contains("Bottom")) plotname4+="Bjets";
      }
      if(plotname.Contains("Vis")) plotname4+="PhaseSpace"; 
      plotname4+="/"+plotname2;

      if(PS=="Hadron"&&(plotname.Contains("Lep")||plotname.Contains("Bottom"))) plotname4+="Gen";
      if(verbose>1) std::cout << plotname << "->" << plotname4 << std::endl;
      TH1F* loadcentral2=0;
      if(verbose>0){
	std::cout << "transfer relative uncertainties to central value from: " << plotname4 << " in " 
		  << groupSpace+"RecentAnalysisRun8TeV/combinedDiffXSecSigMcatnloSummer12PF.root" << std::endl;
      }       
      loadcentral2=getTheoryPrediction(plotname4,groupSpace+"RecentAnalysisRun8TeV/combinedDiffXSecSigMcatnloSummer12PF.root");

      TH1F* central=0;    
      TH1F* central2=0;
      TH1F* ErrorUp=0;
      TH1F* ErrorDown=0;
      if(smoothcurves){
	// replace low statistic parts with fitted curve
	central  =useFittedFunctions(loadcentral , model,plotname,        verbose);
	central2 =useFittedFunctions(loadcentral2, model,plotname,        verbose);      
	ErrorUp  =useFittedFunctions(loadErrorUp,  model,plotname+"_Up",  verbose);
	ErrorDown=useFittedFunctions(loadErrorDown,model,plotname+"_Down",verbose);
	// smoothing 1
	if(errorSmoothFactor){
	  central  ->Smooth(errorSmoothFactor);
	  central2 ->Smooth(errorSmoothFactor);
	  ErrorUp  ->Smooth(errorSmoothFactor);
	  ErrorDown->Smooth(errorSmoothFactor);
	}
	// rebinning
	if(errorRebinFactor){ 
	  central  ->Rebin(errorRebinFactor);
	  central2 ->Rebin(errorRebinFactor);
	  ErrorUp  ->Rebin(errorRebinFactor);
	  ErrorDown->Rebin(errorRebinFactor);
	}
      }
      // variable rebinning for binned curve
      if(!smoothcurves){
	// check if maxima of x-axes of loaded plots and new binning is the same
	double newMax         =binning_[plotname3][binning_[plotname3].size()-1];
	double oldMaxcentral  =loadcentral  ->GetBinLowEdge(load->GetNbinsX()+1);
	double oldMaxcentral2 =loadcentral2 ->GetBinLowEdge(loadcentral2->GetNbinsX()+1);
	double oldMaxErrorUp  =loadErrorUp  ->GetBinLowEdge(load->GetNbinsX()+1);
	double oldMaxErrorDown=loadErrorDown->GetBinLowEdge(load->GetNbinsX()+1);
	if(newMax<=oldMaxcentral){
	  central=(TH1F*)(loadcentral->Clone());
	  reBinTH1F(*central, binning_[plotname3], verbose);
	}
	// central curve- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " central theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  central =reBinTH1FIrregularNewBinning(loadcentral , binning_[plotname3] , plotname3, false);
	}
	if(newMax<=oldMaxcentral2){
	  central2=(TH1F*)(loadcentral2->Clone());
	  reBinTH1F(*central2, binning_[plotname3], verbose);
	}
	// central curve- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " central2 theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  central2=reBinTH1FIrregularNewBinning(loadcentral2, binning_[plotname3] , plotname3, false);
	}	
	if(newMax<=oldMaxErrorUp){
	  ErrorUp=(TH1F*)(loadErrorUp->Clone());
	  reBinTH1F(*ErrorUp, binning_[plotname3], verbose);
	}
	// error up- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " ErrorUp theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  ErrorUp=reBinTH1FIrregularNewBinning(loadErrorUp, binning_[plotname3] , plotname3, false);
	}
	if(newMax<=oldMaxErrorDown){
	  ErrorDown=(TH1F*)(loadErrorDown->Clone());
	  reBinTH1F(*ErrorDown, binning_[plotname3], verbose);
	}
	// error down- if not: use different rebinning function and ignore bin content in the missing range
	else{
	  if(verbose>0){
	    std::cout << "WARNING: in DrawTheoryCurve: maximum of loaded "<< model;
	    std::cout << " ErrorDown theory curve is smaller then the maximum of the requested binning! Will ignore the missing bin content!" << std::endl;
	  }
	  ErrorDown=reBinTH1FIrregularNewBinning(loadErrorDown, binning_[plotname3] , plotname3, false);
	}
      }
      // normalize to area
      if(normalize){
	central  ->Scale(1./(central  ->Integral(0,central  ->GetNbinsX()+1)));
	central2 ->Scale(1./(central2 ->Integral(0,central2 ->GetNbinsX()+1)));
	ErrorUp  ->Scale(1./(central  ->Integral(0,central  ->GetNbinsX()+1)));
	ErrorDown->Scale(1./(central  ->Integral(0,central  ->GetNbinsX()+1)));
      }
      // divide by binwidth
      if(smoothcurves){
	central  ->Scale(1.0/central  ->GetBinWidth(1));
	central2 ->Scale(1.0/central  ->GetBinWidth(1));
	ErrorUp  ->Scale(1.0/ErrorUp  ->GetBinWidth(1));
	ErrorDown->Scale(1.0/ErrorDown->GetBinWidth(1));
      }
      else{
	central  =divideByBinwidth(central  , false);
	central2 =divideByBinwidth(central2 , false);
	ErrorUp  =divideByBinwidth(ErrorUp  , false);
	ErrorDown=divideByBinwidth(ErrorDown, false);
      }
      // smoothing 2
      if(smoothcurves&&errorSmoothFactor){
	central  ->Smooth(errorSmoothFactor);
	central2 ->Smooth(errorSmoothFactor);
	ErrorUp  ->Smooth(errorSmoothFactor);
	ErrorDown->Smooth(errorSmoothFactor);
      }
      // create errorbands
      const int constNBINS = central2->GetNbinsX();
      TGraphAsymmErrors *errorBands = new TGraphAsymmErrors(constNBINS);      
      TString errorBandName=name;
      errorBandName.ReplaceAll("MC@NLO2","MC@NLO");
      errorBandName+="errorBand";
      if(verbose>0) std::cout << "name of error bands: " << errorBandName << std::endl;
      errorBands->SetName(errorBandName);
      errorBands->SetTitle(errorBandName);
      // loop bins
      for(Int_t iBin=1; iBin<=constNBINS; iBin++){
	Double_t centralVal = central  ->GetBinContent(iBin);
	Double_t maxVal     = ErrorUp  ->GetBinContent(iBin);
	Double_t minVal     = ErrorDown->GetBinContent(iBin);
	// transfer relative uncertainty from parton level studies to central value for specific phase space
	Double_t centralValue = central2 ->GetBinContent(iBin);
	Double_t maxValue     = centralValue * (1+(maxVal-centralVal)/centralVal);
	Double_t minValue     = centralValue * (1+(minVal-centralVal)/centralVal);
	if(verbose>1&&(plotname.Contains("Lep")||plotname.Contains("Bottom"))){
	  std::cout << std::endl << plotname << ", bin " << iBin << "/" << central->GetNbinsX() << std::endl;
	  std::cout << "central: " << centralValue << std::endl;
	  std::cout << "up     : " << maxValue     << std::endl;
	  std::cout << "down   : " << minValue     << std::endl;
	  std::cout << "c: " << centralVal << std::endl;
	  std::cout << "u: " << maxVal << " ( = " << 1+(maxVal-centralVal)/centralVal << " )" << std::endl;
	  std::cout << "d: " << minVal << " ( = " << 1+(minVal-centralVal)/centralVal << " )" << std::endl;
	}
	// take care of points outside x-plotrange
	if((rangeLow!=-1.&&central2->GetBinCenter(iBin)<rangeLow)||(rangeHigh!=-1.&&central2->GetBinCenter(iBin)>rangeHigh)){
	  errorBands->SetPoint      (iBin, central2->GetBinCenter(iBin), 0.);
	  //std::cout << "errorBands->SetPoint      (iBin, central2->GetBinCenter(iBin), 0.)" << std::endl;
	  errorBands->SetPointEXlow (iBin, 0.);
	  errorBands->SetPointEXhigh(iBin, 0.);
	}
	else {
	  //std::cout << "errorBands->SetPoint      (iBin, central2->GetBinCenter(iBin), centralValue  )" << std::endl;
	  errorBands->SetPoint(iBin, central2->GetBinCenter(iBin), centralValue  );
	  if(smoothcurves){
	    errorBands->SetPointEXlow (iBin, central2->GetXaxis()->GetBinLowEdge(iBin) );
	    errorBands->SetPointEXhigh(iBin, central2->GetXaxis()->GetBinUpEdge (iBin) );
	  }
	  else{
	    errorBands->SetPointEXlow (iBin, 0.5*central2->GetBinWidth(iBin));
	    errorBands->SetPointEXhigh(iBin, 0.5*central2->GetBinWidth(iBin));
	  }
	}
	// symmetrize errors 
	//double difference=std::abs((maxValue-minValue)*0.5);
	//if((maxValue>centralValue&&minValue>centralValue)||(maxValue<centralValue&&minValue<centralValue)){
	//  difference=(std::abs((maxValue-centralValue))+std::abs((minValue-centralValue)))*0.5;
	//}
	if((maxValue>centralValue&&minValue>centralValue)||(maxValue<centralValue&&minValue<centralValue)){
	  double larger = (minValue > maxValue) ? minValue : maxValue;
	  maxValue=larger;
	  minValue=centralValue - (larger - centralValue); 
	}
	errorBands->SetPointEYhigh(iBin, maxValue-centralValue);
	errorBands->SetPointEYlow (iBin, centralValue-minValue);
	// per bin output
	if(verbose>1){
	  std::cout << "bin " << iBin << std::endl;
	  std::cout << "central TH1: " << std::endl;
	  std::cout << "(<x>,y)= (" << central2->GetBinCenter(iBin) << "," << central2->GetBinContent(iBin) << ")" << std::endl;
	  std::cout << "error input:" << std::endl;
	  std::cout << "central, min, max: " << central2->GetBinContent(iBin);
	  std::cout << ", " << ErrorDown->GetBinContent(iBin) << ", ";
	  std::cout <<  ErrorUp->GetBinContent(iBin) << std::endl;
	  std::cout << "error band" << std::endl;
	  std::cout << "(x,y): (" << errorBands->GetX()[iBin] << "," << errorBands->GetY()[iBin] << ")" << std::endl;
	  std::cout << "error xlow:  " << errorBands->GetErrorXlow(iBin)  << std::endl;
	  std::cout << "error xhigh: " << errorBands->GetErrorXhigh(iBin) << std::endl;
	  std::cout << "error ylow:  " << errorBands->GetErrorYlow(iBin)  << std::endl;
	  std::cout << "error yhigh: " << errorBands->GetErrorYhigh(iBin) << std::endl;
	}
      }
      // plotstyle for error bands
      errorBands->SetFillColor(kGray);
      errorBands->SetFillStyle(1001); // NB: explicitly needed, otherwise filling invisible due to default "0"
      errorBands->SetLineColor  (result->GetLineColor());
      errorBands->SetLineWidth  (result->GetLineWidth());
      errorBands->SetMarkerColor(result->GetMarkerColor());
      // draw error bands
      if(smoothcurves) errorBands->Draw("e3 same");
      else errorBands->Draw("e2 same");
    }    
    // --- 
    //draw central value
    // --- 
    if(!drawOnlyErrors){ 
      result->SetName( TString(result->GetName() ).ReplaceAll("Ngenjets", "Njets"));
      result->SetTitle(TString(result->GetTitle()).ReplaceAll("Ngenjets", "Njets"));
      if(smoothcurves) result->Draw("hist c same"); 
      else DrawSteps(result, "same");//result->Draw("hist same");
      if(verbose>0) std::cout << "theory curve drawn" << std::endl;
    }
  }

  void errorWeightedMeanCombination(const TH1& plotMu, const TH1& plotEl, TH1& plotCombination, unsigned int verbose=0){
    // this function does a simple error weighted mean combination 
    // for all bins of "plotMu" and "plotEl" seperately, saving the
    // result in "plotCombination"
    // verbose: level of output
    // modified quantities: plotCombination
    
    // loop bins
    for(int bin=1; bin<=plotCombination.GetNbinsX(); ++bin){
      if(verbose>1){
	std::cout << ", bin" << bin << ":" << std::endl;
	std::cout << "(" << plotCombination.GetBinLowEdge(bin) << " .... " << plotCombination.GetBinLowEdge(bin+1) << ")" << std::endl;
      }
      // consider only non-empty bins
      if(plotMu.GetBinContent(bin)!=0&&plotEl.GetBinContent(bin)!=0){
	double xSecMu     =plotMu.GetBinContent(bin);
	double xSecEl     =plotEl.GetBinContent(bin);
	double xSecErrorMu=plotMu.GetBinError  (bin);
	double xSecErrorEl=plotEl.GetBinError  (bin);
	double combinedxSec=(xSecMu/(xSecErrorMu*xSecErrorMu)+xSecEl/(xSecErrorEl*xSecErrorEl));
	combinedxSec/=(1/(xSecErrorMu*xSecErrorMu)+1/(xSecErrorEl*xSecErrorEl));
	double combinedxSecError=sqrt(1/(1/(xSecErrorMu*xSecErrorMu)+1/(xSecErrorEl*xSecErrorEl)));
	plotCombination.SetBinContent(bin, combinedxSec     );
	plotCombination.SetBinError  (bin, combinedxSecError);
	if(verbose>1){
	  std::cout << "muon:     " << xSecMu       << " +/- " << xSecErrorMu       << std::endl;
	  std::cout << "elec:     " << xSecEl       << " +/- " << xSecErrorEl       << std::endl;
	  std::cout << "combined: " << combinedxSec << " +/- " << combinedxSecError << std::endl;
	}
      }
    }
  }


    
    // ===============================================================
    // ===============================================================
    

  double regParameter(TString variable, TString decayChannel, int verbose=0, bool fullPS=false, bool tau=false, bool hadronPS=false, TString closureTestSpecifier="", bool probSel=false)
  {
    // this function returns k/value for SVD 
    // unfolding for the corresponding variable
    // modified quantities: NONE
    // used functions: makeVariableBinning
    // used enumerators: NONE
    // variable: name of variable
    // verbose: level of output
    
    // NB: only combined channel is supported
    // NB: only fullPS is supported for top based quantities
    // NB: only visiblePS is supported for b-jet/lepton based quantities
    // NB: all parameters have been derived for the doubleKinFit+prob>0.02 setup
    // NB: for all not supported cases k=Nbins is used!

    double k=-1.0;
    // Default unfolding with number of bins (see makeVariableBinning)
    // SVD unfolding at the moment NOT USED, values are only placeholders
    std::map<TString, std::vector<double> > binning_ = makeVariableBinning();
    if(binning_.count(variable)>0) k = binning_[variable].size()-1;
    
    // Unfolding with optimal tau
    if(tau&&probSel){
      // optimized parameters for each PS, final state, selection and closure test configuration
      if(decayChannel.Contains("combined")){ 
	// top & ttbar:  visible PS (hadron+parton level) not relevant
	if(fullPS){
	  if(closureTestSpecifier==""){
	    if     (variable.Contains("topPtLead")   ) k = MadSpin ?  6.84 : 6.87;
	    else if(variable.Contains("topPtSubLead")) k = MadSpin ?  4.01 : 3.84;
	    else if(variable.Contains("topPtTtbarSys"))k = MadSpin ?  4.74 : 4.86;
	    else if(variable.Contains("topPt")       ) k = MadSpin ?  7.57 : 7.56;
	    else if(variable.Contains("topY" )       ) k = MadSpin ?  9.66 : 9.59;
	    else if(variable.Contains("ttbarPt")     ) k = MadSpin ?  5.66 : 5.85;
	    else if(variable.Contains("ttbarY")      ) k = MadSpin ?  7.79 : 7.68;
	    else if(variable.Contains("ttbarMass")   ) k = MadSpin ?  4.60 : 4.55;
	    else if(variable.Contains("ttbarDelPhi" )) k = MadSpin ? 10.70 : 10.54;
	    else if(variable.Contains("ttbarPhiStar")) k = MadSpin ? 10.62 : 10.42;
	  }
	  else if(closureTestSpecifier=="NoDistort"){
	    if     (variable.Contains("topPtLead")   ) k = MadSpin ?  6.10 : 6.15;
	    else if(variable.Contains("topPtSubLead")) k = MadSpin ?  3.51 : 3.35;
	    else if(variable.Contains("topPtTtbarSys"))k = MadSpin ?  4.06 : 4.15;
	    else if(variable.Contains("topPt")       ) k = MadSpin ?  6.70 : 6.68;
	    else if(variable.Contains("topY" )       ) k = MadSpin ? 10.14 : 10.05;
	    else if(variable.Contains("ttbarPt")     ) k = MadSpin ?  5.79 : 5.98;
	    else if(variable.Contains("ttbarY")      ) k = MadSpin ?  7.94 : 7.80;
	    else if(variable.Contains("ttbarMass")   ) k = MadSpin ?  4.85 : 4.80;
	    else if(variable.Contains("ttbarDelPhi" )) k = MadSpin ? 11.10 : 10.88;
	    else if(variable.Contains("ttbarPhiStar")) k = MadSpin ? 10.92 : 10.67;
	  }
	  else if(closureTestSpecifier=="topPtUp"){
	    if     (variable.Contains("topPtLead")   ) k = MadSpin ?  8.17 : 8.13;
	    else if(variable.Contains("topPtSubLead")) k = MadSpin ?  5.33 : 5.09;
	    else if(variable.Contains("topPtTtbarSys"))k = MadSpin ?  6.09 : 6.23;
	    else if(variable.Contains("topPt")       ) k = MadSpin ?  9.53 : 9.42;
	    else if(variable.Contains("topY" )       ) k = MadSpin ? 10.20 : 10.20;
	    else if(variable.Contains("ttbarPt")     ) k = MadSpin ?  6.02 : 6.19;
	    else if(variable.Contains("ttbarY")      ) k = MadSpin ?  8.08 : 7.95;
	    else if(variable.Contains("ttbarMass")   ) k = MadSpin ?  6.51 : 6.44;
	    else if(variable.Contains("ttbarDelPhi" )) k = MadSpin ? 11.10 : 10.92;
	    else if(variable.Contains("ttbarPhiStar")) k = MadSpin ? 11.01 : 10.80;
	  }
	  else if(closureTestSpecifier=="topPtDown"){
	    if     (variable.Contains("topPtLead")   ) k = MadSpin ?  4.81 : 4.90;
	    else if(variable.Contains("topPtSubLead")) k = MadSpin ?  2.78 : 2.65;
	    else if(variable.Contains("topPtTtbarSys"))k = MadSpin ?  3.21 : 3.28;
	    else if(variable.Contains("topPt")       ) k = MadSpin ?  5.38 : 5.39;
	    else if(variable.Contains("topY" )       ) k = MadSpin ?  9.98 : 9.88;
	    else if(variable.Contains("ttbarPt")     ) k = MadSpin ?  5.62 : 5.78;
	    else if(variable.Contains("ttbarY")      ) k = MadSpin ?  7.80 : 7.66;
	    else if(variable.Contains("ttbarMass")   ) k = MadSpin ?  3.87 : 3.81;
	    else if(variable.Contains("ttbarDelPhi" )) k = MadSpin ? 11.12 : 10.88;
	    else if(variable.Contains("ttbarPhiStar")) k = MadSpin ? 10.81 : 10.47;
	  }
	  else if(closureTestSpecifier=="data"){
	    if     (variable.Contains("topPtLead")   ) k = MadSpin ?  7.20 : 7.18;
	    else if(variable.Contains("topPtSubLead")) k = MadSpin ?  4.33 : 4.13;
	    else if(variable.Contains("topPtTtbarSys"))k = MadSpin ?  4.96 : 5.08;
	    else if(variable.Contains("topPt")       ) k = MadSpin ?  8.02 : 7.98;
	    else if(variable.Contains("topY" )       ) k = MadSpin ? 10.22 : 10.12;
	    else if(variable.Contains("ttbarPt")     ) k = MadSpin ?  5.91 : 6.08;
	    else if(variable.Contains("ttbarY")      ) k = MadSpin ?  8.01 : 7.88;
	    else if(variable.Contains("ttbarMass")   ) k = MadSpin ?  5.72 : 5.69;
	    else if(variable.Contains("ttbarDelPhi" )) k = MadSpin ? 11.07 : 10.77;
	    else if(variable.Contains("ttbarPhiStar")) k = MadSpin ? 10.96 : 10.75;
	  }
	  else if(closureTestSpecifier=="ttbarMassUp"){
	    if     (variable.Contains("topPtLead")   ) k = 5.78;
	    else if(variable.Contains("topPtSubLead")) k = 3.21;
	    else if(variable.Contains("topPtTtbarSys"))k = 3.94;
	    else if(variable.Contains("topPt")       ) k = 6.36;
	    else if(variable.Contains("topY" )       ) k = 9.91;
	    else if(variable.Contains("ttbarPt")     ) k = 5.88;
	    else if(variable.Contains("ttbarY")      ) k = 7.77;
	    else if(variable.Contains("ttbarMass")   ) k = 4.54;
	    else if(variable.Contains("ttbarDelPhi" )) k = 10.92;
	    else if(variable.Contains("ttbarPhiStar")) k = 10.62;
	  }
	  else if(closureTestSpecifier=="ttbarMassDown"){
	    if     (variable.Contains("topPtLead")   ) k = 9.15;
	    else if(variable.Contains("topPtSubLead")) k = 5.37;
	    else if(variable.Contains("topPtTtbarSys"))k = 6.97;
	    else if(variable.Contains("topPt")       ) k = 11.18;
	    else if(variable.Contains("topY" )       ) k = 11.12;
	    else if(variable.Contains("ttbarPt")     ) k = 6.56;
	    else if(variable.Contains("ttbarY")      ) k = 8.01;
	    else if(variable.Contains("ttbarMass")   ) k = 8.35;
	    else if(variable.Contains("ttbarDelPhi" )) k = 10.92;
	    else if(variable.Contains("ttbarPhiStar")) k = 11.07;
	  }
	  else if(closureTestSpecifier=="1000"){
	    if     (variable.Contains("topPtLead")   ) k = MadSpin ?  5.21 : 5.32;
	    else if(variable.Contains("topPtSubLead")) k = MadSpin ?  2.56 : 2.42;
	    else if(variable.Contains("topPtTtbarSys"))k = MadSpin ?  3.10 : 3.17;
	    else if(variable.Contains("topPt")       ) k = MadSpin ?  5.46 : 5.48;
	    else if(variable.Contains("topY" )       ) k = MadSpin ? 10.01 : 9.91;
	    else if(variable.Contains("ttbarPt")     ) k = MadSpin ?  5.73 : 5.90;
	    else if(variable.Contains("ttbarY")      ) k = MadSpin ?  7.86 : 7.73;
	    else if(variable.Contains("ttbarMass")   ) k = MadSpin ?  3.91 : 3.84;
	    else if(variable.Contains("ttbarDelPhi" )) k = MadSpin ? 11.03 : 10.85;
	    else if(variable.Contains("ttbarPhiStar")) k = MadSpin ? 10.85 : 10.60;
	  }
	}
	// lepton&b(jet): full PS and parton level PS  not relevant
	else if(hadronPS){
	  if(closureTestSpecifier==""){
	    if(     variable.Contains("lepPt")       ) k = MadSpin ?  5.62 : 5.65 ;
	    else if(variable.Contains("lepEta")      ) k = MadSpin ?  3.67 : 3.30 ;
	    else if(variable.Contains("bqPt")        ) k = MadSpin ? 12.94 : 12.96;
	    else if(variable.Contains("bqEta")       ) k = MadSpin ?  6.65 : 6.62 ;
	    else if(variable.Contains("bbbarMass")   ) k = MadSpin ?  1.31 : 1.49 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = MadSpin ?  3.41 : 4.93 ;
	    else if(variable.Contains("lbMass")      ) k = MadSpin ? 11.95 : 11.96;
	    else if(variable.Contains("Njets")       ) k = MadSpin ?  1.99 : 2.352;
	    else if(variable.Contains("rho"  )       ) k = MadSpin ?  4.35 : 4.65 ; 
	  }
	  else if(closureTestSpecifier=="NoDistort"){
	    if(     variable.Contains("lepPt")       ) k = MadSpin ?  5.59 : 5.62 ;
	    else if(variable.Contains("lepEta")      ) k = MadSpin ?  3.69 : 3.34 ;
	    else if(variable.Contains("bqPt")        ) k = MadSpin ? 12.21 : 12.35;
	    else if(variable.Contains("bqEta")       ) k = MadSpin ?  6.62 : 6.48 ;
	    else if(variable.Contains("bbbarMass")   ) k = MadSpin ?  1.32 : 1.50 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = MadSpin ?  3.51 : 5.08 ;
	    else if(variable.Contains("lbMass")      ) k = MadSpin ? 12.29 : 12.30;
	    else if(variable.Contains("Njets")       ) k = MadSpin ?  1.87 : 2.22 ;
	    else if(variable.Contains("rho"  )       ) k = MadSpin ?  4.48 : 4.80 ; 
	  }
	  else if(closureTestSpecifier=="topPtDown"){
	    if(     variable.Contains("lepPt")       ) k = MadSpin ?  5.21 : 5.24 ;
	    else if(variable.Contains("lepEta")      ) k = MadSpin ?  3.56 : 3.20 ;
	    else if(variable.Contains("bqPt")        ) k = MadSpin ? 10.63 : 10.64;
	    else if(variable.Contains("bqEta")       ) k = MadSpin ?  6.34 : 6.31 ;
	    else if(variable.Contains("bbbarMass")   ) k = MadSpin ?  1.12 : 1.35 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = MadSpin ?  3.25 : 4.81 ;
	    else if(variable.Contains("lbMass")      ) k = MadSpin ? 12.13 : 12.14;
	    else if(variable.Contains("Njets")       ) k = MadSpin ?  1.78 : 2.11 ;
	    else if(variable.Contains("rho"  )       ) k = MadSpin ?  4.12 : 4.44 ; 
	  }
	  else if(closureTestSpecifier=="topPtUp"){
	    if(     variable.Contains("lepPt")       ) k = MadSpin ?  6.02 : 6.09 ;
	    else if(variable.Contains("lepEta")      ) k = MadSpin ?  3.83 : 3.46 ;
	    else if(variable.Contains("bqPt")        ) k = MadSpin ? 14.75 : 14.82;
	    else if(variable.Contains("bqEta")       ) k = MadSpin ?  6.79 : 6.76 ;
	    else if(variable.Contains("bbbarMass")   ) k = MadSpin ?  1.66 : 1.78 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = MadSpin ?  3.88 : 5.49 ;
	    else if(variable.Contains("lbMass")      ) k = MadSpin ? 12.46 : 12.47;
	    else if(variable.Contains("Njets")       ) k = MadSpin ?  1.98 : 2.355;
	    else if(variable.Contains("rho"  )       ) k = MadSpin ?  5.00 : 5.30 ; 
	  }
	  else if(closureTestSpecifier=="data"){
	    if(     variable.Contains("lepPt")       ) k = MadSpin ?  5.80 : 5.83 ;
	    else if(variable.Contains("lepEta")      ) k = MadSpin ?  3.77 : 3.41 ;
	    else if(variable.Contains("bqPt")        ) k = MadSpin ? 13.49 : 13.51;
	    else if(variable.Contains("bqEta")       ) k = MadSpin ?  6.73 : 6.65 ;
	    else if(variable.Contains("bbbarMass")   ) k = MadSpin ?  1.48 : 1.63 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = MadSpin ?  3.70 : 5.31 ;
	    else if(variable.Contains("lbMass")      ) k = MadSpin ? 12.36 : 12.42;
	    else if(variable.Contains("Njets")       ) k = MadSpin ?  1.92 : 2.29 ;
	    else if(variable.Contains("rho"  )       ) k = MadSpin ?  4.74 : 5.06 ; 
	  }
	  else if(closureTestSpecifier=="ttbarMassUp"){
	    if(     variable.Contains("lepPt")       ) k = 5.44 ;
	    else if(variable.Contains("lepEta")      ) k = 3.32 ;
	    else if(variable.Contains("bqPt")        ) k = 11.71;
	    else if(variable.Contains("bqEta")       ) k = 6.48 ;
	    else if(variable.Contains("bbbarMass")   ) k = 1.46 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = 5.02 ;
	    else if(variable.Contains("lbMass")      ) k = 12.22;
	    else if(variable.Contains("Njets")       ) k = 2.18 ;
	    else if(variable.Contains("rho"  )       ) k = 4.68 ; 
	  }
	  else if(closureTestSpecifier=="ttbarMassDown"){
	    if(     variable.Contains("lepPt")       ) k = 7.15 ;
	    else if(variable.Contains("lepEta")      ) k = 3.50 ;
	    else if(variable.Contains("bqPt")        ) k = 17.62;
	    else if(variable.Contains("bqEta")       ) k = 6.65 ;
	    else if(variable.Contains("bbbarMass")   ) k = 1.95 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = 5.60 ;
	    else if(variable.Contains("lbMass")      ) k = 12.92;
	    else if(variable.Contains("Njets")       ) k = 2.52 ;
	    else if(variable.Contains("rho"  )       ) k = 5.76 ; 
	  }
	  else if(closureTestSpecifier=="1000"){
	    if(     variable.Contains("lepPt")       ) k = MadSpin ?  5.50 : 5.47 ;
	    else if(variable.Contains("lepEta")      ) k = MadSpin ?  3.67 : 3.32 ;
	    else if(variable.Contains("bqPt")        ) k = MadSpin ? 11.23 : 11.22;
	    else if(variable.Contains("bqEta")       ) k = MadSpin ?  6.44 : 6.41 ;
	    else if(variable.Contains("bbbarMass")   ) k = MadSpin ?  1.07 : 1.29 ;
	    else if(variable.Contains("bbbarPt"  )   ) k = MadSpin ?  3.32 : 4.83 ;
	    else if(variable.Contains("lbMass")      ) k = MadSpin ? 12.43 : 12.17;
	    else if(variable.Contains("Njets")       ) k = MadSpin ?  1.85 : 2.20 ;
	    else if(variable.Contains("rho"  )       ) k = MadSpin ?  4.34 : 4.66 ; 
	  }
	}
      }
    }
    else{
      // SVD integer unfolding parameters
      // (NOT USED at the moment, values are only placeholders)
      if     (variable.Contains("topPt")    ) k =  8;
      else if(variable.Contains("topY")     ) k = 10;
      else if(variable.Contains("ttbarPt")  ) k =  6;
      else if(variable.Contains("ttbarY")   ) k = 10;
      else if(variable.Contains("ttbarMass")) k =  7;
      else if(variable.Contains("ttbarPhiStar")) k = 8;
      else if(variable.Contains("ttbarDelPhi" )) k = 4;
      else if(variable.Contains("lepPt")    ) k =  6;
      else if(variable.Contains("lepEta")   ) k = 12;
      else if(variable.Contains("bqPt")     ) k =  4;
      else if(variable.Contains("bqEta")    ) k =  8;
      else if(variable.Contains("bbbarPt"  )) k =  6;
      else if(variable.Contains("bbbarMass")) k =  9;
      else if(variable.Contains("lbMass"   )) k =  3;
      else if(variable.Contains("Njets"    )) k =  5;
      else if(variable.Contains("rho"      )) k =  4;
    }
    // output
    if(verbose>1){
      if(tau) std::cout << "tau";
      else std::cout << "k";
      std::cout << "(" << variable << ") = " << k << std::endl;
    }
    // check result
    if(k<0){
      std::cout << "ERROR in regParameter:" << std::endl; 
      std::cout << "invalid k=" << k << " for variable " << variable << std::endl;
      std::cout << "unknown variable " << variable << " or decayChannel " << decayChannel << std::endl;
    }
    return k;
  }
 
  double getValue(TString fileName, TString canvName, TString plotName, int bin){
    // get file
    TFile* file = TFile::Open(fileName);
    if(!file||(file->IsZombie())){ 
      std::cout << "WARNING in getValue: file " << fileName << " not found" << std::endl;
      return 0.;
    }
    // get canvas
    TCanvas* canv=(TCanvas*)file->Get(canvName);
    if(!canv){ 
      std::cout << "WARNING in getValue: canvas " << canvName << " in file " << fileName << " not found" << std::endl;
      return 0.;
    }
    // get plot
    TH1F* plot  = (TH1F*)canv->GetPrimitive(plotName);
    if(!plot){ 
      std::cout << "WARNING in getValue: plot " << plotName << " within canvas " << canvName << " in file " << fileName << " not found" << std::endl;
      return 0.;
    }
    // get uncertainty value
    if(plot->GetNbinsX()+1<bin||bin<0){
      std::cout << "WARNING in getValue: bin " << bin << " in plot " << plotName << " within canvas " << canvName << " in file " << fileName << " not found (" << plot->GetNbinsX() << " bins)" << std::endl;
      return 0.;
    }
    if(plot->GetBinContent(bin)<0){
      std::cout << "WARNING in getValue: bin " << bin << " in plot " << plotName << " within canvas " << canvName << " in file " << fileName << " has a negative value" << std::endl;
    }
    std::cout << plot->GetBinContent(bin) << std::endl;
    return plot->GetBinContent(bin);
  }

  bool considerInTotalError(unsigned int kSys=sysNo){
    // this function returns false for uncertainties which are not 
    // considered in the total systematic uncertainty and true otherwise
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: systematicVariation

    bool result=true;
    if     (kSys==sysVBosonScaleUp||kSys==sysVBosonScaleDown) result=false;          
    else if(kSys==sysVBosonMatchUp||kSys==sysVBosonMatchDown) result=false;      
    else if(kSys==sysTopMassUp2   ||kSys==sysTopMassDown2   ) result=false;
    else if(kSys==sysTopMassUp3   ||kSys==sysTopMassDown3   ) result=false;
    else if(kSys==sysTopMassUp4   ||kSys==sysTopMassDown4   ) result=false;
    else if(kSys==sysBRUp         ||kSys==sysBRDown         ) result=false;
    else if(kSys==sysUnf                                    ) result=false;
    else if(kSys==sysMad                                    ) result=false;
    else if(kSys==sysGenMCatNLO||kSys==sysGenPowheg||kSys==sysGenPowhegHerwig) result=false;         
    else if(kSys>=ENDOFSYSENUM                              ) result=false;
    //TEST FIXME CAREFUT ATTENTION - the following line is only for testing
    //if(kSys!=sysBtagSFUp&&kSys!=sysBtagSFDown&&kSys!=sysBtagSFShapeUpPt65&&kSys!=sysBtagSFShapeDownPt65&&kSys!=sysBtagSFShapeUpEta0p7&&kSys!=sysBtagSFShapeDownEta0p7) result=false;

    // std::cout << sysLabel(kSys) << ": " << result << std::endl; 
    return result;
  }

  TString pseudoDataFileName(TString closureTestSpecifier, std::string decayChannel = "combined"){
    // this function returns the name of the pseudo data  
    // rootfile for the corresponding closureTestSpecifier
    // modified quantities: NONE
    // used functions: NONE
    // used enumerators: NONE
    // closureTestSpecifier: = \"NoDistort\", \"topPtUp\", \"topPtDown\", \"ttbarMassUp\", \"ttbarMassDown\", \"data\" or \"1000\"

    TString out="";
    // decay channel label
    if(     decayChannel.compare("muon"    )==0) out+="muon"    ;
    else if(decayChannel.compare("electron")==0) out+="electron";
    //else if(decayChannel.compare("combined")==0) out+="combined";
    else{
      std::cout << "ERROR in pseudoDataFile, unknown input decayChannel=" << decayChannel << std::endl;
      exit(0);
    }
    // closure test label
    if      (closureTestSpecifier.Contains("NoDistort"    )) out+="PseudoData19712pb8TeV.root";
    else if (closureTestSpecifier.Contains("topPtUp"      )) out+="PseudoData19712pbReweightedtopPtUp8TeV.root";
    else if (closureTestSpecifier.Contains("topPtDown"    )) out+="PseudoData19712pbReweightedtopPtDown8TeV.root";
    else if (closureTestSpecifier.Contains("ttbarMassUp"  )) out+="PseudoData19712pbReweightedttbarMassUp8TeV.root";
    else if (closureTestSpecifier.Contains("ttbarMassDown")) out+="PseudoData19712pbReweightedttbarMassDown8TeV.root";
    else if (closureTestSpecifier.Contains("data"         )) out+="PseudoData19712pbReweighteddata8TeV.root";
    else if (closureTestSpecifier.Contains("1000"         )) out+="PseudoData19712pbandM1000W100Zprime8TeV.root";
    else if (closureTestSpecifier.Contains("topMass161p5" )) out+="PseudoData19712pbTopMass161p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass163p5" )) out+="PseudoData19712pbTopMass163p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass166p5" )) out+="PseudoData19712pbTopMass166p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass169p5" )) out+="PseudoData19712pbTopMass169p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass175p5" )) out+="PseudoData19712pbTopMass175p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass178p5" )) out+="PseudoData19712pbTopMass178p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass181p5" )) out+="PseudoData19712pbTopMass181p58TeV.root";
    else if (closureTestSpecifier.Contains("topMass184p5" )) out+="PseudoData19712pbTopMass184p58TeV.root";
    
    return out;
  }
    TString pseudoDatalabel(TString closureTestSpecifier){
      // this function return
      // unfolding for the corresponding variable
      // modified quantities: NONE
      // used functions: NONE
      // used enumerators: NONE
      // closureTestSpecifier: = \"NoDistort\", \"topPtUp\", \"topPtDown\", \"ttbarMassUp\", \"ttbarMassDown\", \"data\" or \"1000\"

      TString out="";
      // closure test label
      if      (closureTestSpecifier.Contains("NoDistort"    )) out+="pseudo data=sum of default MC";
      else if (closureTestSpecifier.Contains("topPtUp"      )) out+="pseudo data=softer p_{T}(top)";
      else if (closureTestSpecifier.Contains("topPtDown"    )) out+="pseudo data=harder p_{T}(top)";
      else if (closureTestSpecifier.Contains("ttbarMassUp"  )) out+="pseudo data=harder m(t#bar{t})";
      else if (closureTestSpecifier.Contains("ttbarMassDown")) out+="pseudo data=softer m(t#bar{t})";
      else if (closureTestSpecifier.Contains("data"         )) out+="pseudo data=p_{T}(top)#rightarrow data";
      else if (closureTestSpecifier.Contains("1000"         )) out+="pseudo data=default MC + 1TeV Z'";
      return out;
    }
   
    // ===============================================================
    // ===============================================================

    template <class T> bool SortVectorPrescription(const T& a, const T& b){
	// This function is used to defined
	// a sorting presciption to be used 
	// with std::sort
	return a < b;
    }	

    // ===============================================================
    // ===============================================================
  
    
    

#ifdef DILEPTON_MACRO
}
#endif

#endif
