#include "basicFunctions.h"

void analyzeTopDiffXSecMonitoring(double luminosity = 19712,
				  bool save = true, int verbose=0,
				  TString inputFolderName= "RecentAnalysisRun8TeV_doubleKinFit",
				  //TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/muonDiffXSecData2012ABCDAll.root",
				  //TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/elecDiffXSecData2012ABCDAll.root",
				  TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/elecDiffXSecData2012ABCDAll.root:/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/muonDiffXSecData2012ABCDAll.root",
				  const std::string decayChannel = "combined", 
				  bool withRatioPlot = true, bool extrapolate=true, bool hadron=false, TString addSel="ProbSel")
{ 
  // ============================
  //  Set Root Style
  // ============================
		
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.SetStripDecimals(true);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(2);

  // ============================
  //  Name Conventions
  // ============================

  // a) enumerator sample convention (as defined in basicFunctions.h)

  //  0: kSig      1: kBkg      2: kZjets    3: kWjets
  //  4: kQCD      5: kSTop     6: kDiBos    7: kData
  //  8: kQCDEM1   9: kQCDEM2  10: kQCDEM3  11: kQCDBCE1  12: kQCDBCE2  13: kQCDBCE3  
  // 14: kWW      15: kWZ      16: kZZ      
  // 17: kSTops   18: kSATops  19: kSTopt   20: kSATopt   21: kSToptW   22: kSAToptW
  // 18: ENDOFSAMPLEENUM

  // b) file name convention (implemented in basicFunctions.h)

  // decayChannel+"DiffXSec"+sampleName+systematicVariation+MCProduction+"PF.root"

  // decayChannel      = "electron", "muon", "combined"
  // sampleName        = "Sig", "Bkg", Wjets", "Zjets", "WW", "WZ", "ZZ", "VV", "SingleTopSchannel", 
  //                     "SingleTopTchannel", "SingleTopTWchannel", "QCD"
  // MCProductionCycle = "Summer11","Fall11"
	
  // ============================
  //  Options
  // ============================

  // a) options directly entered when calling function
  // luminosity: [/pb]

  TString lumi = getTStringFromInt(roundToInt((luminosity), false));  

  // save:    save plots?
  // verbose: set detail level of output 
  //          0: no output 
  //          1: std output
  //          2: output for debugging

  // b) options to be configured only once
  // choose if you want to set QCD artificially to 0 to avoid problems with large SF for single events
  bool setQCDtoZero=true;
  //if(addSel=="ProbSel") setQCDtoZero=false;
  // scale ttbar component to measured inclusive xSec	
  bool scaleToMeasured=true;
  // add some shape plots comparing some different ttbar MC generators
  bool compareTTsample=false;
  // produce additional efficiency and acceptance control plots 
  bool effA=true;
  // produce additional control plots for probability selection efficiency
  bool diffProbEff=true;
  if(decayChannel!="combined") effA=false;
  // draw MadGraph+Pythia Theory Shape uncertainties (defined in makeUncertaintyBands within basicFunctions.h)
  bool ttbarTheoryBands=false;
  if(decayChannel!="combined") ttbarTheoryBands=false;
  // activate plots added by hand from tree (using e.g.addDistributions.C/addHistograms.C)
  // will not work with raw .root files!
  bool extra=false;
  // get the .root files from the following folder:
  TString inputFolder = "/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName;
  // see if its 2011 or 2012 data from luminosity
  TString dataSample = "2012";
  // save all plots into the following folder
  TString outputFolder = "./diffXSecFromSignal/plots/"+decayChannel+"/";
  if(dataSample!="") outputFolder+=dataSample+"/";
  outputFolder += (withRatioPlot ? "monitoring/" : "monitoring/withoutRatioPlots/") ;
  // save all plots within a root file named:
  TString outputFileName="diffXSecTopSemi";
  if(decayChannel=="muon"    ) outputFileName+="Mu";
  if(decayChannel=="electron") outputFileName+="Elec";
  if(decayChannel=="combined") outputFileName+="Lep";
  // NOTE: reco plots are the same - whatever PS is chosen but its nice
  //       to have them in the same output file like the cross sections
  // choose phase space
  TString PS="";
  // a) for full PS use extrapolate=true;
  if(!extrapolate) PS="PhaseSpace";
  // b) for restricted phase space:
  // b1) parton PS: hadron=false
  // b2) hadron PS: hadron=true
  TString LV="Parton";
  if(!extrapolate&&hadron) LV="Hadron";
  //if(decayChannel=="electron"||decayChannel=="muon") outputFileName+=dataSample;
  outputFileName+=LV+PS+".root";
  // choose name of the output .pdf file
  TString pdfName="differentialXSecMonitoring"+lumi+"pb";
  // chose wheter SSV instead of CSV btagging control plots should be shown 
  // (WITHOUT ERRORBANDS..)
  bool SSV=false;
  if(SSV) std::cout << std::endl << " for all \"...Tagged/...\" plots the SSV algorithm is chosen!" << std::endl;
  // adjust luminosity and data files for combined control plots
  double luminosityEl=0;
  double luminosityMu=0;
  TString dataFileEl="";
  TString dataFileMu="";
  if(decayChannel=="combined"){
    luminosityEl=constLumiElec;
    luminosityMu=constLumiMuon;
    if(!dataFile.Contains(":")){
      std::cout << "wrong input filenames, should be dataFileEl:dataFileMu, but is ";
      std::cout << dataFile << std::endl;
      exit(0);
    }
    dataFileEl=getStringEntry(dataFile,1 , ":");
    dataFileMu=getStringEntry(dataFile,42, ":");
  }
  else{



  }
  // addSel: xSec from prob selection step?
  //if(!extrapolate) addSel="";
  TString probComplement = addSel=="ProbSel" ? "" : "ProbSel";
  //         0: sysNo                                                       
  //         1: sysLumiUp                   2: sysLumiDown                 
  //         3: sysPUUp                     4: sysPUDown                   
  //         5: sysJESUp                    6: sysJESDown                  
  //         7: sysJERUp                    8: sysJERDown                  
  //         9: sysLepEffSFNormUp          10: sysLepEffSFNormDown         
  //        11: sysLepEffSFShapeUpEta      12: sysLepEffSFShapeDownEta     
  //        13: sysLepEffSFShapeUpPt       14: sysLepEffSFShapeDownPt      
  //        15: sysBtagSFUp                16: sysBtagSFDown               
  //        17: sysBtagSFShapeUpPt65       18: sysBtagSFShapeDownPt65      
  //        19: sysBtagSFShapeUpEta0p7     20: sysBtagSFShapeDownEta0p7    
  //        21: sysMisTagSFUp              22: sysMisTagSFDown             
  //        23: sysTopScaleUp              24: sysTopScaleDown             
  //        25: sysVBosonScaleUp           26: sysVBosonScaleDown          
  //        27: sysSingleTopScaleUp        28: sysSingleTopScaleDown       
  //        29: sysTopMatchUp              30: sysTopMatchDown             
  //        31: sysVBosonMatchUp           32: sysVBosonMatchDown  
  //        33: sysTopMassUp               34: sysTopMassDown  
  //        35: sysTopMassUp2              36: sysTopMassDown2
  //        37: sysTopMassUp3              38: sysTopMassDown3
  //        39: sysTopMassUp4              40: sysTopMassDown4
  //        41: sysQCDUp                   42: sysQCDDown                  
  //        43: sysSTopUp                  44: sysSTopDown                 
  //        45: sysDiBosUp                 46: sysDiBosDown 
  //        47: sysVjetsUp                 48: sysVjetsDown
  //        49: sysBRUp                    50: sysBRDown              
  //        51: sysPDFUp                   52: sysPDFDown 
  //        53: sysUnf,                    54: sysMad,
  //        55: sysHadUp                   56: sysHadDown                  
  //        57: sysGenMCatNLO              58: sysGenPowheg  
  //        59: sysGenPowhegHerwig         60: ENDOFSYSENUM         

  int systematicVariation=sysNo; // MadGraph: sysNo, topPt-reweigthing: sysTest, Powheg+Pythia: sysGenPowheg/sysTestPowheg, Powheg+Herwig: sysGenPowhegHerwig/sysTestPowhegHerwig, McatNLO: sysGenMCatNLO/sysTestMCatNLO
  // use different ttbar MC ("Madgraph", "Powheg", "PowhegHerwig", "McatNLO"), also used for generator uncertainties
  TString ttbarMC="Madgraph";
  if     (systematicVariation==sysGenMCatNLO||systematicVariation==sysTestMCatNLO            ) ttbarMC="Mcatnlo";
  else if(systematicVariation==sysGenPowheg ||systematicVariation==sysTestPowheg             ) ttbarMC="Powheg";
  else if(systematicVariation==sysGenPowhegHerwig ||systematicVariation==sysTestPowhegHerwig ) ttbarMC="PowhegHerwig";
  if(ttbarMC!="Madgraph"){
    pdfName+=ttbarMC;
    outputFolder+=ttbarMC;
    effA=false;
    compareTTsample=false;
    ttbarTheoryBands=false;
  }
  
  // ============================
  //     choose plots
  // ============================

  // a) list plots you would like to see ("folder/plotName") - same as in .root files (for 1D and 2D)

  TString plots1D[ ] = { 
    // (I) preselection
    // (ii) jet monitoring
    //"tightJetKinematicsPreSel/n"  ,
    //"tightJetKinematicsPreSel/en" ,
    //"tightJetKinematicsPreSel/pt" ,
    //"tightJetKinematicsPreSel/eta",
    //"tightJetKinematicsPreSel/phi",
    // (II) before btagging
    // (ii) jet monitoring
    "tightJetKinematics/n"  ,
    //"tightJetKinematicsNjets1/n",
    "tightJetKinematics/en" ,
    "tightJetKinematics/pt" ,
    "tightJetKinematics/eta",
    "tightJetKinematics/phi",
    "tightJetKinematics/ht",
    "tightJetQuality/n_"    ,
    "tightJetQuality/charge",
    "tightJetQuality/nhf"   ,
    "tightJetQuality/nef"   ,
    "tightJetQuality/chf"   ,
    "tightJetQuality/cef"   ,
    "tightJetQuality/ncp"   ,
    "tightLead_0_JetKinematics/en" ,
    "tightLead_0_JetKinematics/pt" ,
    "tightLead_0_JetKinematics/eta",
    "tightLead_0_JetKinematics/phi",
    "tightLead_1_JetKinematics/en" ,
    "tightLead_1_JetKinematics/pt" ,
    "tightLead_1_JetKinematics/eta",
    "tightLead_1_JetKinematics/phi",
    "tightLead_2_JetKinematics/en" ,
    "tightLead_2_JetKinematics/pt" ,
    "tightLead_2_JetKinematics/eta",
    "tightLead_2_JetKinematics/phi",
    "tightLead_3_JetKinematics/en" ,
    "tightLead_3_JetKinematics/pt" ,
    "tightLead_3_JetKinematics/eta",
    "tightLead_3_JetKinematics/phi",
    // (iii) btag monitoring
    "tightJetQuality/btagTrkCntHighPurity",
    "tightJetQuality/btagTrkCntHighEff_"  ,
    "tightJetQuality/btagSimpleSecVtxHighEff_",
    "tightJetQuality/btagSimpleSecVtxHighPur_",
    "tightJetQuality/btagCombSecVtx_",
    "tightJetQuality/btagCombSecVtxMVA_",
    "tightJetQuality/btagJetBProbability_",
    "tightJetQuality/btagJetProbability_",
    "tightJetQuality/btagSoftMuon_",
    "tightJetQuality/btagSoftMuonByPt_",
    "tightJetQuality/btagSoftMuonByIP3d_",
    "bottomJetKinematics/n",
    // (iv) MET monitoring
    "analyzeMETMuon/metEt"   ,
    "analyzeMETMuon/metSumEt",
    // (v) Vertices and pileup
    "PUControlDistributionsBeforeBtagging/pileup",
    "PUControlDistributionsBeforeBtagging/pileup_reweighted",
    "PUControlDistributionsBeforeBtagging/pileup_reweighted_up",
    "PUControlDistributionsBeforeBtagging/pileup_reweighted_down",
    //"PUControlDistributionsBeforeBtagging/npvertex",
    "PUControlDistributionsBeforeBtagging/npvertex_reweighted",
    "PUControlDistributionsBeforeBtagging/npvertex_reweighted_up",
    "PUControlDistributionsBeforeBtagging/npvertex_reweighted_down",
    // (III) after btagging 
    // (ii) jet monitoring
    "tightJetKinematicsTagged/n"  ,
    "tightJetKinematicsTagged/pt" ,
    "tightJetKinematicsTagged/eta",
    "tightJetKinematicsTagged/phi",
    "tightJetKinematicsTagged/ht",
    "tightLead_0_JetKinematicsTagged/pt" ,
    "tightLead_0_JetKinematicsTagged/eta",
    "tightLead_1_JetKinematicsTagged/pt" ,
    "tightLead_1_JetKinematicsTagged/eta",
    "tightLead_2_JetKinematicsTagged/pt" ,
    "tightLead_2_JetKinematicsTagged/eta",
    "tightLead_3_JetKinematicsTagged/pt" ,
    "tightLead_3_JetKinematicsTagged/eta",
    "bottomLead_0_JetKinematicsTagged/pt",
    "bottomLead_1_JetKinematicsTagged/pt",
    "bottomLead_0_JetKinematicsTagged/eta",
    "bottomLead_1_JetKinematicsTagged/eta",
    // (iii) btag monitoring    
    "tightJetQualityTagged/btagCombSecVtx_",
    "bottomJetKinematicsTagged/n"          ,
    // (iv) MET monitoring
    "analyzeMETMuonTagged/metEt"   ,
    "analyzeMETMuonTagged/metPhi"  ,
    "analyzeMETMuonTagged/metSumEt",
    // (v) Vertices and pileup
    "PUControlDistributionsAfterBtagging/pileup",
    "PUControlDistributionsAfterBtagging/pileup_reweighted",
    "PUControlDistributionsAfterBtagging/pileup_reweighted_up",
    "PUControlDistributionsAfterBtagging/pileup_reweighted_down",
    //"PUControlDistributionsAfterBtagging/npvertex",
    "PUControlDistributionsAfterBtagging/npvertex_reweighted",
    "PUControlDistributionsAfterBtagging/npvertex_reweighted_up", 
    "PUControlDistributionsAfterBtagging/npvertex_reweighted_down", 
    // (III) after kinematic fit 
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarAngle",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtTtbarSys",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topY",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topMass",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarY",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarMass",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarHT",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarPhiStar",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lepPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lepEta",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lightqPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lightqEta",   
    "analyzeTopRecoKinematicsKinFit"+addSel+"/bqPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/bqEta",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lepEtaPlus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lepEtaMinus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topEtaPlus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topEtaMinus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lepYPlus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lepYMinus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topYPlus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topYMinus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topYHad",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topYLep",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topYLead",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topYSubLead",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/leadqPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/bbbarPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/bbbarY",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/bbbarMass", 
    "compositedKinematicsKinFit/MWFitJJ",
    "compositedKinematicsKinFit/MJJ",
    "compositedKinematicsKinFit/mHbb",
    "compositedKinematicsKinFit/mbb",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/lbMass",
    // complementary plots - with or without probability selection to have both versions for sure
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarAngle",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topPtTtbarSys",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topY",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topMass",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarY",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarMass",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarHT",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarPhiStar",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/ttbarDelPhi",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lepPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lepEta",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lightqPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lightqEta",   
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/bqPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/bqEta",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lepEtaPlus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lepEtaMinus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topEtaPlus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topEtaMinus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lepYPlus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/lepYMinus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topYPlus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topYMinus",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topYHad",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topYLep",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topYLead",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/topYSubLead",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/leadqPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/bbbarPt",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/bbbarY",
    "analyzeTopRecoKinematicsKinFit"+probComplement+"/bbbarMass",
    // more KinFit plots
    "compositedKinematicsKinFit"+addSel+"/Nbjets",
    "compositedKinematicsKinFit"+addSel+"/Njets",
    "compositedKinematicsKinFit"+addSel+"/rhos",
    "compositedKinematicsKinFit"+addSel+"/leadNonttjetPt",
    "compositedKinematicsKinFit"+addSel+"/leadNonttjetY",
    "compositedKinematicsKinFit"+addSel+"/leadNonttjetEta",
    "analyzeTopRecoKinematicsKinFit/prob", 
    "analyzeTopRecoKinematicsKinFit/chi2",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarDelPhi",
    // gen distributions
    "analyzeTopRecoKinematicsKinFit"+addSel+"/decayChannel", 
    "analyzeTopPartonLevelKinematics/ttbarMass",
    "analyzeTopPartonLevelKinematicsPhaseSpace/ttbarMass",
    "analyzeTopHadronLevelKinematicsLeptonPhaseSpace/lepPtGen",
    // kinfit object shifts
    "compositedKinematicsKinFit/shiftLqPt",
    "compositedKinematicsKinFit/shiftLqEta",
    "compositedKinematicsKinFit/shiftLqPhi",
    "compositedKinematicsKinFit/shiftBqPt",
    "compositedKinematicsKinFit/shiftBqEta",
    "compositedKinematicsKinFit/shiftBqPhi",
    "compositedKinematicsKinFit/shiftLepPt",
    "compositedKinematicsKinFit/shiftLepEta",
    "compositedKinematicsKinFit/shiftLepPhi",
    "compositedKinematicsKinFit/shiftNuPt",
    "compositedKinematicsKinFit/shiftNuEta",
    "compositedKinematicsKinFit/shiftNuPhi",
    // additional top pt plots with splitted contributions
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtPlus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtMinus",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtHad",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtLep",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtLead",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPtSubLead",
    // neutrino plots
    "analyzeTopRecoKinematicsKinFit"+addSel+"/neutrinoPt",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/neutrinoEta",
  };

  TString plots1DLeptons[] = {
    // control plots: common for both lepton types and combined
    // "Lepton" will be replace later by "Electron" or "Muon":
    // e.g: tightLeptonKinematics/n --> tightMuonKinematics/n
    //      tightLeptonKinematics/n --> tightElectronKinematics/n
    // (I) preselection
    //"testIsoLeptonQuality/relIso",
    // (II) before btagging
    "tightLeptonKinematics/n",
    "tightLeptonKinematics/en",
    "tightLeptonKinematics/pt",
    "tightLeptonQuality/relIso",
    // (III) after btagging 
    "tightLeptonQualityTagged/relIso",
    "tightLeptonKinematicsTagged/pt",
   };

  TString plots1Dmu[ ] = { 
    // control plots: muon only
    // (I) preselection
    // (II) before btagging
    "tightMuonKinematics/eta",
    "tightMuonKinematics/y",
    "tightMuonKinematics/phi",
    "tightMuonQuality/nHit",
    "tightMuonQuality/chi2",
    "tightMuonQuality/dB",
    "tightMuonQuality/dz",
    "tightMuonQuality/ecalEn",
    "tightMuonQuality/hcalEn",
    "tightMuonQuality/matches",
   // (III) after btagging 
    "tightMuonKinematicsTagged/eta",
    "tightMuonKinematicsTagged/phi",			 
  };

  TString plots1De[ ] = { 
    // control plots: electron only
    // (I) preselection
    // (II) before btagging
    "tightElectronKinematics/eta",
    "tightElectronKinematics/phi",
    "tightElectronQuality/etaSC",
    "tightElectronQuality/dB",
    "tightElectronQuality/simpleEleId70cIso", 
    "tightElectronQuality/nHitsInner",
    "tightElectronQuality/convDcot",
    "tightElectronQuality/convDist",
    // (III) after btagging 
    "tightElectronKinematicsTagged/eta",
    "tightElectronKinematicsTagged/phi",
  };

  TString plots1Deff[ ] = { 
    // gen distributions for efficiency and acceptance control plots
    "analyzeTopPartonLevelKinematics/topPt",
    "analyzeTopPartonLevelKinematicsPhaseSpace/topPt",
    "analyzeTopPartonLevelKinematics/topY",
    "analyzeTopPartonLevelKinematicsPhaseSpace/topY",
    "analyzeTopPartonLevelKinematics/ttbarPt",
    "analyzeTopPartonLevelKinematicsPhaseSpace/ttbarPt",
    "analyzeTopPartonLevelKinematics/ttbarY",
    "analyzeTopPartonLevelKinematicsPhaseSpace/ttbarY",
    "analyzeTopPartonLevelKinematics/lepPt",
    "analyzeTopPartonLevelKinematicsPhaseSpace/lepPt",
    "analyzeTopPartonLevelKinematics/lepEta",
    "analyzeTopPartonLevelKinematicsPhaseSpace/lepEta",
    "analyzeTopPartonLevelKinematics/bqPt",
    "analyzeTopPartonLevelKinematicsPhaseSpace/bqPt",
    "analyzeTopPartonLevelKinematics/bqEta",
    "analyzeTopPartonLevelKinematicsPhaseSpace/bqEta",
    "analyzeTopPartonLevelKinematics/leadqPt",
    "analyzeTopPartonLevelKinematicsPhaseSpace/leadqPt",
    "analyzeTopPartonLevelKinematics/leadqEta",
    "analyzeTopPartonLevelKinematicsPhaseSpace/leadqEta",
    // reco distributions for efficiency and acceptance control plots
    //"analyzeTopRecoKinematicsKinFit"+addSel+"/leadqPt", // FIXME
    //"analyzeTopRecoKinematicsKinFit"+addSel+"/leadqEta",// FIXME
  };

  TString plots1Dextra[ ] = { 
    // control plots: topPt in topY slices
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_Yless04",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_Ygreater04less10",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_Ygreater10",
    // control plots: topPt in ttbarMass slices
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_ttbarMassless420",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_ttbarMassgreater420less580",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_ttbarMassgreater580",
    // control plots: topPt in ttbarPt slices
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_ttbarPtless40",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_ttbarPtgreater40less130",
    "analyzeTopRecoKinematicsKinFit"+addSel+"/topPt_ttbarPtgreater130",
  };
	
  TString plots2D[ ] = {};	

  // ========================================================
  //  b) list plot axes style
  //  1D: "x-axis title"/"y-axis title"/log/rebin-factor
  //  log = 0 or 1 for linear or logarithmic axis 
  // ========================================================
      
  TString axisLabel1D[ ] = { 
    // (I) preselection
    // (ii) jet monitoring
    //"N_{jets}/Events/0/1",
    //"E(jets)/jets/1/1",
    //"p_{t}(jets)/jets/1/1",
    //"#eta(jets)/jets/0/5",
    //"#phi(jets)/jets/0/10",
    // (II) before btagging
    // (ii) jet monitoring
    "N_{jets};Events;1;1",
    //"N_{jets};Events;1;1",
    "E(jets) [GeV]; jets;1;1",
    "p_{T} #left[GeV#right];Jets;1;1",
    "#eta(jets);Jets;0;5",
    "#phi(jets);Jets;0;10",
    "H_{T} #left[GeV#right];Events;0;50",
    "N(jet constituents);Jets;0;10",
    "jet charge;Jets;0;10"         ,
    "neutral hadron fraction (jets);Jets;1;1"         ,
    "neutral electromagnetic fraction (jets);Jets;0;2",
    "charged hadron fraction (jets);Jets;0;1"         ,
    "charged electromagnetic fraction (jets);Jets;1;1",
    "N_{charged particles} (jets);Jets;0;2"           ,
    "E(lead 1^{st} jet) [GeV];Events;1;2",
    "p_{t}(lead 1^{st} jet) #left[GeV#right];Events;0;5",
    "#eta(lead 1^{st} jet);Events;0;5",
    "#phi(lead 1^{st} jet);Events;0;10",
    "E(lead 2^{nd} jet) [GeV];Events;1;2",
    "p_{t}(lead 2^{nd} jet) #left[GeV#right];Events;0;5",
    "#eta(lead 2^{nd} jet);Events;0;5",
    "#phi(lead 2^{nd} jet);Events;0;10",
    "E(lead 3^{rd} jet) [GeV];Events;1;2",
    "p_{t}(lead 3^{rd} jet) #left[GeV#right];Events;0;5",
    "#eta(lead 3^{rd} jet);Events;0;5",
    "#phi(lead 3^{rd} jet);Events;0;10",
    "E(lead 4^{th} jet) [GeV];Events;1;2",
    "p_{t}(lead 4^{th} jet) #left[GeV#right];Events;0;5",
    "#eta(lead 4^{th} jet);Events;0;5",
    "#phi(lead 4^{th} jet);Events;0;10",
    // (iii) btag monitoring
    "b-discr.(TCHP);Jets;0;2"        ,
    "b-discr.(TCHE);Jets;0;2"	,
    "b-discr.(SSV HEff);Jets;1;2"	,
    "b-discr.(SSV HPur);Jets;0;2"	,
    "b-discr.(CSV);Jets;0;2"	,
    "b-discr.(CSVMVA);Jets;0;2"	,
    "b-discr.(JetBProb);Jets;0;5"	,
    "b-discr.(JetProb);Jets;0;5"	,
    "b-discr.(soft#mu);Jets;0;10"	,
    "b-discr.(soft#muPt);Jets;0;10"  ,                  
    "b-discr.(soft#muIP3d);Jets;0;10",
    "N_{b-jets};Events;1;1"      ,
    // (iv) MET monitoring 
    "#slash{E}_{T} #left[GeV#right];Events;0;10",
    "#SigmaE_{T} [GeV];Events;0;50"  ,
    // (v) Vertices and pileup
    "Number of PU Events;Events;1;1",
    "Number of PU Events (Reweighted);Events;1;1",
    "Number of PU Events (Reweighted sysUp);Events;1;1",
    "Number of PU Events (Reweighted sysDown);Events;1;1",
    "Number of Vertices;Events;1;10",
    //"Number of Vertices (Reweighted);Events;1;10",
    "Number of Vertices (Reweighted sysUp);Events;1;10",
    "Number of Vertices (Reweighted sysDown);Events;1;10",
    // (III) after btagging 
    // (ii) jet monitoring
    "N_{jets};Events;1;1",
    "p_{T} #left[GeV#right];Jets;1;2",
    "#eta(jets);Jets;0;5" ,
    "#phi(jets);Jets;0;10",
    "H_{T} [GeV];Events;0;100",
    "p_{t}(lead 1^{st} jet) #left[GeV#right];Events;1;5",
    "#eta(lead 1^{st} jet);Events;0;5" ,
    "p_{t}(lead 2^{nd} jet) #left[GeV#right];Events;1;5",
    "#eta(lead 2^{nd} jet);Events;0;5" ,
    "p_{t}(lead 3^{rd} jet) #left[GeV#right];Events;1;5",
    "#eta(lead 3^{rd} jet);Events;0;5" ,
    "p_{t}(lead 4^{th} jet) #left[GeV#right];Events;1;5",
    "#eta(lead 4^{th} jet);Events;0;5" ,
    "p_{t}(lead 1^{st} b-tagged jet) #left[GeV#right];Events;1;5",
    "p_{t}(lead 2^{nd} b-tagged jet) #left[GeV#right];Events;1;5",
    "#eta(lead 1^{st} b-tagged jet);Events;0;5" ,
    "#eta(lead 2^{nd} b-tagged jet);Events;0;5" ,
    // (iii) btag monitoring    
    "b-discr.(CSV);Jets;0;2",
    "N_{b-jets};Events;1;1" ,
    // (iv) MET monitoring 
    "#slash{E}_{T} #left[GeV#right];Events;0;20",
    "#phi(#slash{E}_{T});Events;0;5",
    "#SigmaE_{T} [GeV];Events;0;30",
    // (v) Vertices and pileup
    "Number of PU Events;Events;0;1",
    "Number of PU Events (Reweighted);Events;0;1",
    "Number of PU Events (Reweighted sysUp);Events;0;1",
    "Number of PU Events (Reweighted sysDown);Events;0;1",
    "Number of Vertices;Events;0;10",
    //"Number of Vertices (Reweighted);Events;0;10",
    "Number of Vertices (Reweighted sysUp);Events;0;10", 
    "Number of Vertices (Reweighted sysDown);Events;0;10", 
    // (III) after kinematic fit 
    "#angle(t,#bar{t});Events;0;15",
    "p_{T}^{t} #left[GeV#right];Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (t#bar{t} system);Events;0;20",
    "y^{t};Top quarks;0;1",
    "m^{t};Top quarks;0;10",
    "p_{T}^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;20",
    "y^{t#bar{t}};Top-quark pairs;0;1",
    "m^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;50",
    "H_{T}^{t#bar{t}}=#Sigma(E_{T}(jets)) #left[GeV#right];#frac{dN}{dH_{T}^{t#bar{t}}};0;20",
    "#Phi^{*}(t,#bar{t});Events;0;10",
    "p_{T}^{l} #left[GeV#right];N^{l};0;10",    
    "#eta^{l};Leptons;0;1",
    "p_{T}^{q} #left[GeV#right];tt jets;0;20",    
    "#eta^{q};tt jets;0;1",
    "p_{T}^{b} #left[GeV#right];b-jets;0;20",    
    "#eta^{b};b-jets;0;1",
    "#eta^{l+};Events;0;1",
    "#eta^{l-};Events;0;1",
    "#eta^{t+};Events;0;1",
    "#eta^{t-};Events;0;1",
    "y^{l+};Events;0;1",
    "y^{l-};Events;0;1",
    "y^{t+};Events;0;1",
    "y^{t-};Events;0;1",
    "y^{t+};Events;0;1",
    "y^{t-};Events;0;1",
    "y^{Lead.t};Events;0;1",
    "y^{NLead.t};Events;0;1",
    "p_{T}^{leading jet} #left[GeV#right];Events;0;4",
    "p_{T}^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;20",  
    "y^{b#bar{b}}(assigned to t#bar{t} system);Events;0;1",
    "m^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;20",
    "m^{jj} (KinFit W-assignment) #left[GeV#right];Events;1;5",
    "m^{jj} #left[GeV#right];Events;1;10",
    "m^{bb} (KinFit non t#bar{t} b-jets) #left[GeV#right];4 b-tag events;1;50",
    "m^{bb} #left[GeV#right];permutations;1;10",
    "m^{lb} #left[GeV#right];events;1;5",
    // complementary plots - with or without probability selection to have both versions for sure
    "#angle(t,#bar{t});Events;0;15",
    "p_{T}^{t} #left[GeV#right];Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (t#bar{t} system);Events;0;20",
    "y^{t};Top quarks;0;1",
    "m^{t};Top quarks;0;10",
    "p_{T}^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;20",
    "y^{t#bar{t}};Top-quark pairs;0;1",
    "m^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;50",
    "H_{T}(t#bar{t})=#Sigma(E_{T}(jets)) #left[GeV#right];#frac{dN}{dH_{T}(t#bar{t})};0;20",
    "#Phi^{*}(t,#bar{t});Events;0;10",
    "#Delta#Phi(t,#bar{t});Events;0;4",
    "p_{T}^{l} #left[GeV#right];N^{l};0;10",  // keep synchronized with lepPt (Tagged) in axisLabel1DLeptons for tagged/kinfit ratio
    "#eta^{l};Leptons;0;1",
    "p_{T}^{q} #left[GeV#right];tt jets;0;20",    
    "#eta^{q};tt jets;0;1",
    "p_{T}^{b} #left[GeV#right];b-jets;0;20",    
    "#eta^{b};b-jets;0;1",
    "#eta^{l+};Events;0;1",
    "#eta^{l-};Events;0;1",
    "#eta^{t+};Events;0;1",
    "#eta^{t-};Events;0;1",
    "y^{l+};Events;0;1",
    "y^{l-};Events;0;1",
    "y^{t+};Events;0;1",
    "y^{t-};Events;0;1",
    "y^{t+};Events;0;1",
    "y^{t-};Events;0;1",
    "y^{Lead.t};Events;0;1",
    "y^{NLead.t};Events;0;1",
    "p_{T}^{leading jet} #left[GeV#right];Events;0;4",
    "p_{T}^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;20",  
    "y^{b#bar{b}}(assigned to t#bar{t} system);Events;0;1",
    "m^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;20",
    // more KinFit plots
    "N_{bjets};Events;1;1",
    "N_{jets};Events;1;1",
    "#rho_{S}=2m_{0}/m(t#bar{t}+j);Events;1;5",    
    "p_{T}^{leading non t#bar{t}-jet} #left[GeV#right];Events;1;50",
    "y^{leading non t#bar{t}-jet};Events;0;4",
    "#eta^{leading non t#bar{t}-jet};Events;0;4",
    "#chi^{2}-probability;Events;1;25", 
    "#chi^{2};Events;1;10",
    "#Delta#Phi(t,#bar{t});Events;0;4",
    // gen distributions
    "t#bar{t} other decay channel;events;0;1",
    "m^{t#bar{t}} #left[GeV#right] parton all truth;Events;1;1",
    "m^{t#bar{t}} #left[GeV#right] parton lv PS parton truth;Events;1;1",
    "p_{T}^{lep} #left[GeV#right] hadron lv PS truth;Events;1;1",
    // kinfit object shifts
    "#Delta p_{T}^{light jets} #left[GeV#right];Events;0;20",
    "#Delta #eta^{light jets};Events;0;5",
    "#Delta #phi^{light jets};Events;0;10",
    "#Delta p_{T}^{b jets} #left[GeV#right];Events;0;20",
    "#Delta #eta^{b jets};Events;0;5",
    "#Delta #phi^{b jets};Events;0;5",
    "#Delta p_{T}^{lepton} #left[GeV#right];Events;0;20",
    "#Delta #eta^{lepton};Events;0;10",
    "#Delta #phi^{lepton};Events;0;5",
    "#Delta p_{T}^{neutrino} #left[GeV#right];Events;0;20",
    "#Delta #eta^{neutrino};Events;0;20",
    "#Delta #phi^{neutrino};Events;0;10",
    // additional top pt plots with splitted contributions
    "p_{T}^{t+} #left[GeV#right];Events;0;20",
    "p_{T}^{t-} #left[GeV#right];Events;0;20",
    "p_{T}^{had.t} #left[GeV#right];Events;0;20",
    "p_{T}^{lep.t} #left[GeV#right];Events;0;20",
    "p_{T}^{Lead.t} #left[GeV#right];Events;0;20",
    "p_{T}^{NLead.t} #left[GeV#right];Events;0;20",
    // neutrino plots
    "p_{T}^{neutrino} #left[GeV#right];Events;0;10",
    "#eta^{neutrino};Events;0;2",
  };

  TString axisLabel1DLeptons[ ] = {
    // (I) preselection
    //"PF relIso N-1;leptons;0;1",
    // (II) before btagging
    "N;Leptons;0;1"   ,
    "E;Events;0;2"    ,
    "p_{T} [GeV];Leptons;0;1",   
    "PF relIso;Leptons;0;1",
    // (III) after btagging 
    "PF relIso;Leptons;0;1",
    "p^{l}_{T} [GeV];Leptons;0;2",
  };

  TString axisLabel1Dmu[ ] = {
    // (I) preselection
    // (II) before btagging
    "#eta;muons;0;5" ,
    "y;muons;0;5",	
    "#phi;muons;0;5" ,
    "N_{hits}(inner tracker #mu);muons;0;1",
    "#chi^{2} (global trackfit #mu);muons;1;1",
    "d_{xy} (#mu wrt. beamspot) [cm];muons;0;1",
    "d_{z} [cm];muons;0;10",
    "E_{Ecal};muons;1;1",
    "E_{Hcal};muons;1;1",
    "N_{matched #mu segments};muons;0;1",
    // (III) after btagging 
    "#eta;muons;0;1",
    "#phi;muons;0;10",
  };

  TString axisLabel1De[ ] = {
    // (I) preselection
    // (II) before btagging
    "#eta;electrons;0;5",
    "#phi;electrons;0;5",
    "#eta(S.C.);electrons;0;1",
    "d_{xy} (e wrt. beamspot) [cm];electrons;0;1",
    "simpleEleId70cIso;electrons;0;1", 
    "nHitsInner(conv);electrons;0;1",
    "convCot;electrons;0;5",
    "convDist;electrons;0;5",
    // (III) after btagging 
    "#eta;electrons;0;1",
    "#phi;electrons;0;10",
  };

  TString axisLabel1Deff[ ] = { 
    // gen distributions for efficiency and acceptance control plots
    "parton truth p_{T}^{t} #left[GeV#right];Top quarks;0;20",
    "parton truth PS p_{T}^{t} #left[GeV#right];Top quarks;0;20",
    "parton truth y^{t};Top quarks;0;1",
    "parton truth PS y^{t};Top quarks;0;1",
    "parton truth p_{T}^{t#bar{t}} #left[GeV#right];Top quarks;0;20",
    "parton truth PS p_{T}^{t#bar{t}} #left[GeV#right];Top quarks;0;20",
    "parton truth y^{t#bar{t}};Top quarks;0;1",
    "parton truth PS y^{t#bar{t}};Top quarks;0;1",
    "parton truth p_{T}^{l} #left[GeV#right];Events;0;10",
    "parton truth PS p_{T}^{l} #left[GeV#right];Events;0;10",
    "parton truth #eta^{l};Events;0;1",
    "parton truth PS #eta^{l};Events;0;1",
    "parton truth p_{T}^{b} #left[GeV#right];b-quarks;0;20",
    "parton truth PS p_{T}^{b} #left[GeV#right];b-quarks;0;20",
    "parton truth #eta^{b};b-quarks;0;1",
    "parton truth PS #eta^{b};b-quarks;0;1",
    "parton truth p_{T}^{lead q} #left[GeV#right];Events;0;20",
    "parton truth PS p_{T}^{lead q} #left[GeV#right];Events;0;20",
    "parton truth #eta^{lead q};Events;0;1",
    "parton truth PS #eta^{lead q};Events;0;1",
    // reco distributions for efficiency and acceptance control plots
    //"p_{T}^{lead q} #left[GeV#right];Events;0;20",// FIXME
    //"#eta^{lead q};Events;0;1",// FIXME
  };

  TString axisLabel1Dextra[ ] = { 
    // control plots: topPt in topY slices
    "p_{T}^{t} #left[GeV#right] (y^{t}<0.4);Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (0.4<y^{t}<1.0);Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (y^{t}>1.0);Top quarks;0;20",
    // control plots: topPt in ttbarMass slices
    "p_{T}^{t} #left[GeV#right] (m^{t#bar{t}}<420);Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (420<m^{t#bar{t}}<580);Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (m^{t#bar{t}}>580);Top quarks;0;20",
    // control plots: topPt in ttbarPt slices
    "p_{T}^{t} #left[GeV#right] (p_{T}^{t#bar{t}}<40);Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (40<p_{T}^{t#bar{t}}<130);Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (p_{T}^{t#bar{t}}>130);Top quarks;0;20",
  };

  // 2D: "x-axis title"/"y-axis title"
  TString axisLabel2D[ ] = {};

  // ======================================================================
  //  Count number of plots + cross-check to number of axis-label defined 
  // ======================================================================
  unsigned int N1Dplots = sizeof(plots1D)/sizeof(TString) + sizeof(plots1DLeptons)/sizeof(TString);
  if(decayChannel=="electron") N1Dplots += sizeof(plots1De  )/sizeof(TString);
  if(decayChannel=="muon"    ) N1Dplots += sizeof(plots1Dmu )/sizeof(TString);
  if(effA                    ) N1Dplots += sizeof(plots1Deff)/sizeof(TString);
  if(extra                   ) N1Dplots += sizeof(plots1Dextra)/sizeof(TString);
  unsigned int N2Dplots = sizeof(plots2D)/sizeof(TString);
  // get number fo axis labels and check if it corresponds to number of plots
  unsigned int Naxislabels = sizeof(axisLabel1D)/sizeof(TString) + sizeof(axisLabel1DLeptons)/sizeof(TString);
  if(decayChannel=="electron") Naxislabels += sizeof(axisLabel1De  )/sizeof(TString);
  if(decayChannel=="muon"    ) Naxislabels += sizeof(axisLabel1Dmu )/sizeof(TString);
  if(effA                    ) Naxislabels += sizeof(axisLabel1Deff)/sizeof(TString);
  if(extra                   ) Naxislabels += sizeof(axisLabel1Dextra)/sizeof(TString);
  if(N1Dplots != Naxislabels){
    std::cout << "ERROR - 1D plots: Number of plots and axis label do not correspond .... Exiting macro!" << std::endl;
    exit(1);
  }
  if(N2Dplots != sizeof(axisLabel2D)/sizeof(TString)){
    std::cout << "ERROR - 2D plots: Number of plots and axis label do not correspond .... Exiting macro!" << std::endl;
    exit (1);
  }
  // run automatically in batch mode if there are many canvas
  if((N1Dplots+N2Dplots)>15) gROOT->SetBatch();
	
  // ===================================
  //  Open the standard analysis files
  // ===================================
  if(verbose>1) std::cout << "get files" << std::endl;
  std::map<unsigned int, TFile*> files_, filesMu_, filesEl_;
  if(decayChannel!="combined") files_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFile, decayChannel, ttbarMC);
  else{
    filesMu_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFileMu, "muon"    , ttbarMC);
    filesEl_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFileEl, "electron", ttbarMC);
  }

  // =====================
  //  Loading histos
  // =====================
  if(verbose>1) std::cout << "get histos" << std::endl;
  // collect all plot names in vector (first 1D, then 2D)
  std::vector<TString> plotList_, plotListEl_, plotListMu_;
  plotList_.insert(plotList_.begin(), plots1D,        plots1D        + sizeof(plots1D       )/sizeof(TString));
  plotList_.insert(plotList_.end(),   plots1DLeptons, plots1DLeptons + sizeof(plots1DLeptons)/sizeof(TString));
  if(decayChannel=="electron") plotList_.insert(plotList_.end(), plots1De,  plots1De  + sizeof(plots1De )/sizeof(TString));
  if(decayChannel=="muon"    ) plotList_.insert(plotList_.end(), plots1Dmu, plots1Dmu + sizeof(plots1Dmu)/sizeof(TString)); 
  if(effA) plotList_.insert(plotList_.end(), plots1Deff,     plots1Deff     + sizeof(plots1Deff    )/sizeof(TString));
  if(extra                   ) plotList_.insert(plotList_.end(), plots1Dextra, plots1Dextra + sizeof(plots1Dextra)/sizeof(TString)); 
  plotList_.insert( plotList_.end(), plots2D, plots2D + sizeof(plots2D)/sizeof(TString) );	

  if(decayChannel == "combined"){
    plotListEl_ = plotList_;
    plotListMu_ = plotList_;
  }

  // Renaming: for the naming used to save the plots and later to access them internally
  if(verbose>1) std::cout << "internal renaming" << std::endl;
  if(decayChannel != "combined"){
    for (std::vector<TString>::iterator iter = plotList_.begin(); iter != plotList_.end(); iter++){
      // Rename lepton plots according to chosen decay channel
      if (iter->Contains("Lepton")&&!iter->Contains("Hadron")){     
	if (decayChannel=="electron")  iter->ReplaceAll("Lepton","Electron");
	else if (decayChannel=="muon") iter->ReplaceAll("Lepton","Muon");    
      }
      // Rename KinFitProbSel plots
      if (iter->Contains("KinFitProbSel")&&iter->Contains("composited")) iter->ReplaceAll("KinFitProbSel","ProbSel"); 
    }
  }
  else{
    // electron plots
    for (std::vector<TString>::iterator iter = plotListEl_.begin(); iter != plotListEl_.end(); iter++){
      // Rename lepton plots to electron label and change pt to et
      if (iter->Contains("Lepton")&&!iter->Contains("Hadron")) iter->ReplaceAll("Lepton","Electron");
      if (iter->Contains("tightElectronKinematics/pt")) iter->ReplaceAll("tightElectronKinematics/pt","tightElectronKinematics/et");
      // Rename KinFitProbSel plots
      if (iter->Contains("KinFitProbSel")&&iter->Contains("composited")) iter->ReplaceAll("KinFitProbSel","ProbSel"); 
    }
    // muon plots    
    for (std::vector<TString>::iterator iter = plotListMu_.begin(); iter != plotListMu_.end(); iter++){
      // Rename lepton plots to muon label
      if (iter->Contains("Lepton")&&!iter->Contains("Hadron")) iter->ReplaceAll("Lepton","Muon");
      // Rename KinFitProbSel plots
      if (iter->Contains("KinFitProbSel")&&iter->Contains("composited")) iter->ReplaceAll("KinFitProbSel","ProbSel"); 
    }
  }

  // container for all histos (1D&2D)
  // example for acess: histo_["plotName"][sampleNr]
  std::map< TString, std::map <unsigned int, TH1F*> > histo_, histoEl_, histoMu_;
  std::map< TString, TH1F* > histoErrorBand_, histoErrorBandPU_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2_, histo2El_, histo2Mu_;
  // a) Save all histos from plotList_ that exist in files_ into histo_ and histo2_
  // b) Count total number of plots as Nplots
  // c) Special action required for control plots regarding PU and nVertex:
  //    - Compare the reweighted histos in MC to the same but unweighted histo in data
  //    - Strip off appendices from histo-names as redundant part of name to get the according plots from data
  if(verbose>0) std::cout << std::endl;
  int Nplots=0;
  std::vector<TString> vecRedundantPartOfNameInData;
  vecRedundantPartOfNameInData.push_back("_reweighted_up");
  vecRedundantPartOfNameInData.push_back("_reweighted_down");
  vecRedundantPartOfNameInData.push_back("_reweighted");
  // do NOT save plots in open rootfiles
  gROOT->cd();
  if(decayChannel!="combined"){
    getAllPlots(files_, plotList_, histo_, histo2_, N1Dplots, Nplots, verbose, decayChannel, &vecRedundantPartOfNameInData, SSV);
    closeStdTopAnalysisFiles(files_);
  }
  else{
    getAllPlots(filesEl_, plotListEl_, histoEl_, histo2El_, N1Dplots, Nplots, verbose, "electron", &vecRedundantPartOfNameInData, SSV);
    getAllPlots(filesMu_, plotListMu_, histoMu_, histo2Mu_, N1Dplots, Nplots, verbose, "muon"    , &vecRedundantPartOfNameInData, SSV);
    closeStdTopAnalysisFiles(filesEl_);
    closeStdTopAnalysisFiles(filesMu_);
  }

  // ==========================================
  //  Lumiweighting for choosen luminosity
  // ==========================================
  // scale every histo in histo_ and histo2_ to the corresponding luminosity
  // NOTE: luminosity [/pb]
  if(decayChannel!="combined") scaleByLuminosity(plotList_, histo_, histo2_, N1Dplots, luminosity, verbose, systematicVariation, decayChannel);
  else{
    scaleByLuminosity(plotListMu_, histoMu_, histo2Mu_, N1Dplots, luminosityMu, verbose, systematicVariation, "muon"    ,ttbarMC);
    scaleByLuminosity(plotListEl_, histoEl_, histo2El_, N1Dplots, luminosityEl, verbose, systematicVariation, "electron",ttbarMC);
  }

  // =======================================================
  //  Add single top channels and DiBoson contributions
  // =======================================================
  // for SingleTop and DiBoson samples:
  // every plot plotList_ existing in sTop/ diBoson .root file
  // will be combined and saved in the histo_ and histo2_ map
  // reCreate: reCreate combined plots if they are already existing
  bool reCreate=false;
  if(decayChannel!="combined")  AddSingleTopAndDiBoson(plotList_, histo_, histo2_, N1Dplots, verbose, reCreate, decayChannel);
  else{
    AddSingleTopAndDiBoson(plotListMu_, histoMu_, histo2Mu_, N1Dplots, verbose, reCreate, "muon"    ); 
    AddSingleTopAndDiBoson(plotListEl_, histoEl_, histo2El_, N1Dplots, verbose, reCreate, "electron");
  }

  // =======================================================
  //  Add decay channels for combined control plot
  // =======================================================
  if(decayChannel=="combined"){
    // loop samples
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      if(verbose>1) std::cout << sampleLabel(sample, decayChannel) << std::endl;
      // loop plots
      for(unsigned int plot=0; plot<plotList_.size(); ++plot){
	//if(sample==kData&&plotList_[plot].Contains("topPt")&&plotList_[plot].Contains("Lead")) verbose=2;
	if(verbose>1) std::cout << "comb: " << plotList_[plot] << ", mu: " << plotListMu_[plot] << ", el: " << plotListEl_[plot] << std::endl;
	// a) 1D
	if((plot<N1Dplots)&&(histoMu_.count(plotListMu_[plot])>0)&&(histoMu_[plotListMu_[plot]].count(sample)>0)&&(histoEl_.count(plotListEl_[plot])>0)&&(histoEl_[plotListEl_[plot]].count(sample)>0)){ 
	  //if(plotListMu_[plot].Contains("shift")){
	  //  std::cout << "muon: " << plotListMu_[plot] << ", sample " << sampleLabel(sample, "muon") << " (" << histoMu_[plotListMu_[plot]][sample]->GetNbinsX() << " bins, [" << histoMu_[plotListMu_[plot]][sample]->GetBinLowEdge(1) << ".." << histoMu_[plotListMu_[plot]][sample]->GetBinLowEdge(histoMu_[plotListMu_[plot]][sample]->GetNbinsX()+1) << "]" << std::endl;
	  //  std::cout << "electron: " << plotListEl_[plot] << ", sample " << sampleLabel(sample, "electron") << " (" << histoEl_[plotListEl_[plot]][sample]->GetNbinsX() << " bins, [" << histoEl_[plotListEl_[plot]][sample]->GetBinLowEdge(1) << ".." << histoEl_[plotListEl_[plot]][sample]->GetBinLowEdge(histoEl_[plotListEl_[plot]][sample]->GetNbinsX()+1) << "]" << std::endl; 
	  //}
	  if(verbose>1) std::cout << "-> 1D" << std::endl;
	  histo_[plotList_[plot]][sample]=     (TH1F*)(histoMu_[plotListMu_[plot]][sample]->Clone());
	  if (plotListMu_[plot]=="tightMuonKinematics/pt"){
	    histo_[plotList_[plot]][sample]->Add((TH1F*)(histoEl_["tightElectronKinematics/et"][sample]->Clone()));
	  }
	  else histo_[plotList_[plot]][sample]->Add((TH1F*)(histoEl_[plotListEl_[plot]][sample]->Clone()));
	}
	// b) 2D
	else if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	  if(verbose>1) std::cout << "-> 2D" << std::endl;
	  histo2_[plotList_[plot]][sample]=     (TH2F*)(histo2Mu_[plotListMu_[plot]][sample]->Clone());
	  histo2_[plotList_[plot]][sample]->Add((TH2F*)(histo2El_[plotListEl_[plot]][sample]->Clone()));
	}
	else{
	  if(verbose>1){
	    std::cout << "WARNING: " << plotList_[plot] << "NOT FOUND" << std::endl;
	    if(!histoMu_.count(plotListMu_[plot])>0) std::cout << "in histoMu" << std::endl;
	    if( histoMu_.count(plotListMu_[plot])>0&&!histoMu_[plotListMu_[plot]].count(sample)>0) std::cout << "in mu sample " << sampleLabel(sample, "muon"    ) << std::endl;
	    if(!histoEl_.count(plotListEl_[plot])>0) std::cout << "in histoEl" << std::endl;
	    if( histoEl_.count(plotListEl_[plot])>0&&!histoEl_[plotListEl_[plot]].count(sample)>0) std::cout << "in el sample " << sampleLabel(sample, "electron") << std::endl;
	  }
	}
	//if(sample==kData&&plotList_[plot].Contains("topPt")&&plotList_[plot].Contains("Lead")) verbose=0;
      } 
    }
  }

  // ============================
  //  Configure histograms
  // ============================
  // needs: plotList_, histo_, histo2_, N1Dplots, axisLabel_, axisLabel1D, axisLabel2D
  std::vector<TString> axisLabel_;
  axisLabel_.insert(axisLabel_.begin(), axisLabel1D,        axisLabel1D        + sizeof(axisLabel1D)/sizeof(TString));
  axisLabel_.insert(axisLabel_.end(),   axisLabel1DLeptons, axisLabel1DLeptons + sizeof(axisLabel1DLeptons)/sizeof(TString));
  if(decayChannel=="electron") axisLabel_.insert( axisLabel_.end(), axisLabel1De , axisLabel1De +sizeof(axisLabel1De )/sizeof(TString) );
  if(decayChannel=="muon"    ) axisLabel_.insert( axisLabel_.end(), axisLabel1Dmu, axisLabel1Dmu+sizeof(axisLabel1Dmu)/sizeof(TString) );
  if(effA) axisLabel_.insert(axisLabel_.end(), axisLabel1Deff,     axisLabel1Deff     + sizeof(axisLabel1Deff)/sizeof(TString));
  if(extra) axisLabel_.insert(axisLabel_.end(), axisLabel1Dextra,     axisLabel1Dextra     + sizeof(axisLabel1Dextra)/sizeof(TString));
  axisLabel_.insert( axisLabel_.end()  , axisLabel2D, axisLabel2D + sizeof(axisLabel2D)/sizeof(TString) );
  if(verbose>1){
    std::cout << "(plot, x Axis label, y Axis label, log scale?, rebinning factor):" << std::endl;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      std::cout << plotList_[plot];
      std::cout << ": " << getStringEntry(axisLabel_[plot],1,";");
      std::cout << ", " << getStringEntry(axisLabel_[plot],2,";");
      if(plot<N1Dplots){
	std::cout << ", " << getStringEntry(axisLabel_[plot],3,";");
	std::cout << ", " << getStringEntry(axisLabel_[plot],4,";");
      }
      std::cout << std::endl;
    }
  }

  // loop samples
  for(unsigned int sample=kSig; sample<=kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // a) 1D
      if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){ 
	// Default configuration
	histogramStyle(*histo_[plotList_[plot]][sample], sample, true);
	// Special configurations
	if(getStringEntry(plotList_[plot], 2)=="PartonJetDRall") histo_[plotList_[plot]][sample] -> SetNdivisions(816);
	if(plotList_[plot].Contains("tightJetKinematics") &&plotList_[plot].Contains("/n")) histo_[plotList_[plot]][sample]->SetNdivisions(010); 
	if(plotList_[plot].Contains("bottomJetKinematics")&&plotList_[plot].Contains("/n")) histo_[plotList_[plot]][sample]->SetNdivisions(010);
	if(plotList_[plot].Contains("Njets")&&plotList_[plot].Contains("compositedKinematics")) histo_[plotList_[plot]][sample]->SetNdivisions(010);
	// set QCD to 0
	if(setQCDtoZero&&sample==kQCD&&(plotList_[plot].Contains("Tagged")||plotList_[plot].Contains("AfterBtagging")||plotList_[plot].Contains("analyzeTopReco")||plotList_[plot].Contains("compositedKinematicsKinFit"))) histo_[plotList_[plot]][sample]->Scale(0.);
      }
      // b) 2D
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	histStyle2D( *histo2_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), getStringEntry(axisLabel_[plot],1,";"), getStringEntry(axisLabel_[plot],2,";"));
      }
    }
  }
  
  // ============================
  //  Event composition
  // ============================
  double shiftxSec=0.;
  std::map< TString, std::map <unsigned int, double> > events_;
  // a) get event numbers
  std::vector<TString> selection_;
  selection_.push_back("tightJetKinematics/n"      );
  selection_.push_back("tightJetKinematicsTagged/n");
  selection_.push_back("compositedKinematicsKinFit/shiftLepPt");  
  selection_.push_back("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarMass");  
  unsigned int MCBG=ENDOFSAMPLEENUM;
  events_[selection_[0]][MCBG]=0;
  events_[selection_[1]][MCBG]=0; 
  events_[selection_[2]][MCBG]=0;
  events_[selection_[3]][MCBG]=0;
  // loop pretagged/tagged
  for(unsigned int step=0; step<selection_.size(); ++step){
    // loop samples
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      if(verbose>1) std::cout << selection_[step] << " in " << sampleLabel(sample, decayChannel) << std::endl;
      // save number
      events_[selection_[step]][sample]=histo_[selection_[step]][sample]->Integral(0,histo_[selection_[step]][sample]->GetNbinsX()+1);
      // add non ttbar MC
      if(sample>kSig&&sample<kData) events_[selection_[step]][MCBG]+=events_[selection_[step]][sample];
    }
  }
  // b) print composition (only if ratio is also drawn)
  double xSec=0;
  double NttbarFull  = histo_["analyzeTopPartonLevelKinematics/ttbarMass"          ][kSig]->Integral(0, histo_["analyzeTopPartonLevelKinematics/ttbarMass"          ][kSig]->GetNbinsX()+1);
  double NttbarPSHad = histo_["analyzeTopHadronLevelKinematicsLeptonPhaseSpace/lepPtGen"][kSig]->Integral(0, histo_["analyzeTopHadronLevelKinematicsLeptonPhaseSpace/lepPtGen"][kSig]->GetNbinsX()+1);
  double NttbarPS    = histo_["analyzeTopPartonLevelKinematicsPhaseSpace/ttbarMass"][kSig]->Integral(0, histo_["analyzeTopPartonLevelKinematicsPhaseSpace/ttbarMass"][kSig]->GetNbinsX()+1);
  double BR=BRPDG(systematicVariation);
  double A   =NttbarPS   /NttbarFull;
  double AHad=NttbarPSHad/NttbarFull;
  if(verbose>2) std::cout << "NallttbarGenSG=" << NttbarFull << std::endl;
  // loop pretagged/tagged
  for(unsigned int step=0; step<selection_.size(); ++step){    
    // print label
    if(verbose>=0&&withRatioPlot){
      switch (step){
      case 0 : std::cout << std::endl << " Event composition ---- pre-tagged, derived from: "+selection_[step]          << std::endl; break;
      case 1 : std::cout << std::endl << " Event composition ---- after b-tagging, derived from: "+selection_[step]     << std::endl; break;
      case 2 : std::cout << std::endl << " Event composition ---- after kinematic fit (before prob selection), derived from: "+selection_[step] << std::endl; break;
      case 3 : std::cout << std::endl << " Event composition ---- after kinematic fit, derived from: "+selection_[step] << std::endl; break;
      default: break;
      }
    }
    // all MC events
    double NAllMC=events_[selection_[step]][kSig]+events_[selection_[step]][MCBG];
    double NData=events_[selection_[step]][kData];
    if(verbose>2) std::cout << "data/MC=" << NData/NAllMC << std::endl;
    // print yield and composition
    if(verbose>=0&&withRatioPlot){
      std::cout << " Events observed in data: " << NData << std::endl;
      std::cout << " Events expected from MC: " << NAllMC << std::endl;
      std::cout << " -> data/MC ratio: "           << NData/NAllMC << std::endl << std::endl;
      std::cout << " Expected event composition:"   << std::endl; 
      std::cout << " ttbar SG:   " << std::setprecision(4) << std::fixed << events_[selection_[step]][kSig  ] / NAllMC << std::endl;
      std::cout << " ttbar BG:   " << std::setprecision(4) << std::fixed << events_[selection_[step]][kBkg  ] / NAllMC << std::endl;
      std::cout << " W + Jets:   " << std::setprecision(4) << std::fixed << events_[selection_[step]][kWjets] / NAllMC << std::endl; 
      std::cout << " Z + Jets:   " << std::setprecision(4) << std::fixed << events_[selection_[step]][kZjets] / NAllMC << std::endl;
      std::cout << " QCD:        " << std::setprecision(4) << std::fixed << events_[selection_[step]][kQCD  ] / NAllMC;
      if(setQCDtoZero&&events_[selection_[step]][kQCD]==0.0)  std::cout  << " (artificially set to 0)";
      std::cout << std::endl;
      std::cout << " Single Top: " << std::setprecision(4) << std::fixed << events_[selection_[step]][kSTop ] / NAllMC << std::endl;
    }
    // investigate single top composition
    if(decayChannel != "combined"){
      double Nschannel= histo_[selection_[step]][kSTops  ]->Integral(0,histo_[selection_[step]][kSTops  ]->GetNbinsX()+1);
      Nschannel+=       histo_[selection_[step]][kSATops ]->Integral(0,histo_[selection_[step]][kSATops ]->GetNbinsX()+1);
      double Ntchannel= histo_[selection_[step]][kSTopt  ]->Integral(0,histo_[selection_[step]][kSTopt  ]->GetNbinsX()+1);
      Ntchannel+=       histo_[selection_[step]][kSATopt ]->Integral(0,histo_[selection_[step]][kSATopt ]->GetNbinsX()+1);
      double NtWchannel=histo_[selection_[step]][kSToptW ]->Integral(0,histo_[selection_[step]][kSToptW ]->GetNbinsX()+1);
      NtWchannel+=      histo_[selection_[step]][kSAToptW]->Integral(0,histo_[selection_[step]][kSAToptW]->GetNbinsX()+1);
      double NallCh=events_[selection_[step]][kSTop];
      if(withRatioPlot){
	std::cout << "             (s: " << std::setprecision(2) << std::fixed << Nschannel/NallCh;
	std::cout << ", t: "  << std::setprecision(2) << std::fixed << Ntchannel/NallCh;
	std::cout << ", tW: " << std::setprecision(2) << std::fixed << NtWchannel/NallCh << ")" << std::endl;
      }
    }
    if(verbose>=0&&withRatioPlot) std::cout << " DiBoson:    " << std::setprecision(4) << std::fixed << events_[selection_[step]][kDiBos] / NAllMC << std::endl;
    double NnonTtbarBG=NAllMC-events_[selection_[step]][kSig]-events_[selection_[step]][kBkg];
    double ttbarSigFrac=events_[selection_[step]][kSig]/(events_[selection_[step]][kSig]+events_[selection_[step]][kBkg]);
    double eff   =events_[selection_[step]][kSig]/NttbarPS;
    double effHad=events_[selection_[step]][kSig]/NttbarPSHad;
    double luminosity2=luminosity;
    if(decayChannel=="combined") luminosity2=luminosityMu+luminosityEl;
    xSec=(NData-NnonTtbarBG)*ttbarSigFrac/(A*eff*BR*luminosity2);
    if(selection_[step].Contains("shift")) shiftxSec=xSec;
    if(verbose>=0&&withRatioPlot){
      std::cout << std::endl;
      std::cout << "    inclusive cross section: " << std::setprecision(2) << std::fixed << xSec;
      std::cout << " ( +" << std::setprecision(2) << std::fixed << 100*(xSec-ttbarCrossSection)/ttbarCrossSection;
      std::cout << "% wrt. NNLO+NNLL)" << std::endl;
      std::cout << "      N(nonttBG):    " << std::setprecision(3) << std::fixed << NnonTtbarBG << std::endl;
      std::cout << "      tt SG frac:    " << std::setprecision(3) << std::fixed << ttbarSigFrac << std::endl;
      std::cout << "      BR:            " << std::setprecision(3) << std::fixed << BR  << std::endl;
      std::cout << "      efficiency:    " << std::setprecision(3) << std::fixed << eff << std::endl;
      std::cout << "      acceptance:    " << std::setprecision(3) << std::fixed << A   << std::endl;
      std::cout << "hadLV efficiency:    " << std::setprecision(3) << std::fixed << effHad << std::endl;
      std::cout << "hadLV acceptance:    " << std::setprecision(3) << std::fixed << AHad   << std::endl;
      std::cout << "      eff. lumi[pb]: " << std::setprecision(3) << std::fixed << luminosity2 << std::endl;
      std::cout << "Nrec, sel.:    " << std::setprecision(3) << std::fixed << events_[selection_[step]][kSig] << std::endl;
      std::cout << "Ngen, all:    " << std::setprecision(3) << std::fixed << NttbarFull << std::endl;
      std::cout << "Ngen, vis:    " << std::setprecision(3) << std::fixed << NttbarPSHad << std::endl;
      std::cout << "Nrec/Ngen,vis:    " << std::setprecision(3) << std::fixed << events_[selection_[step]][kSig]/NttbarPSHad << std::endl;
      std::cout << "Nrec/Ngen,all:    " << std::setprecision(3) << std::fixed << events_[selection_[step]][kSig]/NttbarFull  << std::endl;

      std::cout << std::endl;
      if(step==selection_.size()-1&&scaleToMeasured){
	std::cout << " scaling ttbar MC wrt. sigma(ttbar)=";
	std::cout << std::setprecision(2) << std::fixed << xSec << " pb" << std::endl;
	for(unsigned int step2=0; step2<selection_.size(); ++step2){
	  std::cout << "   - " << selection_[step2] << std::endl;
	  double NscaledAllMC=events_[selection_[step2]][kSig]+events_[selection_[step2]][MCBG]+(xSec/ttbarCrossSection-1)*(events_[selection_[step2]][kSig]+events_[selection_[step2]][kBkg]);
	  std::cout << "      Events in scaled MC : " << NscaledAllMC << std::endl;
	  std::cout << "      ratio data/scaled MC: " << events_[selection_[step2]][kData]/NscaledAllMC << std::endl << std::endl;
	}
      }
    }
  }
  std::cout << std::endl;
  if(!(verbose>=0&&withRatioPlot)) std::cout << " The event composition is only printed when running the monitoring macro using the option 'withRatioPlot=true'" << std::endl;  
	
  // ============================
  //  Rebinning 1D histograms
  // ============================
  // loop samples
  for(unsigned int sample=kSig; sample<=kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      TString plotName=plotList_[plot];
      // check if plot exists and is 1D
      if((plot<N1Dplots)&&(plotExists(histo_, plotName, sample))){
	// equidistant binning
	double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4,";")).c_str());
	if(reBinFactor>1){
	  equalReBinTH1(reBinFactor, histo_, plotName, sample);
	}
      }
    }
  }

  // ======================================
  //  scale ttbar component to measured xSec
  // ======================================
  if(scaleToMeasured){
    if(verbose>2) std::cout << "data/NNLO+NNLL=" << xSec/ttbarCrossSection << std::endl;
    for(unsigned int sample=kSig; sample<=kBkg; ++sample){
      for(unsigned int plot=0; plot<plotList_.size(); ++plot){
	if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	  if(plotList_[plot].Contains("shift")) histo_[plotList_[plot]][sample]->Scale(shiftxSec/ttbarCrossSection);
	  else histo_[plotList_[plot]][sample]->Scale(xSec/ttbarCrossSection);
	}
      }
    }
  }

  // ===============================================================
  //  Errors for uncertainty bands from ttbar Xsec and luminosity
  // ===============================================================
  if(verbose>0) std::cout << std::endl << " Start calculating error bands for 1D plots .... ";
  if(!SSV&&!scaleToMeasured) makeUncertaintyBands(histo_, histoErrorBand_, plotList_, N1Dplots);
  else if(ttbarTheoryBands) makeTheoryUncertaintyBands(histo_, histoErrorBand_, plotList_, plotListEl_, plotListMu_, Nplots, N1Dplots, inputFolder, dataFileMu, dataFileEl, vecRedundantPartOfNameInData, luminosityMu, luminosityEl);
  if(verbose>0) std::cout << " .... Finished." << std::endl; 

  // ========================================================
  //  Create Legends
  // ========================================================
  unsigned int Nlegends=0; 
  TLegend *leg  = new TLegend(); 
  TLegend *leg0 = new TLegend(0.05, 0.15, 1.05, 0.9);
  TLegend *leg1 = new TLegend(0.05, 0.15, 1.05, 0.9);
  Nlegends+=3;
  legendStyle(*leg,"");
  legendStyle(*leg0,"After Selection, Before b-Tagging");
  legendStyle(*leg1,"After Selection and b-Tagging");

  // Set legend entries
  // - First entry: data
  // - Following entries: MC samples (if existing)
  // - Last entry: uncertainty band

  // Data:
  TString lumilabel = Form("%2.1f fb^{-1}",luminosity/1000);
  leg ->AddEntry(histo_[plotList_[0]][kData], sampleLabel(kData,decayChannel),"P");
  leg0->AddEntry(histo_[plotList_[0]][kData], sampleLabel(kData,decayChannel)+", "+lumilabel,"PL");
  leg1->AddEntry(histo_[plotList_[0]][kData], sampleLabel(kData,decayChannel)+", "+lumilabel,"PL");
 
  // MC samples (add only if sample exists in at least one plot, then quit plot-loop to avoid duplication of entries)
  TLegend *legTheo=(TLegend*)leg->Clone();
  legTheo->SetHeader("#splitline{after KinFit reco,}{#splitline{BG subtracted}{  }}");
  for(unsigned int sample=kSig; sample<kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size()-1; ++plot){  // <plotList_.size()-1, because last entry is for data (see above)
      if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)){
	TString leglabel=sampleLabel(sample,decayChannel);
	if((systematicVariation==sysTest)&&(sample==kSig||sample==kBkg)) leglabel="weighted "+leglabel;
	leg ->AddEntry(histo_[plotList_[plot]][sample], leglabel, "F");
	leg0->AddEntry(histo_[plotList_[plot]][sample], leglabel, "F");
	leg1->AddEntry(histo_[plotList_[plot]][sample], leglabel, "F");
	break;
      }
    }
  }

  // Uncertainty band in legend
  if(histoErrorBand_.size() > 0 && plotList_.size() > 0){
    TString unclabel = ttbarTheoryBands ? "t#bar{t} model unc." : "Norm. Uncertainty";
    //TString unclabel = "JES unc.";
    leg ->AddEntry(histoErrorBand_[plotList_[0]], unclabel, "F");
    leg0->AddEntry(histoErrorBand_[plotList_[0]], unclabel, "F");
    leg1->AddEntry(histoErrorBand_[plotList_[0]], unclabel, "F");
  }
  //std::cout << "last:           "  << std::setprecision(4) << std::fixed << 1./sqrt(1250./lumiweight(kSig, luminosity, sysNo          , "muon")) << std::endl;
  //std::cout << "first:           "  << std::setprecision(4) << std::fixed << 1./sqrt(360./lumiweight(kSig, luminosity, sysNo          , "muon")) << std::endl;
  //std::cout << "sysTopMatchUp:   " << 1./sqrt(2500./lumiweight(kSig, luminosity, sysTopMatchUp  , "muon")) << std::endl;;
  //std::cout << "sysTopMatchDown: " << 1./sqrt(2500./lumiweight(kSig, luminosity, sysTopMatchDown, "muon")) << std::endl;;
  //std::cout << "sysTopScaleUp:   " << 1./sqrt(2500./lumiweight(kSig, luminosity, sysTopScaleUp  , "muon")) << std::endl;;
  //std::cout << "sysTopScaleDown: " << 1./sqrt(2500./lumiweight(kSig, luminosity, sysTopScaleDown, "muon")) << std::endl;;

  // ===================================
  //  Change 1D plots into stack plots
  // ===================================
  // loop plots -> all 1D plots will become stacked plots
  if(verbose>1) std::cout << std::endl;
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    createStackPlot(plotList_, histo_, plot, N1Dplots, verbose, decayChannel);
  }
  if(verbose>1) std::cout << std::endl;
  if(verbose>2) std::cout << "reweighted difference: " << histo_["analyzeTopRecoKinematicsKinFit"+addSel+"/topPt"][kSig]->Integral(0,histo_["analyzeTopRecoKinematicsKinFit"+addSel+"/topPt"][kSig]->GetNbinsX()+1)-histo_["analyzeTopRecoKinematicsKinFit"+addSel+"/topPt"][kData]->Integral(0,histo_["analyzeTopRecoKinematicsKinFit"+addSel+"/topPt"][kData]->GetNbinsX()+1) << std::endl;

  // differential probability cut efficiency
  std::map<TString, std::vector<double> > errEff_;
  bool foundProbEff=false;
  TLegend* legProbEff= new TLegend(0.8, 0.8, 0.9, 0.88);
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    if(plotList_[plot].Contains("ProbSel")&&diffProbEff){
      //naming
      TString nameOri=plotList_[plot];
      TString nameEff=nameOri;
      nameEff.ReplaceAll("ProbSel","ProbEff");
      TString namePur=nameOri;
      namePur.ReplaceAll("ProbSel","");
      if(histo_[namePur][kSig]&&histo_[plotList_[plot]][kSig]){
	plotList_.push_back(nameEff);
	TString newLabel=getStringEntry(axisLabel_[plot],1,";")+";eff(prob>0.02);"+getStringEntry(axisLabel_[plot],3,";")+";"+getStringEntry(axisLabel_[plot],4,";");
	// debug output
	if(verbose>1){
	  std::cout << std::endl << plotList_[plot] << std::endl;
	  std::cout << "-> " << nameEff << std::endl;
	  std::cout << axisLabel_[plot] << std::endl;
	  std::cout << "-> " << newLabel << std::endl;
	}
	axisLabel_.push_back(newLabel);
	// MC histo
	histo_[nameEff][kSig]        =(TH1F*)histo_[plotList_[plot]][kSig]->Clone();
	histo_[nameEff][kSig]->Divide((TH1F*)histo_[namePur][kSig]->Clone());
	// MC unc: consider only statistical uncertainty of ttbar events
	double NdenomGes=0;
	std::string lepTag= decayChannel=="combined" ? "muon" : decayChannel;
	double lumiTag= decayChannel=="combined" ? (luminosityEl+luminosityMu)/2 : luminosity;
	for(int bin=0; bin<=histo_[nameEff][kSig ]->GetNbinsX(); ++bin){
	  double ratio=histo_[nameEff][kSig]->GetBinContent(bin);
	  double Ndenom=histo_[namePur][kSig]->GetBinContent(bin)-histo_[namePur][kZjets]->GetBinContent(bin);
	  Ndenom/=lumiweight(kSig, lumiTag, systematicVariation, lepTag);
	  if(scaleToMeasured) Ndenom/=(xSec/ttbarCrossSection);
	  if(verbose>1&&plotList_[plot].Contains("ttbarMass")){
	    std::cout << "ttbarMass bin #" << bin << ": " << Ndenom << std::endl;
	    NdenomGes+=Ndenom;
	    if(bin==histo_[nameEff][kSig]->GetNbinsX()){
	      std::cout << "SUM=" << NdenomGes << std::endl;
	      std::cout << "lumiweight=" << lumiweight(kSig, lumiTag, systematicVariation, lepTag) << std::endl;
	    }
	  }
	  double temp_err=sqrt((ratio*(1-ratio))/Ndenom);
	  histo_[nameEff][kSig]->SetBinError(bin,temp_err);
	}
	// MC histo style
	histo_[nameEff][kSig ]->SetLineColor(kRed);
	histo_[nameEff][kSig ]->SetFillColor(kWhite);
	// data histo
	histo_[nameEff][kData]        =(TH1F*)histo_[plotList_[plot]][kData]->Clone();
	DivideTwoYields(histo_[nameEff][kData], (TH1F*)histo_[namePur][kData]->Clone());
	N1Dplots++;
	// uncertainty for MC/data ratio 
	std::vector<double> tempErr_;
	for(int bin=1; bin<=histo_[nameEff][kData]->GetNbinsX(); ++bin){
	  double b =histo_[nameEff][kSig ]->GetBinContent(bin);
	  double Db=histo_[nameEff][kSig ]->GetBinError  (bin);
	  double a =histo_[nameEff][kData]->GetBinContent(bin);
	  double Da=histo_[nameEff][kData]->GetBinError  (bin);
	  tempErr_.push_back(sqrt((a*a*Db*Db)/(b*b*b*b)+(Da*Da)/(b*b)));
	}
	errEff_[nameEff]=tempErr_;
	// legend
	if(!foundProbEff){
	  legendStyle(*legProbEff,"");
	  legProbEff->AddEntry(histo_[nameEff][kSig ], "all MC", "L");
	  legProbEff->AddEntry(histo_[nameEff][kData], "data"  , "P");
	  foundProbEff=true;
	}
	if(plotList_[plot].Contains("lepPt")){
	  // lepPt ratio before/ after KinFit
	  TString nameBef="tightLeptonKinematicsTagged/pt";
	  TString nameEff2=nameEff;
	  nameEff2.ReplaceAll("ProbEff","KinFitEff");
	  if(histo_.count(nameBef)>0){
	    // A MC histos
	    // adapt MC bin range
	    histo_[nameEff2][kSig]=cutOffTH1((TH1F*)histo_[namePur][kSig ]->Clone(namePur+"2"), histo_[nameBef][kSig]->GetBinLowEdge(histo_[nameBef][kSig]->GetNbinsX()+1));
	    TH1F* tempSG=cutOffTH1((TH1F*)histo_[nameBef][kSig]->Clone("tempSG"), histo_[nameEff2][kSig]->GetBinLowEdge(histo_[nameEff2][kSig]->GetNbinsX()+1));
	    // MC ratio
	    histo_[nameEff2][kSig]->Divide(cutOffTH1(tempSG, histo_[nameEff2][kSig]->GetBinLowEdge(histo_[nameEff2][kSig]->GetNbinsX()+1)));
	    // MC uncertainty: 
	    for(int bin=0; bin<=histo_[nameEff2][kSig ]->GetNbinsX(); ++bin){
	      double ratio=histo_[nameEff2][kSig]->GetBinContent(bin);
	      double Ndenom=histo_[nameBef][kSig]->GetBinContent(bin)-histo_[nameBef][kZjets]->GetBinContent(bin);
	      Ndenom/=lumiweight(kSig, lumiTag, systematicVariation, lepTag);
	      if(scaleToMeasured) Ndenom/=(xSec/ttbarCrossSection);
	      double temp_err=sqrt((ratio*(1-ratio))/Ndenom);
	      histo_[nameEff2][kSig]->SetBinError(bin,temp_err);
	    }
	    // MC histo style
	    histo_[nameEff2][kSig]->SetLineColor(kRed);
	    histo_[nameEff2][kSig]->SetFillColor(kWhite);
	    // B data histo
	    TH1F* tempDa=(TH1F*)histo_[nameBef][kData]->Clone("tempDa");
	    // adapt data bin range
	    histo_[nameEff2][kData]=cutOffTH1((TH1F*)histo_[namePur][kData]->Clone(), tempDa->GetBinLowEdge(tempDa->GetNbinsX()+1));
	    // data ratio
	    DivideTwoYields(histo_[nameEff2][kData], cutOffTH1(tempDa, histo_[nameEff2][kData]->GetBinLowEdge(histo_[nameEff2][kData]->GetNbinsX()+1)));
	    // uncertainty for MC/data ratio 
	    std::vector<double> tempErr2_;
	    for(int bin=1; bin<=histo_[nameEff2][kData]->GetNbinsX(); ++bin){
	      double b =histo_[nameEff2][kSig ]->GetBinContent(bin);
	      double Db=histo_[nameEff2][kSig ]->GetBinError  (bin);
	      double a =histo_[nameEff2][kData]->GetBinContent(bin);
	      double Da=histo_[nameEff2][kData]->GetBinError  (bin);
	      tempErr2_.push_back(sqrt((a*a*Db*Db)/(b*b*b*b)+(Da*Da)/(b*b)));
	    }
	    errEff_[nameEff2]=tempErr2_;
	    // add new plot to list
	    plotList_.push_back(nameEff2);
	    TString newLabel2=getStringEntry(axisLabel_[plot],1,";")+";eff(kinFit);"+getStringEntry(axisLabel_[plot],3,";")+";"+getStringEntry(axisLabel_[plot],4,";");
	    axisLabel_.push_back(newLabel2);
	    N1Dplots++;
	  }
	}
      }
    }
  }
  
  // ============================
  // introduce some BG subtracted
  // shape plots
  // ============================
  // loop all plots
  if(compareTTsample){
    bool firstBGsub=true;
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      // search for leading non ttbar jet
      if((plotList_[plot].Contains("compositedKinematicsKinFit/leadNontt")||plotList_[plot].Contains("compositedKinematicsKinFit/Njets"))&&!plotList_[plot].Contains("BGSubNorm")){
	// create new Entry
	TString  newName=plotList_[plot];
	newName+="BGSubNorm";
	plotList_ .insert(plotList_ .begin()+plot+1, newName         );
	axisLabel_.insert(axisLabel_.begin()+plot+1, axisLabel_[plot]);
	axisLabel_[plot+1].ReplaceAll("events", "norm. events");
	N1Dplots++;
	// BG subtracted data and pure signal prediction
	histo_[newName][kSig]=(TH1F*)histo_[plotList_[plot]][kSig]->Clone(newName);
	histo_[newName][kSig]->Add(histo_[plotList_[plot]][kBkg],-1);
	histo_[newName][kData]=(TH1F*)histo_[plotList_[plot]][kData]->Clone(newName);
	histo_[newName][kData]->Add(histo_[plotList_[plot]][kBkg],-1);
	// shape normalization
	histo_[newName][kSig ]->Scale(1./histo_[newName][kSig ]->Integral(0,histo_[newName][kSig ]->GetNbinsX()+1));
	histo_[newName][kData]->Scale(1./histo_[newName][kData]->Integral(0,histo_[newName][kData]->GetNbinsX()+1)); 
	// adapt style 
	histo_[newName][kSig]->SetLineColor(histo_[newName][kSig]->GetFillColor());
	histo_[newName][kSig]->SetLineWidth(2);
	histo_[newName][kSig]->SetFillStyle(0);
	// create other theory prediction plots
	unsigned int koneMca= (decayChannel!="electron" ? kSigMca : kBkgMca);
	unsigned int ktwoMca= (decayChannel=="combined" ? kBkgMca : 999);
	unsigned int konePow= (decayChannel!="electron" ? kSigPow : kBkgPow);
	unsigned int ktwoPow= (decayChannel=="combined" ? kBkgPow : 999);
	files_[kSigMca]=TFile::Open(inputFolder+"/"+TopFilename(kSigMca, systematicVariation, "muon"    ));
	files_[kBkgMca]=TFile::Open(inputFolder+"/"+TopFilename(kSigMca, systematicVariation, "electron"));
	files_[kSigPow]=TFile::Open(inputFolder+"/"+TopFilename(kSigPow, systematicVariation, "muon"    ));
	files_[kBkgPow]=TFile::Open(inputFolder+"/"+TopFilename(kSigPow, systematicVariation, "electron"));
	if(files_[kSigMca]&&files_[kBkgMca]){
	  histo_[newName][kSigMca]=(TH1F*)(files_[koneMca]->Get(plotList_[plot])->Clone(newName));
	  if(ktwoMca!=999) histo_[newName][kSigMca]->Add((TH1F*)(files_[ktwoMca]->Get(plotList_[plot])->Clone(newName)));
	}
	if(files_[kSigPow]&&files_[kBkgPow]){
	  histo_[newName][kSigPow]=(TH1F*)(files_[konePow]->Get(plotList_[plot])->Clone(newName));
	  if(ktwoPow!=999) histo_[newName][kSigPow]->Add((TH1F*)(files_[ktwoPow]->Get(plotList_[plot])->Clone(newName)));
	}
	double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4,";")).c_str());
	if(reBinFactor>1){
	  if(files_[kSigPow]&&files_[kBkgPow]) equalReBinTH1(reBinFactor, histo_, newName, kSigPow);
	  if(files_[kSigMca]&&files_[kBkgMca]) equalReBinTH1(reBinFactor, histo_, newName, kSigMca);
	}
	if(files_[kSigMca]&&files_[kBkgMca]){
	  histo_[newName][kSigMca]->Scale(1./histo_[newName][kSigMca]->Integral(0,histo_[newName][kSigMca]->GetNbinsX()+1));
	  histo_[newName][kSigMca]->SetLineColor(constMcatnloColor);
	  histo_[newName][kSigMca]->SetMarkerColor(constMcatnloColor);
	  histo_[newName][kSigMca]->SetLineWidth(2);
	  histo_[newName][kSigMca]->SetLineStyle(constNnloStyle);
	  histo_[newName][kSigMca]->SetFillStyle(0);
	}
	if(files_[kSigPow]&&files_[kBkgPow]){
	  histo_[newName][kSigPow]->Scale(1./histo_[newName][kSigPow]->Integral(0,histo_[newName][kSigPow]->GetNbinsX()+1));
	  histo_[newName][kSigPow]->SetLineColor(constPowhegColor );
	  histo_[newName][kSigPow]->SetMarkerColor(constPowhegColor );
	  histo_[newName][kSigPow]->SetLineWidth(2);
	  histo_[newName][kSigPow]->SetLineStyle(constPowhegStyle);
	  histo_[newName][kSigPow]->SetFillStyle(0);
	}
	if(firstBGsub){
	  legTheo->AddEntry(histo_[newName][kSig], "MadGraph+Pythia", "L");
	  if(files_[kSigPow]&&files_[kBkgPow]) legTheo->AddEntry(histo_[newName][kSigPow], "POWHEG+Pythia"  , "L");
	  if(files_[kSigMca]&&files_[kBkgMca]) legTheo->AddEntry(histo_[newName][kSigMca], "MC@NLO+Herwig"  , "L");
	  firstBGsub=false;
	}
      }
    }
  }

  // =================================
  // introduce some acceptance and 
  // efficiency generator comparisons
  // =================================
  TLegend *legTheo2= new TLegend();
  legendStyle(*legTheo2,"t#bar{t} SG simulations");
  if(effA&&compareTTsample){
    // get POWHEG and MC@NLO ttbar signal files
    files_[kSigMca]=TFile::Open(inputFolder+"/"+TopFilename(kSigMca, systematicVariation, "muon"    ));
    files_[kBkgMca]=TFile::Open(inputFolder+"/"+TopFilename(kSigMca, systematicVariation, "electron"));
    files_[kSigPow]=TFile::Open(inputFolder+"/"+TopFilename(kSigPow, systematicVariation, "muon"    ));
    files_[kBkgPow]=TFile::Open(inputFolder+"/"+TopFilename(kSigPow, systematicVariation, "electron"));
    // create other theory prediction plots
    unsigned int koneMca= (decayChannel!="electron" ? kSigMca : kBkgMca);
    unsigned int ktwoMca= (decayChannel=="combined" ? kBkgMca : 999);
    unsigned int konePow= (decayChannel!="electron" ? kSigPow : kBkgPow);
    unsigned int ktwoPow= (decayChannel=="combined" ? kBkgPow : 999);
    // collect interesting variables
    std::vector<TString> variables_;
    TString effVariables[] = {"topPt" , "topY"  , "ttbarPt", "ttbarY", "ttbarMass", "lepPt", "lepEta", "bqPt", "bqEta", "leadqPt", "leadqEta"};
    variables_.insert(variables_.end(), effVariables, effVariables + sizeof(effVariables)/sizeof(TString));
    for(int var=0; var<int(variables_.size()); ++var){
      TString name=variables_[var];
      TString genall="analyzeTopPartonLevelKinematics/"+name;
      TString genPS ="analyzeTopPartonLevelKinematicsPhaseSpace/"+name;
      TString reco  ="analyzeTopRecoKinematicsKinFit"+addSel+"/"+name;
      // std::cout << "a)" << genall << std::endl;
      // std::cout << "b)" << genPS  << std::endl;
      // std::cout << "c)" << reco   << std::endl;
      TString ext="BGfree";
      // get ttbar SG MadGraph MC
      //std::cout << "get MadGraph" << std::endl;
      //std::cout << "a)" << std::endl;
      histo_[genall+ext][kSig]=(TH1F*)(filesEl_[kSig]->Get(genall)->Clone(genall+ext));
      histo_[genall+ext][kSig]->Sumw2();
      histo_[genall+ext][kSig]->Add((TH1F*)filesMu_[kSig]->Get(genall)->Clone(genall+ext));
      //histo_[genall+ext ][kSig]=(TH1F*)histo_[genall ][kSig]->Clone(genall+ext);
      //std::cout << "b)" << std::endl;
      histo_[genPS+ext][kSig]=(TH1F*)(filesEl_[kSig]->Get(genPS)->Clone(genPS+ext));
      histo_[genPS+ext][kSig]->Sumw2();
      histo_[genPS+ext][kSig]->Add((TH1F*)filesMu_[kSig]->Get(genPS)->Clone(genPS+ext));
      //histo_[genPS+ext ][kSig]=(TH1F*)histo_[genPS ][kSig]->Clone(genPS+ext);
      //std::cout << "c)" << std::endl;
      histo_[reco+ext][kSig]=(TH1F*)(filesEl_[kSig]->Get(reco)->Clone(reco+ext));
      histo_[reco+ext][kSig]->Sumw2();
      histo_[reco+ext][kSig]->Add((TH1F*)filesMu_[kSig]->Get(reco)->Clone(reco+ext));
      //histo_[reco+ext  ][kSig]=(TH1F*)histo_[reco  ][kSig]->Clone(reco+ext);
      // get MC@NLO SG MC 
      //std::cout << "get MC@NLO" << std::endl;
      if(files_[kSigMca]&&files_[kBkgMca]){
	histo_[genall+ext][kSigMca]=(TH1F*)(files_[koneMca]->Get(genall)->Clone(genall+ext));
	histo_[genPS+ext ][kSigMca]=(TH1F*)(files_[koneMca]->Get(genPS )->Clone(genPS+ext ));
	histo_[reco+ext  ][kSigMca]=(TH1F*)(files_[koneMca]->Get(reco  )->Clone(reco+ext  ));
	histo_[genall+ext][kSigMca]->Sumw2();
	histo_[genPS+ext ][kSigMca]->Sumw2();
	histo_[reco+ext  ][kSigMca]->Sumw2();
	if(ktwoMca!=999){
	  // add both channels for combined decayChannel
	  histo_[genall+ext][kSigMca]->Add((TH1F*)(files_[ktwoMca]->Get(genall)->Clone(genall+ext)));
	  histo_[genPS+ext ][kSigMca]->Add((TH1F*)(files_[ktwoMca]->Get(genPS )->Clone(genPS+ext )));
	  histo_[reco+ext  ][kSigMca]->Add((TH1F*)(files_[ktwoMca]->Get(reco  )->Clone(reco+ext  )));  
	}
      }
      // get POWHEG SG MC 
      //std::cout << "get POWHEG" << std::endl;
      if(files_[kSigPow]&&files_[kBkgPow]){
	histo_[genall+ext][kSigPow]=(TH1F*)(files_[konePow]->Get(genall)->Clone(genall+ext));
	histo_[genPS+ext ][kSigPow]=(TH1F*)(files_[konePow]->Get(genPS )->Clone(genPS+ext ));
	histo_[reco+ext  ][kSigPow]=(TH1F*)(files_[konePow]->Get(reco  )->Clone(reco+ext  ));
	histo_[genall+ext][kSigPow]->Sumw2();
	histo_[genPS+ext ][kSigPow]->Sumw2();
	histo_[reco+ext  ][kSigPow]->Sumw2();
	if(ktwoPow!=999){
	  // add both channels for combined decayChannel
	  histo_[genall+ext][kSigPow]->Add((TH1F*)(files_[ktwoPow]->Get(genall)->Clone(genall+ext)));
	  histo_[genPS+ext ][kSigPow]->Add((TH1F*)(files_[ktwoPow]->Get(genPS )->Clone(genPS+ext )));
	  histo_[reco+ext  ][kSigPow]->Add((TH1F*)(files_[ktwoPow]->Get(reco  )->Clone(reco+ext  )));  
	}
      }
      // Rebinning
      // WARNING: adjust to factors listed in corresponding axislabels above
      //          to have same binning for all MCs and PSs
      double reBinFactor = ( name.Contains("Pt") ? 20 : (name.Contains("ttbarMass") ? 50 : 1 ) );
      if(reBinFactor>1){
	equalReBinTH1(reBinFactor, histo_, genall+ext, kSig);
	if(files_[kSigPow]&&files_[kBkgPow]) equalReBinTH1(reBinFactor, histo_, genall+ext, kSigPow);
	if(files_[kSigMca]&&files_[kBkgMca]) equalReBinTH1(reBinFactor, histo_, genall+ext, kSigMca);
	equalReBinTH1(reBinFactor, histo_, genPS+ext , kSig);
	if(files_[kSigPow]&&files_[kBkgPow]) equalReBinTH1(reBinFactor, histo_, genPS+ext , kSigPow);
	if(files_[kSigMca]&&files_[kBkgMca]) equalReBinTH1(reBinFactor, histo_, genPS+ext , kSigMca);
	equalReBinTH1(reBinFactor, histo_, reco+ext  , kSig);
	if(files_[kSigPow]&&files_[kBkgPow]) equalReBinTH1(reBinFactor, histo_, reco+ext  , kSigPow);
	if(files_[kSigMca]&&files_[kBkgMca]) equalReBinTH1(reBinFactor, histo_, reco+ext  , kSigMca);
      }
      // Calculate efficiency / acceptance
      TString eff="efficiencyOnly/"+name;
      TString Acc="acceptanceOnly/"+name;
      // MadGraph
      histo_[eff][kSig]= (TH1F*)histo_[reco+ext  ][kSig]->Clone(eff);
      histo_[eff][kSig]->Divide(histo_[genPS+ext ][kSig]);
      histo_[Acc][kSig]= (TH1F*)histo_[genPS+ext ][kSig]->Clone(Acc);
      histo_[Acc][kSig]->Divide(histo_[genall+ext][kSig]);  
      // MC@NLO
      histo_[eff][kSigMca]= (TH1F*)histo_[reco+ext  ][kSigMca]->Clone(eff);
      histo_[eff][kSigMca]->Divide(histo_[genPS+ext ][kSigMca]);
      histo_[Acc][kSigMca]= (TH1F*)histo_[genPS+ext ][kSigMca]->Clone(Acc);
      histo_[Acc][kSigMca]->Divide(histo_[genall+ext][kSigMca]);  
      // POWHEG
      histo_[eff][kSigPow]= (TH1F*)histo_[reco+ext  ][kSigPow]->Clone(eff);
      histo_[eff][kSigPow]->Divide(histo_[genPS+ext ][kSigPow]);
      histo_[Acc][kSigPow]= (TH1F*)histo_[genPS+ext ][kSigPow]->Clone(Acc);
      histo_[Acc][kSigPow]->Divide(histo_[genall+ext][kSigPow]);  
      // add efficiency and acceptance to plotList
      plotList_ .insert(plotList_ .end(), eff);
      plotList_ .insert(plotList_ .end(), Acc);
      // add corresponding axis labels
      TString effLabel [] = {"p_{T}^{t} [GeV]" ,"y^{t}", "p_{T}^{t#bar{t}} [GeV]", "y^{t#bar{t}}", "m^{t#bar{t}} [GeV]",
      "p_{T}^{l} [GeV]", "#eta^{l}", "p_{T}^{b} [GeV]", "#eta^{b}", "p_{T}^{lead q} [GeV]", "#eta^{lead q}"};
      axisLabel_.insert(axisLabel_.end(), effLabel[var]+";eff = N_{reco} / N_{gen}^{parton PS};0;1"   );
      axisLabel_.insert(axisLabel_.end(), effLabel[var]+";A = N_{gen}^{parton PS} / N_{gen}^{all};0;1");
      // increase 1 D plots
      N1Dplots++; 
      N1Dplots++;
      // choose style
      // MadGraph
      histo_[eff][kSig]->SetLineColor(constMadgraphColor);
      histo_[eff][kSig]->SetMarkerColor(constMadgraphColor);
      histo_[eff][kSig]->SetLineWidth(2);
      histo_[eff][kSig]->SetLineStyle(1);
      histo_[eff][kSig]->SetFillStyle(0);
      histo_[Acc][kSig]->SetLineColor(constMadgraphColor);
      histo_[Acc][kSig]->SetMarkerColor(constMadgraphColor);
      histo_[Acc][kSig]->SetLineWidth(2);
      histo_[Acc][kSig]->SetLineStyle(1);
      histo_[Acc][kSig]->SetFillStyle(0);
      // MC@NLO
      if(files_[kSigMca]&&files_[kBkgMca]){
	histo_[eff][kSigMca]->SetLineColor(constMcatnloColor);
	histo_[eff][kSigMca]->SetMarkerColor(constMcatnloColor);
	histo_[eff][kSigMca]->SetLineWidth(2);
	histo_[eff][kSigMca]->SetLineStyle(constMcatnloStyle);
	histo_[eff][kSigMca]->SetFillStyle(0);
	histo_[Acc][kSigMca]->SetLineColor(constMcatnloColor);
	histo_[Acc][kSigMca]->SetMarkerColor(constMcatnloColor);
	histo_[Acc][kSigMca]->SetLineWidth(2);
	histo_[Acc][kSigMca]->SetLineStyle(constMcatnloStyle);
	histo_[Acc][kSigMca]->SetFillStyle(0);
      }
      // POWHEG
      if(files_[kSigPow]&&files_[kBkgPow]){
	histo_[eff][kSigPow]->SetLineColor(constPowhegColor );
	histo_[eff][kSigPow]->SetMarkerColor(constPowhegColor );
	histo_[eff][kSigPow]->SetLineWidth(2);
	histo_[eff][kSigPow]->SetLineStyle(constPowhegStyle);
	histo_[eff][kSigPow]->SetFillStyle(0);
	histo_[Acc][kSigPow]->SetLineColor(constPowhegColor );
	histo_[Acc][kSigPow]->SetMarkerColor(constPowhegColor );
	histo_[Acc][kSigPow]->SetLineWidth(2);
	histo_[Acc][kSigPow]->SetLineStyle(constPowhegStyle);
	histo_[Acc][kSigPow]->SetFillStyle(0);
      }
      // create small legend
      if(var==0){
	legTheo2->AddEntry(histo_[eff][kSig], "MadGraph", "L");
	if(files_[kSigPow]&&files_[kBkgPow]) legTheo2->AddEntry(histo_[eff][kSigPow], "POWHEG"  , "L");
	if(files_[kSigMca]&&files_[kBkgMca]) legTheo2->AddEntry(histo_[eff][kSigMca], "MC@NLO"  , "L");
      }
    }
  }

  // ==========================
  // compare tt BG composition
  // ==========================
  if(compareTTsample){
    // normalization
    bool area=true;
    bool SGnorm=true;
    // names
    TString loadName="analyzeTopRecoKinematicsKinFit"+addSel+"/decayChannel";
    TString newName=loadName+"GeneratorComparison";
    // get MadGraph, MC@NLO and Powheg ttbar Bkg files 
    files_[kSigMca+50]=TFile::Open(inputFolder+"/"+TopFilename(kBkgMca, systematicVariation, "muon"    ));
    files_[kBkgMca+50]=TFile::Open(inputFolder+"/"+TopFilename(kBkgMca, systematicVariation, "electron"));
    files_[kSigPow+50]=TFile::Open(inputFolder+"/"+TopFilename(kBkgPow, systematicVariation, "muon"    ));
    files_[kBkgPow+50]=TFile::Open(inputFolder+"/"+TopFilename(kBkgPow, systematicVariation, "electron"));
    files_[kSig+50   ]=TFile::Open(inputFolder+"/"+TopFilename(kBkg   , systematicVariation, "muon"    ));
    files_[kBkg+50   ]=TFile::Open(inputFolder+"/"+TopFilename(kBkg   , systematicVariation, "electron"));   
    // get MadGraph, MC@NLO and Powheg ttbar Sig files 
    files_[kSigMca+100]=TFile::Open(inputFolder+"/"+TopFilename(kSigMca, systematicVariation, "muon"    ));
    files_[kBkgMca+100]=TFile::Open(inputFolder+"/"+TopFilename(kSigMca, systematicVariation, "electron"));
    files_[kSigPow+100]=TFile::Open(inputFolder+"/"+TopFilename(kSigPow, systematicVariation, "muon"    ));
    files_[kBkgPow+100]=TFile::Open(inputFolder+"/"+TopFilename(kSigPow, systematicVariation, "electron"));
    files_[kSig+100   ]=TFile::Open(inputFolder+"/"+TopFilename(kSig   , systematicVariation, "muon"    ));
    files_[kBkg+100   ]=TFile::Open(inputFolder+"/"+TopFilename(kSig   , systematicVariation, "electron")); 
    // choose files for correct final state
    std::string one= decayChannel!="electron" ? "muon": "electron";
    std::string two= decayChannel=="combined" ? "electron": "";
    unsigned int koneMca= (decayChannel!="electron" ? kSigMca+50 : kBkgMca+50);
    unsigned int ktwoMca= (decayChannel=="combined" ? kBkgMca+50 : 999);
    unsigned int konePow= (decayChannel!="electron" ? kSigPow+50 : kBkgPow+50);
    unsigned int ktwoPow= (decayChannel=="combined" ? kBkgPow+50 : 999);   
    unsigned int koneMad= (decayChannel!="electron" ? kSig+50 : kBkg+50);
    unsigned int ktwoMad= (decayChannel=="combined" ? kBkg+50 : 999);
    // plot
    // a) MC@NLO
    if(files_[kSigMca+50]&&files_[kBkgMca+50]){
      histo_[newName][kSigMca]=(TH1F*)(files_[koneMca]->Get(loadName)->Clone(newName));
      histo_[newName][kSigMca]->Scale(lumiweight(kBkgMca, luminosity, systematicVariation, one));
      if(ktwoMca!=999){
	TH1F* temp=(TH1F*)(files_[ktwoMca]->Get(loadName)->Clone(newName));
	temp->Scale(lumiweight(kBkgMca, luminosity, systematicVariation, two));
	histo_[newName][kSigMca]->Add((TH1F*)temp->Clone());
      }
    }
    // b) Powheg
    if(files_[kSigPow+50]&&files_[kBkgPow+50]){
      histo_[newName][kSigPow]=(TH1F*)(files_[konePow]->Get(loadName)->Clone(newName));
      histo_[newName][kSigPow]->Scale(lumiweight(kBkgPow, luminosity, systematicVariation, one));
      if(ktwoPow!=999){
	TH1F* temp=(TH1F*)(files_[ktwoPow]->Get(loadName)->Clone(newName));
	temp->Scale(lumiweight(kBkgPow, luminosity, systematicVariation, two));
	histo_[newName][kSigPow]->Add((TH1F*)temp->Clone());
      }
    }
    // c) MadGraph
    if(files_[kSig+50]&&files_[kBkg+50]){
      histo_[newName][kSig]=(TH1F*)(files_[koneMad]->Get(loadName)->Clone(newName));
      histo_[newName][kSig]->Scale(lumiweight(kBkg, luminosity, systematicVariation, one));
      if(ktwoMad!=999){
	TH1F* temp=(TH1F*)(files_[ktwoMad]->Get(loadName)->Clone(newName));
	temp->Scale(lumiweight(kBkg, luminosity, systematicVariation, two));
	histo_[newName][kSig]->Add((TH1F*)temp->Clone());
      }
    }
    // adapt style 
    // MadGraph
    histo_[newName][kSig]->SetLineColor(constMadgraphColor);
    histo_[newName][kSig]->SetMarkerColor(constMadgraphColor);
    histo_[newName][kSig]->SetLineWidth(2);
    histo_[newName][kSig]->SetLineStyle(1);
    histo_[newName][kSig]->SetFillStyle(0);
    // MC@NLO
    if(files_[kSigMca+50]&&files_[kBkgMca+50]){
      histo_[newName][kSigMca]->SetLineColor(constMcatnloColor);
      histo_[newName][kSigMca]->SetMarkerColor(constMcatnloColor);
      histo_[newName][kSigMca]->SetLineWidth(2);
      histo_[newName][kSigMca]->SetLineStyle(constMcatnloStyle);
      histo_[newName][kSigMca]->SetFillStyle(0);
    }
    // POWHEG
    if(files_[kSigPow+50]&&files_[kBkgPow+50]){
      histo_[newName][kSigPow]->SetLineColor(constPowhegColor );
      histo_[newName][kSigPow]->SetMarkerColor(constPowhegColor );
      histo_[newName][kSigPow]->SetLineWidth(2);
      histo_[newName][kSigPow]->SetLineStyle(constPowhegStyle);
      histo_[newName][kSigPow]->SetFillStyle(0);
    }
    // normalization
    // a) to unity
    if(area){
      histo_[newName+"Fraction"][kSig   ]=(TH1F*)histo_[newName][kSig   ]->Clone(newName+"Fraction");
      histo_[newName+"Fraction"][kSigPow]=(TH1F*)histo_[newName][kSigPow]->Clone(newName+"Fraction");
      histo_[newName+"Fraction"][kSigMca]=(TH1F*)histo_[newName][kSigMca]->Clone(newName+"Fraction");
      histo_[newName+"Fraction"][kSig   ]->Scale(1./(histo_[newName+"Fraction"][kSig   ]->Integral(0,histo_[newName+"Fraction"][kSig   ]->GetNbinsX()+1)));
      histo_[newName+"Fraction"][kSigPow]->Scale(1./(histo_[newName+"Fraction"][kSigPow]->Integral(0,histo_[newName+"Fraction"][kSigPow]->GetNbinsX()+1)));
      histo_[newName+"Fraction"][kSigMca]->Scale(1./(histo_[newName+"Fraction"][kSigMca]->Integral(0,histo_[newName+"Fraction"][kSigMca]->GetNbinsX()+1)));
    }
    // b) to ttbar N(ttbar SG)
    if(SGnorm){
      // a) MC@NLO
      if(files_[kSigMca+100]&&files_[kBkgMca+100]){
	TH1F* temp1=(TH1F*)(files_[koneMca+50]->Get(loadName)->Clone("temp1"));
	double NMcaSG=temp1->Integral(0,temp1->GetNbinsX()+1)*lumiweight(kSigMca, luminosity, systematicVariation, one);
	if(ktwoMca!=999){
	  TH1F* temp2=(TH1F*)(files_[ktwoMca+50]->Get(loadName)->Clone("temp2"));
	  NMcaSG+=temp2->Integral(0,temp2->GetNbinsX()+1)*lumiweight(kSigMca, luminosity, systematicVariation,two);
	}
	std::cout<< "N(ttSG)/N(ttBG)[MC@NLO]=" << NMcaSG/histo_[newName][kSigMca]->Integral(0,histo_[newName][kSigMca]->GetNbinsX()+1) << std::endl;
	histo_[newName+"WRTSG"][kSigMca]=(TH1F*)histo_[newName][kSigMca]->Clone(newName+"WRTSG");
	histo_[newName+"WRTSG"][kSigMca]->Scale(1./NMcaSG);
      }
      // b) Powheg
      if(files_[kSigPow+100]&&files_[kBkgPow+100]){
	TH1F* temp1=(TH1F*)(files_[konePow+50]->Get(loadName)->Clone("temp1"));
	double NPowSG=temp1->Integral(0,temp1->GetNbinsX()+1)*lumiweight(kSigPow, luminosity, systematicVariation, one);
	if(ktwoPow!=999){
	  TH1F* temp2=(TH1F*)(files_[ktwoPow+50]->Get(loadName)->Clone("temp2"));
	  NPowSG+=temp2->Integral(0,temp2->GetNbinsX()+1)*lumiweight(kSigPow, luminosity, systematicVariation, two);
	}
	std::cout<< "N(ttSG)/N(ttBG)[POWHEG]=" << NPowSG/histo_[newName][kSigPow]->Integral(0,histo_[newName][kSigPow]->GetNbinsX()+1) << std::endl;
	histo_[newName+"WRTSG"][kSigPow]=(TH1F*)histo_[newName][kSigPow]->Clone(newName+"WRTSG");
	histo_[newName+"WRTSG"][kSigPow]->Scale(1./NPowSG);
      }
      // c) MadGraph
      if(files_[kSig+100]&&files_[kBkg+100]){
	TH1F* temp1=(TH1F*)(files_[koneMad+50]->Get(loadName)->Clone("temp1"));
	double NMadSG=temp1->Integral(0,temp1->GetNbinsX()+1)*lumiweight(kSig, luminosity, systematicVariation, one);	
	if(ktwoMad!=999){
	  TH1F* temp2=(TH1F*)(files_[ktwoMad+50]->Get(loadName)->Clone("temp2"));
	  NMadSG+=temp2->Integral(0,temp2->GetNbinsX()+1)*lumiweight(kSig, luminosity, systematicVariation, two);
	}
	std::cout<< "N(ttSG)/N(ttBG)[MADGRA]=" << NMadSG/histo_[newName][kSig]->Integral(0,histo_[newName][kSig]->GetNbinsX()+1) << std::endl;
	histo_[newName+"WRTSG"][kSig]=(TH1F*)histo_[newName][kSig]->Clone(newName+"WRTSG");
	histo_[newName+"WRTSG"][kSig]->Scale(1./NMadSG);
      }
    }

    // new plot to plotlist
    plotList_ .insert(plotList_ .end(), newName );
    if(area  ) plotList_ .insert(plotList_ .end(), newName+"Fraction" );
    if(SGnorm) plotList_ .insert(plotList_ .end(), newName+"WRTSG"    );
    
    // add corresponding axis labels
    axisLabel_.insert(axisLabel_.end(), "t#bar{t} other final state;N(t#bar{t} other);0;1");
    if(area  ) axisLabel_.insert(axisLabel_.end(), "t#bar{t} other final state;fraction (t#bar{t} other);0;1");
    if(SGnorm) axisLabel_.insert(axisLabel_.end(), "t#bar{t} other final state;N(t#bar{t} other)/N(t#bar{t} SG);0;1");
    
    // increase 1 D plots
    N1Dplots++; 
    if(area  ) N1Dplots++;
    if(SGnorm) N1Dplots++;
  }
  
  // add shape normalized PU plots
  std::map<TString, std::vector<double> > errPU_;
  TLegend* legPU=(TLegend*)leg->Clone("legPU");
  std::vector<TString> PU_;
  PU_.push_back("PUControlDistributionsBeforeBtagging/npvertex_reweighted"     );
  PU_.push_back("PUControlDistributionsAfterBtagging/npvertex_reweighted"      );
  PU_.push_back("PUControlDistributionsBeforeBtagging/npvertex_reweighted_up"  );
  PU_.push_back("PUControlDistributionsAfterBtagging/npvertex_reweighted_up"   );
  PU_.push_back("PUControlDistributionsBeforeBtagging/npvertex_reweighted_down");
  PU_.push_back("PUControlDistributionsAfterBtagging/npvertex_reweighted_down" );
  for(unsigned int pu=0; pu<PU_.size(); ++pu){
    TString name=PU_[pu];
    int PUix =std::find(plotList_.begin(), plotList_.end(), name)-plotList_.begin();
    TString oldLabel=plotList_[PUix];
    TString newLabel=oldLabel+"Norm";  
    if(verbose>2) std::cout << oldLabel << "->" << newLabel << std::endl;
    // add to plotlist
    plotList_.insert(plotList_.begin()+PUix+1,  newLabel);
    N1Dplots++;
    // add corresponding axis labels
    TString newAxisLabel=axisLabel_[PUix];
    newAxisLabel.ReplaceAll("Events", "rel. Events");
    axisLabel_.insert(axisLabel_.begin()+PUix+1, newAxisLabel);
    // produce normalized plot
    if(verbose>2) std::cout << histo_[oldLabel][kSig]->GetNbinsX() << std::endl;
    double integral=histo_[oldLabel][kSig]->Integral(0, histo_[oldLabel][kSig]->GetNbinsX()+1);
    if(verbose>2) std::cout << "int(MC)=" << integral << std::endl;
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      histo_[newLabel][sample]=(TH1F*)(histo_[oldLabel][sample]->Clone(newLabel));
      if(sample<kData) histo_[newLabel][sample]->Scale(1./integral);
      else{
	for(int bin=0; bin<=(int)histo_[oldLabel][kData]->GetNbinsX()+1; ++bin){
	  histo_[newLabel][kData]->SetBinContent(bin, (histo_[oldLabel][kData]->GetBinContent(bin)/histo_[oldLabel][kData]->Integral(0, histo_[oldLabel][kData]->GetNbinsX()+1)));
	  histo_[newLabel][kData]->SetBinError  (bin, (histo_[oldLabel][kData]->GetBinError  (bin)/histo_[oldLabel][kData]->Integral(0, histo_[oldLabel][kData]->GetNbinsX()+1)));
	}
      }
      if(verbose>2) std::cout << "rescaled int(" << sampleLabel(sample, "combined") << ")=" <<  histo_[newLabel][sample]->Integral(0, histo_[newLabel][sample]->GetNbinsX()+1) << std::endl;
    }
    // create PU uncertainty band
    if(name.Contains("down")){
      TString namedown= name;
      TString nameup  = name.ReplaceAll("down", "up");
      TString namestd = name.ReplaceAll("_up" , ""  );
      TString nameErr = namestd+"Norm";
      if(verbose>2){
	std::cout << "namedown: " << namedown<< std::endl;
	std::cout << "nameup  : " << nameup  << std::endl;
	std::cout << "namestd : " << namestd << std::endl;
	std::cout << "nameErr : " << nameErr << std::endl;
      }
      histoErrorBandPU_[nameErr]=(TH1F*)histo_[nameErr][kSig]->Clone(namestd+"Err");
      // loop bins
      std::vector<double> temp_;
      for(int bin=1; bin<=(int)histoErrorBandPU_[nameErr]->GetNbinsX(); ++bin){
	// get error
	double uncUp=std::abs(histo_[namestd][kSig]->GetBinContent(bin)-histo_[nameup  ][kSig]->GetBinContent(bin))/std::abs(histo_[namestd][kSig]->GetBinContent(bin));
	double uncDn=std::abs(histo_[namestd][kSig]->GetBinContent(bin)-histo_[namedown][kSig]->GetBinContent(bin))/std::abs(histo_[namestd][kSig]->GetBinContent(bin));
	double uncSymm=(uncUp+uncDn)/2.;
	if(std::abs(histo_[namestd][kSig]->GetBinContent(bin)==0.)) uncSymm=0.;
	if(verbose>2) std::cout << "unc(bin" << bin << "): " << uncSymm << std::endl;
	// set errorband
	histoErrorBandPU_[nameErr]->SetBinError(bin, uncSymm*histoErrorBandPU_[nameErr]->GetBinContent(bin));
	if(verbose>2) std::cout << "bin" << bin << ": " << histoErrorBandPU_[nameErr]->GetBinContent(bin) << "+/-" << histoErrorBandPU_[nameErr]->GetBinError(bin) << std::endl;
	// keep error ratio
	temp_.push_back(uncSymm*(histo_[namestd][kData]->GetBinContent(bin)/histo_[namestd][kSig]->GetBinContent(bin)));
      }
      errPU_[nameErr]=temp_;
      if(name.Contains("Before")) legPU ->AddEntry(histoErrorBandPU_[nameErr], "pile up unc.", "F");
    }
  }
  
  // =====================
  //  Create canvas
  // =====================
  std::vector<TCanvas*> plotCanvas_;
  // a) create canvas for all plots + legends
  for(unsigned int sample=0; sample<N1Dplots+N2Dplots+Nlegends; sample++){
    char canvname[10];
    sprintf(canvname,"canv%i",sample);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    //canvasStyle(*plotCanvas_[sample]);
  }

  // =====================
  //  Do the printing
  // =====================
  //
  // a) for plots
  //	  
  double xUp=-999;
  double xDn=-999;
  double min = 0;
  int canvasNumber=0;
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    bool first=true;
    // open canvas and set title corresponding to plotname in .root file
    plotCanvas_[canvasNumber]->cd(0);
    plotCanvas_[canvasNumber]->SetTitle(getStringEntry(plotList_[plot], 2)+getStringEntry(plotList_[plot], 1));	    
    // loop samples
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      // a1) for 1D event yields, efficiency and cross section plots (existing)
      if((plot<N1Dplots)||(plot>=N1Dplots+N2Dplots)){
	// check if plot is existing
       	if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)&&histo_[plotList_[plot]][sample]){
	  if(!(plotList_[plot].Contains("shift")&&(plotList_[plot].Contains("Eta")||plotList_[plot].Contains("Phi"))&&!(plotList_[plot].Contains("Nu")))) histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	  if(verbose>0){
	    std::cout << " Plotting "    << plotList_[plot];
	    if(decayChannel!="combined") std::cout << " from sample " << sampleLabel(sample,decayChannel);
	    else std::cout << " from sample " << sampleLabel(sample, "muon") << " and " << sampleLabel(sample, "electron");
	    std::cout << " to canvas "   << canvasNumber << " ( " << plotCanvas_[canvasNumber]->GetTitle() << " )";
	    std::cout << std::endl;
	  }
	  // first plot
	  if(first){
	    xUp=histo_[plotList_[plot]][sample]->GetXaxis()->GetXmax();
	    xDn=histo_[plotList_[plot]][sample]->GetXaxis()->GetXmin();
	    // plot to take maximum from
	    TString maxFrom=plotList_[plot];
	    if(plotList_[plot].Contains("_JetKinematicsTagged/pt")) maxFrom="tightLead_3_JetKinematicsTagged/pt";
	    if(plotList_[plot].Contains("_JetKinematics/pt"))       maxFrom="tightLead_3_JetKinematics/pt";
	    if(!histo_[maxFrom].count(sample)>0||!histo_[maxFrom].count(kData)>0) maxFrom=plotList_[plot];
	    if(verbose>1) std:: cout << " Take max for histo " << plotList_[plot] << " from histo " << maxFrom << std::endl;
	    // min / max
	    double max = 1.3*histo_[maxFrom][sample]->GetMaximum();
	    // if data file exists
	    if(histo_.count(maxFrom)>0&&histo_[maxFrom].count(kData)>0&&histo_[plotList_[plot]][kData]){
	      // and has larger maximum
	      if(max < 1.3*histo_[maxFrom][kData]->GetMaximum()){
		// take this maximum
		max = 1.3*histo_[maxFrom][kData]->GetMaximum();
		if(verbose>1) std:: cout << " Take max from data! " << std::endl;
	      }
	    }
	    //if(plotList_[plot].Contains("npvertex_reweighted")&&plotList_[plot].Contains("Norm")) std::cout << "test int(" << sampleLabel(sample,"combined") << ")=" << histo_[plotList_[plot]][sample]->Integral(0, histo_[plotList_[plot]][sample]->GetNbinsX()) << std::endl;
	    min = 0;
	    // log plots
	    if(getStringEntry(axisLabel_[plot],3,";")=="1"){
	      plotCanvas_[canvasNumber]->SetLogy(1);
	      min=1;
	      max=exp(1.3*(std::log(max)-std::log(min))+std::log(min));
	      if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit")&&plotList_[plot].Contains("/prob")) min=0.1; 
	      if(plotList_[plot].Contains("_JetKinematicsTagged/pt")) min=0.1;
	      if(plotList_[plot].Contains("mHbb")) min=0.1;
	      if(plotList_[plot].Contains("BGSubNorm")){ max=5.0; min=0.0001; }
	      if(plotList_[plot].Contains("reweighted")&&plotList_[plot].Contains("Norm")){min=0.0001; max=99.;}
	    }
	    // get nicer int values if maximum is large enough
	    if(max>3) max = (double)roundToInt(max);
	    if(plotList_[plot].Contains("btagSimpleSecVtx"))max*=0.8;
	    // Set x-axis range for special plots
	    if(getStringEntry(plotList_[plot], 2)=="nHit"  ){xDn=10;xUp=30;}
	    if(getStringEntry(plotList_[plot], 2)=="chi2"  ){ 
	      if(getStringEntry(plotList_[plot], 1).Contains("analyzeTopRecoKinematicsKinFit")){xDn=0;xUp=100;}
	      else{xDn=0;xUp=10;}
	    }
	    if(getStringEntry(plotList_[plot], 2)=="dB"                                                                ){xDn=0.  ;xUp=0.02;}
	    if(plotList_[plot].Contains("relIso")&&!plotList_[plot].Contains("test")                                   ){xDn=0.  ;xUp=0.15;}
	    if(plotList_[plot].Contains("tightJetKinematics/n")||plotList_[plot].Contains("tightJetKinematicsTagged/n")){xDn=3.5 ;xUp=10.5;}
	    if(plotList_[plot].Contains("tightJetKinematicsNjets1/n")                                                  ){xDn=1.5 ;xUp=10.5;}
	    if(plotList_[plot].Contains("_JetKinematics/en")                                                           ){xDn=0   ;xUp=500; }
	    if(plotList_[plot].Contains("_JetKinematics/pt")                                                           ){xDn=0   ;xUp=300; }
	    if(plotList_[plot].Contains("analyzeMETMuon/metSumEt")                                                     ){xDn=200 ;xUp=1400;}
	    if(plotList_[plot].Contains("analyzeMETMuon/metEt")                                                        ){xDn=0   ;xUp=300; }
	    if(plotList_[plot].Contains("bottomJetKinematics")&&plotList_[plot].Contains("/n")                         ){xDn=0   ;xUp=5;   }
	    if(plotList_[plot].Contains("JetKinematicsTagged")&&plotList_[plot].Contains("pt")                         ){xDn=30  ;xUp=350; }
	    if(plotList_[plot].Contains("LeptonKinematicsTagged")&&plotList_[plot].Contains("pt")                      ){xDn=30  ;xUp=150; }
	    if(plotList_[plot].Contains("btagSimpleSecVtx")                                                            ){xDn=-1  ;xUp=7;   }
	    if(plotList_[plot].Contains("btagCombSecVtx")                                                              ){xDn=0.  ;xUp=1.;  }
	    if(plotList_[plot].Contains("lepPt" )                                                                      ){xDn=0.  ;xUp=250.;}
	    if(plotList_[plot].Contains("lepEta")                                                                      ){xDn=-2.5;xUp=2.5; }
	    if((plotList_[plot].Contains("Y")||plotList_[plot].Contains("Eta"))&&(plotList_[plot].Contains("lep")||plotList_[plot].Contains("top"))&&(plotList_[plot].Contains("Plus")||plotList_[plot].Contains("Minus"))){xDn=-2.5;xUp=2.5;}
	    if(plotList_[plot].Contains("bqPt"  )                                                                      ){xDn=0.    ;xUp=250.; }
	    if(plotList_[plot].Contains("bqEta" )                                                                      ){xDn=-2.5  ;xUp=2.5;  }
	    if(plotList_[plot].Contains("topPt" )                                                                      ){xDn=0.    ;xUp=400.; }
	    if(plotList_[plot].Contains("topY"  )                                                                      ){xDn=-2.5  ;xUp=2.5;  }
	    if(plotList_[plot].Contains("ttbarY")                                                                      ){xDn=-2.5  ;xUp=2.5;  }	
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/leadqPt")                            ){xDn=0     ;xUp=400;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetPt")                                   ){xDn=30    ;xUp=400;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetEta")                                  ){xDn=-2.4  ;xUp=2.4;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetY")                                    ){xDn=-2.4  ;xUp=2.4;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/MJJ")                                              ){xDn=0     ;xUp=800;  }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/bbbarMass")                          ){xDn=0     ;xUp=800;  }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/bbbarPt")                            ){xDn=0     ;xUp=400;  }

	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftLqPt")                                        ){xDn=-80    ;xUp=80;    }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftLqEta")                                       ){xDn=-0.02  ;xUp=0.02;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftLqPhi")                                       ){xDn=-0.02  ;xUp=0.02;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftBqPt")                                        ){xDn=-30    ;xUp=30;    }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftBqEta")                                       ){xDn=-0.01  ;xUp=0.01;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftBqPhi")                                       ){xDn=-0.005 ;xUp=0.005; }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftLepPt")                                       ){xDn=-0.4   ;xUp=0.4;   }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftLepEta")                                      ){xDn=-0.000001 ;xUp=0.000001; }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftLepPhi")                                      ){xDn=-0.0000005;xUp=0.0000005;}
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftNuPt")                                        ){xDn=-60    ;xUp=60;    }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftNuEta")                                       ){xDn=-4.    ;xUp=4.;    }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/shiftNuPhi")                                       ){xDn=-0.25  ;xUp=0.25;  }
	    if(plotList_[plot].Contains("compositedKinematicsKinFit/Njets"     )                                       ){xDn=4.5    ;xUp=9.5;   }
	    if(getStringEntry(plotList_[plot], 2).Contains("topMass")                                                  ){xDn=100.   ;xUp=500.;  }

	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFitProbSel/topMass")                               ){xDn=100.   ;xUp=250.;  }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFitProbEff/ttbarMass")                             ){xDn=300.   ;xUp=1500.; }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFitProbEff/ttbarHT")                               ){xDn=50.    ;xUp=600.;  }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFitProbEff/leadqPt")                               ){xDn=0.     ;xUp=300.;  }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFitProbEff/bbbarPt")                               ){xDn=0.     ;xUp=300.;  }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFitProbEff/bbbarMass")                             ){xDn=0.     ;xUp=700.;  }
	    if(plotList_[plot].Contains("npvertex_reweighted")                                                         ){xDn=0.     ;xUp=40.;   }
	    // FIXME: intermediate by hand fix of the vertex distribution
	    // needed because analyzeTopDiffXSec_cfg.py has removed no-PU weight when running with option reduced so that they are not included in the plot...
	    if(plotList_[plot].Contains("PUControlDistributionsAfterBtagging/npvertex_reweighted")&&!plotList_[plot].Contains("Norm")){
	      for(unsigned int testsample=kSig; testsample<=kData; ++testsample){
		histo_[plotList_[plot]][testsample ]->Scale(histo_["tightJetKinematicsTagged/n"][testsample]->Integral(0,100)/histo_[plotList_[plot]][testsample ]->Integral(xDn,xUp));
	      }
	    }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarMass")                          ){xDn=0.     ;xUp=1600.; }
	    if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/neutrinoPt")                         ){xDn=0.     ;xUp=400. ; }
	    // adjust range
	    if(xUp!=histo_[plotList_[plot]][sample]->GetXaxis()->GetXmax()||xDn!=histo_[plotList_[plot]][sample]->GetXaxis()->GetXmin()) histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(xDn,xUp-0.000000001);
	    // adjust labels if overlapping because of too many large numbers
	    //if(xUp>1600) histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(false);
	    if(xUp>1600) histo_[plotList_[plot]][sample]->GetXaxis()->SetNdivisions(506);
	    // Special y-range for paper control plots
	    if (decayChannel == "combined"){
	      bool paperPlot=false;
	      if(plotList_[plot].Contains("tightJetKinematicsTagged/n"))    {min=1.0; max=1.0E06; paperPlot=true;}
	      if(plotList_[plot].Contains("bottomJetKinematics/n"))         {min=1.0; max=3.0E06; paperPlot=true;}
	      if(plotList_[plot].Contains("tightJetKinematicsTagged/pt"))   {min=1.0; max=1.0E06; paperPlot=true;}
	      if(plotList_[plot].Contains("tightLeptonKinematicsTagged/pt")){min=0.0; max=10.0E03; paperPlot=true;}
	      if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/topPt")&&!plotList_[plot].Contains("_Y")){
		min=0; max=7.0E03; paperPlot=true;  //max=2.5E04;
		if(getStringEntry(axisLabel_[plot],2,";").Contains("Events")) max*=0.5; // splitted topPt control plots
	      }
	      if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/topY"))    {min=0; max=2.5E03; paperPlot=true;} //max=1.5E04;}
	      if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarPt")) {min=0; max=6.0E03; paperPlot=true;} //max=2.2E04;}
	      if(plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarY"))  {min=0; max=1.6E03; paperPlot=true;} //max=0.8E04;}
	      if(getStringEntry(plotList_[plot], 2)=="chi2"&&getStringEntry(plotList_[plot], 1).Contains("analyzeTopRecoKinematicsKinFit")) max=10000000; 
	      if(getStringEntry(plotList_[plot], 2)=="prob") max=50000000000; 
	      if(paperPlot){
		max*=(luminosity/5000.); // lumi dependence
		if(dataSample=="2012"){
		  max*=(234./170.); // sqrt(s) dependence
		  max*=(xSec/234.); // measured cross section dependence
		}
		if(addSel=="ProbSel")  max*=0.6;
	      }
	    }
	    // axis style
	    TString titleY=getStringEntry(axisLabel_[plot],2,";");
	    if((plotList_[plot].Contains("tightLeptonKinematicsTagged/pt")||
		plotList_[plot].Contains("tightJetKinematicsTagged/pt"   )||
		plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"))&&
	       !(plotList_[plot].Contains("Eff"))){
	      double width=histo_[plotList_[plot]][sample]->GetBinWidth(1);
	      TString argument2=" / %";
	      TString precMain="2";
	      TString precSub=getTStringFromInt(getRelevantDigits(width));
	      //TString precSub="1";
	      argument2+=precMain+"."+precSub;
	      if(plotList_[plot].Contains("pt")||(plotList_[plot].Contains("Pt"))||(plotList_[plot].Contains("Pt"))||(plotList_[plot].Contains("Mass"))||(plotList_[plot].Contains("HT"))) argument2+="f GeV";
	      else argument2+="f";
	      const char * argument=argument2.Data();
	      if(verbose>1) std::cout <<  plotList_[plot] << ":" << width << " -> " << argument << " -> " << Form(argument,width) << std::endl;
	      titleY += Form(argument,width);
	    }
	    // avoid drawing 0 at y axis
	    if(min==0.&&max>0.001) min+=0.001;
	    if(min==1.&&max>1.001) min+=0.001;
	    axesStyle(*histo_[plotList_[plot]][sample], getStringEntry(axisLabel_[plot],1,";"), titleY, min, max);
	    if(!(plotList_[plot].Contains("shift")&&(plotList_[plot].Contains("Eta")||plotList_[plot].Contains("Phi"))&&(!plotList_[plot].Contains("Nu")))) histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	    if(max<100) histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
	    else histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(false);
	    // draw histos (as stack)
	    TString drawLabel="hist X0";
	    if(plotList_[plot].Contains("ProbEff")||plotList_[plot].Contains("KinFitEff")){
	      histo_[plotList_[plot]][sample]->SetMarkerColor(histo_[plotList_[plot]][sample]->GetLineColor());
	      histo_[plotList_[plot]][sample]->SetMarkerSize(0.01);
	      histo_[plotList_[plot]][sample]->SetFillColor(histo_[plotList_[plot]][sample]->GetLineColor());
	      histo_[plotList_[plot]][sample]->SetFillStyle(1001);
	      gStyle->SetErrorX(0.5);  
	      drawLabel="e2";
	    } 
	    if(plotList_[plot].Contains("shift")&&(plotList_[plot].Contains("Eta")||plotList_[plot].Contains("Phi"))&&(!plotList_[plot].Contains("Nu"))){
	      histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(false);
	    }
	    histo_[plotList_[plot]][sample]->Draw(drawLabel);
	    histo_[plotList_[plot]][ENDOFSAMPLEENUM] = (TH1F*)(histo_[plotList_[plot]][sample]->Clone());
	    if(!(plotList_[plot].Contains("shift")&&(plotList_[plot].Contains("Eta")||plotList_[plot].Contains("Phi"))&&!(plotList_[plot].Contains("Nu")))) histo_[plotList_[plot]][ENDOFSAMPLEENUM]->GetXaxis()->SetNoExponent(true);
	    // draw other theory predictions for BG subtracted plots
	    if(plotList_[plot].Contains("BGSubNorm")){
	      if(histo_[plotList_[plot]].count(kSigMca)>0) histo_[plotList_[plot]][kSigMca]->Draw("hist X0 same");
	      if(histo_[plotList_[plot]].count(kSigPow)>0) histo_[plotList_[plot]][kSigPow]->Draw("hist X0 same");
	      // draw small legend
	      legTheo->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.2);
	      legTheo->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() +0.05 - 0.03 * leg->GetNRows());
	      legTheo->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
	      legTheo->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength());
	      if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetPt")||plotList_[plot].Contains("compositedKinematicsKinFit/Njets"))legTheo->Draw("same");
	    }
	    // draw other theory predictions for efficiency and acceptance study in same plot
	    if(plotList_[plot].Contains("efficiencyOnly")||plotList_[plot].Contains("acceptanceOnly")){
	      if(histo_[plotList_[plot]].count(kSigMca)>0) histo_[plotList_[plot]][kSigMca]->Draw("hist same");
	      if(histo_[plotList_[plot]].count(kSigPow)>0) histo_[plotList_[plot]][kSigPow]->Draw("hist same");
	      // draw small legend
	      legTheo2->SetX1NDC(0.73);
	      legTheo2->SetY1NDC(0.69);
	      legTheo2->SetX2NDC(0.91);
	      legTheo2->SetY2NDC(0.88);
	      legTheo2->Draw("same");
	      // turn on grid
	      plotCanvas_[canvasNumber]->SetGrid(0,1);
	    }
	    // draw other theory predictions for ttbar other composition in same plot
	    if(plotList_[plot].Contains("decayChannel")&&plotList_[plot].Contains("GeneratorComparison")){
	      if(histo_[plotList_[plot]].count(kSigMca)>0) histo_[plotList_[plot]][kSigMca]->Draw("hist same");
	      if(histo_[plotList_[plot]].count(kSigPow)>0) histo_[plotList_[plot]][kSigPow]->Draw("hist same");
	      // turn on grid
	      plotCanvas_[canvasNumber]->SetGrid(0,1);
	      // theory ratios wrt MadGraph
	      if(withRatioPlot){
		std::vector<double> zeroerr_;
		for(int bin=1; bin<=histo_[plotList_[plot]][kSig]->GetNbinsX(); ++bin) zeroerr_.push_back(0);
		int rval1 = drawRatio(histo_[plotList_[plot]][kSigMca], histo_[plotList_[plot]][kSig], 0.5, 1.19, myStyle, verbose, zeroerr_, "simulation", "MadGraph", "hist"     , constMcatnloColor, false, 0.1);
		int rval2 = drawRatio(histo_[plotList_[plot]][kSigPow], histo_[plotList_[plot]][kSig], 0.5, 1.19, myStyle, verbose, zeroerr_, "simulation", "MadGraph", "hist same", constPowhegColor , false, 0.1);
		if (rval1!=0||rval2!=0) std::cout << " Problem occured when creating ratio plot for " << plotList_[plot] << std::endl;
	      }
	    }
	  }
	  // draw other plots into same canvas 
	  else{ 
	    // draw data as points
	    if(sample==kData) histo_[plotList_[plot]][sample]->Draw("p e X0 same");
	    // draw others as histo (stack)
	    else histo_[plotList_[plot]][sample]->Draw("hist X0 same");
	  }
	  first=false;
	  // at the end:
	  if((sample==kData)&&(histo_.count(plotList_[plot])>0)&&histo_[plotList_[plot]][sample]){
	    if(histoErrorBandPU_.count(plotList_[plot])>0){
	      if(verbose>2) std::cout << "draw PU unc. band for " << plotList_[plot] << std::endl;
	      // configure style of and draw uncertainty bands
	      histoErrorBandPU_[plotList_[plot]]->SetMarkerStyle(0);
	      histoErrorBandPU_[plotList_[plot]]->SetFillColor(1);
	      histoErrorBandPU_[plotList_[plot]]->SetFillStyle(3004);
	      gStyle->SetErrorX(0.5);  
	      histoErrorBandPU_[plotList_[plot]]->Draw("E2 SAME");
	    }
	    else if(histoErrorBand_.count(plotList_[plot])>0){
	      // configure style of and draw uncertainty bands
	      histoErrorBand_[plotList_[plot]]->SetMarkerStyle(0);
	      histoErrorBand_[plotList_[plot]]->SetFillColor(1);
	      histoErrorBand_[plotList_[plot]]->SetFillStyle(3004);
	      gStyle->SetErrorX(0.5);  
	      histoErrorBand_[plotList_[plot]]->Draw("E2 SAME");
	    }
	    // redraw axis 
	    histo_[plotList_[plot]][ENDOFSAMPLEENUM]->Draw("axis X0 same");
	    // draw label indicating the analysis cuts applied
	    if((unsigned int)canvasNumber<plotCanvas_.size()-Nlegends){
	      if(withRatioPlot&&(!(plotList_[plot].Contains("bottomJetKinematics/n")||plotList_[plot].Contains("tightJetKinematicsTagged/n")||plotList_[plot].Contains("tightJetKinematicsTagged/pt")||plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/topPt")||plotList_[plot].Contains("tightLeptonKinematicsTagged/pt")||plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/topY")||plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarY")||plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarPt")||plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit"+addSel+"/ttbarMass")||plotList_[plot].Contains("analyzeTopRecoKinematicsKinFit/prob")||plotList_[plot].Contains("ProbEff")||plotList_[plot].Contains("KinFitEff")))){
		// draw cut label
		TString cutLabel="1 lepton, #geq4 Jets";
		if(decayChannel=="muon"    ) cutLabel.ReplaceAll("lepton","#mu");
		if(decayChannel=="electron") cutLabel.ReplaceAll("lepton","e");
		if(decayChannel=="combined") cutLabel.ReplaceAll("lepton","e/#mu");
		if(plotList_[plot].Contains("Tagged")||plotList_[plot].Contains("AfterBtagging")) cutLabel+=", #geq2 b-tags";
		else if(plotList_[plot].Contains("KinFit")) cutLabel+=", #geq2 b-tags, KinFit";
		if(plotList_[plot].Contains("ProbSel")) cutLabel.ReplaceAll("KinFit","prob>0.02");
		double positionX=xUp+0.045*(xUp-xDn)*(gStyle->GetCanvasDefW()/600.);
		double positionY=min;
		//std::cout << plotList_[plot] << ":" << positionY << std::endl;
		//std::cout << plotList_[plot] << ": " << xUp << "+0.03*(" << xUp << "-" << xDn << ")=" << positionX << std::endl;
		TLatex *sellabel = new TLatex(positionX,positionY,cutLabel);
		sellabel->SetTextAlign(11);
		sellabel->SetTextAngle(90);
		sellabel->SetTextSize(0.04);
		sellabel->Draw("same");
	      }
	      double x1=1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25;
	      if(systematicVariation==sysTest) x1-=0.05;
	      leg->SetX1NDC(x1+0.015);
	      leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() -0.05 - 0.03 * leg->GetNRows());
	      leg->SetX2NDC(1.03- gStyle->GetPadRightMargin() - gStyle->GetTickLength());
	      leg->SetY2NDC(1.0-gStyle->GetPadTopMargin()-0.8*gStyle->GetTickLength());
	      leg->SetTextSize(0.035);

	      legPU->SetX1NDC(x1+0.015);
	      legPU->SetY1NDC(1.0 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() -0.05 - 0.03 * legPU->GetNRows());
	      legPU->SetX2NDC(1.03- gStyle->GetPadRightMargin() - gStyle->GetTickLength());
	      legPU->SetY2NDC(1.0-gStyle->GetPadTopMargin()-0.8*gStyle->GetTickLength());
	      legPU->SetTextSize(0.035);

	      if(!plotList_[plot].Contains("KinFitEff")&&!plotList_[plot].Contains("ProbEff")&&!plotList_[plot].Contains("BGSubNorm")&&!(plotList_[plot].Contains("efficiency")||plotList_[plot].Contains("acceptance"))&&!(plotList_[plot].Contains("efficiency")||plotList_[plot].Contains("acceptance"))&&!(plotList_[plot].Contains("reweightedNorm"))) leg->Draw("SAME");
	      else if(plotList_[plot].Contains("reweightedNorm")) legPU->Draw("SAME");
	      // add labels for decay channel, luminosity, energy and CMS preliminary (if applicable)
	      if      (decayChannel=="muon"    ) DrawDecayChLabel("#mu + Jets");
	      else if (decayChannel=="electron") DrawDecayChLabel("e + Jets");
	      else DrawDecayChLabel("e/#mu + Jets Combined");	      
	      DrawCMSLabels(true,luminosity); 
	      // draw data/MC ratio
	      if((histo_[plotList_[plot]].count(kSig)>0) && withRatioPlot && !plotList_[plot].Contains("PartonLevel")){
		if(plotList_[plot].Contains("BGSubNorm")){
		  std::vector<TH1F*>ratiohists_;
		  // get theory
		  ratiohists_.push_back( histo_[plotList_[plot]][kSig   ] );
		  if(histo_[plotList_[plot]].count(kSigMca)>0) ratiohists_.push_back( histo_[plotList_[plot]][kSigMca] );
		  if(histo_[plotList_[plot]].count(kSigPow)>0) ratiohists_.push_back( histo_[plotList_[plot]][kSigPow] );
		  double theomax=1.8;
		  double theomin=0.2;
		  double theoxmax=-1;
		  double theoxmin=-1;
		  if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetPt" )){ theoxmin=30. ; theoxmax=400.; }
		  if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetEta")){ theoxmin=-2.4; theoxmax=2.4 ; }
		  if(plotList_[plot].Contains("compositedKinematicsKinFit/leadNonttjetY"  )){ theoxmin=-2.4; theoxmax=2.4 ; }
		  if(plotList_[plot].Contains("compositedKinematicsKinFit/Njets"          )){ theoxmin=4   ; theoxmax=9   ; }
		  // create ratio canvas
		  std::vector<TCanvas*> plotCanvasRatio_;
		  plotCanvasRatio_.push_back(drawFinalResultRatio(histo_[plotList_[plot]][kData], theomin, theomax, myStyle, 0, ratiohists_, plotCanvas_[canvasNumber], theoxmin, theoxmax));
		  plotCanvasRatio_[0]->Draw();
		  plotCanvasRatio_[0]->Update();
		}
		else{
		  double ratMin=0.1;
		  double ratMax=1.9;
		  int ndivisions=505;
		  if(plotList_[plot].Contains("shift")){ratMin=0.88; ratMax=1.12;}
		  if(plotList_[plot].Contains("TopRecoKinematics")){ratMin=0.49; ratMax=1.49; ndivisions=405;}
		  if(getStringEntry(plotList_[plot], 2)=="prob"){ratMin=0.75; ratMax=1.25;}
		  //if(plotList_[plot].Contains("topMass")){ratMin=0.7; ratMax=2.0;}
		  TString ratioLabelNominator  ="N_{data}";
		  TString ratioLabelDenominator="N_{MC}";
		  std::vector<double> err_=std::vector<double>(0);
		  std::vector<double> err2_;
		  if(plotList_[plot].Contains("reweightedNorm")) err_=errPU_[plotList_[plot]];
		  else if(plotList_[plot].Contains("KinFitEff")||plotList_[plot].Contains("ProbEff")){legProbEff->Draw("same"); ratMin=0.79; ratMax=1.21; ratioLabelNominator  ="eff_{data}"; ratioLabelDenominator="eff_{MC}"; if(errEff_.count(plotList_[plot])>0){err_=errEff_[plotList_[plot]];}}
		  
		  int rval =0;
		  if(!plotList_[plot].Contains("PU")||plotList_[plot].Contains("reweightedNorm")){
		    rval=drawRatio(histo_[plotList_[plot]][kData], histo_[plotList_[plot]][kSig], ratMin, ratMax, myStyle, verbose, err_, ratioLabelNominator, ratioLabelDenominator, "p e X0", kBlack, true, 0.8, ndivisions);
		    if (rval!=0) std::cout << " Problem occured when creating ratio plot for " << plotList_[plot] << std::endl;
		  }
		  // PU shape errorband
		  if(histoErrorBandPU_.count(plotList_[plot])>0&&rval==0){
		    if(verbose>2) std::cout << "draw PU unc. band in ratio for " << plotList_[plot] << std::endl;
		    for(int bin=1; bin<=histoErrorBandPU_[plotList_[plot]]->GetNbinsX(); ++bin){
		      err2_.push_back(histoErrorBandPU_[plotList_[plot]]->GetBinError(bin)/histoErrorBandPU_[plotList_[plot]]->GetBinContent(bin));
		    }
		    int rval2 = drawRatio(histoErrorBandPU_[plotList_[plot]], histoErrorBandPU_[plotList_[plot]], ratMin, ratMax, myStyle, verbose, err2_, ratioLabelNominator, ratioLabelDenominator, "e2 same", kBlack, true, 0.8, ndivisions);
		    if (rval2!=0) std::cout << " Problem occured when creating errorband ratio plot for " << plotList_[plot] << std::endl;		    
		  }
		  // MC shape arrorband
		  else if(histoErrorBand_.count(plotList_[plot])>0&&rval==0){
		    for(int bin=1; bin<=histoErrorBand_[plotList_[plot]]->GetNbinsX(); ++bin){
		      err2_.push_back(histoErrorBand_[plotList_[plot]]->GetBinError(bin)/histoErrorBand_[plotList_[plot]]->GetBinContent(bin));
		    }
		    int rval2 = drawRatio(histoErrorBand_[plotList_[plot]], histoErrorBand_[plotList_[plot]], ratMin, ratMax, myStyle, verbose, err2_, ratioLabelNominator, ratioLabelDenominator, "e2 same", kBlack, true, 0.8, ndivisions);
		    if (rval2!=0) std::cout << " Problem occured when creating errorband ratio plot for " << plotList_[plot] << std::endl;		    
		  }
		}
	      }
	    }
	  }
	  // efficiency and acceptance ratios: theory ratios wrt MadGraph
	  if(withRatioPlot&&(plotList_[plot].Contains("efficiency")||plotList_[plot].Contains("acceptance"))){
	    std::vector<double> zeroerr_;
	    for(int bin=0; bin<histo_[plotList_[plot]][kSig]->GetNbinsX(); ++bin) zeroerr_.push_back(0);
	    int rval1 = drawRatio(histo_[plotList_[plot]][kSigMca], histo_[plotList_[plot]][kSig], 0.7, 1.3, myStyle, verbose, zeroerr_, "simulation", "MadGraph", "hist"     , constMcatnloColor, false, 0.1);
	    int rval2 = drawRatio(histo_[plotList_[plot]][kSigPow], histo_[plotList_[plot]][kSig], 0.7, 1.3, myStyle, verbose, zeroerr_, "simulation", "MadGraph", "hist same", constPowhegColor , false, 0.1);
	    if (rval1!=0||rval2!=0) std::cout << " Problem occured when creating ratio plot for " << plotList_[plot] << std::endl;
	  }
      	} // end: plot existing?
      } // end: 1D plots
      
      //  a2) for 2D plots
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)){
	// new Canvas for every plot
	plotCanvas_[canvasNumber]->cd(0);
	plotCanvas_[canvasNumber]->SetRightMargin(0.15);
	if(plotList_[plot].Contains("shift")&&(plotList_[plot].Contains("Eta")||plotList_[plot].Contains("Phi"))&&!(plotList_[plot].Contains("Nu"))) plotCanvas_[canvasNumber]->SetRightMargin(0.30);
	plotCanvas_[canvasNumber]->SetTitle(getStringEntry(plotList_[plot],2)+getStringEntry(plotList_[plot],1)+getTStringFromInt(sample)); 
	if(verbose>1){
	  std::cout << " Plotting "    << plotList_[plot];
	  std::cout << " from sample " << sampleLabel(sample,decayChannel);
	  std::cout << " to canvas "   << canvasNumber  << " ( " << plotCanvas_[canvasNumber]->GetTitle() << " )";
	  std::cout << std::endl;
	}
	++canvasNumber;
	// draw histo
	histo2_[plotList_[plot]][sample]->Draw("colz");
	// print correlation factor
	DrawLabel(Form("Correlation: %f",histo2_[plotList_[plot]][sample]->GetCorrelationFactor()), 0.35, 0.92, 0.75, 0.99, 0.7);
      }
    }
    // for 1D hists: next canvas
    if((plot<N1Dplots)&&(histo_.count(plotList_[plot])>0)) ++canvasNumber;
  }

  // b) for legends
  //
  // b1) before btagging
  //  
  plotCanvas_[canvasNumber]->cd(0);
  plotCanvas_[canvasNumber]->SetTitle("legendMonitoringBeforeBtagging");
  leg0->Draw("");
  ++canvasNumber;
  //
  // b2) after btagging
  //
  plotCanvas_[canvasNumber]->cd(0);
  plotCanvas_[canvasNumber]->SetTitle("legendMonitoringAfterBtagging");
  leg1->Draw("");
  ++canvasNumber;
	
  // ==============
  //  Saving plots
  // ==============
  //std::cout << "test0" << std::endl;
  if(save){
    // eps and png
    if(verbose==0) gErrorIgnoreLevel=kWarning;
    saveCanvas(plotCanvas_, outputFolder, pdfName, true, false);
    unsigned int idxstop=plotCanvas_.size();
    for(unsigned int idx=0; idx<plotCanvas_.size(); idx++){
      TString title=(TString)(plotCanvas_[idx])->GetTitle();
      if(!title.Contains("canv")){
	TString outputFolder2=outputFolder;
	if(title.Contains("efficiencyOnly")||title.Contains("acceptanceOnly")){
	  outputFolder2.ReplaceAll("monitoring", "effAndAcc");
	}
	plotCanvas_[idx]->Print(outputFolder2+title+".eps"); 
	plotCanvas_[idx]->Print(outputFolder2+title+".png");
      }
      // kinfitshift as pdf
      if(title.Contains("shift")){
	TString pdfname=outputFolder+"KinFitObjectShift.pdf";
	if(title.Contains("shiftLqPt" )) pdfname+="(";
	if(title.Contains("shiftNuPhi")) pdfname+=")";
	plotCanvas_[idx]->Print(pdfname); 
      }
      // prob efficiency as pdf
      if(title.Contains("ProbEff")){
	TString pdfname=outputFolder+"ProbEff.pdf";
	if(title.Contains("ttbarAngle")) pdfname+="(";
	if(title.Contains("bbbarMass" )){ pdfname+=")"; idxstop=idx; }
	if(idx<=idxstop) plotCanvas_[idx]->Print(pdfname); 
      }
    }
    // root file
    if(verbose>0) std::cout << "will save into outputfile named " << outputFileName << std::endl;
    for(unsigned int idx=0; idx<plotCanvas_.size(); ++idx){
      saveToRootFile(outputFileName, plotCanvas_[idx], true, verbose, "monitoring");
    }
  }
  //std::cout << "test" << std::endl;
  // delete pointer
  //closeStdTopAnalysisFiles(files_);
  //std::cout << "done" << std::endl;
}
