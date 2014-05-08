#include "basicFunctions.h"
#include "../../unfolding/TopSVDFunctions.h" 
#include "../../unfolding/TopSVDFunctions.C" 

void analyzeRegularizationTest(TString test="mix", double luminosity = 19712.,
			       bool save = true, int systematicVariation=sysNo, unsigned int verbose=0,
			       TString inputFolderName="RecentAnalysisRun8TeV_doubleKinFit",
			       TString dataFile= groupSpace+AnalysisFolder+"/elecDiffXSecData2012ABCDAll.root:"+groupSpace+AnalysisFolder+"/muonDiffXSecData2012ABCDAll.root",
			       std::string decayChannel = "combined", bool SVDunfold=true, bool extrapolate=true, bool hadron=false, 
			       bool redetermineopttau =false, TString closureTestSpecifier="", TString addSel="ProbSel")
{
  // iterative unfolding?
  bool iterative=false;
  // regularization mode and scanned parameters
  // test : mix, k, tau
  std::vector <double > tau_;
  std::vector <int >    tauColor_;
  std::vector <int >    tauMarker_;
  std::vector <TH1F*>   tauxSecNorm_;
  std::vector <TGraphAsymmErrors*> tauxSecNormM_;
  std::map < TString, double> stdTau_;
  int kindex;
  // ============================  
  //  Set ROOT Style
  // ============================
 
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gROOT->ForceStyle();

  // ============================
  //  Name Conventions
  // ============================
  // a) enumerator sample convention (as defined in basicFunctions.h)
  //    kSig,    kBkg,    kZjets,  kWjets, 
  //    kQCD,    kSTop,   kDiBos,  kData, 
  //    kQCDEM1, kQCDEM2, kQCDEM3, kQCDBCE1, kQCDBCE2, kQCDBCE3,  
  //    kWW,     kWZ,     kZZ, 
  //    kSTops,  kSATops, kSTopt,  kSATopt,  kSToptW , kSAToptW};

  // b) file name convention (implemented in basicFunctions.h)

  // decayChannel+"DiffXSec"+sampleName+systematicVariation+MCProduction+"PF.root"

  // decayChannel      = "elec", "muon"
  // sampleName        = "Sig", "Bkg", Wjets", "Zjets", "WW", "WZ", "ZZ", "VV", "SingleTopSchannel", 
  //                     "SingleTopTchannel", "SingleTopTWchannel", "QCD"
  // MCProductionCycle = "Summer11","Fall11"

  // ============================
  //  Options
  // ============================
  //
  // a) options directly entered when calling function
  //    systematicVariation: which systematic shift do you want to make? from basicFunctions.h:  
  //
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
  //        35: sysQCDUp                   36: sysQCDDown                  
  //        37: sysSTopUp                  38: sysSTopDown                 
  //        39: sysDiBosUp                 40: sysDiBosDown                
  //        41: sysPDFUp                   42: sysPDFDown                  
  //        43: sysHadUp                   44: sysHadDown                  
  //        45: sysGenMCatNLO              46: sysGenPowheg  
  //        47: sysGenPowhegHerwig         48: ENDOFSYSENUM
  
  // systematic variation for which you want to print the
  // inclusive cross section and the jet permutation overview
  // testMe=ENDOFSYSENUM means none
  int testMe=ENDOFSYSENUM;//ENDOFSYSENUM
  // name of the cross section quantity you want to print
  // detailed information for
  // testQuantity="" means all, testQuantity="NOTEST" means none
  TString testQuantity="NOTEST";
  // closureTestSpecifier = \"NoDistort\", \"topPtUp\", \"topPtDown\", \"ttbarMassUp\", \"ttbarMassDown\", \"data\" or \"1000\"
  // errorbands for yield plots
  bool errorbands=false;
  // addSel: xSec from prob selection step
  //TString addSel="ProbSel";
  //if(!extrapolate) addSel="";
  // for systematic variations done in central signal sample (e.g. PU, effSF, btagSF):
  // use MC rootfiles from subfolder /Prob and ignore addSel plotname extension
  TString addSelData="";
  if(addSel=="ProbSel"&&((systematicVariation>=sysPUUp&&systematicVariation<=sysPUDown)||(systematicVariation>=sysLepEffSFNormUp&&systematicVariation<=sysMisTagSFDown))){
    inputFolderName+="/Prob";
    addSelData=addSel;
  }
  // luminosity uncertainties
  if(systematicVariation==sysLumiUp  )      luminosity*=1.044;
  else if(systematicVariation==sysLumiDown) luminosity*=0.956;
  // verbose: set detail level of output 
  // 0: no output, 1: std output 2: more output 3: output for debugging
  // data file: relative path of .root file
  // save: save plots?
  // SVDunfold: use SVD instead of bin to bin unfolding
  // choose phase space
  // -> process only visible hadron level PS and extrapolated parton level PS 
  if(!extrapolate) hadron=true;
  else hadron=false;
  TString PS="";
  // a) for full PS use extrapolate=true;
  if(!extrapolate) PS="PhaseSpace";
  // b) for restricted phase space:
  // b1) parton PS: hadron=false
  // b2) hadron PS: hadron=true
  TString LV="Parton";
  if(!extrapolate&&hadron) LV="Hadron";
  if(verbose>1){
    if(extrapolate) std::cout << "full Phase Space will be used!" << std::endl; 
    else std::cout << LV << " level Phase Space will be used!" << std::endl; 
  }
  // compare SVD and BBB results
  bool compare=false;
  if(!SVDunfold) compare=false;
  // luminosity: [/pb]
  TString lumi = getTStringFromInt(roundToInt((luminosity), false));  
  // b) options to be configured only once
  // get the .root files from the following folder:
  TString inputFolder = groupSpace+inputFolderName;
  // see if its 2010 or 2011 data from luminosity
  TString dataSample="2012";
  // for closure test if desired
  TString closureLabel = "";
  TString zprime="";
  double  zPrimeLumiWeightIni=1.;
  TString zPrimeLumiWeightStr="";
  if(closureTestSpecifier!=""){
    if      (closureTestSpecifier.Contains("NoDistort")) closureLabel = "PseudoData"+closureTestSpecifier;
    else if (closureTestSpecifier.Contains("topPt"   )||
	     closureTestSpecifier.Contains("ttbarMass")) closureLabel = "PseudoDataReweight"+closureTestSpecifier;
    else if (closureTestSpecifier.Contains("data"     )) closureLabel = "PseudoDataReweighttopPt"+closureTestSpecifier;
    else if (closureTestSpecifier.Contains("1000"     )){
      // identify exact closure test
      zprime="1000";
      if      (closureTestSpecifier.Contains("x0p03")) {zPrimeLumiWeightIni=0.03; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.03)";}
      else if (closureTestSpecifier.Contains("x0p1" )) {zPrimeLumiWeightIni=0.1 ; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.1)" ;}
      else if (closureTestSpecifier.Contains("x0p25")) {zPrimeLumiWeightIni=0.25; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.25)";}
      else if (closureTestSpecifier.Contains("x0p5" )) {zPrimeLumiWeightIni=0.5 ; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.5)" ;}
      else if (closureTestSpecifier.Contains("x2"   )) {zPrimeLumiWeightIni=2.  ; zPrimeLumiWeightStr=" (#sigma_{Z'}x2)"   ;}
      else if (closureTestSpecifier.Contains("x4"   )) {zPrimeLumiWeightIni=4.  ; zPrimeLumiWeightStr=" (#sigma_{Z'}x4)"   ;}
      closureLabel = "PseudoDataZprime"+closureTestSpecifier+"GeV";
    }
    else{
      std::cout << "ERROR: unknown closureTestSpecifier=" << closureTestSpecifier << std::endl;
      exit(0);
    }
  }
  if(closureTestSpecifier!=""&&!dataFile.Contains("PseudoData")){
    dataFile=inputFolder+"/pseudodata/"+pseudoDataFileName(closureTestSpecifier, "electron")+":"+inputFolder+"/pseudodata/"+pseudoDataFileName(closureTestSpecifier, "muon");
    // closure test pseudo data is created from ttbar simulation without madspin
    dataFile.ReplaceAll("MadSpin", "");
  }
  // save all plots into the following foldre
  TString outputFolder = "./diffXSecFromSignal/plots/"+decayChannel+"/";
  if(dataSample!="") outputFolder+=dataSample+"/";
  // save all plots within a root file named:
  TString outputFileName="diffXSecTopSemi";
  if(decayChannel=="muon"    ) outputFileName+="Mu";
  if(decayChannel=="electron") outputFileName+="Elec";
  if(decayChannel=="combined") outputFileName+="Lep";
  outputFileName+="RegularizationTest";
  outputFileName+=closureLabel+LV+PS+".root";
  // choose name of the output .pdf file
  TString pdfName="kinFitHypothesis"+lumi+"pb";
  // choose if you want to set QCD artificially to 0 to avoid problems with large SF for single events
  bool setQCDtoZero=true;
  if(setQCDtoZero&&verbose>1) std::cout << "ATTENTION: qcd will artificially be set to 0 (for everything BUT the inclusive cross section)!"; 
  // redetermine optimal tau
  if(!SVDunfold || systematicVariation!=sysNo) redetermineopttau =false; // never determine tau for syst. var or if not SVD unf. is used
  if(redetermineopttau){
    if(verbose>1) std::cout << "ATTENTION: optimal tau for SVD unfolding will be determined! this takes a while"; 
    save=false;
  }
  // adjust luminosity and data files for combined control plots
  double luminosityEl=0;
  double luminosityMu=0;
  TString dataFileEl="";
  TString dataFileMu="";
  if(decayChannel=="combined"){
    luminosityEl=constLumiElec;
    luminosityMu=constLumiMuon;
    if(systematicVariation==sysLumiUp ){      luminosityEl*=(1.0+globalLumiUncertainty); luminosityMu*=(1.0+globalLumiUncertainty);}
    else if(systematicVariation==sysLumiDown){luminosityEl*=(1.0-globalLumiUncertainty); luminosityMu*=(1.0-globalLumiUncertainty);}
    if(!dataFile.Contains(":")){
      std::cout << "wrong input filenames, should be dataFileEl:dataFileMu, but is";
      std::cout << dataFile << std::endl;
      exit(0);
    }
    dataFileEl=getStringEntry(dataFile,1 , ":");
    dataFileMu=getStringEntry(dataFile,42, ":");
  }
  double luminosity2=luminosity;
  if(decayChannel=="combined") luminosity2= ( constLumiElec + constLumiMuon );
  if(systematicVariation==sysLumiUp  )      luminosity2*=(1.0+globalLumiUncertainty);
  else if(systematicVariation==sysLumiDown) luminosity2*=(1.0-globalLumiUncertainty);
  // output
  if(verbose>1||closureTestSpecifier!="") {
    std::cout << "dataFile " << dataFile << std::endl;
    std::cout << "closureTestSpecifier " << closureTestSpecifier << std::endl;
    std::cout << "outputFileName " << outputFileName << std::endl;
  }

  // ---
  //    create list of systematics to be ignored- the std input will be used instead
  // ---
  int systematicVariationMod=systematicVariation;
  std::vector<int> ignoreSys_;
  //ignoreSys_.push_back(sysGenMCatNLO);
  // exclude JES and JER
  //for(int sys=sysJESUp     ; sys<=sysJERDown    ; ++sys) ignoreSys_.push_back(sys);
  // exclude Scale matching and top mass 
  //for(int sys=sysTopScaleUp; sys<=sysTopMassDown; ++sys) ignoreSys_.push_back(sys);
  // exclude Hadronization
  for(int sys=sysHadUp     ; sys<=sysHadDown    ; ++sys) ignoreSys_.push_back(sys);
  // exclude PDF
  //for(int sys=sysPDFUp     ; sys<=sysPDFDown    ; ++sys) ignoreSys_.push_back(sys);
  // FIXME: exclude all
  //for(int sys=3     ; sys<=47    ; ++sys) ignoreSys_.push_back(sys);
  // use std variable for loading plots in case of listed systematics
  for(unsigned int i=0; i<ignoreSys_.size(); ++i){
    if(systematicVariation==ignoreSys_[i]) systematicVariationMod=sysNo;
  }
  // use different ttbar MC ("Madgraph", "Powheg", "PowhegHerwig", "McatNLO"), also used for generator uncertainties
  TString ttbarMC="Madgraph";
  if(systematicVariationMod==sysGenMCatNLO          ) ttbarMC="Mcatnlo";
  else if(systematicVariationMod==sysGenPowheg      ) ttbarMC="Powheg";
  else if(systematicVariationMod==sysGenPowhegHerwig) ttbarMC="PowhegHerwig";
  TString ttbarMC2=ttbarMC;  
  if(systematicVariation==sysHadUp||systematicVariation==sysHadDown) ttbarMC2="Powheg"; // this is the old hadronization uncertainty will be overwritten from 'bothDecayChannelCombination'
  // normalization of differential normalized cross sections
  // true: use integral of unfolded differential cross section ignoring UF/OF bins
  // false: use all events and inclusive eff*A (no BR!) before unfolding
  //bool normToIntegral=true;
  //if(extrapolate)  normToIntegral=false;
  // choose plot input folder corresponding to systematic Variation  
  TString sysInputFolderExtension="";
  TString sysInputGenFolderExtension="";

  // create list of variables you would like to create the efficiency / cross section for (centrally defined in basicFunctions.h)
  std::vector<TString> xSecVariables_, xSecLabel_; 
  // a) top and ttbar quantities
  TString xSecVariablesKinFit2[] = {"topPt"};
  TString xSecLabelKinFit2[]     = {"p_{T}^{t}/[GeV]"};
  if(!hadron){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFit2, xSecVariablesKinFit2 + sizeof(xSecVariablesKinFit2)/sizeof(TString));
    xSecLabel_    .insert(xSecLabel_.end(),     xSecLabelKinFit2,     xSecLabelKinFit2     + sizeof(xSecLabelKinFit2    )/sizeof(TString));
  }
  // b) lepton and b-jet quantities
  if(hadron||!extrapolate){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalState, xSecVariablesFinalState + sizeof(xSecVariablesFinalState)/sizeof(TString));
    xSecLabel_    .insert(xSecLabel_.end()    , xSecLabelFinalState    , xSecLabelFinalState     + sizeof(xSecLabelFinalState    )/sizeof(TString));
  }
  // check
  if(xSecVariables_.size()!=xSecLabel_.size()){
    std::cout << "ERROR: size of xSecVariables_ and xSecLabel_ is unequal!" << std::endl;
    std::cout << "xSecVariables_" << std::endl;
    for(unsigned int ix=0; ix<xSecVariables_.size(); ++ix){
     std::cout << xSecVariables_[ix] << std::endl; 
    }
    std::cout << "xSecLabel_" << std::endl;
    for(unsigned int ix=0; ix<xSecLabel_.size(); ++ix){
     std::cout << xSecLabel_[ix] << std::endl; 
    }
    exit(0);
  }


  // get folder prefix for systematics without extra rootfile
  switch (systematicVariationMod)
  {
    case sysPUUp                     : sysInputFolderExtension="PUup";   sysInputGenFolderExtension=sysInputFolderExtension; break;
    case sysPUDown                   : sysInputFolderExtension="PUdown"; sysInputGenFolderExtension=sysInputFolderExtension; break;
    case sysLepEffSFNormUp           : sysInputFolderExtension="EffSFNormUpSys";            break;
    case sysLepEffSFNormDown         : sysInputFolderExtension="EffSFNormDownSys";          break;
    case sysLepEffSFShapeUpEta       : sysInputFolderExtension="EffSFShapeUpEta";           break;
    case sysLepEffSFShapeDownEta     : sysInputFolderExtension="EffSFShapeDownEta";         break;
    case sysLepEffSFShapeUpPt        : sysInputFolderExtension="EffSFShapeUpPt";            break;
    case sysLepEffSFShapeDownPt      : sysInputFolderExtension="EffSFShapeDownPt";          break;
 // case sysLepEffSFShapeUpPt40      : sysInputFolderExtension="EffSFShapeUpPt40";          break;
 // case sysLepEffSFShapeDownPt40    : sysInputFolderExtension="EffSFShapeDownPt40";        break;
 // case sysLepEffSFNormUpStat       : sysInputFolderExtension="EffSFNormUpStat";           break;
 // case sysLepEffSFNormDownStat     : sysInputFolderExtension="EffSFNormDownStat";         break;
    case sysBtagSFUp                 : sysInputFolderExtension="BtagSFup";                  break;
    case sysBtagSFDown               : sysInputFolderExtension="BtagSFdown";                break;
 // case sysBtagSFShapeUpPt100       : sysInputFolderExtension="BTagSFShapeUpPt100";        break;
 // case sysBtagSFShapeDownPt100     : sysInputFolderExtension="BTagSFShapeDownPt100";      break;
 // case sysBtagSFShapeUpEta1p2      : sysInputFolderExtension="BTagSFShapeUpEta1p2";       break;
 // case sysBtagSFShapeDownEta1p2    : sysInputFolderExtension="BTagSFShapeDownEta1p2";     break;
    case sysBtagSFShapeUpPt65        : sysInputFolderExtension="BTagSFShapeUpPt65";         break;
    case sysBtagSFShapeDownPt65      : sysInputFolderExtension="BTagSFShapeDownPt65";       break;
    case sysBtagSFShapeUpEta0p7      : sysInputFolderExtension="BTagSFShapeUpEta0p7";       break;
    case sysBtagSFShapeDownEta0p7    : sysInputFolderExtension="BTagSFShapeDownEta0p7";     break;
 // case sysBTagSFFullShapeUpPt65    : sysInputFolderExtension="BTagSFFullShapeUpPt65";     break;
 // case sysBTagSFFullShapeDownPt65  : sysInputFolderExtension="BTagSFFullShapeDownPt65";   break;
 // case sysBTagSFFullShapeUpEta0p7  : sysInputFolderExtension="BTagSFFullShapeUpEta0p7";   break;
 // case sysBTagSFFullShapeDownEta0p7: sysInputFolderExtension="BTagSFFullShapeDownEta0p7"; break;
    case sysMisTagSFUp               : sysInputFolderExtension="MisTagSFup";                break;
    case sysMisTagSFDown             : sysInputFolderExtension="MisTagSFdown";              break;
    default: break;
  }
  TString sysInputFolderExtensionRaw=sysInputFolderExtension;
  // additional (control plot folders): NoWeight, OnlyPUWeight, NoBtagSFWeight;

  // add folder extensions for xSec from different selection step if renamed default folder is used
  if(addSel!="") sysInputFolderExtension+=addSel;

  // choose correct input folder for b-quark or b-jet
  TString recPartonBpath= "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension;
  //if(addSel!="") recPartonBpath+=addSel;
  TString recHadronBpath= "analyzeTopRecoKinematicsBjets" +sysInputFolderExtension;
  TString genPartonBpath= "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension;
  TString genHadronBpath= "analyzeTopHadronLevelKinematicsBjetsPhaseSpace"+sysInputGenFolderExtension;
  TString recBpath = ( (!extrapolate&&hadron) ? recHadronBpath : recPartonBpath );
  //TString genBpath = ( (!extrapolate&&hadron) ? genHadronBpath : genPartonBpath );
  TString genBpath = genHadronBpath;
  TString recBlabel = ( (!extrapolate&&hadron) ? "Rec" : "" );
  //TString genBlabel = ( (!extrapolate&&hadron) ? "Gen" : "" );
  TString genBlabel = "Gen";

  // choose correct input folder for status 1 or 3 lepton
  TString recPartonLeppath= "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension;
  //if(addSel!="") recPartonLeppath+=addSel;
  TString recHadronLeppath= "analyzeTopRecoKinematicsLepton" +sysInputFolderExtension;
  TString genPartonLeppath= "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension;
  TString genHadronLeppath= "analyzeTopHadronLevelKinematicsLeptonPhaseSpace"+sysInputGenFolderExtension;
  TString recLeppath  = ( (!extrapolate&&hadron) ? recHadronLeppath : recPartonLeppath );
  //TString genLeppath  = ( (!extrapolate&&hadron) ? genHadronLeppath : genPartonLeppath );
  TString genLeppath=genHadronLeppath;
  TString recLeplabel = ( (!extrapolate&&hadron) ? "Rec" : "" );
  //TString genLeplabel = ( (!extrapolate&&hadron) ? "Gen" : "" );
  TString genLeplabel =  "Gen";

  // choose correct input folder for mixed object analyzer
  TString recMixpath= "compositedKinematics";
  // warning: naming different than for std analyzer
  if(addSel==""||(addSel=="ProbSel"&&inputFolder.Contains("Prob")))  recMixpath+="KinFit";
  else if(addSel=="ProbSel") recMixpath+=addSel;
  recMixpath+=sysInputFolderExtensionRaw;
  TString genMixpath= "compositedHadronGenPhaseSpace";
  // FIXME1: no sys weight gen folder existing for mixed object analyzer
  genMixpath+=sysInputGenFolderExtension;
  // FIXME2: no sys weight reco folders existing for mixed object analyzer in Non-ttbar SG samples outside the .../Prob/ subfolder
  
  std::vector<TString> vecRedundantPartOfNameInNonTTbarSG_;
  //if(sysInputFolderExtensionRaw!=""&&!(addSel=="ProbSel"&&inputFolder.Contains("Prob"))){
  //  vecRedundantPartOfNameInNonTTbarSG_.push_back("Njets/"+sysInputFolderExtensionRaw);
  //  vecRedundantPartOfNameInNonTTbarSG_.push_back("rhos/"+sysInputFolderExtensionRaw);
  //}
  // debug
  //std::cout << genMixpath << std::endl;
  //std::cout << recMixpath << std::endl;
  //for(unsigned int i=0; i<vecRedundantPartOfNameInNonTTbarSG_.size(); ++i){
  //  std::cout << vecRedundantPartOfNameInNonTTbarSG_[i] << std::endl;
  //}

  // ====================
  //     choose plots
  // ====================
  // a) list plots you would like to see ("folder/plotName") - same as in .root files (for 1D and 2D)
  TString plots1D[ ] = { 
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt"    , // REC
    "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/topPt", // GEN
  };

  TString plots2D[ ] = { // reco - gen Match correlation plots (ttbar signal only)
    // b) response matrix top quantities
    "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/topPt_", // Response Matrix
  };

  // b) list plot axes style
  // 1D: "x-axis title"/"y-axis title"/log/rebin-factor
  // log = 0 or 1 for linear or logarithmic axis 

  TString axisLabel1D[ ] = { 
    xSecLabelName("topPt")+"/#frac{dN}{dp_{T}^{t}} #left[GeV^{-1}#right]/0/1",
    xSecLabelName("topPt")+" parton truth/events/0/1",
  };
  
  
  // 2D: "x-axis title"/"y-axis title"
  TString axisLabel2D[ ] = {// reco - gen Match correlation plots (ttbar signal only)
    // b) response matrix Top quantities
    xSecLabelName("topPt")+" gen/"+xSecLabelName("topPt")+" reco",
  };

  // count # plots
  unsigned int N1Dplots        = sizeof(plots1D)/sizeof(TString);
  unsigned int N2Dplots        = sizeof(plots2D)/sizeof(TString);

  // check if all axis labels exist
  if(N1Dplots        != sizeof(axisLabel1D)/sizeof(TString))       { std::cout << "ERROR: some 1D plots or axis label are missing"          << std::endl; exit(1); }

  if(N2Dplots        != sizeof(axisLabel2D)/sizeof(TString))       { std::cout << "ERROR: some 2D plots or axis label are missing"          << std::endl; exit(1); }
  
  // run automatically in batch mode if there are many canvas
  if((N1Dplots+N2Dplots)>15) gROOT->SetBatch();

  // ===================================
  //  Open standard analysis files
  // ===================================
  std::map<unsigned int, TFile*> files_, filesMu_, filesEl_;
  // open files depending on the decay channel
  if(decayChannel!="combined"){
    files_ = getStdTopAnalysisFiles(inputFolder, systematicVariationMod, dataFile, decayChannel, ttbarMC);
  }
  else{
    filesMu_ = getStdTopAnalysisFiles(inputFolder, systematicVariationMod, dataFileMu, "muon"    , ttbarMC);
    filesEl_ = getStdTopAnalysisFiles(inputFolder, systematicVariationMod, dataFileEl, "electron", ttbarMC);
  }
  // =====================
  //  Loading histos
  // =====================
  // collect all plot names in vector (first 1D, then 2D)
  std::vector<TString> plotList_;
  plotList_.insert(plotList_.begin(), plots1D, plots1D + sizeof(plots1D)/sizeof(TString));
  plotList_.insert(plotList_.end(),   plots2D, plots2D + sizeof(plots2D)/sizeof(TString));
  
  // remove irrelevant plots for systemtic variation to speed up

  // container for all histos (1D&2D)
  // example for acess: histo_["plotName"][sampleNr]
  std::map< TString, std::map <unsigned int, TH1F*> > histo_, histoEl_, histoMu_;
  std::map< TString, TH1F* > histoErrorBand_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2_, histo2El_, histo2Mu_;
  // total # plots 
  int Nplots=0;
  // save all histos from plotList_ that exist in files_ into 
  // histo_ and histo2_ and count total # of plots as Nplots
  if(verbose>1) std::cout << std::endl;
  std::vector<TString> vecRedundantPartOfNameInData;
  vecRedundantPartOfNameInData.push_back(sysInputFolderExtension);
  vecRedundantPartOfNameInData[0].ReplaceAll(addSel, "");
  if(decayChannel!="combined") getAllPlots(files_, plotList_, histo_, histo2_, N1Dplots, Nplots, verbose-1, decayChannel, &vecRedundantPartOfNameInData, false, addSelData, &vecRedundantPartOfNameInNonTTbarSG_);
  else{
    getAllPlots(filesEl_, plotList_, histoEl_, histo2El_, N1Dplots, Nplots, verbose-1, "electron", &vecRedundantPartOfNameInData, false, addSelData, &vecRedundantPartOfNameInNonTTbarSG_);
    getAllPlots(filesMu_, plotList_, histoMu_, histo2Mu_, N1Dplots, Nplots, verbose-1, "muon"    , &vecRedundantPartOfNameInData, false, addSelData, &vecRedundantPartOfNameInNonTTbarSG_);
  }
  //exit(0); // debug exit - use with high verbosity to see which plots are loaded
  // ===============================================
  // take care of rec level b quark plots from MC BG
  // ===============================================
  // this is necessary because signal and background do have different input folders
  if(!extrapolate&&hadron){
    if(decayChannel!="combined"){
      bool bqPtExistsInAnySample =false;
      bool bqEtaExistsInAnySample=false;
      bool bbbarPtExistsInAnySample =false;
      bool bbbarMassExistsInAnySample=false;
      bool lbMassExistsInAnySample=false;
      bool lepPtExistsInAnySample =false;
      bool lepEtaExistsInAnySample=false;
      // get input path/folder
      TString path=recPartonBpath;
      TString path2=recPartonBpath;
      TString pathL=recPartonLeppath;
      TString pathL2=recPartonLeppath;
      // loop samples
      for(unsigned int sample=kBkg; sample<=kSAToptW; ++sample){
	// get input path/folder
	path=recPartonBpath;
	path2=recPartonBpath;
	pathL=recPartonLeppath;
	pathL2=recPartonLeppath;
	if(sample==kData){
	  if (vecRedundantPartOfNameInData.size() != 0){
	    std::vector<TString>::iterator iter;
	    for (iter = (vecRedundantPartOfNameInData.begin()); iter != (vecRedundantPartOfNameInData.end()); iter++){
	      path.ReplaceAll( (*iter), ""); 
	      pathL.ReplaceAll((*iter), "");
	    }
	  }
	}
	else{
	   path.ReplaceAll (addSelData, ""); 
	   pathL.ReplaceAll(addSelData, "");
	}
	if(verbose>2){
	  std::cout << std::endl << sampleLabel(sample, decayChannel) << std::endl;
	  std::cout << path+"/bqPt"    << " -> " << path2+"/bqPt"    << std::endl;
	  std::cout << path+"/bqEta"   << " -> " << path2+"/bqEta"   << std::endl;
	  std::cout << path+"/bbbarPt"    << " -> " << path2+"/bbbarPt"    << std::endl;
	  std::cout << path+"/bbbarMass"  << " -> " << path2+"/bbbarMass"  << std::endl;
	  std::cout << path+"/lbMass"  << " -> " << path2+"/lbMass"  << std::endl;
	  std::cout << pathL+"/lepPt"  << " -> " << pathL2+"/lepPt"  << std::endl;
	  std::cout << pathL+"/lepEta" << " -> " << pathL2+"/lepEta" << std::endl;
	}
	// create plot container
	TH1* targetPlotBqPt  =0;
	TH1* targetPlotBqEta =0;
	TH1* targetPlotBbbarPt  =0;
	TH1* targetPlotBbbarMass=0;
	TH1* targetPlotLBMass=0;
	TH1* targetPlotLepPt =0;
	TH1* targetPlotLepEta=0;
	// check if file exists
	if(files_[sample]){
	  // get plot
	  files_[sample]->GetObject(path+"/bqPt" , targetPlotBqPt );
	  files_[sample]->GetObject(path+"/bqEta", targetPlotBqEta);
	  files_[sample]->GetObject(path+"/bbbarPt"  , targetPlotBbbarPt  );
	  files_[sample]->GetObject(path+"/bbbarMass", targetPlotBbbarMass);
	  files_[sample]->GetObject(path+"/lbMass", targetPlotLBMass);
	  files_[sample]->GetObject(pathL+"/lepPt" , targetPlotLepPt );
	  files_[sample]->GetObject(pathL+"/lepEta", targetPlotLepEta);
	  // Check existence of plot
	  if(targetPlotBqPt){ 
	    histo_[path2+"/bqPt" ][sample]=(TH1F*)targetPlotBqPt ->Clone(path+"/bqPt"); 
	    bqPtExistsInAnySample =true;
	  }
	  if(targetPlotBqEta){
	    histo_[path2+"/bqEta"][sample]=(TH1F*)targetPlotBqEta->Clone(path+"/bqEta");
	    bqEtaExistsInAnySample=true;
	  }
	  if(targetPlotBbbarPt){ 
	    histo_[path2+"/bbbarPt"  ][sample]=(TH1F*)targetPlotBbbarPt->Clone(path+"/bbbarPt");
	    bbbarPtExistsInAnySample=true;
	  }
	  if(targetPlotBbbarMass){
	    histo_[path2+"/bbbarMass"][sample]=(TH1F*)targetPlotBbbarMass->Clone(path+"/bbbarMass");
	    bbbarMassExistsInAnySample=true;
	  }
	  if(targetPlotLBMass){
	    histo_[path2+"/lbMass"][sample]=(TH1F*)targetPlotLBMass->Clone(path+"/lbMass");
	    lbMassExistsInAnySample=true;
	  }
	  if(targetPlotLepPt){ 
	    histo_[pathL2+"/lepPt" ][sample]=(TH1F*)targetPlotLepPt ->Clone(pathL+"/lepPt" ); 
	    lepPtExistsInAnySample =true;
	  }
	  if(targetPlotLepEta){
	    histo_[pathL2+"/lepEta"][sample]=(TH1F*)targetPlotLepEta->Clone(pathL+"/lepEta");
	    lepEtaExistsInAnySample=true;
	  }
	  if(verbose>1){
	    std::cout << sampleLabel(sample,decayChannel) << ":" << std::endl;
	    std::cout << "targetPlotBqPt  : " <<  targetPlotBqPt   << std::endl;
	    std::cout << "targetPlotBqEta : " <<  targetPlotBqEta	<< std::endl;
	    std::cout << "targetPlotBbbarPt  : " << targetPlotBbbarPt   << std::endl;
	    std::cout << "targetPlotBbbarMass: " << targetPlotBbbarMass << std::endl;
	    std::cout << "targetPlotLBMass: " <<  targetPlotLBMass  << std::endl;
	    std::cout << "targetPlotLepPt : " <<  targetPlotLepPt  << std::endl;
	    std::cout << "targetPlotLepEta: " <<  targetPlotLepEta << std::endl;
	  }
	}
	// delete plot container
	delete targetPlotBqPt;
	delete targetPlotBqEta;
	delete targetPlotBbbarPt;
	delete targetPlotBbbarMass;
	delete targetPlotLBMass;
	delete targetPlotLepPt;
	delete targetPlotLepEta;
      }
      // check if plot exists at all
      if(!bqPtExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bqPt for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!bqEtaExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bqEta for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!bbbarPtExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bbbarPt for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!bbbarMassExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bbbarMass for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!lbMassExistsInAnySample){
	std::cout << "no plot found with label " << path+"/lbMass for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!lepPtExistsInAnySample){
	std::cout << "no plot found with label " << pathL+"/lepPt for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!lepEtaExistsInAnySample){
	std::cout << "no plot found with label " << pathL+"/lepEta for non ttbar signal" << std::endl;
	exit(0);
      }
    }
    else{
      bool bqPtExistsInAnySample =false;
      bool bqEtaExistsInAnySample=false;
      bool bbbarPtExistsInAnySample =false;
      bool bbbarMassExistsInAnySample=false;
      bool lbMassExistsInAnySample=false;
      bool lepPtExistsInAnySample =false;
      bool lepEtaExistsInAnySample=false;
      // get input path/folder
      TString path=recPartonBpath;
      TString path2=recPartonBpath;
      TString pathL=recPartonLeppath;
      TString pathL2=recPartonLeppath;
      // loop samples
      for(unsigned int sample=kBkg; sample<=kSAToptW; ++sample){
      path=recPartonBpath;
      path2=recPartonBpath;
      pathL=recPartonLeppath;
      pathL2=recPartonLeppath;
	if(sample==kData){
	  if (vecRedundantPartOfNameInData.size() != 0){
	    std::vector<TString>::iterator iter;
	    for (iter = (vecRedundantPartOfNameInData.begin()); iter != (vecRedundantPartOfNameInData.end()); iter++){
	      path.ReplaceAll( (*iter), ""); 
	      pathL.ReplaceAll((*iter), "");
	    }
	  }
	}
	else{
	  path.ReplaceAll (addSelData, ""); 
	  pathL.ReplaceAll(addSelData, "");
	}
	if(verbose>2){
	  std::cout << std::endl << sampleLabel(sample,decayChannel) << std::endl;
	  std::cout << path+"/bqPt"    << " -> " << path2+"/bqPt"    << std::endl;
	  std::cout << path+"/bqEta"   << " -> " << path2+"/bqEta"   << std::endl;
	  std::cout << path+"/bbbarPt"    << " -> " << path2+"/bbbarPt"    << std::endl;
	  std::cout << path+"/bbbarMass"  << " -> " << path2+"/bbbarMass"  << std::endl;
	  std::cout << path+"/lbMass"  << " -> " << path2+"/lbMass"  << std::endl;
	  std::cout << pathL+"/lepPt"  << " -> " << pathL2+"/lepPt"  << std::endl;
	  std::cout << pathL+"/lepEta" << " -> " << pathL2+"/lepEta" << std::endl;
	}
	// create plot container
	TH1* targetPlotBqPt  =0;
	TH1* targetPlotBqEta =0;
	TH1* targetPlotBbbarPt  =0;
	TH1* targetPlotBbbarMass=0;
	TH1* targetPlotLBMass=0;
	TH1* targetPlotLepPt =0;
	TH1* targetPlotLepEta=0;
	TH1* targetMuPlotBqPt  =0;
	TH1* targetMuPlotBqEta =0;
	TH1* targetMuPlotBbbarPt  =0;
	TH1* targetMuPlotBbbarMass=0;
	TH1* targetMuPlotLBMass=0;
	TH1* targetMuPlotLepPt =0;
	TH1* targetMuPlotLepEta=0;
	// check if file exists
	if(filesMu_[sample]&&filesEl_[sample]){
	  if(verbose>2) std::cout << "-> file for exists" << std::endl;
	  // get plot
	  filesEl_[sample]->GetObject(path+"/bqPt"   , targetPlotBqPt  );
	  filesEl_[sample]->GetObject(path+"/bqEta"  , targetPlotBqEta );
	  filesEl_[sample]->GetObject(path+"/bbbarPt"  , targetPlotBbbarPt  );
	  filesEl_[sample]->GetObject(path+"/bbbarMass", targetPlotBbbarMass);
	  filesEl_[sample]->GetObject(path+"/lbMass", targetPlotLBMass);
	  filesEl_[sample]->GetObject(pathL+"/lepPt" , targetPlotLepPt );
	  filesEl_[sample]->GetObject(pathL+"/lepEta", targetPlotLepEta);
	  filesMu_[sample]->GetObject(path+"/bqPt"   , targetMuPlotBqPt  );
	  filesMu_[sample]->GetObject(path+"/bqEta"  , targetMuPlotBqEta );
	  filesMu_[sample]->GetObject(path+"/bbbarPt"  , targetMuPlotBbbarPt  );
	  filesMu_[sample]->GetObject(path+"/bbbarMass", targetMuPlotBbbarMass); 
	  filesMu_[sample]->GetObject(path+"/lbMass" , targetMuPlotLBMass);
	  filesMu_[sample]->GetObject(pathL+"/lepPt" , targetMuPlotLepPt );
	  filesMu_[sample]->GetObject(pathL+"/lepEta", targetMuPlotLepEta);
	  // Check existence of plot
	  if(targetPlotBqPt&&targetMuPlotBqPt){ 
	    //std::cout << "canv there for" << sampleLabel(sample, decayChannel) << std::endl;
	    histoEl_[path2+"/bqPt" ][sample]=(TH1F*)targetPlotBqPt  ->Clone(path+"/bqPt"); 
	    histoMu_[path2+"/bqPt" ][sample]=(TH1F*)targetMuPlotBqPt->Clone(path+"/bqPt"); 
	    bqPtExistsInAnySample =true;
	  }
	  if(targetPlotBqEta&&targetMuPlotBqEta){
	    histoEl_[path2+"/bqEta"][sample]=(TH1F*)targetPlotBqEta  ->Clone(path+"/bqEta");
	    histoMu_[path2+"/bqEta"][sample]=(TH1F*)targetMuPlotBqEta->Clone(path+"/bqEta");
	    bqEtaExistsInAnySample=true;
	  }

	  if(targetPlotBbbarPt){ 
	    histoEl_[path2+"/bbbarPt"  ][sample]=(TH1F*)targetPlotBbbarPt  ->Clone(path+"/bbbarPt");
	    histoMu_[path2+"/bbbarPt"  ][sample]=(TH1F*)targetMuPlotBbbarPt->Clone(path+"/bbbarPt");
	    bbbarPtExistsInAnySample=true;
	  }
	  if(targetPlotBbbarMass){
	    histoEl_[path2+"/bbbarMass"][sample]=(TH1F*)targetPlotBbbarMass  ->Clone(path+"/bbbarMass");
	    histoMu_[path2+"/bbbarMass"][sample]=(TH1F*)targetMuPlotBbbarMass->Clone(path+"/bbbarMass");
	    bbbarMassExistsInAnySample=true;
	  }
	  if(targetPlotLBMass){
	    histoEl_[path2+"/lbMass"][sample]=(TH1F*)targetPlotLBMass  ->Clone(path+"/lbMass");
	    histoMu_[path2+"/lbMass"][sample]=(TH1F*)targetMuPlotLBMass->Clone(path+"/lbxMass");
	    lbMassExistsInAnySample=true;
	  }
	  if(targetPlotLepPt&&targetMuPlotLepPt){ 
	    histoEl_[pathL2+"/lepPt" ][sample]=(TH1F*)targetPlotLepPt  ->Clone(pathL+"/lepPt" );
	    histoMu_[pathL2+"/lepPt" ][sample]=(TH1F*)targetMuPlotLepPt->Clone(pathL+"/lepPt" );  
	    lepPtExistsInAnySample =true;
	  }
	  if(targetPlotLepEta&&targetMuPlotLepEta){
	    histoEl_[pathL2+"/lepEta"][sample]=(TH1F*)targetPlotLepEta  ->Clone(pathL+"/lepEta");
	    histoMu_[pathL2+"/lepEta"][sample]=(TH1F*)targetMuPlotLepEta->Clone(pathL+"/lepEta");
	    lepEtaExistsInAnySample=true;
	  }
	  if(verbose>1){
	    std::cout << sampleLabel(sample,decayChannel) << ":" << std::endl;
	    std::cout << "targetPlotBqPt  : " <<  targetPlotBqPt   << std::endl;
	    std::cout << "targetPlotBqEta : " <<  targetPlotBqEta	<< std::endl;
	    std::cout << "targetPlotBbbarPt  : " << targetPlotBbbarPt	<< std::endl;
	    std::cout << "targetPlotBbbarMass: " << targetPlotBbbarMass << std::endl;
	    std::cout << "targetPlotLBMass: " <<  targetPlotLBMass << std::endl;
	    std::cout << "targetPlotLepPt : " <<  targetPlotLepPt  << std::endl;
	    std::cout << "targetPlotLepEta: " <<  targetPlotLepEta << std::endl;
	  }
	}
	// delete plot container
	delete targetPlotBqPt;
	delete targetPlotBqEta;
	delete targetPlotBbbarPt;
	delete targetPlotBbbarMass;
	delete targetPlotLBMass;
	delete targetPlotLepPt;
	delete targetPlotLepEta;
	delete targetMuPlotBqPt;
	delete targetMuPlotBqEta;
	delete targetMuPlotBbbarPt;
	delete targetMuPlotBbbarMass;
	delete targetMuPlotLBMass;
	delete targetMuPlotLepPt;
	delete targetMuPlotLepEta;
      }
      // check if plot exists at all
      if(!bqPtExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bqPt for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!bqEtaExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bqEta for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!bbbarPtExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bbbarPt for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!bbbarMassExistsInAnySample){
	std::cout << "no plot found with label " << path+"/bbbarMass for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!lbMassExistsInAnySample){
	std::cout << "no plot found with label " << path+"/lbMass for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!lepPtExistsInAnySample){
	std::cout << "no plot found with label " << pathL+"/lepPt for non ttbar signal" << std::endl;
	exit(0);
      }
      if(!lepEtaExistsInAnySample){
	std::cout << "no plot found with label " << pathL+"/lepEta for non ttbar signal" << std::endl;
	exit(0);
      }
    }
  }

  // ====================================================================================
  //  Renaming of PS specific plot names to parton level PS names to unify naming
  // ====================================================================================
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    TString newName=plotList_[plot];
    // ---
    //    a) for parton level FULL PS configuration
    // ---
    // change gen level full PS label to gen level PS label
    if(plotList_[plot].Contains("analyzePartonLevel")&&!plotList_[plot].Contains("PhaseSpace")){
      newName.ReplaceAll("analyzePartonLevelKinematics", "analyzePartonLevelKinematicsPhaseSpace");
    }
    // ---
    //    b) for hadron level PS configuration
    // ---
    //  b1) change b-jet to b-quark label
    // for rec plots
    if(plotList_[plot].Contains(recHadronBpath)){ 
      newName.ReplaceAll(recHadronBpath, recPartonBpath);
      newName.ReplaceAll("Rec", "");
      newName.ReplaceAll("Topo", "TopReco");
    }
    // for gen plots
    if(plotList_[plot].Contains(genHadronBpath)){ 
      newName.ReplaceAll(genHadronBpath, genPartonBpath);
      if(!plotList_[plot].Contains("_")) newName.ReplaceAll("Gen", "");
    }
    //  b2) change lepton status 1 to lepton status 3 label
    // for rec plots
    if(plotList_[plot].Contains(recHadronLeppath)){ 
      newName.ReplaceAll(recHadronLeppath, recPartonLeppath);
      newName.ReplaceAll("Rec", "");
      newName.ReplaceAll("Topo", "TopReco");
    }
    // for mixed object analyzer
    if(plotList_[plot].Contains("composited")){
      newName.ReplaceAll("composited", "analyzeTop");
      newName.ReplaceAll("Ngenjets", "Njets");
      newName.ReplaceAll("rhosGen" , "rhos" );
      if(sysInputGenFolderExtension!=""&&newName.Contains("Gen"+PS)&&!newName.Contains(sysInputGenFolderExtension)) newName.ReplaceAll("Gen"+PS, "Gen"+PS+sysInputGenFolderExtension);
      if(plotList_[plot].Contains("Gen")) newName.ReplaceAll("Gen", "LevelKinematics");
      else{
	newName.ReplaceAll("analyzeTop", "analyzeTopReco"  );
	if(!newName.Contains("KinFit")) newName.ReplaceAll("Kinematics", "KinematicsKinFit");
      }
      TString syAdd=sysInputFolderExtension;
      syAdd.ReplaceAll(addSel, "");
      if(syAdd!=""&&!newName.Contains(syAdd)) newName.ReplaceAll("KinFit","KinFit"+syAdd);
      if(addSel.Contains("ProbSel")&&!newName.Contains("ProbSel")) newName.ReplaceAll("/", "ProbSel/");
	if(addSel.Contains("ProbSel")&&newName.Contains("ProbSel")&&(newName.Contains("Hadron")||newName.Contains("Parton"))) newName.ReplaceAll("ProbSel", "");
    }
    // for all gen plots
    if(plotList_[plot].Contains(genHadronLeppath)){ 
      newName.ReplaceAll(genHadronLeppath, genPartonLeppath);
      if(!plotList_[plot].Contains("_")) newName.ReplaceAll("Gen", "");
    }
    // b3) change hadron PS gen level label to parton PS gen level label
    if(newName.Contains(LV)&&newName.Contains("PhaseSpace")) newName.ReplaceAll(LV, "Parton");
    // b4) remove phase space form name
    if(newName.Contains("PhaseSpace")) newName.ReplaceAll("PhaseSpace", "");
    // check if replacement is necessary
    if(newName!=plotList_[plot]){
      if(verbose>=1||(plotList_[plot].Contains("composited")&&systematicVariation==testMe)) std::cout << plotList_[plot] << " -> " << newName << " for " << std::endl;
      // loop samples
      for(unsigned int sample=kSig; sample<=kData; ++sample){
	// check if plot exists 1D
	if((decayChannel!="combined"&&plotExists(histo_, plotList_[plot], sample))||(decayChannel=="combined"&&plotExists(histoEl_, plotList_[plot], sample)&&plotExists(histoMu_, plotList_[plot], sample))){
	  if(verbose>1) std::cout << sampleLabel(sample,decayChannel) << "(1D)" << std::endl;
	  // replace plot entry
	  if(decayChannel!="combined"){
	    histo_[newName][sample]=(TH1F*)histo_[plotList_[plot]][sample]->Clone(); 
	    histo_[plotList_[plot]].erase(sample);
	  }
	  else{
	    histoEl_[newName][sample]=(TH1F*)histoEl_[plotList_[plot]][sample]->Clone(); 
	    histoEl_[plotList_[plot]].erase(sample);
	    histoMu_[newName][sample]=(TH1F*)histoMu_[plotList_[plot]][sample]->Clone(); 
	    histoMu_[plotList_[plot]].erase(sample);
	  }
	}
	// check if plot exists 2D
	else if((decayChannel!="combined"&&plotExists(histo2_, plotList_[plot], sample))||(decayChannel=="combined"&&plotExists(histo2Mu_, plotList_[plot], sample)&&plotExists(histo2El_, plotList_[plot], sample))){
	  if(verbose>1) std::cout << sampleLabel(sample,decayChannel) << "(2D)" << std::endl;
	  // replace plot entry
	  if(decayChannel!="combined"){
	    histo2_[newName][sample]=(TH2F*)histo2_[plotList_[plot]][sample]->Clone(); 
	    histo2_[plotList_[plot]].erase(sample);
	  }
	  else{
	    histo2El_[newName][sample]=(TH2F*)histo2El_[plotList_[plot]][sample]->Clone(); 
	    histo2El_[plotList_[plot]].erase(sample); 
	    histo2Mu_[newName][sample]=(TH2F*)histo2Mu_[plotList_[plot]][sample]->Clone(); 
	    histo2Mu_[plotList_[plot]].erase(sample);
	  }
	}
	// finally delete the whole name entry
	// attention: some full PS plots still needed!!!
	if(sample==kSAToptW&&!extrapolate) histo_.erase(plotList_[plot]);
      }
      // replace plot list entry
      plotList_[plot]=newName;
    }
  }

  // ==========================================
  //  Lumiweighting for choosen luminosity
  // ==========================================
  // scale every histo in histo_ and histo2_ to the corresponding luminosity
  // Additionally the mu eff SF is applied
  // NOTE: luminosity [/pb]
  if(decayChannel!="combined") scaleByLuminosity(plotList_, histo_, histo2_, N1Dplots, luminosity, verbose-1, systematicVariationMod, decayChannel, ttbarMC2);
  else{
    scaleByLuminosity(plotList_, histoMu_, histo2Mu_, N1Dplots, luminosityMu, verbose-1, systematicVariationMod, "muon"    , ttbarMC2);
    scaleByLuminosity(plotList_, histoEl_, histo2El_, N1Dplots, luminosityEl, verbose-1, systematicVariationMod, "electron", ttbarMC2);
  }
  
  // ========================================================
  //  Add single top channels and DiBoson contributions
  // ========================================================
  // for SingleTop and DiBoson samples:
  // every plot plotList_ existing in sTop/ diBoson .root file
  // will be combined and saved in the histo_ and histo2_ map
  // reCreate: recreate combined plots if they are already existing
  bool reCreate=false;
  if(decayChannel!="combined") AddSingleTopAndDiBoson(plotList_, histo_, histo2_, N1Dplots, verbose-1, reCreate, decayChannel);
  else{
    AddSingleTopAndDiBoson(plotList_, histoMu_, histo2Mu_, N1Dplots, verbose-1, reCreate, "muon"    );
    AddSingleTopAndDiBoson(plotList_, histoEl_, histo2El_, N1Dplots, verbose-1, reCreate, "electron");
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
	if(verbose>1) std::cout << plotList_[plot] << " : " ;
	// a) 1D
	if((plot<N1Dplots)&&(histoMu_.count(plotList_[plot])>0)&&(histoMu_[plotList_[plot]].count(sample)>0)&&(histoEl_.count(plotList_[plot])>0)&&(histoEl_[plotList_[plot]].count(sample)>0)){ 
	  if(verbose>1) std::cout << "1D" << std::endl;
	  histo_[plotList_[plot]][sample]=     (TH1F*)(histoMu_[plotList_[plot]][sample]->Clone());
	  histo_[plotList_[plot]][sample]->Add((TH1F*)(histoEl_[plotList_[plot]][sample]->Clone()));
	}
	// b) 2D
	else if((plot>=N1Dplots)&&(histo2Mu_.count(plotList_[plot])>0)&&(histo2Mu_[plotList_[plot]].count(sample)>0)&&(histo2El_.count(plotList_[plot])>0)&&(histo2El_[plotList_[plot]].count(sample)>0)){
	  if(verbose>1) std::cout << "2D" << std::endl;
	  histo2_[plotList_[plot]][sample]=     (TH2F*)(histo2Mu_[plotList_[plot]][sample]->Clone());
	  histo2_[plotList_[plot]][sample]->Add((TH2F*)(histo2El_[plotList_[plot]][sample]->Clone()));
	}
	else{
	  if(verbose>1) std::cout << "NOT FOUND" << std::endl;
	}
      }   
    }
  }

  // ============================
  //  Configure histograms
  // ============================
  // needs: plotList_, histo_, histo2_, N1Dplots, axisLabel_, axisLabel1D, axisLabel2D
  std::vector<TString> axisLabel_;

  axisLabel_.insert(axisLabel_.begin(), axisLabel1D, axisLabel1D + sizeof(axisLabel1D)/sizeof(TString));
  axisLabel_.insert(axisLabel_.end(),   axisLabel2D, axisLabel2D + sizeof(axisLabel2D)/sizeof(TString));

  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    if(decayChannel=="electron"&&axisLabel_[plot].Contains("#mu")) axisLabel_[plot].ReplaceAll("#mu", "e"  );
    if(decayChannel=="combined"&&axisLabel_[plot].Contains("#mu")) axisLabel_[plot].ReplaceAll("#mu", "lep");
    if(verbose>1){
      std::cout << "(plot, x Axis label , y Axis label , log scale?, rebinning factor):" << std::endl;
      std::cout << plotList_[plot] << ": " << getStringEntry(axisLabel_[plot],1);
      std::cout << " , " << getStringEntry(axisLabel_[plot], 2);
      if(plot<N1Dplots){
	std::cout<< " , " << getStringEntry(axisLabel_[plot], 3);
	std::cout<< " , " << getStringEntry(axisLabel_[plot], 4);
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
	// default
	histogramStyle( *histo_[plotList_[plot]][sample], sample, true);
	// special configurations
	if(getStringEntry(plotList_[plot], 2)=="PartonJetDRall")histo_[plotList_[plot]][sample]->SetNdivisions(816);
	// set QCD to 0
	if(setQCDtoZero&&sample==kQCD) histo_[plotList_[plot]][sample]->Scale(0.);
      }
      // b) 2D
      if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)) histStyle2D( *histo2_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2));
    }
  }
  // =================================================
  //  Clone relevant plots for SVD unfolding xSec
  // =================================================
  int kAllMC=kSAToptW+1;
  for(unsigned int var=0; var<xSecVariables_.size(); ++var){
    TString variable=xSecVariables_[var];
    TString varname="analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable;
    if(!plotExists(histo_, varname, kBkg)) std::cout << " ERROR - Variable " << varname << " can not be found for MC BG, check name, also for separate loading for lep/bq quantities in hadron level PS" << std::endl;
    // ttbar BG yield for signal fraction
    histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kBkg]=(TH1F*)histo_[varname][kBkg]->Clone(variable);
    // create combined BG reco plot
    histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kAllMC]=(TH1F*)histo_[varname][kBkg]->Clone(variable);
    for(int bgsample=kZjets; bgsample<=kDiBos; ++bgsample){
      if(!histo_[varname][bgsample]){
	  std::cout << "missing plot: analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable << " for sample " << sampleLabel(bgsample,decayChannel) << std::endl;
	exit(1);
      }
      if(verbose>1) std::cout << sampleLabel(bgsample,decayChannel) << std::endl;
      histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kAllMC]->Add((TH1F*)histo_[varname][bgsample]->Clone(variable));
    }
    // data event yield
    histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData]=(TH1F*)histo_[varname][kData]->Clone(variable);
    // signal reco plot
    histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kSig]=(TH1F*)histo_[varname][kSig]->Clone(variable);
    // signal gen plot
    if(!histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable][kSig]){
      std::cout << "missing plot: " << "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable << " for ttbar signal" << std::endl;
    }
    else histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable][kSig]=(TH1F*)histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable][kSig]->Clone(variable);
    // response matrix plot
    if(!histo2_[varname+"_"][kSig]){
      std::cout << "missing plot: " << varname+"_" << " for ttbar signal" << std::endl;
    }
    else histo2_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable+"_"][kSig]=(TH2F*)histo2_[varname+"_"][kSig]->Clone(variable);
  }
  
  // ============================
  //  Rebinning 1D histogramme
  // ============================
  // create variable bin edges
  std::map<TString, std::vector<double> > binning_ = makeVariableBinning(false);
  // loop samples
  for(unsigned int sample=kSig; sample<=kData; ++sample){
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      TString plotName=plotList_[plot];
      // check if plot exists and is 1D
      if((plot<N1Dplots)&&(plotExists(histo_, plotName, sample))){
	// equidistant binning
	double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4)).c_str());
	if(reBinFactor>1&&binning_.count(plotName)==0){
	  equalReBinTH1(reBinFactor, histo_, plotName, sample);
	}
	// variable binning
	else{
	  // check existence of e.g. binning_["topPt"]
	  if(binning_.count(getStringEntry(plotName, 2))>0){
	    reBinTH1F(*histo_[plotName][sample], binning_[getStringEntry(plotName, 2)], verbose-1);
	    // remove additional bins with 0 width
	    double *newBinLowEdges = new double[binning_[getStringEntry(plotName, 2)].size()];
	    newBinLowEdges = &binning_[getStringEntry(plotName, 2)].front();
	    histo_[plotName][sample]= (TH1F*)histo_[plotName][sample]->Rebin(binning_[getStringEntry(plotName, 2)].size()-1, histo_[plotName][sample]->GetName(), newBinLowEdges); 
	    if(verbose>1){
	      for(int i=1; i<= histo_[plotName][sample]->GetNbinsX()+1; i++){	       
		std::cout << "bin " << i << ": " << histo_[plotName][sample]->GetBinContent(i) << " / ";
		std::cout << histo_[plotName][sample]->GetBinWidth(i) << std::endl;
	      }
	    }
	    // divide by binwidth
	    histo_[plotName][sample] = divideByBinwidth(histo_[plotName][sample], verbose-1);
	    if(verbose>1){
	      for(int i=1; i<= histo_[plotName][sample]->GetNbinsX()+1; i++){
		std::cout << "bin " << i << "= " << histo_[plotName][sample]->GetBinContent(i) << std::endl;
	      }
	    }
	    // change y axis label
	    if(!axisLabel_[plot].Contains("binwidth")){
	      //axisLabel_[plot].ReplaceAll(getStringEntry(axisLabel_[plot],2), getStringEntry(axisLabel_[plot],2)+" per binwidth");
	      // check if unit is GeV
	      if(getStringEntry(axisLabel_[plot],1).Contains("GeV")){
		//axisLabel_[plot].ReplaceAll(getStringEntry(axisLabel_[plot],2), getStringEntry(axisLabel_[plot],2)+" [GeV^{-1}]");
	      }
	    }
	  }
	}
      }
    }
  }
  
  // ===================================
  //  MC efficiency determination
  // ===================================
  // count # of efficiency histos
  unsigned int NMCeff=0;
  // loop all variables
  for(unsigned int number=0; number<xSecVariables_.size(); ++number){
    TString variable=xSecVariables_[number];    
    TString efficiency="efficiency/"+variable;
    // check if gen and reco plots are available
    if(extrapolate&&plotExists(histo_, "analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable, kSig)&&plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable, kSig)){
      // std::cout << "found gen and reco" << std::endl;
      // get reco plot
      histo_[efficiency][kSig]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      // divide by gen plot
      histo_[efficiency][kSig]->Divide((TH1F*)(histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable][kSig]->Clone()));
      // add plot to list of plots
      plotList_.push_back(efficiency);
      // add axis configuration
      unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
      TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
      axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#epsilon #times A (MC t#bar{t}#rightarrow#mu)/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e"  );
      if(decayChannel=="combined") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "lep");
      // set binomial errors
      if(verbose>1) std::cout << "       eff, events, width, sqrt(eff*(1.-eff)/events" << std::endl;
      for(int bin=1; bin<=histo_[efficiency][kSig]->GetNbinsX(); ++bin){
	double eff=histo_[efficiency][kSig]->GetBinContent(bin); 
	double N=histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable][kSig]->GetBinContent(bin);
	double width=histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable][kSig]->GetBinWidth(bin);
	N*=width*effSFAB(systematicVariationMod,decayChannel);
	if(decayChannel!="combined"){
	  if(ttbarMC2=="Powheg") N/=(lumiweight(kSigPow, luminosity, systematicVariationMod, decayChannel));
	  else N/=(lumiweight(kSig, luminosity, systematicVariationMod, decayChannel));
	}
	if(verbose>1){
	  std::cout << "bin " << bin << ": " << eff << ", " << N << ", " << width;
	  std::cout << ", " <<  sqrt(eff*(1.-eff)/N) << std::endl;
	}
	histo_[efficiency][kSig]->SetBinError(bin, sqrt(eff*(1.-eff)/N));
	// set error to 0 for bins with width 0 or 0 entries
	if(N==0) histo_[efficiency][kSig]->SetBinError(bin, 0);
      }
      // change histo style
      histogramStyle(*histo_[efficiency][kSig], kSig, false, 0.6);
      ++NMCeff;
    }
  }

  // ===============================================================
  //  Differential cross section (extrapolated or PS) determination
  // ===============================================================
  // count # of cross section histos
  unsigned int NXSec=0;
  // loop all variables
  for(unsigned int number=0; number<xSecVariables_.size(); ++number){
    TString variable=xSecVariables_[number];
    TString label =getStringEntry(xSecLabel_[number],1);
    TString label2=getStringEntry(xSecLabel_[number],2);
    if(label2=="[GeV]")  label2="#left[pb GeV^{-1}#right]";
    else if(label2==" ") label2="[pb]";
    TString xSec ="xSec/";
    if(compare) xSec+="BBB";
    xSec+=variable;
    bool calculateXSec=true;
    // loop samples 
    for(unsigned int sample=kSig; sample<=kData; ++sample){
      // check if reco MC plot and data event yield plots are available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable, sample)){
	calculateXSec=false;
      }
    }
    // check if efficiency plot is available
    if(!plotExists(histo_, "efficiency/"+variable, kSig)){
      calculateXSec=false;
    }
    if(calculateXSec&&(!SVDunfold||compare)){
      // a) differential XSec from data
      // get data plot
      TString name=TString(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetName())+"kData";
      histo_[xSec][kData]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->Clone(name));
      // subtract BG(MC)
      //histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg  ]->Clone()), -1); // other ttbar BG is accessed by a signal fraction (see below) 
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos]->Clone()), -1);
      histo_[xSec][kData]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ]->Clone()), -1);
      // correct for fraction of ttbar events that are signal
      histo_[xSec+"AllTtBar"][kSig]=     (TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      histo_[xSec+"AllTtBar"][kSig]->Add((TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg]->Clone()));
      histo_[xSec+"SignalFraction"][kSig]=(TH1F*)histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone();
      histo_[xSec+"SignalFraction"][kSig]->Divide(histo_[xSec+"AllTtBar"][kSig]);
      histo_[xSec][kData]->Multiply(histo_[xSec+"SignalFraction"][kSig]);
      // apply efficiency correction
      DivideYieldByEfficiencyAndLumi(histo_[xSec][kData], ((TH1F*)(histo_["efficiency/"+variable][kSig])), luminosity, 0);
      // b) differential XSec from Signal(MC prediction)
      histo_[xSec][kSig]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      DivideYieldByEfficiencyAndLumi(histo_[xSec][kSig], (TH1F*)(histo_["efficiency/"+variable][kSig]), luminosity, 0);
      // add plot to list of plots
      plotList_.push_back(xSec);
      // add axis configuration
      unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
      TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
      recoAxisLabel.ReplaceAll("KinFit ","");
      //      axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{d#sigma}{d"+label+"} "+label2+" (t#bar{t}#rightarrow #mu prompt)/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{d#sigma}{d"+label+"} "+label2+"/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e"  );
      if(decayChannel=="combined") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "lep");
      // configure xSec plot histo style
      histogramStyle(*histo_[xSec][kData], kData, false);
      histogramStyle(*histo_[xSec][kSig ], kSig , false);
      if(compare){
	histo_[xSec][kData]->SetLineColor  (38);
	histo_[xSec][kData]->SetMarkerColor(38);
	histo_[xSec][kData]->SetMarkerStyle(33);
	histo_[xSec][kSig ]->SetLineColor  (kOrange);
	histo_[xSec][kSig ]->SetMarkerColor(kOrange);
      }
      // restrict axis
      setXAxisRange(histo_[xSec][kData], variable);
      setXAxisRange(histo_[xSec][kSig ], variable);
      if(!compare)++NXSec;
    }
  }
  
  // ---
  //    differential normalized cross section (phase space) determination
  // ---
  if(!SVDunfold||compare){
    // loop all variables
    for(unsigned int number=0; number<xSecVariables_.size(); ++number){
      TString variable=xSecVariables_[number];
      TString label =getStringEntry(xSecLabel_[number],1);
      TString label2=getStringEntry(xSecLabel_[number],2);
      if(label2=="[GeV]")  label2="#left[GeV^{-1}#right]";
      else if(label2==" ") label2="";
      TString xSec ="xSecNorm/";
      if(compare) xSec+="BBB";
      xSec+=variable;
      bool calculateXSec=true;
      // loop samples 
      for(unsigned int sample=kSig; sample<=kData; ++sample){
	// check if reco MC plot and data event yield plots are available
	if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable, sample)){
	  calculateXSec=false;
	}
      }
      // check if efficiency plot is available
      if(!plotExists(histo_, "efficiency/"+variable, kSig)){
	calculateXSec=false;
      }
      if(calculateXSec){
	// a) differential norm XSec from data
	// get data plot
	TString name=TString(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetName())+"kData";
	if(compare) name+="BBB";
	histo_[xSec][kData]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->Clone(name));
	// total number of data events for this variable
	double NdataVariable=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData],verbose-1);
	// total number of BG events for this variable
	double NBGVariable=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets],verbose-1);
	// NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg],verbose-1);
	// BG from other ttbar events is treated in the signal fraction
	NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ],verbose-1);
	NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos],verbose-1);
	NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ],verbose-1);
	NBGVariable+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets],verbose-1);
	// ttbar signal fraction
	double NttbarAll=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig],verbose-1);
	NttbarAll+=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg],verbose-1);
	double sigFracVariable=getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig],verbose-1);
	sigFracVariable/=NttbarAll;
	// inclusive efficiency (* A) for this variable
	double effIncl= getInclusiveXSec(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kSig],verbose-1);
	effIncl/= getInclusiveXSec(histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/"+variable][kSig],verbose-1);
	// loop bins
	for(int bin =0; bin<=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetNbinsX()+1; ++bin){
	  // for inclusive values of xSecResult, effA, BR, Ndata and NBG 
	  // which are used in the calculation are taken from section
	  // "total inclusive cross section (ttbar->X) determination" above
	  // remark: event numbers are taken from yield plots-> *binwidth
	  // get data entries in bin

	  // get binwidth
	  double binwidth=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinWidth(bin);
	  double NdataBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinContent(bin);
	  NdataBin*=binwidth;
	  // get BG in bin
	  double NBGBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets]->GetBinContent(bin);
	  //NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg]->GetBinContent(bin);
	  // BG from other ttbar events is treated in the signal fraction
	  NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets]->GetBinContent(bin);
	  NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ]->GetBinContent(bin);
	  NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos]->GetBinContent(bin);
	  NBGBin+=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ]->GetBinContent(bin);
	  NBGBin*=binwidth;
	  // ttbar signal fraction
	  double NttbarAllBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->GetBinContent(bin);
	  NttbarAllBin      +=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg]->GetBinContent(bin);
	  double sigFracBin=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->GetBinContent(bin);
	  sigFracBin/=NttbarAllBin;
	  if(!isfinite(sigFracBin)||isnan(sigFracBin)) sigFracBin=1;
	  // calculate Ni = number of ttbar only events in bin from data
	  double Ni= NdataBin-NBGBin;
	  // calculate Ntot = sum_i(Ni) number of signal only events from data for this variable
	  double Ntot =NdataVariable-NBGVariable;
	  // get eff*A in bin
	  double effABin=(histo_["efficiency/"+variable][kSig])->GetBinContent(bin);
	  // calculate cross section for this bin
	  double xSecBin = (Ni/Ntot) * (sigFracBin/sigFracVariable) * (effIncl/effABin) * (1/binwidth);
	  // calculate complete error
	  double Nierror=sqrt(NdataBin); // the error of NBG is treated as systematic variation
	  double sumErrorNjSquare = 0; // sum(j,i!=j) [ errorNj*errorNj ]
	  for(int bin2=1; bin2<=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetNbinsX()+1; ++bin2){
	    if(bin2!=bin){    
	      double errorNj = sqrt( histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinContent(bin2) 
				     * (histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetBinWidth(bin2)) );   
	      // remark: Nj=NdataBinj-NBGBinj, as above, the error of NBGBinj is treated 
	      // as systematic variation and therefore not incorporated here
	      sumErrorNjSquare += errorNj*errorNj;
	    }
	  }
	  // calculate complete error
	  double xSecBinError = sqrt( ((Ntot - Ni)*(Ntot - Ni)*Nierror*Nierror) + (Ni*Ni*sumErrorNjSquare) );
	  xSecBinError/= (Ntot*Ntot);
	  xSecBinError*= (effIncl/effABin) * (sigFracBin/sigFracVariable) * (1/binwidth);
	  // check for ambiguities
	  if(!isfinite(xSecBin)||isnan(xSecBin)){
	    if(verbose>1){ 
	      std::cout << "bin " << bin << " in ";
	      std::cout << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable;
	      std::cout << "gives no valid cross section!" << std::endl;
	      std::cout << "!will set this bin to 0!" << std::endl;
	    }
	    xSecBin=0;	 
	    xSecBinError=0;
	  }
	  // set underflow/overflow to 0
	  // because of infinite binwidth
	  if(bin==0||bin==histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->GetNbinsX()+1){
	    xSecBinError=0;
	    xSecBin=0;
	  }
	  // check calculation for one example bin
	  TString testVar  = "ttbarPt";
	  TString testVar2 = "ttbarMass";
	  int testBin=0;
	  if(verbose>0&&xSecBin!=0&&(variable==testVar||variable==testVar2)&&bin==testBin){
	    std::cout << std::endl << "differential normalized xSec ";
	    std::cout << variable+" bin " << bin << ":" << std::endl;
	    std::cout << "binwidth:             " << std::setprecision(3) << std::fixed << binwidth << std::endl;
	    std::cout << "Ndata(bin):           " << std::setprecision(3) << std::fixed << NdataBin << std::endl;
	    std::cout << "NBG(bin):             " << std::setprecision(3) << std::fixed << NBGBin << std::endl;
	    std::cout << "eff*A(bin):           " << std::setprecision(3) << std::fixed << effABin << std::endl;
	    std::cout << "ttbar sig frac(bin):  " << std::setprecision(3) << std::fixed << sigFracBin << std::endl;
	    std::cout << "Ndata(incl):          " << std::setprecision(3) << std::fixed << NdataVariable << std::endl;
	    std::cout << "NBG(incl):            " << std::setprecision(3) << std::fixed << NBGVariable << std::endl;
	    std::cout << "eff*A(incl):          " << std::setprecision(3) << std::fixed << effIncl << std::endl;
	    std::cout << "ttbar sig frac(incl): " << std::setprecision(3) << std::fixed << sigFracVariable << std::endl;
	    std::cout << "result:               " << std::setprecision(3) << std::fixed << xSecBin << std::endl;
	  }
	  // apply cross section and error to plot
	  histo_[xSec][kData]->SetBinContent(bin, xSecBin     );
	  histo_[xSec][kData]->SetBinError  (bin, xSecBinError);
	  //saveToRootFile(outputFileName, histo_[xSec][kData], true, verbose-1,"xSec/NormXSec");
	}
      }
      if(verbose>1){
	std::cout << "integral diff.norm.data distribution (" << variable << "): ";
	std::cout << getInclusiveXSec(histo_[xSec][kData],verbose-1) << std::endl;
      }
      // b) differential normalized XSec from Signal(MC prediction)
      histo_[xSec][kSig]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig]->Clone());
      DivideYieldByEfficiencyAndLumi(histo_[xSec][kSig], (TH1F*)(histo_["efficiency/"+variable][kSig]), luminosity, 0);
      double XSecInclTheoPS= getInclusiveXSec(histo_[xSec][kSig], verbose-1);
      // exclude underflow and overflow bins because they are negligible and treated wrong
      XSecInclTheoPS-=histo_[xSec][kSig]->GetBinContent(0);
      XSecInclTheoPS-=histo_[xSec][kSig]->GetBinContent(histo_[xSec][kSig]->GetNbinsX()+1);
      histo_[xSec][kSig]->Scale(1/(XSecInclTheoPS));
      if(verbose>1){
	std::cout << "integral diff.norm.MC distribution (" << variable << "): ";
	std::cout << getInclusiveXSec(histo_[xSec][kSig],verbose-1) << std::endl;
      }
      // add plot to list of plots
      plotList_.push_back(xSec);
      // add axis configuration
      unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
      TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
      recoAxisLabel.ReplaceAll("KinFit ","");
      //    axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{1}{#sigma}"+" #frac{d#sigma}{d"+label+"} "+label2+" (t#bar{t}#rightarrow #mu prompt"+")/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{1}{#sigma}"+" #frac{d#sigma}{d"+label+"} "+label2+"/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
      if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e");
      else if(decayChannel=="combined") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "l");
      // configure xSec plot histo style
      histogramStyle(*histo_[xSec][kData], kData, false);
      histogramStyle(*histo_[xSec][kSig ], kSig , false);
      if(compare){
	histo_[xSec][kData]->SetLineColor  (38);
	histo_[xSec][kData]->SetMarkerColor(38);
	histo_[xSec][kData]->SetMarkerStyle(33);
	histo_[xSec][kSig ]->SetLineColor  (kOrange);
	histo_[xSec][kSig ]->SetMarkerColor(kOrange);
      }
      // restrict axis
      setXAxisRange(histo_[xSec][kData],variable);
      setXAxisRange(histo_[xSec][kSig ],variable);
      if(!compare) ++NXSec;
    }
  }
  
  // =======================================
  //  Use SVD unfolding for cross sections
  // =======================================
  if(SVDunfold){
    // loop all variables
    for(unsigned int number=0; number<xSecVariables_.size(); ++number){
      // prepare style for xSec plots
      TString variable=xSecVariables_[number];
      bool verboseTest = variable.Contains(testQuantity) ? true : false;
      TString labelNorm =getStringEntry(xSecLabel_[number],1);
      TString label2Norm=getStringEntry(xSecLabel_[number],2);
      if(label2Norm=="[GeV]")  label2Norm="#left[GeV^{-1}#right]";
      else if(label2Norm==" ") label2Norm="";
      TString label =getStringEntry(xSecLabel_[number],1);
      TString label2=getStringEntry(xSecLabel_[number],2);
      if(label2=="[GeV]")  label2="#left[pb GeV^{-1}#right]";
      else if(label2==" ") label2="[pb]";
      TString xSec ="xSec/"+variable;
      TString xSecNorm ="xSecNorm/"+variable;

      // calculate cross section within the chosen PS
      // used for normalization of diff. norm. cross sections
      // NOTE: these cross sections do not include BR corrections!
      TString recpath="analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable;
      TString genpath="analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable;
      // data 
      double Ndataplain= histo_[recpath][kData]->Integral(0, histo_[recpath][kData]->GetNbinsX()+1);
      // inclusive eff*acceptance (PS)
      double NgenSigPS= histo_[genpath][kSig]->Integral(0, histo_[genpath][kSig]->GetNbinsX()+1);
      double NrecSig  = histo_[recpath][kSig]->Integral(0, histo_[recpath][kSig]->GetNbinsX()+1);
      double effAPSforNorm=NrecSig/NgenSigPS;
      // ttbar signal fraction
      double NrecttbarBG = histo_[recpath][kBkg]->Integral(0, histo_[recpath][kBkg]->GetNbinsX()+1);
      double signalFraction=NrecSig/ (NrecSig + NrecttbarBG);
      // non ttbar BG 
      double NrecMCBG = histo_[recpath][kAllMC]->Integral(0, histo_[recpath][kAllMC]->GetNbinsX()+1);
      NrecMCBG-=NrecttbarBG;
      // cross section 
      double xSecPSforNorm= (Ndataplain - NrecMCBG)*signalFraction / (effAPSforNorm * luminosity);
      if(verboseTest||verbose>1){
	std::cout << "calculate incl. xSec PS (without unfolding)" << std::endl;
	std::cout << "Ndataplain :" << Ndataplain << std::endl;
	std::cout << "NrecMCBG :"   << NrecMCBG   << std::endl;
	std::cout << "signalFraction :" << signalFraction << std::endl;
	std::cout << "effAPSforNorm: " << effAPSforNorm << std::endl;
	std::cout << "luminosity: " << luminosity << std::endl;
	std::cout << "xSecPSforNorm: " << xSecPSforNorm << std::endl;
      }
      if(verboseTest||verbose>0){
	std::cout << "inclusive cross section (PS) for normalization: " << xSecPSforNorm << std::endl;
      }
      // ==========================================
      //  Check availablility of input histos
      // ==========================================
      bool calculateXSec=true;
      // check if BG reco plot is available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable, kAllMC)){
	std::cout << "warning: plot " << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable << " not found for BG" << std::endl;
	calculateXSec=false;
      }
      // check if ttbar BG reco plot is available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable, kBkg)){
	std::cout << "warning: plot " << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable << " not found for BG" << std::endl;
	calculateXSec=false;
      }
      // check if data event yield is available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable, kData)){
	std::cout << "warning: plot " << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable << " not found for data" << std::endl;
	calculateXSec=false;
      }
      // check if signal reco plot is available
      if(!plotExists(histo_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable, kSig)){
	std::cout << "warning: plot " << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable << " not found for signal" << std::endl;
	calculateXSec=false;
      }
      // check if signal gen level plot is available
      if(!plotExists(histo_, "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable, kSig)){
	std::cout << "warning: plot " << "analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable << " not found for signal" << std::endl;
	calculateXSec=false;
      }
      // check if response matrix plot is available
      if(!plotExists(histo2_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable+"_", kSig)){
	std::cout << "warning: plot " << "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable+"_" << " not found for signal" << std::endl;
	calculateXSec=false;
      }
      // =====================
      //  Configure labels
      // =====================
      // particle
      TString particle="";
      variable.Contains("ttbar") ? particle="t#bar{t}" : (variable.Contains("lep") ? particle="lepton" : (variable.Contains("top") ? particle="t" : (variable.Contains("bq") ? particle="b" : ( variable.Contains("bbbar") ? particle="b#bar{b}" : particle="unknown") ) ) );
    
      TString particleTex=particle;
      particle.ReplaceAll("#","");
      particle.ReplaceAll("{","");
      particle.ReplaceAll("}","");
      // quantity
      TString quantity="";
      variable.Contains("Mass") ? quantity="mass" : (variable.Contains("Pt") ? quantity="p_{t}" : (variable.Contains("Eta") ? quantity="#eta" : (variable.Contains("Y") ? quantity="y" : quantity="unknown") ) );
      TString quantityTex=xSecLabelName(variable);
      quantity.ReplaceAll("#","");
      quantity.ReplaceAll("{","");
      quantity.ReplaceAll("}","");
      // decay channel
      TString channelTex="#mu + Jets";
      if(decayChannel=="electron") channelTex="e + Jets";
      if(decayChannel=="combined") channelTex="e/#mu + Jets";
      TString special="";
      // ==========================================
      //  Configure the unfolding details
      // ==========================================
      // list of systematics
      TString sysList=sysLabel(systematicVariation);
      //TString sysList[]={sysLabel(systematicVariation)};
      // number of systematic samples to unfold 
      int numSys=0;
      // Regularization parameter
      // Regularization Modus 
      //         0 means: Default setting. Same as 2
      //         1 means: Bin by Bin Unfolding
      //         2 means: SVD Unfolding, minimal regularization.
      //              Regularization will be done with a k value
      //              and the k value will be set to the number of
      //              bins (not counting OF). Parameter 'regPar'
      //              is ignored.
      //         3 means: SVD unfolding, regularization according
      //              to the parameter 'regPar'. Whether this one
      //              is taken as k value or tau value, can be 
      //              specified with the second digit of this
      //              parameter. (default)
      //         4 means: SVD unfolding according to the file
      //              given in 'regFile'. Whether this one
      //              is taken as k value or tau value, can be 
      //              specified with the second digit of this
      //              parameter.Parameter 'regPar'
      //              is ignored.
      int regMode=3;
      // continuous or discrete unfolding?
      //         0 means: Default setting, same as 1
      //         1 means: interpret 'regPar' as k value (default)
      //         2 means: interpret 'regPar' as tau value
      int unfoldWithParameter=2;
      if(systematicVariation!=sysNo) regMode=3;
      //    PRE-WEIGHTING
      //         0 means: Default value, same as 1
      //         1 means: No preweighting of MC is performed (default)
      //         2 means: MC is reweighted to unfolded data (1 iteration)
      //         i means: MC is reweighted to unfolded data ((i-1) iteration), i<=9 (up to 8 iterations)
      int unfPreWeighting=iterative ? 6 : 0;      
      TString unfPreWeightingStr = "";
      if(unfPreWeighting>1) unfPreWeightingStr+="PreWeight"+getTStringFromInt(unfPreWeighting-1)+"Iter";
      // produce scan plots? 
      //         0 means: Default setting, same as 2
      //         1 means: no scan for optimal tau is perforemd,
      //         2 means: tau scan  (default).
      //              A scan for the optimal tau parameter is
      //              performed. The scan is performed around
      //              a "center value" to be specified in 'regpar'.
      //              Whether this one is taken
      //              as k value or tau value, can be 
      //              specified with the second digit of this
      //              parameter.
      //              Note: The scan may take a while!
      int scan =1;
      if(redetermineopttau) scan=2;
      else scan=1;
      bool scanPlots=false;
      if(regMode<3) scan =1;
      if(scan==2) scanPlots=true;
      // print unfolding eps and ps plots?
      // 0 means: no
      // 1 means: only ps file containing all plots
      // 2 means: both ps and individual eps
      int printUnfPlots = 2;
      
      // output files
      TString rootFile="";
      TString psFile="";
      TString epsFile="";
      TString regFile="";
      //         0 means: Default value, same as 4
      //         1 means: no plotting at all
      //         2 means: standard plots
      //         3 means: standard plots + k scan plots
      //         4 means: standard plots + k scan plots + tau scan plots (default)
      int plotting=1;
      if(redetermineopttau) plotting=4;
      if((save||redetermineopttau)&&systematicVariation==sysNo){ 
	// output file: content
	if(plotting==1) plotting+=1;
	if(scan==2&&scanPlots) plotting+=2; // k and tau scan plots
	// output files: labels
	rootFile=outputFolder+"unfolding/diffXSecUnfoldTopSemi";
	if(decayChannel=="muon"    ) rootFile+="Mu";
	else if(decayChannel=="electron") rootFile+="Elec";
	else if(decayChannel=="combined") rootFile+="Lep";
	rootFile+=closureLabel+"RegTest"+variable+LV+PS+unfPreWeightingStr+".root";
	psFile =outputFolder+"unfolding/unfolding"+closureLabel+"RegTest"+variable+LV+PS+unfPreWeightingStr;
	epsFile=outputFolder+"unfolding/unfolding"+closureLabel+"RegTest"+variable+LV+PS+unfPreWeightingStr+".eps";
	if(scan==2) psFile+="Scan";
	psFile+=".ps";
	//if(regMode>2) regFile=outputFolder+"unfolding/optimalSVDRegularization.txt";
	// if ps and/or eps files not desired, set back to ""
	if(printUnfPlots!=1 && printUnfPlots!=2) psFile="";
	if(printUnfPlots!=2 )                    epsFile="";
      }
      
      TString txtfile=outputFolder+"unfolding/diffXSecUnfoldDetailsTopSemi";
      if(decayChannel=="muon"    ) txtfile+="Mu";
      else if(decayChannel=="electron") txtfile+="Elec";
      else if(decayChannel=="combined") txtfile+="Lep";
      txtfile+=closureLabel+"RegTest"+dataSample+LV+PS+unfPreWeightingStr+".txt";
      // save unfolding plots in rootfile?
      //         0 means: Default value, same as 1
      //         1 means: no root file will be written (default)
      //         2 means: standard plots to root file  
      int doRootFile=(rootFile=="" ? 1 : 2);
      // Combine all options in one Steering Parameter
      TString steering="";
      TString steering2=""; // k paramter unfolding version
      // Steering options in parameter 'steering' 
      //     (1) REGMODE, see above
      //         2 means: SVD Unfolding, minimal regularization.
      //              Regularization will be done with a k value
      //              and the k value will be set to the number of
      //              bins (not counting OF). Parameter 'regPar'
      //              is ignored.
      //         3 means: SVD unfolding, regularization according
      //              to the parameter 'regPar'. Whether this one
      //              is taken as k value or tau value, can be 
      //              specified with the second digit of this
      //              parameter. (default)
      steering=getTStringFromInt(regMode)+steering;
      steering2=getTStringFromInt(regMode)+steering2; 
      //     (2) REGULARIZATION PARAMETER, see above
      steering=getTStringFromInt(unfoldWithParameter)+steering; //2 means: interpret 'regPar' as tau value
      steering2=getTStringFromInt(1)+steering2; //1 means: interpret 'regPar' as k value
      //     (3) SCAN 
      steering=getTStringFromInt(scan)+steering;
      steering2=getTStringFromInt(scan)+steering2;
      //     (4) OUTPUT PS CONTENT, see above
      steering=getTStringFromInt(plotting)+steering;
      steering2=getTStringFromInt(plotting)+steering2;
      //     (5) ROOT FILE 
      steering=getTStringFromInt(doRootFile)+steering;
      steering2=getTStringFromInt(doRootFile)+steering2;
      //     (6) TEXT FILE
      //         0 means: Default value, same as 2
      //         1 means: no text file is written
      //         2 means: text file with histogram entries is written (default)
      //         3 means: text file is written if (syst==sysNo)
      int  doTextFile = 1;
      steering=getTStringFromInt(doTextFile)+steering;
      steering2=getTStringFromInt(doTextFile)+steering2;
      //     (7) VERBOSITY
      //         0 means: Default value, same as 2
      //         1 means: no output at all
      //         2 means: standard output (default)
      //         3 means: debugging output
      int verbosity=0;
      if(verbose==0) verbosity=1;
      if(verbose==1) verbosity=2;
      if(redetermineopttau) verbosity=2;
      if(verbose>1) verbosity=3;
      steering=getTStringFromInt(verbosity)+steering;
      steering2=getTStringFromInt(verbosity)+steering2;
      //     (8)  SCANPOINTS
      //          0 means: Default value, same as 3
      //          1 means: 5 scan points
      //          2 means: 25 scan points
      //          3 means: 125 scan points (default)
      //          4 means: 625 scan points
      int scanpoints= (scan==2 ? 3 : 0);
      scanpoints=1; // FIXME: fast tauscan results
      steering=getTStringFromInt(scanpoints)+steering;
      steering2=getTStringFromInt(scanpoints)+steering2;
      //     (9)  SCANRANGE
      //          0 means: Default value, same as 2
      //          1 means: Tau+ / Tau- = 100  
      //          2 means: Tau+ / Tau- = 10000
      //          3 means: Tau+ / Tau- = 1000000
      int scanrange = (scan==2 ? 3 : 0);
      steering=getTStringFromInt(scanrange)+steering;
      steering2=getTStringFromInt(scanrange)+steering2;
      //     (10) LOWER SIDE BIN
      //         0 means: Default value, same as 3
      //         1 means: Regard as regular bin (not encouraged!)
      //         2 means: Regard as regular bin, keep bin content fixed to MC (not encouraged!)
      //         3 means: Cut away on Rec Level, unfold to Gen Level (default)
      //         4 means: Cut away on Rec Level, unfold to Gen Level, keep bin content fixed to MC 
      //         5 means: Cut away on Rec Level, ignore on Gen Level (not encouraged!) 
      int lowsidebin=0;
      steering=getTStringFromInt(lowsidebin)+steering;
      steering2=getTStringFromInt(lowsidebin)+steering2;
      //     (11) UPPER SIDE BIN 
      //         0 means: Default value, same as 3
      //         1 means: Regard as regular bin (not encouraged!)
      //         2 means: Regard as regular bin, keep bin content fixed to MC (not encouraged!)
      //         3 means: Cut away on Rec Level, unfold to Gen Level (default)
      //         4 means: Cut away on Rec Level, unfold to Gen Level, keep bin content fixed to MC 
      //         5 means: Cut away on Rec Level, ignore on Gen Level (not encouraged!) 
      int upsidebin=0;
      steering=getTStringFromInt(upsidebin)+steering;
      steering2=getTStringFromInt(upsidebin)+steering2;
      //     (12) ORIENTATION OF RESPONSE MATRIX
      //          0 means: Default value, same as 2
      //          1 means: Do not transpose input response matrix during rebinning
      //          2 means: Do transpose input response matrix during rebinning (default)
      int matrixori=2;
      steering=getTStringFromInt(matrixori)+steering;
      steering2=getTStringFromInt(matrixori)+steering2;
      //    (13) NORMALIZATION (13. digit from right)
      //         0 means: Default value, same as 1
      //         1 means: Extrinsic Normalization (default)
      //             This means, the normalization is done with a global inclusive cross
      //             section which is calculated from the parameters totalDataEvents, 
      //             totalBgrEvents, totalTtBgrEvents, totalRecEvents and totalGenEvents.
      //             Analytic Error Propagation.
      //         2 means: Intrinsic Normalization.
      //             Each unfolded distribution is normalized with its integral. 
      //             Analytic Error Propagation.
      //         3 means: Intrinsic Normalization.
      //             Each unfolded distribution is normalized with its integral. 
      //             Error Propagation on the basis of pseudo experiments
      int normalizeUnfPlot=3;
      steering=getTStringFromInt(normalizeUnfPlot)+steering;
      steering2=getTStringFromInt(normalizeUnfPlot)+steering2;
      // NEW: !!!
      //    (14) CLOSURE TEST (14. digit from right)
      //         0 means: Default value, same as 1
      //         1 means: No closure test, unfold real data (default)
      //         2 means: Closure test, unfold the reconstructed MC 
      int closureTest=0;
      steering=getTStringFromInt(closureTest)+steering;
      steering2=getTStringFromInt(closureTest)+steering2;
      //    (15) PRE-WEIGHTING (15. digit from right)
      //         0 means: Default value, same as 1
      //         1 means: No preweighting of MC is performed (default)
      //         2 means: MC is reweighted to unfolded data (1 iteration)
      //         i means: MC is reweighted to unfolded data ((i-1) iteration), i<=9 (up to 8 iterations)
      steering=getTStringFromInt(unfPreWeighting)+steering;
      steering2=getTStringFromInt(unfPreWeighting)+steering2;
      //    (16) BACKGROUND REDUCTION (16. digit from right)
      //         0 means: Default value, same as 1
      //         1 means: Be non-forgiving
      //         2 means: If Background>Data, set Data to zero.
      int unfBkgReduction=0;
      steering=getTStringFromInt(unfBkgReduction)+steering;
      steering2=getTStringFromInt(unfBkgReduction)+steering2;
      //    (17) CONTROL PLOT STYLE (17. digit from right)
      //         0 means: Default value, same as 1
      //         1 means: default style (e.g. show all bins incl. UF/OF, show vertical lines for tau scan etc.)
      //         2 means: Uni HH style (e.g. show only bins of measurement, no vertical lines etc.)
      int unfPlotStyle=2;
      steering=getTStringFromInt(unfPlotStyle)+steering;
      steering2=getTStringFromInt(unfPlotStyle)+steering2;

      // ==============
      //  Get binning
      // ==============

      // a) filter relevant bins
      std::vector<double> relevantBins_(binning_[variable]);
      // b) calculate number of considered bins 
      int unfoldbins=relevantBins_.size()-1; // NB: N(bins)=NbinEdges-1 
      // c) refill bin edges to array
      double bins[relevantBins_.size()];
      for(unsigned int bin=0; bin<relevantBins_.size(); ++bin){
	bins[bin]=relevantBins_[bin]; 
	if(verbose>1) std::cout << "bin " << bin << ": " << bins[bin]<< std::endl;
      }
      if(verboseTest||verbose>1) std::cout << "bins used for unfolding: " << unfoldbins << std::endl;

      double * NdataTot = new double[1];
      NdataTot[0] = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData]->Integral(0,histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData]->GetNbinsX()+1);
      double * NBGTot= new double[1];
      NBGTot[0] = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kAllMC]->Integral(0,histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kAllMC]->GetNbinsX()+1);
      double * NttBG= new double[1];
      NttBG[0] = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kBkg]->Integral(0,histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kBkg]->GetNbinsX()+1);
      double * NttSG= new double[1];
      NttSG[0] = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kSig]->Integral(0,histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kSig]->GetNbinsX()+1);
      double * NttSGgen= new double[1];
      NttSGgen[0] = histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable][kSig]->Integral(0,histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable][kSig]->GetNbinsX()+1);
      
      // ============================
      //  Use unfolding machine
      // ============================
      if(calculateXSec){
	for(int bin=1; bin<=histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetNbinsX(); ++bin){
	  double NumSig  = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig  ]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSig  ]->GetBinWidth(bin);
	  double NumBkg  = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg  ]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kBkg  ]->GetBinWidth(bin);     
	  double NumQCD  = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kQCD  ]->GetBinWidth(bin);
	  double NumVV   = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kDiBos]->GetBinWidth(bin);
	  double NumStop = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kSTop ]->GetBinWidth(bin);
	  double NumWjet = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kWjets]->GetBinWidth(bin);
	  double NumZjet = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets]->GetBinContent(bin)*histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kZjets]->GetBinWidth(bin);
	  double dataVal = histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinContent(bin)* histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinWidth(bin);
	  double allBGVal=NumStop+NumWjet+NumZjet+NumVV+NumQCD;
	  if(allBGVal>dataVal){
	    std::cout << "WARNING: data=" << dataVal << "<BG=" << allBGVal << " in bin " << bin;
	    std::cout << "("<< histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinLowEdge(bin) << ".." << histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinLowEdge(bin+1) << ")"<< std::endl;
	    std::cout << " for " << variable << "(" << sysLabel(systematicVariation) << ") !!!" << std::endl;
	    std::cout << "(now fixed by setting data=MC)" << std::endl;
	    int binInRaw= histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData]->FindBin(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinCenter(bin));
	    histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData]->SetBinContent(binInRaw, 1.01*(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData]->GetBinContent(binInRaw)+(allBGVal-dataVal)));
	  }
	  if(allBGVal>dataVal||verboseTest||verbose>3){
	    // N(events) output for each bin and variable
	    std::cout << variable << " (bin" << bin << ")= (" <<  histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinLowEdge(bin) << "," << histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable ][kData]->GetBinLowEdge(bin+1) << ")" << std::endl;
	    std::cout << "     N(data): "   << dataVal << std::endl; 
	    std::cout << "     N(BG all): " << allBGVal << std::endl;
	    std::cout << "     N(Sig): "    << NumSig  << std::endl;
	    std::cout << "     N(Bkg): "    << NumBkg  << std::endl;
	    std::cout << "     N(Stop): "   << NumStop << std::endl;
	    std::cout << "     N(Wjets): "  << NumWjet << std::endl;
	    std::cout << "     N(Zjets): "  << NumZjet << std::endl;
	    std::cout << "     N(VV): "     << NumVV   << std::endl;
	    std::cout << "     N(QCD): "    << NumQCD  << std::endl;
 	  }
	}
	
	// collect reg Parameters to calculate cross section for
	double regPar=regParameter(variable, decayChannel, verbose, extrapolate, true, hadron, closureTestSpecifier, (addSel=="ProbSel" ? true : false) );
	stdTau_[variable]=regPar;
	// tau parameters
	if(test=="tau"){
	  tau_.push_back(3.);	
	  tau_.push_back(5.);
	}
	else if(test=="mix"){
	  tau_.push_back(0.*regPar);	  
	  tau_.push_back(0.1*regPar);
	}
	tau_.push_back(regPar);
	int stdindex=(int)tau_.size(); // used as indicator entries >=kindex == k unfolding
	if(test=="tau"){
	  tau_.push_back(10.);
	  tau_.push_back(15.);
	  tau_.push_back(25.);
	}
        else if(test=="mix"){
	  tau_.push_back(5.0*regPar);
	  tau_.push_back(10.*regPar);
	  tau_.push_back(50.*regPar);
	}
	kindex=(int)tau_.size(); // used as indicator entries >=kindex == k unfolding
	// k parameters
	if(test=="k"){
	  tau_.push_back(1); 
	  tau_.push_back(2); 
	}
	if(test!="tau"){
	  tau_.push_back(3);
	  tau_.push_back(4);
	  tau_.push_back(5);
	  tau_.push_back(6);
	}
        if(test=="k"){
          tau_.push_back(7);
          tau_.push_back(8);
        }
	// colors
	tauColor_.push_back(kCyan+1 );
	tauColor_.push_back(kBlue   );
	//tauColor_.push_back(kBlack  );
	tauColor_.push_back(kGreen+1);
	tauColor_.push_back(kGreen+3);
	tauColor_.push_back(kAzure+8);
	tauColor_.push_back(kOrange+3 ); // k parameter
	tauColor_.push_back(kRed      ); // k parameter
	tauColor_.push_back(kMagenta+2); // k parameter
	tauColor_.push_back(kOrange+8); // k parameter
	// marker styles
	tauMarker_.push_back(33);
	tauMarker_.push_back(22);	
	//tauMarker_.push_back(29);
	tauMarker_.push_back(23);
	tauMarker_.push_back(34);
	tauMarker_.push_back(20);
	tauMarker_.push_back(24); // k parameter
	tauMarker_.push_back(30); // k parameter
	tauMarker_.push_back(32); // k parameter
	tauMarker_.push_back(25); // k parameter
	// default reg parameter: black star
	int defCol [1] ={kBlack};
	int defMar [1] ={29};
	tauColor_ .insert(tauColor_ .begin()+stdindex-1, defCol, defCol+1);
        tauMarker_.insert(tauMarker_.begin()+stdindex-1, defMar, defMar+1);
	
	for(int tau=0; tau<(int)tau_.size(); ++tau){
	  // to avoid overwriting, add the used tau parameter to all output files
	  TString tauVal= tau<kindex ? "Tau"+getTStringFromDouble(tau_[tau],2) : "K"+getTStringFromDouble(tau_[tau],0);;
	  tauVal.ReplaceAll(".","p");
	  TString rootFileMod=rootFile;
	  rootFileMod.ReplaceAll(".root", tauVal+".root");
	  TString psFileMod=psFile;
	  psFileMod.ReplaceAll(".ps", tauVal+".ps");
	  TString epsFileMod=epsFile;
	  epsFileMod.ReplaceAll(".eps", tauVal+".eps");
	  TString txtfileMod=txtfile;
	  txtfileMod.ReplaceAll(".txt", tauVal+".txt");

	  TH1D* unfoldedData    =new TH1D();
	  TH1D* unfoldedDataNorm=new TH1D();
	  TopSVDFunctions::SVD_Unfold(
	  // ==============
	  //  HISTOS
	  // ==============				    
	  // Data Input (RAW Data including the background)
	  (TH1D*)histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kData],                  
	  // Background (will be substracted from data)
	  (TH1D*)histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kAllMC],
	  // ttbar background only (used to calculate a ttbar signal 
	  // fraction instead of subtracting the yield which depends 
	  // on the inclusive ttbar cross section.) 
	  // Note: if 0 pointer is handed over 
	  (TH1D*)histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kBkg],
	  // Generated ttbar SG MC
	  (TH1D*)histo_["analyzeTopPartonLevelKinematics"+PS+sysInputGenFolderExtension+"/raw"+variable][kSig],
	  // Reconstructed ttbar SG MC
	  (TH1D*)histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable][kSig],
	  // Response Matrix 
	  (TH2D*)histo2_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/raw"+variable+"_"][kSig],
	  // Returned: Unfolded Distribution              
	  unfoldedData,
	  // Returned: Normalized Unfolded Distribution              
	  unfoldedDataNorm,
	  // ==============
	  //  NUMBERS
	  // ==============	
	  // For Normalization: Total Data Events
	  // If set to zero, will be taken from the integral of dataInputHist
	  NdataTot,
	  // For Normalization: Total Reconstructed all BG MC Events
	  // If set to NULL, will be taken from the integral of bgrInputHist
	  NBGTot,
	  // For Normalization: Total Reconstructed ttbar BG MC Events
	  // If set to NULL, will be taken from the integral of ttbgrInputHist
	  NttBG,
	  // For Normalization: Total Reconstructed ttbar SG MC Events
	  // If set to NULL, will be taken from the integral of recInputHist
	  NttSG,
	  // For Normalization: Total Generated ttbar SG MC Events
	  // If set to NULL, will be taken from the integral of genInputHist
	  NttSGgen,
	  // ==============
	  //  BINS
	  // ==============	
	  // Binning for the unfolding
	  bins, 
	  // Number of bins for unfolding (not counting OF bins !!!)
	  unfoldbins, 
	  // ==============
	  //  OPTIONS
	  // ==============
	  // Regularization parameter
	  tau_[tau],
	  // steering parameter (as defined above)
	  tau<kindex ? steering : steering2, // difference: information whether to use k or tau regularization
	  // Specify the number of systematic samples to unfold 
	  numSys,
	  // =====================
	  //  TEX(T) AND LABELS
	  // =====================                         
	  // Specify Name for the Channel ("mumu", "emu", "ee" ...)
	  (TString)decayChannel,  
	  // Specify Name for the Physics Object ("Top", "Jets", "Leptons")      
	  particle,
	  // Specify Name for the Quantity ("Eta", "Pt", or "Mass");      
	  quantity,
	  // Specify Name for special run of unfolding
	  special, 
	  // Array of Names for the different systematics
	  // if you run only over the nominal sample, provide NULL
	  sysList,
	  // Nicely formatted name for the channel
	  channelTex,
	  // Nicely formatted name for the physics object 
	  particleTex,
	  // Nicely formatted name for ithe quantity
	  quantityTex,
	  // Nicely formatted name indicating some special condition 
	  special,
	  // Array of Names for the different systematics
	  // if you run only over the nominal sample, provide NULL
	  sysList,
	  // ==============
	  //  OUTPUT
	  // ==============  
	  // If specified, plots will be saved in ROOT File
	  rootFileMod,
	  // If specified, plots will be saved in PS File
	  psFileMod,
	  // If specified, plots will be saved in EPS Files 
	  epsFileMod,
	  // The optimal Reg Parameters will be written to this file.
	  // The file will NOT be overwritten, rather the result will be appended.
	  // The following data will be saved in this order: 
	  // (1) the optimal tau, (2) the two nearest k values,
	  // (3) the k value from the d value method
	  // (4) the number of bins including side bins
	  txtfileMod,
	  // If specified, optimal Reg Parameters will be 
	  // written to this file. The file will NOT be 
	  // overwritten, rather the result will be appended
	  // The following data will be saved in this order: 
	  // (1) the optimal tau, (2) the two nearest k values,
	  // (3) the k value from the d value method
	  // (4) the number of bins including side bins
	  regFile
	  );
	  // ========================================================
	  //  Remaining steps for cross section calculation
	  // ========================================================
	  // use unfolded event yield as input
	  histo_[xSec][kData+42   ]=(TH1F*)(unfoldedData->Clone(variable+"kData"));
	  histo_[xSec][kData+42+42]=(TH1F*)(unfoldedDataNorm->Clone(variable+"kData"));
	  // use reco yield plot clone to get correct and complete binning 
	  histo_[xSec    ][kData]=(TH1F*)(histo_["analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable][kData]->Clone(variable+"kData"));
	  histo_[xSec    ][kData]->Reset("icesm");
	  histo_[xSecNorm][kData]=(TH1F*)(histo_[xSec][kData]->Clone());
	  // remove additional bins with 0 width
	  double *newBinLowEdges = new double[binning_[variable].size()];
	  newBinLowEdges = &binning_[variable].front();
	  histo_[xSec    ][kData] = (TH1F*) histo_[xSec    ][kData]->Rebin(binning_[variable].size()-1, histo_[xSec    ][kData]->GetName(), newBinLowEdges); 
	  histo_[xSecNorm][kData] = (TH1F*) histo_[xSecNorm][kData]->Rebin(binning_[variable].size()-1, histo_[xSecNorm][kData]->GetName(), newBinLowEdges); 
	  // remove additional OF/UF bins from SVD unfolded plots
	  if(verbose>1){ 
	    std::cout << std::endl << variable << ":" << std::endl;
	    for(int bin=0; bin<=histo_[xSec][kData+42]->GetNbinsX()+1; ++bin){
	      std::cout << "bin " << bin;
	      std::cout << " (" << histo_[xSec][kData+42]->GetBinLowEdge(bin) << " .. ";
	      //std::cout << (bin==histo_[xSec][kData+42]->GetNbinsX()+1 ? 999999 : histo_[xSec][kData+42]->GetBinLowEdge(bin+1));
	      std::cout << histo_[xSec][kData+42]->GetBinLowEdge(bin+1);
	      std::cout << "): " << histo_[xSec][kData+42]->GetBinContent(bin) << ", Norm: " << histo_[xSec][kData+42+42]->GetBinContent(bin) << std::endl;
	    }
	  }
	  for(int binSVD=0; binSVD<=histo_[xSec][kData+42]->GetNbinsX()+1; ++binSVD){
	    if(verbose>1) std::cout << "bin " << binSVD;
	    double value=histo_[xSec][kData+42]->GetBinContent(binSVD);
	    double error=histo_[xSec][kData+42]->GetBinError(binSVD);
	    double valueNorm=histo_[xSec][kData+42+42]->GetBinContent(binSVD);
	    double errorNorm=histo_[xSec][kData+42+42]->GetBinError(binSVD);
	    double loweredgeSVD=histo_[xSec][kData+42]->GetBinLowEdge(binSVD);
	    double higheredgeSVD=histo_[xSec][kData+42]->GetBinLowEdge(binSVD+1);
	    // search corresponding final bin
	    bool found=false;
	    for(int bin=1; bin<=histo_[xSec][kData]->GetNbinsX(); ++bin){
	      double loweredge=histo_[xSec][kData]->GetBinLowEdge(bin);
	      double higheredge=histo_[xSec][kData]->GetBinLowEdge(bin+1);
	      if(!found){
		// normal bins
		if(loweredgeSVD>=loweredge&&higheredgeSVD<=higheredge){
		  found=true; 
		  if(verbose>1) std::cout << "-> bin " << bin << std::endl;
		  histo_[xSec][kData]->SetBinContent(bin, histo_[xSec][kData]->GetBinContent(bin)+value);
		  histo_[xSec][kData]->SetBinError(bin, sqrt(histo_[xSec][kData]->GetBinError(bin)*histo_[xSec][kData]->GetBinError(bin)+error*error));
		  histo_[xSecNorm][kData]->SetBinContent(bin, histo_[xSecNorm][kData]->GetBinContent(bin)+valueNorm);
		  histo_[xSecNorm][kData]->SetBinError(bin, sqrt(histo_[xSecNorm][kData]->GetBinError(bin)*histo_[xSecNorm][kData]->GetBinError(bin)+errorNorm*errorNorm));
		}
	      }
	    }
	    // SVD bins outside plot range -> OF/UF
	    if(!found){
	      if(binSVD==0||(loweredgeSVD<=histo_[xSec][kData]->GetBinLowEdge(1)&&higheredgeSVD<=histo_[xSec][kData]->GetBinLowEdge(1))){
		found=true; 
		int bin=0;
		if(verbose>1) std::cout << "-> bin " << bin << std::endl;
		histo_[xSec][kData]->SetBinContent(bin, histo_[xSec][kData]->GetBinContent(bin)+value);
		histo_[xSec][kData]->SetBinError(  bin, sqrt(histo_[xSec][kData]->GetBinError(bin)*histo_[xSec][kData]->GetBinError(bin)+error*error));
		histo_[xSecNorm][kData]->SetBinContent(bin, histo_[xSecNorm][kData]->GetBinContent(bin)+valueNorm);
		histo_[xSecNorm][kData]->SetBinError(  bin, sqrt(histo_[xSecNorm][kData]->GetBinError(bin)*histo_[xSecNorm][kData]->GetBinError(bin)+errorNorm*errorNorm));
	      }
	      else if(binSVD==histo_[xSec][kData+42]->GetNbinsX()+1||(loweredgeSVD>=histo_[xSec][kData]->GetBinLowEdge(histo_[xSec][kData]->GetNbinsX()+1)&&higheredgeSVD>=histo_[xSec][kData]->GetBinLowEdge(histo_[xSec][kData]->GetNbinsX()+1))){
		found=true; 
		int bin=histo_[xSec][kData]->GetNbinsX()+1;
		if(verbose>1) std::cout << "-> bin " << bin << std::endl;
		histo_[xSec][kData]->SetBinContent(bin, histo_[xSec][kData]->GetBinContent(bin)+value);
		histo_[xSec][kData]->SetBinError(  bin, sqrt(histo_[xSec][kData]->GetBinError(bin)*histo_[xSec][kData]->GetBinError(bin)+error*error));
		histo_[xSecNorm][kData]->SetBinContent(bin, histo_[xSecNorm][kData]->GetBinContent(bin)+valueNorm);
		histo_[xSecNorm][kData]->SetBinError(  bin, sqrt(histo_[xSecNorm][kData]->GetBinError(bin)*histo_[xSecNorm][kData]->GetBinError(bin)+errorNorm*errorNorm));
	      }
	    }
	    // make sure there is always a corresponding bin found!
	    if(!found){ 
	      std::cout << std::endl << "ERROR: bin " << binSVD << " (" << loweredgeSVD << " .. " << higheredgeSVD << ") from SVD plot does not match any of the final bins:" << std::endl;
	      for(int i=0; i<(int)binning_[variable].size(); ++i){
		std::cout << binning_[variable][i] << " ";
	      }
	      std::cout << std::endl;
	      exit(0);
	    }
	  }
	  if(verbose>1){
	    // std::cout << std::endl << "->" << std::endl;
	    for(int bin=0; bin<=histo_[xSec][kData]->GetNbinsX()+1; ++bin){
	      std::cout << "bin " << bin;
	      std::cout << " (" << histo_[xSec][kData]->GetBinLowEdge(bin) << " .. ";
	      //std::cout << (bin==histo_[xSec][kData]->GetNbinsX()+1 ? 999999 : histo_[xSec][kData]->GetBinLowEdge(bin+1));
	      std::cout << histo_[xSec][kData]->GetBinLowEdge(bin+1);
	      std::cout << "): " << histo_[xSec][kData]->GetBinContent(bin) << std::endl;
	    }
	  }
	  histo_[xSec][kData]->SetTitle(variable);
	  // divide by binwidth
	  histo_[xSec][kData] = divideByBinwidth(histo_[xSec][kData], verbose-1);
	  // divide by luminosity 
	  histo_[xSec][kData]->Scale(1./(luminosity2));
	  // BR correction
	  if(extrapolate) histo_[xSec][kData]->Scale(1./BRPDG(systematicVariation));
	  // Normalization -> not needed anymore, done in unfolding setup
	  // NB: exclude underflow and overflow bins because they are negligible and treated wrong
	  //histo_[xSecNorm][kData]=(TH1F*)histo_[xSec][kData]->Clone();
	  //double inclXSecPS =getInclusiveXSec(histo_[xSec][kData],verbose-1);
	  //inclXSecPS-=histo_[xSec][kData]->GetBinContent(0);
	  //inclXSecPS-=histo_[xSec][kData]->GetBinContent(histo_[xSec][kData]->GetNbinsX()+1);
	  //if(normToIntegral==true){
	  //  histo_[xSecNorm][kData]->Scale(1./inclXSecPS);
	  //  if(verbose>1) std::cout << "normalizing to integral" << std::endl;
	  //}
	  //else{
	  //  histo_[xSecNorm][kData]->Scale(1./xSecPSforNorm);
	  //  if(verbose>1){
	  //    std::cout << "normalizing to calculated incl. xSec:" << std::endl;
	  //    std::cout << xSecPSforNorm << std::endl;
	  //  }
	  //}
	  // divide Normalized plot by binwidth and set title
	  histo_[xSecNorm][kData]->SetTitle(variable);
	  histo_[xSecNorm][kData] = divideByBinwidth(histo_[xSecNorm][kData], verbose-1);
	  if(verboseTest||verbose>1){ 
	    std::cout << std::endl << variable << ":" << std::endl;
	    for(int bin=0; bin<=histo_[xSecNorm][kData]->GetNbinsX()+1; ++bin){
	      std::cout << "bin " << bin;
	      std::cout << " (" << histo_[xSecNorm][kData]->GetBinLowEdge(bin) << " .. ";
	      std::cout << (bin==histo_[xSecNorm][kData]->GetNbinsX()+1 ? 999999 : histo_[xSecNorm][kData]->GetBinLowEdge(bin+1));
	      std::cout << "): " << histo_[xSecNorm][kData]->GetBinContent(bin) << std::endl;
	    }
	  }
	  if(verboseTest||verbose>0){
	    std::cout << std::endl << variable << std::endl;
	    std::cout << "data preunfolded inclusive abs: " << xSecPSforNorm << std::endl;
	    std::cout << "data unfolded sum abs: "  << getInclusiveXSec(histo_[xSec    ][kData],0) << std::endl;
	    std::cout << "data norm sum: " << getInclusiveXSec(histo_[xSecNorm][kData],0) << std::endl;
	  }
	  // =====================
	  //  Styling issues
	  // =====================
	  // add plot to list of plots - only done once
	  if(tau==0){
	    plotList_.push_back(xSec);
	    plotList_.push_back(xSecNorm);
	    // add axis configuration
	    unsigned int positionOfRecoAxisLabel = positionInVector(plotList_, "analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/"+variable);
	    TString recoAxisLabel =axisLabel_[positionOfRecoAxisLabel];
	    recoAxisLabel.ReplaceAll("KinFit ","");
	    axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{d#sigma}{d"+label+"} "+label2+"/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
	    if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e"  );
	    if(decayChannel=="combined") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "lep");
	    axisLabel_.push_back(""+getStringEntry(recoAxisLabel,1)+"/"+"#frac{1}{#sigma}"+" #frac{d#sigma}{d"+labelNorm+"} "+label2Norm+"/"+getStringEntry(recoAxisLabel,3)+"/"+getStringEntry(recoAxisLabel,4));
	    if(decayChannel=="electron") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "e"  );
	    if(decayChannel=="combined") axisLabel_[axisLabel_.size()-1].ReplaceAll("#mu", "lep");
	    NXSec=NXSec+2;
	  }
	  // configure xSec plot histo style
	  histogramStyle(*histo_[xSecNorm][kData], kData, false);
	  // thicker error bars for comparison
	  if(compare){
	    histo_[xSecNorm][kData]->SetLineWidth(5);
	    histo_[xSec    ][kData]->SetLineWidth(3);
	  }
	  histo_[xSecNorm][kData]->SetMarkerSize (0.1);	  
	  TString tauName=(xSecNorm+Form("%2.2f", tau_[tau]));
	  tauxSecNorm_.push_back((TH1F*)histo_[xSecNorm][kData]->Clone(tauName));
	  tauxSecNorm_[tau]->SetLineStyle  (1);
	  tauxSecNorm_[tau]->SetLineWidth  (1);
	  tauxSecNorm_[tau]->SetLineColor  (tauColor_[tau]);
	  tauxSecNorm_[tau]->SetMarkerColor(tauColor_[tau]);
	  tauxSecNorm_[tau]->SetMarkerStyle(tauMarker_[tau]);
	  tauxSecNorm_[tau]->SetMarkerSize (1.0);
	  TGraphAsymmErrors* tempxSec=new TGraphAsymmErrors(tauxSecNorm_[tau]->GetNbinsX());
	  double binfac=double(tau)/double(tau_.size()+1);
	  for(int bin=0; bin<tauxSecNorm_[tau]->GetNbinsX(); ++bin){
	    double pointXValue=tauxSecNorm_[tau]->GetBinLowEdge(bin+1)+(1.0/double(tau_.size())+binfac)*tauxSecNorm_[tau]->GetBinWidth(bin+1);
	    tempxSec->SetPoint(bin, pointXValue, tauxSecNorm_[tau]->GetBinContent(bin+1));
	  }
	  tempxSec->SetLineStyle  (1);
	  tempxSec->SetLineWidth  (1);
	  tempxSec->SetLineColor  (tauColor_[tau]);
	  tempxSec->SetMarkerColor(tauColor_[tau]);
	  tempxSec->SetMarkerStyle(tauMarker_[tau]);
	  tempxSec->SetMarkerSize (0.8);
	  tauxSecNormM_.push_back((TGraphAsymmErrors*)tempxSec->Clone(tauName+"asym"));
	}

	// ========================================================
	//  Calculate differential XSec from Signal(MC prediction)
	// ========================================================
	// get PS gen event yield/binwidth plots
	histo_[xSec][kSig]=(TH1F*)(histo_["analyzeTopPartonLevelKinematics"+sysInputGenFolderExtension+"/"+variable][kSig]->Clone());
	// divide by lumi
	histo_[xSec][kSig]->Scale(1./luminosity2);
	// BR correction
	if(extrapolate) histo_[xSec][kSig]->Scale(1./BRPDG(systematicVariation));
	// Normalization
	histo_[xSecNorm][kSig]=(TH1F*)(histo_[xSec][kSig]->Clone());
	// NB: exclude underflow and overflow bins because they are mostly negligible and we want to nornmalise to the visible range
	histo_[xSecNorm][kSig]->SetBinContent(0                                    , 0.);
	histo_[xSecNorm][kSig]->SetBinContent(histo_[xSecNorm][kSig]->GetNbinsX()+1, 0.);
	double XSecInclTheoPS= getInclusiveXSec(histo_[xSec][kSig],verbose-1);
	histo_[xSecNorm][kSig]->Scale(1./(XSecInclTheoPS));
	// style
	histogramStyle(*histo_[xSec    ][kSig ], kSig , false);
	histogramStyle(*histo_[xSecNorm][kSig ], kSig , false);
	if(compare){
	  histo_[xSec    ][kSig]->SetLineWidth(3);
	  histo_[xSecNorm][kSig]->SetLineWidth(4);
	}
	setXAxisRange(histo_[xSec    ][kSig],variable);
	setXAxisRange(histo_[xSecNorm][kSig],variable);
	if(verbose>1){
	  std::cout << std::endl << xSec << std::endl;
	  std::cout << "mc abs for normalization: "   << XSecInclTheoPS << std::endl << std::endl;
	  std::cout << "UF: "   << histo_[xSec][kSig]->GetBinContent(0) << std::endl;
	  std::cout << "OF: "   << histo_[xSec][kSig]->GetBinContent(histo_[xSec][kSig]->GetNbinsX()+1) << std::endl;
	  std::cout << "mc abs:"   << getInclusiveXSec(histo_[xSec    ][kSig],2) << std::endl;
	  std::cout << "mc norm: " << getInclusiveXSec(histo_[xSecNorm][kSig],2) << std::endl;
	}
	// =================================================
	//  Additional histos for closure tests
	// =================================================
	if(systematicVariation==sysNo&&variable!="inclusive"&&closureTestSpecifier!=""&&!closureTestSpecifier.Contains("NoDistort")){
	  // get ttbar signal files (with reweighting applied or without for zprime)
	  TString rewfold=groupSpace+inputFolderName+"/";
	  if(zprime=="") rewfold+="ttbarReweight/";
	  TString muReweighted=rewfold+TopFilename(kSig, sysNo, "muon"    );
	  TString elReweighted=rewfold+TopFilename(kSig, sysNo, "electron");
	  // name extension for reweighted files
	  if(zprime==""){
	    muReweighted.ReplaceAll("Summer", "SysDistort"+closureTestSpecifier+"Summer");
	    elReweighted.ReplaceAll("Summer", "SysDistort"+closureTestSpecifier+"Summer");
	  }								\
	  // closure test pseudo data is created from ttbar simulation without madspin
	  muReweighted.ReplaceAll("MadSpin", "");
	  elReweighted.ReplaceAll("MadSpin", "");
	  TFile* mufile = new (TFile)(muReweighted);
	  TFile* elfile = new (TFile)(elReweighted);
	  // get plot
	  TString partonPlot="analyzeTop"+LV+"LevelKinematics"+PS+"/"+variable;
	  if(hadron){
	    // composited objects analyzer
	    if(variable.Contains("rhos")||variable.Contains("Njets")){
	      partonPlot.ReplaceAll("analyzeTop"     , "composited");
	      partonPlot.ReplaceAll("Njets"          , "Ngenjets"  );
	      partonPlot.ReplaceAll("rhos"           , "rhosGen"   );
	      partonPlot.ReplaceAll("LevelKinematics", "Gen"       );
	    }
	    // bjets analyzer
	    else if(variable.Contains("bq")||variable.Contains("lb")||variable.Contains("bbbar")){
	      partonPlot.ReplaceAll("PhaseSpace", "BjetsPhaseSpace");
	      partonPlot+="Gen";
	    }
	    // lepton analyzer
	    else if(variable.Contains("lep")){
	      partonPlot.ReplaceAll("PhaseSpace", "LeptonPhaseSpace");
	      partonPlot+="Gen";
	    }
	  }
	  histo_["reweighted"+variable+"Mu"][kSig]=(TH1F*)(mufile->Get(partonPlot)->Clone("rewmu"+variable));
	  histo_["reweighted"+variable+"El"][kSig]=(TH1F*)(elfile->Get(partonPlot)->Clone("rewel"+variable));
	  // add zprime for zprime closure test
	  if(zprime!=""){
	    // zprime files
	    TString muzprime=rewfold+"zprime/"+TopFilename(kSig, sysNo, "muon"    );
	    TString elzprime=rewfold+"zprime/"+TopFilename(kSig, sysNo, "electron");
	    // name for chosen mass
	    TString massextension="";
	    if(zprime=="1000") massextension="ZprimeM1000W100";
	    muzprime.ReplaceAll("Sig", massextension+"Sig");
	    elzprime.ReplaceAll("Sig", massextension+"Sig");
	    muzprime.ReplaceAll("MadSpin", "");
	    elzprime.ReplaceAll("MadSpin", "");
	    // get files
	    TFile* zprimemufile = new (TFile)(muzprime);
	    TFile* zprimeelfile = new (TFile)(elzprime);
	    // get plots
	    if(zprimemufile) histo_["zprime"+variable+"Mu"][kSig]=(TH1F*)(zprimemufile->Get(partonPlot)->Clone("zpmu"+variable));
	    if(zprimeelfile) histo_["zprime"+variable+"El"][kSig]=(TH1F*)(zprimeelfile->Get(partonPlot)->Clone("zpel"+variable));
	    // apply lumiweights
	    double zPrimeLumiWeight=zPrimeLumiWeightIni;
	    if     (zprime=="1000") zPrimeLumiWeight=(zPrimeLumiWeight*5*luminosity)/104043;
	    if(histo_["zprime"+variable+"Mu"].count(kSig)>0) histo_["zprime"+variable+"Mu"][kSig]->Scale(zPrimeLumiWeight);
	    if(histo_["zprime"+variable+"El"].count(kSig)>0) histo_["zprime"+variable+"El"][kSig]->Scale(zPrimeLumiWeight);
	    histo_["reweighted"+variable+"Mu"][kSig]->Scale(lumiweight(kSig, constLumiMuon, 0, "muon"    ,false));
	    histo_["reweighted"+variable+"El"][kSig]->Scale(lumiweight(kSig, constLumiElec, 0, "electron",false));
	    // add zprime to signal
	    if(histo_["zprime"+variable+"Mu"].count(kSig)>0) histo_["reweighted"+variable+"Mu"][kSig]->Add((TH1F*)histo_["zprime"+variable+"Mu"][kSig]->Clone());
	    if(histo_["zprime"+variable+"El"].count(kSig)>0) histo_["reweighted"+variable+"El"][kSig]->Add((TH1F*)histo_["zprime"+variable+"El"][kSig]->Clone());
	  }
	  // add channels
	  histo_[xSecNorm+"closure"][kSig]=(TH1F*)(histo_["reweighted"+variable+"Mu"][kSig]->Clone(xSecNorm+"closure"));
	  histo_[xSecNorm+"closure"][kSig]->Add(   histo_["reweighted"+variable+"El"][kSig]);
	  // apply standard rebinning
	  reBinTH1F(*histo_[xSecNorm+"closure"][kSig], binning_[variable], verbose-1);
	  // ensure that under- and overflow bins are empty
	  histo_[xSecNorm+"closure"][kSig]->SetBinContent(0                                              , 0.);
	  histo_[xSecNorm+"closure"][kSig]->SetBinContent(histo_[xSecNorm+"closure"][kSig]->GetNbinsX()+1, 0.);
	  // normalization to unit area
	  histo_[xSecNorm+"closure"][kSig]->Scale(1./(histo_[xSecNorm+"closure"][kSig]->Integral(0, histo_[xSecNorm+"closure"][kSig]->GetNbinsX()+1)));
	  // divide by binwidth
	  histo_[xSecNorm+"closure"][kSig] = divideByBinwidth(histo_[xSecNorm+"closure"][kSig], verbose-1);
	  // set style
	  histogramStyle(*histo_[xSecNorm+"closure"][kSig], kSig, false, 1.2, kMagenta);
	  histo_[xSecNorm+"closure"][kSig]->SetLineWidth(4);
	  histo_[xSecNorm+"closure"][kSig]->SetLineStyle(2);
	}
      }
    }
  }
  
  // ===============================================================
  //  Errors for uncertainty bands from ttbar Xsec and luminosity
  // ===============================================================
  
  if(verbose>0) std::cout << std::endl << " Start calculating error bands for 1D plots .... ";
  if(errorbands) makeUncertaintyBands(histo_, histoErrorBand_, plotList_, N1Dplots);
  if(verbose>0) std::cout << " .... Finished." << std::endl; 

  // =========================================
  //  Create one legend for all 1D histos
  // =========================================
  unsigned int Nlegends=0;  
  TLegend *leg  = new TLegend();  
  TLegend *leg0 = new TLegend(0.05, 0.15, 1.05, 0.9);
  Nlegends+=2;
  leg ->SetFillStyle(0);
  leg ->SetBorderSize(0);
  leg ->SetTextSize(0.03);
  leg ->SetTextAlign(12);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->SetHeader("Kin. Fit (after selection)");
  // fill in contributing sample
  // data is to be first entry
  TString lumilabel = Form("%3.2f fb^{-1}",luminosity/1000);
  TString dataLabel=sampleLabel(kData,decayChannel);
  if(closureTestSpecifier!=""){
    dataLabel.ReplaceAll("Data", "Pseudo Data");
    dataLabel.ReplaceAll("data", "pseudo data");
  }
  leg ->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], dataLabel,"P");
  leg0->AddEntry(histo_[plotList_[plotList_.size()-1]][kData], dataLabel+", "+lumilabel,"P");
  // now loop over MC samples
  for(unsigned int sample=kSig; sample<kData; ++sample){
    // check if sampe exists in at least one plot
    bool exit=false;
    // loop plots
    for(unsigned int plot=0; plot<plotList_.size()-1; ++plot){  // <plotList_.size()-1, because last entry is for data (see above)
      // if found: add entry to legend
      if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0)&&(!exit)){
	leg ->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	leg0->AddEntry(histo_[plotList_[plot]][sample], sampleLabel(sample,decayChannel), "F");
	exit=true;
      }
    }
  }
  // add entry for uncertainty to legends
  if (errorbands && histoErrorBand_.size() > 0 && plotList_.size() > 0){
      leg ->AddEntry(histoErrorBand_[plotList_[0]],"Uncertainty","F");
      leg0->AddEntry(histoErrorBand_[plotList_[0]],"Uncertainty","F");
  }
  // resize legends
  leg->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.20);
  leg->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.03*leg->GetNRows());
  leg->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
  leg->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength());

  leg0->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.20);
  leg0->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - 0.03*leg->GetNRows());
  leg0->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
  leg0->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength());

  // create in canvas legend for residuals
  TLegend *legPull = new TLegend(0.6, 0.5, 0.96, 0.9);
  legPull->SetFillStyle(0);
  legPull->SetBorderSize(0);
  legPull->SetTextSize(0.04);
  //legPull->SetEntrySeparation(0.45);
  legPull->AddEntry((TObject*)0, "#bf{rel. variation fit}", "");
  if(histo_.count("analyzeHypoKinFitLepton/leptonPullPtKinFitRec")>0&&histo_["analyzeHypoKinFitLepton/leptonPullPtKinFitRec"].count(kSig)>0){
    legPull->AddEntry(histo_["analyzeHypoKinFitLepton/leptonPullPtKinFitRec"     ][kSig], "(reco - fit) / reco", "L");
  }
  legPull->AddEntry((TObject*)0, "#bf{residuals}", "");
  if(histo_.count("analyzeHypoKinFitLepton/leptonPullPtRecPartonTruth")>0&&histo_["analyzeHypoKinFitLepton/leptonPullPtRecPartonTruth"].count(kSig)>0){
    legPull->AddEntry(histo_["analyzeHypoKinFitLepton/leptonPullPtRecPartonTruth"][kSig], "#splitline{before KinFit}{reco wrt. parton}"      , "L");
  }  
  if(histo_.count("analyzeHypoKinFitLepton/leptonPullPtKinFitPartonTruth")>0&&histo_["analyzeHypoKinFitLepton/leptonPullPtKinFitPartonTruth"].count(kSig)>0){
    legPull->AddEntry(histo_["analyzeHypoKinFitLepton/leptonPullPtKinFitPartonTruth"][kSig], "#splitline{after KinFit}{fitted wrt. parton}"       , "L");
  }

  // create a legend for tau closure tests
  TLegend *legTau = new TLegend(0.47, 0.81, 0.94, 0.87);
  legendStyle(*legTau, "", 0.03, 12);
  TLegend *legTau2 = new TLegend( 0.65, 0.41, 0.93, 0.81);
  TString dataLabel2= (closureTestSpecifier=="" ? "data" : "pseudo data");
  legendStyle(*legTau2, "#font[22]{unfolded "+dataLabel2+"}", 0.03, 12);
  if(closureTestSpecifier!="") legTau->AddEntry(histo_["xSecNorm/topPtclosure"][kSig], "gen truth of pseudo data", "L");
  unsigned int stdtau=positionInVector(tau_, stdTau_["topPt"]);
  legTau->AddEntry(histo_["xSecNorm/topPt"][kSig], "simulation used for unfolding", "L");
  for(int tau=0; tau<(int)tau_.size(); ++tau){
    if(tau<kindex){
      if(tau!=(int)stdtau) legTau2->AddEntry(tauxSecNorm_[tau], TString(Form("#tau=%2.1f",tau_[tau]))+TString(Form("#scale[0.85]{ (=%2.2f*#tau_{d})}", tau_[tau]/tau_[stdtau])), "LP");
      else legTau2->AddEntry(tauxSecNorm_[stdtau], Form("#tau_{d}=%2.1f (default)", tau_[stdtau]), "LP");
    }
    else{
      TString hbw= double(tauxSecNorm_[tau]->GetNbinsX())==2*tau_[tau] ? " (=0.5*N_{bins})" : "";
      legTau2->AddEntry(tauxSecNorm_[tau], TString(Form("k=%2.0f", tau_[tau]))+hbw, "LP");
    }
  }

  // ---
  //     print information about plots, canvas and legends
  // ---
  if(verbose>1){
    std::cout << std::endl << "# different plots requested: 1D + 2D = ";
    std::cout << N1Dplots + N2Dplots << " ( " << N1Dplots << " + " << N2Dplots << " )" << std::endl;
    std::cout << "1D plots for different samples will be drawn ";
    std::cout << "into same canvas as stacked plot" << std::endl;
    std::cout << "every 2D plot will be drawn into a seperate canvas" << std::endl;
    std::cout << "#Canvas for plots: " << Nplots << std::endl;
    std::cout << "#legends: " << Nlegends << std::endl;
    std::cout << "every legend will be drawn into a seperate canvas" << std::endl;
    std::cout << "#Canvas = #plots + #legends: " << Nplots+Nlegends << std::endl<< std::endl;
  }

  // ---
  //    create canvas
  // ---
  std::vector<TCanvas*> plotCanvas_;
  // entries added when printing is done

  // ---
  //    change 1D plots into stack plots
  // ---
  // loop plots 
  if(verbose>2) std::cout << std::endl;
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    // show kinfit permutation plot as relative plot for signal
    // and delete this plot for all other samples 
    if(plotList_[plot].Contains("qAssignment")||plotList_[plot].Contains("decayChannel")){
      TString permutationLabel=plotList_[plot];
      if(verbose>2) std::cout << "adapting plot " << permutationLabel << std::endl;
      for(unsigned int sample=kSig; sample<=kData; ++sample){
	if(histo_[permutationLabel].count(sample)>0){
	  if((sample!=kSig&&plotList_[plot].Contains("qAssignment"))||(sample!=kBkg&&plotList_[plot].Contains("decayChannel"))){
	    histo_[permutationLabel].erase(sample);
	    if(verbose>2) std::cout << "erase plot for sample " << sampleLabel(sample, decayChannel) << std::endl;
	  }
}      
      }
      if((histo_[permutationLabel].count(kSig)>0&&plotList_[plot].Contains("qAssignment"))||(histo_[permutationLabel].count(kBkg)>0&&plotList_[plot].Contains("decayChannel"))){ 
	if(verbose>2) std::cout << "scale to area" << std::endl;
	unsigned int sample = (plotList_[plot].Contains("qAssignment") ? kSig : kBkg);
	histo_[permutationLabel][sample]->Scale(1./histo_[permutationLabel][sample]->Integral(0, histo_[permutationLabel][sample]->GetNbinsX()));
	axisLabel_[plot].ReplaceAll("events","relative fraction of events");
	if(plotList_[plot].Contains("qAssignment")){
	  histogramStyle( *histo_[permutationLabel][kSig], kSig, true);
	  // print permutations
	  if(systematicVariation==testMe||systematicVariation==sysNo){
	    std::cout << std::endl << "permutations wrt jet parton matching: " << std::endl;
	    std::cout << "-------------------------------------" << std::endl;
	    std::cout << "ok     : " << histo_[permutationLabel][kSig]->GetBinContent(1) << std::endl;
	    std::cout << "(correct jet-quark assignment)" << std::endl;
	    std::cout << "b<->b  : " << histo_[permutationLabel][kSig]->GetBinContent(2) << std::endl;
	    std::cout << "(hadronic and leptonic b-jet interchanged)" << std::endl;
	    double others=0.;
	    for(int bin=3; bin<=7; ++bin){
	      others+=histo_[permutationLabel][kSig]->GetBinContent(bin);
	    }
	    std::cout << "q<->q  : " << others << std::endl;
	    std::cout << "(2 or more other quarks interchanged)" << std::endl;
	    std::cout << "jmis   : " << histo_[permutationLabel][kSig]->GetBinContent(8) << std::endl;
	    std::cout << "(jet(s) missing in leading 5) " << std::endl;
	    std::cout << "wrongj : " << histo_[permutationLabel][kSig]->GetBinContent(9) << std::endl;
	    std::cout << "(wrong 4 out of the 5 leading chosen)" << std::endl;
	    std::cout << "nomatch: " << histo_[permutationLabel][kSig]->GetBinContent(10) << std::endl;
	    std::cout << "(no unambigious matching due to e.g. jet splitting(ISR/FSR), jet merging,";
	    std::cout << " jet(s) out of acceptance or jets close by)" << std::endl << std::endl;
	  }
	}
      }
    }
    // pull plots as line, signal only, no stack, relative 
    else if(plotList_[plot].Contains("analyzeHypoKinFitLepton")||plotList_[plot].Contains("analyzeHypoKinFitMET")||plotList_[plot].Contains("analyzeHypoKinFitJets")){
      TString pullLabel=plotList_[plot];
      for(unsigned int sample=kBkg; sample<=kData; ++sample){
	if(histo_[pullLabel].count(sample)>0)  histo_[pullLabel].erase(sample);
      }
      if(histo_[pullLabel].count(kSig)>0){ 
	histo_[pullLabel][kSig]->Scale(1./histo_[pullLabel][kSig]->Integral(0, histo_[pullLabel][kSig]->GetNbinsX()));
	axisLabel_[plot].ReplaceAll("Frequency","event fraction");
	histogramStyle( *histo_[pullLabel][kSig], kSig, false);
	histo_[pullLabel][kSig]->SetLineWidth(3);
	if(pullLabel.Contains("RecPartonTruth")) histo_[pullLabel][kSig]->SetLineColor(kBlue );
	if(pullLabel.Contains("KinFitRec"     )) histo_[pullLabel][kSig]->SetLineColor(kBlack);
      }
    }
    // otherwise: all 1D plots will become stacked plots
    else createStackPlot(plotList_, histo_, plot, N1Dplots, verbose-1, decayChannel);
  }
  if(verbose>2) std::cout << std::endl;
  
  // =====================
  //  Do the printing
  // =====================
  // a) for plots
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      bool first=true;
      // loop samples
      unsigned int sample=0;
      for(sample=kSig; sample<=kData; ++sample){
	  // a1) for 1D event yields, efficiency and cross section plots (existing)
	  if((plot<N1Dplots)||(plot>=N1Dplots+N2Dplots)){
	      // check if plot is existing
	      // draw BBB xSec control plots in same canvas with SVD result
	      if(((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample)>0))&&!plotList_[plot].Contains("BBB")){
		  histo_[plotList_[plot]][sample]->SetName(((TString)histo_[plotList_[plot]][sample]->GetName()).ReplaceAll("Gen",""));
		  // draw all pull distributions in same canvas if RecPartonTruth pull is called
		  if((!plotList_[plot].Contains("Pull")||plotList_[plot].Contains("RecPartonTruth"))){
		      if(verbose>1){
			  unsigned int activeCanv =(plotCanvas_.size()==0 ? 1 : plotCanvas_.size()-1);
			  std::cout << " Plotting "    << plotList_[plot];
			  std::cout << " from sample " << sampleLabel(sample,decayChannel);
			  std::cout << " to canvas "   << activeCanv << std::endl;
		      }
		      // for efficiency plots: draw grid
		      if(getStringEntry(plotList_[plot], 1)=="efficiency") plotCanvas_[plotCanvas_.size()-1]->SetGrid(1,1);
		      // adjust x range
		      setXAxisRange(histo_[plotList_[plot]][sample], getStringEntry(plotList_[plot], 2));
		      if(!(plotList_[plot].Contains("xSec"))&&(getStringEntry(plotList_[plot], 2)=="topY"   ||
							       getStringEntry(plotList_[plot], 2)=="topYHad"|| 
							       getStringEntry(plotList_[plot], 2)=="topYLep")){histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(-3,3);}
		      if(getStringEntry(plotList_[plot], 2)=="PartonJetDRall") histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(0,4);
		      if(getStringEntry(plotList_[plot], 2)=="decayChannel"  ) histo_[plotList_[plot]][sample]->GetXaxis()->SetRange(5,histo_[plotList_[plot]][sample]->GetNbinsX());
		      if(first){ 
			  // create canvas and set titel corresponding to plotname in .root file
			  addCanvas(plotCanvas_);
			  plotCanvas_[plotCanvas_.size()-1]->cd(0);
			  TString canvtitle=getStringEntry(plotList_[plot], 2)+getStringEntry(plotList_[plot], 1);
			  canvtitle.ReplaceAll(addSel,"");
			  //std::cout << canvtitle << std::endl;
			  plotCanvas_[plotCanvas_.size()-1]->SetTitle(canvtitle);
			  // min / max
			  double max = 1.3*histo_[plotList_[plot]][sample]->GetMaximum();
			  // if data file exists
			  if(histo_[plotList_[plot]].count(kData)>0){
			    // and has larger maximum
			    if(max < 1.3*histo_[plotList_[plot]][kData]->GetMaximum()){
			      // take this maximum
			      max = 1.3*histo_[plotList_[plot]][kData]->GetMaximum();
			      if(verbose>1) std:: cout << "take max from data! " << std::endl;
			    }
			  }
			  double min = 0;
			  // log plots
			  if(getStringEntry(axisLabel_[plot],3)=="1"){
			      plotCanvas_[plotCanvas_.size()-1]->SetLogy(1);
			      min=1;
			      max=exp(1.3*(std::log(max)-std::log(min))+std::log(min));
			      if(plotList_[plot]=="analyzeTopRecoKinematicsKinFit"+sysInputFolderExtension+"/prob") min=0.1; 
			      if(plotList_[plot].Contains("ttbarMass")) min=0.1;
			      if(plotList_[plot].Contains("ttbarMass")&&plotList_[plot].Contains("xSec")){
				  min=0.0001;
				  max=1.2*exp(1.3*(std::log(max)-std::log(min))+std::log(min));
				  if(plotList_[plot].Contains("Norm")){
				      min=0.000001;
				      max=0.06;
				  }
			      }
			  }
			  // get nicer int values if maximum is large enough
			  if(max>3) max = (double)roundToInt(max);
			  // style for pull plots
			  if(plotList_[plot].Contains("RecPartonTruth")){
			      min=0;
			      //max=1;
			      TString title=(TString)plotCanvas_[plotCanvas_.size()-1]->GetTitle();
			      title.ReplaceAll("RecPartonTruth","");
			      plotCanvas_[plotCanvas_.size()-1]->SetTitle(title);
			      histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
			      plotCanvas_[plotCanvas_.size()-1]->SetGrid(1,1);
			  }
			  if(plotList_[plot].Contains("qAssignment")||plotList_[plot].Contains("decayChannel")) plotCanvas_[plotCanvas_.size()-1]->SetGrid(1,1);
			  // axis style
			  axesStyle(*histo_[plotList_[plot]][sample], getStringEntry(axisLabel_[plot],1), getStringEntry(axisLabel_[plot],2), min, max); 
			  histo_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
			  if(max>1&&max<100) histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
			  else histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(false);
			  if(plotList_[plot].Contains("qAssignment")||plotList_[plot].Contains("decayChannel")) histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
			  if(getStringEntry(plotList_[plot], 1).Contains("xSec")) histo_[plotList_[plot]][sample]->GetYaxis()->SetTitleOffset(1.6);
			  // draw efficiency plots as line
			  if(getStringEntry(plotList_[plot], 1)=="efficiency") histo_[plotList_[plot]][sample]->Draw("p e");
			  // draw pull plots as line into same canvas with extra legend
			  else if(plotList_[plot].Contains("RecPartonTruth")){
			      histo_[plotList_[plot]][sample]->GetXaxis()->SetRangeUser(-1.5,1.5);
			      histo_[plotList_[plot]][sample]->GetXaxis()->SetTitle(((TString)(histo_[plotList_[plot]][sample]->GetXaxis()->GetTitle())).ReplaceAll(":","/"));
			      histo_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);	
			      histo_[plotList_[plot]][sample]->Draw("hist");
			      TString afterKinFit=plotList_[plot];
			      afterKinFit.ReplaceAll("RecPartonTruth","KinFitPartonTruth");
			      TString kinFitShift=plotList_[plot];
			      kinFitShift.ReplaceAll("RecPartonTruth","KinFitRec"        );
			      if(histo_.count(afterKinFit)>0){ 
				  histo_[plotList_[plot]][sample]->SetMaximum(1.3*histo_[afterKinFit][sample]->GetMaximum());
				  histo_[afterKinFit][sample]->Draw("hist same");
			      }
			      if(histo_.count(kinFitShift)>0) histo_[kinFitShift][sample]->Draw("hist same");
			      legPull->Draw("same");
			  }
			  else if(plotList_[plot].Contains("xSecNorm")){
			    TH1F* tempAxis=(TH1F*)histo_[plotList_[plot]][sample]->Clone();
			    tempAxis->GetXaxis()->SetLabelSize(0);
			    tempAxis->GetXaxis()->SetTitle("");	
			    tempAxis->Draw("AXIS");
			    histo_[plotList_[plot]][sample]->Draw("hist same");
			  }
			  // others as histo (stack)
			  else histo_[plotList_[plot]][sample]->Draw("hist");
			  histo_[plotList_[plot]][42] = (TH1F*)(histo_[plotList_[plot]][sample]->Clone());
			  // draw BBB comparison results in same canvas
			  if(compare&&plotList_[plot].Contains("xSec")){
			      TString BBBlabel=plotList_[plot];
			      BBBlabel.ReplaceAll("xSec/","xSec/BBB");
			      BBBlabel.ReplaceAll("xSecNorm/","xSecNorm/BBB");
			      if(plotExists(histo_, BBBlabel, sample)) histo_[BBBlabel][sample]->Draw("hist same");
			  }
		      }
		      // draw other plots into same canvas
		      else{ 
			  // draw data as points
			  if(sample==kData){ 
			    if(!plotList_[plot].Contains("xSecNorm")) histo_[plotList_[plot]][sample]->Draw("p e X0 same");
			      // draw BBB comparison results in same canvas
			      if(compare&&plotList_[plot].Contains("xSec")){
				  TString BBBlabel=plotList_[plot];
				  BBBlabel.ReplaceAll("xSec/","xSec/BBB");
				  BBBlabel.ReplaceAll("xSecNorm/","xSecNorm/BBB");
				  if(plotExists(histo_, BBBlabel, sample)) histo_[BBBlabel][sample]->Draw("p e X0 same");
			      }
			      if(plotList_[plot].Contains("xSecNorm")){
				// draw theory truth curve
				if(closureTestSpecifier!="") histo_["xSecNorm/topPtclosure"][kSig]->Draw("hist same");
				// draw all tau parameter results into the same plot
				for(int xSecnum=0; xSecnum<(int)tauxSecNorm_.size(); ++xSecnum){
				  //tauxSecNorm_[xSecnum]->Draw("hist same");
				  tauxSecNormM_[xSecnum]->Draw("p same");
				}
				// draw corresponding legend
				legTau ->Draw("same");
				legTau2->Draw("same");
			      }
			  }
			  else{
			      // draw efficiencies as points
			      if(getStringEntry(plotList_[plot], 1)=="efficiency") histo_[plotList_[plot]][sample]->Draw("p e same");
			      // draw others as histo (stack)
			      else histo_[plotList_[plot]][sample]->Draw("hist same");
			  }
		      }
		      first=false;
		  }
		  first=false;	
		  // draw uncertainty bands, add legend and labels and re-draw axis
		  if(sample==kData){
		    // configure style of and draw uncertainty bands
		    if (!plotList_[plot].Contains("xSec")){
		      if(errorbands){
			histoErrorBand_[plotList_[plot]]->SetMarkerStyle(0);
			histoErrorBand_[plotList_[plot]]->SetFillColor(1);
			histoErrorBand_[plotList_[plot]]->SetFillStyle(3004);
			gStyle->SetErrorX(0.5);  
			histoErrorBand_[plotList_[plot]]->Draw("E2 SAME");
		      }
		    }
		    // draw legend for recoYield plots
		    TString tempTitle = plotCanvas_[plotCanvas_.size()-1]->GetTitle();
		    if(tempTitle.Contains("analyzeTopRecoKinematicsKinFit")){
		      if(!plotList_[plot].Contains("qAssignment")&&!plotList_[plot].Contains("decayChannel")) leg->Draw("SAME");
		    }	 
		    // labels
		    TString plotType=getStringEntry(plotList_[plot], 1);
		    if(plotType.Contains("xSec")||plotType.Contains("Reco")){
		      if (decayChannel=="muon")         DrawDecayChLabel("#mu + Jets");
		      else if(decayChannel=="electron") DrawDecayChLabel("e + Jets");
		      else if(decayChannel=="combined") DrawDecayChLabel("e/#mu + Jets");
		      DrawCMSLabels(prelim, luminosity, 0.04, (closureTestSpecifier!="" ? true : false), false, false);
		    }
		    if (!plotList_[plot].Contains("xSec")){
		      // redraw axis
		      histo_[plotList_[plot]][42]->Draw("axis same");
		    }
		    else{
		      // label indicating the type of closure test performed
		      if(closureTestSpecifier!=""){
			double positionX=525;
			double positionY=0.;
			TLatex *closurelabel = new TLatex(positionX,positionY,pseudoDatalabel(closureTestSpecifier));
			closurelabel->SetTextAlign(11);
			closurelabel->SetTextAngle(90);
			closurelabel->SetTextSize(0.04);
			closurelabel->Draw("same");
		      }
		      // ratio for different taus
		      std::vector<double> zeroerr_;
		      for(int bin=1; bin<=tauxSecNorm_[stdtau]->GetNbinsX(); ++bin) zeroerr_.push_back(0);
		      if(closureTestSpecifier=="") histo_["xSecNorm/topPtclosure"][kSig]=histo_[plotList_[plot]][42];
		      double ratioMin=closureTestSpecifier=="" ? 0.5 : 0.7;
		      double ratioMax=closureTestSpecifier=="" ? 1.25 : 1.35;
		      drawRatio(histo_[plotList_[plot]][42], histo_[plotList_[plot]][42], ratioMin, ratioMax, myStyle, verbose, zeroerr_, (closureTestSpecifier=="" ? "prediction" : "unfolded"), (closureTestSpecifier=="" ? "unfolded" : "gen truth"), "AXIS", histo_["xSecNorm/topPtclosure"][kSig]->GetLineColor() , false, 0.1);
		      drawRatio(histo_["xSecNorm/topPtclosure"][kSig], histo_["xSecNorm/topPtclosure"][kSig], 0.8, 1.35, myStyle, verbose, zeroerr_, "unfolded", "gen truth", "hist same", histo_["xSecNorm/topPtclosure"][kSig]->GetLineColor() , false, 0.1);
		      for(int xSecnum=0; xSecnum<(int)tauxSecNorm_.size(); ++xSecnum){
			TGraphAsymmErrors* tempRatio=(TGraphAsymmErrors*)tauxSecNormM_[xSecnum]->Clone();
			for(int bin=1; bin<=tauxSecNorm_[xSecnum]->GetNbinsX(); ++bin){
			  tempRatio->SetPoint(bin-1, tempRatio->GetX()[bin-1], ((tempRatio->GetY()[bin-1])/(histo_["xSecNorm/topPtclosure"][kSig]->GetBinContent(bin))));
			}			
			tempRatio->Draw("p same");
			//drawRatio(tauxSecNorm_[xSecnum], histo_["xSecNorm/topPtclosure"][kSig], 0.8, 1.2, myStyle, verbose, zeroerr_, "unfolded", "gen truth", "hist same", tauxSecNorm_[xSecnum]->GetLineColor() , false, 0.1);
		      }
		      for(int bin=2; bin<=tauxSecNorm_[stdtau]->GetNbinsX(); ++bin){
			double xBorder=tauxSecNorm_[stdtau]->GetBinLowEdge(bin);
			drawLine(xBorder, ratioMin, xBorder, ratioMax, kBlack, 1, 3);
		      }
		    }
		  }
	      } 
	  }
	  // a3) for 2D plots
	  if((plot>=N1Dplots)&&(histo2_.count(plotList_[plot])>0)&&(histo2_[plotList_[plot]].count(sample)>0)&&	  
	     // draw gen-reco correlation only for signal 
	     (!plotList_[plot].Contains("_"   )||sample==kSig)){
	      // new Canvas for every plot
	      addCanvas(plotCanvas_);
	      plotCanvas_[plotCanvas_.size()-1]->cd(0);
	      plotCanvas_[plotCanvas_.size()-1]->SetRightMargin(myStyle.GetPadRightMargin());
	      TString canvname=getStringEntry(plotList_[plot],2)+getStringEntry(plotList_[plot],1)+getTStringFromInt(sample);
	      canvname.ReplaceAll(addSel,"");
	      // special naming of hadron level quanatities
	      if(canvname.Contains("Ngenjets")) canvname.ReplaceAll("Ngenjets","Njets");
	      plotCanvas_[plotCanvas_.size()-1]->SetTitle(canvname);
	      if(verbose>1){
		  std::cout << " Plotting "    << plotList_[plot];
		  std::cout << " from sample " << sampleLabel(sample,decayChannel);
		  std::cout << " to canvas "   << plotCanvas_.size()-1  << " ( ";
		  std::cout << plotCanvas_[plotCanvas_.size()-1]->GetTitle() << " )"  << std::endl;
	      }
	      // draw histo
	      plotCanvas_[plotCanvas_.size()-1]->SetRightMargin(0.15);
	      histo2_[plotList_[plot]][sample]->GetXaxis()->SetNoExponent(true);
	      histo2_[plotList_[plot]][sample]->GetYaxis()->SetNoExponent(true);
	      histo2_[plotList_[plot]][sample]->Draw("colz");
	      // print correlation factor
	      double d = histo2_[plotList_[plot]][sample]->GetCorrelationFactor();
	      char correlation[20];
	      sprintf(correlation, "%f", d);
	      TString corr = (TString)correlation;
	      DrawLabel("correlation: "+corr, 0.35, 0.92, 0.75, 0.99, 0.7);
	  }
      } // end of loop over samples
  } // end of loop over plots
  
  // b) for legends
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle("legend");
  leg0->Draw("");

  // ==============
  //  Saving
  // ==============
  if(save){
    TString iterext=iterative ? "Iter" : "";
    // pdf and eps only for standard analysis without variation
    if(systematicVariation==sysNo){
      if(verbose>0) std::cout << "will save all plots as png/eps/pdf" << std::endl;
      if(verbose<=1) gErrorIgnoreLevel=kWarning;
      // a) as pdf
      saveCanvas(plotCanvas_, outputFolder, pdfName, true, false);
      // b) as eps and png
      for(unsigned int idx=0; idx<plotCanvas_.size(); idx++){
	TString saveToFolder=outputFolder;
	TString title=(plotCanvas_[idx])->GetTitle();	
	saveToFolder+="regularizationTest/";
	if(!title.Contains("canv")){
	  // add additional label that indicates PS for all relevant plots
	  TString universalplotLabel="";
	  if(title.Contains("0")||(!title.Contains("analyzeTopRecoKinematicsKinFit")&&!title.Contains("legend"))){
	    if(extrapolate) universalplotLabel="FullPS";
	    else universalplotLabel=LV+"LvPS";
	  }
	  if(title.Contains("xSec")){
	    if(!SVDunfold) universalplotLabel+="BBB";
	    if(compare   ) universalplotLabel+="BBBcomparison";
	  }
	  // do the saving
	  TString canvtit=plotCanvas_[idx]->GetTitle();
	  if(!(canvtit.Contains("legend")||(canvtit.Contains("xSec")&&!canvtit.Contains("Norm")))){
	    plotCanvas_[idx]->Print(saveToFolder+canvtit+(TString)(universalplotLabel+closureLabel+test+iterext+".eps")); 
	    plotCanvas_[idx]->Print(saveToFolder+canvtit+(TString)(universalplotLabel+closureLabel+test+iterext+".png"));
	  }
	}
      }
      if(verbose>0) std::cout << "done" << std::endl;
    }
    // c) root file
    std::vector< std::pair< TCanvas*,TString > > saveCanvas_;
    if(verbose>0) std::cout << "will save all plots in outputfile named " << outputFileName << std::endl;
    for(unsigned int idx=0; idx<plotCanvas_.size(); ++idx){
      // get correct folder
      TString title=(plotCanvas_[idx])->GetTitle();
      TString outputfolder="";
      TString possibleFolderNames[5]={"efficiency", "analyzeTopRecoKinematicsKinFit", "xSec", "analyzeTopPartonLevelKinematics", "analyzeHypoKinFit"};
      for(unsigned int name=0; name<sizeof(possibleFolderNames)/sizeof(TString); ++name){
	// add another subfoder indicating the systematic variation
	if(title.Contains(possibleFolderNames[name])){ 
	  outputfolder=possibleFolderNames[name]+"/"+sysLabel(systematicVariation);
	  plotCanvas_[idx]->SetTitle(title.ReplaceAll(possibleFolderNames[name]+iterext,""));
	}
      }
      // save only canvas from selected subfolders or legends
      if(outputfolder!=""||title.Contains("legend")){ 
	// collect information in saveCanvas_
	saveCanvas_.push_back(std::make_pair((TCanvas*)plotCanvas_[idx]->Clone(),outputfolder));
      }
    }
    // save all collected plots
    saveToRootFileAll(outputFileName, saveCanvas_, true, verbose-1);
    if(verbose>0) std::cout << "done" << std::endl; 
  }
  
  // delete pointer
  delete leg0;
  closeStdTopAnalysisFiles(files_);
}
