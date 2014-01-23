#include "basicFunctions.h"

void bothDecayChannelsCombination(double luminosity=19712, bool save=true, unsigned int verbose=0,
				  TString inputFolderName="RecentAnalysisRun8TeV_doubleKinFit",
				  bool pTPlotsLog=false, bool extrapolate=true, bool hadron=false, bool addCrossCheckVariables=false, 
				  bool combinedEventYields=true, TString closureTestSpecifier="" , bool smoothcurves=false){

  // closureTestSpecifier = \"NoDistort\", \"topPtUp\", \"topPtDown\", \"ttbarMassUp\", \"ttbarMassDown\", \"data\" or \"1000\"
  // run automatically in batch mode
  gROOT->SetBatch();
    
  // ============================
  //  Set Root Style
  // ============================
  
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.cd();
  TGaxis::SetMaxDigits(2);
  gROOT->SetStyle("HHStyle");
  
  // ============================
  //  Setup
  // ============================
  
  // decay channels
  enum channel {kMuon, kElectron, kCombined};

  // Define variables for combination (centrally defined in basicFunctions.h)
  std::vector<TString> xSecVariables_;
  // a) top and ttbar quantities
  if(!hadron){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFit, xSecVariablesKinFit + sizeof(xSecVariablesKinFit)/sizeof(TString));
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFitNorm, xSecVariablesKinFitNorm + sizeof(xSecVariablesKinFitNorm)/sizeof(TString));
  }
  // b) lepton and b-jet quantities
  if(hadron||!extrapolate){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalState    , xSecVariablesFinalState     + sizeof(xSecVariablesFinalState    )/sizeof(TString));
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalStateNorm, xSecVariablesFinalStateNorm + sizeof(xSecVariablesFinalStateNorm)/sizeof(TString));
  }
  // c) cross check variables (only available for parton level cross-sections)
  if (addCrossCheckVariables && !hadron){
    xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVar    , xSecVariablesCCVar     + sizeof(xSecVariablesCCVar    )/sizeof(TString));
    xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVarNorm, xSecVariablesCCVarNorm + sizeof(xSecVariablesCCVarNorm)/sizeof(TString));
  }
  xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesIncl,      xSecVariablesIncl      + sizeof(xSecVariablesIncl)/sizeof(TString)     );


  // ---
  //    create list of systematics to be treated as uncorrelated between electron and muon channel
  // ---
  // available systematics:
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

  std::vector<unsigned int> uncorrSys_;
  // trigger 
  //uncorrSys_.push_back(sysTriggerEffSFJetNormUp	);
  //uncorrSys_.push_back(sysTriggerEffSFJetNormDown	);
  //uncorrSys_.push_back(sysTriggerEffSFJetShapeUp	);
  //uncorrSys_.push_back(sysTriggerEffSFJetShapeDown);  
  // lepton eff
  uncorrSys_.push_back(sysLepEffSFNormUp       );        
  uncorrSys_.push_back(sysLepEffSFNormDown     );        
  uncorrSys_.push_back(sysLepEffSFShapeUpEta   );    
  uncorrSys_.push_back(sysLepEffSFShapeDownEta );
  uncorrSys_.push_back(sysLepEffSFShapeUpPt    );     
  uncorrSys_.push_back(sysLepEffSFShapeDownPt  );
  //uncorrSys_.push_back(sysLepEffSFShapeUpPt40  );   
  //uncorrSys_.push_back(sysLepEffSFShapeDownPt40);
  //uncorrSys_.push_back(sysLepEffSFNormUpStat   );    
  //uncorrSys_.push_back(sysLepEffSFNormUpStat   );
 
  // Label for datasample
  TString dataSample="2012";
  
  // Choose phase space
    // a) for full PS use extrapolatetrue;
  TString PS = extrapolate ? "" : "PhaseSpace";
  // b) for restricted phase space:
  // b1) parton PS: hadron=false
  // b2) hadron PS: hadron=true
  TString LV="Parton";
  if(!extrapolate&&hadron) LV="Hadron";
  if(verbose>1){
    if(extrapolate) std::cout << "full Phase Space will be used!" << std::endl; 
    else std::cout << LV << " level Phase Space will be used!" << std::endl; 
  }
  TString universalplotLabel = extrapolate ? "FullPS" : LV+"LvPS";
  
  // Choose additional theory curves to plot
  // Version 1: data, MadGraph, MC@NLO, POWHEG NNLO (Kidonakis)
  bool DrawSmoothMadgraph    = false;
  bool DrawMCAtNLOPlot       = true;
  bool DrawPOWHEGPYTHIAPlot  = true;
  bool DrawPOWHEGHERWIGPlot  = true;
  bool DrawNNLOPlot          = true;
  bool DrawMCFMPlot          = false;
  // smooth instead of binned theory curves
  if(smoothcurves){
    DrawSmoothMadgraph   = true;
    DrawPOWHEGPYTHIAPlot = true;
    DrawPOWHEGHERWIGPlot = true;
    DrawMCAtNLOPlot      = true;
    // FIXME: no smooth mc@nlo curve for hadron PS at the moment
    //if(hadron) DrawMCAtNLOPlot = false;
    // FIXME: -> use the one running the sample itself
  }	
  // use large 40M ttbar sample?
  bool largeMGfile=false;
  // add kidonakis plots for full PS
  //if(extrapolate&&smoothcurves) DrawNNLOPlot=true;

  //if(extrapolate==false&&hadron==false){
  //  DrawSmoothMadgraph = true;
  //  DrawMCAtNLOPlot = true;
  //  DrawPOWHEGPYTHIAPlot = true;
  //}
  // use only mcatnlo plots from external studies 
  // useOnlyExternalMCatNLOfile=false: means error bands from external studies
  //                            but central values from std file
  bool useOnlyExternalMCatNLOfile=false;
  // want to draw prediction with SC in addition?
  bool SC=false; 
  // GOSSIE quick fix: cut of m(ttbar) below 345 GeV
  bool cutTtbarMass=true;
  const double constMassRangeLow  =  345.0;
  const double constMassRangeHigh = 1600.0;
  
  // up to which systematic variation should be looped? Default: all
  unsigned int sysEnd = ENDOFSYSENUM;
  // for closure tests:
  if(closureTestSpecifier!=""){
    // no systematics
    sysEnd = 1;
    // no theory curve besides MadGraph+Pythia
    DrawMCAtNLOPlot       = false;
    DrawPOWHEGPYTHIAPlot  = false;
    DrawPOWHEGHERWIGPlot  = false;
    DrawNNLOPlot          = false;
  }
  
  TString dataLabel    = "Data";
  TString closureLabel = "";
  TString dataLabelrew = "#splitline{Reco.&Unfolded}{Reweighted t#bar{t}}";
  // Closure test with reweighted parton level shape or additional z prime signal
  // will plot additionally the modified diff. norm. xSec on parton level
  bool reweightClosure=false;
  TString zprime="";
  double  zPrimeLumiWeightIni=1.;
  TString zPrimeLumiWeightStr="";
  if(closureTestSpecifier!=""){
    if      (closureTestSpecifier.Contains("NoDistort")) {dataLabel="#splitline{Reco.&Unfolded}{Pseudo Data}"; reweightClosure=true; closureLabel = "PseudoData"+closureTestSpecifier;             } 
    else if (closureTestSpecifier.Contains("topPt"   )||
	     closureTestSpecifier.Contains("ttbarMass")) {dataLabel=dataLabelrew; reweightClosure=true; closureLabel = "PseudoDataReweight"+closureTestSpecifier;     }
    else if (closureTestSpecifier.Contains("data"     )) {dataLabel=dataLabelrew; reweightClosure=true; closureLabel = "PseudoDataReweighttopPt"+closureTestSpecifier;}
    else if (closureTestSpecifier.Contains("1000"     )) {
      closureLabel = "PseudoDataZprime"+closureTestSpecifier+"GeV"; 
      // identify exact closure test
      zprime="1000";
      if      (closureTestSpecifier.Contains("x0p03")) {zPrimeLumiWeightIni=0.03; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.03)";}
      else if (closureTestSpecifier.Contains("x0p1" )) {zPrimeLumiWeightIni=0.1 ; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.1)" ;}
      else if (closureTestSpecifier.Contains("x0p25")) {zPrimeLumiWeightIni=0.25; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.25)";}
      else if (closureTestSpecifier.Contains("x0p5" )) {zPrimeLumiWeightIni=0.5 ; zPrimeLumiWeightStr=" (#sigma_{Z'}x0.5)" ;}
      else if (closureTestSpecifier.Contains("x2"   )) {zPrimeLumiWeightIni=2.  ; zPrimeLumiWeightStr=" (#sigma_{Z'}x2)"   ;}
      else if (closureTestSpecifier.Contains("x4"   )) {zPrimeLumiWeightIni=4.  ; zPrimeLumiWeightStr=" (#sigma_{Z'}x4)"   ;}
      // adjust label in plot: data -> pseudo data
      dataLabel = "#splitline{Reco. and Unf.}{t#bar{t} + "+zprime+" GeV Z'"+zPrimeLumiWeightStr+"}";
    }
    else{
      std::cout << "ERROR: unknown closureTestSpecifier=" << closureTestSpecifier << std::endl;
      exit(0);
    }
  }
  //  Define muon and electron input rootfiles
  std::map<unsigned int, TFile*> files_;
  if(!combinedEventYields){
    files_[kMuon    ] = new (TFile)("diffXSecTopSemiMu"  +closureLabel+LV+PS+".root");
    files_[kElectron] = new (TFile)("diffXSecTopSemiElec"+closureLabel+LV+PS+".root");
  }
  else files_[kCombined] = new (TFile)("diffXSecTopSemiLep"+closureLabel+LV+PS+".root");
  // output for chosen input
  if(verbose>0){
    if(combinedEventYields) std::cout << "using COMBINED EVENT YIELD as input!!!" << std::endl;
    else std::cout << "using MU and EL SEPARATE files as input!!!" << std::endl;
  }
  
  // define folders where XSec plots are stored
  TString xSecFolder = "xSec";
  
  // create plot container for combined e+#mu plots
  std::map< TString, std::map <unsigned int, TH1F*> >    histo_;
  std::map< TString, std::map <unsigned int, TCanvas*> > canvas_;
  
  // =============================================================================
  //  Combination for all kinematic variables x systematic variations available
  // =============================================================================
  
  // loop systematic variations
  for(unsigned int sys=sysNo; sys<sysEnd; ++sys){
    //if(sys>sysNo) break; // FIXME: shortcut- process only results without systematic variation 
    TString subfolder=sysLabel(sys);
    // loop variables
    int countNorm=0;
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      // get canvas
      TCanvas* canvasMu   = combinedEventYields ? 0 : (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]));
      TCanvas* canvasEl   = combinedEventYields ? 0 : (TCanvas*)(files_[kElectron]->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]));
      TCanvas* canvasTheo = combinedEventYields ? 0 : (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i]));
      // if already combined at yield level:
      TCanvas* canvasComb     = combinedEventYields ? (TCanvas*)(files_[kCombined]->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]))       : 0;
      TCanvas* canvasTheoComb = combinedEventYields ? (TCanvas*)(files_[kCombined]->Get(xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i])) : 0;
      bool foundCanv=true;
      if(!(canvasMu&&canvasEl&&canvasTheo)&&!combinedEventYields){
	std::cout << std::endl << " WARNING in bothDecayChannelsCombination.C! " << std::endl;
	if(!canvasMu||!canvasEl) std::cout << " canvas " << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " not found in ";
	if(!canvasMu) std::cout << "diffXSecTopSemiMu"+closureLabel+LV+PS+".root (muon file)"     << std::endl;
	if(!canvasEl) std::cout << "diffXSecTopSemiElec"+closureLabel+LV+PS+".root (electron file)" << std::endl;
	if(!canvasTheo){ 
	  std::cout << " theory canvas " << xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i] << " not found in diffXSecTopSemiMu"+closureLabel+LV+PS+".root" << std::endl;
	}
	// use std file instead
	std::cout << "will use the std. files instead -> related uncertainty will be 0!!!!!!!" << std::endl;
	foundCanv=false;
      }
      if(!(canvasComb&&canvasTheoComb)){
	if(!canvasComb) std::cout << " canvas " << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " not found in diffXSecTopSemiLep"+closureLabel+LV+PS+".root (event yield combined file)"     << std::endl;
	if(!canvasTheoComb) std::cout << " theory canvas " << xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i] << " not found in diffXSecTopSemiLep"+closureLabel+LV+PS+".root (event yield combined file)" << std::endl;
	foundCanv=false;
      }
      if(foundCanv){
	// get data plots for all systematics
	TString plotName = xSecVariables_[i];
	plotName.ReplaceAll("Norm","");
	TString plotNameTheo=plotName;
	plotNameTheo.ReplaceAll("inclusive","xSec/inclusiveTheory");
	//std::cout << xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i] << "/" << plotNameTheo << std::endl;
	TH1F* plotMu   = combinedEventYields ? 0 : (TH1F*)canvasMu  ->GetPrimitive(plotName+"kData");
	TH1F* plotEl   = combinedEventYields ? 0 : (TH1F*)canvasEl  ->GetPrimitive(plotName+"kData");
	// if already combined at yield level:
	TH1F* plotComb = combinedEventYields ? (TH1F*)canvasComb    ->GetPrimitive(plotName+"kData") : 0;
	// theory
	TH1F* plotTheo = combinedEventYields ? (TH1F*)canvasTheoComb->GetPrimitive(plotNameTheo) : (TH1F*)canvasTheo->GetPrimitive(plotNameTheo);
	if(!plotTheo){
	  if(verbose>0) std::cout << "INFO: can not find " << plotNameTheo << std::endl;
	  plotNameTheo.ReplaceAll("Theory","");
	  plotTheo = combinedEventYields ? (TH1F*)canvasTheoComb->GetPrimitive(plotNameTheo) : (TH1F*)canvasTheo->GetPrimitive(plotNameTheo);
	}
	if(!plotTheo){
	  if(verbose>0) std::cout << "INFO: can not find " << plotNameTheo << std::endl;
	  plotNameTheo.ReplaceAll("xSec/","");
	  plotTheo = combinedEventYields ? (TH1F*)canvasTheoComb->GetPrimitive(plotNameTheo) : (TH1F*)canvasTheo->GetPrimitive(plotNameTheo);
	}
	if(!plotTheo){
	  if(verbose>0) std::cout << "INFO: can not find " << plotNameTheo << std::endl;
	  // take care of differing naming convention for hadron level plots
	  if(hadron) plotTheo = combinedEventYields ? (TH1F*)canvasTheoComb->GetPrimitive(plotName+"Gen") : (TH1F*)canvasTheo->GetPrimitive(plotName+"Gen");
	  // if theory is not found, load parton level theory plots directly
	  else plotTheo = combinedEventYields ? (TH1F*)canvasTheoComb->GetPrimitive("analyzeTop"+LV+"LevelKinematics"+PS+"/"+plotName) : (TH1F*)canvasTheo->GetPrimitive("analyzeTop"+LV+"LevelKinematics"+PS+"/"+plotName);
	  if(verbose>0) std::cout << "INFO: can not find " << "analyzeTop"+LV+"LevelKinematics"+PS+"/"+plotName << std::endl;
	}
	if(!plotTheo&&(plotNameTheo.Contains("Njets")||plotNameTheo.Contains("rhos"))){
	  if(verbose>0) std::cout << "INFO: can not find " << plotNameTheo << std::endl;
	  // take care of differing naming conventions of mixed object analyzer plots
	  if(hadron){
	    TString plotNameTheoAlter=plotNameTheo;
	    if(plotNameTheo.Contains("Njets")) plotNameTheoAlter.ReplaceAll("Njets"  , "Ngenjets"); 	    
	    if(plotNameTheo.Contains("rhos" )) plotNameTheoAlter.ReplaceAll("rhos"   , "rhosGen" );
	      plotTheo = combinedEventYields ? (TH1F*)canvasTheoComb->GetPrimitive(plotNameTheoAlter) : (TH1F*)canvasTheo->GetPrimitive(plotNameTheoAlter); 
	  }
	}
	if(((plotMu&&plotEl)||plotComb)&&plotTheo){ 
	  plotTheo->SetName(plotName);
	  if(verbose>1){
	    if(!combinedEventYields) std::cout << "plot "+plotName+"kData in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " for both channels found!" << std::endl;
	    else std::cout << "plot "+plotName+"kData in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " for combined event yield found!" << std::endl;
	  }

	  if(verbose>0) std::cout << "loading plots done..." << std::endl;
	  // Store information if variable should be drawn normalized
	  bool normalize=xSecVariables_[i].Contains("Norm");
	  
	  // combine the results
	  TH1F* plotCombination= combinedEventYields ? (TH1F*)(plotComb->Clone()) : (TH1F*)(plotMu->Clone());
	  if(!combinedEventYields) plotCombination->Reset("ICESM");
	  if(xSecVariables_[i].Contains("ttbarY")) plotCombination->GetXaxis()->SetRangeUser(-2.49,2.49);
	  if(verbose>1) std::cout << std::endl << xSecVariables_[i] << ", " << sysLabel(sys) << ":" << std::endl;
	  plotCombination->SetLineWidth(3.5);
	  plotCombination->SetMarkerSize(0.8);
	  // check correlation of systematic variation
	  bool correlated=true;
	  for(unsigned int j =0; j<uncorrSys_.size(); ++j){
	    if(uncorrSys_[j]==sys&&sys!=sysNo) correlated=false;
	  }
	  // (i) hadronization uncertainty: use Powheg+Pythia vs Mc@nlo+Herwig 
	  if (sys==sysHadUp||sys==sysHadDown){
	    TH1F* plotMuPow=0;
	    TH1F* plotElPow=0;
	    TH1F* plotMuMca=0;
	    TH1F* plotElMca=0;
	    TH1F* plotCombPow=0;
	    TH1F* plotCombMca=0;
	    if(!combinedEventYields){
	      errorWeightedMeanCombination(*plotMu, * plotEl, *plotCombination, verbose); // use value from MC
	      // A) get histos: Powheg and MC@NLO
	      // canvas
	      TCanvas* canvasMuPow = (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+sysLabel(sysGenPowheg)+"/"+xSecVariables_[i]));
	      TCanvas* canvasElPow = (TCanvas*)(files_[kElectron]->Get(xSecFolder+"/"+sysLabel(sysGenPowheg)+"/"+xSecVariables_[i]));
	      TCanvas* canvasMuMca = (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+sysLabel(sysGenMCatNLO)+"/"+xSecVariables_[i]));
	      TCanvas* canvasElMca = (TCanvas*)(files_[kElectron]->Get(xSecFolder+"/"+sysLabel(sysGenMCatNLO)+"/"+xSecVariables_[i]));
	      // plots
	      if(canvasMuPow&&canvasElPow){
		plotMuPow = (TH1F*)canvasMuPow->GetPrimitive(plotName+"kData");
		plotElPow = (TH1F*)canvasElPow->GetPrimitive(plotName+"kData");  
	      }
	      if(canvasMuMca&&canvasElMca){
		plotMuMca = (TH1F*)canvasMuMca->GetPrimitive(plotName+"kData");
		plotElMca = (TH1F*)canvasElMca->GetPrimitive(plotName+"kData");  
	      }
	      // check if all plots are available
	      if(plotMuPow&&plotElPow&&plotMuMca&&plotElMca){
		// B) combine channels for Powheg and Mcatnlo
		plotCombPow=(TH1F*)(plotMu->Clone());
		plotCombMca=(TH1F*)(plotMu->Clone());
		errorWeightedMeanCombination(*plotMuPow, *plotElPow, *plotCombPow, verbose);
		errorWeightedMeanCombination(*plotMuMca, *plotElMca, *plotCombMca, verbose);
	      }
	    }
	    else{
	      // for combination on yield level: get combined plots
	      TCanvas* canvasCombPow = (TCanvas*)(files_[kCombined]->Get(xSecFolder+"/"+sysLabel(sysGenPowheg )+"/"+xSecVariables_[i]));
	      TCanvas* canvasCombMca = (TCanvas*)(files_[kCombined]->Get(xSecFolder+"/"+sysLabel(sysGenMCatNLO)+"/"+xSecVariables_[i]));
	      if(canvasCombPow) plotCombPow=(TH1F*)canvasCombPow->GetPrimitive(plotName+"kData");
	      if(canvasCombMca) plotCombMca=(TH1F*)canvasCombMca->GetPrimitive(plotName+"kData");
	    }
	    if(plotCombPow&&plotCombMca){
	      // C) transfer absolute uncertainty from Powheg vs. MCatNlo to systematic shifted plot (hadronization)
	      // loop all bins
	      for(int bin=1; bin<=plotCombination->GetNbinsX(); ++bin){
		// central combined value for this bin
		double stdValue = histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin);
		// calculate absolute difference between Powheg and MCatNLO
		double hadUnc = plotCombPow->GetBinContent(bin)-plotCombMca->GetBinContent(bin);
		// take care of sign: hadDown = -hadUp
		if(sys==sysHadDown) hadUnc*=-1;
		// FIXME MARTIN: ignore Hadronization uncertainty in 1st bin because it comes from the difference of the finite and non-finite top mass
		if(xSecVariables_[i].Contains("ttbarMass")&&bin==1) hadUnc=0;
		// printout
		if(verbose>1){
		  std::cout << "bin " << bin << " MadGraph    xSec value: " << stdValue << std::endl;
		  if(plotElPow) std::cout << "      "      << " Powheg el   xSec value: " <<  plotElPow->GetBinContent(bin)   << std::endl;
		  if(plotMuPow) std::cout << "      "      << " Powheg mu   xSec value: " <<  plotMuPow->GetBinContent(bin)   << std::endl;
		  std::cout << "      "      << " Powheg comb xSec value: " <<  plotCombPow->GetBinContent(bin) << std::endl;
		  if(plotElMca) std::cout << "      "      << " mc@nlo el   xSec value: " <<  plotElMca->GetBinContent(bin)   << std::endl;
		  if(plotMuMca) std::cout << "      "      << " mc@nlo mu   xSec value: " <<  plotMuMca->GetBinContent(bin)   << std::endl;
		  std::cout << "      "      << " mc@nlo comb xSec value: " <<  plotCombMca->GetBinContent(bin) << std::endl;
		  std::cout << "bin " << bin << " xSec value: " << stdValue << std::endl;
		  std::cout << "  old (Powheg vs. MadGraph) uncertainty: " << plotCombination->GetBinContent(bin)-stdValue << std::endl;
		  std::cout << "  new (Powheg vs. MC@nlo  ) uncertainty: " << hadUnc << std::endl;
		}
		// save this absolute error wrt. central combined value as final hadronization uncertainty shifted plot
		plotCombination->SetBinContent(bin, stdValue+hadUnc);
	      }
	    }
	    else{ 
	      std::cout << " ERROR: At least one plot not found for hadronization uncertainty (cf. following histo-pointers): " << std::endl;
	      if(!combinedEventYields){
	      std::cout << " Powheg, muon:" << plotMuPow << std::endl;
	      std::cout << " Powheg, elec:" << plotElPow << std::endl;
	      std::cout << " MC@NLO, muon:" << plotMuMca << std::endl;
	      std::cout << " MC@NLO, elec:" << plotElMca << std::endl;
	      }
	      else{
		std::cout << " Powheg, combined:" << plotCombPow << std::endl;
		std::cout << " MC@NLO, combined:" << plotCombMca << std::endl;
	      }
	    }
	  }
	  // (ii) uncorrelated uncertainties 
	  // FIXME: for combination on event yield level all systematics are treated as correlated so far
	  else if(!correlated&&!combinedEventYields){
	    // get standard xSec histos
	    // a) canvas
	    TCanvas* canvasMuStd = (TCanvas*)(files_[kMuon    ]->Get(xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i]));
	    TCanvas* canvasElStd = (TCanvas*)(files_[kElectron]->Get(xSecFolder+"/"+sysLabel(sysNo)+"/"+xSecVariables_[i]));
	    // b) plots
	    if(canvasMuStd&&canvasElStd){
	      TH1F* plotMuStd = (TH1F*)canvasMuStd->GetPrimitive(plotName+"kData");
	      TH1F* plotElStd = (TH1F*)canvasElStd->GetPrimitive(plotName+"kData");    
	      // keep one decay channel constant and vary the other one
	      TH1F* plotMuConstElecVaried=(TH1F*)(plotMu->Clone());
	      errorWeightedMeanCombination(*plotMuStd, * plotEl   , *plotMuConstElecVaried, verbose);
	      TH1F* plotMuVariedElecConst=(TH1F*)(plotMu->Clone());
	      errorWeightedMeanCombination(*plotMu   , * plotElStd, *plotMuVariedElecConst, verbose);
	      // loop all bins
	      for(int bin=1; bin<=plotCombination->GetNbinsX(); ++bin){
		// central combined value for this bin
		double stdValue = histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin);
		// calculate relative shift for both channels wrt. central combined value
		double shiftEl=plotMuConstElecVaried->GetBinContent(bin)-stdValue;
		double shiftMu=plotMuVariedElecConst->GetBinContent(bin)-stdValue;
		// add seperate variations in quadrature
		double totalShift=sqrt(shiftEl*shiftEl+shiftMu*shiftMu);
		// save this absolute error wrt. central combined value as final systematic shifted plot
		plotCombination->SetBinContent(bin, stdValue+totalShift);
	      }
	    }
	  }
	  // (iii)other correlated uncertainties: just do error weighted mean
	  //      for combination on event yield level: do nothing
	  else if(!combinedEventYields) errorWeightedMeanCombination(*plotMu, * plotEl, *plotCombination, verbose);
	  
	  // =================================================
	  //  Additional histos for closure tests
	  // =================================================
	  if(sys==sysNo&&plotName!="inclusive"&&!closureTestSpecifier.Contains("NoDistort")){
	    // get ttbar signal files (with reweighting applied or without for zprime)
	    TString rewfold="/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/";
	    if(reweightClosure) rewfold+="ttbarReweight/";
	    TString muReweighted=rewfold+TopFilename(kSig, sysNo, "muon"    );
	    TString elReweighted=rewfold+TopFilename(kSig, sysNo, "electron");
	    // name extension for reweighted files
	    if(reweightClosure){
	      muReweighted.ReplaceAll("Summer", "SysDistort"+closureTestSpecifier+"Summer");
	      elReweighted.ReplaceAll("Summer", "SysDistort"+closureTestSpecifier+"Summer");
	    }
	    TFile* mufile = new (TFile)(muReweighted);
	    TFile* elfile = new (TFile)(elReweighted);
	    // get plot
	    TString partonPlot="analyzeTop"+LV+"LevelKinematics"+PS+"/"+plotName;
	    if(hadron){
	      // composited objects analyzer
	      if(plotName.Contains("rhos")||plotName.Contains("Njets")){
		partonPlot.ReplaceAll("analyzeTop"     , "composited");
		partonPlot.ReplaceAll("Njets"          , "Ngenjets"  );
		partonPlot.ReplaceAll("rhos"           , "rhosGen"   );
		partonPlot.ReplaceAll("LevelKinematics", "Gen"       );
	      }
	      // bjets analyzer
	      else if(plotName.Contains("bq")||plotName.Contains("lb")||plotName.Contains("bbbar")){
		partonPlot.ReplaceAll("PhaseSpace", "BjetsPhaseSpace");
		partonPlot+="Gen";
	      }
	      // lepton analyzer
	      else if(plotName.Contains("lep")){
		partonPlot.ReplaceAll("PhaseSpace", "LeptonPhaseSpace");
		partonPlot+="Gen";
	      }
	    }
	    histo_["reweighted"+plotName+"Mu"][kSig]=(TH1F*)(mufile->Get(partonPlot)->Clone("rewmu"+plotName));
	    histo_["reweighted"+plotName+"El"][kSig]=(TH1F*)(elfile->Get(partonPlot)->Clone("rewel"+plotName));
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
	      // get files
	      TFile* zprimemufile = new (TFile)(muzprime);
	      TFile* zprimeelfile = new (TFile)(elzprime);
	      // get plots
	      if(zprimemufile) histo_["zprime"+plotName+"Mu"][kSig]=(TH1F*)(zprimemufile->Get(partonPlot)->Clone("zpmu"+plotName));
	      if(zprimeelfile) histo_["zprime"+plotName+"El"][kSig]=(TH1F*)(zprimeelfile->Get(partonPlot)->Clone("zpel"+plotName));
	      // apply lumiweights
	      double zPrimeLumiWeight=zPrimeLumiWeightIni;
	      if     (zprime=="1000") zPrimeLumiWeight=(zPrimeLumiWeight*5*luminosity)/104043;
	      if(histo_["zprime"+plotName+"Mu"].count(kSig)>0) histo_["zprime"+plotName+"Mu"][kSig]->Scale(zPrimeLumiWeight);
	      if(histo_["zprime"+plotName+"El"].count(kSig)>0) histo_["zprime"+plotName+"El"][kSig]->Scale(zPrimeLumiWeight);
	      histo_["reweighted"+plotName+"Mu"][kSig]->Scale(lumiweight(kSig, constLumiMuon, 0, "muon"    ));
	      histo_["reweighted"+plotName+"El"][kSig]->Scale(lumiweight(kSig, constLumiElec, 0, "electron"));
	      // add zprime to signal
	      if(histo_["zprime"+plotName+"Mu"].count(kSig)>0) histo_["reweighted"+plotName+"Mu"][kSig]->Add((TH1F*)histo_["zprime"+plotName+"Mu"][kSig]->Clone());
	      if(histo_["zprime"+plotName+"El"].count(kSig)>0) histo_["reweighted"+plotName+"El"][kSig]->Add((TH1F*)histo_["zprime"+plotName+"El"][kSig]->Clone());
	    }
	    // add channels
	    histo_["closure"+plotName][kSig]=(TH1F*)(histo_["reweighted"+plotName+"Mu"][kSig]->Clone("closure"+plotName));
	    histo_["closure"+plotName][kSig]->Add(   histo_["reweighted"+plotName+"El"][kSig]);
	    // apply standard rebinning
	    std::map<TString, std::vector<double> > binning_ = makeVariableBinning(addCrossCheckVariables);
	    reBinTH1F(*histo_["closure"+plotName][kSig], binning_[plotName], verbose-1);
	    // ensure that under- and overflow bins are empty
	    histo_["closure"+plotName][kSig]->SetBinContent(0                                              , 0.);
	    histo_["closure"+plotName][kSig]->SetBinContent(histo_["closure"+plotName][kSig]->GetNbinsX()+1, 0.);
	    // normalization to unit area
	    histo_["closure"+plotName][kSig]->Scale(1./(histo_["closure"+plotName][kSig]->Integral(0, histo_["closure"+plotName][kSig]->GetNbinsX()+1)));
	    // divide by binwidth
	    histo_["closure"+plotName][kSig] = divideByBinwidth(histo_["closure"+plotName][kSig], verbose-1);
	    // set style
	    histogramStyle(*histo_["closure"+plotName][kSig], kSig, false, 1.2, kMagenta);
	    histo_["closure"+plotName][kSig]->SetLineWidth(3);
	    histo_["closure"+plotName][kSig]->SetLineStyle(2);
	  }

	  // =====================
	  //  Adapt plot styles
	  // =====================
	  
	  double max = plotTheo->GetMaximum();
	  
	  max*=1.5;
	  plotTheo->GetXaxis()->SetNoExponent(true);
	  if(max>1&&max<100) plotTheo->GetYaxis()->SetNoExponent(true);
	  else plotTheo->GetYaxis()->SetNoExponent(false);
	  TString yTitle=plotTheo->GetYaxis()->GetTitle();
	  TString xTitle=plotTheo->GetXaxis()->GetTitle();
	  if(verbose>1) std::cout << "original title: " << yTitle  << std::endl;
	  xTitle.ReplaceAll("#mu","l");
	  yTitle.ReplaceAll("#mu","l");
	  if(verbose>1) std::cout << "new title: " << yTitle  << std::endl;
	  plotTheo->GetXaxis()->SetTitle(xTitle);
	  plotTheo->GetYaxis()->SetTitle(yTitle);
	  // Create canvases
	  TCanvas* combicanvas = new TCanvas("combicanvasStd", "combicanvasStd",  600, 600);
	  // get style from old canvas
	  combicanvas ->SetLogy(combinedEventYields ? canvasComb->GetLogy() : canvasMu->GetLogy());
	  if(normalize){
	    if (pTPlotsLog && xSecVariables_[i].Contains("Pt") ){
	      plotTheo->SetMinimum(0.0001);
	      if      (plotName.Contains("topPt")) plotTheo->SetMaximum(0.02);
	      else if (plotName=="ttbarPt"       ) plotTheo->SetMaximum(0.07);
	      else if (plotName=="lepPt"         ) plotTheo->SetMaximum(0.07);
	      else if (plotName=="bqPt"          ) plotTheo->SetMaximum(0.02);
	      else if (plotName=="bbbarPt"       ) plotTheo->SetMaximum(0.02);
	      combicanvas->SetLogy(1);
	    }
	    else{
	      if      (plotName.Contains("topPt")) plotTheo->SetMaximum(0.01);
	      else if (plotName=="ttbarPt"       ) plotTheo->SetMaximum(0.025);
	      else if (plotName=="lepPt"         ) plotTheo->SetMaximum(0.03);
	      else if (plotName=="bqPt"          ) plotTheo->SetMaximum(0.018);
	      else if (plotName=="bbbarPt"       ) plotTheo->SetMaximum(0.01);
	    }
	    // bq Pt should alwaye be log
	    // if (!pTPlotsLog && xSecVariables_[i].Contains("bqPt") ){
	    //  plotTheo->SetMinimum(0.00001);
	    //  plotTheo->SetMaximum(0.1);
	    //  combicanvas->SetLogy(1); 
	    //}
	    // adjust max
	    if(plotName.Contains("lepEta")||plotName=="bqEta"||plotName.Contains("topY")||plotName=="ttbarY"||plotName.Contains("rhos")) plotTheo->GetYaxis()->SetNoExponent(true);
	    
	    if      (plotName.Contains("topY")  ) plotTheo->SetMaximum(0.7);
	    else if (plotName=="ttbarY"         ) plotTheo->SetMaximum(0.8);
	    else if (plotName.Contains("lepEta")) plotTheo->SetMaximum(0.6);
	    else if (plotName=="bqEta"          ) plotTheo->SetMaximum(0.6);
	    else if (plotName.Contains("rhos")  ) plotTheo->SetMaximum(4.0);
	  }
	  
	  // activate canvas
	  gStyle->SetEndErrorSize(8);
	  
	  combicanvas->cd(0);
	  combicanvas->SetBottomMargin(myStyle.GetPadBottomMargin());
	  combicanvas->SetLeftMargin(myStyle.GetPadLeftMargin());
	  gStyle->SetEndErrorSize(8);
	  
	  // =====================
	  //  Theory Curves
	  // =====================
	  
	  bool DrawMCAtNLOPlot2      = DrawMCAtNLOPlot;
	  bool DrawPOWHEGPYTHIAPlot2 = DrawPOWHEGPYTHIAPlot;
	  bool DrawPOWHEGHERWIGPlot2 = DrawPOWHEGHERWIGPlot;
	  bool smoothcurves2         = smoothcurves;
	  bool DrawSmoothMadgraph2   = DrawSmoothMadgraph;
	  
	  // ============================
	  //    !!!! PLEASE READ !!!!
	  // ============================
	  
	  // NO theory CURVES for inclusive xSec
 	  if(plotName=="inclusive"){
 	    DrawMCAtNLOPlot2=false;
	    DrawPOWHEGPYTHIAPlot2=false;
	    DrawPOWHEGHERWIGPlot2=false;
	    DrawSmoothMadgraph2=false;
	    smoothcurves2=false;
	  }
	  
	  // Adjust drawing parameters for theory curves
	  if (normalize){
	    smoothcurves2=smoothcurves;
	  }
	  else{
	    // no theory curves for absolute xSecs at the moment
	    DrawMCAtNLOPlot2      = false;
	    DrawPOWHEGPYTHIAPlot2 = false;
	    DrawPOWHEGHERWIGPlot2 = false;
	    smoothcurves2         = false;  
	    DrawSmoothMadgraph2   = false;
	  }
	  
	  // a) draw axis of binned MadGraph curve
	  if(xSecVariables_[i].Contains("ttbarY")) plotTheo->GetXaxis()->SetRangeUser(-2.49,2.49);
	  plotTheo->Draw("AXIS");
	  //if(xSecVariables_[i].Contains("Norm")){
	  
	  // b) binned MadGraph curve
	  //plotTheo->Draw("hist same"); // from analyzeHypothesisKinFit.C
	  
	  // get labels for folder extensions
	  TString hadLevelExtend = "";
	  TString hadLevelPlotExtend = "";
	  if(hadron){
	    if     (xSecVariables_[i].Contains("lep")){ hadLevelExtend="Lepton"; hadLevelPlotExtend="Gen"; }
	    else if(xSecVariables_[i].Contains("bq" )||xSecVariables_[i].Contains("lbMass")||xSecVariables_[i].Contains("bbbar")){ hadLevelExtend="Bjets" ; hadLevelPlotExtend="Gen"; }
	  }
	  // b1) create binned MADGRAPH theory curve (std sample w.o. SC)
	  // load it from combined file
	  TString plotNameMadgraph="analyzeTop"+LV+"LevelKinematics"+hadLevelExtend+PS+"/"+plotName+hadLevelPlotExtend;
	  if(xSecVariables_[i].Contains("Njets"    )) plotNameMadgraph="composited"+LV+"Gen"+PS+"/Ngenjets";
	  if(xSecVariables_[i].Contains("rhos"     )) plotNameMadgraph="composited"+LV+"Gen"+PS+"/rhosGen" ;
	  plotNameMadgraph.ReplaceAll("Norm","");
	  TString MGcombFile="/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/"+TopFilename(kSig, 0, "muon").ReplaceAll("muon", "combined");
	  TString plotNameMadgraph2=plotNameMadgraph;
	  if(xSecVariables_[i].Contains("inclusive")){
	    if(extrapolate) plotNameMadgraph2.ReplaceAll("inclusive", "topPt");
	    else            plotNameMadgraph2="analyzeTopHadronLevelKinematicsLeptonPhaseSpace/lepEtaGen";
	  }
	  TH1F* plotTheo2 = getTheoryPrediction(plotNameMadgraph2, MGcombFile);
	  // inclusive cross section
	  if(xSecVariables_[i]=="inclusive"){
	    // get events in PS from top pt
	    double NGenPS=0.5*plotTheo2->Integral(0,plotTheo2->GetNbinsX()+1);
	    NGenPS*=lumiweight(kSig, 0.5*(constLumiElec+constLumiMuon), sysNo, "muon");
	    // get BR
	    TH1F* plotTheo22 = getTheoryPrediction("analyzeTopPartonLevelKinematics/topPt", MGcombFile);
	    double NGen=0.5*plotTheo22->Integral(0,plotTheo22->GetNbinsX()+1);
	    NGen*=lumiweight(kSig, 0.5*(constLumiElec+constLumiMuon), sysNo, "muon");
	    double BR=NGen*1.0/(ttbarCrossSection*0.5*(constLumiElec+constLumiMuon));
	    // get predicted inclusive cross section in chosen PS 
	    double inclxSec=NGenPS/(BR*0.5*(constLumiElec+constLumiMuon));
	    // fill histo
	    plotTheo2 = new TH1F( plotTheo->GetName(), plotTheo->GetTitle(), 1, 0., 1.0);
	    plotTheo2->SetBinContent(1, inclxSec );
	  }
	  else{
	    // other quantities
	    std::map<TString, std::vector<double> > binning_ = makeVariableBinning(addCrossCheckVariables);
	    reBinTH1F(*plotTheo2, binning_[plotName], verbose-1);
	    // Normalization absolute cross sections
	    if(!normalize){
	      // luminosity
	      double luminosity2=luminosity;
	      if(combinedEventYields) luminosity2= ( constLumiElec + constLumiMuon );
	      plotTheo2->Scale(1./(luminosity2));
	      // muon and electron channel are added in the root file
	      plotTheo2->Scale(1./(2));
	      // event weight (for signal it does not matter if muon or electron)
	      plotTheo2->Scale(lumiweight(kSig, luminosity2, sysNo, "muon"));
	      // large sample correction factor for number of events
	      // lumiweight is for v1 but large sample is v1+v2
	      //if(largeMGfile) plotTheo2->Scale(3697693./(59613991.+3697693.)); // large file not used
	      // BR
	      if(extrapolate) plotTheo2->Scale(1./BRPDG(sysNo));
	      if(verbose>1) std::cout << "area from abs diff MC plot: " << plotTheo2->Integral(0,plotTheo2->GetNbinsX()+1) << std::endl;
	    }
	    // divide by binwidth
	    plotTheo2=divideByBinwidth(plotTheo2, verbose-1);
	    // Normalization normalized cross sections
	    double XSecInclTheoPS= getInclusiveXSec(plotTheo2,verbose-1);
	    if(normalize) plotTheo2->Scale(1./(XSecInclTheoPS));
	  }
	  // styling
	  histogramStyle( *plotTheo2, kSig, false);
	  plotTheo2->SetTitle(plotTheo->GetTitle());
	  plotTheo2->SetName (plotTheo->GetName() );

	  // b2) create binned MADGRAPH theory curve (large sample including SC)
	  // load it from combined file
	  TString MGcombFile2="/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/combinedDiffXSecSigSummer12PFLarge.root";
	  if(!largeMGfile) MGcombFile2.ReplaceAll("Large","");
	  TH1F* plotTheo3 = (xSecVariables_[i]=="inclusive") ? new TH1F( TString(plotTheo->GetName())+"incl", TString(plotTheo->GetTitle())+"incl", 1, 0., 1.0) : getTheoryPrediction(plotNameMadgraph2, MGcombFile2);
	  // inclusive cross section
	  if(xSecVariables_[i]=="inclusive"){
	    // get events in PS from top pt
	    double NGenPS=0.5*plotTheo3->Integral(0,plotTheo3->GetNbinsX()+1);
	    NGenPS*=lumiweight(kSig, 0.5*(constLumiElec+constLumiMuon), sysNo, "muon");
	    // get BR
	    TH1F* plotTheo33 = getTheoryPrediction("analyzeTopPartonLevelKinematics/topPt", MGcombFile);
	    double NGen=0.5*plotTheo33->Integral(0,plotTheo33->GetNbinsX()+1);
	    NGen*=lumiweight(kSig, 0.5*(constLumiElec+constLumiMuon), sysNo, "muon");
	    double BR=NGen*1.0/(ttbarCrossSection*0.5*(constLumiElec+constLumiMuon));
	    // get predicted inclusive cross section in chosen PS 
	    double inclxSec=NGenPS/(BR*0.5*(constLumiElec+constLumiMuon));
	    // fill histo
	    plotTheo3 = new TH1F( plotTheo->GetName(), plotTheo->GetTitle(), 1, 0., 1.0);
	    plotTheo3->SetBinContent(1, inclxSec );
	  }
	  else{
	    // other quantities
	    std::map<TString, std::vector<double> > binning_ = makeVariableBinning(addCrossCheckVariables);
	    reBinTH1F(*plotTheo3, binning_[plotName], verbose-1);
	    // Normalization absolute cross sections
	    if(!normalize){
	      // luminosity
	      double luminosity2=luminosity;
	      if(combinedEventYields) luminosity2= ( constLumiElec + constLumiMuon );
	      plotTheo3->Scale(1./(luminosity2));
	      // muon and electron channel are added in the root file
	      plotTheo3->Scale(1./(2));
	      // event weight (for signal it does not matter if muon or electron)
	      plotTheo3->Scale(lumiweight(kSig, luminosity2, sysNo, "muon"));
	      // large sample correction factor for number of events
	      // lumiweight is for v1 but large sample is v1+v2
	      //if(largeMGfile) plotTheo3->Scale(3697693./(59613991.+3697693.)); // large file not used
	      // BR
	      if(extrapolate) plotTheo3->Scale(1./BRPDG(sysNo));
	      if(verbose>1) std::cout << "area from abs diff MC plot: " << plotTheo3->Integral(0,plotTheo3->GetNbinsX()+1) << std::endl;
	    }
	    // divide by binwidth
	    plotTheo3=divideByBinwidth(plotTheo3, verbose-1);
	    // Normalization normalized cross sections
	    double XSecInclTheoPS= getInclusiveXSec(plotTheo3,verbose-1);
	    if(normalize) plotTheo3->Scale(1./(XSecInclTheoPS));
	  }
	  // styling
	  histogramStyle( *plotTheo3, kSig, false);
	  plotTheo3->SetLineStyle(2);
	  plotTheo3->SetLineColor(kMagenta);
	  plotTheo3->SetTitle(TString(plotTheo->GetTitle())+" SC");
	  plotTheo3->SetName (TString(plotTheo->GetName()) +" SC");
	  // c) MCatNLO
	  int smoothFactor=0;
	  int rebinFactor=0;
	  int errorRebinFactor=0; 
	  int errorSmoothFactor=0;
	  
	  TString plotNameMCAtNLO="";
	  double rangeLow  = -1.0;
	  double rangeHigh = -1.0;
	  TString PSlabel = (extrapolate) ? "" : "Vis";
	  // adjust options for standalone external file
	  if     (xSecVariables_[i].Contains("topPt"    )){ smoothFactor = 10; rebinFactor = 10; errorRebinFactor = 10; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"TopPt"  ;}
	  else if(xSecVariables_[i].Contains("topY"     )){ smoothFactor =  5; rebinFactor =  1; errorRebinFactor =  1; errorSmoothFactor =  5; plotNameMCAtNLO="h"+PSlabel+"TopY"   ;}
	  else if(xSecVariables_[i].Contains("ttbarPt"  )){ smoothFactor = 10; rebinFactor =  1; errorRebinFactor =  1; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"TTbarPt";}
	  else if(xSecVariables_[i].Contains("ttbarY"   )){ smoothFactor =  5; rebinFactor =  1; errorRebinFactor =  5; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"TTbarY" ;}
	  else if(xSecVariables_[i].Contains("ttbarMass")){ smoothFactor =  5; rebinFactor =  1; errorRebinFactor =  1; errorSmoothFactor =  5; plotNameMCAtNLO="h"+PSlabel+"TTbarM" ;
	      if(cutTtbarMass){rangeLow=constMassRangeLow; rangeHigh=constMassRangeHigh;}}
	  else if(xSecVariables_[i].Contains("ttbarDelPhi")){ smoothFactor = 1; rebinFactor = 1; errorRebinFactor =  1; errorSmoothFactor = 1; plotNameMCAtNLO="h"+PSlabel+"TTbarDelPhi" ;}
	  else if(xSecVariables_[i].Contains("ttbarPhiStar")){ smoothFactor = 1; rebinFactor = 1; errorRebinFactor =  1; errorSmoothFactor = 1; plotNameMCAtNLO="h"+PSlabel+"TTbarPhiStar" ;}
	  else if(xSecVariables_[i].Contains("lepPt"    )){ smoothFactor = 10; rebinFactor =  1; errorRebinFactor =  0; errorSmoothFactor =  2; plotNameMCAtNLO="h"+PSlabel+"LepPt"  ;}
	  else if(xSecVariables_[i].Contains("lepEta"   )){ smoothFactor = 10; rebinFactor = 20; errorRebinFactor = 20; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"LepEta" ;}
	  else if(xSecVariables_[i].Contains("bqPt"     )){ smoothFactor = 10; rebinFactor =  2; errorRebinFactor =  5; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"BottomPt"  ;}
	  else if(xSecVariables_[i].Contains("bqEta"    )){ smoothFactor = 10; rebinFactor =  2; errorRebinFactor =  5; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"BottomEta" ;
}
	  else if(xSecVariables_[i].Contains("Njets"    )){ smoothFactor = 1; rebinFactor =  1; errorRebinFactor =  1; errorSmoothFactor = 1; plotNameMCAtNLO="h"+PSlabel+"Njets" ;}
	  else if(xSecVariables_[i].Contains("rhos"     )){ smoothFactor = 1; rebinFactor =  1; errorRebinFactor =  1; errorSmoothFactor = 1; plotNameMCAtNLO="h"+PSlabel+"rhos"  ;}
	  else if(xSecVariables_[i].Contains("bbbarPt"  )){ smoothFactor = 10; rebinFactor =  2; errorRebinFactor =  5; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"BBbarPt"  ;}
	  else if(xSecVariables_[i].Contains("bbbarMass")){ smoothFactor = 10; rebinFactor =  2; errorRebinFactor =  5; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"BBbarMass";}
	  else if(xSecVariables_[i].Contains("lbMass"   )){ smoothFactor = 10; rebinFactor =  2; errorRebinFactor =  5; errorSmoothFactor = 10; plotNameMCAtNLO="h"+PSlabel+"lbMass";}
	  else if(DrawMCAtNLOPlot2){
	    std::cout << " ERROR - Unknown variable " << xSecVariables_[i] << std::endl;
	    // close file and delete pointer
	    closeStdTopAnalysisFiles(files_);
	    exit(-3);
	  }
	  // choose mcatnlo file
	  // -> file from standalone external studies for errorBands
	  // -> std simulation file else
	  TString filename="/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/combinedDiffXSecSigMcatnloSummer12PF.root";
	  TString errorBandFilename="/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/ttbarNtupleCteq6m.root";
	  // error bands for MCatNLO curves
	  // plotname for std simulation file = madgraph-plotname from analyzer structure
	  TString plotNameMCAtNLO2="analyzeTop"+LV+"LevelKinematics"+hadLevelExtend+PS+"/"+plotName+hadLevelPlotExtend;
	  if(xSecVariables_[i].Contains("Njets")) plotNameMCAtNLO2="composited"+LV+"Gen"+PS+"/Ngenjets"; 
	  if(xSecVariables_[i].Contains("rhos" )) plotNameMCAtNLO2="composited"+LV+"Gen"+PS+"/rhosGen"; 
	  plotNameMCAtNLO2.ReplaceAll("Norm","");
	  plotNameMCAtNLO2.ReplaceAll("Minus","");
	  plotNameMCAtNLO2.ReplaceAll("Plus","");
	  // for bbbar and lead/sublead quantities
	  // -> no error bands
	  bool errorbands=true; 
 	  if(xSecVariables_[i].Contains("lbMass")||xSecVariables_[i].Contains("bbbar")||xSecVariables_[i].Contains("Lead")||xSecVariables_[i].Contains("PhiStar")||xSecVariables_[i].Contains("DelPhi")||xSecVariables_[i].Contains("topPtTtbarSys")||xSecVariables_[i].Contains("Njets")||xSecVariables_[i].Contains("rhos")){
	    errorbands=false;
	    // check if only external file should be used
 	    // -> then: no mcatnlo curve
 	    if(useOnlyExternalMCatNLOfile) DrawMCAtNLOPlot2=false;
 	  }
	  else{
	    // check if only external file should be used
	    // -> then: use external file and plotname
	    if(useOnlyExternalMCatNLOfile){
	      plotNameMCAtNLO2=plotNameMCAtNLO;
	      filename=errorBandFilename;
	    }
	  }
	  // c1) draw MC@NLO errorband
	  if (DrawMCAtNLOPlot2&&errorbands) DrawTheoryCurve(errorBandFilename, plotNameMCAtNLO, normalize, smoothFactor, rebinFactor, constMcatnloColor, constMcatnloStyle, rangeLow, rangeHigh, errorBandFilename, errorRebinFactor, errorSmoothFactor, verbose-1, true, false, "mcatnlo", smoothcurves2, LV);
	  // b3) draw binned madgraph theory curves
	  plotTheo2->Draw("hist same");
	  if(SC) plotTheo3->Draw("hist same"); 
	  // g1) draw NLO+NNLL curve for topMass (normalized)
	  if(extrapolate && DrawNNLOPlot && (xSecVariables_[i].Contains("ttbarMassNorm"))){
	    TString plotname=xSecVariables_[i];	   
	    plotname.ReplaceAll("Norm", "");
	    TFile  *file = new TFile("/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/AhrensNNLL8TeV.root");
	    TH1F   *newPlotNNLOHisto = (TH1F*)(file->Get(plotname)->Clone(plotname+"nnlo"));
	    if(newPlotNNLOHisto){
	      newPlotNNLOHisto->SetLineColor(constNnloColor2);
	      newPlotNNLOHisto->SetMarkerColor(constNnloColor2);
	      newPlotNNLOHisto->SetLineStyle(constNnloStyle2);
	      newPlotNNLOHisto->SetLineWidth(2);
	      newPlotNNLOHisto->SetMarkerStyle(7);
	      DrawSteps(newPlotNNLOHisto ,"same");
	    }
	  }
	  // c2) draw MC@NLO central curve
	  if (DrawMCAtNLOPlot2) DrawTheoryCurve(filename, plotNameMCAtNLO2, normalize, smoothFactor, rebinFactor, constMcatnloColor, constMcatnloStyle, rangeLow, rangeHigh, false, errorRebinFactor, errorSmoothFactor, verbose-1, false, false, "mcatnlo", smoothcurves2, LV);
	  // d) POWHEG 
	  // configure configuration
	  smoothFactor=0;
	  rebinFactor=0;
	  TString plotNamePOWHEG="analyzeTop"+LV+"LevelKinematics"+hadLevelExtend+PS+"/"+xSecVariables_[i]+hadLevelPlotExtend;
	  if(xSecVariables_[i].Contains("Njets")) plotNamePOWHEG="composited"+LV+"Gen"+PS+"/Ngenjets";
	  if(xSecVariables_[i].Contains("rhos" )) plotNamePOWHEG="composited"+LV+"Gen"+PS+"/rhosGen";
	  plotNamePOWHEG.ReplaceAll("Norm","");
	  plotNamePOWHEG.ReplaceAll("Minus","");	
	  plotNamePOWHEG.ReplaceAll("Plus","");
	  rangeLow  = -1.0;
	  rangeHigh = -1.0;
	  if     (xSecVariables_[i].Contains("topPt"    )){ smoothFactor = 10; rebinFactor=10; }
	  else if(xSecVariables_[i].Contains("topY"     )){ smoothFactor =  3; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("ttbarPt"  )){ smoothFactor = 10; rebinFactor= 2; }
	  else if(xSecVariables_[i].Contains("ttbarY"   )){ smoothFactor =  5; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("ttbarMass")){ smoothFactor =  5; rebinFactor= 1; if(cutTtbarMass){rangeLow=constMassRangeLow; rangeHigh=constMassRangeHigh;}}
	  else if(xSecVariables_[i].Contains("ttbarDelPhi" )){ smoothFactor = 1; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("ttbarPhiStar")){ smoothFactor = 1; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("lepPt"    )){ smoothFactor =  2; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("lepEta"   )){ smoothFactor = 10; rebinFactor= 2; }
	  else if(xSecVariables_[i].Contains("bqPt"     )){ smoothFactor = 10; rebinFactor= 2; }
	  else if(xSecVariables_[i].Contains("bqEta"    )){ smoothFactor = 10; rebinFactor= 2; }
	  else if(xSecVariables_[i].Contains("Njets"    )){ smoothFactor = 1 ; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("rhos"     )){ smoothFactor = 1 ; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("bbbarPt"  )){ smoothFactor = 1 ; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("bbbarMass")){ smoothFactor = 1 ; rebinFactor= 1; }
	  else if(xSecVariables_[i].Contains("lbMass"   )){ smoothFactor = 1 ; rebinFactor= 1; }
	  // d1) draw curve POWHEG+PYTHIA
	  if(DrawPOWHEGPYTHIAPlot2) DrawTheoryCurve("/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/combinedDiffXSecSigPowhegSummer12PF.root", plotNamePOWHEG, normalize, smoothFactor, rebinFactor, constPowhegColor, constPowhegStyle, -1./*rangeLow*/, -1./*rangeHigh*/, false, 1., 1., verbose-1, false, false, "powheg", smoothcurves2, LV);
	  
	  // d2) draw curve POWHEG+HERWIG
	  if(DrawPOWHEGHERWIGPlot2) DrawTheoryCurve("/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/combinedDiffXSecSigPowhegHerwigSummer12PF.root", plotNamePOWHEG, normalize, smoothFactor, rebinFactor, constPowhegColor2, constPowhegStyle2, -1./*rangeLow*/, -1./*rangeHigh*/, false, 1., 1., verbose-1, false, false, "powhegherwig", smoothcurves2, LV);

	  // e/f) reweighted histos/ttbar+zprime for closure test
	  if(closureTestSpecifier!=""&&!closureTestSpecifier.Contains("NoDistort")&&sys==sysNo&&plotName!="inclusive"){
	      histo_["closure"+plotName][kSig]->Draw("hist same");
	  }
	  // g) draw NNLO curve for topPt (normalized) and topY (normalized) and/or MCFM curves
	  if(extrapolate && (DrawNNLOPlot || DrawMCFMPlot)){ 
	    
	    TString plotname=xSecVariables_[i];	   
	    plotname.ReplaceAll("Norm", "");
	    std::map<TString, std::vector<double> > binning_ = makeVariableBinning(addCrossCheckVariables);

	    // g2) draw NNLO curve for topPt (normalized) and topY (normalized)
	    if(DrawNNLOPlot&&(xSecVariables_[i].Contains("topPtNorm")||xSecVariables_[i].Contains("topYNorm"))){
	      TFile  *file = new TFile("/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/kidonakisApproxNNLO8TeV.root");
	      TH1F   *newPlotNNLOHisto = (TH1F*)(file->Get(plotname)->Clone(plotname+"nnlo"));
	      if(newPlotNNLOHisto){
		newPlotNNLOHisto->GetXaxis()->SetRange(0, newPlotNNLOHisto->GetNbinsX()-1);
		newPlotNNLOHisto->SetLineColor(constNnloColor);
		newPlotNNLOHisto->SetLineStyle(constNnloStyle);
		newPlotNNLOHisto->SetLineWidth(2);
		newPlotNNLOHisto->SetMarkerStyle(7);
		newPlotNNLOHisto->SetMarkerColor(constNnloColor);
		DrawSteps(newPlotNNLOHisto ,"same");
	      }
	    }	  
	    // g3) MCFM curves

	    if (DrawMCFMPlot&&xSecVariables_[i].Contains("Norm")){

	      TString mcfmplotname = "";
	      if      (plotname=="topPt")     mcfmplotname = "pt_t";
	      else if (plotname=="topY")      mcfmplotname = "y_t";
	      else if (plotname=="ttbarPt")   mcfmplotname = "pt_tt";
	      else if (plotname=="ttbarY")    mcfmplotname = "y_tt";
	      else if (plotname=="ttbarMass") mcfmplotname = "m_tt";
	      else if (plotname=="lepPt")     mcfmplotname = "pt_l";
	      else if (plotname=="lepEta")    mcfmplotname = "eta_l";

	      if (mcfmplotname!=""){
		TFile *fileMCFM1 = new TFile("/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/diffCrossSections_normalized_tt_bbh_todk_MSTW200_172_172_ful_central.root","READ");
		TFile *fileMCFM2 = new TFile("/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/diffCrossSections_normalized_tt_bbh_todk_MSTW200_172_172_ful_charCon.root","READ");
		
		if (fileMCFM1 || fileMCFM1){
		  TH1F *plotMCFM1= (TH1F*)(fileMCFM1->Get(mcfmplotname)->Clone(plotname+"MCFM1"));
		  TH1F *plotMCFM2= (TH1F*)(fileMCFM2->Get(mcfmplotname)->Clone(plotname+"MCFM2"));
		  
		  if (plotMCFM1){
		    TH1F *newPlotMCFM1 = reBinTH1FIrregularNewBinning(plotMCFM1,binning_[plotname],mcfmplotname,true);
		    
		    newPlotMCFM1->SetName(plotname+"MCFM1"); 
		    newPlotMCFM1->GetXaxis()->SetRange(0, newPlotMCFM1->GetNbinsX()-1);
		    newPlotMCFM1->SetLineColor(kMagenta);
		    newPlotMCFM1->SetLineWidth(2);
		    newPlotMCFM1->Draw("same");
		  }
		  
		  if (plotMCFM2){
		    TH1F *newPlotMCFM2 = reBinTH1FIrregularNewBinning(plotMCFM2,binning_[plotname],mcfmplotname,true);
		    
		    newPlotMCFM2->SetName(plotname+"MCFM2"); 
		    newPlotMCFM2->GetXaxis()->SetRange(0, newPlotMCFM2->GetNbinsX()-1);
		    newPlotMCFM2->SetLineColor(kBlue);
		    newPlotMCFM2->SetLineWidth(1);
		    newPlotMCFM2->Draw("same");
		  }
		}
	      }	
	      else std::cout << " Unknown variable for MCFM plot - skipping this part for " << plotname << std::endl;
	    }
	  }

	  // g) smooth MadGraph Theory curve
	  smoothFactor=0;
	  rebinFactor=0;
	  rangeLow  = -1.0;
	  rangeHigh = -1.0;
	  if     (xSecVariables_[i].Contains("topPt"    )){ smoothFactor = (largeMGfile ? 1 : 10); rebinFactor = (largeMGfile ? 2 : 10); }
	  else if(xSecVariables_[i].Contains("topY"     )){ smoothFactor = (largeMGfile ? 1 : 5 ); rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("ttbarPt"  )){ smoothFactor = (largeMGfile ? 1 : 5 ); rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("ttbarY"   )){ smoothFactor = (largeMGfile ? 0 : 2 ); rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("ttbarMass"  )){ smoothFactor = (largeMGfile ? 1 : 10);}
	  else if(xSecVariables_[i].Contains("ttbarDelPhi" )){ smoothFactor = 1; rebinFactor = 1; }
	  else if(xSecVariables_[i].Contains("ttbarPhiStar")){ smoothFactor = 1; rebinFactor = 1; }
	  else if(xSecVariables_[i].Contains("lepPt"    )){ smoothFactor = 0; rebinFactor =  0; }
	  else if(xSecVariables_[i].Contains("lepEta"   )){ smoothFactor = (largeMGfile ? 0 : 4); rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("bqPt"     )){ smoothFactor = 0; rebinFactor =  0; }
	  else if(xSecVariables_[i].Contains("bqEta"    )){ smoothFactor = 2; rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("Njets"    )){ smoothFactor = 1; rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("rhos"     )){ smoothFactor = 1; rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("bbbarPt"  )){ smoothFactor = 1; rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("bbbarMass")){ smoothFactor = 1; rebinFactor =  1; }
	  else if(xSecVariables_[i].Contains("lbMass"   )){ smoothFactor = 1; rebinFactor =  1; }
	  if(largeMGfile) MGcombFile="/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName+"/combinedDiffXSecSigSummer12PFLarge.root";
	  //std::cout << plotNameMadgraph << std::endl;
	  if(DrawSmoothMadgraph2) DrawTheoryCurve(MGcombFile, plotNameMadgraph, normalize, smoothFactor, rebinFactor, kRed+1, 1, rangeLow, rangeHigh, false, 1., 1., verbose-1, false, false, "madgraph", DrawSmoothMadgraph2, LV);
	  
	  // ==============
	  //  Legend
	  // ==============
	  
	  TLegend *leg = new TLegend();
	  double legTxtSize         = 0.04;
	  //double x1Subtrahend       = 0.30;
	  //double y1SubtrahendFactor = 0.05;
	  if(reweightClosure || zprime!="") {
	    legTxtSize         = 0.03;
	    //x1Subtrahend       = 0.35;
	    //y1SubtrahendFactor = 0.10;
	  }
	  legendStyle(*leg,"",legTxtSize);
	  

	  // a) Legend - Closure test(s)
	  if(closureTestSpecifier!=""&&sys==sysNo&&plotName!="inclusive"){
	    if(!closureTestSpecifier.Contains("NoDistort")){
	      TString closureentry="";	      
	      if(zprime!="")           closureentry="t#bar{t} + "+zprime+" GeV Z'"+zPrimeLumiWeightStr;
	      else if(reweightClosure){
		closureentry="#splitline{Reweighted t#bar{t}}{#scale[0.75]{(MadGraph+Pythia, "+closureTestSpecifier+")}}";   
		// label cosmetics
		if(closureTestSpecifier=="data") closureentry.ReplaceAll(closureTestSpecifier, "topPt#rightarrow"+closureTestSpecifier);
		if(closureentry.Contains("Up"  )){ closureentry.ReplaceAll(closureTestSpecifier, (closureentry.Contains("ttbarMass") ? "harder " : "softer ")+closureTestSpecifier); closureentry.ReplaceAll("Up"  , "");}
		if(closureentry.Contains("Down")){ closureentry.ReplaceAll(closureTestSpecifier, (closureentry.Contains("ttbarMass") ? "softer " : "harder ")+closureTestSpecifier); closureentry.ReplaceAll("Down", "");}
		closureentry.ReplaceAll("topPt"    , "p_{T}^{t}"   );
		closureentry.ReplaceAll("ttbarMass", "m^{t#bar{t}}");
	      }
	      leg->AddEntry(histo_["closure"+plotName][kSig], closureentry, "L");
	    }
	    leg->SetX1NDC(0.587);
	    leg->SetY1NDC(0.6  );
	    leg->SetX2NDC(0.820);
	    leg->SetY2NDC(0.867);
	  }
	  else {
            leg->SetX1NDC(0.587);
            leg->SetY1NDC(0.6  );
            leg->SetX2NDC(0.820);
            leg->SetY2NDC(0.867);
	  }
	  // b) Legend - Data label
	  leg->AddEntry(plotCombination, dataLabel, "LP");
	  
	  // c1) Legend - Theory prediction - MADGRAPH no SC
	  TString nameMADGRAPHcurve=plotName;
	  TString legLabelMADGRAPH=constMadGraphPythiaLabel;
	  if(reweightClosure || zprime!="") legLabelMADGRAPH="#splitline{t#bar{t} "+constMadGraphPythiaLabel+"}{Used for Unfolding}";
	  nameMADGRAPHcurve.ReplaceAll("Norm","");
	  TH1F* madgraphcurve =(TH1F*)combicanvas->GetPrimitive(nameMADGRAPHcurve);
	  if(madgraphcurve) leg->AddEntry(madgraphcurve, legLabelMADGRAPH, "L" );
	  
 	  // c2) Legend - Theory prediction - MADGRAPH SC
	  TString nameMADGRAPHcurve2=plotName+" SC";
	  TString legLabelMADGRAPH2=constMadGraphPythiaLabelSC;
	  nameMADGRAPHcurve2.ReplaceAll("Norm","");
	  TH1F* madgraphcurve2 =(TH1F*)combicanvas->GetPrimitive(nameMADGRAPHcurve2);
	  if(madgraphcurve2) leg->AddEntry(madgraphcurve2, legLabelMADGRAPH2, "L" );

	  // d) Legend - Theory prediction - MCatNLO
	  if(DrawMCAtNLOPlot2){
	    TString curveName=xSecVariables_[i];
	    curveName.ReplaceAll("Norm","");
	    curveName+="MC@NLOerrorBand";
	    TGraphAsymmErrors *mcatnlocurve =(TGraphAsymmErrors*)combicanvas->GetPrimitive(curveName);
	    if(mcatnlocurve){
	      mcatnlocurve->SetLineStyle(constMcatnloStyle);
	      leg->AddEntry(mcatnlocurve, constMcatnloHerwigLabel, "FL");
	      if(verbose>0) std::cout << curveName << " found!" << std::endl;
	    }	    
	    else{
	      if(verbose>0){
		std::cout << "curve " << curveName << " not found. ";
		curveName.ReplaceAll("errorBand","");
		std::cout << "Will add " << curveName << " instead to the legend" << std::endl;
	      }
	      curveName.ReplaceAll("errorBand","");
	      TH1F* mcatnlocurve2 =(TH1F*)combicanvas->GetPrimitive(curveName);
	      if(mcatnlocurve2){
		leg->AddEntry(mcatnlocurve2, constMcatnloHerwigLabel, "L");	
		if(verbose>0) std::cout << "found!" << std::endl;
	      }
	    }
	  }
	  
	  // e1) Legend - Theory prediction - POWHEG+PYTHIA
	  if(DrawPOWHEGPYTHIAPlot2){
	    TString curveName=xSecVariables_[i];
	    curveName.ReplaceAll("Norm","");
	    curveName+="POWHEG";
	    TH1F* powhegcurve =(TH1F*)combicanvas->GetPrimitive(curveName);
	    if(powhegcurve){
	      leg->AddEntry(powhegcurve, constPowhegPythiaLabel, "L");
	      //std::cout << "found!" << std::endl;
	    }
	  }
	  // e2) Legend - Theory prediction - POWHEG+HERWIG
	  if(DrawPOWHEGHERWIGPlot2){
	    TString curveName=xSecVariables_[i];
	    curveName.ReplaceAll("Norm","");
	    curveName+="POWHEGHERWIG";
	    TH1F* powhegherwigcurve =(TH1F*)combicanvas->GetPrimitive(curveName);
	    if(powhegherwigcurve){
	      leg->AddEntry(powhegherwigcurve, constPowhegHerwigLabel, "L");
	      //std::cout << "found!" << std::endl;
	    }
	  }

	  // f) Legend - Theory prediction - NNLO
	  if(DrawNNLOPlot){
	    TH1F* nnlocurve = 0;
	    if      (xSecVariables_[i].Contains("topPtNorm")){
	      nnlocurve =(TH1F*)combicanvas->GetPrimitive("topPtnnlo");
	    }
	    else if (xSecVariables_[i].Contains("topYNorm" )){
	      nnlocurve =(TH1F*)combicanvas->GetPrimitive("topYnnlo"); 
	    }
	    else if (xSecVariables_[i].Contains("ttbarMassNorm" )){
	      nnlocurve =(TH1F*)combicanvas->GetPrimitive("ttbarMassnnlo");
	    }
	    if(nnlocurve){
	      if (xSecVariables_[i].Contains("ttbarMassNorm" )) leg->AddEntry(nnlocurve, constNloNNLLLabelAhrens, "L");
	      else if (xSecVariables_[i].Contains("topPtNorm")||xSecVariables_[i].Contains("topYNorm" )){
		leg->AddEntry(nnlocurve, constNnloLabelKidonakis, "L");
	      }
	    }
	  }

	  // g) Legend - Theory prediction - MCFM
	  if(DrawNNLOPlot){
	    TH1F* mcfmcurve = 0;
	    TString curveName=xSecVariables_[i];
	    curveName.ReplaceAll("Norm","");
	    curveName+="MCFM1";
	    mcfmcurve =(TH1F*)combicanvas->GetPrimitive(curveName);
	    if(mcfmcurve) leg->AddEntry(mcfmcurve, "MCFM", "L");
	  }
	  
	  
	  // h) Plotting, legend and additional labels

	  if(closureTestSpecifier=="") plotCombination->Draw("e X0 same"); 
	  else  plotCombination->Draw("p X0 same"); 
	  leg->Draw("same");  
	  gPad->RedrawAxis(); 
	  DrawCMSLabels(prelim, luminosity, 0.04, (closureTestSpecifier!="" ? true : false), false, false);
	  DrawDecayChLabel("e/#mu + Jets Combined");
	  if(DrawNNLOPlot&&extrapolate){
	    if (xSecVariables_[i].Contains("topPtNorm")) DrawLabel("(arXiv:1205.3453)", leg->GetX1NDC()+0.06, leg->GetY1NDC()-0.025, leg->GetX2NDC(), leg->GetY1NDC(), 12, 0.025);
	    if (xSecVariables_[i].Contains("topYNorm" )) DrawLabel("(arXiv:1205.3453)", leg->GetX1NDC()+0.06, leg->GetY1NDC()-0.025, leg->GetX2NDC(), leg->GetY1NDC(), 12, 0.025);
	    if (xSecVariables_[i].Contains("ttbarMassNorm")) DrawLabel("(arXiv:1003.5827)", leg->GetX1NDC()+0.06, leg->GetY1NDC()-0.025, leg->GetX2NDC(), leg->GetY1NDC(), 12, 0.025);
	  }

	  histo_[xSecVariables_[i]][sys]=(TH1F*)(plotCombination->Clone());
	  
	  // save combined e+mu plot for systematic error calculation afterwards
	  canvas_[xSecVariables_[i]][sys]=(TCanvas*)(combicanvas->Clone());
	  canvas_[xSecVariables_[i]][sys]->SetTitle(xSecVariables_[i]);
	  canvas_[xSecVariables_[i]][sys]->SetName(xSecVariables_[i]);
	  if(save&&sys==sysNo){
	    if(verbose>0) std::cout << "saving" << std::endl;
	    int initialIgnoreLevel=gErrorIgnoreLevel;
	    if(verbose==0) gErrorIgnoreLevel=kWarning;
	    TString outfolder = closureTestSpecifier=="" ? "xSec/" : "closureTest/";
	    combicanvas->Print("./diffXSecFromSignal/plots/combined/"+dataSample+"/"+outfolder+"xSec"+closureLabel+xSecVariables_[i]+universalplotLabel+".eps"); 
	    combicanvas->Print("./diffXSecFromSignal/plots/combined/"+dataSample+"/"+outfolder+"xSec"+closureLabel+xSecVariables_[i]+universalplotLabel+".png");
	    // pdf for all relevant closureTest plots
	    if(closureTestSpecifier!=""){
	      if(verbose>1) std::cout << "saving for closure test" << std::endl;
	      TString closureTestPdfName ="./diffXSecFromSignal/plots/combined/2012/closureTest/xSecs"+closureLabel+LV+PS;
	      int Nnorm=(xSecVariables_.size())/2;
	      if(xSecVariables_[i].Contains("Norm")&&canvas_.count(xSecVariables_[i])>0&&canvas_[xSecVariables_[i]].count(sysNo)>0&&canvas_[xSecVariables_[i]][sysNo]){
		++countNorm;
		if(verbose>1) std::cout << "normalized xSec " << xSecVariables_[i] << ", "<< countNorm << "/" << Nnorm <<  std::endl;	
		if     (countNorm==1)     combicanvas->Print(closureTestPdfName+".pdf(", "pdf");	
		else if(countNorm==Nnorm) combicanvas->Print(closureTestPdfName+".pdf)", "pdf");
		else                      combicanvas->Print(closureTestPdfName+".pdf" , "pdf");
	      }
	    }
	    gErrorIgnoreLevel=initialIgnoreLevel;
	  }
	  // close Canvas
	  combicanvas->Close();
	  //delete combicanvas;
	  delete leg;
	  delete combicanvas;
	}
	else{ 
	  std::cout << " ERROR in bothDecayChannelsCombination.C! " << std::endl;
	  std::cout << " plot " << xSecVariables_[i];
	  if((plotMu&&plotEl)||plotComb) std::cout << "kData";
	  std::cout << " not found in ";
	  std::cout << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " for:" << std::endl;
	  if(!combinedEventYields){
	    if(!plotMu)   std::cout << "muon channel data"     << std::endl;
	    if(!plotEl)   std::cout << "electron channel data" << std::endl;
	  }
	  else{
	    if(!plotComb) std::cout << "combined event yield channel data" << std::endl;
	  }
	  if(!plotTheo) std::cout << "theory"                << std::endl;
	  // close file and delete pointer
	  closeStdTopAnalysisFiles(files_);
	  exit(-2);
	}
      }
    }
  }
  
  // =================================================
  //  Save combined mu+e plots to ROOT-file
  // =================================================

  // saving in rootfile
  // loop variables
  for(unsigned int i=0; i<xSecVariables_.size(); ++i){
    // loop systematic variations
    for(unsigned int sys=sysNo; sys<sysEnd; ++sys){
      if(histo_[xSecVariables_[i]][sys]){
	saveToRootFile("diffXSecTopSemiLep"+closureLabel+LV+PS+".root", canvas_[xSecVariables_[i]][sys], true, verbose,""+xSecFolder+"/"+sysLabel(sys));
      }
    }
  }
  
  // close file and delete pointer
  closeStdTopAnalysisFiles(files_);
}
