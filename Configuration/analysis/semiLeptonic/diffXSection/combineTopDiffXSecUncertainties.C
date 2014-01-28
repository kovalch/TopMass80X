#include "basicFunctions.h"
#include <numeric>

void combineTopDiffXSecUncertainties(double luminosity=19712., bool save=false, unsigned int verbose=0, TString decayChannel="combined", bool extrapolate=true, bool hadron=false, bool addCrossCheckVariables=false, TString closureTestSpecifier="", bool useBCC=false){

  // ============================
  //  Systematic Variations:
  // ============================
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
 
  // ============================
  //  Set Root Style
  // ============================

  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gStyle->SetErrorX(0.5); // needed to keep the functions alive drawn with DrawSteps
  gROOT->SetBatch(); // you dont want to see all canvas

  // ============================
  //  Parameter Configuration
  // ============================
  // name quantity for which you want to see a detailed uncertainty printout
  TString testVar="topPtNorm";//"topPtNorm";
  // if true: for uncertainties with different versions like eff. SF (norm., eta+pt shape) take only maximum of those
  bool takeMaxOfNormAndShape=true;

  // verbose: set detail level of output 
  // 0: no output, 1: std output 2: output for debugging
  // save: enable saving
  // take care that prescaling of muon channel for full 2011 datset was taken into account
  if(luminosity==4980&&decayChannel=="muon"    ) luminosity=constLumiMuon;
  if(luminosity==4955&&decayChannel=="electron") luminosity=constLumiElec;
  // choose phase space
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
  TString universalplotLabel = extrapolate ? "FullPS" : LV+"LvPS";  
  // dataSample: see if its "2010" or "2011" data
  TString dataSample="2012";
  // for closure test if desired
  TString closureLabel = "";
  if(closureTestSpecifier!=""){
    if      (closureTestSpecifier.Contains("NoDistort")) closureLabel = "PseudoData"+closureTestSpecifier;
    else if (closureTestSpecifier.Contains("topPt"   )||
	     closureTestSpecifier.Contains("ttbarMass")) closureLabel = "PseudoDataReweight"+closureTestSpecifier;
    else if (closureTestSpecifier.Contains("data"     )) closureLabel = "PseudoDataReweighttopPt"+closureTestSpecifier;
    else if (closureTestSpecifier.Contains("1000"     )) closureLabel = "PseudoDataZprime"+closureTestSpecifier+"GeV";
    else{
      std::cout << "ERROR: unknown closureTestSpecifier=" << closureTestSpecifier << std::endl;
      exit(0);
    }
  }
  // outputFile: target rootfile
  // NOTE: this must be identical with TString outputFileName 
  // in analyzeHypothesisKinFit.C
  TString outputFile="diffXSecTopSemi";
  if(decayChannel=="muon"    ) outputFile+="Mu"+closureLabel;
  if(decayChannel=="electron") outputFile+="Elec"+closureLabel;
  if(decayChannel=="combined") outputFile+="Lep"+closureLabel;
  outputFile+=LV+PS+".root";
  // define folder where XSec plots are stored
  TString xSecFolder = "xSec";
  // save all plots into the following folder
  TString outputFolder = "./diffXSecFromSignal/plots/"+decayChannel+"/";
  if(dataSample!="") outputFolder+=dataSample;
  if(extrapolate==false&&hadron==false) useBCC=true;
  // No cross-check variables for hadron phase space
  if (hadron) addCrossCheckVariables=false;

  // ==========================================
  //  Basic Printout and Variable Definitions
  // ==========================================
  
  // Define variables for combination (centrally defined in basicFunctions.h)
  // NOTE: these must be identical to those defined in xSecVariables_ in analyzeHypothesisKinFit.C

  std::vector<TString> xSecVariables_;
  int NormxSecs=0;
  // a) top and ttbar quantities
  if(!hadron){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFit    , xSecVariablesKinFit     + sizeof(xSecVariablesKinFit    )/sizeof(TString));
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFitNorm, xSecVariablesKinFitNorm + sizeof(xSecVariablesKinFitNorm)/sizeof(TString));
    NormxSecs+=sizeof(xSecVariablesKinFitNorm)/sizeof(TString);
  }
  // b) lepton and b-jet quantities
  if(hadron||!extrapolate){
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalState    , xSecVariablesFinalState     + sizeof(xSecVariablesFinalState    )/sizeof(TString));
    xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalStateNorm, xSecVariablesFinalStateNorm + sizeof(xSecVariablesFinalStateNorm)/sizeof(TString));
    NormxSecs+=sizeof(xSecVariablesFinalStateNorm)/sizeof(TString);
  }
  // c) cross check variables presently only available for parton level cross-sections
  if (addCrossCheckVariables && !hadron){
    xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVar,     xSecVariablesCCVar     + sizeof(xSecVariablesCCVar    )/sizeof(TString)    );
    xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVarNorm, xSecVariablesCCVarNorm + sizeof(xSecVariablesCCVarNorm)/sizeof(TString));
    NormxSecs+=sizeof(xSecVariablesCCVarNorm)/sizeof(TString);
  }
  xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesIncl,      xSecVariablesIncl      + sizeof(xSecVariablesIncl)/sizeof(TString)     );
  // chose min/max value[%] for relative uncertainty plots
  double errMax=15.0;
  double errMin= 0.0;
  // container for plots, uncertainties and uncertainty distributions
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TH1F*> > relativeUncertainties_;
  std::map< TString, std::map <unsigned int, TH1F*> > relUncertDistributions_;  // 1st: xSec variable, 2nd: type of uncertainty, 3rd: uncertainty values 
  std::map< TString, std::vector<double> > binning_ = makeVariableBinning(addCrossCheckVariables);
  std::map< TString, TH1F*> statUncertaintyDistributions_;
  std::map< TString, TH1F*> sysUncertaintyDistributions_;
  std::map< TString, TH1F*> totalUncertaintyDistributions_;
  // create container to indicate that plots have been found 
  // and therefore systematics will be calculated
  // The following uncertainties are excluded
  std::map<TString, std::map<unsigned int, bool> > calculateError_, considerError_;
  // create container for combined Errors
  std::map<TString, TGraphAsymmErrors*> totalErrors_;
  // map from new index to old index (enum) of syst. unc.
  std::map<int, int> mapNew2OldSysIndex_;

  // parameter printout

  if(verbose>0){
    std::cout << std::endl;
    std::cout << " Executing combineTopDiffXSecUncertainties with " << dataSample << " data" << std::endl;
    std::cout << std::endl;
    std::cout << " Target file: " << outputFile << std::endl;
    std::cout << " Target variables for cross sections:" << std::endl;

    // loop variables
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      std::cout << xSecVariables_[i]; 
      if(i<xSecVariables_.size()-1) std::cout << ", ";
      else std::cout << std::endl;
    }
    std::cout << " Target folder containing cross section plots: " << xSecFolder << std::endl;
    if(save) std::cout << " Final plots will be saved in " << outputFile << " and as .eps in " << outputFolder << std::endl;
  }
  
  // ============================
  //  Bin Center Corrections
  // ============================
  
  std::map<TString, std::vector<double> > correctedCenters_;
  std::map<TString, std::vector<double> > corrCenterErrors_;
  
  if(useBCC){
    // open file
    TString fileNameBCC = "/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/diffXSecBCC8TeV.root";
    std::cout << " Loading bin center corrections from  file " << fileNameBCC << std::endl;
    TFile* BCCfile=TFile::Open(fileNameBCC, "OPEN");
    if(!BCCfile||BCCfile->IsZombie()){
      std::cout << " Corrupt or missing file " << fileNameBCC << std::endl;
      std::cout << " Aborting execution of macro! " << std::endl;
      exit(0);
    }
    // loop variables
    for(unsigned int var=0; var<xSecVariables_.size(); ++var){
      TString candidate=xSecVariables_[var];
      candidate.ReplaceAll("Norm","");
      if(verbose>1) std::cout << "searching for candidate " << "BCC/"+LV+PS+"/"+candidate+"/Template of Vector class" << std::endl;
      // get object 
      TVectorT<double>* temp;
      BCCfile->GetObject("BCC/"+LV+PS+"/"+candidate+"/Template of Vector class", temp);
      if(temp){
	if(verbose>1) std::cout << "-> found!" << endl;
	// convert it to vector
	std::vector<double> tempvec= std::vector<double>( temp->GetMatrixArray(), temp->GetMatrixArray() + temp->GetNrows() );
	// save it in map
	correctedCenters_[candidate]=tempvec;
	// FIXME: as proper errors for x-BCCs do not exist at the moment, use value itself as dummy...
	corrCenterErrors_[candidate]=tempvec;
      }
      else if(!(candidate.Contains("inclusive"))){
	// draw warning if BCC are missing for a specific variable
	std::cout << "WARNING: BCC for "+candidate+" not found" << endl;
	std::cout << "-> will use the the geometrical bin center" << endl;
	// get geometrical mean from bin borders in makeVariableBinning (basicFunctions.h)
	std::vector<double> tempvec;
	// debug output
	if(verbose>1){
	  std::cout << "searching "+candidate+" in binning_" << std::endl;
	  std::cout << "found? " << binning_.count(candidate) << std::endl;
	}
	// loop all bins
	for(int bin=0; bin<(int)binning_[candidate].size()-1; ++bin){
	  tempvec.push_back(0.5*(binning_[candidate][bin]+binning_[candidate][bin+1]));
	}
	// save centers in map
	correctedCenters_[candidate]=tempvec;
	// FIXME: as proper errors for x-BCCs do not exist at the moment, use value itself as dummy...
	corrCenterErrors_[candidate]=tempvec;
      }
    }
    // Output of results
    for (std::map<TString, std::vector<double> >::iterator iter1 = correctedCenters_.begin(); iter1 != correctedCenters_.end(); iter1++ ){
      std::cout << iter1->first << ": " << std::endl;
      int count=0;
      for (std::vector<double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++){
	++count;
	std::cout << "Bin (" << binning_[iter1->first][count-1] << " .. " <<  binning_[iter1->first][count] << " ): " << (*iter2) << std::endl;
      }
    }
    if(verbose>1){
      for (std::map<TString, std::vector<double> >::iterator iter1 = corrCenterErrors_.begin(); iter1 != corrCenterErrors_.end(); iter1++ ){
	std::cout << iter1->first << " errors : ";
	for (std::vector<double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) std::cout << (*iter2) << " ";
	std::cout << std::endl;
      }
    }

    std::cout << " .... Loading bin center corrections finished." << std::endl;
  }

  // ============================
  //  Open rootfile
  // ============================

  TFile* file = TFile::Open(outputFile, "UPDATE");
  if(verbose>1) std::cout << " Opened file "+outputFile << std::endl;

  // Proceed only if file exists and is not broken  
  
  if(!file || file->IsZombie()){
    std::cout << "ERROR: Target file does not exist or is broken!" << std::endl;
  }
  else{

    if(verbose>0) std::cout << " Target file found!" << std::endl;

    // ===================================
    //  Load hadronization uncertainties
    // ===================================
//     // ATTENTION: values are loaded but not used at the moment!!!
//     // container for hadronization uncertainties
//     std::map< TString, TH1F* > hadUnc_;
//     // get file 
//     TFile* hadfile=TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/ttbarNtupleCteq6mHadronUncert.root", "OPEN");
//     if(!hadfile||hadfile->IsZombie()){
//       std::cout << " Corrupt or missing file " << "/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/ttbarNtupleCteq6mHadronUncert.root" << std::endl;
//       std::cout << " Aborting execution of macro! " << std::endl;
//       exit(0);
//     }

//     // loop quantities 
//     for(unsigned int var=0; var<xSecVariables_.size(); ++var){
//       // exclude inclusive cross section
//       // FIXME: exclude bquark quantities for the moment
//       if(!xSecVariables_[var].Contains("inclusive")&&!xSecVariables_[var].Contains("bq")){
// 	TString name=xSecVariables_[var];       
// 	name.ReplaceAll("Norm","");
// 	name.ReplaceAll("Plus","");
// 	name.ReplaceAll("Minus","");
// 	// get plot 
// 	hadUnc_[xSecVariables_[var]]=(TH1F*)(hadfile->Get(name)->Clone());
// 	//calculateError_[xSecVariables_[var]][sysHadUp  ]=true; // done above
// 	//calculateError_[xSecVariables_[var]][sysHadDown]=true; // done above
//       }
//     }
    
    // ===============================================================
    //  Part A: get all data xSec with all systematic variations
    // ===============================================================

    if(verbose>0) std::cout << std::endl << " Part A: Loading Plots" << std::endl;
    // loop systematic variations
    for(unsigned int sys=sysNo; sys<=ENDOFSYSENUM; ++sys){
      TString subfolder= sys==ENDOFSYSENUM ? "sysHadronizationOld" : sysLabel(sys);      
      // loop variables
      for(unsigned int i=0; i<xSecVariables_.size(); ++i){
	// get canvas       
	TCanvas* canvas = (TCanvas*)file->Get(xSecFolder+"/"+subfolder+"/"+xSecVariables_[i]);	
	if(canvas){
	  // get data plot for all systematics
	  TString plotName = xSecVariables_[i];
	  plotName.ReplaceAll("Norm","");
	  TH1F* plot= (TH1F*)canvas->GetPrimitive(plotName+"kData");
	  if(plot){ 
	    if(verbose>1) std::cout << " Plot "+plotName+"kData found in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << std::endl;
	    // go to root directory, keep plot when closing rootfile
	    gROOT->cd();	    
	    histo_[xSecVariables_[i]][sys]=(TH1F*)(plot->Clone());
	    calculateError_[xSecVariables_[i]][sys]=true;
	    considerError_ [xSecVariables_[i]][sys]=calculateError_[xSecVariables_[i]][sys];
	    if(considerError_[xSecVariables_[i]][sys]) considerError_[xSecVariables_[i]][sys]=considerInTotalError(sys);
	    // no pure generator uncertainties considered for final errors
	    //if(sys==sysGenPowheg||sys==sysGenPowhegHerwig||sys==sysGenMCatNLO) considerError_[xSecVariables_[i]][sys]=false;
	    // no hadronisation uncertainties considered for INCLUSIVE cross section
	    if((sys==sysHadUp||sys==sysHadDown) && xSecVariables_[i]=="inclusive") considerError_[xSecVariables_[i]][sys]=false;
	    // sysHadUp   -> Hadronization: POWHEG+PYTHIA cs. POWHEG+HERWIG
	    // sysHadDown -> Generator: MC@NLO+HERWIG vs. POWHEG+HERWIG
	    if(sys==sysHadDown) considerError_[xSecVariables_[i]][sys]=false;
	    if(sys==sysHadUp  ) considerError_[xSecVariables_[i]][sys]=false;
	    // sysHadronizationOld=MC@NLO+HERWIG vs. POWHEG+PYTHIA
	    if(sys==ENDOFSYSENUM-1){ calculateError_[xSecVariables_[i]][ENDOFSYSENUM]=true; considerError_[xSecVariables_[i]][ENDOFSYSENUM]=true; }
	  }
	  else{
	    if(verbose>1) std::cout << "ERROR: Plot " << plotName +"kData not found in "+ xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << std::endl;
	  }
	  
	  // get MC prediction plot without systematic variation
	  if(sys==sysNo){
	    if (xSecVariables_[i]=="inclusive"&&decayChannel!="combined") plotName = "xSec/"+plotName+"Theory";
	    TH1F* MCplot= (TH1F*)canvas->GetPrimitive(plotName);	    
	    if(MCplot){ 
	      if(verbose>1) std::cout << "MC plot "+plotName+" in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " found!" << std::endl;
	      // go to root directory, keep plot when closing rootfile
	      gROOT->cd();
	      histo_[xSecVariables_[i]+"MC"][sys]=(TH1F*)(MCplot->Clone());
	      if(!histo_[xSecVariables_[i]+"MC"][sys]) std::cout << plotName << " buggy" << std::endl;
	      // Prepare plots for uncertainty distributions	     
	      statUncertaintyDistributions_[xSecVariables_[i]] =(TH1F*)(MCplot->Clone());
	      statUncertaintyDistributions_[xSecVariables_[i]] ->Reset("ICESM");
	      sysUncertaintyDistributions_[xSecVariables_[i]]  =(TH1F*)(MCplot->Clone());
	      sysUncertaintyDistributions_[xSecVariables_[i]]  ->Reset("ICESM");
	      totalUncertaintyDistributions_[xSecVariables_[i]]=(TH1F*)(MCplot->Clone());
	      totalUncertaintyDistributions_[xSecVariables_[i]]->Reset("ICESM");
	    }
	    else if(verbose>1) std::cout << "MC plot '"+plotName+"' in "+xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " not found!" << std::endl;
	  }
	}	
	else if(verbose>1) std::cout << "ERROR: canvas " << xSecFolder+"/"+subfolder+"/"+xSecVariables_[i] << " not found!" << std::endl;
      }
    }

    // ===============================================================
    //  Part B: calculate combined systematic+statistic error
    // ===============================================================

    if(verbose>0) std::cout << std::endl << " Part B: calculate systematic errors for all variables and bins" << std::endl;

    // loop variables
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
    int verbose2=verbose;
    if(xSecVariables_[i]=="inclusive"||xSecVariables_[i]==testVar) verbose2+=1;
      if(verbose2>0) std::cout << std::endl << "variable: " << xSecVariables_[i] << std::endl;
      // check if any plot of the chosen variable has been found
      if(calculateError_.count(xSecVariables_[i])>0&&calculateError_[xSecVariables_[i]][sysNo]==true){
	// define object to save uncertainties
	int Nbins=histo_[xSecVariables_[i]][sysNo]->GetNbinsX();
	TGraphAsymmErrors* combinedErrors = new TGraphAsymmErrors(Nbins);
	// get plot without variation
	TH1F* noSysPlot=(TH1F*)histo_[xSecVariables_[i]][sysNo]->Clone();
	// loop bins
	for(int bin=1; bin<=noSysPlot->GetNbinsX(); ++bin){
	  double totalSystematicError = 0;
	  double stdBinXSecValue      = histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin);
	  double binEdgeDown          = histo_[xSecVariables_[i]][sysNo]->GetBinLowEdge(bin);
	  double binEdgeUp            = histo_[xSecVariables_[i]][sysNo]->GetBinLowEdge(bin+1);
	  double MCpredBinVar         = (decayChannel=="combined") ? histo_[xSecVariables_[i]+"MC"][sysNo]->GetBinContent(bin) : -1.0;
	  // jump to bin without vanishing bin content
	  while(stdBinXSecValue==0&&bin<Nbins){
	    if(verbose>1) std::cout << "bin #" << bin << "/" << Nbins << " has bin content == 0, will be skipped!" << std::endl;
	    ++bin;
	    stdBinXSecValue=histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin);
	    binEdgeDown          = histo_[xSecVariables_[i]][sysNo]->GetBinLowEdge(bin);
	    binEdgeUp            = histo_[xSecVariables_[i]][sysNo]->GetBinLowEdge(bin+1);
	  }
	  // if last bin is empty -> stop looping, next variable
	  if(bin==histo_[xSecVariables_[i]][sysNo]->GetNbinsX()&&stdBinXSecValue==0){
	    if(verbose>1){
	      std::cout << "bin #" << bin << "/" << Nbins << " has bin content == 0, will be skipped!" << std::endl;
	      std::cout << "last bin reached" << std::endl;
	    }
	  }
	  // else: continue with calculation
	  else{
	    if(verbose2>0){ 
	      std::cout << std::endl;
	      if(xSecVariables_[i].Contains("inclusive")) std::cout << "channel: " << decayChannel << std::endl;
	      if(verbose>0) std::cout << xSecVariables_[i] << " bin #" << bin << "/" << Nbins << std::endl;
	      if(verbose>0) std::cout << "( range " << binEdgeDown << " .. " << binEdgeUp << ")" << std::endl;
	      std::cout << "std value[pb/binwidth]: " << std::endl;
	      std::cout << "a) data   : " << stdBinXSecValue << std::endl;
	      if(verbose>0) std::cout << "b) MC pred: " << MCpredBinVar << std::endl;     
	    }
	    // create plot that indicates all relative uncertainties
	    // one plot for every variable in every bin
	    unsigned int nSysTypes = 0;
	    unsigned int nSysCnt   = 0;
	    bool setNewLabel       = 0;
	    TString label="";
	    for(unsigned int sys=sysNo+1; sys<=ENDOFSYSENUM; ++sys){
	      label = sys==ENDOFSYSENUM ? "sysHadronizationOld" : sysLabel(sys);
	      //if(sys==sysNo+1) std::cout << std::endl;
	      //std::cout << label << ": ";
	      // Systemtatic uncertainties are symmetrized later:
	      // --> take only one into account if 'up' and 'down' variations exist
	      // --> take in any case into account if not distinguished between 'up' and 'down' 
	      if      (  label.Contains("Up")   || label.Contains("up"))    nSysTypes++;
	      else if (label.Contains("Hadronization")||!(label.Contains("Down") || label.Contains("down"))) nSysTypes++;
	      //std::cout << nSysTypes << std::endl;
	    }
	    nSysTypes++;; // for MC@NLO+Herwig vs. Powheg+Pythia uncertainty
	    TH1F* relSysPlot = new TH1F("relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i], "relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i], nSysTypes+3, 0.5, 0.5+nSysTypes+3); // +3 for total syst., total stat. and total
	    relSysPlot->GetXaxis()->SetLabelSize(0.025);
	    // loop systematic variations
	    for(unsigned int sys=sysNo+1; sys<=ENDOFSYSENUM+1; ++sys){ // +1 for MC@NLO+Herwig vs. Powheg+Pythia uncertainty
	      // set labels for relative uncertainties for bin & variable in map relativeUncertainties_
	      label = sys<ENDOFSYSENUM ? sysLabel(sys).ReplaceAll("sys","") : "HadronizationOld";
	      // up/down variations: keep only up variations (symmetrisation already done above)
	      if      (!(label.Contains("Hadronization"))&&  label.Contains("Up"  )){nSysCnt++; setNewLabel=1; label = label.ReplaceAll("Up",""); mapNew2OldSysIndex_[nSysCnt]=sys;}
	      else if (!(label.Contains("Hadronization"))&&  label.Contains("up"  )){nSysCnt++; setNewLabel=1; label = label.ReplaceAll("up",""); mapNew2OldSysIndex_[nSysCnt]=sys;}
	      // other (non- up/down variations)
	      else if (  label.Contains("Hadronization") ||!(label.Contains("Down") || label.Contains("down"))){nSysCnt++; setNewLabel=1; mapNew2OldSysIndex_[nSysCnt]=sys;}
	      else setNewLabel=0;
	      if (setNewLabel){
		if(label.Contains("Hadronization")){
		  if (label.Contains("Up")){
		    label="Hadronization";
		  }
		  else if (label.Contains("Down")){
		    label="NLO Generator";
		  }
		}
		relSysPlot->GetXaxis()->SetBinLabel(nSysCnt,label);
		setNewLabel=0;
	      }
	      
	      // create plot that indicates the relative systematic uncertainty
	      double sysDiff=0;
	      if(verbose2>0){
		TString tempi = (sys==ENDOFSYSENUM ? "sysHadronizationOld" : (sys==sysHadDown ? "sysGenerator" : ( sys==sysHadUp ? "sysHadronization" : sysLabel(sys) ) ) );
		std::cout << tempi;
		// check if chosen systematic variation should be considered
		if(calculateError_[xSecVariables_[i]].count(sys)>0&&considerError_[xSecVariables_[i]].count(sys)>0){
		  if(considerError_[xSecVariables_[i]][sys]==true&&calculateError_[xSecVariables_[i]][sys]==true) std::cout << "(considered): ";
		  else if (considerError_[xSecVariables_[i]][sys]==false&&calculateError_[xSecVariables_[i]][sys]==true) std::cout << "(NOT considered): ";
		  else  std::cout << "(NOT calculated): ";
		}
	      }
	      // check if chosen systematic variation of chosen variable has been found and std value is nonzero
	      if(calculateError_[xSecVariables_[i]].count(sys)<=0||stdBinXSecValue==0){
		if(verbose2>0){
		  if(stdBinXSecValue!=0) std::cout << ": not found" << std::endl;
		  else std::cout << ": std value is 0" << std::endl;
		}
	      }
	      else{
		// check if variable has to be considerered for total systematic error
		if(calculateError_[xSecVariables_[i]][sys]==true){
		  double sysBinXSecValue=sys<ENDOFSYSENUM ? histo_[xSecVariables_[i]][sys]->GetBinContent(bin) : stdBinXSecValue;
		  // FIXME TEST MARTIN: chose some fixed test value for empty values
		  //if(sysBinXSecValue==stdBinXSecValue){
		  //  double unc =0.08;
		  //  if(sysLabel(sys).Contains("up")||sysLabel(sys).Contains("Up")) sysBinXSecValue=stdBinXSecValue*(1.+unc);
		  //  else sysBinXSecValue=stdBinXSecValue*(1.-unc);
		  //}
		  sysDiff=std::abs(sysBinXSecValue-stdBinXSecValue);
		  // hadron lv PS lepton and b-jet PDF uncertainties
		  if(!extrapolate&&hadron&&(sys==sysPDFUp||sys==sysPDFDown)&&(xSecVariables_[i].Contains("lep")||xSecVariables_[i].Contains("bq"))){
		    if(verbose>1) std::cout << "load unc PDF for " << xSecVariables_[i] << " bin " << bin << std::endl;
		    TString fileName="/afs/naf.desy.de/group/cms/scratch/tophh/tmp/2012/diffXSecTopSemi";
		    if(decayChannel=="muon"    ) fileName+="Mu";
		    else if(decayChannel=="electron") fileName+="Elec";
		    else if(decayChannel=="combined") fileName+="Lep";
		    fileName+="PartonPhaseSpace.root";   
		    TString canvName="relativeUncertainties/"+xSecVariables_[i]+"/relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i];
		    TString plotName="relSysPlotBin"+getTStringFromInt(bin)+xSecVariables_[i];
		    int sysBin=23; // this is the bin in the symmetrised histo plot
		    sysDiff=stdBinXSecValue*getValue(fileName, canvName, plotName, sysBin)/100.;
		  }
		  // adjust NEW hadronization uncertainty by hand
		  // sysHadUp: as powheg+herwig vs powheg+pythia
		  // sysHadDown: as powheg+herwig vs mcatnlo+herwig
		  if(sys==sysHadUp&&calculateError_[xSecVariables_[i]][sysGenPowheg        ]==true&&calculateError_[xSecVariables_[i]][sysGenPowhegHerwig]==true){
		    sysDiff=std::abs(histo_[xSecVariables_[i]][sysGenPowheg      ]->GetBinContent(bin)-histo_[xSecVariables_[i]][sysGenPowhegHerwig]->GetBinContent(bin)); 
		  }
		  if(sys==sysHadDown&&calculateError_[xSecVariables_[i]][sysGenPowhegHerwig]==true&&calculateError_[xSecVariables_[i]][sysGenMCatNLO     ]==true){
		    sysDiff=std::abs(histo_[xSecVariables_[i]][sysGenPowhegHerwig]->GetBinContent(bin)-histo_[xSecVariables_[i]][sysGenMCatNLO     ]->GetBinContent(bin)); 
		  }		  
		  // adjust OLD hadronization uncertainty by hand
		  // powheg+pythia vs mcatnlo+herwig
		  if(sys==ENDOFSYSENUM&&calculateError_[xSecVariables_[i]][sysGenPowheg        ]==true&&calculateError_[xSecVariables_[i]][sysGenMCatNLO]==true){
		    sysDiff=std::abs(histo_[xSecVariables_[i]][sysGenPowheg      ]->GetBinContent(bin)-histo_[xSecVariables_[i]][sysGenMCatNLO]->GetBinContent(bin)); 
		  }		  
		  // rescale top mass to 0.9 GeV - use rescale factor "SF_TopMass*Uncertainty" defined in basicFunctions.h
		  if      (sys==sysTopMassUp   ) sysDiff *= SF_TopMassUpUncertainty  ;   
		  else if (sys==sysTopMassUp2  ) sysDiff *= (3./6. )*SF_TopMassUpUncertainty;  
		  else if (sys==sysTopMassUp3  ) sysDiff *= (3./9. )*SF_TopMassUpUncertainty;  
		  else if (sys==sysTopMassUp4  ) sysDiff *= (3./12.)*SF_TopMassUpUncertainty; 
		  else if (sys==sysTopMassDown ) sysDiff *= SF_TopMassDownUncertainty; 
		  else if (sys==sysTopMassDown2) sysDiff *= (3./6. )*SF_TopMassDownUncertainty; 
		  else if (sys==sysTopMassDown3) sysDiff *= (3./9. )*SF_TopMassDownUncertainty; 
		  else if (sys==sysTopMassDown4) sysDiff *= (3./11.)*SF_TopMassDownUncertainty; 
		}

		// save relative systematic uncertainties for bin & variable
		// a) non- symmetrized errors (for MC generator based uncertainties and unfolding)
		if (sys==sysGenPowheg || sys==sysGenPowhegHerwig || sys == sysGenMCatNLO || sys==sysUnf || sys==sysHadUp || sys==sysHadDown || sys==ENDOFSYSENUM) relSysPlot->Fill(nSysCnt, 100.0*sysDiff/stdBinXSecValue); 
		// b) weight 0.5 due to error symmetrization
		else relSysPlot->Fill(nSysCnt, 100.0*0.5*sysDiff/stdBinXSecValue);
		// print single systematic uncertainty absolut and relative for bin & variable
		if(verbose2>0) {
		  std::cout << sysDiff << " ( = " << 100*sysDiff/stdBinXSecValue << "% )" << std::endl;
		  if(sys>sysNo && sys <sysHadUp && sys%2==0){
		    std::cout << "Symmetrised:             ----- " << relSysPlot->GetBinContent(nSysCnt)/100.0*stdBinXSecValue << " ( = " << relSysPlot->GetBinContent(nSysCnt)<< "% )" << std::endl;
		  }
		}
	      }
	      // for last systematic 
	      if(sys==ENDOFSYSENUM){
		// calculate total systematic uncertainty
		double maxLepEffSF=0;
		double maxBtagSF=0;
		for (unsigned int n=1; n<=nSysTypes; n++){
		  // exclude uncertainties from total uncertainties and put name in brakets if it is not intended to be considered
		  if (considerError_[xSecVariables_[i]][ mapNew2OldSysIndex_[n] ]==false){
		    relSysPlot->GetXaxis()->SetBinLabel(n, TString("(")+relSysPlot->GetXaxis()->GetBinLabel(n)+")");

		    if(verbose>0) std::cout << " Uncertainties when unfolding with "<< relSysPlot->GetXaxis()->GetBinLabel(n) <<" is excluded from total systematic uncertainty. Idx= "<< n << std::endl;
		  }
		  // add errors quadratically for other uncertainties
		  else {
		    TString tempSysLabel = (TString) relSysPlot->GetXaxis()->GetBinLabel(n);
		    double tempErr = stdBinXSecValue*relSysPlot->GetBinContent(n)/100.0;
		    // for uncertainties with different versions like eff. SF (norm., eta+pt shape) take only maximum of those
		    if(takeMaxOfNormAndShape){
		      if(tempSysLabel=="LepEffSFNorm" || tempSysLabel=="LepEffSFShapeEta" || tempSysLabel=="LepEffSFShapePt"){
			maxLepEffSF=TMath::Max(maxLepEffSF, tempErr);
			if(verbose>0) std::cout << "   Search for the maximum of LepEffSF; now: "<< tempSysLabel << "; with err= " << tempErr << "; maxLepEffSF= " << maxLepEffSF <<std::endl;
			if(tempSysLabel=="LepEffSFShapePt") tempErr=maxLepEffSF;
			else                                tempErr=0.;
		      }
		      else if(tempSysLabel=="BtagSF" || tempSysLabel=="BtagSFShapePt65" || tempSysLabel=="BtagSFShapeEta0p7"){
			maxBtagSF=TMath::Max(maxBtagSF, tempErr);
			if(verbose>0) std::cout << "   Search for the maximum of BtagSF; now: "<< tempSysLabel << "; with err= " << tempErr << "; maxBtagSF= " << maxBtagSF <<std::endl;
			if(tempSysLabel=="BtagSFShapeEta0p7") tempErr=maxBtagSF;
			else                                  tempErr=0.;
		      }
		    }
		    // add finally quadratically!!!
		    if(verbose>0) std::cout << "The following uncertainty is added: "<< tempSysLabel << "; Index= " << n << "; Unc.= " <<tempErr << std::endl;
		    totalSystematicError += pow(tempErr,2);
		  }
		}
		totalSystematicError = sqrt(totalSystematicError);
		// go to root directory keep plot when closing rootfile
		gROOT->cd();
		// finally save relative uncertainties in map relativeUncertainties_
		relativeUncertainties_[xSecVariables_[i]][bin]=(TH1F*)relSysPlot->Clone();		
		// get statistical and total uncertainties
		double statErrorBinVar     = noSysPlot->GetBinError(bin);
		double combinedErrorBinVar = sqrt(totalSystematicError*totalSystematicError + statErrorBinVar*statErrorBinVar);
		// print statistical, systematic and total uncertainties, absolut and relative for bin & variable
		if(verbose2>0){
		  std::cout << "total statistic uncertainty:  " << " +/- " << statErrorBinVar      << " (" << 100*statErrorBinVar/stdBinXSecValue      << "%)" << std::endl;
		  std::cout << "total systematic uncertainty: " << " +/- " << totalSystematicError << " (" << 100*totalSystematicError/stdBinXSecValue << "%)" << std::endl;
		  std::cout << "total uncertainty:            " << " +/- " << combinedErrorBinVar  << " (" << 100*combinedErrorBinVar/stdBinXSecValue  << "%)" << std::endl;
		}
		
		// print MC prediction value and difference in std variations
		double xSecDiff=stdBinXSecValue-MCpredBinVar;
		double absError=combinedErrorBinVar;
		double relativeDifference=xSecDiff/absError; 
		if(verbose>0){
		  std::cout << "difference data-MC = " << xSecDiff << " pb" << std::endl;
		  std::cout << " (" << (xSecDiff/stdBinXSecValue)*100 << "%)" << std::endl;
		  std::cout << " (" << relativeDifference << " std variations)" << std::endl;
		}
		// save relative uncertainties for bin & variable in map relativeUncertainties_
		// a) statistic uncertainty
		relativeUncertainties_[xSecVariables_[i]][bin]   -> SetBinContent(nSysTypes+1, 100*statErrorBinVar/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]   -> GetXaxis()->SetBinLabel(nSysTypes+1, "statistical");
		// b) total systematic uncertainty
		relativeUncertainties_[xSecVariables_[i]][bin]   -> SetBinContent(nSysTypes+2, 100*totalSystematicError/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]   -> GetXaxis()->SetBinLabel(nSysTypes+2, "total systematic");
		// c) total stat+systematic uncertainty
		relativeUncertainties_[xSecVariables_[i]][bin]   -> SetBinContent(nSysTypes+3, 100*combinedErrorBinVar/stdBinXSecValue);
		relativeUncertainties_[xSecVariables_[i]][bin]   -> GetXaxis()->SetBinLabel(nSysTypes+3, "total");
		// d) differential uncertainty distributions
		statUncertaintyDistributions_[xSecVariables_[i]] -> SetBinContent(bin,100*statErrorBinVar/stdBinXSecValue);
		sysUncertaintyDistributions_[xSecVariables_[i]]  -> SetBinContent(bin,100*totalSystematicError/stdBinXSecValue);
		totalUncertaintyDistributions_[xSecVariables_[i]]-> SetBinContent(bin,100*combinedErrorBinVar/stdBinXSecValue);		
		// set combined errors for final xSec plot
		double pointXValue = histo_[xSecVariables_[i]][sysNo]->GetBinCenter(bin);
		double pointXError = 0;
		if(useBCC&&decayChannel=="combined"&&xSecVariables_[i]!="inclusive"){
		  TString plotName=xSecVariables_[i];
		  plotName.ReplaceAll("Norm","");
		  if(verbose>1){ 
		    std::cout << "try to access BCC result map for " << plotName << " (bin " << bin << ")" << std::endl;
		    for(std::map<TString, std::vector<double> >::const_iterator val=correctedCenters_.begin(); val!=correctedCenters_.end(); ++val){
		      std::cout << "found entry " << val->first << " (size " << val->second.size() << ", ";
		      std::cout << "center (bin " << bin << ") " << ((int)val->second.size() >= bin ? val->second.at(bin-1) : -1.) << ")" << std::endl;
		    }
		  }
		  if(correctedCenters_.find(plotName)!=correctedCenters_.end()) pointXValue = correctedCenters_[plotName].at(bin-1);
		  else{
		    std::cout << "ERROR: no entry for " << plotName << " found in correctedCenters_= b.getMapWithCorrectedCentersInX()" << std::endl;
		    exit(0);
		  }
		  if(corrCenterErrors_.find(plotName)!=corrCenterErrors_.end()) pointXError = corrCenterErrors_[plotName].at(bin-1);
		  else{
		    std::cout << "ERROR: no entry for " << plotName << " found in corrCenterErrors_=b.getMapWithCenterErrorsInX()" << std::endl;
		      exit(0);
		  }
		  pointXError=0;
		}
		combinedErrors->SetPoint(0, 0, -1000);
		combinedErrors->SetPoint(bin, pointXValue, histo_[xSecVariables_[i]][sysNo]->GetBinContent(bin));
		combinedErrors->SetPointError(bin, pointXError, pointXError, combinedErrorBinVar, combinedErrorBinVar);
		// define style for relative error plots, uncertainties are only plotted positively after implementing the symmetrization
		histogramStyle(*relativeUncertainties_[xSecVariables_[i]][bin], kSig, true, 2.0, kBlack); 
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->LabelsOption("v");
		relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetLabelSize(0.04);
		double histMax = double(roundToInt(1.1*relativeUncertainties_[xSecVariables_[i]][bin]->GetMaximum()));
		relativeUncertainties_[xSecVariables_[i]][bin]->SetMaximum(histMax);
		relativeUncertainties_[xSecVariables_[i]][bin]->SetMinimum(errMin );
		relativeUncertainties_[xSecVariables_[i]][bin]->SetStats(kFALSE);
		relativeUncertainties_[xSecVariables_[i]][bin]->GetYaxis()->SetTitle("Relative Uncertainty (symmetrized) [%]");
		relativeUncertainties_[xSecVariables_[i]][bin]->GetYaxis()->SetTitleOffset(1.2);
		// save asymmetric errors in map totalErrors_
		combinedErrors->SetLineWidth(3.5);
		combinedErrors->SetMarkerSize(0.8);
		combinedErrors->SetLineColor(histo_[xSecVariables_[i]][sysNo]->GetLineColor());
		totalErrors_[xSecVariables_[i]]=(TGraphAsymmErrors*)(combinedErrors->Clone());
		whipEmptyBinsAway(totalErrors_[xSecVariables_[i]], verbose);
		// remove ugly first point in m(ttbar)
		//if(xSecVariables_[i].Contains("ttbarMass")){
		//totalErrors_[xSecVariables_[i]]->SetPoint(1, 0, -1000);
		//}
	      }
	    }
	    // std::cout << mcatnloIdx  << std::endl;
	    // std::cout << powhegIdx   << std::endl;
	    delete relSysPlot;
	  }
	}
      }
    }
    // ======================================================================
    //  Create plots to show distribution of uncertainties per xSec variable
    // ======================================================================

    // prepare complete uncertainty table

    std::map<TString, double> minMedian;
    std::map<TString, double> maxMedian;
    std::map<TString, double> aveMedian;

    TString uncList = "";
    for (int uncIdx=1; uncIdx<=ENDOFSYSENUM;){

      TString arrayLabelIds[] = {"Up","up","Down","down"};
      TString upTypeLabel   = uncIdx  >=ENDOFSYSENUM+1 ? "" : ( uncIdx  ==ENDOFSYSENUM ? "sysHadronizationOld" : (uncIdx  ==sysHadDown ? "sysGenerator" : ( uncIdx  ==sysHadUp ? "sysHadronization" : sysLabel(uncIdx  ) )) );
      TString downTypeLabel = uncIdx+1>=ENDOFSYSENUM+1 ? "" : ( uncIdx+1==ENDOFSYSENUM ? "sysHadronizationOld" : (uncIdx+1==sysHadDown ? "sysGenerator" : ( uncIdx+1==sysHadUp ? "sysHadronization" : sysLabel(uncIdx+1) )) ); // +1 to get next entry
      
      for (unsigned int j=0; j<(sizeof(arrayLabelIds)/sizeof(arrayLabelIds[0])); j++){
	upTypeLabel.ReplaceAll(arrayLabelIds[j],"");
	downTypeLabel.ReplaceAll(arrayLabelIds[j],"");	    
      }
      
      TString label = upTypeLabel;
      std::cout << label << std::endl;
      (upTypeLabel==downTypeLabel) ? uncIdx+=2 : uncIdx++;

      uncList+=" & ";
      if(uncIdx!=ENDOFSYSENUM&&(!considerInTotalError(uncIdx)||uncIdx==sysHadDown||uncIdx==sysHadUp)) uncList+= "%";
      uncList+=label;
    }
    uncList.ReplaceAll("sys","");
    if(save) writeToFile(uncList, outputFolder+"/uncertaintyDistributionsOverview/uncertaintyTable_"+decayChannel+"_"+universalplotLabel+".txt", false);

    for(unsigned int i=0; i<xSecVariables_.size(); ++i){

      // Distributions of statistical, systematic, and total uncertainties

      std::map<TString, TH1F*>::const_iterator statIter  = statUncertaintyDistributions_.find(xSecVariables_[i]);
      std::map<TString, TH1F*>::const_iterator sysIter   = sysUncertaintyDistributions_.find(xSecVariables_[i]);
      std::map<TString, TH1F*>::const_iterator totalIter = totalUncertaintyDistributions_.find(xSecVariables_[i]);

      //if (statIter  == statUncertaintyDistributions_.end()) {std::cout << "1" << std::endl;}
      //if (sysIter   == sysUncertaintyDistributions_.end())  {std::cout << "2" << std::endl;}
      //if (totalIter == totalUncertaintyDistributions_.end()){std::cout << "3" << std::endl;}

      if (statIter == statUncertaintyDistributions_.end() || sysIter == sysUncertaintyDistributions_.end() || totalIter == totalUncertaintyDistributions_.end() ){
      	std::cout << " Variable '" << xSecVariables_[i] << "' not found in statistical, systematic or total uncertainty distribution." << std::endl;
      }
      else{
	std::cout << " Building statistical, systematic or total uncertainty distributions for '" << xSecVariables_[i] << "'." << std::endl;

	TCanvas *canvasUncertaintyDistributions = new TCanvas("canvasUncertaintyDistributions","canvasUncertaintyDistributions",800,600);
	canvasUncertaintyDistributions->cd();

	totalUncertaintyDistributions_[xSecVariables_[i]]->SetLineColor(38);
	totalUncertaintyDistributions_[xSecVariables_[i]]->SetMinimum(errMin);
	totalUncertaintyDistributions_[xSecVariables_[i]]->SetMaximum(errMax);
	totalUncertaintyDistributions_[xSecVariables_[i]]->GetYaxis()->SetTitle("Relative Uncertainty [%]");
	
	statUncertaintyDistributions_[xSecVariables_[i]] ->SetLineColor(2);
	sysUncertaintyDistributions_ [xSecVariables_[i]] ->SetLineColor(kOrange-3);

	TLegend *uncDistLegend = new TLegend();
	legendStyle(*uncDistLegend,"");
	uncDistLegend->AddEntry(totalUncertaintyDistributions_[xSecVariables_[i]],"Total Uncertainty"      ,"L");
	uncDistLegend->AddEntry(statUncertaintyDistributions_ [xSecVariables_[i]],"Statistical Uncertainty","L");
	uncDistLegend->AddEntry(sysUncertaintyDistributions_  [xSecVariables_[i]],"Systematic Uncertainty" ,"L");
	uncDistLegend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.35);
	uncDistLegend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength() - (double)uncDistLegend->GetNRows()*0.04);
	uncDistLegend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
	uncDistLegend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()   - gStyle->GetTickLength());

	canvasUncertaintyDistributions->cd();
	canvasUncertaintyDistributions->SetTitle(xSecVariables_[i]);
	totalUncertaintyDistributions_[xSecVariables_[i]]->Draw();
	statUncertaintyDistributions_ [xSecVariables_[i]]->Draw("SAME");
	sysUncertaintyDistributions_  [xSecVariables_[i]]->Draw("SAME");
	uncDistLegend                                    ->Draw("SAME");

	int initialIgnoreLevel=gErrorIgnoreLevel;
	if(verbose==0) gErrorIgnoreLevel=kWarning;
	// a) save to rootFile
	if(save) saveToRootFile(outputFile, canvasUncertaintyDistributions, true, verbose, "uncertaintyDistributionsOverview");
	// b) save as eps and png
	if(save){
	  TString saveName=outputFolder+"/uncertaintyDistributionsOverview/relativeUncertainties"+xSecVariables_[i];
	  if(decayChannel=="combined") saveName+="Combined";
	  saveName+=universalplotLabel;
	  canvasUncertaintyDistributions->Print(saveName+".eps");
	  canvasUncertaintyDistributions->Print(saveName+".png");
	}
	gErrorIgnoreLevel=initialIgnoreLevel;

	delete canvasUncertaintyDistributions; canvasUncertaintyDistributions=NULL;
      }

      // Distribution for each uncertainty

      std::map<TString, std::map<unsigned int, TH1F*> >::const_iterator outerIter = relativeUncertainties_.find(xSecVariables_[i]);

      if (outerIter == relativeUncertainties_.end()){
	std::cout << " Variable '" << xSecVariables_[i] << "' not found to build uncertainty distributions." << std::endl;
      }
      else{
	std::cout << " Building single uncertainty distributions for '" << xSecVariables_[i] << "'." << std::endl;

	TString plotNameForBinning = xSecVariables_[i];
	plotNameForBinning.ReplaceAll("Norm","");       
	int NBINS = binning_[plotNameForBinning].size()-1;

	//std::map<unsigned int, TH1F*>::const_iterator innerIter = outerIter->second.begin();
	//int NSYSTYPES   = innerIter->second->GetNbinsX() - 3;	
	int histoBinIdx = 0;

	// prepare uncertainty table with minimum and maximum

	if(universalplotLabel!="PartonLvPS"&&save) writeToFile("Type of unc. & average Relative unc. & min Relative unc. & max Relative unc. \\\\", outputFolder+"/uncertaintyDistributionsOverview/uncertaintyTable_"+decayChannel+"_"+universalplotLabel+"_minmax.txt", false);
	uncList = "";
	uncList+=xSecVariables_[i];

	for (int uncIdx=1; uncIdx<=ENDOFSYSENUM;){ // loop starts at 1 to skip sysNo   

	  histoBinIdx++;

	  TString arrayLabelIds[] = {"Up","up","Down","down"};
	  TString upTypeLabel   = uncIdx  >=ENDOFSYSENUM+1 ? "" : ( uncIdx  ==ENDOFSYSENUM ? "sysHadronizationOld" : (uncIdx  ==sysHadDown ? "sysGenerator" : ( uncIdx  ==sysHadUp ? "sysHadronization" : sysLabel(uncIdx  ) )) );
	  TString downTypeLabel = uncIdx+1>=ENDOFSYSENUM+1 ? "" : ( uncIdx+1==ENDOFSYSENUM ? "sysHadronizationOld" : (uncIdx+1==sysHadDown ? "sysGenerator" : ( uncIdx+1==sysHadUp ? "sysHadronization" : sysLabel(uncIdx+1) )) ); // +1 to get next entry

	  for (unsigned int j=0; j<(sizeof(arrayLabelIds)/sizeof(arrayLabelIds[0])); j++){
	    upTypeLabel.ReplaceAll(arrayLabelIds[j],"");
	    downTypeLabel.ReplaceAll(arrayLabelIds[j],"");	    
	  }
	  
	  TString label = upTypeLabel;
	  if(uncIdx!=ENDOFSYSENUM&&(!considerInTotalError(uncIdx)||uncIdx==sysHadDown||uncIdx==sysHadUp)) label= TString("%")+label;
	  (upTypeLabel==downTypeLabel) ? uncIdx+=2 : uncIdx++;
	  	  
	  //std::cout << histoBinIdx << " " << uncIdx << " " << upTypeLabel << " " << downTypeLabel << " " << label << std::endl;

	  TCanvas *canvasUncertaintyDistributions = new TCanvas("canvasUncertaintyDistributions","canvasUncertaintyDistributions",800,600);
	  std::cout << label;
	  if (label.Contains("HadronizationOld")||(calculateError_[xSecVariables_[i]][2*histoBinIdx] && label != "sysNo")){
	    std::cout << " -> considered" << std::endl;
	    TH1F* tempResult = new TH1F(xSecVariables_[i]+"_"+label,xSecVariables_[i]+"_"+label,NBINS,0.5,NBINS+0.5);
	    tempResult->GetXaxis()->SetTitle("Bin Number ("+xSecVariables_[i]+")");
	    tempResult->GetYaxis()->SetTitle(label+" Relative Uncertainty [%]"); 
	    tempResult->SetFillColor(38);
	    tempResult->SetNdivisions(10,"X");
	    
	    int binCounter = 1;
	    
	    std::vector<float> vecValues;

	    for (std::map<unsigned int, TH1F*>::const_iterator histoIter = outerIter->second.begin(); histoIter != outerIter->second.end();){
	      
	      // jump over leading empty bins
	      
		if (histo_[xSecVariables_[i]][sysNo]->GetBinContent(binCounter) == 0)
		  binCounter++;
		
		else{
		  double value = histoIter->second->GetBinContent(histoBinIdx); 
		  tempResult->Fill(binCounter,value);
		  vecValues.push_back(value);
		  
		  binCounter++;
		  histoIter++;
		}
		
		if (binCounter>NBINS) break;  // to avoid endless loop for empty histogrammes
	      }    

	    tempResult->SetMinimum(0);
	    tempResult->SetMaximum(((int)(tempResult->GetMaximum()/5)+1)*5);
	      
	    relUncertDistributions_[xSecVariables_[i]][histoBinIdx] = (TH1F*)tempResult->Clone();
	    canvasUncertaintyDistributions->cd();
	    tempResult->Draw();
	    
	    std::sort(vecValues.begin(),vecValues.end(),SortVectorPrescription<float>);
	    
	    size_t vecSize = vecValues.size();
	    
	    float mean   = (std::accumulate(vecValues.begin(),vecValues.end(),0.0))/vecSize;
	    float median = ( vecValues.size() % 2 != 0 ) ? vecValues[vecSize/2] : (vecValues[vecSize/2-1] + vecValues[vecSize/2]) / 2;

	    // calculate minimum and maximum of medians for uncertainty table (in the plotted variables)
	    if((((xSecVariables_[i].Contains("bq") || xSecVariables_[i].Contains("lep")) && universalplotLabel=="HadronLvPS") ||
	    (   !(xSecVariables_[i].Contains("bq") || xSecVariables_[i].Contains("lep")) && universalplotLabel=="FullPS"    )   ) 
	    &&    xSecVariables_[i].Contains("Norm")){
	      if(!minMedian[label])minMedian[label]=1000.;
	      if(median < minMedian[label])minMedian[label] = median;
	      if(!maxMedian[label])maxMedian[label]=0.;
	      if(median > maxMedian[label])maxMedian[label] = median;
	      aveMedian[label]+=median;
	    }
	    // fill uncertainty table with minimum and maximum (for last plot)
	    if(i==xSecVariables_.size()-1 && universalplotLabel!="PartonLvPS"){
	      TString uncTable2 = label;
	      uncTable2+=fillspaceT(label, 20);
	      uncTable2+=Form(" & %3.1f",aveMedian[label]/NormxSecs);
	      uncTable2+=fillspace(minMedian[label], 3);
	      uncTable2+=Form(" & %3.1f",minMedian[label]); 
	      uncTable2+="\\%"; 
	      uncTable2+=fillspace(maxMedian[label], 3);
	      uncTable2+=Form(" & %3.1f",maxMedian[label]); 
	      uncTable2+="\\% \\\\";
	      uncTable2.ReplaceAll("sys","");
	      if(save) writeToFile(uncTable2, outputFolder+"/uncertaintyDistributionsOverview/uncertaintyTable_"+decayChannel+"_"+universalplotLabel+"_minmax.txt", true);
	    }

	    // prepare line for the complete uncertainty table

	    uncList+=" & ";
	    uncList+=Form("%3.2f",median);

	    DrawLabel((TString)(Form("Minimum: %3.1f",(*vecValues.begin()))) + "%",  
		      gStyle->GetPadLeftMargin() + gStyle->GetTickLength(),
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05,
			gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.25, 
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength(), 12, 0.03);
	    
	    DrawLabel((TString)(Form("Maximum: %3.1f",(*(vecValues.end()-1)))) + "%",
			gStyle->GetPadLeftMargin() + gStyle->GetTickLength(),
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10,
			gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.25, 
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05, 12, 0.03);							
	    
	    DrawLabel((TString)(Form("Mean: %3.1f",mean)) + "%",
		      gStyle->GetPadLeftMargin() + gStyle->GetTickLength(),
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.15,
		      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.25, 
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.10, 12, 0.03);
	    
	    DrawLabel((TString)(Form("Median: %3.1f",median)) + "%",
		      gStyle->GetPadLeftMargin() + gStyle->GetTickLength(),
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.20,
		      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.25, 
		      1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.15, 12, 0.03);
	    
	    vecValues.clear();
	    
	    int initialIgnoreLevel=gErrorIgnoreLevel;
	    if(verbose==0) gErrorIgnoreLevel=kWarning;
	    // a) save to rootFile
	    //if(save) saveToRootFile(outputFile, canvasUncertaintyDistributions, true, verbose, "uncertaintyDistributions");
	      // b) save as eps and png
	    if(save){
	      TString saveName=outputFolder+"/uncertaintyDistributions/relativeUncertainties"+xSecVariables_[i]+"_"+label;
	      if(decayChannel=="combined") saveName+="Combined";
	      saveName+=universalplotLabel;
	      canvasUncertaintyDistributions->Print(saveName+".eps");
	      canvasUncertaintyDistributions->Print(saveName+".png");
	    }
	    gErrorIgnoreLevel=initialIgnoreLevel;
	      
	    delete tempResult; tempResult = NULL;
	  }
	  else{
	    std::cout << std::endl;
	  }
	  delete canvasUncertaintyDistributions; canvasUncertaintyDistributions=NULL;
	}
	// fill complete uncertainty table
	if(save) writeToFile(uncList, outputFolder+"/uncertaintyDistributionsOverview/uncertaintyTable_"+decayChannel+"_"+universalplotLabel+".txt", true);
      }
    }
    std::cout << "uncertainty tables can be found in " << outputFolder+"/uncertaintyDistributionsOverview/uncertaintyTable_"+decayChannel+"_*.txt" << std::endl;
    // close file
    // needed to be able to use the saveToRootFile function, which also opens the file
    file->Close();
    
    // ===============================================================
    //  Part C: draw plot with combined systematic+statistic error
    //          in new cross section plot
    // ===============================================================

    if(save){
      if(verbose>0) std::cout << std::endl << "part C: save xSec plots with stat.+sys. errors and plots with relative uncertainties" << std::endl;
      // loop variables
      for(unsigned int i=0; i<xSecVariables_.size(); ++i){
	// check if total combined errors are available
	if(totalErrors_.count(xSecVariables_[i])>0){
	  if(verbose>0) std::cout << std::endl << "variable: " << xSecVariables_[i] << std::endl;
	  // save rel uncertainties plot for all variables and bins
	  // loop bins
	  for(int bin=0; bin<=histo_[xSecVariables_[i]][sysNo]->GetNbinsX(); ++bin){
	    // make sure that plot exists
	    if(relativeUncertainties_.count(xSecVariables_[i])>0&&relativeUncertainties_[xSecVariables_[i]].count(bin)>0){
	      // create canvas for plot
	      TCanvas* relUnCertaintyCanvas = new TCanvas(relativeUncertainties_[xSecVariables_[i]][bin]->GetName() ,relativeUncertainties_[xSecVariables_[i]][bin]->GetTitle(), 600, 600);
	      //canvasStyle(*relUnCertaintyCanvas);
	      relUnCertaintyCanvas->SetBottomMargin(0.35);
	      relUnCertaintyCanvas->SetLeftMargin(0.20);
	      relUnCertaintyCanvas->SetRightMargin(0.10);
	      relUnCertaintyCanvas->SetTopMargin(0.07);
	      relUnCertaintyCanvas->cd(0);
	      relUnCertaintyCanvas->SetGrid(1,1);
	      relUnCertaintyCanvas->SetTitle(xSecVariables_[i]+"Bin"+getTStringFromInt(bin));
	      // draw plot into canvas (transparent) for having axis AND grid
	      relativeUncertainties_[xSecVariables_[i]][bin]->SetFillColor(10);
	      relativeUncertainties_[xSecVariables_[i]][bin]->SetLineColor(10);
	      relativeUncertainties_[xSecVariables_[i]][bin]->Draw("hist");	      
	      relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->SetTickLength(0.0);
	      // draw axis also on the right side of canvas, uncertainties are only plotted positively after implementing the symmetrization
	      double xPosition = relativeUncertainties_[xSecVariables_[i]][bin]->GetXaxis()->GetXmax();
	      double max   = relativeUncertainties_[xSecVariables_[i]][bin]->GetMaximum();
	      double min       = errMin;
	      TGaxis *axis = new TGaxis(xPosition,min,xPosition,max,min,max,relativeUncertainties_[xSecVariables_[i]][bin]->GetYaxis()->GetNdivisions(),"+L");
	      axis->SetLabelSize(myStyle.GetLabelSize());
	      axis->SetLabelFont(myStyle.GetLabelFont());
	      axis->SetLabelOffset(myStyle.GetLabelOffset());
	      axis->Draw("same");
	      // redraw legend with not considered systematics in grey
	      TH1F* relUnCAxisCopy = (TH1F*)relativeUncertainties_[xSecVariables_[i]][bin]->Clone();
	      relUnCAxisCopy->GetXaxis()->SetLabelColor(kRed);// FIXME
	      // draw every systematic variation with different color
	      TH1F* relUnCertaintyCopy = (TH1F*)relativeUncertainties_[xSecVariables_[i]][bin]->Clone();
	      int colourCounter =    1;
	      int fillStyle     = 1001;
	      unsigned int colour=kBlack;
	      for(int sys=1; sys<=relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX(); ++sys, ++colourCounter){
		colour = (colourCounter%2==0) ? kYellow-7 : 8;
		fillStyle = 1001;
		if (sys<relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX()-2&&considerError_[xSecVariables_[i]][ mapNew2OldSysIndex_[sys] ]==false){
		  fillStyle=3002;
		}
		relUnCertaintyCopy->Reset("ICEMS");
		relUnCertaintyCopy->SetBinContent(sys, relativeUncertainties_[xSecVariables_[i]][bin]->GetBinContent(sys));
		relUnCertaintyCopy->SetFillStyle(fillStyle);
		// single systematic uncertainties	       
		if(sys<=relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX()-3){
		  relUnCertaintyCopy->SetFillColor(colour);
		  relUnCertaintyCopy->DrawCopy("hist same");
		  // delete label entry for considered uncertainties in axis clone
		  if(!(TString(relUnCAxisCopy->GetXaxis()->GetBinLabel(sys)).Contains("("))){
		    //std::cout<< "delete label " << relUnCAxisCopy->GetXaxis()->GetBinLabel(sys) << std::endl;
		    //relUnCAxisCopy->GetXaxis()->SetBinLabel(sys, "");
		  }
		}
		// statistical uncertainty
		else if(sys==relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX()-2){
		  relUnCertaintyCopy->SetFillColor(kOrange+7);
		  relUnCertaintyCopy->DrawCopy("hist same");
		}
		// total systematic uncertainty
		else if(sys==relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX()-1){
		  relUnCertaintyCopy->SetFillColor(kRed+2);
		  relUnCertaintyCopy->DrawCopy("hist same");
		}
		// total uncertainty
		else if(sys==relativeUncertainties_[xSecVariables_[i]][bin]->GetNbinsX()){
		  relUnCertaintyCopy->SetFillColor(38);
		  relUnCertaintyCopy->DrawCopy("hist same");
		  gPad->RedrawAxis("g");
		  relUnCAxisCopy->Draw("AXIS same");
		}
	      }
	      DrawLabel(xSecVariables_[i]+"Bin"+getTStringFromInt(bin), 0.07, 0.90, 0.4, 1.02);
	      // save canvas to file
	      if(save) saveToRootFile( outputFile, relUnCertaintyCanvas, true, verbose, "relativeUncertainties/"+xSecVariables_[i]);
	      // save canvas as eps
	      int initialIgnoreLevel=gErrorIgnoreLevel;
	      if(verbose==0) gErrorIgnoreLevel=kWarning;
	      if(save){
		TString saveNamePDF=outputFolder+"/uncertainties/relativeUncertainties";
		TString saveName=saveNamePDF+relUnCertaintyCanvas->GetTitle();		
		if(decayChannel=="combined"){ saveName+="Combined"; saveNamePDF+="Combined"; }
		saveNamePDF+=LV+PS;
		saveName+=universalplotLabel;
		relUnCertaintyCanvas->Print(saveName+".eps");
		relUnCertaintyCanvas->Print(saveName+".png");
		if(i==0&&bin==1) relUnCertaintyCanvas->Print(saveNamePDF+".pdf(");
		else if (i==xSecVariables_.size()-1&&bin==histo_[xSecVariables_[i]][sysNo]->GetNbinsX()) relUnCertaintyCanvas->Print(saveNamePDF+".pdf)");
		else relUnCertaintyCanvas->Print(saveNamePDF+".pdf");
	      }
	      gErrorIgnoreLevel=initialIgnoreLevel;
	      // delete canvas
	      delete relUnCertaintyCanvas;
	      delete relUnCertaintyCopy;
	      delete axis;
	    }
	  }
	  // open file to get xSec canvas from it
	  file = TFile::Open(outputFile, "UPDATE");
	  // get Canvas
	  TCanvas* canvas = (TCanvas*)(file->Get(xSecFolder+"/sysNo/"+xSecVariables_[i])->Clone());
	  // close file
	  // needed to be able to use the saveToRootFile function 
	  // which also opens the file
	  file->Close();
	  if(canvas){
	    // for inclusive xSecs
	    if(xSecVariables_[i].Contains("inclusive")){
	      TCanvas* canvas2 = new TCanvas(xSecVariables_[i], xSecVariables_[i], 600, 600);
	      //canvasStyle(*canvas2);
	      canvas2->SetLeftMargin(0.05);
	      canvas2->cd();
	      canvas2->SetTitle(xSecVariables_[i]);
	      //canvas2->SetGrid(1,0);
	      // draw xSec 
	      int NxSecBins=3000;
	      double xSecMin=50.;
	      double xSecMax=300.;
	      if(!extrapolate){ 
		xSecMin=0.;
		xSecMax=10.;
	      }	 
	      TH1F* finalInclusiveXSec = new TH1F( "finalInclusive", "finalInclusive", NxSecBins, xSecMin, xSecMax);
	      // x axis as xSec[pb] indicator
	      int xSecBin=roundToInt(histo_["inclusive"][sysNo]->GetBinContent(1)*((double)NxSecBins)/xSecMax)+1;
	      finalInclusiveXSec->SetBinContent(xSecBin, 1);
	      finalInclusiveXSec->SetBinError(xSecBin, 0);
	      histogramStyle(*finalInclusiveXSec, kData , false);
	      TString inclXSeclabel=" #sigma(t#bar{t}#rightarrowX)";
	      if(!extrapolate) inclXSeclabel+="*BR(#mu or e)*A";
	      inclXSeclabel+=" [pb]";
	      axesStyle(*finalInclusiveXSec, inclXSeclabel, " ", 0, 4);
	      finalInclusiveXSec->GetYaxis()->SetNdivisions(0);
	      finalInclusiveXSec->GetXaxis()->SetNdivisions(510);
	      finalInclusiveXSec->GetXaxis()->SetNoExponent(true);
	      finalInclusiveXSec->GetYaxis()->SetTitleSize(0.0);							
	      finalInclusiveXSec->Draw("axis");
	      // draw Theory expectation     
	      double theoryXSecNLL=245.8;
	      double theoryXSecNLO=225.197;
	      double theoryErrorNLOUp  =23.2/157.5*theoryXSecNLO; // FIXME: use relative uncertainty from 7TeV
	      double theoryErrorNLODown=24.4/157.5*theoryXSecNLO; // FIXME: use relative uncertainty from 7TeV
	      double theoryErrorNLLUp  =sqrt(6.2*6.2+6.2*6.2);
	      double theoryErrorNLLDown=sqrt(8.4*8.4+6.4*6.4);
	      TBox* TheoryError = new TBox(theoryXSecNLO-theoryErrorNLODown, 0.0, theoryXSecNLO+theoryErrorNLOUp, 4.0);
	      TBox* TheoryError2= new TBox(theoryXSecNLL-theoryErrorNLLDown, 0.0 ,theoryXSecNLL+theoryErrorNLLUp, 4.0);
	      TheoryError->SetFillColor(kGray);
	      TheoryError2->SetFillColor(kGray+1);
	      if(extrapolate){
		TheoryError->Draw ("same");
		TheoryError2->Draw("same");
	      }
	      // our Analysis result
	      double xSecValue=histo_["inclusive"][sysNo]->GetBinContent(1);
	      double xSecError=histo_["inclusive"][sysNo]->GetBinError(1);
	      double totalErrorUp  = totalErrors_[xSecVariables_[i]]->GetErrorYhigh(1);
	      double totalErrorDown= totalErrors_[xSecVariables_[i]]->GetErrorYlow(1);
	      // 2011 CMS dilepton result - 
	      double cmsLJetsxSecValue=228.4;
	      double cmsLJetsErrorUp  =sqrt(9*9+(29*29+26*26)/2+10.*10.);
	      double cmsLJetsErrorDown=cmsLJetsErrorUp;
	      double cmsLJetsStatError=9;
	      // 2011 CMS 1fb combined result - TOP-11-024 
	      double cmsDiLepxSecValue=227;
	      double cmsDiLepErrorUp  =sqrt(3*3+11*11+10*10);
	      double cmsDiLepErrorDown=cmsDiLepErrorUp;
	      double cmsDiLepStatError=3;
	      if(verbose>0){
		std::cout << "xSec:            " << xSecValue      << std::endl;
		std::cout << "stat. error: +/- " << xSecError      << " (" << xSecError/xSecValue      << ")" << std::endl;
		std::cout << "tot. error up  : " << totalErrorUp   << " (" << totalErrorUp/xSecValue   << ")" << std::endl;	      
		std::cout << "tot. error down: " << totalErrorDown << " (" << totalErrorDown/xSecValue << ")" << std::endl;
	      }
	      totalErrors_[xSecVariables_[i]+"AllInclusive"] = new TGraphAsymmErrors(3);  // number of data point as argument to constructor
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPoint(0, xSecValue,        1);
	      if(extrapolate){
		totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPoint(1, cmsDiLepxSecValue, 2);
		totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPoint(2, cmsLJetsxSecValue, 3);
	      }
	      else{
		totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPoint(1, -1, 2);
		totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPoint(2, -1, 3);
	      }
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPointError(0, totalErrorDown,   totalErrorUp,   0, 0);
	      if(extrapolate){
		totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPointError(1, cmsDiLepErrorDown, cmsDiLepErrorUp, 0, 0);
		totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetPointError(2, cmsLJetsErrorDown, cmsLJetsErrorUp, 0, 0);
	      }
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetLineWidth(3);
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetMarkerSize(1.5);
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetMarkerStyle(20);
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetMarkerColor(kRed+1);
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->SetLineColor(kRed+1);
	      finalInclusiveXSec->Draw("axis same");
	      drawLine(xSecValue-xSecError,               1.0, xSecValue+xSecError,               1.0, kRed+1, 6, 1);
	      if(extrapolate){
		drawLine(cmsDiLepxSecValue-cmsDiLepStatError, 2.0, cmsDiLepxSecValue+cmsDiLepStatError, 2.0, kRed+1, 6, 1);	  
		drawLine(cmsLJetsxSecValue-cmsLJetsStatError, 3.0, cmsLJetsxSecValue+cmsLJetsStatError, 3.0, kRed+1, 6, 1);
	      }
	      totalErrors_[xSecVariables_[i]+"AllInclusive"]->Draw("p same");
	      //Labels
	      if(extrapolate){
		DrawLabel("CMS preliminary 8TeV" , 0.07, 0.70, 0.4, 0.75); 
		DrawLabel("e/#mu + Jets, 2.8/fb"   , 0.07, 0.65, 0.4, 0.70);
		DrawLabel("(CMS-PAS-TOP-12-006)" , 0.07, 0.60, 0.4, 0.65);
		
		DrawLabel("CMS preliminary 8TeV"   , 0.07, 0.50, 0.4, 0.55); 
		DrawLabel("#mu#mu/e#mu/ee, 2.4/fb" , 0.07, 0.45, 0.4, 0.50);
		DrawLabel("(CMS-PAS-TOP-12-007)"   , 0.07, 0.40, 0.4, 0.45);
	      }
	      TString channelLabel="unknown";
	      TString dataLabel=Form(dataSample+" data, %2.1f fb^{-1}",luminosity/1000);
	      if(decayChannel.Contains("mu"  )) channelLabel="#mu + Jets";
	      if(decayChannel.Contains("el"  )) channelLabel="e + Jets";
	      if(decayChannel.Contains("comb")) channelLabel="e/#mu + Jets";
	      if(!extrapolate){
		if(hadron) channelLabel+=" (particle lv PS)";
		else  channelLabel+=" (parton lv PS)";
	      }
	      DrawLabel(channelLabel, 0.07, 0.30, 0.4, 0.35);
	      DrawLabel(dataLabel,    0.07, 0.25, 0.4, 0.30);
	      canvas=(TCanvas*)canvas2->Clone();
	    }
	    // for differential xSecs
	    else{
	      canvas->cd();
	      canvas->SetTitle(xSecVariables_[i]);
	      // unset exp xaxis
	      totalErrors_[xSecVariables_[i]]->GetXaxis()->SetNoExponent(true);
	      double max=totalErrors_[xSecVariables_[i]]->GetMaximum();
	      max*=1.3;
	      if(max>1&&max<100) totalErrors_[xSecVariables_[i]]->GetYaxis()->SetNoExponent(false);
	      else totalErrors_[xSecVariables_[i]]->GetYaxis()->SetNoExponent(true);
	      totalErrors_[xSecVariables_[i]]->SetName("dataTotalError");
	      // for combined cross sections:
	      // change TH1F with statistical errors into 
	      // TGraphAsymmErrors with bin center corrections
	      if(decayChannel=="combined"&&xSecVariables_[i]!="inclusive"){
		TString plotName = xSecVariables_[i];
		std::cout << "shifted stat value for: variable " << plotName << std::endl;
		plotName.ReplaceAll("Norm","");
		// get data plot
		TH1F* dataStat = (TH1F*)canvas->GetPrimitive(plotName+"kData");
		// convert to TGraphAsymmErrors
		TGraphAsymmErrors* statErrors= new TGraphAsymmErrors(dataStat->GetNbinsX());
		statErrors->SetName("dataStatError");
		statErrors->SetLineWidth(3.5) ; // totalErrors_[xSecVariables_[i]]->GetLineWidth());
		statErrors->SetLineColor(totalErrors_[xSecVariables_[i]]->GetLineColor());
		for(int bin=1; bin<=dataStat->GetNbinsX(); ++bin){
		  // exclude over/underflow bins with empty bin content or zero width 
		  if(dataStat->GetBinWidth(bin)!=0&&dataStat->GetBinContent(bin)!=0){
		    if(useBCC && verbose>1){ 
		      std::cout << "try to access BCC result map (2nd) for " << plotName << " (bin " << bin << ")" << std::endl;
		      for(std::map<TString, std::vector<double> >::const_iterator val=correctedCenters_.begin(); val!=correctedCenters_.end(); ++val){
			std::cout << "found entry " << val->first << " (size " << val->second.size() << ", ";
			std::cout << "center (bin " << bin << ") " << ((int)val->second.size() >= bin ? val->second.at(bin-1) : -1.) << ")" << std::endl;
		      }
		    }
		    double pointXValue = (useBCC&&correctedCenters_.find(plotName)!=correctedCenters_.end()) ? correctedCenters_[plotName].at(bin-1): histo_[xSecVariables_[i]][sysNo]->GetBinCenter(bin);
		    double pointXError = (useBCC&&corrCenterErrors_.find(plotName)!=corrCenterErrors_.end()) ? corrCenterErrors_[plotName].at(bin-1): 0;
		    double pointYError = dataStat->GetBinError(bin);
		    statErrors->SetPoint(bin, pointXValue, dataStat->GetBinContent(bin));
		    statErrors->SetPoint(0, 0, -1000);
		    pointXError=0;
		    statErrors->SetPointError(bin, pointXError, pointXError, pointYError, pointYError);
		  }
		}
		whipEmptyBinsAway(statErrors, verbose);
		// print new object
		statErrors->SetMarkerSize(0.01);
		statErrors->Draw("e p same");
		// delete (old) TH1F plot from canvas
		canvas->GetListOfPrimitives()->Remove(canvas->GetPrimitive(plotName+"kData"));
	      }
	      // Draw errors into Canvas
	      totalErrors_[xSecVariables_[i]]->Draw("][ p z same");
	      canvas->SetName(xSecVariables_[i]);
	    }
	    // save Canvas
	    int initialIgnoreLevel=gErrorIgnoreLevel;
	    if(verbose==0) gErrorIgnoreLevel=kWarning;
	    if(save){
	      // a) within rootFile
	      saveToRootFile(outputFile, canvas, true, verbose, "finalXSec");
	      // b) as eps, png and pdf
	      TString saveNamePDF=outputFolder+"/xSec/finalXSec";
	      TString saveName=saveNamePDF+xSecVariables_[i]; 
	      if(decayChannel=="combined"){ saveName+="Combined"; saveNamePDF+="Combined"; }
	      saveName+=universalplotLabel+closureLabel;
	      saveNamePDF+=LV+PS;
	      canvas->Print(saveName+".eps");
	      canvas->Print(saveName+".png");
	      if(i==0) canvas->Print(saveNamePDF+".pdf(");
	      else if (i==xSecVariables_.size()-1) canvas->Print(saveNamePDF+".pdf)");
	      else canvas->Print(saveNamePDF+".pdf");
	      gErrorIgnoreLevel=initialIgnoreLevel;
	      //delete statErrors;
	    }
	  }
	}
      }
    }
    // close file
    //hadfile->Close();
  }

  // close file
  file->Close();

}
