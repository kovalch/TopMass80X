// --------------------------------------------------------------
// This macro calculates the purity, stability and efficiency
// for a user-defined fixed binning.
// (Modified version of the macro findBinning.C)
// --------------------------------------------------------------

#include <sstream>
#include <string>
#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPaletteAxis.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include "TTree.h"
#include "algorithm"
#include <TMath.h>
#include "HHStyle.h"
#include "basicFunctions.h"
#include "../../../../TopUtils/interface/extract_sigma.h"

void purityStabilityEfficiency(TString variable = "topPt", bool save=true, TString lepton="combined", 
			       TString inputFolderName=AnalysisFolder, bool plotAcceptance = false, 
			       bool plotEfficiencyPhaseSpace = false, bool plotEfficiency2 = false, double chi2Max=7.824,//9999.,//7.824,
			       int verbose=1, bool hadron=false, int qAssignment=-1,
			       bool fitGaussRes=true, bool printSeparateRes = false)
{
  // ARGUMENTS of function:
  // variable:       choose variable to plot, e.g.:
  //                 topPt, topY, ttbarPt, ttbarY, ttbarMass, lepPt, lepY
  // save:           save the purity and stability canvas 
  // lepton:         "muon" or "elec" or "combined"
  // plotAcceptance: in addition to purity and stability also acceptance is plotted (if true)
  // plotEfficiencyPhaseSpace: in addition to purity, stability, acceptance also 
  //                           efficiency in restricted phase space (i.e. Acceptance=1) is plotted (if true)
  // hadron:         draw plots for hadron level
  // qAssignment     make a cut on kinFit jet permutation: -1= no cut, 0=ok, 1=bb, 2=blepq, 3=bhadq, 4=bbqlep, 5=bbqhad, 6=bbqq, 7=jmis, 8=wrongj, 9=nomatch
  // fitGaussRes:    fit a Gauss to the residuum plots to obtain the resolution
  // printSeparateRes save the residuum plots separately

  TString path=groupSpace; // folder with all rootfiles
  bool useTree=true; // use default 2D histo or create 2D histo from tree, allows chi2 cuts
  if(variable.Contains("lbMass")) useTree=false; //FIXME
  TString nameExt="";
  double chi2MaxInit=chi2Max;
  // if chi2 sel is applied-plot chi2 efficiency in addition?
  bool drawChi2Eff=false;
  // if chi2 sel is applied-plot pur and stab without chi2 sel in addition?
  bool plotDefToCompare=false;
  if(chi2Max>100) plotDefToCompare=false;
  if(!useTree){//||variable.Contains("rhos")||variable.Contains("Njets")){
    // use chi2/probability selected events
    if(chi2Max<100) nameExt+="ProbSel";
    chi2Max=99999; // can be done only with tree
  }
  if(chi2Max>100) drawChi2Eff=false;
  // legend drawing: drawSepLegend=false -> only separate canvas with legend is drawn
  bool drawSepLegend=false;
  // legend drawing: drawAlwaysLegenddrawAlwaysLegend=false->only for topPt
  bool drawAlwaysLegend=true;
  if(!drawSepLegend) drawAlwaysLegend=false;
  if(lepton.Contains("combined")) lepton="combined";
  // output folder in case of saving the canvases:
  TString outputFolder = "./diffXSecFromSignal/plots/"+lepton+"/2012/binning";
  if(useTree&&chi2Max<100){ 
    plotEfficiencyPhaseSpace = false;
    plotAcceptance = false;
    plotEfficiency2= false;
  }
  // create resolution gen vs. rec plots for choseb binning
  bool createResolutionPlots=true;
  if(!useTree) createResolutionPlots=false;
  
  // shall b-jets be swapped if deltaR is better then?
  bool swapBb = false;
  if(hadron) swapBb=true;
  
  int initialIgnoreLevel=gErrorIgnoreLevel;
  if(verbose<=1) gErrorIgnoreLevel=kWarning;
  
  if(qAssignment>=0) std::cout << "qAssignment=" << qAssignment << std::endl;
  // input file
  TFile* myFile1 = new TFile(path+""+inputFolderName+"/"+TopFilename(kSig, sysNo, (std::string)lepton), "READ");
  TFile* myFile2 = new TFile(path+""+inputFolderName+"/"+TopFilename(kSig, sysNo, (std::string)lepton), "READ");
  
  if(!myFile1){
    std::cout << "ERROR when opening file " << path+""+inputFolderName+"/"+TopFilename(kSig, sysNo, (std::string)lepton) << std::endl;
    exit(0);
  }

  if(lepton=="electron") lepton="elec";
  
  // choose phase space
  // b) for restricted phase space:
  // b1) parton PS: hadron=false
  // b2) hadron PS: hadron=true
  TString LV="Parton";
  if(hadron) LV="Hadron";
  if(verbose>1){
    std::cout << LV << " level Phase Space will be used!" << std::endl; 
  }
  TString universalplotLabel="";
  universalplotLabel=LV+"LvPS";

  // create more control output for debugging?
  bool debug=false;
  
  //choose proper root folder name with histos and tree
  //default for parton level
  TString folderRecoKin          ="analyzeTopRecoKinematicsKinFit";
  TString folderRecoKinGeneral   ="analyzeTopRecoKinematicsKinFit";
  if(nameExt!=""){
    folderRecoKin+=nameExt;
    folderRecoKinGeneral+=nameExt;
  }
  TString folderGenKin           ="analyzeTopPartonLevelKinematics";
  TString folderGenKinPhaseSpace ="analyzeTopPartonLevelKinematicsPhaseSpace";

  // label if variable is gen or reco in folderRecoKin
  TString genExtTree="PartonTruth";
  TString recExtTree="";
  TString genExtHisto="";
  TString recExtHisto="";
  // text draw option for 2d response matrix
  TString textOpt="TEXT";
  
  // if hadron level
  if(hadron){
    if(variable.Contains("bq")||variable.Contains("bbbar")||variable.Contains("lbMass")){
      folderRecoKin          = "analyzeTopRecoKinematicsBjets";
      folderGenKinPhaseSpace = "analyzeTopHadronLevelKinematicsBjetsPhaseSpace";
      folderGenKin           = "analyzeTopHadronLevelKinematicsBjetsPhaseSpace";
      genExtTree="Gen";
      recExtTree="Rec";
      genExtHisto="Gen";
      recExtHisto="Rec";
    }
    else if(variable.Contains("lep")){
      folderRecoKin          = "analyzeTopRecoKinematicsLepton";
      folderGenKinPhaseSpace = "analyzeTopHadronLevelKinematicsLeptonPhaseSpace";
      folderGenKin           = "analyzeTopHadronLevelKinematicsLeptonPhaseSpace";
      genExtTree="Gen";
      recExtTree="Rec";
      genExtHisto="Gen";
      recExtHisto="Rec";
    }
    else if(variable.Contains("Njets")||variable.Contains("rhos")){
      folderRecoKin          = "compositedKinematicsKinFit";
      folderRecoKinGeneral   = "compositedKinematicsKinFit";
      folderGenKinPhaseSpace = "compositedHadronGenPhaseSpace";
      folderGenKin           = "compositedHadronGenPhaseSpace";
      genExtTree="True";
    }
    else{
      std::cout << "--------------------------------------------------------------!" <<std::endl;
      std::cout << "No hadron level plots calculated!!! Only for b-jet and lepton based quantities!" <<std::endl;
      std::cout << "Exit programme!" <<std::endl;
      return;
    }
    plotAcceptance=false;
  }
  // take care of name extension, e.g. due to probability selection
  TString folderRecoKinAlter=folderRecoKin;
  if(nameExt!=""){
    folderRecoKin+=nameExt;
    folderRecoKinGeneral+=nameExt;
    //if(folderRecoKin.Contains(       "composited")) folderRecoKin.ReplaceAll(       "KinFitProbSel", "ProbSel");
    //if(folderRecoKinGeneral.Contains("composited")) folderRecoKinGeneral.ReplaceAll("KinFitProbSel", "ProbSel");
  }

  
  ///---------------------------------------------------------------------------------------------
  // determine here the binning:
  // change x axis range and determine edges of residual histo for resolution studies
  double rangeUserLeft = -1e6;
  double rangeUserRight = 1e6;
  double resEdgeL = -1e6;
  double resEdgeR = 1e6;
  double relResEdgeL = -1e6;
  double relResEdgeR = 1e6;
  if ( variable == "topPt")          { rangeUserLeft = 0         ; rangeUserRight = 500-0.001;  resEdgeL = -500.0;  resEdgeR=400.0;  relResEdgeL = -1.0; relResEdgeR=4.0; }
  else if ( variable == "topY")      { rangeUserLeft = -2.5+0.001; rangeUserRight = 2.5-0.001;  resEdgeL =   -3.0;  resEdgeR=  3.0;  relResEdgeL = -1.5; relResEdgeR=1.5; }
  else if ( variable == "ttbarPt")   { rangeUserLeft = 0         ; rangeUserRight = 300-0.001;  resEdgeL = -300.0;  resEdgeR=300.0;  relResEdgeL = -1.0; relResEdgeR=6.0; }
  else if ( variable == "ttbarY")    { rangeUserLeft = -2.5+0.001; rangeUserRight = 2.5-0.001;  resEdgeL =   -2.0;  resEdgeR=  2.0;  relResEdgeL = -2.0; relResEdgeR=2.0; }
  else if ( variable == "ttbarMass") { rangeUserLeft = 345+0.001 ; rangeUserRight = 1600-0.001; resEdgeL = -400.0;  resEdgeR=400.0;  relResEdgeL = -1.0; relResEdgeR=4.0; }
  else if ( variable == "lepPt" )    { rangeUserLeft = 30        ; rangeUserRight = 400-0.001;  resEdgeL =  -35.0;  resEdgeR= 35.0;  relResEdgeL = -0.2; relResEdgeR=0.2; /*textOpt="";*/}
  else if ( variable == "lepEta")    { rangeUserLeft = -2.1+0.001; rangeUserRight = 2.1-0.001;  resEdgeL =  -0.005; resEdgeR= 0.005; relResEdgeL = -0.5; relResEdgeR=0.5; /*textOpt="";*/}
  else if ( variable == "bqPt")      { rangeUserLeft = 30+0.001  ; rangeUserRight = 350-0.001;  resEdgeL = -200.0;  resEdgeR=200.0;  relResEdgeL = -1.0; relResEdgeR=2.0 ;}
  else if ( variable == "bqEta")     { rangeUserLeft = -2.4+0.001; rangeUserRight = 2.4-0.001;  resEdgeL =   -0.2;  resEdgeR=  0.2;  relResEdgeL = -1.0; relResEdgeR=1.0; }
  else if ( variable == "bbbarPt")   { rangeUserLeft = 0         ; rangeUserRight =1000-0.001;  resEdgeL = -400.0;  resEdgeR=  400;  relResEdgeL = -1.5; relResEdgeR=4.0; }
  else if ( variable == "bbbarMass") { rangeUserLeft = 0         ; rangeUserRight =1500-0.001;  resEdgeL = -400.0;  resEdgeR=  400;  relResEdgeL = -1.5; relResEdgeR=4.0; }
  else if ( variable == "Njets"    ) { rangeUserLeft = -0.5      ; rangeUserRight =14.5      ;  resEdgeL =  0.   ;  resEdgeR=  5.0;  relResEdgeL = -0.5; relResEdgeR=0.5; }
  else if ( variable == "rhos"     ) { rangeUserLeft = 0         ; rangeUserRight =1.1-0.001 ;  resEdgeL = -2.0  ;  resEdgeR=  2.0;  relResEdgeL = -1.5; relResEdgeR=1.5; }


  //else if ( variable == "lepPt"  && lepton == "muon")    { rangeUserLeft = 20  ; rangeUserRight = 400-0.001; }
  //else if ( variable == "lepPt"  && lepton == "elec")    { rangeUserLeft = 30  ; rangeUserRight = 400-0.001; }
  //else if ( variable == "lepEta" && lepton == "muon")    { rangeUserLeft = -2.1; rangeUserRight = 2.1-0.001; }
  //else if ( variable == "lepEta" && lepton == "elec")    { rangeUserLeft = -2.5; rangeUserRight = 2.5-0.001; }
  
  // determine here the binning:
  // Attention: binning ALWAYS should start from left boundary of first and bin of the input file;
  // it ends with the last bin of the input file even if this is not part of the array defined below
  
  std::map<TString, std::vector<double> > binningMap = makeVariableBinning();
  std::vector<double> xBins = binningMap[variable];
  int NxBins = xBins.size()-1;
  
  // pt
  // double xBins[] = {0, 60, 120, 200, 280, 400};
  // y
  // double xBins[] = {-5., -2.5, -1.5, -1., -0.5, 0, 0.5, 1., 1.5, 2.5};
  // ttbar pt
  // double xBins[] = {0, 20, 60, 110, 200};
  // ttbar Y
  // double xBins[] = {-5., -1.3, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.3};
  // ttbarMass
  // double xBins[] = {0, 400, 500, 630, 900};
  // double xBins[] = {0, 345, 410, 480, 580, 750};
  // Korea:
  // double xBins[] = {0, 350, 400, 450, 500, 550,600,700,800};
  // lepPt
  // double xBins[] = {0., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100., 120., 150., 200., 275., 400., 1200.};
  // lepEta
  // double xBins[] = {-5., -2.5, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5, 5.};
  // b-quark Pt
  // double bqPtBins[]  = {0, 30, 60, 95, 140, 200, 400 };
  // b-quark Eta 
  // double bqEtaBins[] = {-5.0, -2.5, -1.5, -1.0, -0.5, 0., 0.5, 1.0, 1.5, 2.5, 5.0};
  //int NxBins = sizeof(xBins)/sizeof(double);

  if(verbose>1) std::cout << "size of all calculated bins, xBins: " << NxBins << std::endl;
  
  ///---------------------------------------------------------------------------------------------
  // histos for resolution studies for each bin of the variable of interest:
  // vector with one histo for each bin of the distribution of the variable of interest
  std::vector<TH1F*> residualHistos_;
  std::vector<TH1F*> relativeResidualHistos_;
  for(int iBin=0; iBin < NxBins; iBin++){
    TString histoName  = "residualHisto"+variable+universalplotLabel+"Bin"+Form("%i",iBin);
    TString histoTitle = "residualHisto"+universalplotLabel+" "+Form("%4.1f",xBins[iBin])+"<"+variable+"<"+Form("%4.1f",xBins[iBin+1]);
    residualHistos_.push_back        ( new TH1F (histoName, histoTitle, 100, resEdgeL, resEdgeR) );
    relativeResidualHistos_.push_back( new TH1F ("relative"+histoName, "relative"+histoTitle, 100, relResEdgeL, relResEdgeR) );
  }
    
  // make a nice style
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.cd();
  TGaxis::SetMaxDigits(3);
  gROOT->SetStyle("HHStyle");
  gROOT->ForceStyle();	

  // ---
  //    get the 2D histogram
  // ---
  // x-Axis should be generated value, y-Axis should be reconstructed value
  TH2F* myHist2d = 0;
  // empty histogram
  if(useTree) myHist2d = new TH2F("correlationMatrix", "correlationMatrix", NxBins, &xBins[0], NxBins, &xBins[0]);
  else{
    TH2* targetPlot;
    TString targetname=folderRecoKin+"/"+variable+"_";
    myFile1->GetObject(targetname, targetPlot);
    if(targetPlot) myHist2d = (TH2F*)(targetPlot->Clone("correlationMatrix"));
    else{
      std::cout << "ERROR: can not open plot " << targetname << " from file" << myFile1->GetName() << std::endl; 
      exit(0);
    }
    // apply rebinning wrt chosen variable test binning
    myHist2d=reBinTH2F(*myHist2d, xBins, verbose);
  }
  // control histo to compare e.g. results for chi2 selection wrt. default
  TH2F* myHist2d2 = (TH2F*)myHist2d->Clone();
  // create response matrix with and without efficiency correction in final binning; will be only filled if useTree==true
  TH2F* responseMatrix  = new TH2F("responseMatrix", "responseMatrix", NxBins, &xBins[0], NxBins, &xBins[0]);
  TH2F* migrationMatrix = new TH2F("migrationMatrix", "migrationMatrix", NxBins, &xBins[0], NxBins, &xBins[0]);
  // create pointer to avoid problems with object persistence  
  TH1F* chi2eff = NULL;
  TH1F* all     = NULL;
  // loading tree
  if (useTree){
    std::cout << "file " << myFile1->GetName() << ", folder " << folderRecoKin << "/tree" << std::endl;
    TTree *tree=(TTree*)(myFile1->Get(folderRecoKin+"/tree"));
    // for hadron level, the specific RecoKin folder does not contain a complete tree with all variables;
    //  add the general tree
    TTree *treeGeneral=(TTree*)(myFile1->Get(folderRecoKinGeneral+"/tree"));
    if(!tree||tree->IsZombie()){
      std::cout << "there seems to be a problem with the chosen tree " << folderRecoKin<<"/tree" << std::endl;
      exit(0);  
    }
    if(!treeGeneral||treeGeneral->IsZombie()){
      std::cout << "there seems to be a problem with the chosen tree " << folderRecoKinGeneral<<"/tree" << std::endl;
      exit(0);  
    }
    if(tree->GetEntries()!=treeGeneral->GetEntries()){
      std::cout << "ERROR!!! Entries of tree (" << tree->GetEntries() << ") != entries of treeGeneral (" << treeGeneral->GetEntries() <<")" <<std::endl;
      exit(0);
    }
    // list relevant tree entries
    std::vector<TString> variable_;
    if(variable.Contains("top")&&!variable.Contains("TtbarSys")){
      variable_.push_back(variable+"Lep");
      variable_.push_back(variable+"Had");
    }
    else if(variable.Contains("bqPt")){
      if(!hadron){
	variable_.push_back("bqPtLep");
	variable_.push_back("bqPtHad");
      }
      else{
	variable_.push_back("bqPt");
	variable_.push_back("bbarqPt");	
      }
    }  
    else if(variable.Contains("bqEta")){
      if(!hadron){
	variable_.push_back("bqEtaLep");
	variable_.push_back("bqEtaHad");
      }
      else{
	variable_.push_back("bqEta");
	variable_.push_back("bbarqEta");	
      }
    }
    else variable_.push_back(variable);
    TString genHistoName=variable+genExtHisto;
    // TString genHistoName2=genHistoName;
    if(genHistoName.Contains("rhos" )) genHistoName="rhosGen";
    if(genHistoName.Contains("Njets")) genHistoName="Ngenjets";
    // container for values read from tree
    std::map< TString, float > value_;
    //std::map< TString, double > value_;
    // container for chi2 cut efficiency
    chi2eff = (TH1F*)(((TH1F*)myFile1->Get(folderGenKinPhaseSpace+"/"+genHistoName))->Clone());
    chi2eff->Scale(0);
    all=(TH1F*)chi2eff->Clone();
  
    // initialize map entries with 0 
    // to avoid problems with the map re-ordering
    // when new entries are added
    value_["weight"]=1;
    value_["chi2"]  =0;
    bool bbSwapBetter  = false;
    value_["qAssignment"]  =-2;
    for(unsigned int i=0; i<variable_.size(); ++i){
      value_[variable_[i]+recExtTree]=0;
      value_[variable_[i]+genExtTree]=0;
    }
    // initialize branches
    tree->SetBranchStatus("*",0);
    treeGeneral->SetBranchStatus("*",0);
    tree->SetBranchStatus("weight",1);
    tree->SetBranchAddress("weight",(&value_["weight"]));
    treeGeneral->SetBranchStatus("qAssignment",1);
    if(nameExt!="ProbSel") treeGeneral->SetBranchStatus("chi2",1);
    if(nameExt!="ProbSel") treeGeneral->SetBranchAddress("chi2",(&value_["chi2"]));
    treeGeneral->SetBranchAddress("qAssignment",(&value_["qAssignment"]));
    if(variable.Contains("bq") && swapBb){
      tree->SetBranchStatus("bbSwapBetter",1);
      tree->SetBranchAddress("bbSwapBetter",(&bbSwapBetter));
    }
      
    for(unsigned int i=0; i<variable_.size(); ++i){
      if(verbose>1) std::cout << "searching branch " << variable_[i]+recExtTree << " in " << tree->GetName() << std::endl;
      // activate branches
      tree->SetBranchStatus(variable_[i]+recExtTree,1);
      tree->SetBranchStatus(variable_[i]+genExtTree,1);
      // save branch values in map
      tree->SetBranchAddress(variable_[i]+recExtTree,(&value_[variable_[i]+recExtTree]));
      tree->SetBranchAddress(variable_[i]+genExtTree,(&value_[variable_[i]+genExtTree]));
    }
    // loop all events to fill plots
    for(unsigned int event=0; event<tree->GetEntries(); ++event){
      if(debug) std::cout << "event " << event+1 << "/" << tree->GetEntries() << std::endl;
      tree->GetEntry(event);
      treeGeneral->GetEntry(event);
      double weight=value_["weight"];
      double chi2= nameExt=="ProbSel" ? 0. : value_["chi2"];
      for(unsigned int i=0; i<variable_.size(); ++i){
	double rec =value_[variable_[i]+recExtTree];
	TString genName=variable_[i]+genExtTree;
	//if(genName.Contains("rhos" )||genName.Contains("Njets")) weight=1.; // FIXME MARTIN: weight is wrong at the moment for composited... analyzer
	if(genName.Contains("rhos" )) genName="rhosTrue";
	if(genName.Contains("Njets")) genName="NjetsTrue";	
	double gen =value_[genName];
	// swap b and bbar if deltaR is better
	if(genName.Contains("bq" )&&swapBb && bbSwapBetter){
	  gen = (i==0) ? value_[variable_[1]+genExtTree] : value_[variable_[0]+genExtTree];
	}
	if(rec==-9999||gen==-9999){ 
	  std::cout << "variable " << variable << " is not filled properly:" << std::endl;
	  std::cout << "rec " <<  variable_[i]+recExtTree << ": " << rec << std::endl;
	  std::cout << "gen " <<  genName                 << ": " << gen << std::endl;
	  exit(0);
	}
	else if(debug){
	  std::cout << variable << ":" << std::endl;
	  std::cout << "rec " <<  variable_[i]+recExtTree << ": " << rec << std::endl;
	  std::cout << "gen " <<  genName                 << ": " << gen << std::endl;
	  std::cout << "chi2: " <<   chi2 << std::endl;
	  std::cout << "qAssignment: " << value_["qAssignment"] << std::endl;
	  std::cout << "weight: " << weight << std::endl;
	}
	// check that reco value is valid
	// (needed as eg not all events in tree have an additional jet and therefore a valid rhos value)
	if(rec>-1000){
	  if(debug) std::cout << "filling: ";
	  // apply chi2 and qAssignment cut
	  if(chi2<chi2Max && (qAssignment==-1 || qAssignment==value_["qAssignment"])){
	    if(debug) std::cout << "responseMatrix myHist2d ";
	    myHist2d       ->Fill(gen, rec, weight);
	    responseMatrix ->Fill(gen, rec, weight);
	    // fill events passing chi2
	    if(nameExt!="ProbSel"){
	      if(debug) std::cout << "chi2eff ";
	      chi2eff->Fill(rec, weight);
	    }
	    // fill residual histo of the right bin
	    if(debug) std::cout << "residuals ";
	    for(int iBin=0; iBin < NxBins; iBin++){
	      if(gen>xBins[iBin] && gen< xBins[iBin+1]){
		residualHistos_.at(iBin)->Fill(rec-gen, weight);
		if(gen!=0) relativeResidualHistos_.at(iBin)->Fill((rec-gen)/gen, weight);
		else       relativeResidualHistos_.at(iBin)->Fill(-1);
		break;
	      }
	    }
	  }
	  if((qAssignment==-1 || qAssignment==value_["qAssignment"])){
	    myHist2d2->Fill(gen, rec, weight);
	    if(debug) std::cout << "myHist2d2 ";
	  }
	  // fill all events
	  all->Fill(rec, weight);
	  if(debug) std::cout << "all " << std::endl;
	}
      }
    }
  }
  else{
    // use loaded plot
    responseMatrix=(TH2F*)myHist2d->Clone();
    TH2* targetPlot;
    TString targetname=folderRecoKinAlter+"/"+variable+"_";
    myFile1->GetObject(targetname, targetPlot);
    if(targetPlot) myHist2d2 = (TH2F*)(targetPlot->Clone("correlationMatrix2"));
    else{
      std::cout << "ERROR: can not open plot " << targetname << " from file" << myFile1->GetName() << std::endl;
      exit(0);
    }
    // apply rebinning wrt chosen variable test binning
    myHist2d2=reBinTH2F(*myHist2d2, xBins, verbose);
  }
  int xbins = myHist2d->GetNbinsX();
  double xmax = myHist2d->GetXaxis()->GetXmax();
  //double xmin = myHist2d->GetXaxis()->GetXmin();
  // draw the obtained histo
  TCanvas* Canv1 = new TCanvas("analyzeTopRecoKinematicsKinFitCorr","analyzeTopRecoKinematicsKinFitCorr",600,600);
  Canv1->cd();
  Canv1->SetRightMargin(0.1);
  myHist2d->SetStats(kFALSE);
  myHist2d->GetXaxis()->SetTitle("gen");
  myHist2d->GetYaxis()->SetTitle("rec");
  //myHist2d->Draw("colz");
  myHist2d->Draw("");
  
  // initialize variables
  int numberOfBins=0;
  vector<double> binvec;

  vector<double> purityvec;
  double purity=0.;
  double puritydenom=0.;

  vector<double> stabilityvec;
  double stability=0.;
  double stabilitydenom=0.;

  double myenum=0.;

  vector<double> purityvec2;
  double purity2=0.;
  double puritydenom2=0.;

  vector<double> stabilityvec2;
  double stability2=0.;
  double stabilitydenom2=0.;

  double myenum2=0.;

  int iSave=0;
  
  // calculate purity and stability
  // loop over bins of one axis of 2D histogram
  // (need to be the same binning for x and y anyway!)
  for(int i=1; i<=xbins; i++)
    {
      // save bin number which will be the lower edge of the final bin
      if(!iSave)iSave=i;
      // loop over bins of the other axis of 2D histogram
      for(int j=0; j<=xbins+1; j++)
	{
	  // denominator for purity (number of reconstructed events)
	  puritydenom +=myHist2d ->GetBinContent(j,i);
	  puritydenom2+=myHist2d2->GetBinContent(j,i);
	  // denominator for stability (number of generated events)
	  stabilitydenom +=myHist2d ->GetBinContent(i,j);
	  stabilitydenom2+=myHist2d2->GetBinContent(i,j);
	}
      // loop over bins starting from lower edge of the final bin
      for(int k=iSave; k<=i; k++)
	{
	  // loop over bins starting from lower edge of the final bin
	  for(int l=iSave; l<=i; l++)
	    {
	      // enumerator
	      myenum +=myHist2d ->GetBinContent(l,k);
	      myenum2+=myHist2d2->GetBinContent(l,k);
	    }
	}
      // calculate purity (avoid dividing by zero)
      if(puritydenom )purity =myenum /puritydenom;
      else purity=0;
      if(puritydenom2)purity2=myenum2/puritydenom2;
      else purity2=0;
      // calculate stability (avoid dividing by zero)
      if(stabilitydenom )stability =myenum /stabilitydenom;
      else stability=0;
      if(stabilitydenom2)stability2=myenum2/stabilitydenom2;
      else stability2=0;
      if( (numberOfBins+1 < NxBins+1 && myHist2d->GetBinCenter(i+1) > xBins[numberOfBins+1]) || i==xbins){
	if(debug) std::cout << "(stability/stability2)=" << stability << " / " << stability2 << std::endl;
	if(debug) std::cout << "(purity/purity2)=" << purity << " / " << purity2 << std::endl;
	// write results to vector
	purityvec.push_back(purity);
	stabilityvec.push_back(stability);
	binvec.push_back(myHist2d->GetBinLowEdge(iSave));
	purityvec2.push_back(purity2);
	stabilityvec2.push_back(stability2);
	// count bins
	numberOfBins++;
	// reset variables
	puritydenom=0.;
	stabilitydenom=0.;
	puritydenom2=0.;
	stabilitydenom2=0.;
	iSave=0;
      }
      myenum=0.;
      myenum2=0.;
    }
  binvec.push_back(xmax);

  // print the binning
  if(verbose>1){
    cout<<numberOfBins<<" bins: {";
    for(int i=0; i<numberOfBins; i++)
      {
	cout<<binvec[i]<<", ";
      }
    cout<<xmax<<"}"<<endl;
  }
  
  // calculate stability & purity histogram with the determined binning
  TH1F* purityhist    = new TH1F("purityhist"   ,"purityhist"   ,numberOfBins, &(binvec.front()));
  TH1F* stabilityhist = new TH1F("stabilityhist","stabilityhist",numberOfBins, &(binvec.front()));
  TH1F* purityhist2    = new TH1F("purityhist2"   ,"purityhist2"   ,numberOfBins, &(binvec.front()));
  TH1F* stabilityhist2 = new TH1F("stabilityhist2","stabilityhist2",numberOfBins, &(binvec.front()));
  
  for(int i=0; i<numberOfBins; i++)
    {
      purityhist   ->SetBinContent (i+1,purityvec[i]    );
      stabilityhist->SetBinContent (i+1,stabilityvec[i] );
      purityhist2   ->SetBinContent(i+1,purityvec2[i]   );
      stabilityhist2->SetBinContent(i+1,stabilityvec2[i]);
    }
  
  //------------------------------------------------------------------------------------
  // ----------- Bin-by-bin (BBB) Efficiency calculation--------------------------------

  // ----------- BBB Efficiency*Acceptance--------------------------------
  TString genHistoName=variable+genExtHisto;
  if(genHistoName.Contains("rhos" )) genHistoName="rhosGen";
  if(genHistoName.Contains("Njets")) genHistoName="Ngenjets";
  // get histogram of generated quantity
  if(debug) std::cout << "genHist: "<< folderGenKin+"/"+genHistoName << std::endl;
  TH1* targetPlotgenHist;
  TString targetnamegenHist=folderGenKin+"/"+genHistoName;
  myFile1->GetObject(targetnamegenHist, targetPlotgenHist);
  if(!targetPlotgenHist){
    std::cout << "ERROR: can not open plot " << targetnamegenHist << " from file" << myFile1->GetName() << std::endl; 
    exit(0);
  }
  TH1F* genHist = (TH1F*)(targetPlotgenHist->Clone());
  //genHist->SetStats(kFALSE);
  TCanvas* Canv2 = new TCanvas(folderGenKin+"2",folderGenKin+"2",600,600);
  Canv2->cd();
  genHist->Draw();
  // rebin histogram of generated quantity
  genHist = (TH1F*)genHist->Rebin(binvec.size()-1, genHist->GetName(), &(binvec.front()));
  // get histogram of reconstructed quantity
  TH1F* effHistBBB = (TH1F*)(((TH1F*)myFile1->Get(folderRecoKin+"/"+variable+recExtHisto))->Clone());
  //effHistBBB->SetStats(kFALSE);
  TCanvas* Canv3 = new TCanvas(folderRecoKin,folderRecoKin,600,600);
  Canv3->cd();
  effHistBBB->Draw();
  // rebin histogram of reconstructed quantity
  effHistBBB = (TH1F*)effHistBBB->Rebin(binvec.size()-1, effHistBBB->GetName(), &(binvec.front()));
  TH1F* effHistBBBPS=(TH1F*)(effHistBBB->Clone());
  // calculate efficiency histogram
  effHistBBB->Divide(genHist);
  
  // ----------- BBB Efficiency (i.e. in in limited phase) --------------------------------
  TH1F* genHistPS = 0;
  TH1F* accHistBBB=0;
  // get histogram of generated quantity
  TH1* targetPlotgenHistPS;
  TString targetnamegenHistPS=folderGenKinPhaseSpace+"/"+genHistoName;
  myFile1->GetObject(targetnamegenHistPS, targetPlotgenHistPS);
  if(!targetPlotgenHistPS){
    std::cout << "ERROR: can not open plot " << targetnamegenHistPS << " from file" << myFile1->GetName() << std::endl; 
    exit(0);
  }
  genHistPS = (TH1F*)(targetPlotgenHist->Clone());
  if(debug) std::cout << "genHistPS: "<< folderGenKinPhaseSpace+"/"+genHistoName << std::endl;
  if(plotEfficiencyPhaseSpace){
  TCanvas* Canv4 = new TCanvas(folderGenKinPhaseSpace,folderGenKinPhaseSpace,600,600);
  Canv4->cd();
  //genHistPS->SetStats(kFALSE);
  genHistPS->Draw();
  // rebin histogram of generated quantity
  genHistPS = (TH1F*)genHistPS->Rebin(binvec.size()-1, genHistPS->GetName(), &(binvec.front()));
  //genHistPS->Draw("");
  // calculate efficiency histogram
  effHistBBBPS->Divide(genHistPS);
    
  // ----------- BBB Acceptance --------------------------------
  accHistBBB=(TH1F*)(effHistBBB->Clone());
  accHistBBB->Divide(effHistBBBPS);
  }
 
  // if 2 effiency plots are supposed to be compared:
  TH1F* effHistBBB2 = 0;
  if(plotEfficiency2){
    // get histogram of generated quantity
    TH1* targetPlotgenHist2;
    TString targetnamegenHist2=folderGenKin+"/"+genHistoName;
    myFile1->GetObject(targetnamegenHist2, targetPlotgenHist2);
    if(!targetPlotgenHist2){
      std::cout << "ERROR: can not open plot " << targetnamegenHist2 << " from file" << myFile1->GetName() << std::endl; 
      exit(0);
    }
    TH1F* genHist2 = (TH1F*)(targetPlotgenHist->Clone());
    if(debug) std::cout << "genHist2: "<< folderGenKin+"/"+genHistoName << std::endl;
    TCanvas* Canv4 = new TCanvas(folderGenKin+"3",folderGenKin+"3",600,600);
    Canv4->cd();
    genHist2->Draw();
    // rebin histogram of generated quantity
    genHist2 = (TH1F*)genHist2->Rebin(binvec.size()-1, genHist2->GetName(), &(binvec.front()));
    effHistBBB2 = (TH1F*)(((TH1F*)myFile2->Get(folderRecoKin+"/"+variable+recExtHisto))->Clone());
    // rebin histogram of reconstructed quantity
    effHistBBB2 = (TH1F*)effHistBBB2->Rebin(binvec.size()-1, effHistBBB2->GetName(), &(binvec.front()));
    // calculate efficiency histogram
    effHistBBB2->Divide(genHist2);
  }
  
  //------------------------------------------------------------------------------------
  // ----------- Pure Efficiency calculation--------------------------------
  
  // create histo with efficiency counting all events stemming from a gen bin that are reconstructed anywhere 
  // wrt. all events initially created in that bin, which are reconstructed or not 
  TH1F* effHistGen   = new TH1F("effHistGen","effHistGen", NxBins, &xBins[0]);
  TH1F* effHistGenPS = new TH1F("effHistGenPS","effHistGenPS", NxBins, &xBins[0]);
  TH1F* accHistGen   = new TH1F("accHist","accHist", NxBins, &xBins[0]);
  
  for(int iGenBin=1; iGenBin <=NxBins; iGenBin++){
    if(debug) std::cout << "genbin #" << iGenBin << "/" << NxBins << std::endl;
    double allRecEvtsFromGenBin=0.;
    double effDenom=-1.;
    double effDenomPS=-1.;
    double effGen=0.;
    double effGenPS=0.;
    double acceptance=0.;
    double effDenomMigrationOnly = responseMatrix->Integral(iGenBin,iGenBin,0,NxBins+1);
    if     (genHistPS&&genHistPS->GetNbinsX()>=iGenBin&&genHistPS->GetBinContent(iGenBin)) effDenomPS=genHistPS->GetBinContent(iGenBin);
    if     (genHist  &&genHist  ->GetNbinsX()>=iGenBin&&genHist  ->GetBinContent(iGenBin)) effDenom  =genHist  ->GetBinContent(iGenBin);
    for(int iRecBin=1; iRecBin <=NxBins; iRecBin++){
      double binContentResponse=0.;
      double binContentMigration=0.;
      if(debug){
	std::cout << "   recbin #" << iRecBin << "/" << NxBins << std::endl;
	std::cout << "   effDenom: " << effDenom << std::endl;
	std::cout << "   effDenomMigrationOnly: " << effDenomMigrationOnly << std::endl;
	std::cout << "   effDenomPS: " << effDenomPS << std::endl;
	std::cout << "   responseMatrix->GetBinContent(iGenBin, iRecBin)=" << responseMatrix->GetBinContent(iGenBin, iRecBin) << std::endl;
      }
      if(effDenom>0)              binContentResponse  = responseMatrix->GetBinContent(iGenBin, iRecBin)/ effDenom;
      if(effDenomMigrationOnly>0) binContentMigration = responseMatrix->GetBinContent(iGenBin, iRecBin)/ effDenomMigrationOnly;
      //binContent=99;
      //  	if(verbose>1) std::cout<<"num= "<<responseMatrix->GetBinContent(iGenBin, iRecBin) <<"; denom= "<<effDenom<<"; binContentResponse= "<< binContent<<std::endl;
      allRecEvtsFromGenBin+=responseMatrix->GetBinContent(iGenBin, iRecBin);
      responseMatrix  ->SetBinContent(iGenBin, iRecBin, binContentResponse*100);  // in %
      migrationMatrix ->SetBinContent(iGenBin, iRecBin, binContentMigration*100); // in %
    }
    if(effDenom>0)   effGen    =(double)allRecEvtsFromGenBin/effDenom;
    if(effDenomPS>0) effGenPS  =(double)allRecEvtsFromGenBin/effDenomPS;
    if(effDenom>0)   acceptance=(double)effDenomPS/effDenom;
    //    if(verbose>1){
    // 	std::cout<<"-------------------------"<<std::endl;
    // 	std::cout<<"EffNum= "<<allRecEvtsFromGenBin<<"; Eff.= "<<effGen <<std::endl;
    // 	std::cout<<"effDenomMigrationOnly= "<<effDenomMigrationOnly <<std::endl;
    //       std::cout<<"acceptance = "<<acceptance <<std::endl;
    //       std::cout<<"effDenom = "<<effDenom<<std::endl;
    //       std::cout<<"effDenomPS = "<<effDenomPS<<std::endl;
    // 	std::cout<<"-------------------------"<<std::endl;
    //    }
    effHistGen  ->SetBinContent(iGenBin, effGen);
    effHistGenPS->SetBinContent(iGenBin, effGenPS);
    accHistGen  ->SetBinContent(iGenBin, acceptance);
  }

  // draw purity, stability, pure efficiency and acceptance
  // canvas for purity, stability, pure efficiency and acceptance
  TCanvas* purstab = new TCanvas("purstab","purstab",600,600);
  // compare pure efficiency and bin-by-bin correction factor
  TCanvas* CanvEffGenBBBcomp = new TCanvas("CanvEffGenBBBcomp","CanvEffGenBBBcomp",600,600);
    
  purityhist->SetTitle("");
  //     TString xtitle = "";
  //     TString xtitleNoUnit = "";
  //     if(variable.Contains("Pt"))xtitle+="p_{T}";
  //     else if(variable.Contains("Y"))xtitle+="y";
  //     else if(variable.Contains("Eta"))xtitle+="#eta";
  //     else if(variable.Contains("Mass"))xtitle+="m";
  //     if(variable.Contains("top"))xtitle+="^{t and #bar{t}}";
  //     else if(variable.Contains("ttbar"))xtitle+="^{t#bar{t}}";
  //     else if(variable.Contains("lep"))xtitle+="^{l}";
  //     else if(variable.Contains("bq"))xtitle+="^{b and #bar{b}}";
  //     xtitleNoUnit = xtitle;
  //     if(variable.Contains("Pt"))xtitle+=" [GeV]";
  //     else if(variable.Contains("Mass"))xtitle+=" [GeV]";
  //     if(variable=="topWAngle")xtitle="Angle(top,W)";
  TString xtitle = xSecLabelName(variable);
  TString xtitleNoUnit = xSecLabelName(variable,true);
  if(variable=="topWAngle")xtitle="Angle(top,W)";
  purityhist->GetXaxis()->SetTitle(xtitle);
  purityhist->GetXaxis()->SetNoExponent(true);
  purityhist->GetYaxis()->SetNoExponent(true);

  //purityhist->GetXaxis()->SetTitleSize  ( 0.07);
  //purityhist->GetXaxis()->SetTitleColor (    1);
  //purityhist->GetXaxis()->SetTitleOffset(  0.95);
  //purityhist->GetXaxis()->SetTitleFont  (   62);
  //purityhist->GetXaxis()->SetLabelSize  ( 0.07);
  //purityhist->GetXaxis()->SetLabelFont  (   62);
  //purityhist->GetXaxis()->SetNdivisions (  505);
  purityhist->GetXaxis()->SetRangeUser  (rangeUserLeft, rangeUserRight);
  
  //purityhist->GetYaxis()->SetTitleSize  ( 0.07);
  //purityhist->GetYaxis()->SetTitleColor (    1);
  //purityhist->GetYaxis()->SetTitleOffset(  1.0);
  //purityhist->GetYaxis()->SetTitleFont  (   62);
  //purityhist->GetYaxis()->SetLabelSize  ( 0.07);
  //purityhist->GetYaxis()->SetLabelFont  (   62);
  //purityhist->GetYaxis()->SetNdivisions (  505);
  purityhist->GetYaxis()->SetRangeUser  (0, 0.2);

  purityhist->SetMinimum(0.);
  //double max=purityhist->GetMaximum();
  //if(stabilityhist->GetMaximum()>max)max=stabilityhist->GetMaximum();
  //purityhist->SetMaximum(1.1*max);
  purityhist->SetMaximum(1.04);
  purityhist->SetLineColor(2);
  purityhist->SetLineWidth(4);
  stabilityhist->SetLineColor(4);
  stabilityhist->SetLineStyle(2);
  stabilityhist->SetLineWidth(4);
  accHistGen->SetLineColor(1);
  accHistGen->SetLineStyle(3);
  accHistGen->SetLineWidth(4);
  accHistGen->GetXaxis()->SetNoExponent(true);
  accHistGen->GetYaxis()->SetNoExponent(true);
  purityhist->SetStats(kFALSE);
  purityhist->GetYaxis()->SetNdivisions(511);
  purityhist2->SetLineColor(8);
  purityhist2->SetLineWidth(4);
  stabilityhist2->SetLineColor(7);
  stabilityhist2->SetLineStyle(2);
  stabilityhist2->SetLineWidth(4);

  //gPad->SetBottomMargin(0.19);
  purstab->cd();
  //purstab->SetGrid(1,1);
  purityhist->Draw();
  stabilityhist->Draw("same");
  if(plotDefToCompare){
    purityhist2   ->Draw("same");
    stabilityhist2->Draw("same");
  }
  if(plotAcceptance)accHistGen->Draw("same");
  if(plotEfficiencyPhaseSpace){
    effHistGenPS->SetLineColor(1);
    effHistGenPS->SetLineStyle(1);
    effHistGenPS->SetLineWidth(4);
    effHistGenPS->Draw("same");
  }
  if(plotEfficiency2){
    effHistBBB2->SetLineColor(1);
    effHistBBB2->SetLineStyle(2);
    effHistBBB2->SetLineWidth(4);
    effHistBBB2->Draw("same");
  }
  // chi2 cut efficiency
  if(nameExt!="ProbSel"&&useTree && chi2Max<100){    
    chi2eff = (TH1F*)chi2eff->Rebin(binvec.size()-1, chi2eff->GetName(), &(binvec.front()));
    all     = (TH1F*)all    ->Rebin(binvec.size()-1, all->GetName()    , &(binvec.front()));
    chi2eff->Divide(all);
    chi2eff->SetLineColor(1);
    if(plotDefToCompare) chi2eff->SetLineStyle(9);
    chi2eff->SetLineWidth(4);
    if(drawChi2Eff) chi2eff->Draw("same");
  }
  if (lepton=="muon") DrawDecayChLabel("#mu + Jets");
  else if(lepton=="elec") DrawDecayChLabel("e + Jets");
  else if(lepton=="combined") DrawDecayChLabel("e/#mu + Jets Combined");
  //double legEdge = 0.4;
  //if(plotAcceptance)legEdge = effHistBBB->GetMinimum();
  TLegend* leg=new TLegend(0.47,0.63,0.67,0.82);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  TString addLabel= (chi2MaxInit<100&&((plotDefToCompare)||(!useTree&&nameExt!=""))) ? " (#chi^{2}<"+getTStringFromDouble(chi2MaxInit)+")" : "";
  if(plotDefToCompare) leg->AddEntry(purityhist2,   "Purity"    ,"l");
  leg->AddEntry(purityhist,   TString("Purity"   )+addLabel ,"l");
  if(plotDefToCompare) leg->AddEntry(stabilityhist2,"Stability" ,"l");
  leg->AddEntry(stabilityhist,TString("Stability")+addLabel ,"l");
  //purityhist->GetYaxis()->SetRangeUser  (0, 0.2);
  //if(plotAcceptance)leg->AddEntry(effHistBBB,"Eff*A full PS Spring11","l");
  if(plotEfficiencyPhaseSpace)leg->AddEntry(effHistGenPS,"Efficiency","l");
  if(plotAcceptance)leg->AddEntry(accHistGen,"Acceptance","l");
  if(plotEfficiency2)leg->AddEntry(effHistBBB2,"eff.*A full PS Summer11","l");
  if(chi2Max<100){
    TString entry="eff (#chi^{2} < ";
    entry+=getTStringFromDouble(chi2Max)+")";
    if(drawChi2Eff) leg->AddEntry(chi2eff,entry,"l");
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  if(drawSepLegend&&(variable == "topPt"||drawAlwaysLegend)) leg->Draw("same");
    
  TLegend* legFull=new TLegend(0.1,0.1,1.,0.9);
  legFull->SetFillStyle(0);
  legFull->SetBorderSize(0);
  legFull->AddEntry(purityhist,   "Purity"    ,"l");
  legFull->AddEntry(stabilityhist,"Stability" ,"l");
  double Nlegfull=2.;
  //purityhist->GetYaxis()->SetRangeUser  (0, 0.2);
  //if(plotAcceptance)leg->AddEntry(effHistBBB,"Eff*A full PS Spring11","l");
  //     if(plotEfficiencyPhaseSpace)legFull->AddEntry(effHistGenPS,"Efficiency","l");
  if(plotEfficiencyPhaseSpace){legFull->AddEntry(effHistGenPS,"#epsilon^{vis. PS}","l"); Nlegfull++; }
  //     if(plotAcceptance)legFull->AddEntry(accHistGen,"Acceptance","l");
  if(plotAcceptance ){ legFull->AddEntry(accHistGen ,"A"                      ,"l"); Nlegfull++; }
  if(plotEfficiency2){ legFull->AddEntry(effHistBBB2,"eff.*A full PS Summer11","l"); Nlegfull++; }
  if(chi2Max<100){
    TString entry="eff (#chi^{2} < ";
    entry+=getTStringFromDouble(chi2Max)+")";
    if(drawChi2Eff){ legFull->AddEntry(chi2eff,entry,"l"); Nlegfull++; }
  }
  legFull->SetTextSize(0.1);
  // canvas for legend
  TCanvas* purstabLeg = new TCanvas("purstabLeg","purstabLeg",600,600);
  legFull->Draw("");
    
  //------------------------------------------------------------------------------
  // draw effGen together with effective efficiencies from bin-by-bin
  // OBSERVATION: acceptance from both methods the same!!! -> overlap
  CanvEffGenBBBcomp->cd();
  effHistBBBPS->GetXaxis()->SetTitle("gen "+xtitle);
  effHistBBBPS->GetXaxis()->SetNoExponent(true);
  effHistBBBPS->GetYaxis()->SetNoExponent(true);
  effHistBBBPS->SetLineColor(kGreen+2);
  effHistBBBPS->SetLineStyle(1);
  effHistBBBPS->SetLineWidth(4);
  effHistBBBPS->GetXaxis()->SetRangeUser  (rangeUserLeft, rangeUserRight);
  effHistBBBPS->GetYaxis()->SetRangeUser  (0., 0.5);
  effHistBBBPS->Draw("");

  if(plotAcceptance){
    accHistBBB->SetLineColor(kGreen+2);
    accHistBBB->SetLineStyle(3);
    accHistBBB->SetLineWidth(4);
    accHistBBB->Draw("same");
  }

  if(plotAcceptance)           accHistGen->Draw("same");
  if(plotEfficiencyPhaseSpace) effHistGenPS->Draw("same");
  
  TLegend* leg2=new TLegend(0.2,0.68,0.67,0.87);
  leg2->SetTextSize(0.03);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  if(plotAcceptance)           leg2->AddEntry(accHistGen,     "A (no migratio)"    ,"l");
  if(plotEfficiencyPhaseSpace) leg2->AddEntry(effHistGenPS,   "eff. (no migration)"    ,"l");
  if(plotAcceptance)           leg2->AddEntry(accHistBBB,     "eff.*C*A (bin-by-bin)","l");
  if(plotEfficiencyPhaseSpace) leg2->AddEntry(effHistBBBPS,   "eff.*C (bin-by-bin)","l");
  leg2->Draw();
    
  //------------------------------------------------------------------------------
  // calculate response matrix as probability matrix: divide number in bin(gen,rec) by number of all
  // initially generated events in the gen bin (i.e. also if they have been lost due to inefficiencies);
  // if restricted phase space shall be plotted, it is calculated for that one
  TCanvas* CanvResponseMatrix = new TCanvas("CanvResponseMatrix","CanvResponseMatrix",600,600);
  // additionally: calculate response matrix that only includes migration effects, not limited acc. or eff.
  TCanvas* CanvMigrationMatrix = new TCanvas("CanvMigrationMatrix","CanvMigrationMatrix",600,600);
  // collection of residual histos for determination of resolution
  //     TCanvas* CanvResidualHistos = new TCanvas("CanvResidualHistos","CanvResidualHistos",1200,1200);
  //     TCanvas* CanvRelativeResidualHistos = new TCanvas("CanvRelativeResidualHistos","CanvRelativeResidualHistos",1200,1200);
  //     int nDivide = (int)sqrt(NxBins) +1;
  //     CanvResidualHistos->Divide(nDivide,nDivide);
  //     CanvRelativeResidualHistos->Divide(nDivide,nDivide);
  int nRows = (int)((NxBins+2)/3); // number of rows for 3 columns
  TCanvas* CanvResidualHistos = new TCanvas("CanvResidualHistos","CanvResidualHistos",840,200*nRows);
  TCanvas* CanvRelativeResidualHistos = new TCanvas("CanvRelativeResidualHistos","CanvRelativeResidualHistos",840,200*nRows);
  CanvResidualHistos->Divide(3,nRows);
  CanvRelativeResidualHistos->Divide(3,nRows);
  // resolution
  TCanvas* CanvResolution         = new TCanvas("CanvResolution","CanvResolution",600,600);
  TCanvas* CanvRelativeResolution = new TCanvas("CanvRelativeResolution","CanvRelativeResolution",600,600);
  TCanvas* CanvBinwidthOverResol  = new TCanvas("CanvBinwidthOverResol","CanvBinwidthOverResol",600,600);
  TH1F* hisRelativeResolRMS=0;
  TH1F* hisRelativeResolFit=0;
  // calculate resolution
  if(createResolutionPlots){
    //----------------------------------------------------
    // draw residuals and calculate resolution
    TH1F* hisResolRMS = new TH1F("hisResolRMS","hisResolRMS", NxBins, &xBins[0]);
    TH1F* hisResolFit = new TH1F("hisResolFit","hisResolFit", NxBins, &xBins[0]);
    hisRelativeResolRMS = new TH1F("hisRelativeResolRMS","hisRelativeResolRMS", NxBins, &xBins[0]);
    hisRelativeResolFit = new TH1F("hisRelativeResolFit","hisRelativeResolFit", NxBins, &xBins[0]);
    TH1F* hisBinwidthOverResolRMS = new TH1F("hisBinwidthOverResolRMS","hisBinwidthOverResolRMS", NxBins, &xBins[0]);
    TH1F* hisBinwidthOverResolFit = new TH1F("hisBinwidthOverResolFit","hisBinwidthOverResolFit", NxBins, &xBins[0]);
    for(int iGenBin=0; iGenBin <NxBins; iGenBin++){
      // absolute resolution	
      CanvResidualHistos->cd(iGenBin+1);
      residualHistos_.at(iGenBin)->GetXaxis()->SetNoExponent(true);
      residualHistos_.at(iGenBin)->GetXaxis()->SetTitle(xtitle+" (rec-gen)");
      residualHistos_.at(iGenBin)->GetYaxis()->SetTitle("N");
      residualHistos_.at(iGenBin)->GetXaxis()->SetTitleSize(0.06);
      residualHistos_.at(iGenBin)->GetXaxis()->SetLabelSize(0.06);
      residualHistos_.at(iGenBin)->GetXaxis()->SetTitleOffset(1.1);
      residualHistos_.at(iGenBin)->GetYaxis()->SetTitleSize(0.06);
      residualHistos_.at(iGenBin)->GetYaxis()->SetLabelSize(0.06);
      residualHistos_.at(iGenBin)->GetYaxis()->SetTitleOffset(1.1);
      residualHistos_.at(iGenBin)->SetLineWidth(2);
      residualHistos_.at(iGenBin)->Draw();
      // fit
      //     TF1 *f1 = new TF1("f1","gaus",-200,200);
      //     f1->SetLineWidth(2);    
      std::vector<double> f1sigma_;
      if(fitGaussRes) f1sigma_=extract_sigma(residualHistos_.at(iGenBin),verbose);
      //     residualHistos_.at(iGenBin)->Fit("f1");
      //     double f1norm  = f1->GetParameter(0);
      //     double f1mean  = f1->GetParameter(1);
      //     double f1sigma = f1->GetParameter(2);
      //residualHistos_.at(iGenBin)->GetFunction("f1");
      TPaveText *label = new TPaveText(0.7, 0.6, 0.9, 0.8, "br NDC");
      label->AddText("Fit:");
      //  label->AddText(Form("x_{0} = %3.2f",f1mean));
      if(fitGaussRes) label->AddText(Form("#sigma = %3.2f",f1sigma_[0]));
      label->SetFillStyle(0);
      label->SetBorderSize(0);
      label->SetTextSize(0.04);
      label->SetTextAlign(12);
      if(fitGaussRes) label->Draw("same");
      TString labelText  = TString("       ")+ Form("%4.1f",xBins[iGenBin])+"<"+xtitle+"<"+Form("%4.1f",xBins[iGenBin+1]);
      DrawLabel(labelText,gStyle->GetPadLeftMargin(),1.0-gStyle->GetPadTopMargin(),1.-gStyle->GetPadRightMargin(),1.,12,0.07);
      hisResolRMS->SetBinContent(iGenBin+1, residualHistos_.at(iGenBin)->GetRMS());
      if(fitGaussRes) hisResolFit->SetBinContent(iGenBin+1, f1sigma_[0]);
      hisBinwidthOverResolRMS->SetBinContent(iGenBin+1, hisBinwidthOverResolRMS->GetBinWidth(iGenBin+1)/(residualHistos_.at(iGenBin)->GetRMS()));
      if(fitGaussRes) hisBinwidthOverResolFit->SetBinContent(iGenBin+1, hisBinwidthOverResolFit->GetBinWidth(iGenBin+1)/f1sigma_[0]);
    
      //------------------
      // relative resolution	
      CanvRelativeResidualHistos->cd(iGenBin+1);
      relativeResidualHistos_.at(iGenBin)->GetXaxis()->SetNoExponent(true);
      relativeResidualHistos_.at(iGenBin)->GetXaxis()->SetTitle(xtitleNoUnit+" (rec-gen)/gen");
      relativeResidualHistos_.at(iGenBin)->GetYaxis()->SetTitle("N");
      relativeResidualHistos_.at(iGenBin)->GetXaxis()->SetTitleSize(0.06);
      relativeResidualHistos_.at(iGenBin)->GetXaxis()->SetLabelSize(0.06);
      relativeResidualHistos_.at(iGenBin)->GetXaxis()->SetTitleOffset(1.1);
      relativeResidualHistos_.at(iGenBin)->GetYaxis()->SetTitleSize(0.06);
      relativeResidualHistos_.at(iGenBin)->GetYaxis()->SetLabelSize(0.06);
      relativeResidualHistos_.at(iGenBin)->GetYaxis()->SetTitleOffset(1.1);
      relativeResidualHistos_.at(iGenBin)->SetLineWidth(2);
      relativeResidualHistos_.at(iGenBin)->Draw();
    
      // fit
      //       TF1 *f1 = new TF1("f1","gaus",-200,200);
      //       f1->SetLineWidth(2);    
      std::vector<double> f1RelativeSigma_;
      if(fitGaussRes) f1RelativeSigma_=extract_sigma(relativeResidualHistos_.at(iGenBin), verbose);
      //     residualHistos_.at(iGenBin)->Fit("f1");
      //     double f1norm  = f1->GetParameter(0);
      //     double f1mean  = f1->GetParameter(1);
      //     double f1sigma = f1->GetParameter(2);
      //residualHistos_.at(iGenBin)->GetFunction("f1");
      TPaveText *label2 = new TPaveText(0.7, 0.6, 0.9, 0.8, "br NDC");
      label2->AddText("Fit:");
      //  label2->AddText(Form("x_{0} = %3.2f",f1mean));
      if(fitGaussRes) label2->AddText(Form("#sigma = %3.2f",f1RelativeSigma_[0]));
      label2->SetFillStyle(0);
      label2->SetBorderSize(0);
      label2->SetTextSize(0.04);
      label2->SetTextAlign(12);
      if(fitGaussRes) label2->Draw("same");
      TString label2Text  = TString("       ")+ Form("%4.1f",xBins[iGenBin])+"<"+xtitle+"<"+Form("%4.1f",xBins[iGenBin+1]);
      DrawLabel(label2Text,gStyle->GetPadLeftMargin(),1.0-gStyle->GetPadTopMargin(),1.-gStyle->GetPadRightMargin(),1.,12,0.07);
    
      // calculate resolution
      hisRelativeResolRMS->SetBinContent(iGenBin+1, relativeResidualHistos_.at(iGenBin)->GetRMS());
      if(fitGaussRes) hisRelativeResolFit->SetBinContent(iGenBin+1, f1RelativeSigma_[0]);
    }

    // draw resolution
    CanvResolution->cd();
    hisResolRMS->GetXaxis()->SetTitle("gen "+xtitle);
    hisResolRMS->GetYaxis()->SetTitle("Resolution of "+xtitle);
    //hisResolRMS->GetYaxis()->SetTitleOffset(1.1);
    hisResolRMS->GetXaxis()->SetNoExponent(true);
    //hisResolRMS->GetYaxis()->SetNoExponent(true);
    hisResolRMS->SetLineColor(1);
    hisResolRMS->SetLineStyle(1);
    hisResolRMS->SetLineWidth(2);
    hisResolRMS->GetXaxis()->SetRangeUser  (rangeUserLeft, rangeUserRight);
    hisResolRMS->GetYaxis()->SetRangeUser  (0., 1.2*hisResolRMS->GetMaximum());
    hisResolRMS->Draw("");
    hisResolFit->SetLineColor(2);
    hisResolFit->SetLineStyle(1);
    hisResolFit->SetLineWidth(2);
    if(fitGaussRes) hisResolFit->Draw("same");
    TLegend* legResol=new TLegend(0.2,0.68,0.67,0.87);
    legResol->SetTextSize(0.03);
    legResol->SetFillStyle(0);
    legResol->SetBorderSize(0);
    legResol->AddEntry(hisResolRMS,   "#sigma_{histo}"    ,"l");
    if(fitGaussRes) legResol->AddEntry(hisResolFit,   "#sigma_{Gauss fit}"    ,"l");
    if(fitGaussRes) legResol->DrawClone("same");
  
    // draw relative resolution
    CanvRelativeResolution->cd();
    hisRelativeResolRMS->GetXaxis()->SetTitle("gen "+xtitle);
    hisRelativeResolRMS->GetYaxis()->SetTitle("Relative Resolution of "+xtitleNoUnit);
    //hisRelativeResolRMS->GetYaxis()->SetTitleOffset(1.1);
    hisRelativeResolRMS->GetXaxis()->SetNoExponent(true);
    //hisRelativeResolRMS->GetYaxis()->SetNoExponent(true);
    hisRelativeResolRMS->SetLineColor(1);
    hisRelativeResolRMS->SetLineStyle(1);
    hisRelativeResolRMS->SetLineWidth(2);
    hisRelativeResolRMS->GetXaxis()->SetRangeUser  (rangeUserLeft, rangeUserRight);
    hisRelativeResolRMS->GetYaxis()->SetRangeUser  (0., 1.2*hisRelativeResolRMS->GetMaximum());
    hisRelativeResolRMS->Draw("");
    hisRelativeResolFit->SetLineColor(2);
    hisRelativeResolFit->SetLineStyle(1);
    hisRelativeResolFit->SetLineWidth(2);
    if(fitGaussRes) hisRelativeResolFit->Draw("same");
    //   TLegend* legResol=new TLegend(0.2,0.68,0.67,0.87);
    //   legResol->SetTextSize(0.03);
    //   legResol->SetFillStyle(0);
    //   legResol->SetBorderSize(0);
    //   legResol->AddEntry(hisResolRMS,   "#sigma_{histo}"    ,"l");
    //   legResol->AddEntry(hisResolFit,   "#sigma_{Gauss fit}"    ,"l");
    if(fitGaussRes) legResol->DrawClone("same");
    // draw binwidth over resolution
    CanvBinwidthOverResol->cd();
    hisBinwidthOverResolRMS->GetXaxis()->SetTitle("gen "+xtitle);
    hisBinwidthOverResolRMS->GetYaxis()->SetTitle("Bin Width / Resolution");
    hisBinwidthOverResolRMS->GetXaxis()->SetNoExponent(true);
    hisBinwidthOverResolRMS->GetYaxis()->SetNoExponent(true);
    hisBinwidthOverResolRMS->SetLineColor(1);
    hisBinwidthOverResolRMS->SetLineStyle(1);
    hisBinwidthOverResolRMS->SetLineWidth(2);
    hisBinwidthOverResolRMS->GetXaxis()->SetRangeUser  (rangeUserLeft, rangeUserRight);
    hisBinwidthOverResolRMS->GetYaxis()->SetRangeUser  (0., 3.*hisBinwidthOverResolRMS->GetMaximum());
    hisBinwidthOverResolRMS->Draw("");
    hisBinwidthOverResolFit->SetLineColor(2);
    hisBinwidthOverResolFit->SetLineStyle(1);
    hisBinwidthOverResolFit->SetLineWidth(2);
    if(fitGaussRes) hisBinwidthOverResolFit->Draw("same");
    //   TLegend* legResol=new TLegend(0.2,0.68,0.67,0.87);
    //   legResol->SetTextSize(0.03);
    //   legResol->SetFillStyle(0);
    //   legResol->SetBorderSize(0);
    //   legResol->AddEntry(hisResolRMS,   "#sigma_{histo}"    ,"l");
    //   legResol->AddEntry(hisResolFit,   "#sigma_{Gauss fit}"    ,"l");
    if(fitGaussRes) legResol->DrawClone("same");
  }

  // draw responseMatrix incl. eff.
  // text format for TH2 labels
  gStyle->SetPaintTextFormat("3.0f");
  CanvResponseMatrix->cd();
  CanvResponseMatrix->SetRightMargin(0.2);
  responseMatrix->SetStats(kFALSE);
  responseMatrix->GetXaxis()->SetTitle("gen "+xtitle);
  responseMatrix->GetYaxis()->SetTitle("rec "+xtitle);
  responseMatrix->GetXaxis()->SetRangeUser(rangeUserLeft,rangeUserRight);
  responseMatrix->GetYaxis()->SetRangeUser(rangeUserLeft,rangeUserRight);
  responseMatrix->GetXaxis()->SetNoExponent(true);
  responseMatrix->GetYaxis()->SetNoExponent(true);
  responseMatrix->GetZaxis()->SetNoExponent(true);
  responseMatrix->GetZaxis()->SetRangeUser(0,30);
  responseMatrix->GetZaxis()->SetTitle("Transition Probability [\%]");
  responseMatrix->GetZaxis()->SetTitleOffset(1.4);
  //responseMatrix->SetMarkerColor(kWhite);
  responseMatrix->SetMarkerSize(1.2);
  responseMatrix->SetTitle("Response Matrix - Including Efficiency");
  responseMatrix->Draw(Form("COLZ %s",textOpt.Data()));
  //TPaletteAxis *palette = (TPaletteAxis*)responseMatrix->GetListOfFunctions()->FindObject("palette");
  //gPad->Update();
  if(plotEfficiencyPhaseSpace) DrawLabel("Full Response Matrix - Migration and Efficiency",gStyle->GetPadLeftMargin(),1.0-gStyle->GetPadTopMargin(),1.-gStyle->GetPadRightMargin(),1.,12,0.03);
  else DrawLabel("Full Response Matrix - Migration and Eff.*Acc.",gStyle->GetPadLeftMargin(),1.0-gStyle->GetPadTopMargin(),1.-gStyle->GetPadRightMargin(),1.,12,0.03);

  // draw migrationMatrix w/o eff.
  CanvMigrationMatrix->cd();
  //CanvMigrationMatrix->SetRightMargin(0.2);
  CanvMigrationMatrix->SetRightMargin(0.13);
  migrationMatrix->SetStats(kFALSE);
  migrationMatrix->GetXaxis()->SetTitle("gen "+xtitle);
  migrationMatrix->GetYaxis()->SetTitle("rec "+xtitle);
  migrationMatrix->GetXaxis()->SetRangeUser(rangeUserLeft,rangeUserRight);
  migrationMatrix->GetYaxis()->SetRangeUser(rangeUserLeft,rangeUserRight);
  migrationMatrix->GetXaxis()->SetNoExponent(true);
  migrationMatrix->GetYaxis()->SetNoExponent(true);
  migrationMatrix->GetZaxis()->SetNoExponent(true);
  migrationMatrix->GetZaxis()->SetRangeUser(0,100);
  migrationMatrix->SetTitle("");
  //migrationMatrix->GetZaxis()->SetTitle("Transition Probability [\%]");
  migrationMatrix->GetZaxis()->SetTitleOffset(1.3);
  //migrationMatrix->SetMarkerColor(kWhite);
  //   migrationMatrix->SetMarkerSize(1.2);
  migrationMatrix->SetMarkerSize(1.4);
  if(variable=="lepPt") migrationMatrix->SetMarkerSize(1.);
  migrationMatrix->SetTitle("Migration Matrix - Transition Probability [\%]");
  migrationMatrix->Draw(Form("COLZ %s",textOpt.Data()));
  DrawLabel("Migration Matrix - Transition Probability [\%]", gStyle->GetPadLeftMargin()-0.03, 1.0-gStyle->GetPadTopMargin(), 1.-gStyle->GetPadRightMargin(), 1., 12,0.04);
   
  if(save){
    if(verbose==0) gErrorIgnoreLevel=kWarning; 
    TString chi="";
    if(chi2Max<100){ 
      chi+="chi";
      chi+=(getTStringFromDouble(chi2Max).ReplaceAll(".", "p"));
    }
    TString qAssStr="";
    if(qAssignment>=0){ 
      qAssStr+="Perm";
      if(qAssignment==0)      qAssStr+="Ok";
      else if(qAssignment==1) qAssStr+="b<->b";
      else if(qAssignment==9) qAssStr+="NoMatch";
      else                    qAssStr+=getTStringFromInt(qAssignment);
    }
    purstab->Print(outputFolder+"/purStabEff_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
    purstab->Print(outputFolder+"/purStabEff_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
    purstabLeg->Print(outputFolder+"/purStabLeg_"+universalplotLabel+"_"+".png");
    purstabLeg->Print(outputFolder+"/purStabLeg_"+universalplotLabel+"_"+".eps");
    //png
    CanvResponseMatrix->Print(outputFolder+"/responseMatrix_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
    CanvMigrationMatrix->Print(outputFolder+"/migrationMatrix_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
    //eps
    CanvResponseMatrix->Print(outputFolder+"/responseMatrix_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
    CanvMigrationMatrix->Print(outputFolder+"/migrationMatrix_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
    if(useTree){
      CanvEffGenBBBcomp->Print(outputFolder+"/effGenBBBcomp_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
      CanvResidualHistos->Print(outputFolder+"/residualHistos_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
      CanvRelativeResidualHistos->Print(outputFolder+"/relativeResidualHistos_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
      CanvResolution->Print(outputFolder+"/resolution_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
      CanvRelativeResolution->Print(outputFolder+"/relativeResolution_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
      //eps
      CanvEffGenBBBcomp->Print(outputFolder+"/effGenBBBcomp_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
      CanvResidualHistos->Print(outputFolder+"/residualHistos_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
      CanvRelativeResidualHistos->Print(outputFolder+"/relativeResidualHistos_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
      CanvResolution->Print(outputFolder+"/resolution_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
      CanvRelativeResolution->Print(outputFolder+"/relativeResolution_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
      // print residuums separately if desired
      if(printSeparateRes){
	for(int iBin =0; iBin < NxBins; iBin++){
	  CanvResidualHistos->cd(iBin+1)->Print(outputFolder+"/residualHistoBin"+Form("%i",iBin)+"_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
	  CanvRelativeResidualHistos->cd(iBin+1)->Print(outputFolder+"/relativeResidualHistoBin"+Form("%i",iBin)+"_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".png");
	  CanvResidualHistos->cd(iBin+1)->Print(outputFolder+"/residualHistoBin"+Form("%i",iBin)+"_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");
	  CanvRelativeResidualHistos->cd(iBin+1)->Print(outputFolder+"/relativeResidualHistoBin"+Form("%i",iBin)+"_"+universalplotLabel+"_"+lepton+"_"+variable+chi+qAssStr+".eps");	}
      }
    }
    if(verbose>1) std::cout<<"Canvas with purity and stability plots is saved in "<<outputFolder<<std::endl;
  }

  // control value for pur/stab
  double minMigr=0.4;
  // print range of migration variables
  double maxStab=0.;
  double maxPur =0.;
  double minStab=1.;
  double minPur =1.;
  double minRes =1.;
  double maxRes =0.;
  int purityMinimumBin=0;
  int purityMaximumBin=0;
  int stabilityMinimumBin=0;
  int stabilityMaximumBin=0;
  int resolutionMinimumBin=0;
  int resolutionMaximumBin=0;
  int Nbins=0;
  int NPurok=0;
  int NStabok=0;
  // loop bins
  for(int bin=1; bin<=purityhist->GetNbinsX(); ++bin){
    // check range
    if(purityhist->GetBinCenter(bin)>=rangeUserLeft&&purityhist->GetBinCenter(bin)<=rangeUserRight){
      // migration in this bin
      double TMPStab=stabilityhist->GetBinContent(bin);
      double TMPPur =purityhist   ->GetBinContent(bin);
      double TMPRes = createResolutionPlots ? hisRelativeResolRMS->GetBinContent(bin) : -1.;
      if(TMPStab<1.&&TMPPur>0.){
	// count all bins used for measurement
	++Nbins;
	// search for min/max values
	if(TMPStab<1.&&TMPStab>maxStab){ maxStab= TMPStab; stabilityMaximumBin=bin;}
	if(TMPPur<1. &&TMPPur>maxPur  ){ maxPur = TMPPur; purityMaximumBin=bin;    }
	if(TMPStab>0.&&TMPStab<minStab){ minStab= TMPStab; stabilityMinimumBin=bin;}
	if(TMPPur>0. &&TMPPur<minPur  ){ minPur = TMPPur; purityMinimumBin=bin;    }
	if(TMPRes>-1){
	  if(TMPRes>0.&&TMPRes<minRes){ minRes= TMPRes; resolutionMinimumBin=bin;}
	  if(TMPRes<1.&&TMPRes>maxRes){ maxRes= TMPRes; resolutionMaximumBin=bin;}
	}
	// check minimal migration condition
	if(TMPStab>minMigr)NStabok++;
	if(TMPPur >minMigr)NPurok++;
      }
    }
  }
  gErrorIgnoreLevel=initialIgnoreLevel;

  // contol output
  std::cout << std::endl << variable << "(" << lepton << "):" << std::endl;
  // purity
  std::cout << "--------------------" << std::endl;
  std::cout << "a) purity" << std::endl;
  // all bins
  if(verbose>0){
    for(int bin=1; bin<=purityhist->GetNbinsX(); ++bin){
      std::cout << std::setprecision(2) << std::fixed << "[" << purityhist->GetBinLowEdge(bin); 
      for(int space=(getBigitFromDouble(purityhist->GetBinLowEdge(purityhist->GetNbinsX()+1))-getBigitFromDouble(purityhist->GetBinLowEdge(bin))); space>0; space--){ std::cout << " "; }
      std::cout << ".." <<  purityhist->GetBinLowEdge(bin+1);
      for(int space=(getBigitFromDouble(purityhist->GetBinLowEdge(purityhist->GetNbinsX()+1))-getBigitFromDouble(purityhist->GetBinLowEdge(bin+1))); space>0; space--){ std::cout << " "; }
      std::cout << "]: ";
      std::cout << std::setprecision(2) << std::fixed << purityhist->GetBinContent(bin);
      if(plotDefToCompare){
	std::cout << " (chi2<" << chi2MaxInit << "), no chi2 cut: ";
	std::cout << std::setprecision(2) << std::fixed << purityhist2->GetBinContent(bin);
      }
      std::cout << std::endl;      
    }
  }
  // min/max value
  std::cout << std::setprecision(2) << std::fixed << minPur  << " (<" << variable; 
  std::cout << ">=" << std::setprecision(0) << std::fixed;
  std::cout << purityhist->GetBinCenter(purityMinimumBin) << ") < purity < ";
  std::cout << std::setprecision(2) << std::fixed << maxPur  << " (<" << variable;
  std::cout << ">=" << std::setprecision(0) << std::fixed;
  std::cout << purityhist->GetBinCenter(purityMaximumBin) << ")" << std::endl;
  std::cout << "( purity > " << std::setprecision(2) << std::fixed;
  std::cout << minMigr << " in " << NPurok << "/" << Nbins << " bins)" << std::endl;
  // stability
  std::cout << "--------------------" << std::endl;
  std::cout << "b) stability" << std::endl;
  // all bins
  if(verbose>0){
    for(int bin=1; bin<=stabilityhist->GetNbinsX(); ++bin){
      std::cout << "[" << stabilityhist->GetBinLowEdge(bin);
      for(int space=(getBigitFromDouble(purityhist->GetBinLowEdge(purityhist->GetNbinsX()+1))-getBigitFromDouble(purityhist->GetBinLowEdge(bin))); space>0; space--){ std::cout << " "; }
      std::cout << ".." <<  stabilityhist->GetBinLowEdge(bin+1);
      for(int space=(getBigitFromDouble(purityhist->GetBinLowEdge(purityhist->GetNbinsX()+1))-getBigitFromDouble(purityhist->GetBinLowEdge(bin+1))); space>0; space--){ std::cout << " "; }
      std::cout << "]: ";
      std::cout << std::setprecision(2) << std::fixed << stabilityhist->GetBinContent(bin);
      if(plotDefToCompare){
	std::cout << " (chi2<" << chi2MaxInit << "), no chi2 cut: ";
	std::cout << std::setprecision(2) << std::fixed << stabilityhist2->GetBinContent(bin);
      }
      std::cout << std::endl;
    }
  }
  // min/max value
  std::cout << std::setprecision(2) << std::fixed << minStab << " (<" << variable;
  std::cout << ">=" << std::setprecision(0) << std::fixed;
  std::cout << stabilityhist->GetBinCenter(stabilityMinimumBin) << ") < stability < ";
  std::cout << std::setprecision(2) << std::fixed << maxStab << " (<" << variable;
  std::cout << ">=" << std::setprecision(0) << std::fixed;
  std::cout << stabilityhist->GetBinCenter(stabilityMaximumBin) << ")"<< std::endl;
  std::cout << "( stability > " << std::setprecision(2) << std::fixed << minMigr;
  std::cout << " in " << NStabok << "/" << Nbins << " bins)" << std::endl<< std::endl;

  // gen-rec correlation
  std::cout << "--------------------" << std::endl;
  std::cout << "c) correlation" << std::endl;
  std::cout << "corr(evtMatrix) : " << myHist2d2      ->GetCorrelationFactor() << std::endl;
  std::cout << "corr(respMatrix): " << responseMatrix ->GetCorrelationFactor() << std::endl;
  std::cout << "corr(migrMatrix): " << migrationMatrix->GetCorrelationFactor() << std::endl;
  // resolution
  if(createResolutionPlots){
  std::cout << "--------------------" << std::endl;
  std::cout << "d) resolution" << std::endl;
  // all bins
  double SumRes   =0.;
  double SumResFit=0.;
  for(int bin=1; bin<=hisRelativeResolRMS->GetNbinsX(); ++bin){
    SumRes+=hisRelativeResolRMS->GetBinContent(bin);
    if(fitGaussRes) SumResFit+= hisRelativeResolFit->GetBinContent(bin);
    if(verbose>0){
      std::cout << "[" << hisRelativeResolRMS->GetBinLowEdge(bin);
      for(int space=(getBigitFromDouble(purityhist->GetBinLowEdge(purityhist->GetNbinsX()+1))-getBigitFromDouble(purityhist->GetBinLowEdge(bin)))  ; space>0; space--){ std::cout << " "; }
      std::cout << ".." <<  hisRelativeResolRMS->GetBinLowEdge(bin+1);
      for(int space=(getBigitFromDouble(purityhist->GetBinLowEdge(purityhist->GetNbinsX()+1))-getBigitFromDouble(purityhist->GetBinLowEdge(bin+1))); space>0; space--){ std::cout << " "; }
      std::cout << "]: ";
      std::cout << std::setprecision(2) << std::fixed << hisRelativeResolRMS->GetBinContent(bin);
      if(fitGaussRes) std::cout << " (1#sigma gaussian fit=" << std::setprecision(2) << std::fixed << hisRelativeResolFit->GetBinContent(bin) << ")";
      std::cout << std::endl;
    }
  }
  // average
  std::cout << "<resolution(histo)>=" << SumRes   /hisRelativeResolRMS->GetNbinsX() << std::endl;
  std::cout << "<resolution(gauss)>=" << SumResFit/hisRelativeResolFit->GetNbinsX() << std::endl;
  // min/max value
  if(verbose>1){
    std::cout << std::setprecision(2) << std::fixed << minRes << " (<" << variable;
    std::cout << ">=" << std::setprecision(0) << std::fixed;
    std::cout << hisRelativeResolRMS->GetBinCenter(resolutionMinimumBin) << ") < resolution < ";
    std::cout << std::setprecision(2) << std::fixed << maxRes << " (<" << variable;
    std::cout << ">=" << std::setprecision(0) << std::fixed;
    std::cout << hisRelativeResolRMS->GetBinCenter(resolutionMaximumBin) << ")"<< std::endl;
  }
  }
}
