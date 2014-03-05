#include "basicFunctions.h"

TH1F* distortPDF(const TH1& hist, TString variation, TString variable, TString inputFolderName, TString phaseSpace, int verbose);
TH1F* distort   (const TH1& hist, TString variation, TString variable, int verbose);
double linSF(const double x, const double xmax, const double a, const double b);

void analyzeTopDiffXSecMCdependency(double luminosity = 12148, std::string decayChannel="electron", bool save=true, int verbose=2, TString inputFolderName="newRecentAnalysisRun8TeV",
				    TString dataFile=  groupSpace+AnalysisFolder+"/elecDiffXSecData2012ABCDAll.root",
				    //TString dataFile= groupSpace+AnalysisFolder+"/muonDiffXSecData2012ABCDAll.root"
                                    bool doPDFuncertainty=true, bool addCrossCheckVariables=false)
{
  // ---
  //     Configuration
  // ---
  // folder with all rootfiles
  TString path=groupSpace;
  // take care that prescaling of muon channel for full 2011 datset was taken into account
  if(luminosity==4980&&decayChannel=="muon"    ) luminosity=constLumiMuon;
  if(luminosity==4955&&decayChannel=="electron") luminosity=constLumiElec;
  // set up common analysis style
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  TGaxis::SetMaxDigits(2);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  
  // output
  int initialgErrorLv=gErrorIgnoreLevel;
  if(verbose==1) gErrorIgnoreLevel=kWarning;
  if(verbose==0) gErrorIgnoreLevel=kFatal;
  
  // define names
  // file name for input rootfile
  TString analysisFileName="";
  TString SampleTag="Summer12";
  if     (decayChannel=="muon"    ) analysisFileName=path+inputFolderName+"/"+TopFilename(kSig, sysNo, "muon"    );
  else if(decayChannel=="electron") analysisFileName=path+inputFolderName+"/"+TopFilename(kSig, sysNo, "electron");
  else{
    std::cout << "ERROR: decay channel " << decayChannel << " is no valid choice, use electron or muon!" << std::endl;
    exit(0);
  }
  
  // file name for output rootfiles
  TString outputFileNameUp   = analysisFileName;
  TString outputFileNameDown = analysisFileName;
  if (doPDFuncertainty) {
    outputFileNameUp.ReplaceAll(inputFolderName, inputFolderName+"/PDFUp");
    outputFileNameDown.ReplaceAll(inputFolderName, inputFolderName+"/PDFDown");
    outputFileNameUp.ReplaceAll(SampleTag+"PF", "PdfVarUp"+SampleTag+"PF"  );
    outputFileNameDown.ReplaceAll(SampleTag+"PF", "PdfVarDown"+SampleTag+"PF");
  } 
  else { 
    outputFileNameUp.ReplaceAll(inputFolderName, inputFolderName+"/MCShapeUp");
    outputFileNameDown.ReplaceAll(inputFolderName, inputFolderName+"/MCShapeDown");
    outputFileNameUp.ReplaceAll("PF", "MCShapeVarUp"+SampleTag+"PF"  );
    outputFileNameDown.ReplaceAll("PF", "MCShapeVarDown"+SampleTag+"PF");
  }
  // name for folder where tree is read from
  TString folder="analyzeTopRecoKinematicsKinFit";
  TString genfolder="analyzeTopPartonLevelKinematics";
  TString genfolder2=genfolder;
  genfolder+="PhaseSpace";
  // name for folder where plots are saved to
  TString savePlotsTo="./diffXSecFromSignal/plots/"+decayChannel+"/2012/shapeReweighting";
  // define variables
  std::vector<TString> variable_;
  TString variable[] ={"topPt", "topY", "ttbarPt", "ttbarMass", "ttbarY"};//, "lepPt", "lepEta"};
  // cross check variables
  variable_.insert( variable_.begin(), variable, variable + sizeof(variable)/sizeof(TString) );
  if (addCrossCheckVariables) variable_.insert(variable_.end(),   xSecVariablesCCVar,  xSecVariablesCCVar + sizeof( xSecVariablesCCVar)/sizeof(TString));
  // b quark variables only for PDF uncertainties
  //TString bvariables[]={"bqPt", "bqEta"};
  //if(doPDFuncertainty) variable_.insert(variable_.end(), bvariables, bvariables+sizeof(bvariables)/sizeof(TString));
  
  // container for values read from tree
  std::map< TString, float > value_;
  // container for original and shape distorted
  // parton Level (phase space) plots
  std::map< TString, TH1F*> plots_;
  std::map< TString, TH1F*> plotsScaled_;
  std::map< TString, TH2F*> plots2D_;
  std::map< TString, TH2F*> plotsScaled2D_;
  // container for final distorted plots
  // lep/had top quantities are filled in the same plot
  std::map< TString, std::map< TString, TH1F*> > finalPlots_;
  std::map< TString, std::map< TString, TH2F*> > finalPlots2D_;
  
  // printout of configuration
  if(verbose>0){
    std::cout << std::endl << "create rootfile with genbased reweighted distibutions" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << "input rootfile: " << analysisFileName << std::endl;
    std::cout << "parton level plots for the calculation of the distortion ";
    std::cout << "SFs are taken from folder: " << genfolder << std::endl;
    std::cout << "used tree for filling distorted reco ";
    std::cout << "level plots: " << folder+"/tree" << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;
    std::cout << "output rootfiles: " << outputFileNameUp << ", " << outputFileNameDown << std::endl;
    std::cout << "reweighted parton level distributions are stored in folder: " << genfolder << std::endl;
    std::cout << "reweighted reconstruction level distributions are stored in folder: " << folder << std::endl;
  }

  // define list of variations
  std::vector<TString> variation_;
  variation_.push_back("Up");
  variation_.push_back("Down");
  
  // ---
  //     Loading
  // ---
  // load input file
  TFile*inFile= new TFile(analysisFileName, "read");
  if(!inFile||inFile->IsZombie()){
    std::cout << "there seems to be a problem with the chosen input file " << analysisFileName << std::endl;
    exit(0);
  }
  // draw a copy of the input file
  // the cross section input plots will 
  // be replaced in part F
  if(save){
    if(verbose>0){
      std::cout << "copy old root file, will only replace shifted plots" << std::endl;
      std::cout << "a) variation up" << std::endl;
    }
    inFile->Cp(analysisFileName, outputFileNameUp);
    if(verbose>0) std::cout << "b) variation down" << std::endl;
    inFile->Cp(analysisFileName, outputFileNameDown);
  }
  // loading trees
  TTree *tree=(TTree*)(inFile->Get(folder+"/tree"));
  if(!tree||tree->IsZombie()){
     std::cout << "there seems to be a problem with the chosen tree " << folder+"/tree" << std::endl;
     exit(0);  
  }

  // initialize map entries with 0 
  // to avoid problems with the map re-ordering
  // when new entries are added
  value_["weight"]=0;
  if(verbose>0) std::cout << "considered variables: ";
  for(unsigned int i=0; i<variable_.size();++i){
    TH1F* hist=(TH1F*)(inFile->Get(folder+"/"+variable_[i]));
    TH2F* hist2D=(TH2F*)(inFile->Get(folder+"/"+variable_[i]+"_"));
    if(verbose>0) std::cout << variable_[i]+" ";
    value_[variable_[i]              ]=0;
    value_[variable_[i]+"PartonTruth"]=0;
    if(variable_[i].Contains("top")||variable_[i].Contains("bq")){
      value_[variable_[i]+"Lep"]=0;
      value_[variable_[i]+"Had"]=0;
      value_[variable_[i]+"LepPartonTruth"]=0;
      value_[variable_[i]+"HadPartonTruth"]=0;
    }
    // initialize plots with correct binning
    // original plots for parton level phase space / parton level extrapolated / reconstruction level / gen - rec correlation plots
    plots_      [variable_[i]+"PartonTruth"    ]=new TH1F( variable_[i]+"PartonTruth"    , variable_[i]+"PartonTruth"    , hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    plots_      [variable_[i]+"PartonTruthFull"]=new TH1F( variable_[i]+"PartonTruthFull", variable_[i]+"PartonTruthFull", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    plots_      [variable_[i]                  ]=new TH1F( variable_[i]                  , variable_[i]                  , hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    plots2D_    [variable_[i]                  ]=new TH2F( variable_[i]                  , variable_[i]                  , hist2D->GetNbinsX(), hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax(), hist2D->GetNbinsX(), hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());

    // loop all variations 
    for(unsigned int var=0; var<variation_.size();++var){
      TString Var=variation_[var];
      // systematic distorted plots for parton level phase space / reconstruction level
      // used for SF determination
      plotsScaled_[variable_[i]+"PartonTruth"+Var    ]=new TH1F( variable_[i]+"PartonTruth"+Var    , variable_[i]+"PartonTruth"+Var    , hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
      plotsScaled_[variable_[i]+"PartonTruthFull"+Var]=new TH1F( variable_[i]+"PartonTruthFull"+Var, variable_[i]+"PartonTruthFull"+Var, hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
      plotsScaled_[variable_[i]              +Var    ]=new TH1F( variable_[i]+Var                  , variable_[i]+Var                  , hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
      plotsScaled2D_[variable_[i]                    ]=new TH2F( variable_[i]+Var                  , variable_[i]+Var                  , hist2D->GetNbinsX(), hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax(), hist2D->GetNbinsX(), hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());
      for(unsigned int j=0; j<variable_.size();++j){
	// final plots: reweighted plots from tree entries
	// (reweighting based on all events on parton level within phase space)
	// reconstruction level plots and parton level plots for reconstructed events
	finalPlots_[variable_[j]+"PartonTruth"][variable_[i]+Var]=new TH1F( variable_[j]+"PartonTruth"+variable_[i]+"SF"+Var, variable_[j]+"PartonTruth"+variable_[i]+"SF"+Var, hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
	finalPlots_[variable_[j]              ][variable_[i]+Var]=new TH1F( variable_[j]              +variable_[i]+"SF"+Var, variable_[j]              +variable_[i]+"SF"+Var, hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
	finalPlots2D_[variable_[j]            ][variable_[i]+Var]=new TH2F( variable_[j]              +variable_[i]+"SF"+Var, variable_[j]              +variable_[i]+"SF"+Var, hist2D->GetNbinsX(), hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax(), hist2D->GetNbinsX(), hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());
      }
    }
  }
  if(verbose>0) std::cout << std::endl;

  // ---
  //    part A: get values for all variables from tree
  // ---
  if(verbose>0) std::cout << std::endl << "part A: get tree entries" << std::endl;
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("weight",1);
  tree->SetBranchAddress("weight",(&value_["weight"]));
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    // take care about splitted leptonic and hadronic top distributions
    if(variable_[i].Contains("top")||variable_[i].Contains("bq")){
      // activate branches
      tree->SetBranchStatus(variable_[i]+"Lep",1);
      tree->SetBranchStatus(variable_[i]+"Had",1);
      tree->SetBranchStatus(variable_[i]+"LepPartonTruth",1);
      tree->SetBranchStatus(variable_[i]+"HadPartonTruth",1);
      // save branch values in map
      tree->SetBranchAddress(variable_[i]+"Lep",(&value_[variable_[i]+"Lep"]));
      tree->SetBranchAddress(variable_[i]+"Had",(&value_[variable_[i]+"Had"]));
      tree->SetBranchAddress(variable_[i]+"LepPartonTruth",(&value_[variable_[i]+"LepPartonTruth"]));
      tree->SetBranchAddress(variable_[i]+"HadPartonTruth",(&value_[variable_[i]+"HadPartonTruth"]));
    }
    // same for ttbar / lepton quantities
    else{
      tree->SetBranchStatus(variable_[i],1);
      tree->SetBranchStatus(variable_[i]+"PartonTruth",1);
      tree->SetBranchAddress(variable_[i],(&value_[variable_[i]]));
      tree->SetBranchAddress(variable_[i]+"PartonTruth",(&value_[variable_[i]+"PartonTruth"]));
    }
  }
  // check if tree entries are ok
  std::vector<TString> uninitialized_;
  int Nevents=0;
  // loop all events
  for(unsigned int event=0; event<tree->GetEntries(); ++event){
    tree->GetEntry(event);
    Nevents++;
    // loop all variables
    for(unsigned int i=0; i<variable_.size();++i){
      // check leptonic/ hadronic top quantities
      if(variable_[i].Contains("top")||variable_[i].Contains("bq")){
	// leptonic top
	if(value_[variable_[i]+"Lep"]==-9999){
	  // check if already an invalid branch entry was found
	  vector<TString>::iterator it=find(uninitialized_.begin(), uninitialized_.end(), variable_[i]+"Lep");
	  if(uninitialized_.size()==0||it!=uninitialized_.end()) uninitialized_.push_back(variable_[i]+"Lep");
	}
	// hadronic top
	if(value_[variable_[i]+"Had"]==-9999){
	  // check if already an invalid branch entry was found
	  vector<TString>::iterator it=find(uninitialized_.begin(), uninitialized_.end(), variable_[i]+"Had");
	  if(uninitialized_.size()==0||it!=uninitialized_.end()) uninitialized_.push_back(variable_[i]+"Had");
	}
      }
      // check ttbar / lepton quantities
      else{ 
	if(value_[variable_[i]]==-9999){
	  // check if already an invalid branch entry was found
	  vector<TString>::iterator it=find(uninitialized_.begin(), uninitialized_.end(), variable_[i]);
	  if(uninitialized_.size()==0||it!=uninitialized_.end()) uninitialized_.push_back(variable_[i]);
	}
      }
    }
  }
  // print list of uninitialized tree entries
  if(uninitialized_.size()>0){
    std::cout << "there are " << uninitialized_.size();
    std::cout << " (partly) uninitialized branches in the tree:" << std::endl;
    for(unsigned int var=0; var<uninitialized_.size();++var){
      std::cout << uninitialized_[var] << std::endl;
    }
    exit(0);
  }
  else{
    if(verbose>0){
      std::cout << "all trees seem to be filled properly" << std::endl;
      std::cout << Nevents << " events in total" << std::endl;
    }
  }

  // ---
  //    part B: get weights for shape distortion
  // ---
  if(verbose>0) std::cout << std::endl << "part B: get weights for shape distortion" << std::endl;
  // B1 get original parton level distribution (phasespace) plots
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    TH1* targetPlot;
    TH1* targetPlot2;
    TH1* targetPlot3;
    inFile->GetObject(genfolder+"/"+variable_[i], targetPlot);
    inFile->GetObject(folder+"/"+variable_[i]   , targetPlot2);
    inFile->GetObject(genfolder2+"/"+variable_[i], targetPlot3);
    // check if plot exits
    if(!targetPlot){ 
      std::cout << "plot " << genfolder+"/"+variable_[i] << "can not ";
      std::cout << "be found in file " << analysisFileName << std::endl;
      exit(0);
    }
    else{
      plots_[variable_[i]+"PartonTruth"]=(TH1F*)(targetPlot->Clone());
      // check for empty plots
      if(plots_[variable_[i]+"PartonTruth"]->GetEntries()==0){
	std::cout << "ERROR: empty plot " << analysisFileName;
	std::cout << "/"+genfolder+"/"+variable_[i] << std::endl;
	exit(0);
      }
    }
    if(!targetPlot2){ 
      std::cout << "plot " << folder+"/"+variable_[i] << "can not ";
      std::cout << "be found in file " << analysisFileName << std::endl;
      exit(0);
    }
    else{
      plots_[variable_[i]]=(TH1F*)(targetPlot2->Clone());
      // check for empty plots
      if(plots_[variable_[i]]->GetEntries()==0){
	std::cout << "ERROR: empty plot " << analysisFileName;
	std::cout << "/"+folder+"/"+variable_[i] << std::endl;
	exit(0);
      }
    }
    if(!targetPlot3){ 
      std::cout << "plot " << genfolder2+"/"+variable_[i] << "can not ";
      std::cout << "be found in file " << analysisFileName << std::endl;
      exit(0);
    }
    else{
      plots_[variable_[i]+"PartonTruthFull"]=(TH1F*)(targetPlot3->Clone());
      // check for empty plots
      if(plots_[variable_[i]+"PartonTruthFull"]->GetEntries()==0){
	std::cout << "ERROR: empty plot " << analysisFileName;
	std::cout << "/"+genfolder2+"/"+variable_[i] << std::endl;
	exit(0);
      }
    }
  }
  // check if all plots are filled properly
  int NeventsPS=plots_["ttbarMassPartonTruth"]->GetEntries();
  if(verbose>0) std::cout << NeventsPS << " events in chosen Phase space" << std::endl;
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    // count number of entries
    int entries=plots_[variable_[i]+"PartonTruth"]->GetEntries();
    // take into account that there are two top quarks
    if(variable_[i].Contains("top")||variable_[i].Contains("bq")) entries*=0.5;
    // compare number of events
    if(entries!=NeventsPS){
	std::cout << "ERROR: wrong number of events in plot " << analysisFileName;
	std::cout << "/"+genfolder+"/"+variable_[i] << " or /ttbarMass" << std::endl;
	exit(0);
    }
  }
  if(verbose>0) std::cout << "original parton level phase space plots loaded successfully" << std::endl;

  // B2 distort distributions
  if(verbose>0) std::cout << "performing a fit and change parameter for Y(ttbar,top) and m(ttbar), otherwise a linear distortion is performed" << std::endl;
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    // loop all variations 
    for(unsigned int var=0; var<variation_.size();++var){
      TString Var=variation_[var];
      if (doPDFuncertainty) {
	// restricted phase space
	plotsScaled_[variable_[i]+"PartonTruth"    +Var]=distortPDF(*plots_[variable_[i]+"PartonTruth"    ], Var, variable_[i], inputFolderName, ""       , verbose-1);
	// full phase space
	TString pdfUncPS="Full";
	plotsScaled_[variable_[i]+"PartonTruthFull"+Var]=distortPDF(*plots_[variable_[i]+"PartonTruthFull"], Var, variable_[i], inputFolderName, pdfUncPS   , verbose-1);
      } 
      else {
	// restricted phase space
	plotsScaled_[variable_[i]+"PartonTruth"+Var    ]=distort   (*plots_[variable_[i]+"PartonTruth"], Var, variable_[i], verbose-1);
	// full phase space
	plotsScaled_[variable_[i]+"PartonTruthFull"+Var]=distort   (*plots_[variable_[i]+"PartonTruthFull"], Var, variable_[i], verbose-1);
      }
    }
  }
  
  // B3 calculate SF
  std::map< TString, TH1F* > SF_;

  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    // loop all variations 
    for(unsigned int var=0; var<variation_.size();++var){
      TString Var=variation_[var];
      // calculate lepton variation SF from full PS
      // others from restricted      // temporary necessary because of different gen and reco cuts
      TString PS="Full";
      //if(!(variable_[i]=="lepPt"||variable_[i]=="lepEta")) PS="";
      if(variable_[i].Contains("Y")) PS="";
      SF_[variable_[i]+Var]=(TH1F*)plotsScaled_[variable_[i]+"PartonTruth"+PS+Var]->Clone();
      SF_[variable_[i]+Var]->Divide(plots_[variable_[i]+"PartonTruth"+PS]);
      // if SF is derived by full PS
      if(PS=="Full"){
	for(int bin=0; bin<=SF_[variable_[i]+Var]->GetNbinsX()+1; ++bin){
	  //std::cout << "bin #" << bin << "(" << SF_[variable_[i]+Var]->GetBinCenter(bin) << "): " << SF_[variable_[i]+Var]->GetBinContent(bin) << std::endl;
	  // replace full PS parton truth plot by limited PS plot
	  // and apply weight from full PS
	  double SFbin=SF_[variable_[i]+Var]->GetBinContent(bin);
	  plotsScaled_[variable_[i]+"PartonTruth"+Var]->SetBinContent(bin, plots_[variable_[i]+"PartonTruth"]->GetBinContent(bin)*SFbin);
	}
      }
//      // Following lines are just to check whether reweighting values from MC@NLO to MadGraph works
//      if (doPDFuncertainty) {
//        // add central values (witouth PDF reweighting) from MC@NLO to the list of plots using MadGraph normalisation (see also ~line 462)
//        plots_[variable_[i]+"PartonTruth"+"MCatNLO"] = getTheoryPrediction(variable_[i], "/afs/naf.desy.de/group/cms/scratch/tophh/" + inputFolderName + "/ttbarNtupleCteq6mPDFuncertOnly.root");
//        plots_[variable_[i]+"PartonTruth"+"MCatNLO"]->Scale(plots_[variable_[i]+"PartonTruth"]->Integral(0,plots_[variable_[i]+"PartonTruth"]->GetNbinsX()+1));
//      }
    }
  }

  // B4 make sure to normalize the distorted plot 
  //    to the same area as the original one
  // -> dont change the absolut amount of events (in Phase space)!
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    // loop all variations 
    for(unsigned int var=0; var<variation_.size();++var){
      TString Var=variation_[var];
      double initArea = plots_[variable_[i]+"PartonTruth"]->Integral(0,plots_[variable_[i]+"PartonTruth"]->GetNbinsX()+1);
      double areaSFUp   = initArea/plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ]->Integral(0,plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ]->GetNbinsX()+1);
      double areaSFDown = initArea/plotsScaled_[variable_[i]+"PartonTruth"+"Down"]->Integral(0,plotsScaled_[variable_[i]+"PartonTruth"+"Down"]->GetNbinsX()+1);
      if (doPDFuncertainty) {
        if(verbose>0) std::cout << "NOTE: not doing any normalization for PDF uncertainties as this would underestimate the maximum variations (normalization should cancel out anyway in reco/parton level ratio, right?!)" << std::endl;
      } else {
        SF_[variable_[i]+"Up"  ]->Scale(areaSFUp);
        plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ]->Scale(areaSFUp);
        SF_[variable_[i]+"Down"]->Scale(areaSFDown);
        plotsScaled_[variable_[i]+"PartonTruth"+"Down"]->Scale(areaSFDown);
      }
      if(verbose>1) std::cout << "area central value    for "   << variable_[i] << ": " << initArea   << std::endl;
      if(verbose>1) std::cout << "area SF(up variation) for "   << variable_[i] << ": " << areaSFUp   << std::endl;
      if(verbose>1) std::cout << "area SF(down variation) for " << variable_[i] << ": " << areaSFDown << std::endl;
    }
  }
  if(verbose>0){
    std::cout << "distortion SF keeps normalization for original ";    
    std::cout << "(phase space) parton level truth distribution" << std::endl;
  }

  // B5 check for negative SF
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    // loop all variations 
    for(unsigned int var=0; var<variation_.size();++var){
      TString Var=variation_[var];
      // loop all bins
      for(int bin=0; bin<=SF_[variable_[i]+Var]->GetNbinsX()+1; ++bin){
	// search for negative SFs
	// but only in non-empty bins
        // NB: skip check for PDF uncertainties 
	if (!doPDFuncertainty&&SF_[variable_[i]+Var]->GetBinContent(bin)<0.&&plots_[variable_[i]+"PartonTruth"]->GetBinContent(bin)>0.){
	  std::cout << "ERROR: distortion SF for variable "+variable_[i]+" variation "+Var+" bin " << bin << " is negative!" << std::endl;
	  exit(1);
	}
      }
    }
  }
    
  // ---
  //    part C: Draw control plots for weighted (phase space) parton level distributions
  // ---
  if(verbose>0) std::cout << std::endl << "part C: create control plots" << std::endl;
  std::vector<TCanvas*> plotCanvas_;
  // C1 create canvas
  for(unsigned int plot=0; plot<variable_.size(); ++plot){
    char canvname[10];
    sprintf(canvname,"canv%i",plot);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    //canvasStyle(*plotCanvas_[sample]);
  }
  // C2 create legend
  TLegend *legPS = new TLegend(0.7, 0.6, 1.0, 0.85);
  legPS->SetFillStyle(0);
  legPS->SetBorderSize(0);
  legPS->SetEntrySeparation(0.45);
  legPS->SetHeader("#splitline{parton level}{(phase space)}");
  legPS->AddEntry(plotsScaled_[variable_[0]+"PartonTruth"+"Up"  ], "up"     , "L");
  legPS->AddEntry(plots_      [variable_[0]+"PartonTruth"       ], "central", "L");
  legPS->AddEntry(plotsScaled_[variable_[0]+"PartonTruth"+"Down"], "down"   , "L");

  // C3 fill scaled and original plot into one canvas
  int canvasNumber=0;
  // loop all variables
  for(unsigned int i=0; i<variable_.size();++i){
    plotCanvas_[canvasNumber]->cd(0);
    plotCanvas_[canvasNumber]->SetTitle(variable_[i]+"PartonLevelPhaseSpace");
    histogramStyle(*plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ], 0, false);
    histogramStyle(*plotsScaled_[variable_[i]+"PartonTruth"+"Down"], 2, false);
    histogramStyle(*plots_      [variable_[i]+"PartonTruth"       ], 0, false);
    plots_      [variable_[i]+"PartonTruth"       ]->SetLineColor(kBlack);
    plotsScaled_[variable_[i]+"PartonTruth"+"Up"]->GetXaxis()->SetNoExponent(true);
    plotsScaled_[variable_[i]+"PartonTruth"+"Up"]->GetXaxis()->SetTitle(variable_[i]);
    plotsScaled_[variable_[i]+"PartonTruth"+"Up"]->GetYaxis()->SetTitle("events");
    // take care of maximum
    double maxDown=plotsScaled_[variable_[i]+"PartonTruth"+"Down"]->GetMaximum();
    double maxUp  =plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ]->GetMaximum();
    double maxStd =plots_      [variable_[i]+"PartonTruth"       ]->GetMaximum();
    double maximum=std::max(std::max(maxUp,maxDown), maxStd);
    plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ]->SetMaximum(1.2*maximum);
    // adapt plot range
    if(variable_[i]=="lepPt") plotsScaled_[variable_[i]+"PartonTruth"+"Up"]->GetXaxis()->SetRangeUser(0,400);
    //draw plots
    plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ]->Draw("");
    plotsScaled_[variable_[i]+"PartonTruth"+"Down"]->Draw("same");
    plots_      [variable_[i]+"PartonTruth"       ]->Draw("same");
    // only for cross check
    //plots_      [variable_[i]+"PartonTruth"+"MCatNLO"]->Draw("same");
    legPS->Draw("same");
    ++canvasNumber;
  }
  //int partonLevelPlots=canvasNumber;

  // ---
  //    part D: create weighted plots
  // ---
  if(verbose>0) std::cout << std::endl << "part D: use parton level based SF to reweight reco histos" << std::endl;
  // loop all events
  for(unsigned int event=0; event<tree->GetEntries(); ++event){
    tree->GetEntry(event);
    double weight=value_["weight"];
    if(weight==0) std::cout << "WARNING: standard event weight (from tree) is 0!" << std::endl;
    // i: loop all variables (for SF)
    for(unsigned int i=0; i<variable_.size();++i){
      if(verbose>3) std::cout << "SF: " << variable_[i] << ", "; 
      // j: loop all variables (for plot to be scaled)
      for(unsigned int j=0; j<variable_.size();++j){
	if(verbose>3) std::cout << "variable: " << variable_[j] << std::endl;
	// loop all variations 
	for(unsigned int var=0; var<variation_.size();++var){
	  TString Var=variation_[var];
	  if(verbose>3) std::cout << "variation: " << Var << std::endl;
	  // get distortion SF (from variable i)
	  double LepShapeWeight=0;
	  double HadShapeWeight=0;
	  // for leptonic and hadronic top: get 2xSF (from variable i)
	  if(variable_[i].Contains("top")||variable_[i].Contains("bq")){
	    LepShapeWeight=SF_[variable_[i]+Var]->GetBinContent(SF_[variable_[i]+Var]->FindBin(value_[variable_[i]+"Lep"+"PartonTruth"]));
	    HadShapeWeight=SF_[variable_[i]+Var]->GetBinContent(SF_[variable_[i]+Var]->FindBin(value_[variable_[i]+"Had"+"PartonTruth"]));
	  }
	  // for other quantities: use one SF only (from variable i)
	  else{
	    LepShapeWeight=SF_[variable_[i]+Var]->GetBinContent(SF_[variable_[i]+Var]->FindBin(value_[variable_[i]+"PartonTruth"]));
	    HadShapeWeight=LepShapeWeight;
	  }
	  if(verbose>3) std::cout << "had weight: " << HadShapeWeight << std::endl;
	  if(verbose>3) std::cout << "lep weight: " << LepShapeWeight << std::endl;
	  if(verbose>3&&i==j&&variable_[i].Contains("Mass")&&HadShapeWeight==0.){
	    std::cout << "variable: " << variable_[i] << std::endl;
	    std::cout << "weight: " << HadShapeWeight << std::endl;
	    std::cout << "gen value " << value_[variable_[i]+"PartonTruth"] << std::endl;
	    std::cout << "rec value " << value_[variable_[i]] << std::endl;
	  }
	  // filling for top quantities: fill leptonic and hadronic top quantities in same plot
	  if(variable_[j].Contains("top")||variable_[i].Contains("bq")){
	    // parton truth plots
	    finalPlots_[variable_[j]+"PartonTruth"][variable_[i]+Var]->Fill(value_[variable_[j]+"Lep"+"PartonTruth"], weight*LepShapeWeight);
	    finalPlots_[variable_[j]+"PartonTruth"][variable_[i]+Var]->Fill(value_[variable_[j]+"Had"+"PartonTruth"], weight*HadShapeWeight);
	    if(verbose>3) std::cout << "lep parton truth: " << value_[variable_[j]+"Lep"+"PartonTruth"] << "*" <<  weight << "*" << LepShapeWeight << std::endl;
	    if(verbose>3) std::cout << "had parton truth: " << value_[variable_[j]+"Had"+"PartonTruth"] << "*" <<  weight << "*" << HadShapeWeight << std::endl;
	    // reconstruction level plots
	    finalPlots_[variable_[j]][variable_[i]+Var]->Fill(value_[variable_[j]+"Lep"], weight*LepShapeWeight);
	    finalPlots_[variable_[j]][variable_[i]+Var]->Fill(value_[variable_[j]+"Had"], weight*HadShapeWeight);
	    if(verbose>3) std::cout << "lep reco: " << value_[variable_[j]+"Lep"] << "*" <<  weight << "*" << LepShapeWeight << std::endl;
	    if(verbose>3) std::cout << "had reco: " << value_[variable_[j]+"Had"] << "*" <<  weight << "*" << HadShapeWeight << std::endl;
	    // correlation plots
	    finalPlots2D_[variable_[j]][variable_[i]+Var]->Fill(value_[variable_[i]+"Lep"+"PartonTruth"], value_[variable_[j]+"Lep"], weight*LepShapeWeight);
	    finalPlots2D_[variable_[j]][variable_[i]+Var]->Fill(value_[variable_[i]+"Had"+"PartonTruth"], value_[variable_[j]+"Had"], weight*HadShapeWeight);
	  }
	  // other quantities: take average of leptonic and hadronic top SF
	  else{
	    // parton truth plots
	    finalPlots_[variable_[j]+"PartonTruth"][variable_[i]+Var]->Fill(value_[variable_[j]+"PartonTruth"], weight*0.5*(LepShapeWeight+HadShapeWeight));
	    if(verbose>3) std::cout << "parton level: " << value_[variable_[j]+"PartonTruth"] << "*" <<  weight << "*" << 0.5*(LepShapeWeight+HadShapeWeight) << std::endl;
	    // reconstruction level plots
	    finalPlots_[variable_[j]][variable_[i]+Var]->Fill(value_[variable_[j]], weight*0.5*(LepShapeWeight+HadShapeWeight));
	    if(verbose>3) std::cout << "reco: " << value_[variable_[j]] << "*" <<  weight << "*" << 0.5*(LepShapeWeight+HadShapeWeight) << std::endl;
	    // correlation plots
	    finalPlots2D_[variable_[j]][variable_[i]+Var]->Fill(value_[variable_[j]+"PartonTruth"], value_[variable_[j]], weight*0.5*(LepShapeWeight+HadShapeWeight));
	  }
	}
      }
    }
  }

  // ---
  //    part E: get BG subtracted data shape plots
  // ---
  // define data reference plot list
  std::vector<TString> plotList_;
  std::vector<TString> axisLabel_;
  for(unsigned int plot=0; plot<variable_.size(); ++plot){
    plotList_.push_back( folder+"/"+variable_[plot] );
    TString rebinFactor="1";
    if(variable_[plot]=="ttbarMass") rebinFactor="20";
    if(variable_[plot]=="ttbarPt"  ) rebinFactor="6";
    if(variable_[plot]=="topPt"    ) rebinFactor="10";
    if(variable_[plot]=="lepPt"    ) rebinFactor="10";
    axisLabel_.push_back(variable_[plot]+" from data/events/0/"+rebinFactor);
    if(verbose>2) std::cout << variable_[plot]+" from data/events/0/"+rebinFactor << std::endl;
    //    plotList_.push_back( genfolder+"/"+variable_[i] );
  }
  unsigned int N1Dplots = plotList_.size();

  // open standard analysis files
  std::map<unsigned int, TFile*> files_ = getStdTopAnalysisFiles(path+inputFolderName, 0, dataFile, decayChannel);
  // define container for data plots
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2_;
  // total # plots 
  int Nplots=0;
  std::vector<TString> vecRedundantPartOfNameInData;
  getAllPlots(files_, plotList_, histo_, histo2_, N1Dplots, Nplots, verbose-1, decayChannel, &vecRedundantPartOfNameInData);
  // apply lumiweighting
  scaleByLuminosity(plotList_, histo_, histo2_, N1Dplots, luminosity, verbose-1, 0, decayChannel);
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    // apply rebinning
    TString plotName=plotList_[plot];
    double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4)).c_str());
    equalReBinTH1(reBinFactor, histo_, plotName, kBkg  );
    equalReBinTH1(reBinFactor, histo_, plotName, kSTop );
    equalReBinTH1(reBinFactor, histo_, plotName, kWjets);
    equalReBinTH1(reBinFactor, histo_, plotName, kZjets);
    equalReBinTH1(reBinFactor, histo_, plotName, kDiBos);
    equalReBinTH1(reBinFactor, histo_, plotName, kQCD  );
    equalReBinTH1(reBinFactor, histo_, plotName, kData );
    // subtract non ttbar prompt lepton BG from data
    histo_[plotName][kData]->Add(histo_[plotName][kBkg  ],-1);
    histo_[plotName][kData]->Add(histo_[plotName][kSTop ],-1);
    histo_[plotName][kData]->Add(histo_[plotName][kWjets],-1);
    histo_[plotName][kData]->Add(histo_[plotName][kZjets],-1);
    histo_[plotName][kData]->Add(histo_[plotName][kDiBos],-1);
    //histo_[plotName][kData]->Add(histo_[plotName][kQCD  ],-1); // neglect QCD for the moment as its spiky and unphysical
    // normalize data to ttbar MC reco yield
    double areaMC = plots_[variable_[plot]]->Integral(0,plots_[variable_[plot]]->GetNbinsX()+1);
    double areaData = histo_[plotName][kData]->Integral(0,histo_[plotName][kData]->GetNbinsX()+1);
    histo_[plotList_[plot]][kData]->Scale(areaMC/areaData);
    // adapt style
    histogramStyle( *histo_[plotName][kData], kData, true);
  }
  
  // ---
  //    part F: draw all reweighted reco level plots
  // ---
  TLegend* legReco= new TLegend(0.65, 0.55, 1.2, 0.85);
  legReco->SetFillStyle(0);
  legReco->SetBorderSize(0);
  legReco->SetEntrySeparation(0.25);
  legReco->SetHeader("reconstruction level");
  legReco->AddEntry(finalPlots_[variable_[0]][variable_[0]+"Up"  ], "up"     , "L");
  legReco->AddEntry(plots_      [variable_[0]]                    , "central", "L");
  legReco->AddEntry(finalPlots_[variable_[0]][variable_[0]+"Down"], "down"   , "L");
  legReco->AddEntry(histo_[plotList_[0]][kData]                   , "data"   ,"PL");
  if(verbose>0) std::cout << std::endl << "part E: collect all final plots" << std::endl;
  // i: loop all variables (for SF)
  for(unsigned int i=0; i<variable_.size();++i){
    // j: loop all variables (for shown plot)
    for(unsigned int j=0; j<variable_.size();++j){
      // fill original/scaled reco plots only into canvas   
      if(i==j){
	// create canvas for every plot
	char canvname[10];
	sprintf(canvname,"canv%i",canvasNumber);
	plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
	plotCanvas_[canvasNumber]->cd(0);
	// set canvas title and headline
	TString title = finalPlots_[variable_[j]][variable_[i]+"Up"]->GetTitle();
	TString PS2="RecoLevel";
	if(title.Contains("PartonTruth")) PS2="PartonLevel";
	finalPlots_[variable_[j]][variable_[i]+"Up"  ]->SetTitle(variable_[j]+PS2+"SFfrom"+variable_[i]+"Up"  );
	finalPlots_[variable_[j]][variable_[i]+"Down"]->SetTitle(variable_[j]+PS2+"SFfrom"+variable_[i]+"Down");
	plotCanvas_[canvasNumber]->SetTitle(variable_[j]+PS2+"SFfrom"+variable_[i]);
	// set style
	histogramStyle(*finalPlots_[variable_[j]][variable_[i]+"Up"  ], 0, false);
	histogramStyle(*finalPlots_[variable_[j]][variable_[i]+"Down"], 2, false);
	histogramStyle(*plots_[variable_[j]], 0, false);
	plots_[variable_[j]]->SetLineColor(kGreen);
	finalPlots_[variable_[j]][variable_[i]+"Up"]->GetXaxis()->SetNoExponent(true);
	finalPlots_[variable_[j]][variable_[i]+"Up"]->GetXaxis()->SetTitle(variable_[j]);
	finalPlots_[variable_[j]][variable_[i]+"Up"]->GetYaxis()->SetTitle("events");
	// clone histos
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]=(TH1F*)finalPlots_[variable_[j]][variable_[i]+"Up"  ]->Clone();
	finalPlots_[variable_[j]][variable_[i]+"DownRebinned"]=(TH1F*)finalPlots_[variable_[j]][variable_[i]+"Down"]->Clone();
	plots_[variable_[j]+"Rebinned"]=(TH1F*)plots_[variable_[j]]->Clone();
	// rebin MC histos (only for this plot, NOT for the rootfile) 
	double reBinFactor = atof(((string)getStringEntry(axisLabel_[i],4)).c_str());
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]->Rebin(reBinFactor);
	finalPlots_[variable_[j]][variable_[i]+"DownRebinned"]->Rebin(reBinFactor);
	plots_[variable_[j]+"Rebinned"]->Rebin(reBinFactor);
	// take care of maximum
	double maxDown=finalPlots_[variable_[j]][variable_[i]+"DownRebinned"]->GetMaximum();
	double maxUp  =finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]->GetMaximum();
	double maxData=0;
	if(i==j&&PS2=="RecoLevel") maxData=histo_[plotList_[i]][kData]->GetMaximum();
	double maxStd =plots_[variable_[j]+"Rebinned"]->GetMaximum();
	double maximum=std::max(std::max(maxUp,maxDown), std::max(maxData, maxStd));
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"]->SetMaximum(1.2*maximum);
	// remaining style issues
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]->GetXaxis()->SetNoExponent(true);
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]->GetXaxis()->SetTitle(variable_[i]);
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]->GetYaxis()->SetTitle("events");
	// adapt plot range
	if(variable_[j]=="lepPt") finalPlots_[variable_[j]][variable_[i]+"UpRebinned"]->GetXaxis()->SetRangeUser(0,400);
	// draw histos
	finalPlots_[variable_[j]][variable_[i]+"UpRebinned"  ]->DrawClone("");
	finalPlots_[variable_[j]][variable_[i]+"DownRebinned"]->DrawClone("same");
	plots_[variable_[j]+"Rebinned"]->DrawClone("same");
	if(i==j&&PS2=="RecoLevel") histo_[plotList_[i]][kData]->DrawClone("p e1 same");
	// change name
	finalPlots_[variable_[j]][variable_[i]+"Down"]->SetName (variable_[j]);
	finalPlots_[variable_[j]][variable_[i]+"Down"]->SetTitle(variable_[j]);
	finalPlots_[variable_[j]][variable_[i]+"Up"  ]->SetName (variable_[j]);
	finalPlots_[variable_[j]][variable_[i]+"Up"  ]->SetTitle(variable_[j]);
	finalPlots2D_[variable_[j]][variable_[i]+"Up"  ]->SetName (variable_[j]+"_");
	finalPlots2D_[variable_[j]][variable_[i]+"Up"  ]->SetTitle(variable_[j]+"_");
	// draw legend
	if(PS2=="RecoLevel") legReco->Draw("same");
	// adapt counter
	++canvasNumber;
      }
    }
  }
  
  // ---
  //    part G: save weighted plots and control plots
  // ---
  bool saveIntoRootfile=true;
  if(save){
    if(verbose>0){
      std::cout << std::endl << "part F1: save rweighted parton level and";
      std::cout << "reconstruction level plots to output folder" << std::endl;
    }
    // F1: save control plots
    for(unsigned int plot=0; plot<plotCanvas_.size(); ++plot){
      plotCanvas_[plot]->Print(savePlotsTo+"/"+plotCanvas_[plot]->GetTitle()+".png");
      plotCanvas_[plot]->Print(savePlotsTo+"/"+plotCanvas_[plot]->GetTitle()+".eps");
    }
    // F2: save scaled parton level distribution 
    // and reco plots in root file for systematic variation
    if(verbose>0) std::cout << std::endl << "part F2: save reweighted reco histos to output files" << std::endl;
    // i: loop all variables (for SF)
    for(unsigned int i=0; i<variable_.size();++i){
      // j: loop all variables (for shown plot)
      for(unsigned int j=0; j<variable_.size();++j){
	// up to now: save only plots where reco plot is weighted with its own gen distribution
	if(saveIntoRootfile&&i==j){
	  // parton level plot
	  // attention: the distorted PS parton level plots used for SF determination are used
	  saveToRootFile(outputFileNameUp  , plotsScaled_[variable_[i]+"PartonTruth"+"Up"  ], true, verbose-1, genfolder);
	  saveToRootFile(outputFileNameDown, plotsScaled_[variable_[i]+"PartonTruth"+"Down"], true, verbose-1, genfolder);
	  // reco plot (weighted with distorted parton level distribution of the same quantity)
	  saveToRootFile(outputFileNameUp  , finalPlots_[variable_[j]][variable_[i]+"Up"  ] , true, verbose-1, folder);
	  saveToRootFile(outputFileNameDown, finalPlots_[variable_[j]][variable_[i]+"Down"] , true, verbose-1, folder);
	  // correlation plot (weighted with distorted parton level distribution of the same quantity)
	  saveToRootFile(outputFileNameUp  , finalPlots2D_[variable_[j]][variable_[i]+"Up"  ] , true, verbose-1, folder);
	  saveToRootFile(outputFileNameDown, finalPlots2D_[variable_[j]][variable_[i]+"Down"] , true, verbose-1, folder);

	}
      }
    }
  }

  // restore gErrorIgnoreLv
  gErrorIgnoreLevel=initialgErrorLv;
}

TH1F* distort(const TH1& hist, TString variation, TString variable, int verbose)
{
  // this function does a linear shape distortion.
  // for variation=="up"/"down" bins that are close in 
  // x to the x value with maximum bin content(xmax)
  // are scaled more/less
  if(verbose>1) std::cout << variable << " (" << variation << ")" << std::endl;
  TH1F* result=(TH1F*)(hist.Clone());
  TH1F* histo =(TH1F*)(hist.Clone());
  if(variable.Contains("lepEta")){
    // fit gaussian
    //    TF1* myFunc = new TF1("myFunc", "gaus", xmin, xmax);
    TF1* myFunc = new TF1("myFunc", "[0]*TMath::Gaus(x, 0, [1])", -5.0, 5.0);
    myFunc->SetParameter(0, 18000); 
    myFunc->SetParameter(1, 1.3); 
    myFunc->SetParLimits(0, 100, 100000);
    myFunc->SetParLimits(1, 0.5 , 10);
    histo->Fit(myFunc,"Q","same");
    // do variation
    // vary shape 
    TF1* myVarFunc=(TF1*)(myFunc->Clone("myVarFunc"));
    if(     variation=="Down") myVarFunc->SetParameter(1, 0.85*myVarFunc->GetParameter(1));
    else if(variation=="Up"  ) myVarFunc->SetParameter(1, 1.5*myVarFunc->GetParameter(1));
    // get SF
    TH1F* fit = (TH1F*)result->Clone();
    fit->Scale(0.);
    fit->Add(myFunc);
    TH1F* fitVar = (TH1F*)result->Clone();
    fitVar->Scale(0.);
    fitVar->Add(myVarFunc);
    TH1F* SFVar=(TH1F*)fitVar->Clone();
    SFVar->Divide(fit);
    // apply SF to original plot
    result->Multiply(SFVar);
    return result;
  }
  if(variable.Contains("ttbarMass")){
    double xmin=200;
    double xmax=1200;
    TF1* myFunc = new TF1("myFunc", "landau", xmin, xmax);
    myFunc->SetParameter(0, 2890); 
    myFunc->SetParameter(1, 420); 
    myFunc->SetParameter(2, 34); 
    myFunc->SetParLimits(0, 1000, 3500);
    myFunc->SetParLimits(1, 300 , 500 );
    myFunc->SetParLimits(2, 10, 80);
    if(verbose>1){
      std::cout << "" << std::endl;
      std::cout << "fitting landau to " << variable << std::endl;
      std::cout << "---------------------------------------" << std::endl;
    }
    //histo->Fit(myFunc,"","same",xmin, xmax);
    histo->Fit(myFunc,"Q","same");
    // vary shape 
    TF1* myVarFunc=(TF1*)(myFunc->Clone("myVarFunc"));
    if(variation=="Down"){
      myVarFunc->SetParameter(0, 0.98*myFunc->GetParameter(0)); 
      myVarFunc->SetParameter(1, 1.02*myFunc->GetParameter(1)); 
      myVarFunc->SetParameter(2, 1.1*myFunc->GetParameter(2));
    }
    else if(variation=="Up"){
      myVarFunc->SetParameter(0, 1.05*myFunc->GetParameter(0)); 
      myVarFunc->SetParameter(1, 0.95*myFunc->GetParameter(1)); 
      myVarFunc->SetParameter(2, 0.8*myFunc->GetParameter(2));
    }
    // derive SF
    TH1F* fit = (TH1F*)result->Clone();
    fit->Scale(0.);
    fit->Add(myFunc);
    TH1F* fitVar = (TH1F*)result->Clone();
    fitVar->Scale(0.);
    fitVar->Add(myVarFunc);
    TH1F* SFVar=(TH1F*)fitVar->Clone();
    SFVar->Divide(fit);
    // apply SF to original plot
    result->Multiply(SFVar);
    // additional scaling for the tail
    for(int bin=1; bin<=result->GetNbinsX(); ++bin){
      if(result->GetBinCenter(bin)>=420){
	double scalefactor=myVarFunc->GetParameter(1)/myFunc->GetParameter(1);
	//	std::cout << variable << "(" << variation << "): " << scalefactor << std::endl;
	result->SetBinContent(bin,result->GetBinContent(bin)*scalefactor*scalefactor);
	if(result->GetBinCenter(bin)>=450) result->SetBinContent(bin,result->GetBinContent(bin)*scalefactor*scalefactor);
	if(result->GetBinCenter(bin)>=650) result->SetBinContent(bin,result->GetBinContent(bin)*scalefactor*scalefactor);
	if(result->GetBinCenter(bin)>=700) result->SetBinContent(bin,result->GetBinContent(bin)*scalefactor*scalefactor);
      }
    }
    return result;
  }
  else if(variable.Contains("Y")){
    // create empty original histogram to have binning at hand
    result->Scale(0.);
    // print out info
    if(verbose>1){
      std::cout << "" << std::endl;
      std::cout << "fitting [0]*TMath::Cos(x*[1])+[2]*TMath::Gaus(x, 0, [3]) to " << variable << std::endl;
      std::cout << "---------------------------------------" << std::endl;
    }
    double xmin=0;
    double xmax=0;
    if(variable.Contains("top")){
      xmin=-2.5;
      xmax=2.5;
    }
    else if(variable.Contains("ttbar")){
      xmin=-2.0;
      xmax=2.0;
    }
    else{
      std::cout << "ERROR: unknown quantity " << variable << std::endl;
      exit(0);
    }
    TF1* myFunc = new TF1("myFunc","[0]*TMath::Cos(x*[1])+[2]*TMath::Gaus(x, 0, [3])", xmin, xmax);
    if(variable.Contains("top")){
      myFunc->SetParameter(0, 4.34485e+03); 
      myFunc->SetParameter(1, 1.33482e+00); 
      myFunc->SetParameter(2, 4.35032e+03); 
      myFunc->SetParameter(3, 9.59915e+01);
      myFunc->SetParLimits(0, 2000, 10000);
      myFunc->SetParLimits(1, 1   , 1.6  );
      myFunc->SetParLimits(2, 3000, 10000);
      myFunc->SetParLimits(3, 0.5 , 1000 ); 
    }
    else if(variable.Contains("ttbar")){
      myFunc->SetParameter(0, 2.02216e+03); 
      myFunc->SetParameter(1, 1.76526e+00); 
      myFunc->SetParameter(2, 3.28315e+03); 
      myFunc->SetParameter(3, 1.92234e+00); 
      myFunc->SetParLimits(0, 0   , 10000);
      myFunc->SetParLimits(1, 0.5 , 2.0  );
      myFunc->SetParLimits(2, 0   , 10000);
      myFunc->SetParLimits(3, 1   , 1000 );
    }
    // do fit for hist in range [xmin,xmax]
    histo->Fit(myFunc,"Q","same",xmin, xmax);
    // vary shape 
    double scaleFactor=0;
    if(variation=="Up") scaleFactor=0.9;
    else if(variation=="Down") scaleFactor=2;
    else{ 
      std::cout << "ERROR: unknown systematic variation" << variation << std::endl;
      exit(0);
    }
    if(variation=="Down"){
      // vary frequency of cosinus function
      myFunc->SetParameter(1, myFunc->GetParameter(1)/(0.6*scaleFactor));
      // vary normalization of gaussian
      myFunc->SetParameter(2, scaleFactor*scaleFactor*scaleFactor*scaleFactor*myFunc->GetParameter(2));
      // vary broughtness of gaussian parameter
      myFunc->SetParameter(3, scaleFactor*scaleFactor*scaleFactor*myFunc->GetParameter(3));
    }
    // vary broughtness of gaussian parameter
    else{
      if(variable.Contains("ttbar")) myFunc->SetParameter(3,0.8*myFunc->GetParameter(3));
      else myFunc->SetParameter(3, scaleFactor*0.7*myFunc->GetParameter(3));
    }
    // separate gaussian and cosinus part
    TF1* mygaus = new TF1("mygaus","[0]*TMath::Gaus(x, 0, [1])",xmin, xmax);
    TF1* mycos  = new TF1("mycos" ,"[0]*TMath::Cos(x*[1])",xmin, xmax);
    mygaus->SetParameter(0, myFunc->GetParameter(2));
    mygaus->SetParameter(1, myFunc->GetParameter(3));
    //    mygaus->SetRange(scaleFactor*xmin, scaleFactor*xmax);
    mycos ->SetParameter(0, myFunc->GetParameter(0));
    mycos ->SetParameter(1, myFunc->GetParameter(1));
    // apply chosen binning to obtained histogram
    result->Add(mycos);
    result->Add(mygaus);
    // correct for offset <0
    TF1* myOffset  = new TF1("myOffset","[0]",scaleFactor*xmin, scaleFactor*xmax);
    // for more narrow distribution: draw only part for y above 0
    if(myFunc->GetMinimum()<0){
      double min=std::abs(result->GetXaxis()->GetBinCenter(result->FindFirstBinAbove(0,1)));
      if(verbose>1){
	std::cout << "min ,0 (findbin, lowedge, min)" << std::endl;
	std::cout << result->FindFirstBinAbove(0,1) << std::endl;
	std::cout << result->GetXaxis()->GetBinCenter(result->FindFirstBinAbove(0,1)) << std::endl;
	std::cout << min << std::endl;
      }
      mygaus->SetRange(-1.*min,min);
      mycos->SetRange(-1.*min,min);
      result->Scale(0.);
      result->Add(mycos);
      result->Add(mygaus);
    }
    // for more brought distribution: remove offset, draw only up to minimum
    else if(myFunc->GetMinimum()!=0){
      double offset=myFunc->GetMinimum();
      double min=std::abs(myFunc->GetX(offset));
      myOffset->SetParameter(0, -1*offset);
      myOffset->SetRange(-1.*min,min);
      if(verbose>1){
	std::cout << "min ,0 (offset, xvalue)" << std::endl;
	std::cout << offset << std::endl;
	std::cout << min << std::endl;
      }
      mygaus->SetRange(-1.*min,min);
      mycos ->SetRange(-1.*min,min);
      result->Scale(0.);
      result->Add(mycos);
      result->Add(mygaus);
      result->Add(myOffset);
    }
    // edit color of fit and extrapolate to whole region of x
    
    //  std::cout << "a = " << myPol->GetParameter(1) << " +/- " << myPol->GetParError(1) << std::endl;
    //std::cout << "b = " << myPol->GetParameter(0) << " +/- " << myPol->GetParError(0) << std::endl;
    //std::cout << "chi2/ndof = " << myPol->GetChisquare() / myPol->GetNDF() << std::endl;
    //std::cout << "probability = " << TMath::Prob(myPol->GetChisquare(),myPol->GetNDF()) << std::endl;

    
    // top Y: [0]*TMath::Cos(x*[1])+[2]*TMath::Gaus(x, 0, [3])
    // range: -2.5 - 2,5
    //   parameter,   value,        uncertainty    boundary
    //1  p0           4.34485e+03   1.81315e+01    2000 - 10000
    //2  p1           1.33482e+00   1.53763e-03    1    - 1.6
    //3  p2           4.35032e+03   1.24700e+01    3000 - 10000
    //4  p3           9.59915e+01   2.29207e+02    0.5  - 1000


    // ttbar Y: [0]*TMath::Cos(x*[1])+[2]*TMath::Gaus(x, 0, [3])
    // range: -2.0 - 2.0
    //   parameter,   value,        uncertainty    boundary
    //1  p0           2.02216e+03   2.79619e+02    0-10000
    //2  p1           1.76526e+00   6.07716e-02    0.5-2
    //3  p2           3.28315e+03   2.92333e+02    0-10000
    //4  p3           1.92234e+00   4.89337e-01    1-1000
  }
  // linear reweighting for other quantities
  else{
    // parameter of linear smearing
    double a=0;
    if     (variation=="Up"  ) a=1.0;
    else if(variation=="Down") a=-1.0;
    else{
      std::cout << "ERROR in function distort: chose variation " << variation << " not known" << std::endl;
      exit(0);
    }   
    // adjust amount of smearing
    // depending on reco level data shape
    if(variable=="topPt"){
      if     (variation=="Up"  ) a*=0.5;
      else if(variation=="Down") a*=10.0;
    }
    else if(variable=="ttbarPt"  ){ 
      if     (variation=="Up"  ) a*=0.5;
      else if(variation=="Down") a*=10.0;
    }
    else if(variable=="lepEta"  ){
      if     (variation=="Up"  ) a*=1.5;
      else if(variation=="Down") a*=4.0;
    }
    else if(variable=="lepPt"   ){
      if     (variation=="Up"  ) a*=0.1;
      else if(variation=="Down") a*=10.0;
    }
    else{
     std::cout << "ERROR in function distort: chose variable " << variable << " not known" << std::endl;
      exit(0);
    }
    // SF control variables
    double SFmin=10000;
    double SFmax=0;
    // find reference point so that 
    // weight is always positive 0 for up 
    // and x_max for down variation
    double xmax= 0.;
    if(variation=="Down") xmax=std::abs(hist.GetBinCenter(hist.GetNbinsX())+1);
    // loop bins
    for(int bin=1; bin<result->GetNbinsX(); ++bin){ 
      // get SF 
      double x= std::abs(hist.GetBinCenter(bin));
      double SF=linSF(x, xmax, a, 0.);
      if(variable=="lepPt"){ 
	if(a<0) SF=pow(SF, 1.4);
	else SF=pow(SF, 0.4);
      }
      else if(variable=="topPt"&&variation=="Up"){
	SF=pow(SF, 0.7);	
      }
     else if(variable=="ttbarPt"){
       if(variation=="Up") SF=pow(SF, 0.7);
       else SF=pow(SF, 0.6);
      }
      // avoid SF=0 
      if(SF==0){
	std::cout << "WARNING: SF for bin " << bin << " is 0" << std::endl;
      	double xPrevious= hist.GetBinCenter(bin-1);
	double xNext=xPrevious;
      	if(xPrevious!=hist.GetNbinsX()) xNext=hist.GetBinCenter(bin+1);
      	// take average SF of surrounding bins
      	SF= 0.5*(linSF(x, xNext, a, 0)+linSF(x, xPrevious, a, 0));
      }
      // search max/min SF
      if(verbose>1){
	if(SF>SFmax) SFmax=SF;
	if(SF<SFmin) SFmin=SF;
      }
      // scale bin content
      result->SetBinContent(bin,SF*hist.GetBinContent(bin));
    }
    // print max/min SF
    if(verbose>1){
      std::cout << "minimum SF: "<< SFmin << std::endl;
      std::cout << "maximum SF: "<< SFmax << std::endl;
    }
  }
  // return distorted plot
  if(verbose>1) std::cout << "return result" << std::endl;
  return result;
}

double linSF(const double x, const double xmax, const double a, const double b){
  double SF=(x-xmax)*a+b;
  if(SF<0){
    std::cout << "ERROR in function SF: negative result!" << std::endl;
    std::cout << "(" << a << "*" << "[" << x << "-" << xmax << "])" << "+" << b << "=" << SF << std::endl;
    exit(0);
  }
  if(a<0) SF=pow(SF, 5.);
  else  SF=pow(SF, 0.5);
  return SF;
}

TH1F* distortPDF(const TH1& hist, TString variation, TString variable, TString inputFolderName, TString phaseSpace, int verbose)
{
  if (verbose>0) std::cout << "Using hard-coded path '"+groupSpace+"CommonFiles' instead of '"+groupSpace+"'" << inputFolderName << " for file with PDF-uncertainties." << std::endl;
  // this function loads the max/min PDF uncertainties as determined externally with MC@NLO for the desired variables and applies it to MadGraph
  TString fileName    = groupSpace+"CommonFiles/ttbarNtupleCteq6mPDFuncertOnly.root";
  TString plotNameVar = (variation == "") ? variable : variable + "_" + variation                                              ;
  TString plotNameNom = variable                                                                                               ;                 
  TString directory   = (phaseSpace == "Full") ? "FullPhaseSpace" : "VisiblePhaseSpace"                                        ; 
  if(verbose>0) std::cout << "distortPDF: using plot " << directory + "/" + plotNameNom << " from file " << fileName << " for estimated PDF uncertainties in plot " << plotNameVar << std::endl;

  TH1F * histPdfNom = getTheoryPrediction(directory + "/" + plotNameNom,fileName) ;  // the central values of MC@NLO (needed to rescale to central values of MadGraph)
  TH1F * histPdfVar = getTheoryPrediction(directory + "/" + plotNameVar,fileName) ;  // the varied  values of MC@NLO (needed to add     to central values of MadGraph)
  TH1F * histPdf    = getTheoryPrediction(directory + "/" + plotNameVar,fileName) ;  // the desired output histogram (PDF variations for MadGraph)

  // check if binning and ranges are consistent
  if ( histPdf->GetNbinsX() != hist.GetNbinsX() || histPdf->GetXaxis()->GetXmin() != hist.GetXaxis()->GetXmin() || histPdf->GetXaxis()->GetXmax() != hist.GetXaxis()->GetXmax() )
  {
    std::cout << "ERROR in distortPDF function: histograms do not have the same binning/ranges" << std::endl ;
    exit(1);
  }

  // scale each bin such that it represents the PDF uncertainty found with MC@NLO but with MadGraph as the default central value
  for (int i=0 ; i <= hist.GetNbinsX()+1; ++i) {
    histPdf->SetBinContent(i,0.) ;
    if (histPdfNom->GetBinContent(i) != 0.) histPdf->SetBinContent(i, histPdfVar->GetBinContent(i) / histPdfNom->GetBinContent(i) * hist.GetBinContent(i) ) ; 
  }

  return histPdf ;

}
