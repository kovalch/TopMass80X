#include "basicFunctions.h"

std::map<TString, std::vector<double> > makeVariableBinningA(int binning=1);
bool PythiaSample (unsigned int s);
bool ATLASSample  (unsigned int s);

// sample list
unsigned int kMad    =0;
unsigned int kPow    =1;
unsigned int kPowHer =2;
unsigned int kMca    =3;
unsigned int kAlp    =4;
unsigned int kPowA   =5;
unsigned int kPowA2  =6;
unsigned int kPowHerA=7;
unsigned int kMcaA   =8;

void ATLASCompTreeSGsamples(bool save = true, int verbose=1, int binning=0, TString outputfolder="./diffXSecFromSignal/plots/combined/2012/comparisonATLAS/"){
//void ATLASCompTreeSGsamples(bool save = true, int verbose=1, int binning=0, TString outputfolder="./ttgencomparison/"){
  // binning= 0:fine binning, 1:ATLAS, 2:CMS
  bool debug  = verbose>0 ? true : false;
  bool debug2 = verbose>1 ? true : false;
  bool excludeATLAS=false;

  // ============================
  //  documentation on how to run
  // ----------------------------
  // run via root -q -b -l ATLASCompTreeSGsamples.C++g
  // a) parameters
  // - choose binning via "binning"= 0:fine binning, 1:ATLAS, 2:CMS
  // - choose output level via "verbose"= 0: minimal, 1: detailed, 2: debug
  // - choose via "save" whether you want to save single plots as eps and in rootfile and all plots in one pdf
  // - choose destination where plots are save via "outputfolder"
  // - change the binning in "makeVariableBinningA"
  // ============================


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
  //  load rootfiles
  // ============================
  std::vector<TFile* > file_;
  //file_.push_back(TFile::Open("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun/combinedDiffXSecSigFall11PFLarge.root"                                        , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/work/i/iasincru/public/TopLHCWG_DiffXSex_CMS/MC_theory_samples/CMSttbarMadGraphZ2Pythia6CTEQ6L1Fall11MCProductionCycle.root", "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/work/i/iasincru/public/TopLHCWG_DiffXSex_CMS/MC_theory_samples/CMSttbarPowhegZ2Pythia6CTEQ6MFall11MCProductionCycle.root"   , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/work/i/iasincru/public/TopLHCWG_DiffXSex_CMS/MC_theory_samples/CMSttbarPowhegAUET2Herwig6CTEQ6MFall11MCProductionCycle.root", "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/work/i/iasincru/public/TopLHCWG_DiffXSex_CMS/MC_theory_samples/CMSttbarMC@NLOHerwig6CTEQ6MFall11MCProductionCycle.root"     , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/user/d/disipio/public/toplhcwg/ntuples_atlas/AlpgenJimmyttbarlnqq.root"                                                     , "Open"));
  //file_.push_back(TFile::Open("/afs/cern.ch/user/d/disipio/public/toplhcwg/ntuples_atlas/AlpGenPythia_P2011_CTEQ5L_ttbarlnqq.root"                                      , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/user/d/disipio/public/toplhcwg/ntuples_atlas/TTbar_PowHeg_Pythia_P2011C.root"                                               , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/user/d/disipio/public/toplhcwg/ntuples_atlas/TTbar_PowHeg_Pythia_AUET2.root"                                                , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/user/d/disipio/public/toplhcwg/ntuples_atlas/TTbar_PowHeg_Jimmy.root"                                                       , "Open"));
  file_.push_back(TFile::Open("/afs/cern.ch/user/d/disipio/public/toplhcwg/ntuples_atlas/T1_McAtNlo_Jimmy.root"                                                         , "Open"));

  // list plots of relevance
  std::vector<TString> plotList_, axisLabel_;
  TString plots1D[ ] = {
    // KinFit plots before prob cut 
    "topPt",
    "topY",
    "ttbarPt",
    "ttbarY",
    "ttbarMass",
    //"decayChannel",
  };
  TString axisLabel1D[ ] = { 
    // KinFit plots before prob cut 
    "p_{T}^{t} #left[GeV#right];Top quarks (norm.);0;20",
    "y^{t};Top quarks (norm.);0;1",
    "p_{T}^{t#bar{t}} #left[GeV#right];t#bar{t} pairs (norm.);0;20",
    "y^{t#bar{t}};t#bar{t} pairs (norm.);0;1",
    "m^{t#bar{t}} #left[GeV#right];t#bar{t} pairs (norm.);0;50",
    //"t#bar{t} decay Channel;relative #Top-quark pairs;0;1",
  };

  plotList_ .insert(plotList_ .begin(), plots1D    , plots1D    + sizeof(plots1D    )/sizeof(TString));
  axisLabel_.insert(axisLabel_.begin(), axisLabel1D, axisLabel1D+ sizeof(axisLabel1D)/sizeof(TString));
  if(plotList_.size() != axisLabel_.size()){
    std::cout << "ERROR - 1D plots: Number of plots and axis label do not correspond .... Exiting macro!" << std::endl;
    exit(1);
  }
  // run automatically in batch mode if there are many canvas
  if(plotList_.size()>15) gROOT->SetBatch();

  // create canvas container
  std::vector<TCanvas*> plotCanvas_, plotCanvas2_;
  // create legend
  TLegend* leg = new TLegend(0.0, 0.0, 1.0, 1.0);
  legendStyle(*leg ,TString("#bf{t#bar{t} simulation, ")+(PHD ? "CMS" : "LHC")+", #sqrt{s}=7 TeV}",(PHD ? 0.05 : 0.03)  );
  TLegend* leg2= new TLegend(0.4, 0.6, 0.85, 0.88);
  legendStyle(*leg2,"#bf{t#bar{t} PYTHIA simulation, #sqrt{s}=7 TeV}");
  TLegend* leg3= new TLegend(0.4, 0.6, 0.85, 0.88);
  legendStyle(*leg3,"#bf{t#bar{t} HERWIG simulation, #sqrt{s}=7 TeV}");
  TLegend* leg4= new TLegend(0.4, 0.7, 0.85, 0.88);
  legendStyle(*leg4,"#bf{t#bar{t} default simulation, #sqrt{s}=7 TeV}");
  TLegend* leg5= new TLegend(0.4, 0.7, 0.85, 0.88);
  legendStyle(*leg5,"#bf{t#bar{t} Powheg simulation, #sqrt{s}=7 TeV}");

  // ============================
  //  get histos from tree
  // ============================
  unsigned int kfirst    =kMad;
  unsigned int klast     =kMcaA;
  unsigned int krelative1 =kPow; 
  TString krelative1lab=excludeATLAS||PHD ? "#scale[0.85]{Powheg+Pythia}" : "#scale[0.85]{#splitline{Powheg+Pythia}{     (CMS)}}";
  unsigned int krelative2 =kPow;
  TString krelative2lab=krelative1lab;
  unsigned int krelative3 =kPowHer;
  TString krelative3lab=excludeATLAS||PHD  ? "#scale[0.85]{Powheg+Herwig}" : "#scale[0.85]{#splitline{Powheg+Herwig}{     (CMS)}}";  
  unsigned int krelative4 =kMad;
  TString krelative4lab="#scale[0.55]{MadGraph+Pythia (CMS)}";   
  unsigned int krelative5 =kPowHer;
  TString krelative5lab=krelative3lab;
  TString treePath ="genTree/tree";
  TString treePathA="tree";
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map<TString, std::vector<double> > binning_=makeVariableBinningA(binning);
  std::map< unsigned int, double> Ntotev_;

  // get template histos
  std::vector<TH1F*> template_, template2_, template3_, template4_, template5_;
  if(debug) std::cout << "get template histos" << std::endl; 
  // loop all plots
  for(int plot=0; plot<(int)plotList_.size(); ++plot){
    TString thname="analyzeTopPartonLevelKinematics/"+plotList_[plot];
    if(debug) std::cout << thname << std::endl;
    //TH1F* temp=new TH1F(plotList_[plot], plotList_[plot], 100000, -5000., 5000.);
    TH1F* temp=new TH1F(plotList_[plot], plotList_[plot], binning_[plotList_[plot]].size()-1, &(binning_[plotList_[plot]][0]));
    //double rebinFactor =atof(((string)getStringEntry(axisLabel_[plot],4,";")).c_str());
    //temp->Rebin(rebinFactor);
    //reBinTH1F(*temp, binning_[plotList_[plot]], 2);
    temp->Reset("icms");
    temp->SetTitle("");
    temp->GetXaxis()->SetTitle(getStringEntry(axisLabel_[plot],1,";"));
    temp->GetYaxis()->SetTitle(getStringEntry(axisLabel_[plot],2,";"));
    temp->GetXaxis()->SetNoExponent(true);
    temp->SetStats(kFALSE);
    temp->SetLineWidth(3);
    temp->SetMarkerSize(1.25);
    //int binMax=temp->GetNbinsX()+1;
    template_ .push_back(temp);
    template2_.push_back(temp);
    template3_.push_back(temp);
    template4_.push_back(temp);
    template5_.push_back(temp);
  }
  if(debug) std::cout << "process trees" << std::endl; 
  // loop all samples
  for(unsigned int sample=kfirst; sample<=klast; ++sample){
    bool CMS = TString(file_.at(sample)->GetName()).Contains("CMS");
    if(debug){
      std::cout << "sample " << sample;
      if(CMS)  std::cout << " (CMS)";
      else  std::cout << " (ATLAS)";
      std::cout << std::endl; 
    }
    // get tree
    TString treePath2=treePath;
    if(!CMS)treePath2=treePathA;
    TTree* tree = (TTree*)(file_[sample]->Get(treePath2));
    if(!tree){
      std::cout << "WARNING: tree " << treePath << " not found in sample " << sample << ", will continue and neglect this one " << std::endl;
      // exit(0);
    }
    // container for values read from tree
    std::map< TString, float  > value_;
    std::map< TString, float > valueA_;
    // initialize map entries with 0 
    value_["weight"   ]=1;
    value_["topPt"    ]=0;
    value_["topbarPt" ]=0;
    value_["topY"     ]=0;
    value_["topbarY"  ]=0;
    value_["ttbarPt"  ]=0;
    value_["ttbarY"   ]=0;
    value_["ttbarMass"]=0;
    //value_["decayChannel"]=0;
    valueA_["weight"   ]=1;
    valueA_["topPt"    ]=0;
    valueA_["topbarPt" ]=0;
    valueA_["topY"     ]=0;
    valueA_["topbarY"  ]=0;
    valueA_["ttbarPt"  ]=0;
    valueA_["ttbarY"   ]=0;
    valueA_["ttbarMass"]=0;
    //valueA_["decayChannel"]=0;
    // initialize branches
    if(tree){
      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("weight"      ,1);
      tree->SetBranchStatus("topPt"       ,1);
      tree->SetBranchStatus("topbarPt"    ,1);
      tree->SetBranchStatus("topY"        ,1);
      tree->SetBranchStatus("topbarY"     ,1);
      tree->SetBranchStatus("ttbarPt"     ,1);
      tree->SetBranchStatus("ttbarY"      ,1);
      tree->SetBranchStatus("ttbarMass"   ,1);
      //tree->SetBranchStatus("decayChannel",1);
      if(CMS){
	tree->SetBranchAddress("weight"      , (&value_["weight"   ]   ));
	tree->SetBranchAddress("topPt"       , (&value_["topPt"    ]   ));
	tree->SetBranchAddress("topbarPt"    , (&value_["topbarPt" ]   ));
	tree->SetBranchAddress("topY"        , (&value_["topY"     ]   ));
	tree->SetBranchAddress("topbarY"     , (&value_["topbarY"  ]   ));
	tree->SetBranchAddress("ttbarPt"     , (&value_["ttbarPt"  ]   ));
	tree->SetBranchAddress("ttbarY"      , (&value_["ttbarY"   ]   ));
	tree->SetBranchAddress("ttbarMass"   , (&value_["ttbarMass"]   ));
	//tree->SetBranchAddress("decayChannel", (&value_["decayChannel"]));
      }
      else{
	tree->SetBranchAddress("weight"      , (&valueA_["weight"   ]   ));
	tree->SetBranchAddress("topPt"       , (&valueA_["topPt"    ]   ));
	tree->SetBranchAddress("topbarPt"    , (&valueA_["topbarPt" ]   ));
	tree->SetBranchAddress("topY"        , (&valueA_["topY"     ]   ));
	tree->SetBranchAddress("topbarY"     , (&valueA_["topbarY"  ]   ));
	tree->SetBranchAddress("ttbarPt"     , (&valueA_["ttbarPt"  ]   ));
	tree->SetBranchAddress("ttbarY"      , (&valueA_["ttbarY"   ]   ));
	tree->SetBranchAddress("ttbarMass"   , (&valueA_["ttbarMass"]   ));
	//tree->SetBranchAddress("decayChannel", (&valueA_["decayChannel"]));
      }
    }

    // initialize histo
    if(debug) std::cout << "initialize histos from template" << std::endl; 
    int color =kRed+7; TString sampleName="MadGraph+Pythia old";
    //if(     sample==kMadOld){ color =kRed+7  ; sampleName="MadGraph+Pythia old" ;}
    if(sample==kMad||sample==kAlp){ color =kRed    ; sample==kMad ? sampleName="MadGraph+Pythia" : sampleName="Alpgen+Herwig"     ; sample==kMad ? sampleName+="(Z2)" : sampleName+="(AUET2)"     ;}
    else if(sample==kPow||sample==kPowA||sample==kPowA2){ 
      color =kGreen  ; sampleName="Powheg+Pythia";
      if(     sample==kPow  ) sampleName+="(Z2)";
      else if(sample==kPowA ) sampleName+="(P11)";
      else if(sample==kPowA2){ sampleName+="(AUET2)"; color=kGreen+3;}
    }
    else if(sample==kPowHer||sample==kPowHerA){ color =kMagenta; sampleName="Powheg+Herwig(AUET2)";} 
    else if(sample==kMca   ||sample==kMcaA   ){ color =kBlue   ; sampleName="MC@NLO+Herwig"; if(sample==kMcaA) sampleName+="(AUET2)";}
    if(CMS&&!excludeATLAS&&!PHD) sampleName+="(CMS)"  ;
    else if(!CMS){
      sampleName+="(ATLAS)";
      //color+=1;      
    }
    for(int plot=0; plot<(int)plotList_.size(); ++plot){
      // initialize result plots from template
      histo_[plotList_[plot]][sample]=(TH1F*)template_[plot]->Clone(plotList_[plot]+getTStringFromInt(sample));
      histo_[plotList_[plot]][sample]->SetLineColor(color);
      histo_[plotList_[plot]][sample]->SetMarkerColor(color);
      if(!CMS){
	histo_[plotList_[plot]][sample]->SetLineStyle(2);
	if(sample==kPowA2) histo_[plotList_[plot]][sample]->SetLineStyle(3);
	 histo_[plotList_[plot]][sample]->SetLineWidth(4.5);
      }
    }
    // Add legend entry
    if(CMS||!excludeATLAS){
      std::cout << "sample: " << sample << std::endl;
      std::cout << "CMS? " <<  CMS << std::endl;
      std::cout << "excludeATLAS? " <<  excludeATLAS << std::endl;
      if(!PHD||!ATLASSample(sample)){
	leg->AddEntry(histo_[plotList_[0]][sample], sampleName, "L");
	if( PythiaSample(sample)) leg2->AddEntry(histo_[plotList_[0]][sample], sampleName, "L");
	if(!PythiaSample(sample)) leg3->AddEntry(histo_[plotList_[0]][sample], sampleName, "L");
	if(sample==kMad||sample==kAlp) leg4->AddEntry(histo_[plotList_[0]][sample], sampleName, "L");
	if(sample==kPow||sample==kPowA||sample==kPowA2||sample==kPowHer||sample==kPowHerA) leg5->AddEntry(histo_[plotList_[0]][sample], sampleName, "L");
      }
    }
    // loop tree
    if(tree&&(CMS||!excludeATLAS)){
      if(debug) std::cout << "fill plots from tree" << std::endl; 
      for(unsigned int event=0; event<tree->GetEntries(); ++event){
	// get event
	tree->GetEntry(event);
	// get relevant quantities
	float weight       = CMS ? value_["weight"      ] : valueA_["weight"      ];
	//float decayChannel = CMS ? value_["decayChannel"] : valueA_["decayChannel"];
	float topPtLep     = CMS ? value_["topPt"       ] : valueA_["topPt"       ];
	float topPtHad     = CMS ? value_["topbarPt"    ] : valueA_["topbarPt"    ];
	float topYLep      = CMS ? value_["topY"        ] : valueA_["topY"        ];
	float topYHad      = CMS ? value_["topbarY"     ] : valueA_["topbarY"     ];
	float ttbarPt      = CMS ? value_["ttbarPt"     ] : valueA_["ttbarPt"     ];
	float ttbarY       = CMS ? value_["ttbarY"      ] : valueA_["ttbarY"      ];
	float ttbarMass    = CMS ? value_["ttbarMass"   ] : valueA_["ttbarMass"   ];
	// debugging output
	if(debug2){
	  std::cout << "event " << event+1 << "/" << tree->GetEntries() << std::endl;
	  std::cout << "weight:       " << weight       << std::endl;
	  //std::cout << "decayChannel: " << decayChannel << std::endl;
	  std::cout << "topPtLep:     " << topPtLep     << std::endl;
	  std::cout << "topPtHad:     " << topPtHad     << std::endl;
	  std::cout << "topYLep:      " << topYLep      << std::endl;
	  std::cout << "topYHad:      " << topYHad      << std::endl;
	  std::cout << "ttbarPt:      " << ttbarPt      << std::endl;
	  std::cout << "ttbarMass:    " << ttbarMass    << std::endl;
	  std::cout << "ttbarY:       " << ttbarY       << std::endl;
	}
	// fill histo for all
	histo_["topPt"     ][sample]->Fill(topPtLep , weight);
	histo_["topPt"     ][sample]->Fill(topPtHad , weight);
	histo_["topY"      ][sample]->Fill(topYLep  , weight);
	histo_["topY"      ][sample]->Fill(topYHad  , weight);
	histo_["ttbarPt"   ][sample]->Fill(ttbarPt  , weight);
	histo_["ttbarY"    ][sample]->Fill(ttbarY   , weight);
	histo_["ttbarMass" ][sample]->Fill(ttbarMass, weight);
      } // end loop event
    } // end if tree
    // save N(events)
    if(sample==krelative1||sample==krelative2||sample==krelative3||sample==krelative4||sample==krelative5) Ntotev_[sample]=histo_[plotList_[0]][sample]->Integral(0.,histo_[plotList_[0]][sample]->GetNbinsX()+1);
    for(int plot=0; plot<(int)plotList_.size(); ++plot){
      // normalize to unity
      histo_[plotList_[plot]][sample]->Scale(1./histo_[plotList_[plot]][sample]->Integral(0.,histo_[plotList_[plot]][sample]->GetNbinsX()+1)); 
    }
  } // end loop samples


  // create label
  TPaveText *label = new TPaveText();
  label -> SetX1NDC(gStyle->GetPadLeftMargin());
  label -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label -> SetY2NDC(1.0);
  label -> SetTextFont(42);
  label -> AddText(excludeATLAS ? "CMS simulation, #sqrt{s} = 7 TeV" : (PHD ?  "private work, simulation at #sqrt{s} = 7 TeV" : "TOPLHCWG Preliminary, #sqrt{s} = 7 TeV"));
  label -> SetFillStyle(0);
  label -> SetBorderSize(0);
  label -> SetTextSize(0.04);
  label -> SetTextAlign(32);

  // ============================
  //  create canvas
  // ============================
  if(debug) std::cout << "create canvas" << std::endl; 
  for(int set=1; set<=5; ++set){
    // loop plots
    for(int plot=0; plot<(int)plotList_.size(); ++plot){
      TString name=plotList_[plot];
      TH1F* temptemplate=0;
      if(set==1) temptemplate=(TH1F*)template_ [plot]->Clone(TString(template_ [plot]->GetName())+"1");
      if(set==2) temptemplate=(TH1F*)template2_[plot]->Clone(TString(template2_[plot]->GetName())+"2");
      if(set==3) temptemplate=(TH1F*)template3_[plot]->Clone(TString(template3_[plot]->GetName())+"3");
      if(set==4) temptemplate=(TH1F*)template4_[plot]->Clone(TString(template4_[plot]->GetName())+"4");
      if(set==5) temptemplate=(TH1F*)template5_[plot]->Clone(TString(template4_[plot]->GetName())+"5");
      TString nameExt= set==2 ? "PYTHIA" : (set==3 ? "HERWIG" : (set==4 ? "MadgraphAlpgen" : ( set==5 ? "Powheg" : "")));
      if(debug) std::cout << "plot " << name << std::endl; 
      addCanvas(plotCanvas_);
      //######################################
      if(plotList_[plot].Contains("topPt")||plotList_[plot].Contains("ttbarPt")){
        plotCanvas_[plotCanvas_.size()-1]->SetLogy(1);
        temptemplate->GetYaxis()->SetRangeUser(0.0001, 1000);
      }
      //######################################
      plotCanvas_[plotCanvas_.size()-1]->cd(0);
      plotCanvas_[plotCanvas_.size()-1]->SetName(name+nameExt);
      plotCanvas_[plotCanvas_.size()-1]->SetTitle(name+nameExt);
      temptemplate->SetMaximum(1.5*histo_[name][kPow]->GetMaximum());
      if(plotCanvas_[plotCanvas_.size()-1]->GetLogy()) temptemplate->SetMaximum(10*temptemplate->GetMaximum());
      if(plotList_[plot].Contains("Y")) temptemplate->SetMaximum(1.2*temptemplate->GetMaximum());
      if(plotList_[plot].Contains("ttbarMass")) temptemplate->SetLabelSize(0.03);
      //temptemplate->GetXaxis()->SetLabelSize(0);
      //temptemplate->GetXaxis()->SetTitleSize(0);
      temptemplate->Draw("AXIS");
      // draw all samples
      for(unsigned int sample=kfirst; sample<=klast; ++sample){
	if(histo_[name].count(sample)>0&&(set==1||(set==2&&PythiaSample(sample))||(set==3&&!PythiaSample(sample))||(set==4&&(sample==kMad||sample==kAlp))||(set==5&&(sample==kPow||sample==kPowA||sample==kPowA2||sample==kPowHer||sample==kPowHerA)))&&(!excludeATLAS||!ATLASSample(sample))){
	  if(debug) std::cout << " - draw sample " << sample << std::endl; 
	  if(!PHD||!ATLASSample(sample)) histo_[name][sample]->Draw("hist same"); // exclude ATLAS samples for PHD
	}
      }
      // legend
      TLegend* templeg=0;
      // - clone correct legend
      if     (set==1) templeg=(TLegend*)(leg ->Clone(TString("templeg")+plotList_[plot]+getTStringFromInt(set)));
      else if(set==2) templeg=(TLegend*)(leg2->Clone(TString("templeg")+plotList_[plot]+getTStringFromInt(set)));
      else if(set==3) templeg=(TLegend*)(leg3->Clone(TString("templeg")+plotList_[plot]+getTStringFromInt(set)));
      else if(set==4) templeg=(TLegend*)(leg4->Clone(TString("templeg")+plotList_[plot]+getTStringFromInt(set)));
      else if(set==5) templeg=(TLegend*)(leg5->Clone(TString("templeg")+plotList_[plot]+getTStringFromInt(set)));
      // - adjust legend position for different quantities
      if(plotList_[plot].Contains("topPt")){
	templeg->SetX1(0.22);
	templeg->SetY1(0.40);
	templeg->SetX2(0.67);
	templeg->SetY2(0.61);
      }
      else if(plotList_[plot].Contains("topY")){
	templeg->SetX1(0.40);
	templeg->SetY1(0.37);
	templeg->SetX2(0.85);
	templeg->SetY2(0.55);
      }
      else if(plotList_[plot].Contains("ttbarPt")){	
	templeg->SetX1(0.46);
	templeg->SetY1(0.63);
	templeg->SetX2(0.91);
	templeg->SetY2(0.87);
      }
      else if(plotList_[plot].Contains("ttbarY")){	
	templeg->SetX1(0.45);
	templeg->SetY1(0.68);
	templeg->SetX2(0.90);
	templeg->SetY2(0.88);
      }
      else if(plotList_[plot].Contains("ttbarMass")){	
	templeg->SetX1(0.44);
	templeg->SetY1(0.60);
	templeg->SetX2(0.89);
	templeg->SetY2(0.86);
      }
      // - draw it
      if(debug) std::cout << " - draw legend" << std::endl; 
      if(!PHD) templeg->Draw("same");
      label->Draw("same");
      if(debug) std::cout << " - draw label" << std::endl; 
      // zero error
      std::vector<double> zeroerr_;
      for(int bin=0; bin<temptemplate->GetNbinsX(); ++bin) zeroerr_.push_back(0);
      // draw ratios
      if(debug) std::cout << " - draw ratios:" << std::endl; 
      //bool first=true;
      unsigned int krelative=krelative1;
      TString krelativelab=krelative1lab;
      TString nominatorLabel= "simulation";
      if(set==2 ){ krelative=krelative2; krelativelab=krelative2lab;}
      if(set==3 ){ krelative=krelative3; krelativelab=krelative3lab;}
      if(set==4 ){ krelative=krelative4; krelativelab=krelative4lab; nominatorLabel= "#scale[0.55]{ATLAS Alpgen+Herwig}";}
      if(set==5 ){ krelative=krelative5; krelativelab=krelative5lab;}
      // ratio y axis min/max valuse
      double min=0.3;
      double max=1.7;
      if(binning>0){min=0.7 ; max=1.3 ;}
      if(set==4)   {min=0.85; max=1.15;}
      if(binning==0&&(name.Contains("topY")||name.Contains("ttbarY"))) {min=0.80; max=1.15;}
      if(binning==0&&(name.Contains("topPt"    ))) {min=0.75; max=1.25;}
      if(binning==0&&(name.Contains("ttbarPt"  ))) {min=0.40; max=1.15;}
      if(binning==0&&(name.Contains("ttbarMass"))) {min=0.40; max=1.35;}
      // axis and labels for ratio plot
      drawRatio(temptemplate, temptemplate, min, max, myStyle, verbose, zeroerr_, krelativelab, nominatorLabel, "AXIS", kWhite);
      // draw errorband for relative MC
      if(histo_[name].count(krelative)>0){
	std::vector<double> err_;
	for(int bin=0; bin<histo_[name][krelative]->GetNbinsX(); ++bin){
	  err_.push_back(sqrt(histo_[name][krelative]->GetBinContent(bin)*1000000)); // NOTE: 1M events are assumed for the statistical error at the moment
	}
	TH1F* ratiotemp =(TH1F*)(histo_[name][krelative]->Clone(TString(histo_[name][krelative]->GetName())+"errup"));
	TH1F* ratiotemp2=(TH1F*)(histo_[name][krelative]->Clone(TString(histo_[name][krelative]->GetName())+"errdn"));
	TH1F* ratiotemp3=(TH1F*)(histo_[name][krelative]->Clone(TString(histo_[name][krelative]->GetName())+"errc" ));
	if(debug) std::cout << "draw uncertainty bands" << std::endl;
	for(int bin=0; bin<=histo_[name][krelative]->GetNbinsX()+1; ++bin){
	  if(debug2) std::cout << "bin: #" << bin << " - ";
	  double Ntotev=Ntotev_[krelative]; 
	  double relUnc=1. / ( sqrt(histo_[name][krelative]->GetBinContent(bin)*Ntotev) );
	  // take care of empty bins
	  if(histo_[name][krelative]->GetBinContent(bin)==0.){
	    relUnc=max-1.0;
	    ratiotemp ->SetBinContent(bin, 0.000001);
	    ratiotemp2->SetBinContent(bin, 0.000001);
	    ratiotemp3->SetBinContent(bin, 0.000001);
	  }
	  ratiotemp ->SetBinContent(bin, (1.+relUnc)*ratiotemp ->GetBinContent(bin)); 	  
	  ratiotemp2->SetBinContent(bin, (1.-relUnc)*ratiotemp2->GetBinContent(bin)); 
	  if(debug2) std::cout << "content: " << histo_[name][krelative]->GetBinContent(bin) << ", unc: " << relUnc << std::endl;
	}
	int ratioColor =kGray;
	int whiteColor=10;
	ratiotemp ->SetFillStyle(1001);
	ratiotemp2->SetFillStyle(1001);
	ratiotemp ->SetLineWidth(1);
	ratiotemp2->SetLineWidth(1);
	ratiotemp ->SetFillColor(ratioColor);
	ratiotemp2->SetFillColor(whiteColor);
	ratiotemp ->SetLineColor(ratioColor);
	ratiotemp2->SetLineColor(whiteColor);
	// central value+1sigma statistics filled in gray
	drawRatio(ratiotemp   , ratiotemp3, min, max, myStyle, verbose-1, zeroerr_, nominatorLabel, krelativelab, "hist same", ratioColor, false, 0.1);
	// central value+1sigma statistics filled in white
	drawRatio(ratiotemp2  , ratiotemp3, min, max, myStyle, verbose-1, zeroerr_, nominatorLabel, krelativelab, "hist same", whiteColor, false, 0.1);
	// redraw axis
	drawRatio(temptemplate, temptemplate, min, max, myStyle, verbose-1, zeroerr_, nominatorLabel, krelativelab, "AXIS same", kWhite);
      }
      for(unsigned int sample=kfirst; sample<=klast; ++sample){
	if((histo_[name].count(sample)>0&&histo_[name].count(krelative)>0)&&((set==1)||(set==2&&PythiaSample(sample))||(set==3&&!PythiaSample(sample))||(set==4&&(sample==kAlp||sample==kMad))||(set==5&&(sample==kPow||sample==kPowA||sample==kPowA2||sample==kPowHer||sample==kPowHerA)))&&(!excludeATLAS||!ATLASSample(sample))){
	  if(debug) std::cout << "   sample " << sample << " / sample " << krelative << std::endl; 
	  //TString opt = first ? "hist" : "hist same";
	  TString opt = "hist same";
	  //if(first) first=false;
	  if(!PHD||!ATLASSample(sample)) drawRatio(histo_[name][sample], histo_[name][krelative], min, max, myStyle, verbose, zeroerr_, nominatorLabel, krelativelab, opt, histo_[name][sample]->GetLineColor());
	}
      }
      //if(set>1) DrawLabel(TString(template_[plot]->GetXaxis()->GetTitle()), 0.2, 0.07, 0.95, 0.3, 32, 0.15);
      temptemplate->Draw("AXIS same");
    } // end for loop plot
  } // end for loop sep

  // legend
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetName ("Legend");
  plotCanvas_[plotCanvas_.size()-1]->SetTitle("Legend");
  leg->Draw("");
  
  // saving
  if(save){
    TString name=outputfolder+"treeATLASCMScomparison";
    if(excludeATLAS) name+="CMSonly";
    TString name2=binning==1 ? "ATLASbinning" : (binning==2 ? "CMSbinning" : "FineBinning");
    if(debug) std::cout << "save plots as pictures" << std::endl;  
    saveCanvas(plotCanvas_, name+name2, "", true, true, true);
    if(debug) std::cout << "save plots in rootfile" << std::endl;
    for(unsigned int i=0; i<plotCanvas_.size(); ++i){
      if(debug) {
	std::cout << i+1 << "/" << plotCanvas_.size();
	std::cout << ": "<< plotCanvas_[i]->GetTitle();
	std::cout << "->" << name+".root" << " (subfolder " << name2 << ")" << std::endl;
      }
      saveToRootFile(name+".root", plotCanvas_[i], true, 0, name2);
    }
  }
}


std::map<TString, std::vector<double> > makeVariableBinningA(int binning)
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
  double topPtBinsA[]= {0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 350.0, 800.};
  double topPtBins[] = {0.0, 60.0, 100.0, 150.0, 200.0, 260.0, 320.0, 400.};
  if(binning==2) bins_.insert( bins_.begin(), topPtBins, topPtBins + sizeof(topPtBins)/sizeof(double) );
  else if(binning ==1) bins_.insert( bins_.begin(), topPtBinsA, topPtBinsA + sizeof(topPtBinsA)/sizeof(double) );
  else{
    double width=20.;
    double min=0.;
    double max=500.;
    for(double bin=min; bin<=max; bin+=width){
      bins_.push_back(bin);
    }
  }
  result["topPt"        ] = bins_;
  bins_.clear();

  // y(top)
  double topYBins[]={-2.5, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.5}; 
  if(binning>0) bins_.insert( bins_.begin(), topYBins, topYBins + sizeof(topYBins)/sizeof(double) );
  else{
    double width=0.1;
    double min=-2.5;
    double max= 2.5;
    for(double bin=min; bin<=max; bin+=width){
      bins_.push_back(bin);
    }
  }
  result["topY"   ]   = bins_;
  bins_.clear();

  // pt(ttbar)
  double ttbarPtBinsA[] = {0.0, 40.0, 170.0, 340.0, 1000.0};
  double ttbarPtBins[]  = {0.0, 20.0, 45.0 , 75.0 , 120.0 , 190.0, 300.0};
  if(binning==2) bins_.insert( bins_.begin(), ttbarPtBins, ttbarPtBins + sizeof(ttbarPtBins)/sizeof(double) );
  else if(binning ==1) bins_.insert( bins_.begin(), ttbarPtBinsA, ttbarPtBinsA + sizeof(ttbarPtBinsA)/sizeof(double) );
  else{
    double width=20.;
    double min=0.;
    double max=500.;
    for(double bin=min; bin<=max; bin+=width){
      bins_.push_back(bin);
    }
  }
  result["ttbarPt"]=bins_;
  bins_.clear();

  // y(ttbar)
  double ttbarYBinsA[] = {-2.5, -1.0, -0.5, 0.0, 0.5, 1.0, 2.5};
  double ttbarYBins[]  = {-2.5, -1.3, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.3, 2.5};
  if(binning==2) bins_.insert( bins_.begin(), ttbarYBins, ttbarYBins + sizeof(ttbarYBins)/sizeof(double) );
  else if(binning ==1) bins_.insert( bins_.begin(), ttbarYBinsA, ttbarYBinsA + sizeof(ttbarYBinsA)/sizeof(double) );
  else{
    double width=0.1;
    double min=-2.5;
    double max= 2.5;
    for(double bin=min; bin<=max; bin+=width){
      bins_.push_back(bin);
    }
  }
  result["ttbarY"]=bins_;
  bins_.clear();

  // m(ttbar)
  double ttbarMassBinsA[] = {345.0, 450.0, 550.0, 700.0, 950.0, 2700.0};
  double ttbarMassBins[]  = {345.0, 400.0, 470.0, 550.0, 650.0, 800.0, 1100.0, 1600.0};
  if(binning==2) bins_.insert( bins_.begin(), ttbarMassBins, ttbarMassBins + sizeof(ttbarMassBins)/sizeof(double) );
  else if(binning ==1) bins_.insert( bins_.begin(), ttbarMassBinsA, ttbarMassBinsA + sizeof(ttbarMassBinsA)/sizeof(double) );
  else{
    double width=20;
    double min=250.;
    double max=2000.;
    for(double bin=min; bin<=max; bin+=width){
      if(bin>1000) bin+=width;
      if(bin>1600) bin+=width;
      bins_.push_back(bin);
    }
  }
  result["ttbarMass"]=bins_;
  bins_.clear();

  return result;
}

bool PythiaSample (unsigned int s)
{
  if(s==kMad  ) return true;
  if(s==kPow  ) return true;
  if(s==kPowA ) return true;
  if(s==kPowA2) return true;
  return false;
}

bool ATLASSample (unsigned int s)
{
  if(s==kAlp    ) return true;
  if(s==kPowA   ) return true;
  if(s==kPowA2  ) return true;
  if(s==kPowHerA) return true;
  if(s==kMcaA   ) return true;
 return false;
}
