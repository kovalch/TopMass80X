#include "basicFunctions.h"

void analyzeFileComparison(bool save = true,  bool usePAPERbinning=true, int verbose=0){
 
  // usePAPERbinning: use paper binning instead of equidistant binning for cross section quantities
  // usePAPERbinning=true / false
  // ============================
  //  Set Root Style
  // ============================
		
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.SetStripDecimals(true);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);

  // ============================
  //  load rootfiles
  // ============================
  std::vector<TFile* > file_;
  TString MS= ""; //MadSpin ? "MadSpin" : ""; -> keep old files for the moment
  TString targetfolder=groupSpace+AnalysisFolder;
  file_.push_back(TFile::Open(targetfolder+"/elecDiffXSecSig"+MS+"Summer12PF.root"                                        , "Open"));
  file_.push_back(TFile::Open(targetfolder+"/muonDiffXSecSig"+MS+"Summer12PF.root"                                        , "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/elecDiffXSecSig"+MS+"TopMassConstraint173p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/muonDiffXSecSig"+MS+"TopMassConstraint173p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/elecDiffXSecSig"+MS+"TopMassConstraint171p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/muonDiffXSecSig"+MS+"TopMassConstraint171p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/elecDiffXSecSig"+MS+"TopMassConstraint174p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/muonDiffXSecSig"+MS+"TopMassConstraint174p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/elecDiffXSecSig"+MS+"TopMassConstraint170p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/muonDiffXSecSig"+MS+"TopMassConstraint170p5Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/elecDiffXSecSig"+MS+"TopMassConstraint"+(MS=="" ? "180p0" : "176p5")+"Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/muonDiffXSecSig"+MS+"TopMassConstraint"+(MS=="" ? "180p0" : "176p5")+"Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/elecDiffXSecSig"+MS+"TopMassConstraint"+(MS=="" ? "160p0" : "168p5")+"Summer12PF.root", "Open"));
  file_.push_back(TFile::Open(targetfolder+"/TopMassConstraint/muonDiffXSecSig"+MS+"TopMassConstraint"+(MS=="" ? "160p0" : "168p5")+"Summer12PF.root", "Open"));

  std::vector<double > massConstraint_;
  massConstraint_.push_back(172.5);
  massConstraint_.push_back(172.5);
  massConstraint_.push_back(173.5);
  massConstraint_.push_back(173.5);
  massConstraint_.push_back(171.5);
  massConstraint_.push_back(171.5);
  massConstraint_.push_back(174.5);
  massConstraint_.push_back(174.5);
  massConstraint_.push_back(170.5);
  massConstraint_.push_back(170.5);
  massConstraint_.push_back(MS=="" ? 180.0 : 176.5);
  massConstraint_.push_back(MS=="" ? 180.0 : 176.5);
  massConstraint_.push_back(MS=="" ? 160.0 : 168.5);
  massConstraint_.push_back(MS=="" ? 160.0 : 168.5);

  // list plots of relevance
  std::vector<TString> plotList_, axisLabel_;
  TString plots1D[ ] = {
    // KinFit plots before prob cut 
    "analyzeTopRecoKinematicsKinFit/prob",
    "analyzeTopRecoKinematicsKinFit/chi2",
    "analyzeTopRecoKinematicsKinFit/ttbarAngle",
    "analyzeTopRecoKinematicsKinFit/topPt",
    "analyzeTopRecoKinematicsKinFit/topPtTtbarSys",
    "analyzeTopRecoKinematicsKinFit/topY",
    "analyzeTopRecoKinematicsKinFit/topMass",
    "analyzeTopRecoKinematicsKinFit/ttbarPt",
    "analyzeTopRecoKinematicsKinFit/ttbarY",
    "analyzeTopRecoKinematicsKinFit/ttbarMass",
    "analyzeTopRecoKinematicsKinFit/ttbarHT",
    "analyzeTopRecoKinematicsKinFit/lepPt",
    "analyzeTopRecoKinematicsKinFit/lepEta",
    "analyzeTopRecoKinematicsKinFit/lightqPt",
    "analyzeTopRecoKinematicsKinFit/lightqEta",   
    "analyzeTopRecoKinematicsKinFit/bqPt",
    "analyzeTopRecoKinematicsKinFit/bqEta",
    "analyzeTopRecoKinematicsKinFit/lepEtaPlus",
    "analyzeTopRecoKinematicsKinFit/lepEtaMinus",
    "analyzeTopRecoKinematicsKinFit/topEtaPlus",
    "analyzeTopRecoKinematicsKinFit/topEtaMinus",
    "analyzeTopRecoKinematicsKinFit/lepYPlus",
    "analyzeTopRecoKinematicsKinFit/lepYMinus",
    "analyzeTopRecoKinematicsKinFit/topYPlus",
    "analyzeTopRecoKinematicsKinFit/topYMinus",
    "analyzeTopRecoKinematicsKinFit/topYHad",
    "analyzeTopRecoKinematicsKinFit/topYLep",
    "analyzeTopRecoKinematicsKinFit/topYLead",
    "analyzeTopRecoKinematicsKinFit/topYSubLead",
    "analyzeTopRecoKinematicsKinFit/leadqPt",
    "analyzeTopRecoKinematicsKinFit/bbbarPt",
    "analyzeTopRecoKinematicsKinFit/bbbarY",
    "analyzeTopRecoKinematicsKinFit/bbbarMass", 
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
    "compositedKinematicsKinFit/rhos",
    "compositedKinematicsKinFit/Njets",
    // KinFit plots after prob cut 
    "analyzeTopRecoKinematicsKinFitProbSel/ttbarAngle",
    "analyzeTopRecoKinematicsKinFitProbSel/topPt",
    "analyzeTopRecoKinematicsKinFitProbSel/topPtTtbarSys",
    "analyzeTopRecoKinematicsKinFitProbSel/topY",
    "analyzeTopRecoKinematicsKinFitProbSel/topMass",
    "analyzeTopRecoKinematicsKinFitProbSel/ttbarPt",
    "analyzeTopRecoKinematicsKinFitProbSel/ttbarY",
    "analyzeTopRecoKinematicsKinFitProbSel/ttbarMass",
    "analyzeTopRecoKinematicsKinFitProbSel/ttbarHT",
    "analyzeTopRecoKinematicsKinFitProbSel/lepPt",
    "analyzeTopRecoKinematicsKinFitProbSel/lepEta",
    "analyzeTopRecoKinematicsKinFitProbSel/lightqPt",
    "analyzeTopRecoKinematicsKinFitProbSel/lightqEta",   
    "analyzeTopRecoKinematicsKinFitProbSel/bqPt",
    "analyzeTopRecoKinematicsKinFitProbSel/bqEta",
    "analyzeTopRecoKinematicsKinFitProbSel/lepEtaPlus",
    "analyzeTopRecoKinematicsKinFitProbSel/lepEtaMinus",
    "analyzeTopRecoKinematicsKinFitProbSel/topEtaPlus",
    "analyzeTopRecoKinematicsKinFitProbSel/topEtaMinus",
    "analyzeTopRecoKinematicsKinFitProbSel/lepYPlus",
    "analyzeTopRecoKinematicsKinFitProbSel/lepYMinus",
    "analyzeTopRecoKinematicsKinFitProbSel/topYPlus",
    "analyzeTopRecoKinematicsKinFitProbSel/topYMinus",
    "analyzeTopRecoKinematicsKinFitProbSel/topYHad",
    "analyzeTopRecoKinematicsKinFitProbSel/topYLep",
    "analyzeTopRecoKinematicsKinFitProbSel/topYLead",
    "analyzeTopRecoKinematicsKinFitProbSel/topYSubLead",
    "analyzeTopRecoKinematicsKinFitProbSel/leadqPt",
    "analyzeTopRecoKinematicsKinFitProbSel/bbbarPt",
    "analyzeTopRecoKinematicsKinFitProbSel/bbbarY",
    "analyzeTopRecoKinematicsKinFitProbSel/bbbarMass", 
    "compositedKinematicsProbSel/rhos",
    "compositedKinematicsProbSel/Njets",
  };
  TString axisLabel1D[ ] = { 
    // KinFit plots before prob cut 
    "probability (best fit hypothesis);Events;1;25", 
    "#chi^{2} (best fit hypothesis);Events;1;20",
    "#angle(t,#bar{t});Events;0;35",
    "p_{T}^{t} #left[GeV#right];Top quarks;0;2",
    "p_{T}^{t} #left[GeV#right] (t#bar{t} system);Events;0;20",
    "y^{t};Top quarks;0;5",
    "m^{t};Top quarks;0;20",
    "p_{T}^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;20",
    "y^{t#bar{t}};Top-quark pairs;0;5",
    "m^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;50",
    "H_{T}^{t#bar{t}}=#Sigma(E_{T}(jets)) #left[GeV#right];#frac{dN}{dH_{T}^{t#bar{t}}};0;20",
    "p_{T}^{l} #left[GeV#right];N^{l};0;10",
    "#eta^{l};Leptons;0;5",
    "p_{T}^{q} #left[GeV#right];light jets;0;50",    
    "#eta^{q};light jets;0;10",
    "p_{T}^{b} #left[GeV#right];b jets;0;20",    
    "#eta^{b};b jets;0;5",
    "#eta^{l+};Events;0;5",
    "#eta^{l-};Events;0;5",
    "#eta^{t+};Events;0;25",
    "#eta^{t-};Events;0;25",
    "y^{l+};Events;0;5",
    "y^{l-};Events;0;5",
    "y^{t+};Events;0;5",
    "y^{t-};Events;0;5",
    "y^{t+};Events;0;5",
    "y^{t-};Events;0;5",
    "y^{Lead.t};Events;0;5",
    "y^{NLead.t};Events;0;5",
    "p_{T}^{leading jet} #left[GeV#right];Events;0;25",
    "p_{T}^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;20",  
    "y^{b#bar{b}}(assigned to t#bar{t} system);Events;0;5",
    "m^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;25",
    "#Delta p_{T}^{light jets} #left[GeV#right];Events;0;100",
    "#Delta #eta^{light jets};Events;0;10",
    "#Delta #phi^{light jets};Events;0;20",
    "#Delta p_{T}^{b jets} #left[GeV#right];Events;0;25",
    "#Delta #eta^{b jets};Events;0;20",
    "#Delta #phi^{b jets};Events;0;20",
    "#Delta p_{T}^{lepton} #left[GeV#right];Events;0;20",
    "#Delta #eta^{lepton};Events;0;20",
    "#Delta #phi^{lepton};Events;0;1",
    "#Delta p_{T}^{neutrino} #left[GeV#right];Events;0;20",
    "#Delta #eta^{neutrino};Events;0;50",
    "#Delta #phi^{neutrino};Events;0;50",
    xSecLabelName("rhos" )+";Events;0;11",
    xSecLabelName("Njets")+";Events;0;1",
    // KinFit plots after prob cut 
    "#angle(t,#bar{t});Events;0;35",
    "p_{T}^{t} #left[GeV#right];Top quarks;0;20",
    "p_{T}^{t} #left[GeV#right] (t#bar{t} system);Events;0;20",
    "y^{t};Top quarks;0;5",
    "m^{t};Top quarks;0;20",
    "p_{T}^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;20",
    "y^{t#bar{t}};Top-quark pairs;0;5",
    "m^{t#bar{t}} #left[GeV#right];Top-quark pairs;0;50",
    "H_{T}^{t#bar{t}}=#Sigma(E_{T}(jets)) #left[GeV#right];#frac{dN}{dH_{T}^{t#bar{t}}};0;20",
    "p_{T}^{l} #left[GeV#right];Leptons;0;10",
    "#eta^{l};Leptons;0;5",
    "p_{T}^{q} #left[GeV#right];light jets;0;50",    
    "#eta^{q};light jets;0;10",
    "p_{T}^{b} #left[GeV#right];b jets;0;20",    
    "#eta^{b};b jets;0;5",
    "#eta^{l+};Events;0;5",
    "#eta^{l-};Events;0;5",
    "#eta^{t+};Events;0;25",
    "#eta^{t-};Events;0;25",
    "y^{l+};Events;0;5",
    "y^{l-};Events;0;5",
    "y^{t+};Events;0;5",
    "y^{t-};Events;0;5",
    "y^{t+};Events;0;5",
    "y^{t-};Events;0;5",
    "y^{Lead.t};Events;0;5",
    "y^{NLead.t};Events;0;5",
    "p_{T}^{leading jet} #left[GeV#right];Events;0;25",
    "p_{T}^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;20",  
    "y^{b#bar{b}}(assigned to t#bar{t} system);Events;0;5",
    "m^{b#bar{b}}(assigned to t#bar{t} system) #left[GeV#right];Events;0;25",
    xSecLabelName("rhos" )+";events;0;22",
    xSecLabelName("Njets")+";events;0;1",
  };
  plotList_ .insert(plotList_ .begin(), plots1D    , plots1D    + sizeof(plots1D    )/sizeof(TString));
  axisLabel_.insert(axisLabel_.begin(), axisLabel1D, axisLabel1D+ sizeof(axisLabel1D)/sizeof(TString));
  if(plotList_.size() != axisLabel_.size()){
    std::cout << "ERROR - 1D plots: Number of plots and axis label do not correspond .... Exiting macro!" << std::endl;
    exit(1);
  }
  // run automatically in batch mode if there are many canvas
  if(plotList_.size()>15) gROOT->SetBatch();

  // canvas container
  std::vector<TCanvas*> plotCanvas_;
  // legend
  TLegend* leg= new TLegend(0.6, 0.55, 0.9, 0.88);
  legendStyle(*leg,"t#bar{t} SG MC (m_{#lower[-0.3]{top}}^{#lower[-0.8]{gen}}=172.5 GeV)");
  
  // ============================
  //  get histos
  // ============================
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map<TString, std::vector<double> > binning_ = makeVariableBinning(false);
  // loop all plots
  for(int plot=0; plot<(int)plotList_.size(); ++plot){
    TString name=plotList_[plot];
    // loop all samples
    for(int sample=0; sample<int(file_.size()); sample=sample+2){
      double mcon=massConstraint_[sample];
      if(mcon!=massConstraint_[sample+1]){
	std::cout << "ERROR: listed mass constraint value for sample " << sample << " and " << sample+1 << " are unequal, the expected format is that muon and elec files are listed after each other in file_ with the same constraint value listed in massConstraint_" << std::endl;
      }
      int kcon=roundToInt(mcon*10);
      // nominal mtop
      histo_[name][kcon]     =(TH1F*)file_[sample  ]->Get(name);
      histo_[name][kcon]->Add((TH1F*)file_[sample+1]->Get(name));
      // rebinning
      if(usePAPERbinning&&binning_.count(getStringEntry(name,2,"/"))>0){
	reBinTH1F(*histo_[name][kcon], binning_[getStringEntry(name, 2, "/")], verbose-1);
      }
      else{
	double reBinFactor = atof(((string)getStringEntry(axisLabel_[plot],4,";")).c_str());
	if(reBinFactor>1) equalReBinTH1(reBinFactor, histo_, name, kcon);
      }
      // legend
      if(plot==0){
	if(sample==0) leg->AddEntry(histo_[name][kcon], Form("m_{#lower[-0.1]{top}}^{constr.}=%2.1f GeV", mcon), "L");
	else{
	  TString legentry="#Delta m_{#lower[-0.1]{top}}^{#lower[-0.7]{constr.}}=#scale[0.5]{ }";
	  mcon > 172.5 ? legentry+="+" : legentry+="#scale[0.95]{#bf{#font[122]{-}}}";
	  legentry+= Form("%2.1f GeV", std::abs(mcon-172.5));
	  leg->AddEntry(histo_[name][kcon], legentry, "L");
	}
      }
      // histostyle
      histo_[name][kcon ]->SetStats(false);
      unsigned int color=kBlack;
      if     (mcon==172.5) color=kBlack  ;
      else if(mcon==170.5) color=kAzure+8;
      else if(mcon==171.5) color=kMagenta;
      else if(mcon==173.5) color=kRed;
      else if(mcon==174.5) color=kBlue;
      else if(mcon==160.0) color=kOrange+3;
      else if(mcon==180.0) color=kGreen+3;
      histo_[name][kcon]->SetLineColor(color);
      histo_[name][kcon]->SetLineWidth(2);
      if     (mcon==172.5) histo_[name][kcon]->SetLineStyle(2);
    }
    // ============================
    //  create canvas
    // ============================
    char canvname[10];
    sprintf(canvname,"canv%i",plot);    
    plotCanvas_.push_back( new TCanvas( canvname, canvname, 600, 600) );
    TString canvtitle=getStringEntry(plotList_[plot], 2)+getStringEntry(plotList_[plot], 1);
    plotCanvas_[plotCanvas_.size()-1]->SetTitle(canvtitle);
    // min/max (y-axis)
    double max=histo_[name][1725]->GetMaximum();
    // loop all samples
    for(int sample=0; sample<int(file_.size()); sample=sample+2){
      double mcon=massConstraint_[sample];
      int kcon=roundToInt(mcon*10);
      double newmax=histo_[name][kcon]->GetMaximum();
      if(max<newmax) max=newmax;
    }
    max=1.6*max;
    double min = 0.;
    // log plots
    if(getStringEntry(axisLabel_[plot],3,";")=="1"){
      plotCanvas_[plot]->SetLogy(1);
      min=1;
      max=exp(1.4*(std::log(max)-std::log(min))+std::log(min));
    }
    if(min==0.&&max>0.001) min+=0.001;
    if(min==1.&&max>1.001) min+=0.001;
    axesStyle(*histo_[name][1725], getStringEntry(axisLabel_[plot],1,";"), getStringEntry(axisLabel_[plot],2,";"), min, max);	  
    histo_[name][1725]->GetXaxis()->SetNoExponent(true);
    if(TString(histo_[name][1725]->GetXaxis()->GetTitle()).Contains("#frac")) histo_[name][1725]->GetXaxis()->SetTitleOffset(1.2*histo_[name][1725]->GetXaxis()->GetTitleOffset());
    // adjust range (x-axis)
    double xUp=histo_[name][1725]->GetXaxis()->GetXmax();
    double xDn=histo_[name][1725]->GetXaxis()->GetXmin();
    if(name.Contains("topPt" )){xDn=0.  ; xUp=500.;}
    if(name.Contains("topMass")){xDn=100.  ; xUp=300.;}
    if(name.Contains("topY"  )){xDn=-2.5; xUp=2.5 ;}
    if(name.Contains("ttbarY")){xDn=-2.5; xUp=2.5 ;}
    if(name.Contains("ttbarMass")){xDn=300.;xUp=1100.;}
    if(name.Contains("ttbarHT"  )){xDn=50. ;xUp=500.;}
    if(name.Contains("lepPt" )){xDn=0.;xUp=200.;}
    if(name.Contains("lepEta")){xDn=-2.1;xUp=2.1;}
    if(name.Contains("leadqPt" )){xDn=0.;xUp=300.;}
    if(name.Contains("leadqEta")){xDn=-3.;xUp=3.;}
    if(name.Contains("bqPt" )){xDn=0. ;xUp=400.;}
    if(name.Contains("lightqPt" )){xDn=0. ;xUp=400.;}


    if(name.Contains("bqEta")){xDn=-2.4;xUp=2.4;}
    if(name.Contains("bbbarPt")){xDn=0.;xUp=500.;}
    if(name.Contains("bbbarY"   )){xDn=-2.5;xUp=2.5;}
    if(name.Contains("bbbarMass")){xDn=0.;xUp=504.;}
    if(name.Contains("shiftLqPt"  )){xDn=-80    ;xUp=80;    }
    if(name.Contains("shiftLqEta" )){xDn=-0.02  ;xUp=0.02;  }
    if(name.Contains("shiftLqPhi" )){xDn=-0.02  ;xUp=0.02;  }
    if(name.Contains("shiftBqPt"  )){xDn=-30    ;xUp=30;    }
    if(name.Contains("shiftBqEta" )){xDn=-0.012 ;xUp=0.012; }
    if(name.Contains("shiftBqPhi" )){xDn=-0.008 ;xUp=0.008; }
    if(name.Contains("shiftLepPt" )){xDn=-0.4   ;xUp=0.4;   }
    if(name.Contains("shiftLepEta")){xDn=-0.000001 ;xUp=0.000001; }
    if(name.Contains("shiftLepPhi")){xDn=-0.0000005;xUp=0.0000005;}
    if(name.Contains("shiftNuPt"  )){xDn=-60    ;xUp=60;    }
    if(name.Contains("shiftNuEta" )){xDn=-4.    ;xUp=4.;    }
    if(name.Contains("shiftNuPhi" )){xDn=-0.3  ;xUp=0.3;  }
//     if(name.Contains("shiftLqEta")){xDn=-0.4;xUp=0.4;}
//     if(name.Contains("shiftLqPhi")){xDn=-0.2;xUp=0.2;}
//     if(name.Contains("shiftBqPt")){xDn=-30.;xUp=30.;}
//     if(name.Contains("shiftBqEta")){xDn=-0.02;xUp=0.02;}
//     if(name.Contains("shiftBqPhi")){xDn=-0.1;xUp=0.1;}
//     if(name.Contains("shiftLepPt")){xDn=-5.;xUp=5.;}
//     if(name.Contains("shiftNuPt" )){xDn=-50.;xUp=50.;}
//     if(name.Contains("shiftNuPhi" )){xDn=-1.;xUp=1.;}
    if(name.Contains("Njets")){xDn=3.51;xUp=9.49;}
    //if(name.Contains("")){xDn=;xUp=;}
    if(xUp!=histo_[name][1725]->GetXaxis()->GetXmax()||xDn!=histo_[name][1725]->GetXaxis()->GetXmin()) histo_[name][1725]->GetXaxis()->SetRangeUser(xDn,xUp-0.000001);
    
    // ============================
    //  plotting
    // ============================
    histo_[name][1725]->Draw("AXIS");
    // loop all samples
    for(int sample=0; sample<int(file_.size()); sample=sample+2){
      double mcon=massConstraint_[sample];
      int kcon=roundToInt(mcon*10);
      histo_[name][kcon]->Draw("hist same");
    }
    histo_[name][1725]->Draw("hist same");
    if(!PHD||(name.Contains("bqPt")||(name.Contains("lepPt")||(name.Contains("topPt")&&!name.Contains("Sys"))||name.Contains("ttbarMass")))) leg->Draw("same");
    // draw cut label
    TString cutLabel="1 lepton, #geq4 Jets";
    if(name.Contains("Tagged")||plotList_[plot].Contains("AfterBtagging")) cutLabel+=", #geq2 b-tags";
    else if(name.Contains("KinFit")) cutLabel+=", #geq2 b-tags, KinFit";
    if(plotList_[plot].Contains("ProbSel")){
      if(name.Contains("KinFit")) cutLabel.ReplaceAll("KinFit","prob>0.02");
      else cutLabel+=", #geq2 b-tags, prob>0.02";
    }
    double positionX=xUp+0.045*(xUp-xDn)*(gStyle->GetCanvasDefW()/600.);
    double positionY=min;
    //double xUpp=histo_[name][1725]->GetBinLowEdge(histo_[name][1725]->FindBin(xUp)+1);
    //if((name.Contains("bbbarPt"))||(name.Contains("bbbarMass"))){
    //  positionX=xUpp+0.045*(xUpp-xDn)*(gStyle->GetCanvasDefW()/600.);
    //  std::cout << name << ": " << xDn << ".." << xUp << " -> " << positionX << std::endl;
    //}
    TLatex *sellabel = new TLatex(positionX,positionY,cutLabel);
    sellabel->SetTextAlign(11);
    sellabel->SetTextAngle(90);
    sellabel->SetTextSize(0.04);
    sellabel->Draw("same");
    DrawDecayChLabel("e/#mu + Jets Combined");
    DrawCMSLabels(prelim, 19714., 0.04, true, false, false);
    // ratio
    std::vector<double> err_;
    for(int bin=1; bin<=histo_[name][1725]->GetNbinsX(); ++bin){
      double a =histo_[name][1725]->GetBinContent(bin);
      err_.push_back(sqrt(a)/(a));
    }
    TString nominator="m_{top} constrain";
    TString denominator="172.5 GeV";
    drawRatio(histo_[name][1725], histo_[name][1725], 0.9, 1.12, myStyle, verbose, err_, denominator, nominator, "e p", kBlack, false, 0.2);
    for(int sample=0; sample<int(file_.size()); sample=sample+2){
      double mcon=massConstraint_[sample];
      int kcon=roundToInt(mcon*10);
      //std::cout << kcon << std::endl;      
      //gStyle->SetErrorX(0.5);
      drawRatio(histo_[name][kcon], histo_[name][1725], 0.9, 1.12, myStyle, verbose, err_, denominator, nominator, "hist same", histo_[name][kcon]->GetLineColor(), false, 0.2);
    }
    drawRatio(histo_[name][1725], histo_[name][1725], 0.9, 1.12, myStyle, verbose, err_, denominator, nominator, "e p same", kBlack, false, 0.2);
  }

  if(save){
    if(verbose==0) gErrorIgnoreLevel=kWarning;
    saveCanvas(plotCanvas_, "diffXSecFromSignal/plots/combined/2012/massConstraintTest/topMassBias", "Compilation", true, true);
  }
}
