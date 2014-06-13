#include "basicFunctions.h"

TH1F* convertToHist(TGraphAsymmErrors* in,  TH1F* templatePlot, double NbinsIn);

void ATLASCMSCOMPARISON(TString quantity="ttbarY"){

  bool cmssim=true;

  // ---
  //    canvas style 
  // ---
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.SetErrorX(0.5);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gStyle->SetEndErrorSize(10);
  gStyle->SetOptFit(0);
  TGaxis::SetMaxDigits(2);

  // ---
  //    top Pt 7 TeV
  // ---
  int Nbins7=7;
  if(     quantity=="topPt"    ) Nbins7=7;
  else if(quantity=="ttbarMass") Nbins7=5;
  else if(quantity=="ttbarY"   ) Nbins7=6;

  // CMS data
  TGraphAsymmErrors* CMSdata7 = new TGraphAsymmErrors(Nbins7);
  if(quantity=="topPt"){
    CMSdata7->SetPoint( 0, 25.0 , 0.004032 ); 
    CMSdata7->SetPoint( 1, 75.0 , 0.006746 ); 
    CMSdata7->SetPoint( 2, 125.0, 0.004737 ); 
    CMSdata7->SetPoint( 3, 175.0, 0.002506 ); 
    CMSdata7->SetPoint( 4, 225.0, 0.001140 ); 
    CMSdata7->SetPoint( 5, 300.0, 0.000334 ); 
    CMSdata7->SetPoint( 6, 575.0, 0.000019 ); 
    
    CMSdata7->SetPointError( 0, 25., 25., ( 5.5/100)*0.004032, ( 5.5/100)*0.004032 );
    CMSdata7->SetPointError( 1, 25., 25., ( 4.1/100)*0.006746, ( 4.1/100)*0.006746 );
    CMSdata7->SetPointError( 2, 25., 25., ( 4.0/100)*0.004737, ( 4.0/100)*0.004737 );
    CMSdata7->SetPointError( 3, 25., 25., ( 5.4/100)*0.002506, ( 5.4/100)*0.002506 );
    CMSdata7->SetPointError( 4, 25., 25., ( 5.6/100)*0.001140, ( 5.6/100)*0.001140 );
    CMSdata7->SetPointError( 5, 50., 50., ( 8.4/100)*0.000334, ( 8.4/100)*0.000334 );
    CMSdata7->SetPointError( 6, 225.,225., (14.1/100)*0.000019, (14.1/100)*0.000019 );
  }
  else if(quantity=="ttbarMass"){
    CMSdata7->SetPoint( 0, 350.0 , 0.002588 ); 
    CMSdata7->SetPoint( 1, 500.0 , 0.002685 ); 
    CMSdata7->SetPoint( 2, 625.0 , 0.000953 ); 
    CMSdata7->SetPoint( 3, 825.0 , 0.000232 ); 
    CMSdata7->SetPoint( 4, 1725.0, 0.000008 ); 
    
    CMSdata7->SetPointError( 0, 100., 100., ( 4.9 /100)*0.002588, (4.9 /100)*0.002588 );
    CMSdata7->SetPointError( 1, 50. , 50. , ( 6.9 /100)*0.002685, (6.9 /100)*0.002685 );
    CMSdata7->SetPointError( 2, 75. , 75. , ( 7.7 /100)*0.000953, (7.7 /100)*0.000953 );
    CMSdata7->SetPointError( 3, 125., 125., ( 14.4/100)*0.000232, (14.4/100)*0.000232 );
    CMSdata7->SetPointError( 4, 775., 775., ( 27.6/100)*0.000008, (27.6/100)*0.000008 );
  }
  else if(quantity=="ttbarY"){
    CMSdata7->SetPoint( 0, -1.75, 0.082140 ); 
    CMSdata7->SetPoint( 1, -0.75, 0.318979 ); 
    CMSdata7->SetPoint( 2, -0.25, 0.427823 ); 
    CMSdata7->SetPoint( 3,  0.25, 0.44591  ); 
    CMSdata7->SetPoint( 4,  0.75, 0.318820 ); 
    CMSdata7->SetPoint( 5,  1.75, 0.080457 );

    CMSdata7->SetPointError( 0, 0.75, 0.75, (8.1/100)*0.082140, (8.1/100)*0.082140 );
    CMSdata7->SetPointError( 1, 0.25, 0.25, (3.3/100)*0.318979, (3.3/100)*0.318979 );
    CMSdata7->SetPointError( 2, 0.25, 0.25, (3.0/100)*0.427823, (3.0/100)*0.427823 );
    CMSdata7->SetPointError( 3, 0.25, 0.25, (3.5/100)*0.44591 , (3.5/100)*0.44591  );
    CMSdata7->SetPointError( 4, 0.25, 0.25, (3.5/100)*0.318820, (3.5/100)*0.318820 );
    CMSdata7->SetPointError( 5, 0.75, 0.75, (5.9/100)*0.080457, (5.9/100)*0.080457 );
  }

  CMSdata7->SetLineWidth(3.);
  CMSdata7->SetMarkerSize(1.2);
  CMSdata7->SetMarkerStyle(24);
  CMSdata7->SetLineStyle(1);
  CMSdata7->SetMarkerColor(kBlue);
  CMSdata7->SetLineColor(kBlue);
  CMSdata7->SetFillStyle(3004);
  CMSdata7->SetFillColor(kBlue);

  // CMS MadGraph+Pythia(Z2*)
  TGraphAsymmErrors* CMSMadGraph7 = new TGraphAsymmErrors(Nbins7);
  if(quantity=="topPt"){
    CMSMadGraph7->SetPoint( 0, 25.0 , 0.003331 ); 
    CMSMadGraph7->SetPoint( 1, 75.0 , 0.006495 ); 
    CMSMadGraph7->SetPoint( 2, 125.0, 0.005077 ); 
    CMSMadGraph7->SetPoint( 3, 175.0, 0.002748 ); 
    CMSMadGraph7->SetPoint( 4, 225.0, 0.001282 ); 
    CMSMadGraph7->SetPoint( 5, 300.0, 0.000413 ); 
    CMSMadGraph7->SetPoint( 6, 575.0, 0.000027 ); 

    CMSMadGraph7->SetPointError( 0, 25., 25., 0., 0. );
    CMSMadGraph7->SetPointError( 1, 25., 25., 0., 0. );
    CMSMadGraph7->SetPointError( 2, 25., 25., 0., 0. );
    CMSMadGraph7->SetPointError( 3, 25., 25., 0., 0. );
    CMSMadGraph7->SetPointError( 4, 25., 25., 0., 0. );
    CMSMadGraph7->SetPointError( 5, 50., 50., 0., 0. );
    CMSMadGraph7->SetPointError( 6, 225., 225., 0., 0. );
  }
  else if(quantity=="ttbarMass"){
    CMSMadGraph7->SetPoint( 0, 350.0 , 0.002541 ); 
    CMSMadGraph7->SetPoint( 1, 500.0 , 0.002759 ); 
    CMSMadGraph7->SetPoint( 2, 625.0 , 0.000988 ); 
    CMSMadGraph7->SetPoint( 3, 825.0 , 0.000222 ); 
    CMSMadGraph7->SetPoint( 4, 1725.0, 0.000008 ); 
    
    CMSMadGraph7->SetPointError( 0, 100., 100., 0., 0. );
    CMSMadGraph7->SetPointError( 1, 50. , 50. , 0., 0. );
    CMSMadGraph7->SetPointError( 2, 75. , 75. , 0., 0. );
    CMSMadGraph7->SetPointError( 3, 125., 125., 0., 0. );
    CMSMadGraph7->SetPointError( 4, 775., 775., 0., 0. );
  }
  else if(quantity=="ttbarY"){
    CMSMadGraph7->SetPoint( 0, -1.75, 0.088374 ); 
    CMSMadGraph7->SetPoint( 1, -0.75, 0.319964 ); 
    CMSMadGraph7->SetPoint( 2, -0.25, 0.414022 ); 
    CMSMadGraph7->SetPoint( 3,  0.25, 0.415108 ); 
    CMSMadGraph7->SetPoint( 4,  0.75, 0.320370 ); 
    CMSMadGraph7->SetPoint( 5,  1.75, 0.088218 );

    CMSMadGraph7->SetPointError( 0, 0.75, 0.75, 0., 0. );
    CMSMadGraph7->SetPointError( 1, 0.25, 0.25, 0., 0. );
    CMSMadGraph7->SetPointError( 2, 0.25, 0.25, 0., 0. );
    CMSMadGraph7->SetPointError( 3, 0.25, 0.25, 0., 0. );
    CMSMadGraph7->SetPointError( 4, 0.25, 0.25, 0., 0. );
    CMSMadGraph7->SetPointError( 5, 0.75, 0.75, 0., 0. ); 
  }

  CMSMadGraph7->SetLineWidth(3.);
  CMSMadGraph7->SetMarkerSize(1.2);
  CMSMadGraph7->SetLineStyle(1);
  CMSMadGraph7->SetMarkerStyle(20);
  CMSMadGraph7->SetMarkerColor(kAzure+6);
  CMSMadGraph7->SetLineColor(kAzure+6);


  // ATLAS data
  TGraphAsymmErrors* ATLASdata7 = new TGraphAsymmErrors(Nbins7);
  if(quantity=="topPt"){
    ATLASdata7->SetPoint( 0, 25.0 , 0.0034  ); 
    ATLASdata7->SetPoint( 1, 75.0 , 0.0067  ); 
    ATLASdata7->SetPoint( 2, 125.0, 0.0052  ); 
    ATLASdata7->SetPoint( 3, 175.0, 0.00266 ); 
    ATLASdata7->SetPoint( 4, 225.0, 0.00114 ); 
    ATLASdata7->SetPoint( 5, 300.0, 0.00033 ); 
    ATLASdata7->SetPoint( 6, 575.0, 0.000018); 

    ATLASdata7->SetPointError( 0, 25., 25., (4.47 /100)*0.0034  , (4.47 /100)*0.0034   );
    ATLASdata7->SetPointError( 1, 25., 25., (1.41 /100)*0.0067  , (1.41 /100)*0.0067   );
    ATLASdata7->SetPointError( 2, 25., 25., (2.83 /100)*0.0052  , (2.83 /100)*0.0052   );
    ATLASdata7->SetPointError( 3, 25., 25., (3.61 /100)*0.00266 , (3.61 /100)*0.00266  );
    ATLASdata7->SetPointError( 4, 25., 25., (3.61 /100)*0.00114 , (3.61 /100)*0.00114  );
    ATLASdata7->SetPointError( 5, 50., 50., (5.83 /100)*0.00033 , (5.83 /100)*0.00033  );
    ATLASdata7->SetPointError( 6, 225., 225., (11.66/100)*0.000018, (11.66/100)*0.000018 );
  }
  else if(quantity=="ttbarMass"){
    ATLASdata7->SetPoint( 0, 350.0 , 0.00250   ); 
    ATLASdata7->SetPoint( 1, 500.0 , 0.00273   ); 
    ATLASdata7->SetPoint( 2, 625.0 , 0.00102   ); 
    ATLASdata7->SetPoint( 3, 825.0 , 0.00023   ); 
    ATLASdata7->SetPoint( 4, 1725.0, 0.0000076 ); 
    
    ATLASdata7->SetPointError( 0, 100., 100., 0.00008   , 0.00008  );
    ATLASdata7->SetPointError( 1, 50. , 50. , 0.00007   , 0.00007  );
    ATLASdata7->SetPointError( 2, 75. , 75. , 0.00004   , 0.00004  );
    ATLASdata7->SetPointError( 3, 125., 125., 0.00001   , 0.00001  );
    ATLASdata7->SetPointError( 4, 775., 775., 0.0000005, 0.0000005 );
  }
  else if(quantity=="ttbarY"){
    ATLASdata7->SetPoint( 0, -1.75, 0.081 ); 
    ATLASdata7->SetPoint( 1, -0.75, 0.321 ); 
    ATLASdata7->SetPoint( 2, -0.25, 0.436 ); 
    ATLASdata7->SetPoint( 3,  0.25, 0.423 ); 
    ATLASdata7->SetPoint( 4,  0.75, 0.321 ); 
    ATLASdata7->SetPoint( 5,  1.75, 0.087 );

    ATLASdata7->SetPointError( 0, 0.75, 0.75, 0.003, 0.003);
    ATLASdata7->SetPointError( 1, 0.25, 0.25, 0.009, 0.009);
    ATLASdata7->SetPointError( 2, 0.25, 0.25, 0.009, 0.009);
    ATLASdata7->SetPointError( 3, 0.25, 0.25, 0.007, 0.007);
    ATLASdata7->SetPointError( 4, 0.25, 0.25, 0.005, 0.005);
    ATLASdata7->SetPointError( 5, 0.75, 0.75, 0.005, 0.005);
  }
  ATLASdata7->SetLineWidth(3.);
  ATLASdata7->SetMarkerSize(1.2);
  ATLASdata7->SetMarkerStyle(22);
  ATLASdata7->SetLineStyle(2);
  ATLASdata7->SetMarkerColor(kRed);
  ATLASdata7->SetLineColor(kRed);
  ATLASdata7->SetFillStyle(3005);
  ATLASdata7->SetFillColor(kRed);

  // ---
  //    dummy plots for axis
  // ---
  // create variable bin edges
  std::map<TString, std::vector<double> > binning_ = makeVariableBinning(false);
  std::vector<double> newTopBins_;
  if(quantity=="topPt"){
    newTopBins_.push_back(0.);
    newTopBins_.push_back(50.);
    newTopBins_.push_back(100.);
    newTopBins_.push_back(150.);
    newTopBins_.push_back(200.);
    newTopBins_.push_back(250.);
    newTopBins_.push_back(350.);
    newTopBins_.push_back(800.);
  }

  else if(quantity=="ttbarMass"){
    newTopBins_.push_back(250.);
    newTopBins_.push_back(450.);
    newTopBins_.push_back(550.);
    newTopBins_.push_back(700.);
    newTopBins_.push_back(950.);
    newTopBins_.push_back(2500.);
  }
  else if(quantity=="ttbarY"   ){
    newTopBins_.push_back(-2.5);
    newTopBins_.push_back(-1.0);
    newTopBins_.push_back(-0.5);
    newTopBins_.push_back(0.  );
    newTopBins_.push_back(0.5 );
    newTopBins_.push_back(1.0 );
    newTopBins_.push_back(2.5 );
  }

  binning_[quantity]=newTopBins_;
  double start=0.;
  double range=800.;
  int Nfinebins=800;
  if(quantity=="topPt"){ Nfinebins=800; start=0.; range=800.;}
  else if(quantity=="ttbarMass"){ Nfinebins=2250; start=250.; range=2500.; }
  else if(quantity=="ttbarY"   ){ Nfinebins=500 ; start=-2.5; range=2.5;}
  TH1F* dummy= new TH1F("","",Nfinebins,start,range);
  reBinTH1F(*dummy, binning_[quantity], 0);
  histogramStyle(*dummy, kSig);
  TString label="p_{T}^{t} [GeV]";
  TString label2="p_{T}^{t}";
  if(quantity=="topPt"){label="p_{T}^{t} [GeV]"; label2="p_{T}^{t}";}
  else if(quantity=="ttbarMass"){ label="m^{t#bar{t}} [GeV]"; label2="m^{t#bar{t}}"; }
  else if(quantity=="ttbarY"   ){ label="y^{t#bar{t}}"; label2=label; }
  dummy->GetXaxis()->SetTitle(label);
  dummy->GetYaxis()->SetTitle(TString("#frac{1}{#sigma} #frac{d#sigma}{")+label2+"}");
  double max=0.01;
  double min=0.;
  if(quantity=="topPt"    ){max=0.01 ; min=0.;}
  if(quantity=="ttbarMass"){max=0.004; min=0.;}
  if(quantity=="ttbarY"   ){max=1.0  ; min=0.;}
  dummy->SetMaximum(max);
  dummy->SetMinimum(min);
  dummy->GetXaxis()->SetNoExponent(true);
  TH1F* A=convertToHist(ATLASdata7, dummy, Nbins7);
  TH1F* C=convertToHist(CMSdata7  , dummy, Nbins7);
  histogramStyle(*A, kSig, true, 1.2, kRed);
  A->SetLineWidth(3.);
  A->SetMarkerSize(1.2);
  A->SetMarkerStyle(22);
  A->SetLineStyle(1);
  A->SetMarkerColor(kRed);
  A->SetLineColor(kRed);
  A->SetFillStyle(3005);
  histogramStyle(*C, kSig, true, 1.2, kBlue);
  C->SetLineWidth(3.);
  C->SetMarkerSize(1.2);
  C->SetMarkerStyle(24);
  C->SetLineStyle(1);
  C->SetMarkerColor(kBlue);
  C->SetLineColor(kBlue);
  C->SetFillStyle(3004);
  TH1F* M=convertToHist(CMSMadGraph7  , dummy, Nbins7);
  histogramStyle(*M, kSig, false, 1.2, kAzure+6);
  M->SetLineWidth(3.);
  M->SetMarkerSize(1.2);
  M->SetMarkerStyle(22);
  M->SetLineStyle(1);
  M->SetMarkerColor(kAzure+6);
  M->SetLineColor(kAzure+6);


  // ---
  //    legend
  // ---
  TLegend *leg0 = new TLegend(0.38, 0.63, 0.95, 0.89);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  leg0->SetHeader("LHC, #sqrt{s}=7TeV (t#bar{t} #rightarrow e/#mu+Jets final state)");
  TLegend *leg1=(TLegend*)leg0->Clone(); 
  if(cmssim) leg0->AddEntry( CMSMadGraph7, "CMS MadGraph+Pythia(Z2*)", "LP");
  leg0->AddEntry( CMSdata7    , "#splitline{CMS data, 5.0 fb^{-1} (Eur.Phys.J.C73(2013)2339}{in ATLAS binning)}", "LP");
  leg0->AddEntry( ATLASdata7  , "ATLAS data, 4.6 fb^{-1}"                         , "LP");
  leg1->AddEntry( ATLASdata7  , "#splitline{ATLAS data, 4.6 fb^{-1}}{(ATLAS-CONF-2013-099)}"   , "FP");
  leg1->AddEntry( CMSdata7    , "#splitline{CMS data, 5.0 fb^{-1} (Eur.Phys.J. C73}{(2013) 2339, redone in ATLAS binning)}", "FP");
  if(cmssim) leg1->AddEntry( CMSMadGraph7, "#splitline{CMS MC (MadGraph+Pythia,}{CTEQ6L1 PDF, Z2* UE tune)}", "LP");

  // ---
  //    privatworklabel
  // ---
  TPaveText *privatworklabel = new TPaveText();
  privatworklabel -> SetX1NDC(gStyle->GetPadLeftMargin());
  privatworklabel -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
  privatworklabel -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
  privatworklabel -> SetY2NDC(1.0);
  privatworklabel -> SetTextFont(42);
  privatworklabel -> AddText("private work");
  privatworklabel->SetFillStyle(0);
  privatworklabel->SetBorderSize(0);
  privatworklabel->SetTextSize(0.04);
  privatworklabel->SetTextAlign(32);

  // canvas
  std::vector<TCanvas*> plotCanvas_;
  // a) linear
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(quantity+" measurement comparison");
  // drawing
  dummy->Draw("axis");
  if(cmssim) CMSMadGraph7->Draw("p e1 same");
  ATLASdata7  ->Draw("p e2 same");
  CMSdata7    ->Draw("p e2 same");
  leg0 ->Draw("same");
  if(!PHD) privatworklabel->Draw("same");
  // b) log scale
  TH1F* dummy2=(TH1F*)dummy->Clone();
  dummy2->SetMinimum(0.00001);
  dummy2->SetMaximum(max*10);
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(quantity+" measurement comparison");
  plotCanvas_[plotCanvas_.size()-1]->SetLogy();
  // drawing
  dummy2->Draw("axis");
  if(cmssim) CMSMadGraph7->Draw("p e1 same");
  CMSdata7    ->Draw("p e2 same");
  ATLASdata7  ->Draw("p e2 same");
  leg0 ->Draw("same");
  if(!PHD) privatworklabel->Draw("same");
  // c) both data with bands
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(quantity+" measurement comparison");
  // drawing
  dummy->Draw("axis");
  //DrawSteps(C, "e2 same");
  //C->Draw("e2 same");
  //DrawSteps(A, "e2 same");
  //A->Draw("e2 same");
  if(cmssim) CMSMadGraph7->Draw("p e1 same");
  CMSdata7    ->Draw("p e2 same");
  ATLASdata7  ->Draw("p e2 same");
  // new pad for log plot
  TPad *rPad = new TPad("rPad","",0.4,0.15,0.95,0.85);
  rPad->SetFillStyle(0);
  rPad->SetFillColor(0);
  rPad->SetBorderSize(0);
  rPad->SetBorderMode(0);
  rPad->SetLogy(0);
  rPad->SetLogx(0);
  rPad->SetTicky(1);
  rPad->Draw("");
  rPad->cd();
  // log plot curves
  rPad->SetLogy();
  dummy2->Draw("axis");
  if(cmssim) CMSMadGraph7->Draw("p e1 same");
  CMSdata7    ->Draw("p e2 same");
  ATLASdata7  ->Draw("p e2 same");
  leg1->Draw("same");
  if(!PHD) privatworklabel->Draw("same");
  // d) both data with ratio
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(quantity+" measurement comparison");
  // drawing
  dummy->Draw("axis");
  if(cmssim) CMSMadGraph7->Draw("p e1 same");
  //DrawSteps(C, "e2 same");
  //C->Draw("e2 same");
  //DrawSteps(A, "e2 same");
  //A->Draw("e2 same");
  CMSdata7    ->Draw("p e2 same");
  ATLASdata7  ->Draw("p e2 same");
  leg1 ->Draw("same");
  if(!PHD) privatworklabel->Draw("same");
  std::vector<double> errA_;
  for(int bin=1; bin<=Nbins7; ++bin){
    errA_.push_back((ATLASdata7->GetY()[bin-1]/CMSdata7->GetY()[bin-1])*((ATLASdata7->GetErrorYhigh(bin-1)/ATLASdata7->GetY()[bin-1])));
  }
  std::vector<double> errC_;
  for(int bin=1; bin<=Nbins7; ++bin){
    errC_.push_back((CMSdata7->GetErrorYhigh(bin-1)/CMSdata7->GetY()[bin-1]));
  }
  TString denomlab  ="data";
  TString nomlab="#scale[0.75]{MadGraph+Pythia}";
  if(cmssim){
    std::vector<double> errM_;
    for(int bin=1; bin<=Nbins7; ++bin){
      errM_.push_back(0.);
      //errM_.push_back((CMSMadGraph7->GetY()[bin-1]/CMSdata7->GetY()[bin-1])*((CMSMadGraph7->GetErrorYhigh(bin-1)/CMSMadGraph7->GetY()[bin-1])));
    }
    drawRatio(dummy, dummy, 0.5, 1.5,  myStyle, 0, errM_, nomlab, denomlab, "axis", 10, false, 0.7);
    //drawRatio(M, M, 0.5, 1.5,  myStyle, 0, errM_, nomlab, denomlab, "hist same" , kAzure+6, false, 0.7);
    drawRatio(C, M, 0.5, 1.5,  myStyle, 0, errC_, nomlab, denomlab, "p e2 same" , kBlue, true, 0.7);
    drawRatio(A, M, 0.5, 1.5,  myStyle, 0, errA_, nomlab, denomlab, "p e2 same" , kRed , true, 0.7);
  }

  //saving
  TString path="./diffXSecFromSignal/plots/combined/2012/comparisonATLAS/";
  plotCanvas_[0]->Print(path+quantity+"ATLASvsCMS7TeV.eps");
  plotCanvas_[0]->Print(path+quantity+"ATLASvsCMS7TeV.png");
  plotCanvas_[1]->Print(path+quantity+"ATLASvsCMS7TeVLog.eps");
  plotCanvas_[1]->Print(path+quantity+"ATLASvsCMS7TeVLog.png");
  plotCanvas_[2]->Print(path+quantity+"ATLASvsCMS7TeVNoratio.eps");
  plotCanvas_[2]->Print(path+quantity+"ATLASvsCMS7TeVNoratio.png");
  plotCanvas_[plotCanvas_.size()-1]->Print(path+quantity+"ATLASvsCMS7TeVratio.eps");
  plotCanvas_[plotCanvas_.size()-1]->Print(path+quantity+"ATLASvsCMS7TeVratio.png");
  plotCanvas_[plotCanvas_.size()-1]->SetTitle(quantity);
  saveToRootFile("ATLASvsCMSDataComparisonPlots.root", plotCanvas_[plotCanvas_.size()-1], true, 0,"");
}

TH1F* convertToHist(TGraphAsymmErrors* in,  TH1F* templatePlot, double NbinsIn){
  // create poutput histo from template plot  
  std::cout << "out: " << templatePlot->GetNbinsX() << ", in: " << NbinsIn << std::endl;
  TH1F* out=(TH1F*)templatePlot->Clone(TString(in->GetTitle())+"Histo");
  if(out->GetNbinsX()==NbinsIn){
    // loop bins
    for(unsigned int bin=1; bin<=(unsigned int)templatePlot->GetNbinsX(); ++bin){
      out->SetBinContent(bin, in->GetY()[bin-1]);
      out->SetBinError(bin, in->GetErrorYhigh(bin-1));
    }
  }
  else{
    std::cout << "ERROR: !!! different number of bins for input TgRaphAsymmErrors and output TH1F* !!!" << std::endl;
  }
  // return the result
  return out;
}
