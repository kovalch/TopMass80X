#include "basicFunctions.h"

void addTAE(TGraphAsymmErrors* in, TGraphAsymmErrors* out);

void getMCcorrectionfactors(TString func = "exp"){

  // some parameters
  bool grey=false;
  bool plotsepfit   =true;
  bool plotfiterrors=false;
  TString optD= plotsepfit ? "" : "0";
  double xmax=500;
  TString Txmax=getTStringFromDouble(xmax);

  // colors
  int color7=kRed-4;
  int color8=kBlue+2;
  int ljets7color=color7;//kRed;//kRed+1;
  int dilep7color=color7;//kOrange+7;
  int ljets8color=color8;//kBlue+2;//kBlue;
  int dilep8color=color8;//kGreen+2;//kAzure+6;
  int fit7color=color7;//kMagenta-4;
  int fit8color=color8;//kTeal+3;
  int colorband=kCyan-7;

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

  // ---
  //    collect all curves
  // ---
  // 7 TeV
  int NbinsLjets7=7;
  TGraphAsymmErrors* SFljets = new TGraphAsymmErrors(NbinsLjets7);
  int NbinsDilep7=5;
  TGraphAsymmErrors* SFdilep = new TGraphAsymmErrors(NbinsDilep7);
  //int Nbins7=NbinsLjets7+NbinsDilep7;
  TGraphAsymmErrors* SF7 = new TGraphAsymmErrors(0);
  // 8 TeV
  int NbinsLjets8=8;
  TGraphAsymmErrors* SFljets8 = new TGraphAsymmErrors(NbinsLjets8);
  int NbinsDilep8=5;
  TGraphAsymmErrors* SFdilep8 = new TGraphAsymmErrors(NbinsDilep8);
  //int Nbins8=NbinsLjets8+NbinsDilep8;
  TGraphAsymmErrors* SF8 = new TGraphAsymmErrors(0);
  // combined
  //int Nbins=NbinsLjets7+NbinsDilep7+NbinsLjets8+NbinsDilep8;
  TGraphAsymmErrors* SF = new TGraphAsymmErrors(0);
  // ---
  //    top Pt data / MC ratio
  // ---
  // a) l+jets 7TeV data points
  //           bin x(BCC)    data  / Madgraph         // BCCNNLO // BCC MG
  SFljets->SetPoint( 0, 28   , 0.004536 / 0.003806 ); //28       // 26.2
  SFljets->SetPoint( 1, 85.6 , 0.006658 / 0.006574 ); //85.6     // 88.8
  SFljets->SetPoint( 2, 125  , 0.004740 / 0.004740 ); //125      // 126.2
  SFljets->SetPoint( 3, 173.6, 0.002501 / 0.002748 ); //173.6    // 173.8
  SFljets->SetPoint( 4, 227.5, 0.001042 / 0.001195 ); //227.5    // 228.8
  SFljets->SetPoint( 5, 287.3, 0.000378 / 0.000454 ); //287.3    // 288.8
  SFljets->SetPoint( 6, 355.8, 0.000120 / 0.000154 ); //355.8    // 356.2
  //                   x errors   rel.err(data) *( data  / Madgraph)
  SFljets->SetPointError( 0, 0., 0., (4.4 /100.)*(0.004536 / 0.003806), (4.4 /100.)*(0.004536 / 0.003806) );
  SFljets->SetPointError( 1, 0., 0., (5.5 /100.)*(0.006658 / 0.006574), (5.5 /100.)*(0.006658 / 0.006574) );
  SFljets->SetPointError( 2, 0., 0., (4.0 /100.)*(0.004740 / 0.004740), (4.0 /100.)*(0.004740 / 0.004740) );
  SFljets->SetPointError( 3, 0., 0., (5.8 /100.)*(0.002501 / 0.002748), (5.8 /100.)*(0.002501 / 0.002748) );
  SFljets->SetPointError( 4, 0., 0., (6.2 /100.)*(0.001042 / 0.001195), (6.2 /100.)*(0.001042 / 0.001195) );
  SFljets->SetPointError( 5, 0., 0., (9.0 /100.)*(0.000378 / 0.000454), (9.0 /100.)*(0.000378 / 0.000454) );
  SFljets->SetPointError( 6, 0., 0., (11.1/100.)*(0.000120 / 0.000154), (11.1/100.)*(0.000120 / 0.000154) );
  //style of ratio
  SFljets->SetLineWidth(3.);
  SFljets->SetMarkerSize(1.5);
  SFljets->SetMarkerStyle(26);
  SFljets->SetMarkerColor(ljets7color);
  SFljets->SetLineColor(ljets7color);

  // b) dilepton 7TeV data points
  //           bin x(BCC)    data  / Madgraph               // BCCNNLO // BCC MG
  SFdilep->SetPoint( 0, 33.7,  (0.00509572 / 0.00453114 )  );// 33.7    // 34 
  SFdilep->SetPoint( 1, 107 ,  (0.00626002 / 0.00600115 )  );// 106     // 107
  SFdilep->SetPoint( 2, 162 ,  (0.00296467 / 0.00321705 )  );// 162     // 163
  SFdilep->SetPoint( 3, 242 ,  (0.000701592/ 0.000931674)  );// 242     // 247
  SFdilep->SetPoint( 4, 343 ,  (0.00012036 / 0.000191065)  );// 343     // 350
  //                   x errors   rel.err(data) *( data  / Madgraph)
  SFdilep->SetPointError( 0, 0., 0., 0.0601381*(0.00509572 / 0.00453114 ), 0.0601381*(0.00509572 / 0.00453114 ) );
  SFdilep->SetPointError( 1, 0., 0., 0.0469906*(0.00626002 / 0.00600115 ), 0.0469906*(0.00626002 / 0.00600115 ) );
  SFdilep->SetPointError( 2, 0., 0., 0.0555114*(0.00296467 / 0.00321705 ), 0.0555114*(0.00296467 / 0.00321705 ) );
  SFdilep->SetPointError( 3, 0., 0., 0.071274* (0.000701592/ 0.000931674), 0.071274* (0.000701592/ 0.000931674) );
  SFdilep->SetPointError( 4, 0., 0., 0.0924826*(0.00012036 / 0.000191065), 0.0924826*(0.00012036 / 0.000191065) );
  //style of ratio
  SFdilep->SetLineWidth(3.);
  SFdilep->SetMarkerSize(1.5);
  SFdilep->SetMarkerStyle(22);
  SFdilep->SetMarkerColor(dilep7color);
  SFdilep->SetLineColor(dilep7color);

  // collect 8 TeV BCC x values for analysis binning
  std::vector<double> xBCCljets_;
  xBCCljets_.push_back( 28  );   //   0.0 ..  60.0
  xBCCljets_.push_back( 86  );   //  60.0 .. 100.0
  xBCCljets_.push_back(125  );   // 100.0 .. 150.0
  xBCCljets_.push_back(173  );   // 150.0 .. 200.0
  xBCCljets_.push_back(227.3);   // 200.0 .. 260.0
  xBCCljets_.push_back(288  );   // 260.0 .. 320.0
  xBCCljets_.push_back(356  );   // 320.0 .. 400.0
  xBCCljets_.push_back(444  );   // 400.0 .. 500.0
  std::vector<double> xBCCdilep_;
  xBCCdilep_.push_back( 29.7);  //   0.0 ..  65.0
  xBCCdilep_.push_back( 99.6);  //  65.0 .. 125.0
  xBCCdilep_.push_back(159.7);  // 125.0 .. 200.0
  xBCCdilep_.push_back(239.1);  // 200.0 .. 290.0
  xBCCdilep_.push_back(336.2);  // 290.0 .. 400.0

  // c) l+jets 8TeV data points
  for(int p=0; p<NbinsLjets8; ++p){
    // get line with all informations
    TString line= readLineFromFile(p+1, "/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/topPtInputForReweighting/diffXSecTopSemiLepPartontopPt.txt");
    // data value
    TString temp = getStringEntry(line, 3 , "&");
    temp.ReplaceAll(" ","");
    double data=atof(temp.Data());
    temp = getStringEntry(line, 2 , "&");
    temp.ReplaceAll(" ","");
    double MC  =atof(temp.Data());
    SFljets8->SetPoint( p,  xBCCljets_.at(p) , data/MC ); 
    temp = getStringEntry(line, 6 , "&");
    double unc=atof(temp.Data());
    SFljets8->SetPointError( p, 0., 0., (unc/100.)*(data/MC), (unc /100.)*(data/MC) );
  }
  whipEmptyBinsAway(SFljets8, 0);

  //style of ratio
  SFljets8->SetLineWidth(3.);
  SFljets8->SetMarkerSize(1.5);
  SFljets8->SetMarkerStyle(24);
  //SFljets8->SetLineStyle(2);
  SFljets8->SetMarkerColor(ljets8color);
  SFljets8->SetLineColor(ljets8color);

  // d) dilepton 8TeV data points
  // MC prediction point (as not in provided table)
  std::vector<double> MCdilep_;
  MCdilep_.push_back(0.00396076 ); 
  MCdilep_.push_back(0.00620269 ); 
  MCdilep_.push_back(0.00336987 ); 
  MCdilep_.push_back(0.00102834 ); 
  MCdilep_.push_back(0.000228163); 
  for(int p=0; p<NbinsDilep8; ++p){
    // get line with all informations
    TString line= readLineFromFile(p+4, "/afs/naf.desy.de/group/cms/scratch/tophh/CommonFiles/topPtInputForReweighting/HypToppTLaTeX.txt");
    // data value
    TString temp = getStringEntry(line, 3 , "&");
    temp.ReplaceAll(" ","");
    double data=atof(temp.Data());
    //temp = getStringEntry(line, 2 , "&");
    //temp.ReplaceAll(" ","");
    double MC  =MCdilep_[p];
    SFdilep8->SetPoint( p,  xBCCdilep_.at(p) , data/MC ); 
    temp = getStringEntry(line, 6 , "&");
    double unc=atof(temp.Data());
    SFdilep8->SetPointError( p, 0., 0., (unc/100.)*(data/MC), (unc /100.)*(data/MC) );
  }

  //style of ratio
  SFdilep8->SetLineWidth(3.);
  SFdilep8->SetMarkerSize(1.5);
  SFdilep8->SetMarkerStyle(20);
  SFdilep8->SetMarkerColor(dilep8color);
  SFdilep8->SetLineColor(dilep8color);

  // e) combined 7 TeV data points
  addTAE(SFdilep , SF7);
  addTAE(SFljets , SF7);
  //style of ratio
  SF7->SetLineWidth(3.);
  SF7->SetMarkerSize(0.1);
  SF7->SetMarkerStyle(20);
  SF7->SetMarkerColor(kWhite);
  SF7->SetLineColor(kWhite);

  // f) combined 8 TeV data points
  addTAE(SFdilep8, SF8);
  addTAE(SFljets8, SF8);
  //style of ratio
  SF8->SetLineWidth(3.);
  SF8->SetMarkerSize(0.1);
  SF8->SetMarkerStyle(20);
  SF8->SetMarkerColor(kWhite);
  SF8->SetLineColor(kWhite);

  // g) combined 7+8TeV data points
  addTAE(SF7, SF);
  addTAE(SF8, SF);
  //style of ratio
  SF->SetLineWidth(3.);
  SF->SetMarkerSize(0.1);
  SF->SetMarkerStyle(20);
  SF->SetMarkerColor(kWhite);
  SF->SetLineColor(kWhite);

  // ---
  //    dummy plots for axis
  // ---
  TH1F* dummy= new TH1F("","",1,0.,xmax);
  histogramStyle(*dummy, kSig);
  dummy->GetXaxis()->SetTitle("p_{T}^{t} [GeV]");
  dummy->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{t}} Ratio: (Data / Simulation)");
  dummy->GetYaxis()->SetTitleOffset(0.9*dummy->GetYaxis()->GetTitleOffset());
  dummy->SetMaximum(1.8);
  dummy->SetMinimum(0.5);

  // ---
  //    legends
  // ---
  double x1=0.29;
  double x2=0.86;
    
  TLegend *leg0 = new TLegend(x1, 0.69, x2, 0.87);
  leg0->SetFillStyle(0);
  leg0->SetTextSize(0.035);
  leg0->SetBorderSize(0);
  leg0->SetHeader("#font[22]{Data / MadGraph+PYTHIA(CTEQ6L1)}");

  TLegend *leg1 = new TLegend(x1, 0.57, x2, 0.69);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->SetBorderSize(0);
  leg1->SetHeader("#font[22]{Fit: exp(a+b#upointx)}");

  if(plotsepfit) leg1->SetY1(leg1->GetY1()-0.2);

  // canvas
  std::vector<TCanvas*> plotCanvas_;
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle("data/MC top Pt ratio");
  // drawing
  dummy->Draw("axis");
  SF->Draw("p e1 same");
  SF7->Draw("p e1 same");
  SF8->Draw("p e1 same");
  SFljets->Draw("p e1 same");
  SFdilep->Draw("p e1 same");
  SFljets8->Draw("p e1 same");
  SFdilep8->Draw("p e1 same");
  // fit polynomial or exponential function
  TString def = "";
  if(func=="pol2")def="[0]*x*x+[1]*x+[2]";
  if(func=="exp" )def="exp([0]+[1]*x)";
  double fitLowEdge=0.;
  double fitHighEdge=xmax;
  // a) to all 8 and 7 TeV points
//   TF1* function=new TF1("function",def,fitLowEdge, fitHighEdge);
//   function->SetLineColor(kMagenta+2);
//   SF->Fit(function,"R","same",fitLowEdge, fitHighEdge);
//   for(int i=0; i<function->GetNumberFreeParameters(); i++){
//     function->SetParameter(i,round(function->GetParameter(i),3));
//   }
//   TString fitEntry="#splitline{}{#splitline{}{#splitline{}{#splitline{combined fit: ";
//   fitEntry+=function->GetExpFormula("p")+",}{                          #chi^{2}/ndof=";
//   fitEntry+=getTStringFromDouble(function->GetChisquare())+"/"+getTStringFromInt(function->GetNDF())+"}}}}";
//   fitEntry.ReplaceAll("+(","");
//   fitEntry.ReplaceAll("))",")");
//   leg0->AddEntry( function, fitEntry, "L");
  // b) to all 7 TeV points
  TF1* function7=new TF1("function7",def,fitLowEdge, fitHighEdge);
  function7->SetLineColor(fit7color);
  function7->SetLineWidth(6);
  function7->SetLineStyle(2);
  SF7->Fit(function7,"R","same",fitLowEdge, fitHighEdge);
  for(int i=0; i<function7->GetNumberFreeParameters(); i++){
    function7->SetParameter(i,round(function7->GetParameter(i),3));
  }
  //TString fitEntry7="fit 7 TeV: ";
  //fitEntry7+=function7->GetExpFormula("p");
  //fitEntry7+=",  #chi^{2}/ndof=";
  //fitEntry7+=getTStringFromDouble(function7->GetChisquare())+"/"+getTStringFromInt(function7->GetNDF());
  //fitEntry7.ReplaceAll("+(","");
  //fitEntry7.ReplaceAll("))",")");
  TString fitEntry7="7 TeV: ";
  if(plotfiterrors) fitEntry7+="              ";
  fitEntry7+="a=";
  fitEntry7+=getTStringFromDouble(function7->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntry7+="#pm";
    fitEntry7+=getTStringFromDouble(function7->GetParError(0) , 3);
  }
  fitEntry7+=", b=";
  fitEntry7+=getTStringFromDouble(function7->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntry7+="#pm";
    fitEntry7+=getTStringFromDouble(function7->GetParError(1) , 5);
  }

  // b1) to l+jets 7 TeV points
  TF1* functionljets7=new TF1("functionljets7",def,fitLowEdge, fitHighEdge);
  functionljets7->SetLineColor(kRed+1);
  functionljets7->SetLineWidth(2);
  SFljets->Fit(functionljets7,"R"+optD,"same",fitLowEdge, fitHighEdge);
  for(int i=0; i<functionljets7->GetNumberFreeParameters(); i++){
    functionljets7->SetParameter(i,round(functionljets7->GetParameter(i),3));
  }
  //TString fitEntryljets7="fit 7 TeV l+jets: ";
  //fitEntryljets7+=functionljets7->GetExpFormula("p");
  //fitEntryljets7+=",  #chi^{2}/ndof=";
  //fitEntryljets7+=getTStringFromDouble(functionljets7->GetChisquare())+"/"+getTStringFromInt(functionljets7->GetNDF());
  //fitEntryljets7.ReplaceAll("+(","");
  //fitEntryljets7.ReplaceAll("))",")");
  TString fitEntryljets7="7 TeV l+jets:     ";
  fitEntryljets7+="a=";
  fitEntryljets7+=getTStringFromDouble(functionljets7->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntryljets7+="#pm";
    fitEntryljets7+=getTStringFromDouble(functionljets7->GetParError(0) , 3);
  }
  fitEntryljets7+=", b=";
  fitEntryljets7+=getTStringFromDouble(functionljets7->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntryljets7+="#pm";
    fitEntryljets7+=getTStringFromDouble(functionljets7->GetParError(1) , 5);
  }

  // b2) to dilepton 7 TeV points
  TF1* functiondilep7=new TF1("functiondilep7",def,fitLowEdge, fitHighEdge);
  functiondilep7->SetLineColor(kOrange+7);
  functiondilep7->SetLineWidth(2);
  SFdilep->Fit(functiondilep7,"R"+optD,"same",fitLowEdge, fitHighEdge);
  for(int i=0; i<functiondilep7->GetNumberFreeParameters(); i++){
    functiondilep7->SetParameter(i,round(functiondilep7->GetParameter(i),3));
  }
  //TString fitEntrydilep7="fit 7 TeV dilepton: ";
  //fitEntrydilep7+=functiondilep7->GetExpFormula("p");
  //fitEntrydilep7+=",  #chi^{2}/ndof=";
  //fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetChisquare())+"/"+getTStringFromInt(functiondilep7->GetNDF());
  //fitEntrydilep7.ReplaceAll("+(","");
  //fitEntrydilep7.ReplaceAll("))",")");
  TString fitEntrydilep7="7 TeV dilepton: ";
  fitEntrydilep7+="a=";
  fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntrydilep7+="#pm";
    fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParError(0) , 3);
  }
  fitEntrydilep7+=", b=";
  if(plotfiterrors){
    fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParameter(1), 5);
    fitEntrydilep7+="#pm";
  }
  fitEntrydilep7+=getTStringFromDouble(functiondilep7->GetParError(1) , 5);

  // c) to all 8 TeV points
  TF1* function8=new TF1("function8",def,fitLowEdge, fitHighEdge);
  function8->SetLineWidth(6);
  function8->SetLineColor(fit8color);
  function8->SetLineStyle(2);
  SF8->Fit(function8,"R","same",fitLowEdge, fitHighEdge);
  for(int i=0; i<function8->GetNumberFreeParameters(); i++){
    function8->SetParameter(i,round(function8->GetParameter(i),3));
  }
  //TString fitEntry8="fit 8 TeV: ";
  //fitEntry8+=function8->GetExpFormula("p");
  //fitEntry8+=",  #chi^{2}/ndof=";
  //fitEntry8+=getTStringFromDouble(function8->GetChisquare())+"/"+getTStringFromInt(function8->GetNDF());
  //fitEntry8.ReplaceAll("+(","");
  //fitEntry8.ReplaceAll("))",")");
  TString fitEntry8="8 TeV: ";
  if(plotfiterrors) fitEntry8+="              ";
  fitEntry8+="a=";
  fitEntry8+=getTStringFromDouble(function8->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntry8+="#pm";
    fitEntry8+=getTStringFromDouble(function8->GetParError(0) , 3);
  }
  fitEntry8+=", b=";
  fitEntry8+=getTStringFromDouble(function8->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntry8+="#pm";
    fitEntry8+=getTStringFromDouble(function8->GetParError(1) , 5);
  }

  // c1) to l+jets 8 TeV points
  TF1* functionljets8=new TF1("functionljets8",def,fitLowEdge, fitHighEdge);
  functionljets8->SetLineColor(kBlue);
  functionljets8->SetLineWidth(2);
  SFljets8->Fit(functionljets8,"R"+optD,"same",fitLowEdge, fitHighEdge);
  for(int i=0; i<functionljets8->GetNumberFreeParameters(); i++){
    functionljets8->SetParameter(i,round(functionljets8->GetParameter(i),3));
  }
  //TString fitEntryljets8="fit 8 TeV l+jets: ";
  //fitEntryljets8+=functionljets8->GetExpFormula("p");
  //fitEntryljets8+=",  #chi^{2}/ndof=";
  //fitEntryljets8+=getTStringFromDouble(functionljets8->GetChisquare())+"/"+getTStringFromInt(functionljets8->GetNDF());
  //fitEntryljets8.ReplaceAll("+(","");
  //fitEntryljets8.ReplaceAll("))",")");
  TString fitEntryljets8="8 TeV l+jets:     ";
  fitEntryljets8+="a=";
  fitEntryljets8+=getTStringFromDouble(functionljets8->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntryljets8+="#pm";
    fitEntryljets8+=getTStringFromDouble(functionljets8->GetParError(0) , 3);
  }
  fitEntryljets8+=", b=";
  fitEntryljets8+=getTStringFromDouble(functionljets8->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntryljets8+="#pm";
    fitEntryljets8+=getTStringFromDouble(functionljets8->GetParError(1) , 5);
  }

  // c2) to dilepton 8 TeV points
  TF1* functiondilep8=new TF1("functiondilep8",def,fitLowEdge, fitHighEdge);
  functiondilep8->SetLineColor(kAzure+6);
  functiondilep8->SetLineWidth(2);
  
  SFdilep8->Fit(functiondilep8,"R"+optD,"same",fitLowEdge, fitHighEdge);
  for(int i=0; i<functiondilep8->GetNumberFreeParameters(); i++){
    functiondilep8->SetParameter(i,round(functiondilep8->GetParameter(i),3));
  }
  //TString fitEntrydilep8="fit 8 TeV dilepton: ";
  //fitEntrydilep8+=functiondilep8->GetExpFormula("p");
  //fitEntrydilep8+=",  #chi^{2}/ndof=";
  //fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetChisquare())+"/"+getTStringFromInt(functiondilep8->GetNDF());
  //fitEntrydilep8.ReplaceAll("+(","");
  //fitEntrydilep8.ReplaceAll("))",")");
  TString fitEntrydilep8="8 TeV dilepton: ";
  fitEntrydilep8+="a=";
  fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParameter(0), 3);
  if(plotfiterrors){
    fitEntrydilep8+="#pm";
    fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParError(0) , 3);
  }
  fitEntrydilep8+=", b=";
  fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParameter(1), 5);
  if(plotfiterrors){
    fitEntrydilep8+="#pm";
    fitEntrydilep8+=getTStringFromDouble(functiondilep8->GetParError(1) , 5);
  }

  // Draw legend
  leg0->AddEntry(SFljets, "7 TeV e/#mu+jets  (TOP-11-013)"   , "P");
  leg0->AddEntry(SFdilep, "7 TeV ee/e#mu/#mu#mu (TOP-11-013)", "P");
  leg0->AddEntry(SFljets8,"8 TeV e/#mu+jets  (TOP-12-028)"   , "P");
  leg0->AddEntry(SFdilep8,"8 TeV ee/e#mu/#mu#mu (TOP-12-028)", "P");
  leg0->Draw("same");
  leg1->AddEntry( function7, fitEntry7, "L");
  if(plotsepfit) leg1->AddEntry( functiondilep7, fitEntrydilep7, "L");
  if(plotsepfit) leg1->AddEntry( functionljets7, fitEntryljets7, "L");
  leg1->AddEntry( function8, fitEntry8, "L");
  if(plotsepfit) leg1->AddEntry( functiondilep8, fitEntrydilep8, "L");
  if(plotsepfit) leg1->AddEntry( functionljets8, fitEntryljets8, "L");
  leg1->Draw("same");
  // Draw cms label
  TPaveText *label = new TPaveText();
  label -> SetX1NDC(gStyle->GetPadLeftMargin());
  label -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label -> SetY2NDC(1.0);
  label -> SetTextFont(42);
  TString CMSlab="";
  if(!PHD) CMSlab+="CMS Preliminary, ";  
  CMSlab+="5.0/19.7 fb^{-1} at #sqrt{s} = 7/8 TeV";
  label -> AddText(CMSlab);
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.04);
  label->SetTextAlign(32);
  label-> Draw("same");
  // BCC label
  double positionX=xmax+0.045*xmax*(gStyle->GetCanvasDefW()/600.);
  double positionY=0.5;
  TLatex *bcclabel = new TLatex(positionX,positionY, " (horizontal BCC wrt. NNLO^{approx}, arXiv:1205.3453)");
  bcclabel->SetTextAlign(11);
  bcclabel->SetTextAngle(90);
  bcclabel->SetTextSize(0.035);
  bcclabel->Draw("same");
  if(grey) plotCanvas_[0]->SetGrayscale();
  //saving
  plotCanvas_[0]->Print("diffXSecFromSignal/plots/combined/2012/xSec/dataVsMadgraph7and8TeV.eps");
  plotCanvas_[0]->Print("diffXSecFromSignal/plots/combined/2012/xSec/dataVsMadgraph7and8TeV.png");

  // ---
  // ERROR band plot
  // ---
  addCanvas(plotCanvas_);
  plotCanvas_[plotCanvas_.size()-1]->cd(0);
  plotCanvas_[plotCanvas_.size()-1]->SetTitle("data/MC top Pt errorband");
  // drawing
  dummy->GetYaxis()->SetTitle("(Data / Simulation) SF (#sqrt{s}=8 TeV)");
  dummy->GetYaxis()->SetRangeUser(0.35, 1.65);
  dummy->SetFillColor(10);
  dummy->SetLineColor(10);
  dummy->Draw("axis");
  // extract parameters
  double a=function8->GetParameter(0);
  double b=function8->GetParameter(1);
  // turning point
  double min=0;
  double max=500;
  double TP=-a/b;
  // get functions for high, low and central
  TF1* centralErr=new TF1("centralErr",def,min, max);
  centralErr->SetParameter(0, a);
  centralErr->SetParameter(1, b);
  TF1* upErr=new TF1("upErr",def,min, max);
  upErr->SetParameter(0, 2*a);
  upErr->SetParameter(1, 2*b);
  TF1* dnErr=new TF1("upErr",def,min, max);
  dnErr->SetParameter(0, 0.);
  dnErr->SetParameter(1, 0.);
  // draw errorbands
  upErr->SetFillStyle(1001);
  dnErr->SetFillStyle(1001);
  upErr->SetLineColor(10);
  dnErr->SetLineColor(10);
  upErr->SetFillColor(colorband);
  upErr->SetRange(min,TP);
  upErr->DrawClone("hist same");
  dnErr->SetFillColor(10);
  dnErr->SetLineColor(10);
  dnErr->SetRange(min,TP);
  dnErr->DrawClone("hist same");
  dnErr->SetFillColor(colorband);
  dnErr->SetLineColor(colorband);
  dnErr->SetRange(TP, max);
  dnErr->DrawClone("hist same");
  upErr->SetFillColor(10);
  upErr->SetLineColor(10);
  upErr->SetRange(TP, max);
  upErr->DrawClone("hist same");
  drawLine(TP, 0.35, TP, 1.05, 10, 2, 1);
  // draw central prediction
  centralErr->SetFillStyle(0);
  centralErr->SetFillColor(0);
  centralErr->SetLineColor(kBlue);
  centralErr->SetLineWidth(6);
  centralErr->SetLineColor(fit8color);
  centralErr->SetLineStyle(2);
  centralErr->Draw("hist same");
  // legend and labels
  dummy->Draw("axis same");
  TPaveText *label2 = new TPaveText();
  label2 -> SetX1NDC(gStyle->GetPadLeftMargin());
  label2 -> SetY1NDC(1.0-gStyle->GetPadTopMargin());
  label2 -> SetX2NDC(1.0-gStyle->GetPadRightMargin());
  label2 -> SetY2NDC(1.0);
  label2 -> SetTextFont(42);
  TString CMSlab2="";
  if(!PHD) CMSlab2+="CMS Preliminary, ";  
  CMSlab2+="19.7 fb^{-1} at #sqrt{s} = 8 TeV";
  label2->AddText(CMSlab2);
  label2->SetFillStyle(0);
  label2->SetBorderSize(0);
  label2->SetTextSize(0.04);
  label2->SetTextAlign(32);
  label2->Draw("same");
  TLegend *leg3 = new TLegend(x1+0.05, 0.7, x2+0.05, 0.85);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.035);
  leg3->SetBorderSize(0);
  leg3->SetHeader("#font[22]{Parametrisation: exp(a+b#upointx)}");
  TString entryErr=fitEntry8;
  entryErr.ReplaceAll("8 TeV: ", "");
  leg3->AddEntry(centralErr, entryErr , "L");
  leg3->AddEntry(dnErr     , "a,b #pm 100%", "F");
  leg3->Draw("same");
  //saving
  if(grey) plotCanvas_[1]->SetGrayscale();
  plotCanvas_[1]->Print("diffXSecFromSignal/plots/combined/2012/xSec/topPtReweighting8TeVunc.eps");
  plotCanvas_[1]->Print("diffXSecFromSignal/plots/combined/2012/xSec/topPtReweighting8TeVunc.png");

}


void addTAE(TGraphAsymmErrors* in, TGraphAsymmErrors* out){
  // this functions adds the points of "in" to out

  // calculate new number of points
  int Nin =in ->GetMaxSize();
  int Nout=out->GetMaxSize();
  int Nsum=Nin+Nout;
  out->Expand(Nsum);

  // add points
  for(int p=0; p<Nin; ++p){
    out->SetPoint     ( Nout+p, in->GetX()[p], in->GetY()[p] );
    out->SetPointError( Nout+p, in->GetEXlow()[p], in->GetEXhigh()[p], in->GetEYlow()[p], in->GetEYhigh()[p] );
  }
}
