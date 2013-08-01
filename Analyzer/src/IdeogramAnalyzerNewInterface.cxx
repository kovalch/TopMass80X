#include "IdeogramAnalyzerNewInterface.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <iomanip>

#include "TCanvas.h"
//#include "TColor.h"
#include "TF2.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

typedef ProgramOptionsReader po;

IdeogramAnalyzerNewInterface::IdeogramAnalyzerNewInterface(const TString& identifier, TTree* tree) :
    MassAnalyzer(identifier, tree),
    fptr_(0),
    combLikelihood_(0),
    channelID_(Helper::channelID())
{
}

void IdeogramAnalyzerNewInterface::Analyze(const TString& cuts, int i, int j) {
  Scan(cuts, i, j, 154, 190, 2, 0.9, 1.1, 0.02);

  double mass = GetValue("mass_mTop_JES").first;
  double JES  = GetValue("JES_mTop_JES" ).first;
  Scan(cuts, i, j, mass-2 , mass+2 , 0.1, JES-0.015, JES+0.015, 0.00075);

  double epsilon = 1e-6;
  mass = GetValue("mass_mTop_JES").first;
  JES  = GetValue("JES_mTop_JES" ).first;
  Scan(cuts, i, j, mass-10, mass+10, 0.1 , 1.-epsilon, 1.+epsilon, epsilon, false);
}


void IdeogramAnalyzerNewInterface::Scan(const TString& cuts, int i, int j, double firstBinMass, double lastBinMass,
			    double resolMass, double firstBinJes, double lastBinJes, double resolJes, bool fit2D)
{
  //*
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.185);
  //*/
  
  /*
  Double_t r[]    = {0., 0.0, 1.0, 1.0, 1.0};
  Double_t g[]    = {0., 0.0, 0.0, 1.0, 1.0};
  Double_t b[]    = {0., 1.0, 0.0, 0.0, 1.0};
  Double_t stop[] = {0., .25, .50, .75, 1.0};
  TColor::CreateGradientColorTable(5, stop, r, g, b, 100);
  //*/
  
  bool blackWhite = false;
  bool syst       = false;
  //bool cmsPrel    = true;
  
  bool debug = false;
  int nDebug = 1;
  int minDebug = 0;
  int maxDebug = 200;
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int binsMass = int((lastBinMass-firstBinMass)/resolMass);
  int binsJes  = int((lastBinJes -firstBinJes )/resolJes );
  
  double pullWidth = 1.;
  
  /*
  if (debug) {
    firstBinMass = 100;
    lastBinMass  = 350;
    binsMass     = 125;
    
    firstBinJes = 0.7;
    lastBinJes  = 1.3;
    binsJes     = 60;
  }
  //*/
  
  if(!combLikelihood_ || !fptr_){
    if(channelID_ == Helper::kAllJets){
      IdeogramCombLikelihoodAllJets *fptrAllJets = new IdeogramCombLikelihoodAllJets();
      fptr_ = fptrAllJets;
      combLikelihood_ = new TF2("combLikelihood",fptrAllJets,&IdeogramCombLikelihoodAllJets::Evaluate, firstBinMass, lastBinMass, firstBinJes, lastBinJes, 7, "IdeogramCombLikelihoodAllJets", "Evaluate");
    }
  }
  //TF2* gausJESConstraint = new TF2("gausJESConstraint", "x*0 +((y-1.)/0.013)**2", firstBinMass, lastBinMass, firstBinJes,lastBinJes);

  //TF1* combBackground = new TF1("combBackground",fptr,&IdeogramCombLikelihood::CrystalBall,150,200,1);

  TF1* fitParabola = new TF1("fitParabola", "abs([1])*(x-[0])^2+[2]");
  fitParabola->SetParNames("mass", "massCurv", "minL");
  fitParabola->SetParLimits(0, firstBinMass, lastBinMass);
  fitParabola->SetParLimits(1, 0.0001, 1000);
  fitParabola->SetLineColor(kRed+1);
  
  TF2* fitParaboloid  = new TF2("fitParaboloid" , "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2 + [5]");
  TF2* systParaboloid = new TF2("systParaboloid", "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2 + [5]");
  fitParaboloid->SetNpx(300);
  fitParaboloid->SetNpy(300);
  fitParaboloid->SetParNames("mass", "massCurv", "jes", "jesCurv", "alpha", "minL");
  fitParaboloid->SetParLimits(0, firstBinMass, lastBinMass);
  fitParaboloid->SetParLimits(1, 0.0001, 1000);
  fitParaboloid->SetParLimits(2, 0.5, 1.5);
  fitParaboloid->SetParLimits(3, 1, 10000000);
  fitParaboloid->SetParLimits(4, 0, 6.3);
  fitParaboloid->SetParameter(4, 1.5);
  if (blackWhite) fitParaboloid->SetLineColor(kRed+1);
  else fitParaboloid->SetLineColor(kWhite);
  fitParaboloid->SetLineWidth(3);

  TF2* null  = new TF2("null" , "0 + 0*x + 0*y");
  TF2* unity = new TF2("unity", "1 + 0*x + 0*y");
  TH2D* hUnity = new TH2D("hUnity","hUnity", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  hUnity->Eval(unity);

  TH2D* eventLikelihood    = new TH2D("eventLikelihood"   , "eventLikelihood"   , binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* logEventLikelihood = new TH2D("logEventLikelihood", "logEventLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* sumLogLikelihood   = new TH2D("sumLogLikelihood"  , "sumLogLikelihood"  , binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* productLikelihood  = new TH2D("productLikelihood" , "productLikelihood" , binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  sumLogLikelihood->Eval(null);
  
  eventLikelihood->SetTitle("L(m_{t}|event)");
  eventLikelihood->SetXTitle("m_{t} [GeV]");
  eventLikelihood->SetYTitle("JES");
  logEventLikelihood->SetTitle("-2#upointln{L(m_{t}|event)}");
  logEventLikelihood->SetXTitle("m_{t} [GeV]");
  logEventLikelihood->SetYTitle("JES");
  sumLogLikelihood->SetTitle("-2#upointln{L(m_{t}|sample)}");
  sumLogLikelihood->SetXTitle("m_{t} [GeV]");
  sumLogLikelihood->SetYTitle("JES");
  sumLogLikelihood->SetZTitle("-2 #Delta log(L)");
  productLikelihood->SetTitle("-2#upointln{L(m_{t}|sample)}");
  productLikelihood->SetXTitle("m_{t} [GeV]");
  productLikelihood->SetYTitle("JES");

  TopEvent    *topEvent    = new TopEvent();
  WeightEvent *weightEvent = new WeightEvent();

  double fitWeight;
  int nEvents = 0;
  double sumWeights = 0.;

  std::cout << "fTree: " << fTree_->GetEntries() << std::endl;

  fTree_->SetBranchAddress("top", &topEvent);
  fTree_->SetBranchAddress("weight", &weightEvent);

  double isFastSim                     = po::GetOption<int   >("fastsim");
  double shapeSystematic               = po::GetOption<double>("shape"  );
  double permutationFractionSystematic = po::GetOption<double>("permu"  );

  TString plotPath("plot/Ideogram/");
  // Build Likelihood
  for (int iEntry = 0, length = fTree_->GetEntries(); iEntry < length; ++iEntry) {
    topEvent->init();
    weightEvent->init();
    fTree_->GetEntry(iEntry);

    //if (event == currentEvent) continue;
    //currentEvent = event;
    ++nEvents;
    if ((debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) || iEntry%1000 == 0) std::cout << iEntry << " - " << topEvent->event() << std::endl;
    
    //std::cout << "eventLikelihood: " << eventLikelihood->GetEntries() << std::endl;
    eventLikelihood->Eval(null);
    eventLikelihood->SetFillColor(0);
    //weight = 0;
    fitWeight = 0;
    //currentWeight = 0;
    
    if (debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) {
    std::cout << std::setiosflags(std::ios::left)
              << std::setw(04) << "n"
              << std::setw(10) << "mt"
              << std::setw(10) << "mW1"
              << std::setw(10) << "mW2"
              << std::setw(12) << "fitProb"
              << std::setw(11) << "dRbb"
      //<< std::setw(11) << "weight"
              << std::endl;
    }
    
    //for (int iComb = 0; iComb < 24; iComb++) {
    //  if (fTree_->GetEntries() < iEntry + iComb + 1) break;
    //  fTree_->GetEntry(iEntry + iComb);
    //  
    //  if (event != currentEvent) break;
    fitWeight += topEvent->fitProb()[0]; //*CombinedWeight;
      
    if (debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) {
      std::cout << std::setw(04) << topEvent->fitProb().size()
		<< std::setw(10) << topEvent->fitTop1()[0].M()
		<< std::setw(10) << topEvent->recoW1()[0].M()
		<< std::setw(10) << topEvent->recoW2()[0].M()
		<< std::setw(12) << topEvent->fitProb()[0]
		<< std::setw(11) << topEvent->fitB1()[0].DeltaR(topEvent->fitB2()[0])
		<< std::endl;
    }
      
    if (topEvent->fitProb()[0] != 0) {
      //bScaleEstimator = 1;

      // Set Likelihood parameters
      combLikelihood_->SetParameters(topEvent->fitProb()[0], topEvent->fitTop1()[0].M(), (topEvent->recoW1()[0].M()+topEvent->recoW2()[0].M())/2.0, 1., shapeSystematic, permutationFractionSystematic, isFastSim);
      // add permutation to event likelihood
      eventLikelihood->Eval(combLikelihood_, "A");
    }
    //}
    
    eventLikelihood->Scale(1./fitWeight);
    sumWeights += fitWeight;
    //if (weight == 0) continue;
    
    logEventLikelihood->Eval(null);

    for (int i = 0; i<=binsMass; i++) {
      for (int j = 0; j<=binsJes; j++) {
        logEventLikelihood->SetBinContent(i, j, -2*TMath::Log(eventLikelihood->GetBinContent(i, j)));
      }
    }

    //TString sEvent("(run=="); sEvent += run; sEvent += " & luminosityBlock=="; 
    //sEvent += luminosityBlock; sEvent += " & event=="; sEvent += event; sEvent += ")";
    
    //TString sEventWeighted = sEvent; sEventWeighted += "*("; sEventWeighted += "CombinedWeight"; sEventWeighted += ")";
    
    //double eventWeight = fTree_->GetEntries(sEventWeighted)/fTree_->GetEntries(sEvent);

    sumLogLikelihood->Add(logEventLikelihood, fitWeight/(pullWidth*pullWidth)); // add weight here
    
    if (debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) {
      TCanvas* eventCanvas = new TCanvas("eventCanvas", "eventCanvas", 1200, 400);
      eventCanvas->Divide(3, 1);
      
      eventCanvas->cd(1);
      eventLikelihood->SetEntries(1);
      eventLikelihood->Draw("COLZ");
      eventLikelihood->SetFillColor(kRed+1);
      /*
      if (weight > 1./16) eventLikelihood->SetFillColor(kGreen);
      if (weight > 1./8) eventLikelihood->SetFillColor(kYellow);
      if (weight > 1./4) eventLikelihood->SetFillColor(kOrange);
      if (weight > 1./2) eventLikelihood->SetFillColor(kRed);
      */
      
      eventCanvas->cd(2);
      eventLikelihood->SetEntries(1);
      logEventLikelihood->Draw("COLZ");
      
      eventCanvas->cd(3);
      eventLikelihood->SetEntries(1);
      sumLogLikelihood->Draw("COLZ");
      
      TString eventPath(plotPath); eventPath += fIdentifier_; eventPath += "_"; eventPath += iEntry; eventPath += "_"; eventPath += topEvent->event(); eventPath += ".eps";
      std::cout << eventPath << std::endl;
      eventCanvas->Print(eventPath);
      
      delete eventCanvas;
    }
  }
  
  ctemp->cd();
  
  std::cout << "Finishing..." << std::endl;
  
  std::cout << "Sum of weights: " << sumWeights << std::endl;
  std::cout << "Total number of events: " << nEvents << std::endl;
  sumLogLikelihood->Scale(nEvents/sumWeights);

  /* JES constraint
  eventLikelihood->Eval(gausJESConstraint);
  for (int i = 0; i<=binsMass; i++) {
    for (int j = 0; j<=binsJes; j++) {
      sumLogLikelihood->SetBinContent(i, j, sumLogLikelihood->GetBinContent(i, j) + eventLikelihood->GetBinContent(i, j));
    }
  }
  //*/
  
  int minBinX;
  int minBinY;
  int minBinZ;
  
  sumLogLikelihood->GetMinimumBin(minBinX, minBinY, minBinZ);
  
  double minMass = sumLogLikelihood->GetXaxis()->GetBinCenter(minBinX);
  double minJes  = sumLogLikelihood->GetYaxis()->GetBinCenter(minBinY);
  
  std::cout << "minMass: " << minMass << ", minJes: " << minJes << std::endl;
  
  /*
  sumLogLikelihood->SetAxisRange(minMass - 3, minMass + 3, "X");
  sumLogLikelihood->SetAxisRange(minJes - 0.03, minJes + 0.03, "Y");
  sumLogLikelihood->SetAxisRange(0, 20, "Z");
  //*/
  
  //sumLogLikelihood->SetMarkerStyle(20);
  //sumLogLikelihood->SetMarkerColor(kRed+1);
  
  std::cout << "Minimum likelihood: " << sumLogLikelihood->GetMinimum(0) << "\tMaximum likelihood (in range): " << sumLogLikelihood->GetMaximum() << std::endl;
  
  Helper* helper = new Helper();
  
  if (fit2D) {  
  
    fitParaboloid->SetParLimits(0, minMass-2*resolMass, minMass+2*resolMass);
    fitParaboloid->SetParameter(0, minMass);
    fitParaboloid->SetParLimits(2, minJes-2*resolJes, minJes+2*resolJes);
    fitParaboloid->SetParameter(2, minJes);
    fitParaboloid->SetParameter(3, 1000000);
    fitParaboloid->SetParLimits(5, sumLogLikelihood->GetMinimum(0)-10., sumLogLikelihood->GetMinimum(0)+10.);
    fitParaboloid->SetParameter(5, sumLogLikelihood->GetMinimum(0));
  
    //fitParaboloid->SetRange(minMass - 1, minJes - 0.01, minMass + 1, minJes + 0.01);
    fitParaboloid->SetRange(minMass - 4, minJes - 0.04, minMass + 4, minJes + 0.04);
    //fitParaboloid->SetRange(minMass - 20, minJes - 0.2, minMass + 20, minJes + 0.2);

    //std::cout << "sumLogLikelihood: " << sumLogLikelihood << std::endl;
    sumLogLikelihood->Fit("fitParaboloid","EMR0");
  
    double semiMajor, semiMinor, alpha = 0;

    double mass = -1;
    double JES  = -1;
    double massError = -1;
    double JESError  = -1;
    if (firstBinMass+1 < fitParaboloid->GetParameter(0) && fitParaboloid->GetParameter(0) < lastBinMass-1) {
      mass = fitParaboloid->GetParameter(0);
      JES  = fitParaboloid->GetParameter(2);
      if (TMath::Sqrt(1/fitParaboloid->GetParameter(1)) < 2*TMath::Sqrt(mass)) {
        semiMajor = TMath::Sqrt(1/fitParaboloid->GetParameter(1));
        semiMinor = TMath::Sqrt(1/fitParaboloid->GetParameter(3));
        alpha     = fitParaboloid->GetParameter(4);

        massError = sqrt(pow(semiMajor * cos(alpha), 2) + pow(semiMinor * sin(alpha), 2));
        JESError  = sqrt(pow(semiMajor * sin(alpha), 2) + pow(semiMinor * cos(alpha), 2));
      }
      else {
        massError = -1;
        JESError  = -1;
      }
    }
    else {
      mass      = -1;
      massError = -1;
      JES       = -1;
      JESError  = -1;
    }

    // Fit again with previous result as range
    double sigmaLevel = 4;
    if (3*massError < resolMass) sigmaLevel = 8;
    fitParaboloid->SetRange(mass - sigmaLevel*massError, JES - sigmaLevel*JESError,
			                mass + sigmaLevel*massError, JES + sigmaLevel*JESError);

    sumLogLikelihood->Fit("fitParaboloid","EMR0");
  
    double contours[3] = {1, 4, 9};
    fitParaboloid->SetContour(3, contours);

    if (firstBinMass+1 < fitParaboloid->GetParameter(0) && fitParaboloid->GetParameter(0) < lastBinMass-1) {
      mass = fitParaboloid->GetParameter(0);
      JES  = fitParaboloid->GetParameter(2);
      if (TMath::Sqrt(1/fitParaboloid->GetParameter(1)) < 2*TMath::Sqrt(mass)) {
        semiMajor = TMath::Sqrt(1/fitParaboloid->GetParameter(1));
        semiMinor = TMath::Sqrt(1/fitParaboloid->GetParameter(3));
        alpha     = fitParaboloid->GetParameter(4);

        massError = sqrt(pow(semiMajor * cos(alpha), 2) + pow(semiMinor * sin(alpha), 2));
        JESError  = sqrt(pow(semiMajor * sin(alpha), 2) + pow(semiMinor * cos(alpha), 2));
      }
      else{
        massError = -1;
        JESError  = -1;
      }
    }
    else {
      mass      = -1;
      JES       = -1;
      massError = -1;
      JESError  = -1;
    }
    SetValue("mass_mTop_JES", mass, massError);
    SetValue("JES_mTop_JES" , JES , JESError );

    fitParaboloid->SetParameter(5, 0);

    // stat+syst ellipsis  
    double mSyst = 1.18;
    double jSyst = 0.012;
  
    double sm2 = massError*massError + mSyst*mSyst;
    double sj2 = JESError *JESError  + jSyst*jSyst;
  
    fitParaboloid->Copy(*systParaboloid);
    systParaboloid->SetRange(mass - 40, JES - 0.4, mass + 40, JES + 0.4);
    systParaboloid->SetParameter(1, 2*cos(2*alpha)/(sm2 - sj2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
    systParaboloid->SetParameter(3, 2*cos(2*alpha)/(sj2 - sm2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
    systParaboloid->SetLineColor(kBlack);
    systParaboloid->SetLineStyle(7);
    //systParaboloid->SetLineWidth(5);
  
    //* Set minL to 0
    sumLogLikelihood->Add(hUnity, -sumLogLikelihood->GetMinimum(0) + 1e-2);
    sumLogLikelihood->SetAxisRange(0, 25, "Z");
    if (blackWhite) sumLogLikelihood->Draw("AXIG");
    else sumLogLikelihood->Draw("COLZ");
    if (syst) systParaboloid->Draw("cont3 same");
    fitParaboloid->Draw("cont3 same");
    //*/
  
    // create legend
    TLegend *leg0 = new TLegend(0.2, 0.15, 0.45, 0.25);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    //leg0->AddEntry((TObject*)0, "1, 2, 3#sigma", "");
    leg0->AddEntry(fitParaboloid, "stat", "L");
    if (syst)    leg0->AddEntry(systParaboloid, "stat + syst", "L");
    leg0->Draw();
  
    helper->DrawCMS();
  
    std::cout << "massError: " << massError << std::endl;
  
    TString path(plotPath); path+= fIdentifier_; path += "_"; path += i; path += "_"; path += j; path += ".eps";
    ctemp->Print(path);
    delete leg0;
  }
  else{

    double leftMargin  = ctemp->GetLeftMargin();
    double rightMargin = ctemp->GetRightMargin();
    ctemp->SetLeftMargin (0.20);
    ctemp->SetRightMargin(0.02);

    //std::cout << "sumLogLikelihood entries: " << sumLogLikelihood->GetEntries() << std::endl;
    TH1D* sumLogLikelihood1D = sumLogLikelihood->ProjectionX("sumLogLikelihood1D", sumLogLikelihood->GetYaxis()->FindBin(1.), sumLogLikelihood->GetYaxis()->FindBin(1.));
    sumLogLikelihood1D->GetYaxis()->SetTitle("-2 #Delta ln(L)");
    sumLogLikelihood1D->GetYaxis()->SetTitleOffset(1.65);

    TF1* unity1D = new TF1("unity1D", "1 + 0*x");
    TH1D* hUnity1D = new TH1D("hUnity1D","hUnity1D", binsMass, firstBinMass, lastBinMass);
    hUnity1D->Eval(unity1D);
    sumLogLikelihood1D->Add(hUnity1D, -sumLogLikelihood1D->GetMinimum() + 1e-2);

    sumLogLikelihood1D->Draw("E");
  
    fitParabola->SetParameter(2, sumLogLikelihood1D->GetMinimum(0));
    fitParabola->SetParameter(1, 100);
    fitParabola->SetRange(minMass - 2, minMass + 2);
  
    //fitParabola->SetRange(sumLogLikelihood1D->GetBinCenter(sumLogLikelihood->GetMinimumBin()) - 3, sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) + 3);

    sumLogLikelihood1D->Fit("fitParabola","EMR");

    //sumLogLikelihood1D->GetXaxis()->SetRangeUser(sumLogLikelihood->GetXaxis()->GetBinLowEdge(1), sumLogLikelihood->GetXaxis()->GetBinLowEdge(sumLogLikelihood->GetNbinsX()+1));
    sumLogLikelihood1D->GetXaxis()->SetRangeUser(minMass - 3, minMass + 3);
 
    double massConstJES      = -1;
    double massConstJESError = -1;
    if (firstBinMass+1 < fitParabola->GetParameter(0) && fitParabola->GetParameter(0) < lastBinMass-1) {
      massConstJES      = fitParabola->GetParameter(0);
      massConstJESError = TMath::Sqrt(1/fitParabola->GetParameter(1));
    }
    else {
      massConstJES      = -1;
      massConstJESError = -1;
    }
    SetValue("mass_mTop", massConstJES, massConstJESError);

    std::cout << "Fixed JES: m_t = ";
    std::cout << massConstJES << " +/- " << massConstJESError << std::endl;

    helper->DrawCMS();
  
    TString path1D(plotPath); path1D+= fIdentifier_; path1D += "_"; path1D += i; path1D += "_"; path1D += j; path1D += "_1D.eps";
    ctemp->Print(path1D);
  
    //sumLogLikelihood1D->Delete();
    ctemp->SetLeftMargin(leftMargin);
    ctemp->SetRightMargin(rightMargin);

    delete sumLogLikelihood1D;
    delete unity1D;
    delete hUnity1D;
  }
  
  delete fitParabola;
  delete fitParaboloid;
  delete systParaboloid;
  delete null;
  delete unity;
  delete ctemp;
  delete eventLikelihood;
  delete logEventLikelihood;
  delete sumLogLikelihood;
  delete productLikelihood;

  delete hUnity;
  delete helper;

  delete topEvent;
  delete weightEvent;

  std::cout << "IdeogramAnalyzerNewInterface done" << std::endl;
}

IdeogramAnalyzerNewInterface::~IdeogramAnalyzerNewInterface()
{
  delete fptr_;
  delete combLikelihood_;
}
