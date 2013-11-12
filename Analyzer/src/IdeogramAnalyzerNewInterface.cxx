#include "IdeogramAnalyzerNewInterface.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <cmath>
#include <iomanip>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>

#include "TCanvas.h"
//#include "TColor.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"

typedef ProgramOptionsReader po;

IdeogramAnalyzerNewInterface::IdeogramAnalyzerNewInterface(const std::string& identifier, TTree* tree) :
    MassAnalyzer(identifier, tree),
    sample_(*(new DataSample())),
    fptr_(0),
    combLikelihood_(0),
    channelID_(Helper::channelID()),
    pullWidth_                    (po::GetOption<double>("pullWidth")),
    isFastSim_                    (po::GetOption<int   >("fastsim"  )),
    shapeSystematic_              (po::GetOption<double>("shape"    )),
    permutationFractionSystematic_(po::GetOption<double>("permu"    )),
    topBranchName_                (po::GetOption<std::string>("topBranchName"))
{
}

void IdeogramAnalyzerNewInterface::Analyze(const std::string& cuts, int i, int j) {
  std::cout << "Starting IdeogramAnalyzer ..." << std::endl;
  time_t start, end;
  time(&start);
  time(&end);

  Scan(cuts, i, j, 154, 190, 2, 0.9, 1.1, 0.02);

  double mass = GetValue("mass_mTop_JES").first;
  double JES  = GetValue("JES_mTop_JES" ).first;
  Scan(cuts, i, j, mass-2 , mass+2 , 0.1, JES-0.015, JES+0.015, 0.00075);

  double epsilon = 1e-6;
  mass = GetValue("mass_mTop_JES").first;
  //JES  = GetValue("JES_mTop_JES" ).first;
  Scan(cuts, i, j, mass-10, mass+10, 0.1 , 1.-(0.5*epsilon), 1.+(0.5*epsilon), epsilon, false);

  time(&end);
  std::cout << "Finished IdeogramAnalyzer in " << difftime(end, start) << " seconds." << std::endl;
}


void IdeogramAnalyzerNewInterface::Scan(const std::string& cuts, int i, int j, double firstBinMass, double lastBinMass,
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

  bool debug = false;
  int nDebug = 1;
  int minDebug = 0;
  int maxDebug = 10;

  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int binsMass = int((lastBinMass-firstBinMass)/resolMass);
  int binsJes  = int((lastBinJes -firstBinJes )/resolJes );

  if(!combLikelihood_ || !fptr_){
    if(channelID_ == Helper::kAllJets){
      IdeogramCombLikelihoodAllJets *fptrAllJets = new IdeogramCombLikelihoodAllJets();
      fptr_ = fptrAllJets;
      combLikelihood_ = new TF2("combLikelihood",fptrAllJets,&IdeogramCombLikelihoodAllJets::Evaluate, firstBinMass, lastBinMass, firstBinJes, lastBinJes, 7, "IdeogramCombLikelihoodAllJets", "Evaluate");
    }
    else {
      IdeogramCombLikelihoodLeptonJets *fptrLeptonJets = new IdeogramCombLikelihoodLeptonJets();
      fptr_ = fptrLeptonJets;
      combLikelihood_ = new TF2("combLikelihood",fptrLeptonJets,&IdeogramCombLikelihoodLeptonJets::Evaluate, firstBinMass, lastBinMass, firstBinJes, lastBinJes, 7, "IdeogramCombLikelihoodLeptonJets", "Evaluate");
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
  TF2* systParaboloid = 0;
  if (syst) systParaboloid = new TF2("systParaboloid", "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2 + [5]");
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
  //TH2D* productLikelihood  = new TH2D("productLikelihood" , "productLikelihood" , binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
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
  //productLikelihood->SetTitle("-2#upointln{L(m_{t}|sample)}");
  //productLikelihood->SetXTitle("m_{t} [GeV]");
  //productLikelihood->SetYTitle("JES");

  double fitWeight;
  int nEvents = sample_.nEvents;
  double sumWeights = 0.;

  std::cout << "nEvents: " << nEvents << std::endl;

  std::string plotPath("plot/Ideogram/");
  {
    // Build Likelihood
    //for (int iEntry = 0, length = sample_->nEvents; iEntry < length; ++iEntry) {
    boost::progress_display progress((int)nEvents, std::cout);
    int iEntry = -1;
    for (const auto& event : sample_.events) {
      ++iEntry;
      ++progress;

      eventLikelihood->Eval(null);
      eventLikelihood->SetFillColor(0);
      //weight = 0;
      fitWeight = 0;
      //currentWeight = 0;

      if (debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) {
        std::cout << std::setiosflags(std::ios::left)
        << std::setw(04) << "n"
        << std::setw(10) << "mt"
        << std::setw(10) << "mW"
        << std::setw(12) << "prob"
        << std::endl;
      }

      for (int iComb = 0, maxComb = event.permutations.size(); iComb < maxComb; ++iComb) {
        double topMass = event.permutations.at(iComb).topMass;
        double wMass   = event.permutations.at(iComb).wMass;
        double prob    = event.permutations.at(iComb).prob;
        int leptonFlavour = event.leptonFlavour;

        fitWeight += prob; //*CombinedWeight;

        if (debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) {
          std::cout << std::setw(04) << iComb
              << std::setw(10) << topMass
              << std::setw(10) << wMass
              << std::setw(12) << prob
              << std::endl;
        }

        if (prob != 0) {
          // Set Likelihood parameters
          // TODO electron channel
          //if(channelID_ == Helper::kAllJets){
          combLikelihood_->SetParameters(prob, topMass, wMass, abs(leptonFlavour), shapeSystematic_, permutationFractionSystematic_, isFastSim_);
          // add permutation to event likelihood
          eventLikelihood->Eval(combLikelihood_, "A");
        }
      }

      eventLikelihood->Scale(1./fitWeight);
      sumWeights += fitWeight;
      //if (weight == 0) continue;

      logEventLikelihood->Eval(null);

      for (int i = 0; i<=binsMass; i++) {
        for (int j = 0; j<=binsJes; j++) {
          logEventLikelihood->SetBinContent(i, j, -2*log(eventLikelihood->GetBinContent(i, j)));
        }
      }

      sumLogLikelihood->Add(logEventLikelihood, fitWeight/(pullWidth_*pullWidth_)); // add weight here

      if (debug && iEntry%nDebug == 0 && iEntry > minDebug && iEntry < maxDebug) {
        TCanvas* eventCanvas = new TCanvas("eventCanvas", "eventCanvas", 1200, 400);
        eventCanvas->Divide(3, 1);

        eventCanvas->cd(1);
        eventLikelihood->SetEntries(1);
        eventLikelihood->Draw("COLZ");
        eventLikelihood->SetFillColor(kRed+1);
        /*
        if (weight > 1./16) eventLikelihood->SetFillColor(kGreen);
        if (weight > 1./8 ) eventLikelihood->SetFillColor(kYellow);
        if (weight > 1./4 ) eventLikelihood->SetFillColor(kOrange);
        if (weight > 1./2 ) eventLikelihood->SetFillColor(kRed);
        */

        eventCanvas->cd(2);
        eventLikelihood->SetEntries(1);
        logEventLikelihood->Draw("COLZ");

        eventCanvas->cd(3);
        eventLikelihood->SetEntries(1);
        sumLogLikelihood->Draw("COLZ");

        std::string eventPath(plotPath); eventPath += HelperFunctions::cleanedName(fIdentifier_); eventPath += std::string("_"); eventPath += boost::lexical_cast<std::string>(iEntry); eventPath += std::string(".eps");
        //std::cout << eventPath << std::endl; //Print does this all by itself ...
        eventCanvas->Print(eventPath.c_str());

        delete eventCanvas;
      } // end if debug
    } // end for
  }
  delete eventLikelihood;
  delete logEventLikelihood;

  ctemp->cd();

  std::cout << "Fitting ..." << std::endl;

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
  double minLike = sumLogLikelihood->GetMinimum(0);

  std::cout << "minMass: " << minMass << ", minJes: " << minJes << std::endl;

  /*
  sumLogLikelihood->SetAxisRange(minMass - 3, minMass + 3, "X");
  sumLogLikelihood->SetAxisRange(minJes - 0.03, minJes + 0.03, "Y");
  sumLogLikelihood->SetAxisRange(0, 20, "Z");
  //*/

  //sumLogLikelihood->SetMarkerStyle(20);
  //sumLogLikelihood->SetMarkerColor(kRed+1);

  std::cout << "Minimum likelihood: " << minLike << "\tMaximum likelihood (in range): " << sumLogLikelihood->GetMaximum() << std::endl;

  Helper* helper = new Helper();

  if (fit2D) {
    fitParaboloid->SetParLimits(0, minMass-2*resolMass, minMass+2*resolMass);
    fitParaboloid->SetParameter(0, minMass);
    fitParaboloid->SetParameter(1, 100);
    fitParaboloid->SetParLimits(2, minJes-2*resolJes, minJes+2*resolJes);
    fitParaboloid->SetParameter(2, minJes);
    fitParaboloid->SetParameter(3, 1000000);
    fitParaboloid->SetParLimits(5, 0.95*minLike, 1.05*minLike);
    fitParaboloid->SetParameter(5, minLike);

    //fitParaboloid->SetRange(minMass - 1, minJes - 0.01, minMass + 1, minJes + 0.01);
    fitParaboloid->SetRange(minMass - 4, minJes - 0.04, minMass + 4, minJes + 0.04);
    //fitParaboloid->SetRange(minMass - 20, minJes - 0.2, minMass + 20, minJes + 0.2);

    //std::cout << "sumLogLikelihood: " << sumLogLikelihood << std::endl;
    sumLogLikelihood->Fit("fitParaboloid","EMR0");

    double semiMajor, semiMinor, alpha;

    double mass = fitParaboloid->GetParameter(0);
    double JES  = -1;
    double massError = -1;
    double JESError  = -1;
    if (firstBinMass+1 < mass && mass < lastBinMass-1) {
      JES  = fitParaboloid->GetParameter(2);
      if (sqrt(1./fitParaboloid->GetParameter(1)) < 2*sqrt(mass)) {
        semiMajor = sqrt(1./fitParaboloid->GetParameter(1));
        semiMinor = sqrt(1./fitParaboloid->GetParameter(3));
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

    mass = fitParaboloid->GetParameter(0);
    if (firstBinMass+1 < mass && mass < lastBinMass-1) {
      JES  = fitParaboloid->GetParameter(2);
      if (sqrt(1./fitParaboloid->GetParameter(1)) < 2*sqrt(mass)) {
        semiMajor = sqrt(1./fitParaboloid->GetParameter(1));
        semiMinor = sqrt(1./fitParaboloid->GetParameter(3));
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

    std::cout << "m_t = " << mass << " +/- " << massError << std::endl;
    std::cout << "JES = " << JES  << " +/- " <<  JESError << std::endl;

    fitParaboloid->SetParameter(5, 0);

    // stat+syst ellipsis
    double mSyst = 1.18;
    double jSyst = 0.012;

    double sm2 = massError*massError + mSyst*mSyst;
    double sj2 = JESError *JESError  + jSyst*jSyst;

    if (syst) {
      fitParaboloid->Copy(*systParaboloid);
      systParaboloid->SetRange(mass - 40, JES - 0.4, mass + 40, JES + 0.4);
      systParaboloid->SetParameter(1, 2*cos(2*alpha)/(sm2 - sj2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
      systParaboloid->SetParameter(3, 2*cos(2*alpha)/(sj2 - sm2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
      systParaboloid->SetLineColor(kBlack);
      systParaboloid->SetLineStyle(7);
      //systParaboloid->SetLineWidth(5);
    }

    //* Set minL to 0
    sumLogLikelihood->Add(hUnity, -minLike + 1e-2);
    sumLogLikelihood->SetAxisRange(0, 25, "Z");
    if (blackWhite) sumLogLikelihood->Draw("AXIG");
    else sumLogLikelihood->Draw("COLZ");
    if (syst) systParaboloid->Draw("cont3 same");
    fitParaboloid->Draw("cont3 same");
    //*/

    // create legend
    //TLegend *leg0 = new TLegend(0.2, 0.15, 0.45, 0.25);
    TLegend *leg0 = new TLegend(0.185, 0.135, 0.435, 0.235);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    //leg0->AddEntry((TObject*)0, "1, 2, 3#sigma", "");
    //leg0->AddEntry(fitParaboloid, "stat", "L");
    leg0->AddEntry(fitParaboloid, "#splitline{1#sigma, 2#sigma, 3#sigma}{contours}", "L");
    if (syst) leg0->AddEntry(systParaboloid, "stat + syst", "L");
    leg0->Draw();

    helper->DrawCMS();

    std::string path(plotPath); path+= HelperFunctions::cleanedName(fIdentifier_); path += "_"; path += boost::lexical_cast<std::string>(i); path += std::string("_"); path += boost::lexical_cast<std::string>(j); path +=  std::string("_"); path += HelperFunctions::cleanedName(topBranchName_); path += std::string(".eps");
    ctemp->Print(path.c_str());
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

    fitParabola->SetParameter(2, 0);
    fitParabola->SetParameter(1, 100);
    fitParabola->SetRange(minMass - 2, minMass + 2);

    //fitParabola->SetRange(sumLogLikelihood1D->GetBinCenter(sumLogLikelihood->GetMinimumBin()) - 3, sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) + 3);

    sumLogLikelihood1D->Fit("fitParabola","EMR");

    //sumLogLikelihood1D->GetXaxis()->SetRangeUser(sumLogLikelihood->GetXaxis()->GetBinLowEdge(1), sumLogLikelihood->GetXaxis()->GetBinLowEdge(sumLogLikelihood->GetNbinsX()+1));
    sumLogLikelihood1D->GetXaxis()->SetRangeUser(minMass - 3, minMass + 3);

    double massConstJES      = fitParabola->GetParameter(0);
    double massConstJESError = -1;
    if (firstBinMass+1 < massConstJES && massConstJES < lastBinMass-1) {
      massConstJESError = sqrt(1./fitParabola->GetParameter(1));
    }
    else {
      massConstJES      = -1;
      massConstJESError = -1;
    }
    SetValue("mass_mTop", massConstJES, massConstJESError);

    std::cout << "Fixed JES: m_t = " << massConstJES << " +/- " << massConstJESError << std::endl;

    helper->DrawCMS();

    std::string path1D(plotPath); path1D+= HelperFunctions::cleanedName(fIdentifier_); path1D += std::string("_"); path1D += boost::lexical_cast<std::string>(i); path1D += std::string("_"); path1D += boost::lexical_cast<std::string>(j); path1D += std::string("_"); path1D += HelperFunctions::cleanedName(topBranchName_); path1D += std::string("_1D.eps");
    ctemp->Print(path1D.c_str());

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
  delete sumLogLikelihood;
  //delete productLikelihood;

  delete hUnity;
  delete helper;

  std::cout << "Fitting done" << std::endl;
}

IdeogramAnalyzerNewInterface::~IdeogramAnalyzerNewInterface()
{
  delete fptr_;
  delete combLikelihood_;
}
