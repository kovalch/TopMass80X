#include "IdeogramAnalyzerMinimizer.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "IdeogramSampleLikelihood.h"
#include "IdeogramEventLikelihood.h"
#include "Helper.h"
#include "CMS_lumi.h"
#include "ProgramOptionsReader.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/progress.hpp>

#include "TArrow.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF2.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"

typedef ProgramOptionsReader po;

//Tree out of Constructor is never used! needs DataSample over SetDataSample()
IdeogramAnalyzerMinimizer::IdeogramAnalyzerMinimizer(const std::string& identifier, TTree* tree) :
    MassAnalyzer(identifier, tree),
    sample_(*(new DataSample())),
    channelID_(Helper::channelID()),
    entries_(0),
    maxPermutations_              (po::GetOption<int   >("analysisConfig.maxPermutations")),
    drawIdeograms_                (po::GetOption<int   >("drawIdeograms")),
    isFastSim_                    (po::GetOption<int   >("fastsim")),
    shapeSystematic_              (po::GetOption<double>("shape"  )),
    shapeSystematic2_             (po::GetOption<double>("shape2" )),
    permutationFractionSystematic_(po::GetOption<double>("permu"  ))
    //topBranchName_                (po::GetOption<std::string>("topBranchName"))
{
}

void IdeogramAnalyzerMinimizer::Analyze(const std::string& cuts, int iBin, int jBin) {
  std::cout << "Starting IdeogramAnalyzer ..." << std::endl;
  time_t start, end;
  time(&start);
  time(&end);

  Scan(cuts, iBin, jBin);
  if (drawIdeograms_ > 0) { // just draw the ideograms
    DrawIdeograms(drawIdeograms_);
  }
  else { // extract results from sample
    NumericalMinimization();
  }
  CleanUp();

  time(&end);
  std::cout << "Finished IdeogramAnalyzer in " << difftime(end, start) << " seconds." << std::endl;
}


void IdeogramAnalyzerMinimizer::Scan(const std::string& cuts, int iBin, int jBin)
{
  int nEvents = sample_.nEvents;

  std::cout << "nEvents: " << nEvents << std::endl;

  // Create functions
  if (eventFunctions_.size() == 0) {
    if (po::GetOption<std::string>("task") == "pe") nEvents *= 1.5;
    boost::progress_display progress((int)nEvents, std::cout);
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
      ++progress;
      
      std::vector<IdeogramCombLikelihood*> permutationFunctions;
      
      if (channelID_ == Helper::kHamburg || channelID_ == Helper::kAllJets) {
        for (int iComb = 0; iComb < maxPermutations_; ++iComb) {
          permutationFunctions.push_back(new IdeogramCombLikelihoodAllJets());
          permutationFunctions.back()->SetActive(false);
        }
      }
      else {
        for (int iComb = 0; iComb < 4; ++iComb) {
          permutationFunctions.push_back(new IdeogramCombLikelihoodLeptonJets());
          permutationFunctions.back()->SetActive(false);
        }
      }

      eventFunctions_.push_back(permutationFunctions);

      if (channelID_ == Helper::kHamburg) {
        permutationFunctions.clear();
        for (int iComb = 0; iComb < 4; ++iComb) {
          permutationFunctions.push_back(new IdeogramCombLikelihoodLeptonJets());
          permutationFunctions.back()->SetActive(false);
        }
        eventFunctions_.push_back(permutationFunctions);
      }
    } // end for
  }
  
  {
    // Set Likelihood parameters
    entries_   = 0;
    int iEvent = 0, iEventLept = 0, iEventJets = 0;
    for (const auto& event : sample_.events) {

      //TODO - negative weights

      int leptonFlavour = event.leptonFlavour;

      int counter = iEvent;
      if(channelID_ == Helper::kHamburg) {
        if(std::abs(leptonFlavour) == 0){
          counter = 2*iEventJets;
          ++iEventJets;
        }
        else{
          counter = 2*iEventLept+1;
          ++iEventLept;
        }
      }

      for (int iComb = 0, maxComb = event.permutations.size(); iComb < maxComb; ++iComb) {
        double topMass = event.permutations.at(iComb).topMass;
        double wMass   = event.permutations.at(iComb).wMass;
        double prob    = event.permutations.at(iComb).prob;
        double weight  = event.weight/fabs(event.weight);
        int bin        = event.permutations.at(iComb).bin;
        
        if (bin != iBin) continue;

        entries_ = entries_ + (int) weight;

        if (prob != 0) {
          eventFunctions_[counter][iComb]->SetFixedParams(prob, topMass, wMass, abs(leptonFlavour), shapeSystematic_, shapeSystematic2_, permutationFractionSystematic_, isFastSim_, weight);
          eventFunctions_[counter][iComb]->SetActive(true);
        }
      }
      ++iEvent;
    } // end for
  }
}

IdeogramAnalyzerMinimizer::~IdeogramAnalyzerMinimizer()
{
}

void IdeogramAnalyzerMinimizer::NumericalMinimization() {
  // create minimizer giving a name and a name (optionally) for the specific
  // algorithm
  // possible choices are: 
  //     minName                  algoName
  // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
  //  Minuit2                     Fumili2
  //  Fumili
  //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
  //                              BFGS2, SteepestDescent
  //  GSLMultiFit
  //   GSLSimAn
  //   Genetic
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type 
  IdeogramSampleLikelihood likelihood;
  likelihood.AddFunctions(eventFunctions_);

  ROOT::Math::Functor f(likelihood, likelihood.NDim());

  min->SetFunction(f);

  std::vector<allowedVariables> toFit;
  if(channelID_ == Helper::kAllJets || channelID_ == Helper::kHamburg) toFit = {kMass, kJES, kFSig, kFCP};
  else toFit = {kMass, kJES};
  IterateVariableCombinations(min, toFit);

  delete min;
}

void IdeogramAnalyzerMinimizer::IterateVariableCombinations(ROOT::Math::Minimizer* min, std::vector<IdeogramAnalyzerMinimizer::allowedVariables> toFit, unsigned int start)
{
  if(!toFit.size()) return;
  // starting point
  static double variable[] = {172.5, 1., po::GetOption<double>("templates.fSig"), po::GetOption<double>("templates.fCP")};
  double step[] = {0.01,0.01,0.001,0.001};

  min->SetFixedVariable(0, "mass", variable[0]);
  min->SetFixedVariable(1, "jes" , variable[1]);
  min->SetFixedVariable(2, "fSig", variable[2]);
  min->SetFixedVariable(3, "fCP" , variable[3]);
  min->SetFixedVariable(4, "jsfc", 0.);

  std::string nameFreeVariables;
  for(unsigned int i = 0; i < toFit.size(); ++i){
    // Set the free variables to be minimized!
    if(toFit[i] == kMass) {
      min->SetVariable(0, "mass", variable[0], step[0]);
      nameFreeVariables += "_mTop";
    }
    else if(toFit[i] == kJES ) {
      min->SetVariable(1, "jes" , variable[1], step[1]);
      nameFreeVariables += "_JES";
    }
    else if(toFit[i] == kFSig) {
      min->SetLimitedVariable(2, "fSig", variable[2], step[2], 0.0, 1.0);
      nameFreeVariables += "_fSig";
    }
    else if(toFit[i] == kFCP ) {
      min->SetLimitedVariable(3, "fCP" , variable[3], step[3], 0.0, 1.0);
      nameFreeVariables += "_fCP";
   }
  }

  // do the minimization
  min->Minimize();
  //print cov matrix
  /*
  std::cout << "covariance matrix:\n";
  for(unsigned int i = 0; i < 4 ; ++i) {
    for(unsigned int j = 0 ; j < 4; ++j) {
      std::cout << min->CovMatrix(i,j) << " ";
    }
    std::cout << '\n';
  }
  */
  for(unsigned int i = 0; i < toFit.size(); ++i){
    // Set the free variables to be minimized!
    if     (toFit[i] == kMass) { SetValue("mass"+nameFreeVariables, min->X()[0], min->Errors()[0]); }
    else if(toFit[i] == kJES ) { SetValue("JES" +nameFreeVariables, min->X()[1], min->Errors()[1]); }
    else if(toFit[i] == kFSig) { SetValue("fSig"+nameFreeVariables, min->X()[2], min->Errors()[2]); }
    else if(toFit[i] == kFCP ) { SetValue("fCP" +nameFreeVariables, min->X()[3], min->Errors()[3]); }
  }
  SetValue("Entries", entries_, sqrt(entries_));
  /*
  for(short i=0; i<4; ++i){
    for(short j=0; j<4; ++j){
      std::cout << min->Correlation(i,j) << " ";
    }
    std::cout << std::endl;
  }
  */
  // DRAW
  if (po::GetOption<bool>("temPlot")) {
    if(channelID_ == Helper::kAllJets){
      if(nameFreeVariables == "_mTop_JES_fSig"){
	PlotResult2(min);
      }
    }
    else{
      if(nameFreeVariables == "_mTop_JES"){
        PlotResult2(min);
      }
    }
  }
  if (po::GetOption<bool>("minPlot")) {
    if(channelID_ == Helper::kAllJets){
      if(nameFreeVariables == "_mTop_JES_fSig_fCP"){
        PlotResult(min, kMass, kJES );
        PlotResult(min, kMass, kFSig);
        PlotResult(min, kMass, kFCP );
        PlotResult(min, kJES , kFSig);
        PlotResult(min, kJES , kFCP );
        PlotResult(min, kFSig, kFCP );
      }
    }
    else{
      if(nameFreeVariables == "_mTop_JES"){
        PlotResult(min);
      }
    }
  }
  
  // do the minimization with JSF constraint
  if (po::GetOption<bool>("constrainJSF")) {
    double s2 = min->Errors()[1] * min->Errors()[1];
    if (po::GetOption<double>("constrainJSFsigma") > -1.) {
      min->SetFixedVariable(4, "jsfc", po::GetOption<double>("constrainJSFsigma"));
    }
    else { // calculate constraint from weight
      min->SetFixedVariable(4, "jsfc", sqrt(s2/po::GetOption<double>("constrainJSFweight") - s2));
    }
    min->Minimize();
    for(unsigned int i = 0; i < toFit.size(); ++i){
      // Set the free variables to be minimized!
      if     (toFit[i] == kMass) { SetValue("mass"+nameFreeVariables+"_jsfc", min->X()[0], min->Errors()[0]); }
      else if(toFit[i] == kJES ) { SetValue("JES" +nameFreeVariables+"_jsfc", min->X()[1], min->Errors()[1]); }
      else if(toFit[i] == kFSig) { SetValue("fSig"+nameFreeVariables+"_jsfc", min->X()[2], min->Errors()[2]); }
      else if(toFit[i] == kFCP ) { SetValue("fCP" +nameFreeVariables+"_jsfc", min->X()[3], min->Errors()[3]); }
    }
    SetValue("Entries", entries_, sqrt(entries_));
  }
  
  if (po::GetOption<bool>("minPlot")) {
    if(nameFreeVariables == "_mTop_JES"){
      PlotResult(min, kMass, kJES, true);
    }
  }

  // do the next combination of variables
  for(unsigned int i = start; i < toFit.size(); ++i){
    std::vector<allowedVariables> copyFit = toFit;
    copyFit.erase(copyFit.begin()+i);
    IterateVariableCombinations(min, copyFit, i);
  }
  return;
}

void IdeogramAnalyzerMinimizer::PlotResult(ROOT::Math::Minimizer* min, IdeogramAnalyzerMinimizer::allowedVariables x, IdeogramAnalyzerMinimizer::allowedVariables y, bool hybrid){
  Helper helper;

  TCanvas canv("canv", "Top mass", 500, 500);
  canv.SetFrameFillColor(kRed);
  canv.SetFrameFillStyle(1001);
  /*
  canv.Update();
  canv.GetFrame()->SetFillColor(21);
  canv.GetFrame()->SetBorderSize(12);
  canv.GetFrame()->SetFillStyle(1001);
  canv.Modified();
  */
  canv.SetLeftMargin (0.20);
  canv.SetRightMargin(0.04);
  canv.cd();

  unsigned int numPoints = 100;
  double* contourxs = new double[numPoints+1];
  double* contourys = new double[numPoints+1];

  int lineColor = kBlack;
  int lineWidth = 1;

  min->SetErrorDef(9.);
  min->Contour(x, y, numPoints, contourxs, contourys);
  contourxs[numPoints] = contourxs[0]; contourys[numPoints] = contourys[0];
  TGraph gr3(numPoints+1, contourxs, contourys);
  gr3.SetNameTitle("gr3","gr3");
  gr3.SetFillColor(kSpring-9);
  gr3.SetLineColor(lineColor);
  gr3.SetLineWidth(lineWidth);

  min->SetErrorDef(4.);
  min->Contour(x, y, numPoints, contourxs, contourys);
  contourxs[numPoints] = contourxs[0]; contourys[numPoints] = contourys[0];
  TGraph gr2(numPoints+1, contourxs, contourys);
  gr2.SetNameTitle("gr2","gr2");	
  gr2.SetFillColor(kAzure+1);
  gr2.SetLineColor(lineColor);
  gr2.SetLineWidth(lineWidth);

  min->SetErrorDef(1.);
  min->Contour(x, y, numPoints, contourxs, contourys);
  contourxs[numPoints] = contourxs[0]; contourys[numPoints] = contourys[0];
  TGraph gr1(numPoints+1, contourxs, contourys);
  gr1.SetNameTitle("gr1","gr1");
  gr1.SetFillColor(kViolet+9);
  gr1.SetLineColor(lineColor);
  gr1.SetLineWidth(lineWidth);
  
  min->SetErrorDef(1.);

  std::string plotNamePostfix("_");
  if     (x == kMass) { gr3.GetXaxis()->SetTitle("m_{t} [GeV]"); plotNamePostfix += "mass_"; }
  else if(x == kJES ) { gr3.GetXaxis()->SetTitle("JSF");         plotNamePostfix += "JES_" ; }
  else if(x == kFSig) { gr3.GetXaxis()->SetTitle("f_{sig}");     plotNamePostfix += "fSig_"; }
  else if(x == kFCP ) { gr3.GetXaxis()->SetTitle("f_{CP}");      plotNamePostfix += "fCP_" ; }
  plotNamePostfix += "vs_";
  if     (y == kMass) { gr3.GetYaxis()->SetTitle("m_{t} [GeV]"); plotNamePostfix += "mass"; }
  else if(y == kJES ) { gr3.GetYaxis()->SetTitle("JSF");         plotNamePostfix += "JES" ; }
  else if(y == kFSig) { gr3.GetYaxis()->SetTitle("f_{sig}");     plotNamePostfix += "fSig"; }
  else if(y == kFCP ) { gr3.GetYaxis()->SetTitle("f_{CP}");      plotNamePostfix += "fCP" ; }
  if (hybrid) plotNamePostfix += "_hyb";
  gr3.GetYaxis()->SetTitleOffset(1.7);

  gr3.Draw("ACF");
  gr3.Draw("C,SAME");
  gr2.Draw("CF,SAME");
  gr2.Draw("C,SAME");
  gr1.Draw("CF,SAME");
  gr1.Draw("C,SAME");

  TGraph gr0(1, &min->X()[x], &min->X()[y]);
  gr0.SetNameTitle("gr0","gr0");
  gr0.SetMarkerColor(kWhite);
  gr0.SetMarkerStyle(2);
  gr0.SetMarkerSize(2);
  gr0.Draw("P,SAME");

  TLegend leg0(0.70, 0.75, 0.93, 0.92);
  leg0.SetFillStyle(1001);
  leg0.SetFillColor(kWhite);
  leg0.SetBorderSize(1);
  leg0.AddEntry(&gr1, "-2#Delta log(L) = 2.30", "F");
  leg0.AddEntry(&gr2, "-2#Delta log(L) = 6.17", "F");
  leg0.AddEntry(&gr3, "-2#Delta log(L) = 11.8", "F");
  leg0.Draw();

  // draw tangents to -2#Delta log(L) = 1 ellipse
  if (channelID_ == Helper::kAllJets) {
    double minX = min->X()[x]-min->Errors()[x];
    double maxX = min->X()[x]+min->Errors()[x];
    double minY = min->X()[y]-min->Correlation(x,y)*min->Errors()[y];
    double maxY = min->X()[y]+min->Correlation(x,y)*min->Errors()[y];
    double yAtXaxis = gr3.GetYaxis()->GetXmin();
    TArrow line1(minX,minY,minX,yAtXaxis,0.04,"-------|>");
    TArrow line2(maxX,maxY,maxX,yAtXaxis,0.04,"-------|>");
    line1.SetAngle(45);
    line2.SetAngle(45);
    line1.SetLineColor(1);
    line2.SetLineColor(1);

    line1.Draw();
    line2.Draw();
  }
  
  //helper.DrawCMS();

  std::string path("plot/Ideogram/"); path+=HelperFunctions::cleanedName(fIdentifier_)+plotNamePostfix;
  canv.Print((path+std::string(".eps" )).c_str(), "eps");
  canv.Print((path+std::string(".root")).c_str(), "root");

  delete[] contourxs;
  delete[] contourys;
}

double IdeogramAnalyzerMinimizer::evalAllJets(double *x, double *p)
{
  IdeogramCombLikelihoodAllJets like;
  like.SetFixedParams(1, x[0], x[1], 0, 0, 0, 0, 0, 1);
  like.SetActive(true);
  double pNew[] = {p[0],p[1],p[2],p[3]};
  double val = like.Evaluate(pNew,0);
  return val;
}

void IdeogramAnalyzerMinimizer::PlotResult2(ROOT::Math::Minimizer* min, IdeogramAnalyzerMinimizer::allowedVariables x, IdeogramAnalyzerMinimizer::allowedVariables y){
  TFile* histFile = TFile::Open("tmpTESThistos.root","RECREATE");
  TH1D* histT = new TH1D("histT","",po::GetOption<double>("templates.maxTopMass")-100,100,po::GetOption<double>("templates.maxTopMass"));
  TH1D* histW = new TH1D("histW","",110,65,120);
  TH2D* hist = new TH2D("hist","",po::GetOption<double>("templates.maxTopMass")-100,100,po::GetOption<double>("templates.maxTopMass"),110,65,120);
  
  TF2* func = new TF2("func",evalAllJets,100,po::GetOption<double>("templates.maxTopMass"),65,120,4);
  func->SetParameters(min->X()[x],min->X()[y], min->X()[kFSig], po::GetOption<double>("templates.fCP"));
  for(int i=1; i<histT->GetNbinsX(); ++i){
    double x = histT->GetBinCenter(i);
    for(int j=1; j<histW->GetNbinsX(); ++j){
      double y = histW->GetBinCenter(j);
      double val = func->Eval(x,y);
      histT->Fill(x,val);
      histW->Fill(y,val);
      hist->Fill(x,y,val);
    } // end for j
  } // end for i
  histFile->cd();
  histT->Write();
  histW->Write();
  hist->Write();
  func->Write();
  histFile->Close();
}

void IdeogramAnalyzerMinimizer::DrawIdeograms(int n) {
  Helper* helper = new Helper();
  helper->SetTDRStyle();
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetHatchesLineWidth(1.5);
  
  for (int i = 0; i < n; ++i) {
    if (!eventFunctions_[i][0]->IsActive()) continue;
    IdeogramEventLikelihood* ideogram = new IdeogramEventLikelihood();
    ideogram->SetFunction(eventFunctions_[i]);
    
    std::cout << i << std::endl;
    for (const auto& permutation : eventFunctions_[i]) {
      if (permutation->IsActive()) {
        std::cout << "prob " << permutation->GetFixedParam(0);
        std::cout << ", mt " << permutation->GetFixedParam(1);
        std::cout << ", mw " << permutation->GetFixedParam(2) << std::endl;
      }
    }
    
    TCanvas* canv = new TCanvas("canv", "ideogram", 500, 500);
    
    TF2* func = new TF2("func", ideogram, &IdeogramEventLikelihood::DoEval, 100, 250, 0.75, 1.25, 9, "IdeogramEventLikelihood", "DoEval");
    func->SetTitle(";m_{t} [GeV];JSF");
    func->SetNpx(200);
    func->SetNpy(200);
    func->Draw("cont1z");
    
    helper->DrawCMS(-1, -1, canv);
    
    canv->Print((std::string("plot/eventLikelihood/")+std::to_string(i)+std::string(".eps")).c_str());
  }
}

// cleanup needed to run pseudo-experiments
void IdeogramAnalyzerMinimizer::CleanUp(){
  for (unsigned int i=0, l=eventFunctions_.size(); i<l; ++i){
    for (unsigned int j=0, m=eventFunctions_[i].size(); j<m; ++j){
      eventFunctions_[i][j]->SetActive(false);
    }
  }
}
