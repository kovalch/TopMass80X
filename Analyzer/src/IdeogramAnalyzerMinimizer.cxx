#include "IdeogramAnalyzerMinimizer.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "IdeogramSampleLikelihood.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/progress.hpp>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"

typedef ProgramOptionsReader po;

IdeogramAnalyzerMinimizer::IdeogramAnalyzerMinimizer(const std::string& identifier, TTree* tree) :
    MassAnalyzer(identifier, tree),
    sample_(*(new DataSample())),
    channelID_(Helper::channelID()),
    entries_(0),
    maxPermutations_              (po::GetOption<int   >("analysisConfig.maxPermutations")),
    isFastSim_                    (po::GetOption<int   >("fastsim")),
    shapeSystematic_              (po::GetOption<double>("shape"  )),
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
  NumericalMinimization();
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
      
      if (channelID_ == Helper::kAllJets) {
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
    } // end for
  }
  
  {
    // Set Likelihood parameters
    entries_   = 0;
    int iEvent = 0;
    for (const auto& event : sample_.events) {
      
      //TODO - negative weights

      for (int iComb = 0, maxComb = event.permutations.size(); iComb < maxComb; ++iComb) {
        double topMass = event.permutations.at(iComb).topMass;
        double wMass   = event.permutations.at(iComb).wMass;
        double prob    = event.permutations.at(iComb).prob;
        int leptonFlavour = event.leptonFlavour;
        int bin        = event.permutations.at(iComb).bin;
        
        if (bin != iBin) continue;

        ++entries_;

        if (prob != 0) {
          eventFunctions_[iEvent][iComb]->SetFixedParams(prob, topMass, wMass, abs(leptonFlavour), shapeSystematic_, permutationFractionSystematic_, isFastSim_);
          eventFunctions_[iEvent][iComb]->SetActive(true);
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
  if(channelID_ == Helper::kAllJets) toFit = {kMass, kJES, kFSig, kFCP};
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
  for(unsigned int i = 0; i < toFit.size(); ++i){
    // Set the free variables to be minimized!
    if     (toFit[i] == kMass) { SetValue("mass"+nameFreeVariables, min->X()[0], min->Errors()[0]); }
    else if(toFit[i] == kJES ) { SetValue("JES" +nameFreeVariables, min->X()[1], min->Errors()[1]); }
    else if(toFit[i] == kFSig) { SetValue("fSig"+nameFreeVariables, min->X()[2], min->Errors()[2]); }
    else if(toFit[i] == kFCP ) { SetValue("fCP" +nameFreeVariables, min->X()[3], min->Errors()[3]); }
  }
  SetValue("Entries", entries_, sqrt(entries_));

  // DRAW
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

  // do the next combination of variables
  for(unsigned int i = start; i < toFit.size(); ++i){
    std::vector<allowedVariables> copyFit = toFit;
    copyFit.erase(copyFit.begin()+i);
    IterateVariableCombinations(min, copyFit, i);
  }
  return;
}

void IdeogramAnalyzerMinimizer::PlotResult(ROOT::Math::Minimizer* min, IdeogramAnalyzerMinimizer::allowedVariables x, IdeogramAnalyzerMinimizer::allowedVariables y){
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
  gr3.SetNameTitle("","");
  gr3.SetFillColor(kSpring-9);
  gr3.SetLineColor(lineColor);
  gr3.SetLineWidth(lineWidth);

  min->SetErrorDef(4.);
  min->Contour(x, y, numPoints, contourxs, contourys);
  contourxs[numPoints] = contourxs[0]; contourys[numPoints] = contourys[0];
  TGraph gr2(numPoints+1, contourxs, contourys);
  gr2.SetFillColor(kAzure+1);
  gr2.SetLineColor(lineColor);
  gr2.SetLineWidth(lineWidth);

  min->SetErrorDef(1.);
  min->Contour(x, y, numPoints, contourxs, contourys);
  contourxs[numPoints] = contourxs[0]; contourys[numPoints] = contourys[0];
  TGraph gr1(numPoints+1, contourxs, contourys);
  gr1.SetFillColor(kViolet+9);
  gr1.SetLineColor(lineColor);
  gr1.SetLineWidth(lineWidth);

  std::string plotNamePostfix("_");
  if     (x == kMass) { gr3.GetXaxis()->SetTitle("m_{t} [GeV]"); plotNamePostfix += "mass_"; }
  else if(x == kJES ) { gr3.GetXaxis()->SetTitle("JES");         plotNamePostfix += "JES_" ; }
  else if(x == kFSig) { gr3.GetXaxis()->SetTitle("f_{sig}");     plotNamePostfix += "fSig_"; }
  else if(x == kFCP ) { gr3.GetXaxis()->SetTitle("f_{CP}");      plotNamePostfix += "fCP_" ; }
  plotNamePostfix += "vs_";
  if     (y == kMass) { gr3.GetYaxis()->SetTitle("m_{t} [GeV]"); plotNamePostfix += "mass"; }
  else if(y == kJES ) { gr3.GetYaxis()->SetTitle("JES");         plotNamePostfix += "JES" ; }
  else if(y == kFSig) { gr3.GetYaxis()->SetTitle("f_{sig}");     plotNamePostfix += "fSig"; }
  else if(y == kFCP ) { gr3.GetYaxis()->SetTitle("f_{CP}");      plotNamePostfix += "fCP" ; }
  gr3.GetYaxis()->SetTitleOffset(1.7);

  gr3.Draw("ACF");
  gr3.Draw("C,SAME");
  gr2.Draw("CF,SAME");
  gr2.Draw("C,SAME");
  gr1.Draw("CF,SAME");
  gr1.Draw("C,SAME");

  TGraph gr0(1, &min->X()[x], &min->X()[y]);
  gr0.SetMarkerColor(kWhite);
  gr0.SetMarkerStyle(2);
  gr0.SetMarkerSize(2);
  gr0.Draw("P,SAME");

  TLegend leg0(0.70, 0.75, 0.93, 0.92);
  leg0.SetFillStyle(1001);
  leg0.SetFillColor(kWhite);
  leg0.SetBorderSize(1);
  leg0.AddEntry(&gr1, "1#sigma contour", "F");
  leg0.AddEntry(&gr2, "2#sigma contour", "F");
  leg0.AddEntry(&gr3, "3#sigma contour", "F");
  leg0.Draw();

  helper.DrawCMS();

  std::string path("plot/Ideogram/"); path+=HelperFunctions::cleanedName(fIdentifier_)+plotNamePostfix;
  canv.Print((path+std::string(".eps" )).c_str(), "eps");
  canv.Print((path+std::string(".root")).c_str(), "root");

  delete[] contourxs;
  delete[] contourys;
}

// cleanup needed to run pseudo-experiments
void IdeogramAnalyzerMinimizer::CleanUp(){
  for (unsigned int i=0, l=eventFunctions_.size(); i<l; ++i){
    for (unsigned int j=0, m=eventFunctions_[i].size(); j<m; ++j){
      eventFunctions_[i][j]->SetActive(false);
    }
  }
}
