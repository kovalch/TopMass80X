#include "IdeogramAnalyzerMinimizer.h"

#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "IdeogramSampleLikelihood.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <cmath>
#include <iomanip>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>

#include "TCanvas.h"
#include "TFrame.h"
//#include "TColor.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>
#include "TGraph.h"

typedef ProgramOptionsReader po;

IdeogramAnalyzerMinimizer::IdeogramAnalyzerMinimizer(const std::string& identifier, TTree* tree) :
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

void IdeogramAnalyzerMinimizer::Analyze(const std::string& cuts, int iBin, int jBin) {
  std::cout << "Starting IdeogramAnalyzer ..." << std::endl;
  time_t start, end;
  time(&start);
  time(&end);

  Scan(cuts, iBin, jBin);
  NumericalMinimization();

  time(&end);
  std::cout << "Finished IdeogramAnalyzer in " << difftime(end, start) << " seconds." << std::endl;
}


void IdeogramAnalyzerMinimizer::Scan(const std::string& cuts, int iBin, int jBin)
{
  int nEvents = sample_.nEvents;

  std::cout << "nEvents: " << nEvents << std::endl;

  std::string plotPath("plot/Ideogram/");
  {
    // Build Likelihood
    boost::progress_display progress((int)nEvents, std::cout);
    for (const auto& event : sample_.events) {
      ++progress;
      
      std::vector<IdeogramCombLikelihood*> permutationFunctions;
      
      //TODO - negative weights

      for (int iComb = 0, maxComb = event.permutations.size(); iComb < maxComb; ++iComb) {
        double topMass = event.permutations.at(iComb).topMass;
        double wMass   = event.permutations.at(iComb).wMass;
        double prob    = event.permutations.at(iComb).prob;
        int leptonFlavour = event.leptonFlavour;
        int bin        = event.permutations.at(iComb).bin;
        
        if (bin != iBin) continue;

        if (prob != 0) {
          if (channelID_ == Helper::kAllJets) {
            permutationFunctions.push_back(new IdeogramCombLikelihoodAllJets());
          }
          else {
            permutationFunctions.push_back(new IdeogramCombLikelihoodLeptonJets());
          }
          permutationFunctions.back()->SetFixedParams(prob, topMass, wMass, abs(leptonFlavour));
        }
      }
      eventFunctions_.push_back(permutationFunctions);
    } // end for
  }
}

IdeogramAnalyzerMinimizer::~IdeogramAnalyzerMinimizer()
{
  delete fptr_;
  delete combLikelihood_;
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

  ROOT::Math::Functor f(likelihood, 2); 
  double step[2] = {0.01,0.01};

  // starting point
  double variable[2] = { 172.5, 1.};

  min->SetFunction(f);

  // Set the free variables to be minimized!
  min->SetVariable(0, "mass", variable[0], step[0]);
  min->SetVariable(1, "jes" , variable[1], step[1]);

  // do the minimization
  min->Minimize(); 
  
  double mass      = min->X()[0];
  double massError = min->Errors()[0];
  double JES       = min->X()[1];
  double JESError  = min->Errors()[1];
  
  std::cout << "Minimum: f(" << mass << "," << JES << "): " 
           << min->MinValue()  << std::endl;

  // DRAW
  if (po::GetOption<bool>("minPlot")) {
    Helper* helper = new Helper();
    
    TCanvas* canv = new TCanvas("canv", "Top mass", 500, 500);
    canv->SetFrameFillColor(kRed);
    canv->SetFrameFillStyle(1001);
    /*
    canv->Update();
    canv->GetFrame()->SetFillColor(21);
    canv->GetFrame()->SetBorderSize(12);
    canv->GetFrame()->SetFillStyle(1001);
    canv->Modified();*/
    
    canv->SetLeftMargin (0.20);
    canv->SetRightMargin(0.04);
    canv->cd();
    
    unsigned int numPoints = 100;
    double* contourxs = new double[numPoints];
    double* contourys = new double[numPoints];

    unsigned int x = 0; unsigned int y = 1;
    
    int lineColor = kBlack;
    int lineWidth = 1;
    
    min->SetErrorDef(9.);
    min->Contour(x, y, numPoints, contourxs, contourys);
    TGraph* gr3 = new TGraph(numPoints, contourxs, contourys);
    gr3->SetFillColor(kSpring-9);
    gr3->SetLineColor(lineColor);
    gr3->SetLineWidth(lineWidth);
    
    min->SetErrorDef(4.);
    min->Contour(x, y, numPoints, contourxs, contourys);
    TGraph* gr2 = new TGraph(numPoints, contourxs, contourys);
    gr2->SetFillColor(kAzure+1);
    gr2->SetLineColor(lineColor);
    gr2->SetLineWidth(lineWidth);
    
    min->SetErrorDef(1.);
    min->Contour(x, y, numPoints, contourxs, contourys);
    TGraph* gr1 = new TGraph(numPoints, contourxs, contourys);
    gr1->SetFillColor(kViolet+9);
    gr1->SetLineColor(lineColor);
    gr1->SetLineWidth(lineWidth);
    
    gr3->SetTitle("-2 #Delta log(L); m_{t} [GeV]; JES");
    gr3->GetYaxis()->SetTitleOffset(1.7);
    
    gr3->Draw("ACF");
    gr3->Draw("C,SAME");
    gr2->Draw("CF,SAME");
    gr2->Draw("C,SAME");
    gr1->Draw("CF,SAME");
    gr1->Draw("C,SAME");
    
    TGraph* gr0 = new TGraph(1, &mass, &JES);
    gr0->SetMarkerColor(kWhite);
    gr0->SetMarkerStyle(2);
    gr0->SetMarkerSize(2);
    gr0->Draw("P,SAME");
    
    TLegend *leg0 = new TLegend(0.70, 0.75, 0.93, 0.92);
    leg0->SetFillStyle(1001);
    leg0->SetFillColor(kWhite);
    leg0->SetBorderSize(1);
    leg0->AddEntry(gr1, "1#sigma contour", "F");
    leg0->AddEntry(gr2, "2#sigma contour", "F");
    leg0->AddEntry(gr3, "3#sigma contour", "F");
    leg0->Draw();

    helper->DrawCMS();
    
    std::string path("plot/Ideogram/"); path+= fIdentifier_; path += std::string(".eps");
    canv->Print(path.c_str());
  }
  
  min->SetFixedVariable(1, "jes", 1.);
  min->Minimize();
  double mass1d      = min->X()[0];
  double mass1dError = min->Errors()[0];
  
  SetValue("mass_mTop_JES", mass, massError);
  SetValue("JES_mTop_JES" , JES , JESError );
  SetValue("mass_mTop", mass1d, mass1dError);
  
  std::cout << "m_t = " << mass << " +/- " << massError << std::endl;
  std::cout << "JES = " << JES  << " +/- " <<  JESError << std::endl;
  std::cout << "m_t(1d) = " << mass1d << " +/- " << mass1dError << std::endl;
}
