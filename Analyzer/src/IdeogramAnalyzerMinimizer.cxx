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
//#include <boost/lexical_cast.hpp>

#include "TCanvas.h"
//#include "TColor.h"
//#include "TError.h"
//#include "TF2.h"
//#include "TFrame.h"
#include "TGraph.h"
//#include "TH1D.h"
//#include "TH2D.h"
#include "TLegend.h"
//#include "TRandom2.h"
//#include "TStyle.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"

typedef ProgramOptionsReader po;

IdeogramAnalyzerMinimizer::IdeogramAnalyzerMinimizer(const std::string& identifier, TTree* tree) :
    MassAnalyzer(identifier, tree),
    sample_(*(new DataSample())),
    //fptr_(0),
    //combLikelihood_(0),
    channelID_(Helper::channelID()),
    isFastSim_                    (po::GetOption<int   >("fastsim"  )),
    shapeSystematic_              (po::GetOption<double>("shape"    )),
    permutationFractionSystematic_(po::GetOption<double>("permu"    ))
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
        permutationFunctions.push_back(new IdeogramCombLikelihoodAllJets());
        permutationFunctions.back()->SetActive(false);
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
  //delete fptr_;
  //delete combLikelihood_;
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

  // all possible result variables go here
  double mass3d = -10, mass3dError = -10;
  double  JES3d = -10,  JES3dError = -10;
  double fSig3d = -10, fSig3dError = -10;

  double mass2d1 = -10, mass2d1Error = -10;
  double  JES2d1 = -10,  JES2d1Error = -10;

  double mass2d2 = -10, mass2d2Error = -10;
  double fSig2d2 = -10, fSig2d2Error = -10;

  double  JES2d3 = -10,  JES2d3Error = -10;
  double fSig2d3 = -10, fSig2d3Error = -10;

  double mass1d = -10, mass1dError = -10;
  double  JES1d = -10,  JES1dError = -10;
  double fSig1d = -10, fSig1dError = -10;
  // create funciton wrapper for minmizer
  // a IMultiGenFunction type 
  IdeogramSampleLikelihood likelihood;
  likelihood.AddFunctions(eventFunctions_);

  ROOT::Math::Functor f(likelihood, likelihood.NDim());
  double step[] = {0.01,0.01,0.001};

  // starting point
  double variable[] = {172.5, 1., 0.70669};

  min->SetFunction(f);

  // Set the free variables to be minimized!
  min->SetVariable(0, "mass", variable[0], step[0]);
  min->SetVariable(1, "jes" , variable[1], step[1]);
  if (channelID_ == Helper::kAllJets) {
    min->SetLimitedVariable(2, "fSig", variable[2], step[2], 0.0, 1.0);
  }
  else {
    min->SetFixedVariable(2, "fSig", 1.0);
  }

  // do the minimization
  min->Minimize();

  if (channelID_ == Helper::kAllJets) {
    mass3d      = min->X()[0];
    mass3dError = min->Errors()[0];
     JES3d      = min->X()[1];
     JES3dError = min->Errors()[1];
    fSig3d      = min->X()[2];
    fSig3dError = min->Errors()[2];
    std::cout << "Minimum: f(" << mass3d << "," << JES3d << "," << fSig3d << "): " << min->MinValue()  << std::endl;
  }
  else{
    mass2d1      = min->X()[0];
    mass2d1Error = min->Errors()[0];
     JES2d1      = min->X()[1];
     JES2d1Error = min->Errors()[1];
    std::cout << "Minimum: f(" << mass2d1 << "," << JES2d1 << "): " << min->MinValue()  << std::endl;
  }

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
    
    double mass, JES;
    if (channelID_ == Helper::kAllJets) {
      mass = mass3d;
       JES =  JES3d;
    }
    else {
      mass = mass2d1;
       JES =  JES2d1;
    }
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
  if (channelID_ == Helper::kAllJets) {
    min->SetVariable(0, "mass", variable[0], step[0]);
    min->SetVariable(1, "jes" , variable[1], step[1]);
    min->SetFixedVariable(2, "fSig", -1.0);
    min->Minimize();

    mass2d1      = min->X()[0];
    mass2d1Error = min->Errors()[0];
     JES2d1      = min->X()[1];
     JES2d1Error = min->Errors()[1];

    min->SetVariable(0, "mass", variable[0], step[0]);
    min->SetFixedVariable(1, "jes", variable[1]);
    min->SetLimitedVariable(2, "fSig", variable[2], step[2], 0.0, 1.0);
    min->Minimize();

    mass2d2      = min->X()[0];
    mass2d2Error = min->Errors()[0];
    fSig2d2      = min->X()[2];
    fSig2d2Error = min->Errors()[2];

    min->SetFixedVariable(0, "mass", variable[0]);
    min->SetVariable(1, "jes", variable[1], step[1]);
    min->SetLimitedVariable(2, "fSig", variable[2], step[2], 0.0, 1.0);
    min->Minimize();

     JES2d3      = min->X()[1];
     JES2d3Error = min->Errors()[1];
    fSig2d3      = min->X()[2];
    fSig2d3Error = min->Errors()[2];

    min->SetFixedVariable(0, "mass", variable[0]);
    min->SetVariable(1, "jes" , variable[1], step[1]);
    min->SetFixedVariable(2, "fSig", -1.0);
    min->Minimize();

    JES1d      = min->X()[1];
    JES1dError = min->Errors()[1];

    min->SetFixedVariable(0, "mass", variable[0]);
    min->SetFixedVariable(1, "jes", variable[1]);
    min->SetLimitedVariable(2, "fSig", variable[2], step[2], 0.0, 1.0);
    min->Minimize();

    fSig1d      = min->X()[2];
    fSig1dError = min->Errors()[2];
  }

  min->SetVariable(0, "mass", variable[0], step[0]);
  min->SetFixedVariable(1, "jes", variable[1]);
  if (channelID_ == Helper::kAllJets) {
    min->SetFixedVariable(2, "fSig", -1.0);
  }
  else {
    min->SetFixedVariable(2, "fSig", 1.0);
  }
  min->Minimize();
  mass1d      = min->X()[0];
  mass1dError = min->Errors()[0];

  SetValue("mass_mTop_JES_fSig", mass3d, mass3dError);
  SetValue("JES_mTop_JES_fSig" ,  JES3d,  JES3dError);
  SetValue("fSig_mTop_JES_fSig", fSig3d, fSig3dError);

  SetValue("mass_mTop_JES", mass2d1, mass2d1Error);
  SetValue("JES_mTop_JES" ,  JES2d1,  JES2d1Error);

  SetValue("mass_mTop_fSig", mass2d2, mass2d2Error);
  SetValue("fSig_mTop_fSig", fSig2d2, fSig2d2Error);

  SetValue("JES_JES_fSig" ,  JES2d3,  JES2d3Error);
  SetValue("fSig_JES_fSig", fSig2d3, fSig2d3Error);
  SetValue("mass_mTop", mass1d, mass1dError);
  SetValue("JES_JES"  ,  JES1d,  JES1dError);
  SetValue("fSig_fSig", fSig1d, fSig1dError);
    
  delete min;
}

// cleanup needed to run pseudo-experiments
void IdeogramAnalyzerMinimizer::CleanUp(){
  for (unsigned int i=0, l=eventFunctions_.size(); i<l; ++i){
    for (unsigned int j=0, m=eventFunctions_[i].size(); j<m; ++j){
      eventFunctions_[i][j]->SetActive(false);
    }
  }
}
