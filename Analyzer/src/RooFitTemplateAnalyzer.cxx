#include "RooFitTemplateAnalyzer.h"

#include "Helper.h"
#include "ProgramOptionsReader.h"

#include <boost/lexical_cast.hpp>

#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"

#include "RooAddPdf.h"
#include "RooBifurGauss.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooLandau.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"

typedef ProgramOptionsReader po;

void RooFitTemplateAnalyzer::Analyze(const std::string& cuts, int i, int j) {
  Scan(cuts, i, j, "mTop_JES_fSig");
  Scan(cuts, i, j, "mTop_fSig");
  Scan(cuts, i, j, "mTop_JES");
  Scan(cuts, i, j, "mTop");
}

void RooFitTemplateAnalyzer::Scan(const std::string& cuts, int i, int j, TString variables)
{
  const int nTemplTypes = 2;
  const int nComboTypes = 3;

  RooRealVar comboTypeVar   = RooRealVar("comboType"     ,"comboType"  ,  0.,   6.,"");
  RooRealVar prob           = RooRealVar("prob"          ,"P(#chi^{2})",  0.,   1.,"");
  RooRealVar MTOP           = RooRealVar("topMass"       ,"m_{t}^{fit}",100., 550.,"GeV");
  RooRealVar meanMW         = RooRealVar("meanWMass"     ,"m_{W}^{rec}",  0.,7000.,"GeV");
  RooRealVar combinedWeight = RooRealVar("combinedWeight","weight"     ,  0., 100.,"");
  
  MTOP.setRange("mTopFitRange",100.,550.);
  
  RooAddPdf* topAdd = 0;
  RooAddPdf*   wAdd = 0;
  
  std::string name = "";

  // set permutation fractions
  workspace->var("fCP")->setVal(2.31859251856803894e-01);
  workspace->var("fWP")->setVal(4.23728721216320992e-03+9.87463630735874176e-04+2.18985587358474731e-01);
  workspace->var("fUN")->setVal(5.43930411338806152e-01);


  // no calibration
  workspace->var("MassOffset"      )->setVal(0.0);
  workspace->var("MassSlopeMass"   )->setVal(0.0);
  workspace->var("MassSlopeJES"    )->setVal(0.0);
  workspace->var("MassSlopeMassJES")->setVal(0.0);

  workspace->var("JESOffset"       )->setVal(0.0);
  workspace->var("JESSlopeMass"    )->setVal(0.0);
  workspace->var("JESSlopeJES"     )->setVal(0.0);
  workspace->var("JESSlopeMassJES" )->setVal(0.0);

  /*
  // test calibration
  workspace->var("MassOffset"      )->setVal( 4.77254e-01);
  workspace->var("MassSlopeMass"   )->setVal( 2.54275e-02);
  workspace->var("MassSlopeJES"    )->setVal(-4.13860e-01);
  workspace->var("MassSlopeMassJES")->setVal(-4.93806e-01);

  workspace->var("JESOffset"       )->setVal(-3.55015e-03);
  workspace->var("JESSlopeMass"    )->setVal(-7.87469e-05);
  workspace->var("JESSlopeJES"     )->setVal(-5.19422e-02);
  workspace->var("JESSlopeMassJES" )->setVal( 4.41256e-03);
  */
  /*
  // parametrization of mean and width independent
  workspace->var("MassOffset"      )->setVal( 4.29624e-01+1.01432e-01);
  workspace->var("MassSlopeMass"   )->setVal( 8.00953e-02);
  workspace->var("MassSlopeJES"    )->setVal(-8.77033e+00);
  workspace->var("MassSlopeMassJES")->setVal(-4.08061e-01);
  
  workspace->var("JESOffset"       )->setVal(-9.38332e-03-1.21815e-03);
  workspace->var("JESSlopeMass"    )->setVal(-2.67342e-04);
  workspace->var("JESSlopeJES"     )->setVal( 1.07114e-01);
  workspace->var("JESSlopeMassJES" )->setVal( 2.72959e-03);

  // parametrization of width as const*mean
  workspace->var("MassOffset"      )->setVal( 8.74877e-01+1.48996e-01);
  workspace->var("MassSlopeMass"   )->setVal( 1.14649e-02);
  workspace->var("MassSlopeJES"    )->setVal(-1.18615e+01);
  workspace->var("MassSlopeMassJES")->setVal(-3.44788e-01);

  workspace->var("JESOffset"       )->setVal(-1.14916e-02-1.45748e-03);
  workspace->var("JESSlopeMass"    )->setVal( 2.42972e-05);
  workspace->var("JESSlopeJES"     )->setVal( 1.19464e-01);
  workspace->var("JESSlopeMassJES" )->setVal( 3.07011e-03);
  */
  /*
  // parametrization of mean and constant width
  workspace->var("MassOffset"      )->setVal( 8.43067e-01);
  workspace->var("MassSlopeMass"   )->setVal( 1.11208e-02);
  workspace->var("MassSlopeJES"    )->setVal(-8.57688e+00);
  workspace->var("MassSlopeMassJES")->setVal(-3.33429e-01);

  workspace->var("JESOffset"       )->setVal(-1.06016e-02);
  workspace->var("JESSlopeMass"    )->setVal( 3.05550e-05);
  workspace->var("JESSlopeJES"     )->setVal( 8.41141e-02);
  workspace->var("JESSlopeMassJES" )->setVal( 3.29863e-03);
  */

  topAdd = (RooAddPdf*)workspace->pdf("mTopPDF");
  wAdd   = (RooAddPdf*)workspace->pdf("mWPDF"  );
  
  RooArgSet varSet = RooArgSet(prob,MTOP,meanMW,"varSet");
  
  RooDataSet data = RooDataSet("data","data",varSet,RooFit::Import(*fTree_),RooFit::WeightVar("prob"));
  
  name = "";
  for(int templType = 0; templType < nTemplTypes; ++templType){
    for(int comboType = 0; comboType < nComboTypes; ++comboType){
      int h = nComboTypes*templType + comboType;
      const unsigned nAlpha = (templType == 0 && comboType == 0) ? 2 : ((templType == 1) ? 3 : 4);
      const unsigned nPar   = 4*nAlpha;
      ////const unsigned nAlpha = (templType == 0 && comboType == 0) ? 2 : 3;
      //const unsigned nAlpha = (templType == 0 && comboType == 0) ? 2 : ((templType == 1) ? 3 : 4);
      //const unsigned n4ParVars = (templType == 0 && comboType != 0) ? 2 : 1;
      //const unsigned nPar   = 4*n4ParVars + nAlpha-n4ParVars;
      for(unsigned p=0; p<nPar; ++p) {
        name = "par_"; name += h; name += "_"; name += p;
        workspace->var(name.c_str())->setConstant(true);
        //workspace->var(name)->removeError();
        //workspace->var(name)->removeAsymError();
      }
    }
  }

  RooRealVar topBKGGammaNorm  = RooRealVar("topBKGGammaNorm" , "topBKGGammaNorm" ,  0.852647/*0.848667*/);
  RooRealVar topBKGGammaGamma = RooRealVar("topBKGGammaGamma", "topBKGGammaGamma",  3.88501/*4.14752*/  );
  RooRealVar topBKGGammaBeta  = RooRealVar("topBKGGammaBeta" , "topBKGGammaBeta" , 33.7309/*33.7793*/   );
  RooRealVar topBKGGammaMu    = RooRealVar("topBKGGammaMu"   , "topBKGGammaMu"   , 93.2179/*89.6645*/   );
  RooGamma   topBKGGamma      = RooGamma  ("topBKGGamma"     , "topBKGGamma"     , MTOP, topBKGGammaGamma, topBKGGammaBeta, topBKGGammaMu);
  topBKGGamma.setNormRange("mTopFitRange");
  
  RooRealVar topBKGLandauMean  = RooRealVar("topBKGLandauMean" , "topBKGLandauMean" , 206.862/*217.114*/ );
  RooRealVar topBKGLandauSigma = RooRealVar("topBKGLandauSigma", "topBKGLandauSigma",  25.7198/*27.6076*/);
  RooLandau  topBKGLandau      = RooLandau ("topBKGLandau"     , "topBKGLandau"     , MTOP, topBKGLandauMean, topBKGLandauSigma);
  topBKGLandau.setNormRange("mTopFitRange");
  
  RooAddPdf topBKG = RooAddPdf("topBKG", "topBKG", topBKGGamma, topBKGLandau, topBKGGammaNorm);
  topBKG.setNormRange("mTopFitRange");
  
  RooRealVar wBKGMean       = RooRealVar("wBKGMean"      , "wBKGMean"      , 91.9298/*86.9853*/ );
  RooRealVar wBKGSigmaLeft  = RooRealVar("wBKGSigmaLeft" , "wBKGSigmaLeft" ,  6.23531/*5.78569*/);
  RooRealVar wBKGSigmaRight = RooRealVar("wBKGSigmaRight", "wBKGSigmaRight",  8.12332/*7.12755*/);
  RooBifurGauss wBKG        = RooBifurGauss("wBKG", "wBKG", meanMW, wBKGMean, wBKGSigmaLeft, wBKGSigmaRight);
  wBKG.setNormRange("mTopFitRange");

  RooRealVar fSig = RooRealVar("fSig", "fSig", 0.539, 0, 1);
  //fSig.setVal(1);
  
  RooAddPdf topModel = RooAddPdf("topModel", "topModel", *topAdd, topBKG, fSig);
  RooAddPdf wModel   = RooAddPdf("wModel"  , "wModel"  , *wAdd  , wBKG  , fSig);
  topModel.setNormRange("mTopFitRange");
  wModel  .setNormRange("mTopFitRange");
  
  RooProdPdf model = RooProdPdf("model", "model", topModel, wModel);
  model.setNormRange("mTopFitRange");
  
  workspace->var("JES" )->setVal(1.0);
  workspace->var("mTop")->setVal(172.5);

  //workspace->var("JES" )->removeError();
  //workspace->var("JES" )->removeAsymError();
  //workspace->var("mTop")->removeError();
  //workspace->var("mTop")->removeAsymError();

  workspace->var("JES" )->setError(0.0);
  workspace->var("mTop")->setError(0.0);
  workspace->var("JES" )->setAsymError(0.0,0.0);
  workspace->var("mTop")->setAsymError(0.0,0.0);

  if( variables.Contains("JES") ) workspace->var("JES" )->setConstant(false);
  else                            workspace->var("JES" )->setConstant(true);
  if( variables.Contains("mTop")) workspace->var("mTop")->setConstant(false);
  else                            workspace->var("mTop")->setConstant(true);
  if(!variables.Contains("fSig")) fSig.setConstant(true);
  else                            fSig.setConstant(false);
  
  RooAbsReal* nll = model.createNLL(data,RooFit::CloneData(false));
  RooMinuit mini(*nll);
  mini.setStrategy(2);
  mini.migrad();
  mini.migrad();
  //mini.improve();
  mini.hesse();
  
  RooFitResult* fitResult = mini.save();
  int fitStatus = fitResult->status();
  delete fitResult;

  std::cout << "FIT STATUS: " << fitStatus << std::endl;

  if(variables.Contains("fSig") && variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      SetValue("fSig_mTop_JES_fSig", fSig.getVal(), fSig.getError());
      SetValue("mass_mTop_JES_fSig", workspace->var("mTop")->getVal(), workspace->var("mTop")->getError());
      SetValue("JES_mTop_JES_fSig" , workspace->var("JES" )->getVal(), workspace->var("JES" )->getError());
    }
    else{
      SetValue("fSig_mTop_JES_fSig", -1., -1.);
      SetValue("mass_mTop_JES_fSig", -1., -1.);
      SetValue("JES_mTop_JES_fSig" , -1., -1.);
    }
  }
  if(variables.Contains("fSig") && !variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      SetValue("mass_mTop_fSig", workspace->var("mTop")->getVal(), workspace->var("mTop")->getError());
      SetValue("fSig_mTop_fSig", fSig.getVal(), fSig.getError());
    }
    else{
      SetValue("mass_mTop_fSig", -1., -1.);
      SetValue("fSig_mTop_fSig", -1., -1.);
    }
  }
  if(!variables.Contains("fSig") && variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      SetValue("mass_mTop_JES", workspace->var("mTop")->getVal(), workspace->var("mTop")->getError());
      SetValue("JES_mTop_JES" , workspace->var("JES" )->getVal(), workspace->var("JES" )->getError());
    }
    else{
      SetValue("mass_mTop_JES", -1., -1.);
      SetValue("JES_mTop_JES" , -1., -1.);
    }
  }
  if(!variables.Contains("fSig") && !variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      SetValue("mass_mTop", workspace->var("mTop")->getVal(), workspace->var("mTop")->getError());
    }
    else{
      SetValue("mass_mTop", -1., -1.);
    }
  }
  
  if(!po::GetOption<std::string>("task").compare("sm")){
    TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 1, 1, 600, 600);
    canvas1->cd();
    RooPlot* frame1 = MTOP.frame(RooFit::Range(100,350));
    data  .statOn (frame1, RooFit::Layout(.5, .9, .9));
    model .paramOn(frame1, RooFit::Layout(.5, .9, .7));
    data   .plotOn(frame1);
    model  .plotOn(frame1, RooFit::LineColor(8));
    topAdd->plotOn(frame1, RooFit::LineColor(kRed) , RooFit::Normalization(   fSig.getVal()));
    topBKG .plotOn(frame1, RooFit::LineColor(kBlue), RooFit::Normalization(1.-fSig.getVal()));
    //workspace->pdf("sig_0")   ->plotOn(frame1, RooFit::LineColor(kRed+1), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fCP")->getVal()*/));
    //workspace->pdf("sig_1")   ->plotOn(frame1, RooFit::LineColor(kRed+2), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fWP")->getVal()*/));
    //workspace->pdf("sig_2")   ->plotOn(frame1, RooFit::LineColor(kRed+3), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fUN")->getVal()*/));
    frame1->Draw();
    std::string path("plot/RooFit/mTop_"); path+= fIdentifier_; path += std::string("_"); path += boost::lexical_cast<std::string>(i); path += std::string("_"); path += boost::lexical_cast<std::string>(j); path += std::string("_"); path += variables; path += std::string(".eps");
    canvas1->Print(path.c_str());
  
    /*
    //TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 641, 1, 600, 600);
    //canvas2->cd();
    //
    //RooPlot* frame2 = fSig.frame();
    ////nll.plotOn(frame2,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_fSig = nll->createProfile(fSig);
    //pll_fSig->plotOn(frame2,RooFit::LineColor(kRed));
    //frame2->Draw();
    //
    //TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 641, 1, 600, 600);
    //canvas3->cd();
    //
    //RooPlot* frame3 = workspace->var("mTop")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
    ////nll.plotOn(frame3,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_mTop = nll->createProfile(*workspace->var("mTop"));
    //pll_mTop->plotOn(frame3,RooFit::LineColor(kRed));
    //frame3->Draw();
    //
    //TCanvas* canvas4 = new TCanvas("canvas4", "canvas4", 641, 1, 600, 600);
    //canvas4->cd();
    //
    //RooPlot* frame4 = workspace->var("JES")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
    ////nll.plotOn(frame4,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_JES = nll->createProfile(*workspace->var("JES"));
    //pll_JES->plotOn(frame4,RooFit::LineColor(kRed));
    //frame4->Draw();
    */

    TCanvas* canvas5 = new TCanvas("canvas5", "canvas5", 641, 1, 600, 600);
    canvas5->cd();
    RooPlot* frame5 = meanMW.frame(RooFit::Range(60,140));
    data .statOn (frame5, RooFit::Layout(.5, .9, .9));
    model.paramOn(frame5, RooFit::Layout(.5, .9, .7));
    data  .plotOn(frame5);
    model .plotOn(frame5, RooFit::LineColor(8));
    wAdd ->plotOn(frame5, RooFit::LineColor(kRed) , RooFit::Normalization(   fSig.getVal()));
    wBKG  .plotOn(frame5, RooFit::LineColor(kBlue), RooFit::Normalization(1.-fSig.getVal()));
    //workspace->pdf("sig_3")   ->plotOn(frame5, RooFit::LineColor(kRed+1), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fCP")->getVal()*/));
    //workspace->pdf("sig_4")   ->plotOn(frame5, RooFit::LineColor(kRed+2), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fWP")->getVal()*/));
    //workspace->pdf("sig_5")   ->plotOn(frame5, RooFit::LineColor(kRed+3), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fUN")->getVal()*/));
    frame5->Draw();
    path = "plot/RooFit/mW_"; path+= fIdentifier_; path += std::string("_"); path += boost::lexical_cast<std::string>(i); path += std::string("_"); path += boost::lexical_cast<std::string>(j); path += std::string("_"); path += variables; path += std::string(".eps");
    canvas5->Print(path.c_str());
  
    delete canvas1;
    delete canvas5;
    delete frame1;
    delete frame5;
  }

  delete nll;
  //delete workspace;

  //tmpFile->Close();
  //delete tmpFile;

  std::cout << "RooFitTemplateAnalyzer done" << std::endl;
}

RooWorkspace* RooFitTemplateAnalyzer::workspace(0);

RooFitTemplateAnalyzer::RooFitTemplateAnalyzer(const std::string& identifier, TTree* tree) : MassAnalyzer(identifier, tree)
{
  if(!workspace){
    std::string workspaceFilePath = std::string(gSystem->pwd()) + std::string("/RooFitCalibration.root");
    TFile * workspaceFile = TFile::Open(workspaceFilePath.c_str());
    gROOT->cd();
    workspace = (RooWorkspace*)((RooWorkspace*)workspaceFile->Get("workspaceMtop"))->Clone();
    workspaceFile->Close();
    delete workspaceFile;
  }
}

RooFitTemplateAnalyzer::~RooFitTemplateAnalyzer()
{
}
