#include "TCanvas.h"
//#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
//#include "TH2F.h"
//#include "TTree.h"
//#include "TMath.h"

#include <iostream>
#include <string>

void deriveShapeUncertainty()
{
  //TFile* fBkgVal7TeV = new TFile("/afs/desy.de/user/e/eschliec/xxl/topMass_7TeV_MixOverToy.root");
  TFile* fBkgVal7TeV = new TFile("/nfs/dust/cms/user/eschliec/TopMass/19/topMass_toy_data_ttbar.root");

  TCanvas* cBkgVal7TeV = (TCanvas*)fBkgVal7TeV->Get("c1_n2");

  TH1F* hTopToy    = (TH1F*)cBkgVal7TeV->GetPrimitive("hTopToy"   );
  TH1F* hTopToyMix = (TH1F*)cBkgVal7TeV->GetPrimitive("hTopToyMix");

  double nRebin = 1;
  hTopToy->GetYaxis()->SetRangeUser(0,0.05*nRebin);

  hTopToy   ->Rebin(nRebin);
  hTopToyMix->Rebin(nRebin);

  TCanvas* cBkgVal7TeVFits = new TCanvas("cBkgVal7TeVFits","",600,600);
  cBkgVal7TeVFits->cd();

  TF1* myGamma = new TF1("myGamma","[0]*TMath::GammaDist(x,[1],[2],[3])", 100, 550);
  myGamma->FixParameter(0,8.30073*nRebin);
  //myGamma->SetParameter(0,8.30073);
  //myGamma->SetParLimits(0,0,100);
  myGamma->SetParameter(1,3.5);
  myGamma->SetParLimits(1,0,10);
  myGamma->FixParameter(2,100);
  //myGamma->SetParameter(2,100);
  //myGamma->SetParLimits(2,80,120);
  myGamma->FixParameter(3,42.3623);
  //myGamma->SetParameter(3,70);
  //myGamma->SetParLimits(3,0,200);

  myGamma->SetLineColor(kGreen+3);
  hTopToy   ->Fit(myGamma,"I L M Q","0", 100, 550);
  hTopToy->Draw();
  myGamma->DrawCopy("same");
  double gammaToy = myGamma->GetParameter(1);

  myGamma->SetLineColor(4);
  hTopToyMix->Fit(myGamma,"I L M Q","0", 100, 550);
  hTopToyMix->Draw("same");
  myGamma->DrawCopy("same");
  double gammaToyMix = myGamma->GetParameter(1);

  double gammaMixOverToy = gammaToyMix/gammaToy;

  std::cout << gammaToyMix << " / " << gammaToy << " = " << gammaMixOverToy << std::endl;

  double gammaDef = 4.15584;
  double betaDef  = 42.3623;
  double muDef    = 100;

  TF1* myGammaRatio = new TF1("myGammaRatio", "[0]*TMath::GammaDist(x,[1],[2],[3])/([4]*TMath::GammaDist(x,[5],[6],[7]))", 100, 550);

  TCanvas* cShapeVariations = new TCanvas("cShapeVariations","",600,600);
  cShapeVariations->cd();

  // shape variation in 7 TeV closure test
  myGammaRatio->SetParameters(1,gammaMixOverToy*gammaDef,muDef,betaDef,1,gammaDef,muDef,betaDef);
  myGammaRatio->SetLineColor(kGreen+3);
  myGammaRatio->DrawCopy();

  // shape variation from statistical uncertainty
  myGammaRatio->SetParameters(1,1.148*gammaDef,muDef,betaDef,1,gammaDef,muDef,betaDef);
  myGammaRatio->SetLineColor(kRed);
  myGammaRatio->DrawCopy("same");
  myGammaRatio->SetParameters(1,0.852*gammaDef,muDef,betaDef,1,gammaDef,muDef,betaDef);
  myGammaRatio->DrawCopy("same");

  // shape variation from MC-loose vs. data-loose
  myGammaRatio->SetParameters(1,1.174*gammaDef,muDef,betaDef,1,gammaDef,muDef,betaDef);
  myGammaRatio->SetLineColor(kBlue);
  myGammaRatio->DrawCopy("same");
}
