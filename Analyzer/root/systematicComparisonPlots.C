//#include "TArrow.h"
//#include "TAxis.h"
#include "TCanvas.h"
//#include "TF1.h"
//#include "TF2.h"
#include "TFile.h"
//#include "TGraphErrors.h"
//#include "TGraph2DErrors.h"
//#include "TH1F.h"
//#include "TLatex.h"
#include "TLegend.h"
//#include "TPaveText.h"
//#include "TMath.h"
//#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
//#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSeqCollection.h"
//#include "TGraphAsymmErrors.h"
#include "TPaveLabel.h"

#include "iostream"

TTree* modifiedTree(TTree * tree);
void savePlotsToDisk(TString format, bool rootFile);
void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color = kBlack);
void DrawCMSPrel();
void DrawCMSSim();

TH2F * bTagEff = 0;
TH2F * cTagEff = 0;
TH2F * lTagEff = 0;

//TString samplePath("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/");
TString samplePath("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/");

void systematicComparisonPlots()
{
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TFile * bTagFile = TFile::Open(samplePath+TString("bTagFile.root"));
  gROOT->cd();
  bTagEff = (TH2F*)bTagFile->Get("histb")->Clone();
  cTagEff = (TH2F*)bTagFile->Get("histc")->Clone();
  lTagEff = (TH2F*)bTagFile->Get("histl")->Clone();
  bTagFile->Close();

  TString plotString[] = {"Z2_F11_172_5"   , "Z2_F11_BJES_Up"      , "Z2_F11_BJES_Down",
                          "Z2_F11_172_5"   , "Z2_F11_BTAG_Up"      , "Z2_F11_BTAG_Down",
                          "Z2_F11_172_5"   , "Z2_F11_JES_Up"       , "Z2_F11_JES_Down",
  			  "Z2_F11_172_5"   , "Z2_F11_JER_Up"       , "Z2_F11_JER_Down",
                          "Z2_F11_172_5"   , "Z2_F11_Scale_Up"     , "Z2_F11_Scale_Down",
                          "Z2_F11_172_5"   , "Z2_F11_Matching_Up"  , "EMPTY",
                          "Z2_F11_FAST"    , "EMPTY"               , "Z2_F11_Matching_Down_FAST",
                          "Z2_F11_P11_FAST", "Z2_F11_P11mpiHi_FAST", "Z2_F11_P11TeV_FAST",
                          "Z2_F11_P11"     , "Z2_F11_P11_NoCR"     , "EMPTY"
                          };
  //TString plotString[] = {"Z2_F11_P11_FAST", "Z2_F11_P11mpiHi_FAST", "Z2_F11_P11TeV_FAST"};

  const int plots = sizeof(plotString)/sizeof(TString)/3;

  TFile * tmpFile = TFile::Open("tmpFileSystematicComparison.root","RECREATE");
  tmpFile->cd();

  // add histograms for Pt leading
  TH1F * hMT = new TH1F("hMT",";m_{t}^{fit} (GeV);arb. units",60,100,250);
  TH1F * hMW = new TH1F("hMW",";m_{W}^{reco} (GeV);arb. units",60,65,125);

  hMT->GetYaxis()->SetTitleOffset(1.6);
  hMW->GetYaxis()->SetTitleOffset(1.6);
  hMT->GetXaxis()->SetTitleOffset(1.1);
  hMW->GetXaxis()->SetTitleOffset(1.1);

  for(int i = 0; i < plots; ++i){
    TString canvName = "cMT_"; canvName += plotString[3*i+(i==6 ? 2 : 1)];
    TCanvas * canvMT = new TCanvas(canvName,canvName,1,1,600,600);
    canvMT->SetMargin(0.12,0.08,0.1,0.1);
    canvName = "cMW_"; canvName += plotString[3*i+(i==6 ? 2 : 1)];
    TCanvas * canvMW = new TCanvas(canvName,canvName,641,1,600,600);
    canvMW->SetMargin(0.12,0.08,0.1,0.1);

    TLegend * leg = new TLegend(0.55,0.65,0.94,0.89);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    for(int j = 0; j < 3; ++j){
      std::cout << "Creating " << i << " " << j << " " << plotString[3*i+j] << std::endl;
      if(plotString[3*i+j] == "EMPTY") continue;
      tmpFile->cd();
      TTree * tree = 0;
      if(tmpFile->Get(TString("tree_")+plotString[3*i+j])) tree = (TTree*)tmpFile->Get(TString("tree_")+plotString[3*i+j]);
      else{
	TString fileName = plotString[3*i+j]+TString("_sig.root");
	TFile * file = TFile::Open(samplePath+TString(fileName));
	TTree * oriTree = (TTree*)file->Get("FullHadTreeWriter/tree");
  
	tmpFile->cd();
	tree = modifiedTree(oriTree);
	file->Close();
      }

      double prob,topMass,meanWMass,combinedWeight,dRbb,ttMass;
      tree->SetBranchAddress("prob", &prob);
      tree->SetBranchAddress("topMass", &topMass);
      tree->SetBranchAddress("meanWMass", &meanWMass);
      tree->SetBranchAddress("combinedWeight", &combinedWeight);
      tree->SetBranchAddress("dRbb", &dRbb);
      tree->SetBranchAddress("ttMass", &ttMass);
      TClonesArray * jets = new TClonesArray("TLorentzVector");
      tree->SetBranchStatus("jets", 1);
      tree->GetBranch("jets")->SetAutoDelete(kFALSE);
      tree->SetBranchAddress("jets", &jets);
      int Njet;
      tree->SetBranchStatus ("Njet", 1);
      tree->SetBranchAddress("Njet", &Njet);
      float bTag[50];
      tree->SetBranchStatus ("bTag_CSV", 1);
      tree->SetBranchAddress("bTag_CSV", &bTag );
      unsigned short comboType;
      tree->SetBranchStatus("comboType", 1);
      tree->SetBranchAddress("comboType", &comboType);
      
      TString histName = TString("hMT_")+plotString[3*i+j];
      hMT->Clone(histName);
      TString exp = TString("topMass>>")+histName;
      canvMT->cd();
      if(j==0) tree->Draw(exp,"combinedWeight","norm hist");
      else     tree->Draw(exp,"combinedWeight","norm hist same");
      TH1F* hist = (TH1F*)gDirectory->Get(histName);
      hist->SetLineWidth(3);
      hist->GetYaxis()->SetRangeUser(0.,0.1);
      std::cout << "mTop (mean/sigma/skewness/kurtosis): " << hist->GetMean() << "/" << hist->GetRMS() << "/" << hist->GetSkewness() << "/" << hist->GetKurtosis() << std::endl;
      if     (j==0) hist->SetLineColor(8);     
      else if(j==1) hist->SetLineColor(kRed);  
      else if(j==2) hist->SetLineColor(kBlue); 
      leg->AddEntry(hist, plotString[3*i+j], "l");
     
      histName = TString("hMW_")+plotString[3*i+j];
      hMW->Clone(histName);
      exp = TString("meanWMass>>")+histName;
      canvMW->cd();
      if(j==0) tree->Draw(exp,"combinedWeight","norm hist");
      else     tree->Draw(exp,"combinedWeight","norm hist same");
      hist = (TH1F*)gDirectory->Get(histName);
      hist->SetLineWidth(3);
      hist->GetYaxis()->SetRangeUser(0.,0.1);
      std::cout << "mW (mean/sigma/skewness/kurtosis): " << hist->GetMean() << "/" << hist->GetRMS() << "/" << hist->GetSkewness() << "/" << hist->GetKurtosis() << std::endl;
      if     (j==0) hist->SetLineColor(8);
      else if(j==1) hist->SetLineColor(kRed);
      else if(j==2) hist->SetLineColor(kBlue);
    }
    canvMT->cd();
    leg->Draw("same");
    DrawCMSSim();
    canvMW->cd();
    leg->Draw("same");
    DrawCMSSim();
  }

  savePlotsToDisk("eps", true);

  //tmpFile->Close();
}

// return the PU weights for the different samples
enum enumForPUWeights {kSummer11, kSummer11Plus05, kSummer11Minus05, kFall11, kFall11Plus05, kFall11Minus05, kFall10};
double calcPUWeight_(enum enumForPUWeights sample, short nPU)
{
  // ----
  // default weight to be used for fall11 samples with S6 PU scenario
  // ----
  // 3.54fb-1 !!! precale weighted !!! 73.5mb
  //double weightsPUFall11[] = {0.0953411, 0.0209421, 0.0389577, 0.355278, 1.3007, 1.8981, 1.95124, 1.81828, 1.63994, 1.55631, 1.54116, 1.44885, 1.24065, 0.915466, 0.580092, 0.350426, 0.209707, 0.113744, 0.0509568, 0.0183767, 0.0054363, 0.00141399, 0.000360394, 0.000106747, 3.96548e-05, 1.65259e-05, 6.57705e-06, 2.32168e-06, 7.15271e-07, 1.99901e-07, 5.67217e-08, 1.86344e-08, 7.44984e-09, 3.25689e-09, 1.39639e-09, 5.58487e-10, 2.00145e-10, 6.60545e-11, 1.89381e-11, 4.82982e-12, 1.12634e-12, 2.21706e-13, 4.06965e-14, 6.1115e-15, 9.20171e-16, 1.05937e-16, 1.38254e-17, 1.14592e-18, 1.13769e-19, 43889.5};
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUFall11[] = {0.0956436, 0.0219907, 0.0713245, 0.735605, 1.85831, 2.19452, 2.04027, 1.76898, 1.58702, 1.53704, 1.45119, 1.24648, 0.904423, 0.552217, 0.319224, 0.179453, 0.0864062, 0.0326249, 0.0095542, 0.00229068, 0.000501167, 0.000121599, 3.80191e-05, 1.3964e-05, 4.89459e-06, 1.48493e-06, 3.82906e-07, 8.94646e-08, 2.18024e-08, 6.46154e-09, 2.31321e-09, 8.61013e-10, 3.0218e-10, 9.44128e-11, 2.57874e-11, 6.17332e-12, 1.27288e-12, 2.34488e-13, 3.65422e-14, 4.94041e-15, 5.96064e-16, 5.92555e-17, 5.36326e-18, 3.87753e-19, 2.74435e-20, 1.45014e-21, 8.48143e-23, 3.07617e-24, 1.3049e-25, 58771.4};
  // ----
  // fall11 weight with S6 PU scenario shifted up by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUFall11Plus05[] = {0.0954358, 0.020682, 0.0484447, 0.47369, 1.50053, 2.01825, 1.99328, 1.80561, 1.61769, 1.54968, 1.51784, 1.3867, 1.12852, 0.773361, 0.467097, 0.277645, 0.157531, 0.0761859, 0.02937, 0.00905285, 0.00233558, 0.000560472, 0.00014654, 4.85625e-05, 1.90022e-05, 7.38414e-06, 2.54485e-06, 7.59523e-07, 2.01995e-07, 5.28236e-08, 1.59273e-08, 5.8689e-09, 2.42327e-09, 9.83602e-10, 3.67109e-10, 1.23672e-10, 3.6646e-11, 9.87495e-12, 2.28817e-12, 4.67308e-13, 8.65061e-14, 1.34004e-14, 1.91938e-15, 2.23009e-16, 2.57593e-17, 2.25592e-18, 2.2207e-19, 1.37666e-20, 1.01363e-21, 50283.5};

  // ----
  // fall11 weight with S6 PU scenario shifted down by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUFall11Minus05[] = {0.0962068, 0.0254305, 0.10942, 1.09375, 2.23819, 2.34231, 2.05248, 1.72372, 1.56859, 1.50191, 1.34308, 1.04368, 0.658236, 0.373059, 0.206629, 0.0997198, 0.0369724, 0.0103136, 0.00227993, 0.000454437, 0.000101297, 3.00008e-05, 1.01979e-05, 3.2025e-06, 8.38403e-07, 1.8544e-07, 3.81727e-08, 8.90286e-09, 2.62892e-09, 8.75387e-10, 2.83582e-10, 8.20241e-11, 2.07262e-11, 4.47805e-12, 8.2373e-13, 1.30017e-13, 1.73391e-14, 2.0282e-15, 1.9709e-16, 1.63193e-17, 1.18443e-18, 6.95731e-20, 3.6548e-21, 1.50639e-22, 5.9703e-24, 1.73528e-25, 5.48351e-27, 1.0555e-28, 2.33406e-30, 73177.8};

  // ----
  // default weight to be used for summer11 samples with S4 PU scenario
  // ----
  // 3.54 fb-1 !!! precale weighted !!! 73.5mb
  //double weightsPUSummer11[] = {0.015788, 0.168999, 0.398126, 0.720134, 1.04936, 1.31397, 1.48381, 1.5613, 1.57478, 1.54969, 1.50736, 1.46098, 1.41641, 1.3734, 1.33783, 1.30218, 1.26742, 1.23224, 1.19459, 1.15701, 1.11803, 1.07132, 1.02626, 0.982213, 0.936123, 0.886997, 0.840387, 0.783362, 0.7366, 0.697965, 0.645112, 0.586587, 0.536603, 0.514893, 0.46368};
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUSummer11[] = {0.022778, 0.231718, 0.517403, 0.888432, 1.23288, 1.47585, 1.59957, 1.62107, 1.57897, 1.50264, 1.41361, 1.32374, 1.23758, 1.15455, 1.07948, 1.00629, 0.936187, 0.868571, 0.802406, 0.739714, 0.679652, 0.618656, 0.562462, 0.510472, 0.460945, 0.413436, 0.370475, 0.326335, 0.289734, 0.259019, 0.225714, 0.193379, 0.166592, 0.150473, 0.127517};
  // ----
  // summer11 weight with S6 PU scenario shifted up by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUSummer11Plus05[] = {0.0181301, 0.190491, 0.439904, 0.780434, 1.11676, 1.37519, 1.52946, 1.58712, 1.58037, 1.53624, 1.47625, 1.4131, 1.35212, 1.29288, 1.24081, 1.18892, 1.13828, 1.08789, 1.03619, 0.985577, 0.93491, 0.879107, 0.82611, 0.775364, 0.724452, 0.67272, 0.624431, 0.570059, 0.524814, 0.486735, 0.440211, 0.391576, 0.350349, 0.32874, 0.289457};

  // ----
  // summer11 weight with S6 PU scenario shifted down by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUSummer11Minus05[] = {0.0285702, 0.28113, 0.606747, 1.0079, 1.35549, 1.57602, 1.66289, 1.64388, 1.56397, 1.45448, 1.33669, 1.22153, 1.11292, 1.01022, 0.917654, 0.829964, 0.748293, 0.672149, 0.600689, 0.535311, 0.475158, 0.417591, 0.366348, 0.320643, 0.279063, 0.241112, 0.208012, 0.176312, 0.150554, 0.12939, 0.108351, 0.0891757, 0.0737799, 0.0639892, 0.0520639};



  switch(sample){
  case kFall11:
    return (nPU < int(sizeof(weightsPUFall11)/sizeof(double)) && nPU > -1) ? weightsPUFall11[nPU] : 0. ;
  case kFall11Plus05:
    return (nPU < int(sizeof(weightsPUFall11Plus05)/sizeof(double)) && nPU > -1) ? weightsPUFall11Plus05[nPU] : 0. ;
  case kFall11Minus05:
    return (nPU < int(sizeof(weightsPUFall11Minus05)/sizeof(double)) && nPU > -1) ? weightsPUFall11Minus05[nPU] : 0. ;
  case kSummer11:
    return (nPU < int(sizeof(weightsPUSummer11)/sizeof(double)) && nPU > -1) ? weightsPUSummer11[nPU] : 0. ;
  case kSummer11Plus05:
    return (nPU < int(sizeof(weightsPUSummer11Plus05)/sizeof(double)) && nPU > -1) ? weightsPUSummer11Plus05[nPU] : 0. ;
  case kSummer11Minus05:
    return (nPU < int(sizeof(weightsPUSummer11Minus05)/sizeof(double)) && nPU > -1) ? weightsPUSummer11Minus05[nPU] : 0. ;
  case kFall10:
    return 1.;
  default:
    return -1.;
  }
}

// calculate the probability of b-tagging one event with 2 b-tags
double eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags){
  double bTaggingEfficiency = 1.;
  double tmp = 1.;

  //std::cout << bTaggingEfficiency << std::endl;

  // subtract probability that no jet is tagged
  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff)
    tmp *= (*eff);
  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis)
    tmp *= (*mis);
  bTaggingEfficiency -= tmp;

  //std::cout << bTaggingEfficiency << std::endl;

  // subtract probability that 1 bJet is tagged
  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
    tmp = 1.-(*eff);
    for(std::vector<double>::const_iterator eff2 = oneMinusBEffies.begin(); eff2 != oneMinusBEffies.end(); ++eff2)
      if(eff != eff2) tmp *= (*eff2);
    for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis)
      tmp *= (*mis);
    bTaggingEfficiency -= tmp;
  }

  //std::cout << bTaggingEfficiency << std::endl;

  // subtract probability that 1 non-bJet is tagged
  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
    tmp = 1.-(*mis);
    for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff)
      tmp *= (*eff);
    for(std::vector<double>::const_iterator mis2 = oneMinusBMistags.begin(); mis2 != oneMinusBMistags.end(); ++mis2)
      if(mis != mis2) tmp *= (*mis2);
    bTaggingEfficiency -= tmp;
  }

  //std::cout << bTaggingEfficiency << std::endl;

  return bTaggingEfficiency;
}

double calcBTagWeight_(int Njet, short * pdgId, TClonesArray * jets)
{
  double bTaggingEfficiency = 0., bTaggingEfficiency_scaled = 0.;
  //double bTaggingEfficiency_scaled_EffUp = 0., bTaggingEfficiency_scaled_EffDown = 0., bTaggingEfficiency_scaled_MisUp = 0., bTaggingEfficiency_scaled_MisDown = 0.;
  double pt, eta, eff, effyScale_pt; //, effVariation_pt, misTagScale_pt, misVariation_pt;

  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0) ; //, oneMinusBEffies_scaled_EffUp(0) , oneMinusBEffies_scaled_EffDown(0) ;
  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0); //, oneMinusBMistags_scaled_MisUp(0), oneMinusBMistags_scaled_MisDown(0);
  for(int i = 0; i < Njet; ++i){
    pt  = ((TLorentzVector*)jets->At(i))->Pt();
    eta = ((TLorentzVector*)jets->At(i))->Eta();

    if(pt > 670.)
      effyScale_pt    = 0.901615*((1.+(0.552628*670.))/(1.+(0.547195*670.)));
    if(pt < 30.)
      effyScale_pt    = 0.901615*((1.+(0.552628*30.))/(1.+(0.547195*30.)));
    else
      effyScale_pt    = 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));
    //effVariation_pt = bTagEffScaleFactor->GetBinError(bTagEffScaleFactor->FindBin(pt));

    if(pdgId[i] == 5 || pdgId[i] == -5){
      eff = bTagEff->GetBinContent(bTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBEffies               .push_back(1.- eff);
      oneMinusBEffies_scaled        .push_back(1.-(eff* effyScale_pt));
      //if(isDefaultSample){
      //  oneMinusBEffies_scaled_EffUp  .push_back(1.-(eff*(effyScale_pt+effVariation_pt)));
      //  oneMinusBEffies_scaled_EffDown.push_back(1.-(eff*(effyScale_pt-effVariation_pt)));
      //}
    }
    else if(pdgId[i] == 4 || pdgId[i] == -4){
      eff = cTagEff->GetBinContent(cTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBMistags               .push_back(1.- eff);
      oneMinusBMistags_scaled        .push_back(1.-(eff* effyScale_pt));
      //if(pt<240){
      //	oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff*(effyScale_pt+(2*effVariation_pt))));
      //	oneMinusBMistags_scaled_MisDown.push_back(1.-(eff*(effyScale_pt-(2*effVariation_pt))));
      //}
      //else{
      //	oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff*(effyScale_pt+effVariation_pt)));
      //	oneMinusBMistags_scaled_MisDown.push_back(1.-(eff*(effyScale_pt-effVariation_pt)));
      //}
    }
    else{
      eff = lTagEff->GetBinContent(lTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBMistags               .push_back(1.- eff);
      oneMinusBMistags_scaled        .push_back(1.-(eff*(((0.948463+(0.00288102*pt))+(-7.98091e-06*(pt*pt)))+(5.50157e-09*(pt*(pt*pt)))) ));
      //if(isDefaultSample){
      //  oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff* ((0.997077+(0.00473953*pt))+(-1.34985e-05*(pt*pt)))+(1.0032e-08*(pt*(pt*pt))) ));
      //  oneMinusBMistags_scaled_MisDown.push_back(1.-(eff* ((0.899715+(0.00102278*pt))+(-2.46335e-06*(pt*pt)))+(9.71143e-10*(pt*(pt*pt))) ));
      //}
    }
  }
  bTaggingEfficiency        = eventBTagProbability_(oneMinusBEffies       , oneMinusBMistags       );
  bTaggingEfficiency_scaled = eventBTagProbability_(oneMinusBEffies_scaled, oneMinusBMistags_scaled);
  //if(isDefaultSample){
  //  bTaggingEfficiency_scaled_EffUp   += eventBTagProbability_(oneMinusBEffies_scaled_EffUp  , oneMinusBMistags_scaled        );
  //  bTaggingEfficiency_scaled_EffDown += eventBTagProbability_(oneMinusBEffies_scaled_EffDown, oneMinusBMistags_scaled        );
  //  bTaggingEfficiency_scaled_MisUp   += eventBTagProbability_(oneMinusBEffies_scaled        , oneMinusBMistags_scaled_MisUp  );
  //  bTaggingEfficiency_scaled_MisDown += eventBTagProbability_(oneMinusBEffies_scaled        , oneMinusBMistags_scaled_MisDown);
  //}
  //std::cout << bTaggingEfficiency << " " << bTaggingEfficiency_scaled << " " << bTaggingEfficiency_scaled/bTaggingEfficiency << std::endl;
  return bTaggingEfficiency_scaled/bTaggingEfficiency;
}

TTree*
modifiedTree(TTree * tree)
{
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("probs", 1);
  tree->SetBranchStatus("dRbb", 1);
  tree->SetBranchStatus("topMasses", 1);
  tree->SetBranchStatus("w1Mass", 1);
  tree->SetBranchStatus("w2Mass", 1);
  tree->SetBranchStatus("comboTypes", 1);
  tree->SetBranchStatus("Njet", 1);
  tree->SetBranchStatus("jets", 1);
  tree->SetBranchStatus("bTag_CSV", 1);
  tree->SetBranchStatus("partonFlavour", 1);
  tree->SetBranchStatus("ttMass", 1);
  //tree->SetBranchStatus("", 1);
  
  const int kMAX = 20000;
  double * probs     = new double[kMAX];
  double * topMasses = new double[kMAX];
  double * w1Masses  = new double[kMAX];
  double * w2Masses  = new double[kMAX];
  float dRbbs;
  double ttMass;
  unsigned short * comboTypes = new unsigned short[kMAX];
  short nPU = -1;

  const int kMAX_(50);
  TString bTagAlgo_ = "CSV";
  int Njet;
  tree->SetBranchStatus ("Njet", 1);
  tree->SetBranchAddress("Njet", &Njet );
  float bTag[kMAX_];
  tree->SetBranchStatus (TString("bTag_")+bTagAlgo_, 1);
  tree->SetBranchAddress(TString("bTag_")+bTagAlgo_, &bTag );
  short pdgId[kMAX_];
  tree->SetBranchStatus ("partonFlavour", 1);       // pdgId  *or*  partonFlavour
  tree->SetBranchAddress("partonFlavour", &pdgId); //  pdgId  *or*  partonFlavour
  TClonesArray * jets = new TClonesArray("TLorentzVector");
  tree->SetBranchStatus("jets", 1);
  tree->GetBranch("jets")->SetAutoDelete(kFALSE);
  tree->SetBranchAddress("jets", &jets);

  double PUWeight = 1.;
  double MCWeight = 1.;
  double BTagWeight = 1.;

  tree->SetBranchAddress("probs",probs);
  tree->SetBranchAddress("dRbb",&dRbbs);
  tree->SetBranchAddress("ttMass",&ttMass);
  tree->SetBranchAddress("topMasses",topMasses);
  tree->SetBranchAddress("w1Mass",w1Masses);
  tree->SetBranchAddress("w2Mass",w2Masses);
  tree->SetBranchAddress("comboTypes",comboTypes);
  //tree->SetBranchAddress("",);

  // event weights
  TString filename = tree->GetCurrentFile()->GetName();
  enumForPUWeights whichSample = kFall10;  
  if(filename.Contains("S11")) whichSample = kSummer11;
  else if(filename.Contains("F11")) whichSample = kFall11;

  if(whichSample == kSummer11){
    tree->SetBranchStatus("nPU", 1);
    tree->SetBranchAddress("nPU", &nPU);
  }
  else if(whichSample == kFall11){
    tree->SetBranchStatus("nPUTru", 1);
    tree->SetBranchAddress("nPUTru", &nPU);
  }
  tree->SetBranchStatus("MCweight", 1);
  tree->SetBranchAddress("MCweight", &MCWeight);

  TString treeName = TString("tree_")+filename;
  treeName.ReplaceAll(samplePath,"");
  treeName.ReplaceAll("_sig.root","");

  double prob,dRbb,topMass,meanWMass,combinedWeight; //w1Mass,w2Mass;
  unsigned short comboType;
  TTree * newTree = new TTree(treeName,treeName);
  newTree->Branch("Njet", &Njet, "Njet/I");
  newTree->Branch("prob", &prob, "prob/D");
  newTree->Branch("topMass", &topMass, "topMass/D");
  newTree->Branch("meanWMass", &meanWMass, "meanWMass/D");
  newTree->Branch("combinedWeight", &combinedWeight, "combinedWeight/D");
  newTree->Branch("jets",&jets, 32000, -1);
  newTree->Branch("dRbb", &dRbb, "dRbb/D");
  newTree->Branch("ttMass", &ttMass, "ttMass/D");
  newTree->Branch("bTag_CSV", bTag, "bTag_CSV[Njet]/F");
  newTree->Branch("comboType", &comboType, "comboType/s");
  //newTree->Branch("w1Mass", &w1Mass, "w1Mass/D");
  //newTree->Branch("w2Mass", &w2Mass, "w2Mass/D");
  //newTree->Branch("", &);

  for(int i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    //if(probs[0] < 0.00) continue;
    if(probs[0] < 0.09) continue;
    if(dRbbs < 1.5) continue;
    if(topMasses[0] < 100.0) continue;
    if(topMasses[0] > 550.0) continue;
    
    prob = probs[0];
    dRbb = dRbbs;
    comboType = comboTypes[0];
    topMass = topMasses[0];
    meanWMass = (w1Masses[0]+w2Masses[0])/2.;
    PUWeight = calcPUWeight_(whichSample, nPU);
    BTagWeight = calcBTagWeight_(Njet, pdgId, jets);
    combinedWeight = PUWeight * MCWeight * BTagWeight; // * probs[0];
    //combinedWeight = 1.;
    //w1Mass = w1Masses[0];
    //w2Mass = w2Masses[0];
    newTree->Fill();
  }
  return newTree;
}

// save plots to disk
void savePlotsToDisk(TString format, bool rootFile)
{
  TSeqCollection * list = gROOT->GetListOfCanvases();
  TString nameOfCanvas;

  TFile * output = TFile::Open("systematicComparison/plots.root", "UPDATE");
  output->cd();

  TString workspaceFolder = output->GetName(); workspaceFolder.ReplaceAll("plots.root","");
  TString folder = workspaceFolder + TString("/plots");
  if(format != "") mkdir(folder, S_IRWXU|S_IRWXG|S_IRWXO);

  int canv = -1;
  do{
    ++canv;
    if(format != "ps" && format != "pdf"){
      nameOfCanvas = TString(list->At(canv)->GetName());
      nameOfCanvas.ReplaceAll("/","_");
    }
    else nameOfCanvas = "plots";
    if(format != ""){
      int errorLevel = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kWarning;
      if     ((format == "ps" || format == "pdf") && list->At(canv) == list->First()) list->At(canv)->Print(folder+TString("/")+nameOfCanvas+TString(".")+format+TString("("));
      else if((format == "ps" || format == "pdf") && list->At(canv) == list->Last())  list->At(canv)->Print(folder+TString("/")+nameOfCanvas+TString(".")+format+TString(")"));
      else                                                                            list->At(canv)->Print(folder+TString("/")+nameOfCanvas+TString(".")+format             );
      gErrorIgnoreLevel = errorLevel;
    }
    if(rootFile) list->At(canv)->Write(nameOfCanvas, TObject::kOverwrite);
  } while(list->At(canv) != list->Last());
  output->Close();
}

void DrawLabel(TString text, const double x1, const double y1, const double x2, Color_t color)
{
  // function to directly draw a label into the active canvas
  double y2 = y1 + 0.05;
  double yOffset = 0.02;
  TPaveLabel *label = new TPaveLabel(x1, y1+yOffset, x2, y2+yOffset, text, "br NDC");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.75);
  label->SetTextAlign(12);
  label->SetTextColor(color);
  label->Draw("same");
}

void DrawCMSPrel() {
  DrawLabel("CMS Preliminary, 3.54 fb^{-1},  #sqrt{s}=7 TeV", 0.11, 0.89, 0.9);
}

void DrawCMSSim() {
  DrawLabel("CMS Simulation Preliminary,  #sqrt{s}=7 TeV", 0.11, 0.89, 0.9);
}

