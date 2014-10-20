#include <iostream>
#include <utility>
#include <vector>

#include "TChain.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TString.h"

double sigmaOfMtop(double mTop = 172.5)
{
  return 245.794 * pow((173.3/mTop),4) * ( 1 - 1.1125 * (mTop-173.3)/173.3 + 0.70778 * pow((mTop-173.3)/173.3,2));
}
double sigmaScaleUpOfMtop(double mTop = 172.5)
{
  return 252.034 * pow((173.3/mTop),4) * ( 1 - 1.11826 * (mTop-173.3)/173.3 + 0.719951 * pow((mTop-173.3)/173.3,2));
}
double sigmaScaleDownOfMtop(double mTop = 172.5)
{
  return 237.375 * pow((173.3/mTop),4) * ( 1 - 1.09562 * (mTop-173.3)/173.3 + 0.677798 * pow((mTop-173.3)/173.3,2));
}
double sigmaPDFUpOfMtop(double mTop = 172.5)
{
  return 251.968 * pow((173.3/mTop),4) * ( 1 - 1.09584 * (mTop-173.3)/173.3 + 0.682769 * pow((mTop-173.3)/173.3,2));
}
double sigmaPDFDownOfMtop(double mTop = 172.5)
{
  return 239.441 * pow((173.3/mTop),4) * ( 1 - 1.12779 * (mTop-173.3)/173.3 + 0.731019 * pow((mTop-173.3)/173.3,2));
}

void mcInfos()
{
  double Lref = 18352.;
  std::string samplePath = "/nfs/dust/cms/user/eschliec/TopMass/2012/Skim_05/";
  std::vector<std::pair<std::string,double> > samples;
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_161_5_sig.root",sigmaOfMtop(161.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_161_5_sig.root",sigmaOfMtop(161.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_161_5_sig.root",sigmaOfMtop(161.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_161_5_sig.root",sigmaOfMtop(161.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_161_5_sig.root",sigmaOfMtop(161.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_163_5_sig.root",sigmaOfMtop(163.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_163_5_sig.root",sigmaOfMtop(163.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_163_5_sig.root",sigmaOfMtop(163.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_163_5_sig.root",sigmaOfMtop(163.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_163_5_sig.root",sigmaOfMtop(163.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_166_5_sig.root",sigmaOfMtop(166.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_166_5_sig.root",sigmaOfMtop(166.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_166_5_sig.root",sigmaOfMtop(166.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_166_5_sig.root",sigmaOfMtop(166.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_166_5_sig.root",sigmaOfMtop(166.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_169_5_sig.root",sigmaOfMtop(169.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_169_5_sig.root",sigmaOfMtop(169.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_169_5_sig.root",sigmaOfMtop(169.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_169_5_sig.root",sigmaOfMtop(169.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_169_5_sig.root",sigmaOfMtop(169.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_175_5_sig.root",sigmaOfMtop(175.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_175_5_sig.root",sigmaOfMtop(175.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_175_5_sig.root",sigmaOfMtop(175.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_175_5_sig.root",sigmaOfMtop(175.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_175_5_sig.root",sigmaOfMtop(175.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_178_5_sig.root",sigmaOfMtop(178.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_178_5_sig.root",sigmaOfMtop(178.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_178_5_sig.root",sigmaOfMtop(178.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_178_5_sig.root",sigmaOfMtop(178.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_178_5_sig.root",sigmaOfMtop(178.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_181_5_sig.root",sigmaOfMtop(181.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_181_5_sig.root",sigmaOfMtop(181.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_181_5_sig.root",sigmaOfMtop(181.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_181_5_sig.root",sigmaOfMtop(181.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_181_5_sig.root",sigmaOfMtop(181.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_184_5_sig.root",sigmaOfMtop(184.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_184_5_sig.root",sigmaOfMtop(184.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_184_5_sig.root",sigmaOfMtop(184.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_184_5_sig.root",sigmaOfMtop(184.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_184_5_sig.root",sigmaOfMtop(184.5)*Lref));

  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_096_172_5_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_098_172_5_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_172_5_sig.root",sigmaScaleDownOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_172_5_sig.root",sigmaPDFDownOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_172_5_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_172_5_sig.root",sigmaPDFUpOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_172_5_sig.root",sigmaScaleUpOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_102_172_5_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_104_172_5_sig.root",sigmaOfMtop(172.5)*Lref));

  //samples.push_back(std::make_pair("Z2_S12_ABS_JES_100_172_5_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_MassiveBinDecay_ABS_JES_100_172_5_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_BJES_Down_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_BJES_Up_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_CORFLAVQUARKJES_Down_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_CORFLAVQUARKJES_Up_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_CORTOTNOFLAVJES_Down_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_CORTOTNOFLAVJES_Up_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_JER_Down_sig.root",sigmaOfMtop(172.5)*Lref));
  //samples.push_back(std::make_pair("Z2_S12_JER_Up_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_Matching_Down_MadSpin_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_Matching_Up_MadSpin_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_Scale_Down_MadSpin_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_Scale_Up_MadSpin_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_P11_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_P11mpiHi_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_P11TeV_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_P11NoCR_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_POWHEG_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_POWHER_sig.root",sigmaOfMtop(172.5)*Lref));
  samples.push_back(std::make_pair("Z2_S12_MCNLO_sig.root",sigmaOfMtop(172.5)*Lref));

  double nData1   = 32072.;
  double nDataAll = 76973.;

  std::cout << "sample: total events, selected events, selection efficiency [%], signal fraction [%], fCP [%], fWP [%], fUN [%]" << std::endl;
  for(auto& sample : samples){
    TChain* chain = new TChain("analyzeKinFit/eventTree");

    chain->Add((samplePath+sample.first).c_str());
    TString draw1 = "top.combinationType[0]";
    TString weight1 = "weight.combinedWeight*top.fitProb[0]"; // temporary divide by BR (/0.456976), due to wrong weights in sample
    TString selection1 = "*(top.fitProb[0] > 0.1 && (pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2)) > 1.5*1.5 && max(max(max(max(max(recoJetIdxB1[0],recoJetIdxW1Prod1[0]),recoJetIdxW1Prod2[0]),recoJetIdxB2[0]),recoJetIdxW2Prod1[0]),recoJetIdxW2Prod2[0])<99)";

    TString drawAll = draw1; drawAll.ReplaceAll("[0]","");
    TString weightAll = weight1; weightAll.ReplaceAll("[0]","");
    TString selectionAll = selection1; selectionAll.ReplaceAll("[0]","");

    chain->Draw(draw1  +TString(">>hist1(20,-10,10)"),weight1  +selection1  , "goff");
    chain->Draw(drawAll+TString(">>histA(20,-10,10)"),weightAll+selectionAll, "goff");
    TH1F *hist1   = (TH1F*)gDirectory->Get("hist1");
    double integral1 = hist1->Integral(0,hist1->GetNbinsX()+1);
    double fCP1 = hist1->GetBinContent(hist1->FindBin(1.0)) / integral1;
    double fWP1 = hist1->Integral(0,hist1->FindBin(0.0)) / integral1;
    double fUN1 = hist1->Integral(hist1->FindBin(2.0),hist1->GetNbinsX()+1) / integral1;
    //integral *= sample.second/sigmaOfMtop(172.5)/Lref;
    std::cout.width(60);
    std::cout << sample.first << " 1. : ";
    std::cout.precision(1);
    std::cout << std::fixed << sample.second << ", " << integral1 << ", ";
    std::cout.precision(3);
    std::cout << integral1/sample.second * 100 << ", ";
    std::cout.precision(2);
    std::cout << integral1/nData1 * 100 << ", " << fCP1 * 100 << ", " << fWP1 * 100 << ", " << fUN1 * 100 << std::endl;

    TH1F *histAll = (TH1F*)gDirectory->Get("histA");
    double integralAll = histAll->Integral(0,histAll->GetNbinsX()+1);
    double fCPAll = histAll->GetBinContent(histAll->FindBin(1.0)) / integralAll;
    double fWPAll = histAll->Integral(0,histAll->FindBin(0.0)) / integralAll;
    double fUNAll = histAll->Integral(histAll->FindBin(2.0),histAll->GetNbinsX()+1) / integralAll;
    //integral *= sample.second/sigmaOfMtop(172.5)/Lref;
    std::cout.width(60);
    std::cout << sample.first << " All: ";
    std::cout.precision(1);
    std::cout << std::fixed << sample.second << ", " << integralAll << ", ";
    std::cout.precision(3);
    //std::cout << integralAll/sample.second * 100 << ", ";
    std::cout << "X.XXX" << ", ";
    std::cout.precision(2);
    std::cout << integralAll/nDataAll * 100 << ", " << fCPAll * 100 << ", " << fWPAll * 100 << ", " << fUNAll * 100 << std::endl;

  }

}
