#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TFitResult.h"

#include "tdrstyle.C"

enum styles          { kDown, kNominal, kUp};
int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
int marker_  [ 3 ] = { 23, 20, 22};

struct ensemble {
  const char* file;
  bool takeLargest;
  double expectedJES;
  double mass;
  double jes;
  ensemble(const char* f, bool t = true, double j = 0)
  : file(f), takeLargest(t), expectedJES(j) {}
};

std::vector<ensemble> ensembles;

//TString sPath("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/topmass_130913_1901c/");
TString sPath("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/topmass_131213_1201/");
//TString sPath("");

double mass[11];
double jes [11];

std::vector<TFile*> files;
std::vector<TTree*> trees;


void ensembleTreeSys()
{
  double totalMassUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  
  ensembles.push_back(ensemble("ensemble_S12_1725_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jesTotNoFlav_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jesTotNoFlav_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jesFlavorQuark_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jesFlavorQuark_up.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jesTotNoFlav_down.root", true, 0.988));
  ////ensembles.push_back(ensemble("ensemble_S12_jesTotNoFlav_up.root", true, 1.012));
  ////ensembles.push_back(ensemble("ensemble_S12_jesFlavorUds_down.root", true, 0.988));
  ////ensembles.push_back(ensemble("ensemble_S12_jesFlavorUds_up.root", true, 1.012));
  //ensembles.push_back(ensemble("ensemble_S12_jer_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jer_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_match_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_match_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_scale_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_scale_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_bjes_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_bjes_up.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_fsig_down.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_fsig_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_pu_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_pu_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_btagSF_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_btagSF_up.root"));
  //ensembles.push_back(ensemble("ensemble_S12_mtagSF_down.root"));
  //ensembles.push_back(ensemble("ensemble_S12_mtagSF_up.root"));
  //////ensembles.push_back(ensemble("ensemble_S12_P11.root"));
  //////ensembles.push_back(ensemble("ensemble_S12_P11.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_shape_095.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_shape_105.root"));
  //////ensembles.push_back(ensemble("ensemble_S12_shape_090.root"));
  //////ensembles.push_back(ensemble("ensemble_S12_shape_110.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+20.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+20.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-20.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-20.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+30.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+30.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-30.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-30.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+40.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_+40.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-40.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_permu_-40.root"));
  ////ensembles.push_back(ensemble("ensemble_S11.root"));
  ////ensembles.push_back(ensemble("ensemble_S11.root"));
  ensembles.push_back(ensemble("ensemble_S12_POWHEG.root"));
  ensembles.push_back(ensemble("ensemble_S12_POWHEG.root"));
  //ensembles.push_back(ensemble("ensemble_S12_POWHER.root"));
  //ensembles.push_back(ensemble("ensemble_S12_POWHER.root"));
  //ensembles.push_back(ensemble("ensemble_S12_MCNLO.root"));
  //ensembles.push_back(ensemble("ensemble_S12_MCNLO.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jet4_02.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jet4_02.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet4_05.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet4_05.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet4_10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet4_10.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jet5_02.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jet5_02.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet5_05.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet5_05.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet5_10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet5_10.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jet6_02.root"));
  //ensembles.push_back(ensemble("ensemble_S12_jet6_02.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet6_05.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet6_05.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet6_10.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_jet6_10.root"));

  //ensembles.push_back(ensemble("ensemble_S12_1725_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_BackgroundModel.root"));
  //ensembles.push_back(ensemble("ensemble_S12_BackgroundModel.root"));
  //ensembles.push_back(ensemble("ensemble_S12_BackgroundModel2.root"));
  //ensembles.push_back(ensemble("ensemble_S12_BackgroundModel2.root"));

  //ensembles.push_back(ensemble("ensemble_S12_P11.root"));
  //ensembles.push_back(ensemble("ensemble_S12_P11_NoCR.root"));
  //ensembles.push_back(ensemble("ensemble_S12_P11_NoCR.root"));
  //ensembles.push_back(ensemble("ensemble_S12_P11mpiHi.root"));
  //ensembles.push_back(ensemble("ensemble_S12_P11TeV.root"));

  //ensembles.push_back(ensemble("ensemble_S12_Z2_FAST.root"));
  //ensembles.push_back(ensemble("ensemble_S12_Z2_Matching_Down_FAST.root"));
  //ensembles.push_back(ensemble("ensemble_S12_Z2_Matching_Down_FAST.root"));
  //ensembles.push_back(ensemble("ensemble_S12_Z2_096_FAST.root"));
  //ensembles.push_back(ensemble("ensemble_S12_Z2_096_FAST.root"));
  //ensembles.push_back(ensemble("ensemble_S12_Z2_104_FAST.root"));
  //ensembles.push_back(ensemble("ensemble_S12_Z2_104_FAST.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_Z2_fragCor_FAST.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_Z2_fragCor_FAST.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_Z2_fragP11_FAST.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_Z2_fragP11_FAST.root"));

  //ensembles.push_back(ensemble("ensemble_S12_1725_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1615_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1615_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1615_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1615_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1615_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1635_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1635_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1635_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1635_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1635_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1665_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1665_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1665_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1665_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1665_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1695_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1695_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1695_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1695_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1695_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1725_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1725_098.root"));
  ////ensembles.push_back(ensemble("ensemble_S12_1725_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1725_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1725_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1755_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1755_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1755_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1755_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1755_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1785_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1785_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1785_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1785_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1785_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1815_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1815_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1815_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1815_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1815_104.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1845_096.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1845_098.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1845_100.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1845_102.root"));
  //ensembles.push_back(ensemble("ensemble_S12_1845_104.root"));

  for (int i = 0; i < (int) ensembles.size(); ++i) {
    files.push_back(new TFile(sPath + ensembles[i].file));
    trees.push_back((TTree*) files[i]->Get("tree"));
  }
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    TF1* gaus = new TF1("gaus", "gaus");
    
    //trees[i]->Fit("gaus", "mass_mTop_JES", "mass_mTop_JES>0 & JES_mTop_JES>0", "EMQ0");
    trees[i]->Fit("gaus", "mass_mTop_JES", "mass_mTop_JES>0 & JES_mTop_JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].mass = gaus->GetParameter(1);
    
    //trees[i]->Fit("gaus", "JES_mTop_JES", "mass_mTop_JES>0 & JES_mTop_JES>0", "EMQ0");
    trees[i]->Fit("gaus", "JES_mTop_JES", "mass_mTop_JES>0 & JES_mTop_JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].jes = gaus->GetParameter(1);
  }
  
  for (int i = 1; i < (int) ensembles.size(); i+=2) {
    std::cout << ensembles[i].file << " / " << ensembles[i+1].file << std::endl;
    
    printf("\t- %4.3f GeV / o %4.3f GeV / + %4.3f GeV \n", ensembles[i].mass, ensembles[0].mass, ensembles[i+1].mass);
    double largestDM = max(abs(ensembles[0].mass-ensembles[i].mass), abs(ensembles[0].mass-ensembles[i+1].mass));
    double meanDM    = (abs(ensembles[0].mass-ensembles[i].mass) + abs(ensembles[0].mass-ensembles[i+1].mass))/2;
    printf("\tLargest uncertainty: %4.2f GeV / Mean uncertainty: %4.2f GeV \n", largestDM, meanDM);
    
    (ensembles[i].takeLargest) ? totalMassUncertainty2 += largestDM*largestDM : totalMassUncertainty2 += meanDM*meanDM;
    
    printf("\t- %4.4f / o %4.4f / + %4.4f \n", ensembles[i].jes, ensembles[0].jes, ensembles[i+1].jes);
    double largestDJ = 0;
    double meanDJ    = 0;
    if (ensembles[i].expectedJES > 0) {
      largestDJ = max(abs(ensembles[i].expectedJES-ensembles[i].jes), abs(ensembles[i+1].expectedJES-ensembles[i+1].jes));
      meanDJ    = (abs(ensembles[i].expectedJES-ensembles[i].jes) + abs(ensembles[i+1].expectedJES-ensembles[i+1].jes))/2;
    }
    else {
      largestDJ = max(abs(ensembles[0].jes-ensembles[i].jes), abs(ensembles[0].jes-ensembles[i+1].jes));
      meanDJ    = (abs(ensembles[0].jes-ensembles[i].jes) + abs(ensembles[0].jes-ensembles[i+1].jes))/2;
    }
    printf("\tLargest uncertainty: %4.3f / Mean uncertainty: %4.3f \n\n", largestDJ, meanDJ);
    
    (ensembles[i].takeLargest) ? totalJESUncertainty2 += largestDJ*largestDJ : totalJESUncertainty2 += meanDJ*meanDJ;
  }
  
  printf("Systematic uncertainty on top mass: %4.2f GeV \n", sqrt(totalMassUncertainty2));
  printf("Systematic uncertainty on JES     : %4.3f \n", sqrt(totalJESUncertainty2));
}


