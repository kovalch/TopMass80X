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
#include "TString.h"

#include "tdrstyle.C"

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 2;

struct ensemble {
  const char* file;
  double size; // EFFECTIVE sample size
  bool correlated;
  bool takeLargest;
  double expectedJES;
  int reference;
  double mass;
  double massWidth;
  double jes;
  double jesWidth;
  ensemble(const char* f, double s, bool c = true, bool t = true, double j = 0, int r = 0)
  : file(f), size(s), correlated(c), takeLargest(t), expectedJES(j), reference(r) {}
};

struct staticUncertainty {
  const char* name;
  double massUncertainty;
  double jesUncertainty;
  staticUncertainty(const char* n, double mu, double ju)
  : name(n), massUncertainty(mu), jesUncertainty(ju) {}
};

TString globalPath("/scratch/hh/current/cms/user/mseidel/topmass_120412_2120cp/");

double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
double genMassN[]     = {1620072, 1633197, 1669034, 1606570, 59613991, 1538301, 1648519, 1665350, 1671859};
//double genMassN[]     = {1620072, 1.5, 1.5, 1.5, 59613991, 1.5, 1.5, 1.5, 1.5};
double maxMCWeight[]  = {1.7, 2.2, 2.2, 2.2, 1.7, 2.2, 2.2, 2.2, 1.7};
double crossSection   = 164.4;
double peLumi         = 5000.;
  

std::vector<TFile*> files;
std::vector<TTree*> trees;


void ensembleTreeSys()
{
  TString sPath = globalPath; sPath += lepton_[channel]; sPath += "/";
  
  std::vector<ensemble> ensembles;
  std::vector<staticUncertainty> staticUncertainties;
  
  double totalMassUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  
  switch(channel) {
    case kElectron:
      staticUncertainties.push_back(staticUncertainty("Calibration", 0.09, 0.001));
      staticUncertainties.push_back(staticUncertainty("PDF", 0.07, 0.001));
      staticUncertainties.push_back(staticUncertainty("UE", 0.244, 0.0011));
      break;
    case kMuon:
      staticUncertainties.push_back(staticUncertainty("Calibration", 0.08, 0.001));
      staticUncertainties.push_back(staticUncertainty("PDF", 0.07, 0.001));
      staticUncertainties.push_back(staticUncertainty("UE", 0.256, 0.0025));
      break;
    case kAll:
      staticUncertainties.push_back(staticUncertainty("Calibration", 0.06, 0.001));
      staticUncertainties.push_back(staticUncertainty("PDF", 0.07, 0.001));
      staticUncertainties.push_back(staticUncertainty("UE", 0.187, 0.0018)); // 0.153 +/- 0.187
      staticUncertainties.push_back(staticUncertainty("Background", 0.126, 0.001));
      staticUncertainties.push_back(staticUncertainty("Statistical", 0.428, 0.003));
      break;
  }
  
  ensembles.push_back(ensemble("Fall11_TTJets1725_1.00/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_flavor:down/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_flavor:up/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jes:down/ensemble.root", 59613991./1.7, true, true, 0.984));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jes:up/ensemble.root", 59613991./1.7, true, true, 1.016));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jer:down/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jer:up/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_matchingup/ensemble.root", 4029823./1.7, false));
  ensembles.push_back(ensemble("Fall11_TTJets1725_matchingup/ensemble.root", 4029823./1.7, false));
  ensembles.push_back(ensemble("Fall11_TTJets1725_scaledown/ensemble.root", 4004587./1.7, false));
  ensembles.push_back(ensemble("Fall11_TTJets1725_scaleup/ensemble.root", 3696269./1.7, false));
  ensembles.push_back(ensemble("Fall11_TTJets1725_P11noCR/ensemble.root", 8841768./1.7, false, true, 0., 12));
  ensembles.push_back(ensemble("Fall11_TTJets1725_P11/ensemble.root", 8621453./1.7, false));
  /*
  ensembles.push_back(ensemble("Fall11_TTJets1725_powheg/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_mcatnlo/ensemble.root"));
  //*/
  //*/
  //*
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeightDown/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeightUp/ensemble.root", 59613991./1.7));
  if (!(channel == kAll)) {
    ensembles.push_back(ensemble("fSig_0.92/ensemble.root", 59613991./1.7));
    ensembles.push_back(ensemble("fSig_0.84/ensemble.root", 59613991./1.7));
  }
  ensembles.push_back(ensemble("bDisc_0.61/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("bDisc_0.75/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_EES_down/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_EES_up/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_MES_down/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_MES_up/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_UNC_0.9/ensemble.root", 59613991./1.7));
  ensembles.push_back(ensemble("Fall11_TTJets1725_UNC_1.1/ensemble.root", 59613991./1.7));
  //*/
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    files.push_back(new TFile(sPath + ensembles[i].file));
    trees.push_back((TTree*) files[i]->Get("tree"));
  }
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    TF1* gaus = new TF1("gaus", "gaus");
    
    trees[i]->Fit("gaus", "mass", "mass>0 & JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].mass = gaus->GetParameter(1);
    ensembles[i].massWidth = gaus->GetParameter(2);
    
    trees[i]->Fit("gaus", "JES", "mass>0 & JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].jes = gaus->GetParameter(1);
    ensembles[i].jesWidth = gaus->GetParameter(2);
  }
  
  for (int i = 1; i < (int) ensembles.size(); i+=2) {
    std::cout << "\n" << ensembles[i].file << " / " << ensembles[i+1].file << std::endl;
    
    printf("\t- %4.3f GeV / o %4.3f GeV / + %4.3f GeV \n", ensembles[i].mass, ensembles[0].mass, ensembles[i+1].mass);
    double largestDM = max(abs(ensembles[ensembles[i].reference].mass-ensembles[i].mass), abs(ensembles[ensembles[i].reference].mass-ensembles[i+1].mass));
    double meanDM    = (abs(ensembles[ensembles[i].reference].mass-ensembles[i].mass) + abs(ensembles[ensembles[i].reference].mass-ensembles[i+1].mass))/2;
    printf("\tLargest uncertainty: %4.2f GeV / Mean uncertainty: %4.2f GeV", largestDM, meanDM);
    
    double statDM = sqrt(pow(ensembles[ensembles[i].reference].massWidth / sqrt(ensembles[ensembles[i].reference].size/(crossSection*peLumi)), 2) + pow(ensembles[i].massWidth / sqrt(ensembles[i].size/(crossSection*peLumi)), 2));
    printf(" / Statistical precision: %4.2f GeV \n", statDM);
    if (statDM > largestDM && !ensembles[i].correlated) largestDM = statDM;
    
    (ensembles[i].takeLargest) ? totalMassUncertainty2 += largestDM*largestDM : totalMassUncertainty2 += meanDM*meanDM;
    
    printf("\t- %4.4f / o %4.4f / + %4.4f \n", ensembles[i].jes, ensembles[0].jes, ensembles[i+1].jes);
    double largestDJ = 0;
    double meanDJ    = 0;
    if (ensembles[i].expectedJES > 0) {
      largestDJ = max(abs(ensembles[i].expectedJES-ensembles[i].jes), abs(ensembles[i+1].expectedJES-ensembles[i+1].jes));
      meanDJ    = (abs(ensembles[i].expectedJES-ensembles[i].jes) + abs(ensembles[i+1].expectedJES-ensembles[i+1].jes))/2;
    }
    else {
      largestDJ = max(abs(ensembles[ensembles[i].reference].jes-ensembles[i].jes), abs(ensembles[ensembles[i].reference].jes-ensembles[i+1].jes));
      meanDJ    = (abs(ensembles[ensembles[i].reference].jes-ensembles[i].jes) + abs(ensembles[ensembles[i].reference].jes-ensembles[i+1].jes))/2;
    }
    printf("\tLargest uncertainty: %4.3f / Mean uncertainty: %4.3f", largestDJ, meanDJ);
    
    double statDJ = sqrt(pow(ensembles[ensembles[i].reference].jesWidth / sqrt(ensembles[ensembles[i].reference].size/(crossSection*peLumi)), 2) + pow(ensembles[i].jesWidth / sqrt(ensembles[i].size/(crossSection*peLumi)), 2));
    printf(" / Statistical precision: %4.3f \n", statDJ);
    if (statDJ > largestDJ && !ensembles[i].correlated) largestDJ = statDJ;
    
    (ensembles[i].takeLargest) ? totalJESUncertainty2 += largestDJ*largestDJ : totalJESUncertainty2 += meanDJ*meanDJ;
  }
  
  for (int i = 0; i < (int) staticUncertainties.size(); ++i) {
    std::cout << staticUncertainties[i].name << std::endl;
    printf("\n\tMass uncertainty: %4.2f / JES uncertainty: %4.3f \n\n", staticUncertainties[i].massUncertainty, staticUncertainties[i].jesUncertainty);
    
    totalMassUncertainty2 += staticUncertainties[i].massUncertainty*staticUncertainties[i].massUncertainty;
    totalJESUncertainty2 += staticUncertainties[i].jesUncertainty*staticUncertainties[i].jesUncertainty;
  }
  
  printf("Systematic uncertainty on top mass: %4.2f GeV \n", sqrt(totalMassUncertainty2));
  printf("Systematic uncertainty on JES     : %4.3f \n", sqrt(totalJESUncertainty2));
}


