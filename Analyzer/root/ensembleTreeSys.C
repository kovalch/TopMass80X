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

//enum sysVariation {kNominal, kFlavorDown, kFlavorUp};
TString sPath("/scratch/hh/current/cms/user/mseidel/topmass_120412_2120cp/all/");

double mass[11];
double jes [11];

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
  double totalMassUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  
  ensembles.push_back(ensemble("Fall11_TTJets1725_1.00/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_flavor:down/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_flavor:up/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jes:down/ensemble.root", true, 0.985));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jes:up/ensemble.root", true, 1.015));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jer:down/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_jer:up/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_matchingup/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_matchingup/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_scaledown/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_scaleup/ensemble.root"));
  /*
  ensembles.push_back(ensemble("Fall11_TTJets1725_METCL1/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_METCL1/ensemble.root"));
  */
  //*
  ensembles.push_back(ensemble("Fall11_TTJets1725_P11/ensemble.root"));
  ensembles.push_back(ensemble("Fall11_TTJets1725_P11/ensemble.root"));
  //*/
  /*
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeightDown/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeightUp/ensemble.root"));
  ensembles.push_back(ensemble("fSig_0.92/ensemble.root"));
  ensembles.push_back(ensemble("fSig_0.92/ensemble.root"));
  ensembles.push_back(ensemble("wbb_fSig_0.96/ensemble.root"));
  ensembles.push_back(ensemble("wbb_fSig_0.96/ensemble.root"));
  ensembles.push_back(ensemble("bDisc_0.61/ensemble.root"));
  ensembles.push_back(ensemble("bDisc_0.75/ensemble.root"));
  //*/
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    files.push_back(new TFile(sPath + ensembles[i].file));
    trees.push_back((TTree*) files[i]->Get("tree"));
  }
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    TF1* gaus = new TF1("gaus", "gaus");
    
    trees[i]->Fit("gaus", "mass", "mass>0 & JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].mass = gaus->GetParameter(1);
    
    trees[i]->Fit("gaus", "JES", "mass>0 & JES>0 & genMass==172.5 & genJES==1", "EMQ0");
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


