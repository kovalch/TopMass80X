#include <vector>
#include <iostream>
#include <sstream> 

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

struct ensemble {
  const char* file;
  bool takeLargest;
  double expectedJES;
  int reference;
  double mass;
  double massConstJES;
  double jes;
  ensemble(const char* f, bool t = true, double j = 0, int r = 0)
  : file(f), takeLargest(t), expectedJES(j), reference(r) {}
};

//TString path("/scratch/hh/dust/naf/cms/user/eschliec/topmass_120522_1500/");
//TString path("/scratch/hh/current/cms/user/eschliec/topmass_120522_1500/");
TString path("/scratch/hh/current/cms/user/eschliec/topmass_120810_1202a/");

double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
double genMassN[]     = {1620072, 1633197, 1669034, 1606570, 59613991, 1538301, 1648519, 1665350, 1671859};
//double genMassN[]     = {1620072, 1.5, 1.5, 1.5, 59613991, 1.5, 1.5, 1.5, 1.5};
double maxMCWeight[]  = {1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7};
double crossSection   = 164.4;
double peLumi         = 3544.844;
  

std::vector<TFile*> files;
std::vector<TTree*> trees;


void ensembleTreeSysPDF()
{
  std::vector<ensemble> ensembles;
  
  double totalMassUncertainty2 = 0;
  double totalMassConstJESUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  
  
  ensembles.push_back(ensemble("ensemble_F11_1725_100.root"));
  // >>> for i in range(4): print "ensembles.push_back(ensemble(\"[" + str(i) + "]/\"));"
  ensembles.push_back(ensemble("ensemble_F11_pdf01_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf01_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf02_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf02_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf03_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf03_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf04_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf04_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf05_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf05_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf06_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf06_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf07_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf07_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf08_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf08_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf09_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf09_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf10_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf10_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf11_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf11_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf12_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf12_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf13_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf13_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf14_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf14_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf15_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf15_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf16_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf16_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf17_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf17_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf18_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf18_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf19_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf19_down.root"));
  ensembles.push_back(ensemble("ensemble_F11_pdf20_up.root  "));
  ensembles.push_back(ensemble("ensemble_F11_pdf20_down.root"));
 
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    files.push_back(new TFile(path + ensembles[i].file));
    trees.push_back((TTree*) files[i]->Get("tree"));
  }
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    TF1* gaus = new TF1("gaus", "gaus");
    
    trees[i]->Fit("gaus", "mass", "mass>0 & JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].mass = gaus->GetParameter(1);
    
    trees[i]->Fit("gaus", "JES", "mass>0 & JES>0 & genMass==172.5 & genJES==1", "EMQ0");
    ensembles[i].jes = gaus->GetParameter(1);

    trees[i]->Fit("gaus", "massConstJES", "massConstJES>0 & genMass==172.5", "EMQ0");
    ensembles[i].massConstJES = gaus->GetParameter(1);   
  }
  
  for (int i = 1; i < (int) ensembles.size(); i+=2) {
    std::cout << ensembles[i].file << " / " << ensembles[i+1].file << std::endl;
    
    printf("\t- %4.3f GeV / o %4.3f GeV / + %4.3f GeV \n", ensembles[i].mass, ensembles[0].mass, ensembles[i+1].mass);
    double meanDM    = (ensembles[i].mass - ensembles[i+1].mass)/2;
    printf("\tSymmetrized uncertainty: %4.2f GeV \n", meanDM);
    
    totalMassUncertainty2 += meanDM*meanDM;
    
    printf("\t- %4.4f / o %4.4f / + %4.4f \n", ensembles[i].jes, ensembles[0].jes, ensembles[i+1].jes);
    double meanDJ    = 0;
    
    meanDJ    = (ensembles[i].jes - ensembles[i+1].jes)/2;
    printf("\tSymmetrized uncertainty: %4.3f \n\n", meanDJ);
    
    totalJESUncertainty2 += meanDJ*meanDJ;

    printf("\t- %4.3f GeV / o %4.3f GeV / + %4.3f GeV \n", ensembles[i].massConstJES, ensembles[0].massConstJES, ensembles[i+1].massConstJES);
    double meanDMConstJES    = (ensembles[i].massConstJES - ensembles[i+1].massConstJES)/2;
    printf("\tSymmetrized uncertainty: %4.2f GeV \n", meanDMConstJES);
    
    totalMassConstJESUncertainty2 += meanDMConstJES*meanDMConstJES;    
  }
  
  printf("Systematic uncertainty on top mass: %4.2f GeV \n", sqrt(totalMassUncertainty2));
  printf("Systematic uncertainty on JES     : %4.3f \n", sqrt(totalJESUncertainty2));
  printf("Systematic uncertainty on mTop(1D): %4.2f GeV \n", sqrt(totalMassConstJESUncertainty2));
}

