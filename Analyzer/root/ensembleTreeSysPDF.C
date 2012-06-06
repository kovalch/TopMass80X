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

enum lepton           { kElectron, kMuon, kAll};
TString lepton_ [3] = { "electron", "muon", "all"};

int channel = 2;

struct ensemble {
  const char* file;
  bool takeLargest;
  double expectedJES;
  int reference;
  double mass;
  double jes;
  ensemble(const char* f, bool t = true, double j = 0, int r = 0)
  : file(f), takeLargest(t), expectedJES(j), reference(r) {}
};

TString path("/scratch/hh/current/cms/user/mseidel/topmass_120412_2120cp/");

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
  TString sPath = path; sPath += lepton_[channel]; sPath += "/";
  
  std::vector<ensemble> ensembles;
  
  double totalMassUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  
  
  ensembles.push_back(ensemble("Fall11_TTJets1725_1.00/ensemble.root"));
  // >>> for i in range(4): print "ensembles.push_back(ensemble(\"muWeight-bWeight-PUWeight-pdfWeights[" + str(i) + "]/ensemble.root\"));"
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[0]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[1]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[2]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[3]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[4]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[5]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[6]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[7]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[8]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[9]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[10]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[11]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[12]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[13]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[14]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[15]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[16]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[17]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[18]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[19]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[20]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[21]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[22]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[23]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[24]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[25]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[26]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[27]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[28]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[29]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[30]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[31]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[32]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[33]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[34]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[35]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[36]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[37]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[38]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[39]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[40]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[41]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[42]/ensemble.root"));
  ensembles.push_back(ensemble("muWeight-bWeight-PUWeight-pdfWeights[43]/ensemble.root"));

  


  
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
    double meanDM    = (ensembles[i].mass - ensembles[i+1].mass)/2;
    printf("\tSymmetrized uncertainty: %4.2f GeV \n", meanDM);
    
    totalMassUncertainty2 += meanDM*meanDM;
    
    printf("\t- %4.4f / o %4.4f / + %4.4f \n", ensembles[i].jes, ensembles[0].jes, ensembles[i+1].jes);
    double meanDJ    = 0;
    
    meanDJ    = (ensembles[i].jes - ensembles[i+1].jes)/2;
    printf("\tSymmetrized uncertainty: %4.3f \n\n", meanDJ);
    
    totalJESUncertainty2 += meanDJ*meanDJ;
  }
  
  printf("Systematic uncertainty on top mass: %4.2f GeV \n", sqrt(totalMassUncertainty2));
  printf("Systematic uncertainty on JES     : %4.3f \n", sqrt(totalJESUncertainty2));
}


