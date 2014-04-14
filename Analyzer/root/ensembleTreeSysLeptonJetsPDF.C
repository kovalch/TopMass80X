#include <vector>
#include <iostream>
#include <sstream> 

#include "TFile.h"
#include "TChain.h"
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
  std::string file;
  bool takeLargest;
  double expectedJES;
  int reference;
  double mass;
  double mass1d;
  double jes;
  TChain* chain;
  ensemble(std::string f, bool t = true, double j = 0, int r = 0)
  : file(f), takeLargest(t), expectedJES(j), reference(r) {}
};

std::string path("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_140305/lepton");
//std::string path("/nfs/dust/cms/user/eschliec/TopMass/topmass_140401_1201c/Z2_S12_ABS_JES_100_172_5");

std::string itos(int number)
{
   stringstream ss;
   ss << number;
   return ss.str();
}

void ensembleTreeSysLeptonJetsPDF()
{
  std::vector<ensemble> ensembles;
  
  TString massName   = "mass_mTop_JES";
  TString  jesName   =  "JES_mTop_JES";
  TString mass1dName = "mass_mTop";

  //TString massName   = "mass_mTop_JES_fSig_fCP";
  //TString  jesName   =  "JES_mTop_JES_fSig_fCP";
  //TString mass1dName = "mass_mTop";

  double massMin =  1000.;
  double massMax = -1000.;
  
  for (int i = 0; i < 147; ++i) {
    ensembles.push_back(ensemble(std::string("/pdf_weight/")+itos(i)+std::string("/job_*.root")));
  }
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {    
    ensembles[i].chain = new TChain("tree");
    ensembles[i].chain->Add((path + ensembles[i].file).c_str());
  }
  
  for (int i = 0; i < (int) ensembles.size(); ++i) {
    TF1* gaus = new TF1("gaus", "gaus");
    
    ensembles[i].chain->Fit("gaus", massName, massName+">0 &"+jesName+">0 & genMass==172.5 & genJES==1", "LEMQ0");
    ensembles[i].mass = gaus->GetParameter(1);
    std::cout << i << " mass: " << ensembles[i].mass << std::endl;
    if (massMin>ensembles[i].mass) massMin = ensembles[i].mass;
    if (massMax<ensembles[i].mass) massMax = ensembles[i].mass;

    ensembles[i].chain->Fit("gaus", jesName, massName+">0 &"+jesName+">0 & genMass==172.5 & genJES==1", "LEMQ0");
    ensembles[i].jes = gaus->GetParameter(1);

    ensembles[i].chain->Fit("gaus", mass1dName, mass1dName+">0 & genMass==172.5 & genJES==1", "LEMQ0");
    ensembles[i].mass1d = gaus->GetParameter(1);
  }
  
  std::cout << "CTEQ6l1: " << 172.57 << std::endl;
  
  // CT10 PDF uncertainty
  double massPlus2    = 0.;
  double massMinus2   = 0.;
  double jesPlus2     = 0.;
  double jesMinus2    = 0.;
  double mass1dPlus2  = 0.;
  double mass1dMinus2 = 0.;
  for (int i = 1; i <= 52; i=i+2) {
    massPlus2  += pow(max(max(ensembles[i].mass-ensembles[0].mass, ensembles[i+1].mass-ensembles[0].mass), 0.), 2);
    massMinus2 += pow(max(max(ensembles[0].mass-ensembles[i].mass, ensembles[0].mass-ensembles[i+1].mass), 0.), 2);
    jesPlus2  += pow(max(max(ensembles[i].jes-ensembles[0].jes, ensembles[i+1].jes-ensembles[0].jes), 0.), 2);
    jesMinus2 += pow(max(max(ensembles[0].jes-ensembles[i].jes, ensembles[0].jes-ensembles[i+1].jes), 0.), 2);
    mass1dPlus2  += pow(max(max(ensembles[i].mass-ensembles[0].mass, ensembles[i+1].mass-ensembles[0].mass), 0.), 2);
    mass1dMinus2 += pow(max(max(ensembles[0].mass-ensembles[i].mass, ensembles[0].mass-ensembles[i+1].mass), 0.), 2);
  }
  // CT10 alpha_s uncertainty
  massPlus2  += pow(ensembles[54].mass-ensembles[0].mass ,2);
  massMinus2 += pow(ensembles[53].mass-ensembles[0].mass ,2);
  jesPlus2  += pow(ensembles[54].jes-ensembles[0].jes ,2);
  jesMinus2 += pow(ensembles[53].jes-ensembles[0].jes ,2);
  mass1dPlus2  += pow(ensembles[54].mass1d-ensembles[0].mass1d ,2);
  mass1dMinus2 += pow(ensembles[53].mass1d-ensembles[0].mass1d ,2);
  
  std::cout << "CT10:  " << ensembles[0].mass << std::endl;
  std::cout << "CT10+: " << sqrt(massPlus2) << std::endl;
  std::cout << "CT10-: " << sqrt(massMinus2) << std::endl;
  double ct10massPlus  = ensembles[0].mass + sqrt(massPlus2);
  double ct10massMinus = ensembles[0].mass - sqrt(massMinus2);
  double ct10jesPlus  = ensembles[0].jes + sqrt(jesPlus2);
  double ct10jesMinus = ensembles[0].jes - sqrt(jesMinus2);
  double ct10mass1dPlus  = ensembles[0].mass1d + sqrt(mass1dPlus2);
  double ct10mass1dMinus = ensembles[0].mass1d - sqrt(mass1dMinus2);
  
  // MSTW PDF uncertainty
  massPlus2  = 0.;
  massMinus2 = 0.;
  jesPlus2     = 0.;
  jesMinus2    = 0.;
  mass1dPlus2  = 0.;
  mass1dMinus2 = 0.;
  for (int i = 56; i <= 95; i=i+2) {
    massPlus2  += pow(max(max(ensembles[i].mass-ensembles[55].mass, ensembles[i+1].mass-ensembles[55].mass), 0.), 2);
    massMinus2 += pow(max(max(ensembles[55].mass-ensembles[i].mass, ensembles[55].mass-ensembles[i+1].mass), 0.), 2);
    jesPlus2  += pow(max(max(ensembles[i].jes-ensembles[55].jes, ensembles[i+1].jes-ensembles[55].jes), 0.), 2);
    jesMinus2 += pow(max(max(ensembles[55].jes-ensembles[i].jes, ensembles[55].jes-ensembles[i+1].jes), 0.), 2);
    mass1dPlus2  += pow(max(max(ensembles[i].mass-ensembles[55].mass, ensembles[i+1].mass-ensembles[55].mass), 0.), 2);
    mass1dMinus2 += pow(max(max(ensembles[55].mass-ensembles[i].mass, ensembles[55].mass-ensembles[i+1].mass), 0.), 2);
  }
  // MSTW alpha_s uncertainty
  massPlus2  += pow(ensembles[97].mass-ensembles[55].mass ,2);
  massMinus2 += pow(ensembles[96].mass-ensembles[55].mass ,2);
  jesPlus2  += pow(ensembles[97].jes-ensembles[55].jes ,2);
  jesMinus2 += pow(ensembles[96].jes-ensembles[55].jes ,2);
  mass1dPlus2  += pow(ensembles[97].mass1d-ensembles[55].mass1d ,2);
  mass1dMinus2 += pow(ensembles[96].mass1d-ensembles[55].mass1d ,2);
  
  std::cout << "MSTW:  " << ensembles[55].mass << std::endl;
  std::cout << "MSTW+: " << sqrt(massPlus2) << std::endl;
  std::cout << "MSTW-: " << sqrt(massMinus2) << std::endl;
  double mstwmassPlus  = ensembles[55].mass + sqrt(massPlus2);
  double mstwmassMinus = ensembles[55].mass - sqrt(massMinus2);
  double mstwjesPlus  = ensembles[55].jes + sqrt(jesPlus2);
  double mstwjesMinus = ensembles[55].jes - sqrt(jesMinus2);
  double mstwmass1dPlus  = ensembles[55].mass1d + sqrt(mass1dPlus2);
  double mstwmass1dMinus = ensembles[55].mass1d - sqrt(mass1dMinus2);
  
  // NNPDF PDF+alpha_s uncertainty
  TH1F* massNNPDF = new TH1F("massNNPDF", "massNNPDF", 500, 172.0, 173.0);
  TH1F* jesNNPDF = new TH1F("jesNNPDF", "jesNNPDF", 100, 0.99, 1.01);
  TH1F* mass1dNNPDF = new TH1F("mass1dNNPDF", "mass1dNNPDF", 100, 172.0, 173.0);
  for (int i = 98; i <= 147; ++i) {
    massNNPDF->Fill(ensembles[i].mass);
    jesNNPDF->Fill(ensembles[i].jes);
    mass1dNNPDF->Fill(ensembles[i].mass1d);
  }
  TF1* gaus = new TF1("gaus", "gaus");
  massNNPDF->Fit("gaus", "LEM");
  std::cout << "NNPDF:    " << gaus->GetParameter(1) << std::endl;
  std::cout << "NNPDF+/-: " << gaus->GetParameter(2) << std::endl;
  double nnpdfmassPlus  = gaus->GetParameter(1) + gaus->GetParameter(2);
  double nnpdfmassMinus = gaus->GetParameter(1) - gaus->GetParameter(2);
  jesNNPDF->Fit("gaus", "LEM");
  double nnpdfjesPlus  = gaus->GetParameter(1) + gaus->GetParameter(2);
  double nnpdfjesMinus = gaus->GetParameter(1) - gaus->GetParameter(2);
  mass1dNNPDF->Fit("gaus", "LEM");
  double nnpdfmass1dPlus  = gaus->GetParameter(1) + gaus->GetParameter(2);
  double nnpdfmass1dMinus = gaus->GetParameter(1) - gaus->GetParameter(2);
  
  
  /*
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

    printf("\t- %4.3f GeV / o %4.3f GeV / + %4.3f GeV \n", ensembles[i].mass1d, ensembles[0].mass1d, ensembles[i+1].mass1d);
    double meanDMConstJES    = (ensembles[i].mass1d - ensembles[i+1].mass1d)/2;
    printf("\tSymmetrized uncertainty: %4.2f GeV \n", meanDMConstJES);
    
    totalmass1dUncertainty2 += meanDMConstJES*meanDMConstJES;    
  }
  
  printf("Systematic uncertainty on top mass: %4.2f GeV \n", sqrt(totalMassUncertainty2));
  printf("Systematic uncertainty on JES     : %4.3f \n", sqrt(totalJESUncertainty2));
  printf("Systematic uncertainty on mTop(1D): %4.2f GeV \n", sqrt(totalmass1dUncertainty2));
  */
  
  printf("Lower envelope mass: %4.3f GeV, upper envelope mass: %4.3f GeV \n", min(min(ct10massMinus,mstwmassMinus),nnpdfmassMinus), max(max(ct10massPlus,mstwmassPlus),nnpdfmassPlus));
  printf("Lower envelope jes: %4.4f GeV, upper envelope jes: %4.4f GeV \n", min(min(ct10jesMinus,mstwjesMinus),nnpdfjesMinus), max(max(ct10jesPlus,mstwjesPlus),nnpdfjesPlus));
  printf("Lower envelope mass1d: %4.3f GeV, upper envelope mass1d: %4.3f GeV \n", min(min(ct10mass1dMinus,mstwmass1dMinus),nnpdfmass1dMinus), max(max(ct10mass1dPlus,mstwmass1dPlus),nnpdfmass1dPlus));
}

