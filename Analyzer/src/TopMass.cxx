#include "TopMass.h"

#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <time.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"

#include "Helper.h"
#include "ProgramOptionsReader.h"
#include "XMLConfigReader.h"

typedef ProgramOptionsReader po;
typedef XMLConfigReader xml;

TopMass::TopMass() :
  fMethod_ (po::GetOption<std::string>("method")),
  fBinning_(po::GetOption<std::string>("binning")),
  fLumi_   (po::GetOption<double     >("lumi"))
{
  // check existence of a temp directory and create one if not available
  TString tempDir(gSystem->Getenv("TMPDIR"));
  if(tempDir.IsNull() || !tempDir.Length()){
    tempDir = gSystem->GetFromPipe("mktemp -d");
    gSystem->Setenv("TMPDIR", tempDir);
  }
  std::cout << "Directory to be used for temporary files: " << tempDir << std::endl;

  // set environment variables needed for LHAPDF
  gSystem->Setenv("LHAPATH", "/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/share/lhapdf/PDFsets");
  
  // Define binning
  
  std::vector<float> vBinning;
  
  if (!po::GetOption<std::string>("binning").compare("deltaThetaHadWHadB")) {
      float xbins[] = {0, TMath::Pi()};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("hadTopPt")) {
      float xbins[] = {0, 50, 75, 100, 125, 150, 200, 400};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("hadTopEta")) {
      float xbins[] = {-5, -2, -1, -0.3, 0.3, 1, 2, 5};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("hadBPt")) {
      float xbins[] = {0, 50, 75, 100, 125, 400};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("hadBEta")) {
      float xbins[] = {-2.5, -1, -0.3, 0.3, 1, 2.5};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("TTBarMass")) {
      float xbins[] = {200, 400, 450, 500, 550, 700 ,1000};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("TTBarPt")) {
      float xbins[] = {0, 20, 30, 40, 50, 100, 250};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("deltaRHadQHadQBar")) {
      float xbins[] = {0.5, 1.25, 1.5, 1.75, 2, 3, 6};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!po::GetOption<std::string>("binning").compare("deltaRHadBLepB")) {
      float xbins[] = {0.5, 1.25, 2, 2.5, 3, 6};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  // Start task
  
  if (!po::GetOption<std::string>("task").compare("pe")) {
    WriteEnsembleTest(vBinning);
  }
  
  else if (!po::GetOption<std::string>("task").compare("sm")) {
    Analysis* analysis = new Analysis(vBinning);
    analysis->Analyze();
  }
}


void TopMass::WriteEnsembleTest(std::vector<float> vBinning) {
  time_t start, end;
  time(&start);
  time(&end);
  
  //if (readCalibration) LoadXML();
  
  int nEnsembles = po::GetOption<int>("number");
  int n = 0;
  
  double genMass, genJES, genfSig;
  int bin;
  std::map<TString, double> values;

  bool treeCreated = false;
  TTree* tree = 0;
  TFile* ensembleFile = new TFile("ensemble.root", "UPDATE");
  ensembleFile->GetObject("tree", tree);
  
  Analysis* analysis = new Analysis(vBinning);
  
  while (difftime(end, start) < po::GetOption<double>("walltime") * 60 && n < nEnsembles) {
    std::cout << "\n- - - - - - - - - - " << n << " - - - - - - - - - -\n" << std::endl;
    
    analysis->Analyze();
    
    const std::map<TString, TH1F*> histograms = analysis->GetH1s();

    if(!treeCreated){

      for(std::map<TString, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = -1;
      }

      if (!tree) {
        ensembleFile->cd();
        tree = new TTree("tree", "tree");
        tree->Branch("genMass", &genMass, "genMass/D");
        tree->Branch("genJES" , &genJES , "genJES/D" );
        tree->Branch("genfSig", &genfSig, "genfSig/D");
        tree->Branch("bin", &bin, "bin/I");

        for(std::map<TString, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          TString leafName = hist->first +TString("/D");
          tree->Branch(hist->first, &values[hist->first], leafName);
        }
      }
      else {
        tree->SetBranchAddress("genMass", &genMass);
        tree->SetBranchAddress("genJES" , &genJES );
        tree->SetBranchAddress("genfSig", &genfSig);
        tree->SetBranchAddress("bin", &bin);

        for(std::map<TString, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          tree->SetBranchAddress(hist->first, &values[hist->first]);
        }
      }
      treeCreated = true;
      std::cout << "branching finished" << std::endl;
    }

    genMass   = po::GetOption<double>("mass");
    genJES    = po::GetOption<double>("jes" );
    genfSig   = po::GetOption<double>("fsig");
    for (unsigned int i = 0; i < vBinning.size()-1; i++) {
      for(std::map<TString, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = hist->second->GetBinContent(i+1);
      }

      bin       = i+1;

      tree->Fill();
    }
    
    std::cout << "fill finished" << std::endl;
    
    ++n;
    time(&end);
  }
  delete analysis;

  tree->GetCurrentFile()->Write();
  
  Helper* helper = new Helper(po::GetOption<std::string>("binning"), vBinning);
  TH1F* hBinning = helper->GetH1("hBinning");
  hBinning->Write();
  
  ensembleFile->Close();
}

bool TopMass::fexists(const char *filename) {
  ifstream ifile(filename);
  return ifile;
}
