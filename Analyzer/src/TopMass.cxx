#include "TopMass.h"

#include <cmath>
//#include <fstream>
#include <map>
#include <time.h>

#include "TFile.h"

#include "Analysis.h"
#include "Helper.h"
#include "ProgramOptionsReader.h"

typedef ProgramOptionsReader po;

TopMass::TopMass() :
  fBinning_(po::GetOption<std::string>("binning")),
  fTask_   (po::GetOption<std::string>("task"))
{
  // set environment variables needed for LHAPDF
  //gSystem->Setenv("LHAPATH", "/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/share/lhapdf/PDFsets");
  
  // Define binning
  
  std::vector<float> vBinning;
  
  if (!strcmp(fBinning_.c_str(),"deltaThetaHadWHadB")) {
    vBinning = {0, float(M_PI)};
  }
  else if (!strcmp(fBinning_.c_str(),"hadTopPt")) {
    vBinning = {0, 50, 75, 100, 125, 150, 200, 400};
  }
  else if (!strcmp(fBinning_.c_str(),"hadTopEta")) {
    vBinning = {-5, -2, -1, -0.3, 0.3, 1, 2, 5};
  }
  else if (!strcmp(fBinning_.c_str(),"hadBPt")) {
    vBinning = {0, 50, 75, 100, 125, 400};
  }
  else if (!strcmp(fBinning_.c_str(),"hadBEta")) {
    vBinning = {-2.5, -1, -0.3, 0.3, 1, 2.5};
  }
  else if (!strcmp(fBinning_.c_str(),"TTBarMass")) {
    vBinning = {200, 400, 450, 500, 550, 700 ,1000};
  }
  else if (!strcmp(fBinning_.c_str(),"TTBarPt")) {
    vBinning = {0, 20, 30, 40, 50, 100, 250};
  }
  else if (!strcmp(fBinning_.c_str(),"deltaRHadQHadQBar")) {
    vBinning = {0.5, 1.25, 1.5, 1.75, 2, 3, 6};
  }
  else if (!strcmp(fBinning_.c_str(),"deltaRHadBLepB")) {
    vBinning = {0.5, 1.25, 2, 2.5, 3, 6};
  }
  else if (!strcmp(fBinning_.c_str(),"topMass")) {
    vBinning = {100, 550};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitTop1[0].M()")) {
    vBinning = {0, 2000};
  }
  else if (!strcmp(fBinning_.c_str(),"bRegTop.fitTop1[0].M()")) {
    vBinning = {100, 550};
  }
  else if (!strcmp(fBinning_.c_str(),"nVertex")) {
    vBinning = {0, 10, 20, 40};
  }
  else{
    std::cerr << "Stopping analysis! Binning " << fBinning_.c_str() << " not defined" <<std::endl;
    exit(0);
  }

  // Start task
  std::cout << "starting task now" << std::endl;

  
  if (!strcmp(fTask_.c_str(),"pe")) {
    WriteEnsembleTest(vBinning);
  }
  
  else if (!strcmp(fTask_.c_str(),"sm")) {
    Analysis* analysis = new Analysis(vBinning);
    analysis->Analyze();
    delete analysis;
  }
}


void TopMass::WriteEnsembleTest(const std::vector<float>& vBinning) {
  time_t start, end;
  time(&start);
  time(&end);
  
  const int nEnsembles = po::GetOption<int>("number");
  int n = 0;
  
  const double genMassRead = po::GetOption<double>("mass");
  const double genJESRead  = po::GetOption<double>("jes" );
  const double genfSigRead = po::GetOption<double>("fsig");

  double genMass, genJES, genfSig;
  int bin;
  std::map<std::string, double> values;

  bool treeCreated = false;
  TTree* tree = 0;
  TFile* ensembleFile = new TFile("ensemble.root", "UPDATE");
  ensembleFile->GetObject("tree", tree);
  
  Analysis* analysis = new Analysis(vBinning);
  
  const double walltime = po::GetOption<double>("walltime");
  while (difftime(end, start) < walltime * 60 && n < nEnsembles) {
    std::cout << "\n- - - - - - - - - - " << n << " - - - - - - - - - -\n" << std::endl;
    
    analysis->Analyze();
    
    const std::map<std::string, TH1F*> histograms = analysis->GetH1s();

    if(!treeCreated){

      for(std::map<std::string, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = -1;
      }

      if (!tree) {
        ensembleFile->cd();
        tree = new TTree("tree", "tree");
        tree->Branch("genMass", &genMass, "genMass/D");
        tree->Branch("genJES" , &genJES , "genJES/D" );
        tree->Branch("genfSig", &genfSig, "genfSig/D");
        tree->Branch("bin", &bin, "bin/I");

        for(std::map<std::string, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          std::string leafName = hist->first+std::string("/D");
          tree->Branch(hist->first.c_str(), &values[hist->first], leafName.c_str());
        }
      }
      else {
        tree->SetBranchAddress("genMass", &genMass);
        tree->SetBranchAddress("genJES" , &genJES );
        tree->SetBranchAddress("genfSig", &genfSig);
        tree->SetBranchAddress("bin", &bin);

        for(std::map<std::string, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          tree->SetBranchAddress(hist->first.c_str(), &values[hist->first]);
        }
      }
      treeCreated = true;
      std::cout << "branching finished" << std::endl;
    }

    genMass = genMassRead;
    genJES  = genJESRead ;
    genfSig = genfSigRead;
    for (unsigned int i = 0; i < vBinning.size()-1; i++) {
      for(std::map<std::string, TH1F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = hist->second->GetBinContent(i+1);
      }
      bin = i+1;
      tree->Fill();
    }
    
    std::cout << "fill finished" << std::endl;
    
    ++n;
    time(&end);
  }
  delete analysis;

  tree->GetCurrentFile()->Write();
  
  Helper* helper = new Helper(fBinning_, vBinning);
  TH1F* hBinning = helper->GetH1("hBinning");
  hBinning->Write();
  
  ensembleFile->Close("R");
}

//bool TopMass::fexists(const char *filename) {
//  ifstream ifile(filename);
//  return ifile;
//}
