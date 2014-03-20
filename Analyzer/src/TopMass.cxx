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
  
  if (!strcmp(fBinning_.c_str(),"top.fitTop1.Pt()")) {
    vBinning = {0, 50, 100, 150, 200, 400};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitTop1.Eta()")) {
    vBinning = {0, 0.5, 1, 2, 4};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitW1.Pt()")) {
    vBinning = {0, 50, 100, 150, 200, 400};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitW1.Eta()")) {
    vBinning = {0, 0.5, 1, 2, 4};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitB1.Pt()")) {
    vBinning = {30, 50, 75, 100, 150, 250};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitB1.Eta()")) {
    vBinning = {0, 0.5, 1, 2.4};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitW1Prod1.Pt()")) {
    vBinning = {30, 50, 75, 100, 150, 250};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitW1Prod1.Eta()")) {
    vBinning = {0, 0.5, 1, 2.4};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitTTBar.M()")) {
    vBinning = {200, 400, 500, 600, 700 ,1000, 1500};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitTTBar.Pt()")) {
    vBinning = {0, 20, 40, 60, 100, 160};
  }
  else if (!strcmp(fBinning_.c_str(),"sqrt(pow(top.fitW1Prod1.Eta()-top.fitW1Prod2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitW1Prod1.Phi()-top.fitW1Prod2.Phi()),2))")) {
    vBinning = {0.5, 1, 2, 3, 6};
  }
  else if (!strcmp(fBinning_.c_str(),"sqrt(pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2))")) {
    vBinning = {0.5, 1, 2, 3, 6};
  }
  else if (!strcmp(fBinning_.c_str(),"top.fitB1.Pt()+top.fitB1.Pt()+top.fitW1Prod1.Pt()+top.fitW1Prod2.Pt()")) {
    vBinning = {100, 200, 300, 400, 500, 600};
  }
  else if (!strcmp(fBinning_.c_str(),"jet.@jet.size()")) {
    vBinning = {3.5, 4.5, 5.5, 6.5};
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
  else if (!strcmp(fBinning_.c_str(),"weight.nVertex")) {
    vBinning = {0, 10, 15, 20, 25, 40};
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
