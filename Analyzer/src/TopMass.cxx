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

#include "Helper.h"
#include "ProgramOptionsReader.h"
#include "XMLConfigReader.h"

typedef ProgramOptionsReader po;
typedef XMLConfigReader xml;

TopMass::TopMass() :
  fMethod_(po::GetOption<std::string>("method")),
  fBins_  (po::GetOption<int        >("bins")),
  fLumi_  (po::GetOption<double     >("lumi"))
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
  
  if (!vm["binning"].as<std::string>().compare("deltaThetaHadWHadB")) {
      float xbins[] = {0, TMath::Pi()};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("hadTopPt")) {
      float xbins[] = {0, 50, 75, 100, 125, 150, 200, 400};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("hadTopEta")) {
      float xbins[] = {-5, -2, -1, -0.3, 0.3, 1, 2, 5};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("hadBPt")) {
      float xbins[] = {0, 50, 75, 100, 125, 400};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("hadBEta")) {
      float xbins[] = {-2.5, -1, -0.3, 0.3, 1, 2.5};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("TTBarMass")) {
      float xbins[] = {200, 400, 450, 500, 550, 700 ,1000};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("TTBarPt")) {
      float xbins[] = {0, 20, 30, 40, 50, 100, 250};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("deltaRHadQHadQBar")) {
      float xbins[] = {0.5, 1.25, 1.5, 1.75, 2, 3, 6};
      vBinning.assign(xbins, xbins + sizeof(xbins) / sizeof(xbins[0]));
  }
  
  else if (!vm["binning"].as<std::string>().compare("deltaRHadBLepB")) {
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
    
    const std::map<TString, TH2F*> histograms = analysis->GetH2s();

    if(!treeCreated){

      for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = -1;
      }

      if (!tree) {
        ensembleFile->cd();
        tree = new TTree("tree", "tree");
        tree->Branch("genMass", &genMass, "genMass/D");
        tree->Branch("genJES" , &genJES , "genJES/D" );
        tree->Branch("genfSig", &genfSig, "genfSig/D");
        tree->Branch("bin", &bin, "bin/I");

        for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          TString leafName = hist->first +TString("/D");
          tree->Branch(hist->first, &values[hist->first], leafName);
        }
      }
      else {
        tree->SetBranchAddress("genMass", &genMass);
        tree->SetBranchAddress("genJES" , &genJES );
        tree->SetBranchAddress("genfSig", &genfSig);
        tree->SetBranchAddress("bin", &bin);

        for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          tree->SetBranchAddress(hist->first, &values[hist->first]);
        }
      }
      treeCreated = true;
      std::cout << "branching finished" << std::endl;
    }
    
    for (int i = 0; i < vBinning.size()-1; i++) {
      // add bin to tree
      genMass   = vm["mass"].as<double>();
      mass      = analysis->GetH1Mass()->GetBinContent(i+1);
      massError = analysis->GetH1MassError()->GetBinContent(i+1);
      massPull  = (mass - genMass)/massError;
      massAlt   = analysis->GetH1MassAlt()->GetBinContent(i+1);
      massAltError = analysis->GetH1MassAltError()->GetBinContent(i+1);
      
      genJES    = vm["jes" ].as<double>();
      JES       = analysis->GetH1JES()->GetBinContent(i+1);
      JESError  = analysis->GetH1JESError()->GetBinContent(i+1);
      JESPull   = (JES - genJES)/JESError;
      
      bin       = i+1;
      
      tree->Fill();
    }

    genMass   = po::GetOption<double>("mass");
    genJES    = po::GetOption<double>("jes" );
    genfSig   = po::GetOption<double>("fsig");
    for (int i = 0; i < vBinning.size()-1; i++) {
      for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = hist->second->GetCellContent(i+1, j+1);
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
  
  Helper* helper = new Helper(vm["binning"].as<std::string>(), vBinning);
  TH1F* hBinning = helper->GetH1("hBinning");
  hBinning->Write();
  
  ensembleFile->Close();
  tempFile = new TFile("temp.root", "RECREATE");
  tempFile->Close();
}

int main(int ac, char** av)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("method,m", po::value<std::string>()->default_value("Ideogram"),
      "Top mass measurement method\n"
      "  GenMatch: \tGaussian fit of correct permutations (MC only)\n"
      "  Ideogram: \tJoint likelihood fit of mt and JES given data sample"
    )
    ("task,t", po::value<std::string>()->default_value("sm"),
      "Task to be done\n"
      "  sm: \tSingle measurement based on input parameters\n"
      "  pe: \tPerform pseudo-experiments, additional settings necessary\n"
      "  diff: \tDifferential top mass, specify binning option"
    )
    ("input,i", po::value<std::string>()->default_value("Summer11_TTJets1725_1.00", "1725_1.00"),
      "Identifier of input file to be analyzed")
    ("lepton,l", po::value<std::string>()->default_value("all"),
      "Lepton+jets channel\n"
      "  electron \t\n"
      "  muon \t\n"
      "  all: \telectron+muon"
    )
    ("binning,b", po::value<std::string>()->default_value("deltaThetaHadWHadB"),
      "Phasespace binning\n"
      "  deltaThetaHadWHadB\n"
      "  hadTopPt\n"
      "  hadBEta"
    )
    ("weight,w", po::value<std::string>()->default_value("muWeight*bWeight*PUWeight"),
      "Event weight used in pseudo-experiments")
    ("fsig,f", po::value<double>()->default_value(1.0), "Signal fraction used in pseudo-experiments")
    ("bdisc,B", po::value<double>()->default_value(0.679), "Threshold for b-jets")
    ("mass,M", po::value<double>()->default_value(172.5), "Input top mass for pseudo-experiments")
    ("jes,J", po::value<double>()->default_value(1.0), "Input JES for pseudo-experiments")
    ("lumi,L", po::value<double>()->default_value(4700.0), "Luminosity for each pseudo-experiment")
    ("number,N", po::value<int>()->default_value(10000), "Number of pseudo-experiments per job")
    ("walltime,W", po::value<int>()->default_value(10), "set walltime limit for pseudo-experiments in minutes")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }
  
  TopMass* top = new TopMass(vm);
}

bool TopMass::fexists(const char *filename) {
  ifstream ifile(filename);
  return ifile;
}
