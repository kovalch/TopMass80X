#include "TopMass.h"

TopMass::TopMass(po::variables_map vm) {
  fMethod = vm["method"].as<std::string>();
  fLumi   = vm["lumi"].as<double>();
  
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
  
  if (!vm["task"].as<std::string>().compare("pe")) {
    WriteEnsembleTest(vm, vBinning);
  }
  
  else if (!vm["task"].as<std::string>().compare("sm")) {
    Analysis* analysis = new Analysis(vm, vBinning);
    analysis->Analyze(vm);
  }
}


void TopMass::WriteEnsembleTest(po::variables_map vm, std::vector<float> vBinning) {
  time_t start, end;
  time(&start);
  time(&end);
  
  //if (readCalibration) LoadXML();
  
  int nEnsembles = vm["number"].as<int>();
  int n = 0;
  
  double genMass, mass, massError, massPull, genJES, JES, JESError, JESPull, massAlt, massAltError;
  int bin;
  TFile* ensembleFile;
  TTree* tree;
  
  TFile* tempFile = new TFile("temp.root", "RECREATE");
  ensembleFile = new TFile("ensemble.root", "UPDATE");
  ensembleFile->GetObject("tree", tree);
  
  if (!tree) {
    tree = new TTree("tree", "tree");
    tree->Branch("genMass", &genMass, "genMass/D");
    tree->Branch("mass", &mass, "mass/D");
    tree->Branch("massError", &massError, "massError/D");
    tree->Branch("massPull", &massPull, "massPull/D");
    tree->Branch("genJES", &genJES, "genJES/D");
    tree->Branch("JES", &JES, "JES/D");
    tree->Branch("JESError", &JESError, "JESError/D");
    tree->Branch("JESPull", &JESPull, "JESPull/D");
    tree->Branch("massAlt", &massAlt, "massAlt/D");
    tree->Branch("massAltError", &massAltError, "massAltError/D");
    tree->Branch("bin", &bin, "bin/I");
  }
  else {
    tree->SetBranchAddress("genMass", &genMass);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("massError", &massError);
    tree->SetBranchAddress("massPull", &massPull);
    tree->SetBranchAddress("genJES", &genJES);
    tree->SetBranchAddress("JES", &JES);
    tree->SetBranchAddress("JESError", &JESError);
    tree->SetBranchAddress("JESPull", &JESPull);
    tree->SetBranchAddress("massAlt", &massAlt);
    tree->SetBranchAddress("massAltError", &massAltError);
    tree->SetBranchAddress("bin", &bin);
  }
  
  tempFile->cd();
  Analysis* analysis = new Analysis(vm, vBinning);
  
  while (difftime(end, start) < vm["walltime"].as<int>() * 60 && n < nEnsembles) {
    std::cout << "\n- - - - - - - - - - " << n << " - - - - - - - - - -\n" << std::endl;
    
    tempFile->cd();
    analysis->Analyze(vm);
    
    ensembleFile->cd();
    
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
    
    std::cout << "fill finished" << std::endl;
    
    n++;
    time(&end);
  }
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
