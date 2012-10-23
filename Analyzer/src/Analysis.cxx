#include "Analysis.h"

Analysis::Analysis(po::variables_map vm, std::vector<float> v)
  : vBinning(v) {
  fMethod = vm["method"].as<std::string>();
  fLumi   = vm["lumi"].as<double>();
  fBinning = vm["binning"].as<std::string>();
  fWeight = vm["weight"].as<std::string>();
  fSig    = vm["fsig"].as<double>();
  fBDisc  = vm["bdisc"].as<double>();
  fIdentifier = vm["input"].as<std::string>();
  fLepton = vm["lepton"].as<std::string>();
  
  if (!strcmp(fLepton, "muon") || !strcmp(fLepton, "all")) {
    fFileMuon = "/scratch/hh/lustre/cms/user/mseidel/";
    fFileMuon += fIdentifier;
    fFileMuon += "_muon/analyzeTop.root";
    
    fTreeTTmu = PrepareTree(fFileMuon);
    fTreeWmu  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root");
    fTreeSTmu = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root");
  }
  
  if (!strcmp(fLepton, "electron") || !strcmp(fLepton, "all")) {
    fFileElectron = "/scratch/hh/lustre/cms/user/mseidel/";
    fFileElectron += fIdentifier;
    fFileElectron += "_electron/analyzeTop.root";
    
    fTreeTTe  = PrepareTree(fFileElectron);
    fTreeWe   = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root");
    fTreeSTe  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root");
  }
  
  CreateHistos();
}

void Analysis::Analyze(po::variables_map vm) {

  std::cout << "Analyze " << fIdentifier << " with method " << fMethod << std::endl;
  
  CreateRandomSubset();
  
  MassAnalyzer* fAnalyzer;
  
  if (!strcmp(fMethod, "GenMatch")) {
    fAnalyzer = new GenMatchAnalyzer(fIdentifier + "_" + fLepton, fTree);
  }
  else if (!strcmp(fMethod, "MVA")) {
    fAnalyzer = new MVAAnalyzer(fIdentifier + "_" + fLepton, fTree);
  }
  else if (!strcmp(fMethod, "Ideogram")) {
    fAnalyzer = new IdeogramAnalyzer(fIdentifier + "_" + fLepton, fTree);
  }
  else {
    return;
  }
   
  Helper* helper = new Helper(fBinning, vBinning);
  helper->SetTDRStyle();
  
  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();
  
  for(int i = 0; i < vBinning.size()-1; i++) {
    // calculate cuts
    TString cuts;
    std::stringstream stream;
    stream << vBinning[i] << " < " << fBinning << " & "
           << fBinning << " < " << vBinning[i+1];
    cuts = stream.str();
    
    if (!strcmp(fMethod, "GenMatch")) {
      cuts += " & target == 1";
    }
    else if (!strcmp(fMethod, "MVA")) {
      cuts += " & mvaDisc > 0";
    }
    else if (!strcmp(fMethod, "Ideogram")) {
      //cuts += " & target == 1";
      //cuts += " & (target == 0 | target == 1)";
      //cuts += " & sumBPt > 50";
      //cuts += " & TTBarPt < 50";
      //cuts += " & jetMultiplicity != 4";
      cuts += " & hitFitProb > 0.2";
      //cuts += " & bProbSSV > 0.3";
      //cuts += " & bProbSSV*hitFitProb > 0.05";
      //cuts += " & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hadBBSSV>1.74 & lepBBSSV>1.74";   // b-TAG
      cuts += " & hadQBCSV<"; cuts += fBDisc;
      cuts += " & hadQBarBCSV<"; cuts += fBDisc;
      cuts += " & hadBBCSV>"; cuts += fBDisc;
      cuts += " & lepBBCSV>"; cuts += fBDisc;   // b-TAG
      //cuts += " & abs(hadQF)!=5 & abs(hadQBarF)!=5 & abs(hadBF)==5 & abs(lepBF)==5";   // b-TRUTH
      //cuts += " & hadQPt>33 & hadBPt>33 & lepBPt>33";
      //cuts += " & run < 168000";
      //cuts += " & nlJetPt/(hadQPt+hadQBarPt)>0.3";
      cuts += " & leptonPt > 30";
      //cuts += " & jetsPt[2] > 50";
      //cuts += " & leptonC == 1";
      //cuts += " & nVertex >= 5";
      //cuts += " & noPtEtaJetPt < 50";
    }
    
    int entries = fTree->GetEntries(cuts);
    std::cout << cuts << std::endl;
    std::cout << entries << std::endl;

    hEntries->SetBinContent(i+1, entries);
    
    fAnalyzer->Analyze(cuts, i, 0);

    hMass         ->SetBinContent (i+1, fAnalyzer->GetMass());
    hMass         ->SetBinError   (i+1, fAnalyzer->GetMassError());
    hMassError    ->SetBinContent (i+1, fAnalyzer->GetMassError());
    hMassSigma    ->SetBinContent (i+1, fAnalyzer->GetMassSigma());
    hMassAlt      ->SetBinContent (i+1, fAnalyzer->GetMassAlt());
    hMassAlt      ->SetBinError   (i+1, fAnalyzer->GetMassAltError());
    hMassAltError ->SetBinContent (i+1, fAnalyzer->GetMassAltError());
    hJES          ->SetBinContent (i+1, fAnalyzer->GetJES());
    hJES          ->SetBinError   (i+1, fAnalyzer->GetJESError());
    hJESError     ->SetBinContent (i+1, fAnalyzer->GetJESError());
    
    std::cout << "Measured mass: " << fAnalyzer->GetMass() << " +/- "
              << fAnalyzer->GetMassError() << " GeV" << std::endl;
    std::cout << "Measured JES: " << fAnalyzer->GetJES() << " +/- "
              << fAnalyzer->GetJESError() << " " << std::endl;
    std::cout << "Measured mass (alt): " << fAnalyzer->GetMassAlt() << " +/- "
              << fAnalyzer->GetMassAltError() << " GeV" << std::endl;
  }
  
  canvas->Clear();
  canvas->Divide(2,2);
  
  canvas->cd(1);
  hEntries->Draw("E1");

  canvas->cd(2);
  hMass->Draw("E1");
  hMass->Fit("pol0");

  canvas->cd(3);
  hMassAlt->Draw("E1");
  hMassAlt->Fit("pol0");
  
  canvas->cd(4);
  hJES->Draw("E1");
  hJES->Fit("pol0");
  
  TString path("plot/"); path += fMethod; path += "_"; path += fIdentifier; path += "_"; path += fBinning; path += ".eps";
  canvas->Print(path);
  
  TString pathr("plot/"); pathr += fMethod; pathr += "_"; pathr += fIdentifier; pathr += "_"; pathr += fBinning; pathr += ".root";
  canvas->Print(pathr);
  
  delete canvas;
  delete fAnalyzer;
  //delete tempFile;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(fBinning, vBinning);

  hEntries = helper->GetH1("Entries");
  hMass = helper->GetH1("Mass");
  hMassError = helper->GetH1("MassError");
  hMassSigma = helper->GetH1("MassSigma");
  hMassAlt = helper->GetH1("MassAlt");
  hMassAltError = helper->GetH1("MassAltError");
  
  hJES = helper->GetH1("JES");
  hJESError = helper->GetH1("JESError");
  
  hMassCalibrated = helper->GetH1("Mass (Calibrated)");
  hMassErrorCalibrated = helper->GetH1("MassError (Calibrated)");
}

void Analysis::CreateRandomSubset() {  
  if (fLumi>0) {
    std::cout << "Create random subset..." << std::endl;
    
    time_t start, end;
    time(&start);
    time(&end);
  
    fTree = new TTree("fTree", "fTree");
  
    fTree->Branch("target", &target, "target/I");
    fTree->Branch("run", &run, "run/I");
    fTree->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/I");
    fTree->Branch("event", &event, "event/I");
    fTree->Branch("combi", &combi, "combi/I");
    fTree->Branch("nVertex", &nVertex, "nVertex/I");
    fTree->Branch("leptonId", &leptonId, "leptonId/I");
    
    fTree->Branch("hadTopMass", &hadTopMass, "hadTopMass/D");
    fTree->Branch("hadWRawMass", &hadWRawMass, "hadWRawMass/D");
    fTree->Branch("leptonPt", &leptonPt, "leptonPt/D");
    fTree->Branch("leptonC", &leptonC, "leptonC/D");
    fTree->Branch("hitFitProb", &hitFitProb, "hitFitProb/D");
    fTree->Branch("deltaThetaHadWHadB", &deltaThetaHadWHadB, "deltaThetaHadWHadB/D");
    fTree->Branch("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar, "deltaThetaHadQHadQBar/D");
    fTree->Branch("PUWeight", &PUWeight, "PUWeight/D");
    fTree->Branch("PUWeightUp", &PUWeightUp, "PUWeightUp/D");
    fTree->Branch("PUWeightDown", &PUWeightDown, "PUWeightDown/D");
    fTree->Branch("muWeight", &muWeight, "muWeight/D");
    fTree->Branch("bWeight", &bWeight, "bWeight/D");
    fTree->Branch("bWeight_bTagSFUp", &bWeight_bTagSFUp, "bWeight_bTagSFUp/D");
    fTree->Branch("bWeight_bTagSFDown", &bWeight_bTagSFDown, "bWeight_bTagSFDown/D");
    fTree->Branch("bWeight_misTagSFUp", &bWeight_misTagSFUp, "bWeight_misTagSFUp/D");
    fTree->Branch("bWeight_misTagSFDown", &bWeight_misTagSFDown, "bWeight_misTagSFDown/D");
    fTree->Branch("MCWeight", &MCWeight, "MCWeight/D");
    fTree->Branch("mcWeight", &mcWeight, "mcWeight/D");
    fTree->Branch("hadQBCSV", &hadQBCSV, "hadQBCSV/D");
    fTree->Branch("hadQBarBCSV", &hadQBarBCSV, "hadQBarBCSV/D");
    fTree->Branch("hadBBCSV", &hadBBCSV, "hadBBCSV/D");
    fTree->Branch("lepBBCSV", &lepBBCSV, "lepBBCSV/D");
    fTree->Branch("lepBBCSV", &lepBBCSV, "lepBBCSV/D");
    fTree->Branch("jetsPt", &jetsPt, "jetsPt[4]/D");
    fTree->Branch("pdfWeights", &pdfWeights, "pdfWeights[44]/D");  
    
    TRandom3* random = new TRandom3(0);
    
    // DATA
    // eventTree->GetEntries("leptonPt>30 & bottomCSVJetMultiplicity > 1 & hitFitProb>0.2 & combi==0")
    // (Long64_t)5144
    double nEventsDataMuon      = 2906.;
    double nEventsDataElectron  = 2268.;
    
    int eventsPEMuon      = random->Poisson(nEventsDataMuon/5000.*fLumi);
    int eventsPEElectron  = random->Poisson(nEventsDataElectron/5000.*fLumi);
    
    fFileMuon = "/scratch/hh/lustre/cms/user/mseidel/";
    fFileMuon += fIdentifier;
    fFileMuon += "_muon/analyzeTop.root";
    
    if (!strcmp(fLepton, "muon") || !strcmp(fLepton, "all")) {
      DrawEvents(fTreeTTmu, eventsPEMuon*fSig);
      DrawEvents(fTreeWmu,  eventsPEMuon*(1.-fSig)*1./4.);
      DrawEvents(fTreeSTmu, eventsPEMuon*(1.-fSig)*3./4.);
    }
    
    if (!strcmp(fLepton, "electron") || !strcmp(fLepton, "all")) {
      DrawEvents(fTreeTTe, eventsPEElectron*fSig);
      DrawEvents(fTreeWe,  eventsPEElectron*(1.-fSig)*1./4.);
      DrawEvents(fTreeSTe, eventsPEElectron*(1.-fSig)*3./4.);
    }
    
    time(&end);
    std::cout << "Created random subset in " << difftime(end, start) << " seconds." << std::endl;
  }
  else {
    fFileMuon = "/scratch/hh/lustre/cms/user/mseidel/";
    fFileMuon += fIdentifier;
    fFileMuon += "_muon/analyzeTop.root";
  
    fFileElectron = "/scratch/hh/lustre/cms/user/mseidel/";
    fFileElectron += fIdentifier;
    fFileElectron += "_electron/analyzeTop.root";
    
    fChain = new TChain("analyzeHitFit/eventTree");
    if (!strcmp(fLepton, "muon") || !strcmp(fLepton, "all")) fChain->Add(fFileMuon);
    if (!strcmp(fLepton, "electron") || !strcmp(fLepton, "all")) fChain->Add(fFileElectron);
    fTree = fChain;
  }
}
  
void Analysis::DrawEvents(TTree* tempTree, double nEventsPE) {
  std::cout << "nEventsPE: " << nEventsPE << std::endl;
  
  mcWeight = 1;
  
  fTree->CopyAddresses(tempTree);
  
  TRandom3* random = new TRandom3(0);
  
  int permsMC  = tempTree->GetEntries("");
  
  std::vector<std::string> vWeight;
  boost::split( vWeight, fWeight, boost::is_any_of("-*"));
  
  double maxMCWeight = 1.; // calculate upper bound for combined MCWeight
  double maxMCWeightB = 1.;
  
  for (int i = 0; i < vWeight.size(); ++i) {
    std::cout << vWeight[i] << std::endl;
    if (strncmp((TString) vWeight[i], "pdfWeights", 10)) maxMCWeight *= tempTree->GetMaximum((TString) vWeight[i]);
    else maxMCWeight *= tempTree->GetMaximum("pdfWeights");
  }
  
  std::cout << "maxMCWeight(" << fWeight << "):" << maxMCWeight  << std::endl;
  
  if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
  if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }
  
  std::cout << "while (eventsDrawn < eventsPE)..." << std::endl;
  
  int eventsDrawn = 0;
  int nAttempts = 0;
  
  boost::progress_display progress((int)nEventsPE, std::cout);
  
  while (eventsDrawn < (int)nEventsPE) {
    int drawn = random->Integer(permsMC);
    tempTree->GetEntry(drawn);
    ++nAttempts;
    //if (combi != 0) tempTree->GetEntry(drawn-combi); // Biases towards events with more combis
    if (combi == 0) {
      double eventWeight = 1.;
      for (int i = 0; i < vWeight.size(); ++i) {
        if (!strcmp((TString) vWeight[i], "MCWeight")) eventWeight *= MCWeight;
        else if (!strcmp((TString) vWeight[i], "muWeight")) eventWeight *= muWeight;
        else if (!strcmp((TString) vWeight[i], "PUWeight")) eventWeight *= PUWeight;
        else if (!strcmp((TString) vWeight[i], "PUWeightUp")) eventWeight *= PUWeightUp;
        else if (!strcmp((TString) vWeight[i], "PUWeightDown")) eventWeight *= PUWeightDown;
        else if (!strcmp((TString) vWeight[i], "bWeight")) eventWeight *= bWeight;
        else if (!strcmp((TString) vWeight[i], "bWeight_bTagSFUp")) eventWeight *= bWeight_bTagSFUp;
        else if (!strcmp((TString) vWeight[i], "bWeight_bTagSFDown")) eventWeight *= bWeight_bTagSFDown;
        else if (!strcmp((TString) vWeight[i], "bWeight_misTagSFUp")) eventWeight *= bWeight_misTagSFUp;
        else if (!strcmp((TString) vWeight[i], "bWeight_misTagSFDown")) eventWeight *= bWeight_misTagSFDown;
        else if (!strncmp((TString) vWeight[i], "pdfWeights", 10)) {
          std::string sub = vWeight[i].substr(11);
          int pdfWeightN = atol(sub.c_str());
          //std::cout << "pdfWeightN: " << pdfWeightN << std::endl;
          eventWeight *= pdfWeights[pdfWeightN];
        }
      }
      
      if (eventWeight > random->Uniform(0, maxMCWeight)) {
        for (int iComb = 0; iComb < 4; iComb++) {
          tempTree->GetEntry(drawn + iComb);
          
          if ((iComb != 0 && combi == 0)) {
            break;
          }
          
          fTree->Fill();
        }
        
        //std::cout << "mcWeight: " << mcWeight << std::endl;
        
        if (mcWeight < 0) {
          eventsDrawn += -1;
          //progress    += -1;
        }
        else {
          ++eventsDrawn;
          ++progress;
        }
      }
    }
  }
  std::cout << eventsDrawn << " events drawn in " << nAttempts << " attempts." << std::endl;
  //tempFile->Write();
  //tempFile = new TFile("tempTree.root", "RECREATE");
  //delete tempTree;
}

TTree* Analysis::PrepareTree(TString file) {
  TChain* chain = new TChain("analyzeHitFit/eventTree");
  chain->Add(file);
  
  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("target", 1);
  chain->SetBranchStatus("hadTopMass", 1);
  chain->SetBranchStatus("hadWRawMass", 1);
  chain->SetBranchStatus("leptonPt", 1);
  chain->SetBranchStatus("leptonC", 1);
  chain->SetBranchStatus("leptonId", 1);
  chain->SetBranchStatus("hitFitProb", 1);
  chain->SetBranchStatus("run", 1);
  chain->SetBranchStatus("luminosityBlock", 1);
  chain->SetBranchStatus("event", 1);
  chain->SetBranchStatus("combi", 1);
  chain->SetBranchStatus("deltaThetaHadWHadB", 1);
  chain->SetBranchStatus("deltaThetaHadQHadQBar", 1);
  chain->SetBranchStatus("PUWeight", 1);
  chain->SetBranchStatus("PUWeightUp", 1);
  chain->SetBranchStatus("PUWeightDown", 1);
  chain->SetBranchStatus("muWeight", 1);
  chain->SetBranchStatus("bWeight", 1);
  chain->SetBranchStatus("bWeight_bTagSFUp", 1);
  chain->SetBranchStatus("bWeight_bTagSFDown", 1);
  chain->SetBranchStatus("bWeight_misTagSFUp", 1);
  chain->SetBranchStatus("bWeight_misTagSFDown", 1);
  chain->SetBranchStatus("MCWeight", 1);
  chain->SetBranchStatus("mcWeight", 1);
  chain->SetBranchStatus("nVertex", 1);
  chain->SetBranchStatus("hadQBCSV", 1);
  chain->SetBranchStatus("hadQBarBCSV", 1);
  chain->SetBranchStatus("hadBBCSV", 1);
  chain->SetBranchStatus("lepBBCSV", 1);
  chain->SetBranchStatus("jetsPt", 1);
  chain->SetBranchStatus("pdfWeights", 1);
  
  TString selection("leptonPt>30 & hitFitProb>0.2");
  selection += " & hadQBCSV<";    selection += fBDisc;
  selection += " & hadQBarBCSV<"; selection += fBDisc;
  selection += " & hadBBCSV>";    selection += fBDisc;
  selection += " & lepBBCSV>";    selection += fBDisc;
  
  std::cout << file << ": " << chain->GetEntries(selection + " & combi==0") << " events" << std::endl;
  
  return chain->CopyTree(selection);
}

TH1F* Analysis::GetH1Mass() {
  return hMass;
}

TH1F* Analysis::GetH1MassError() {
  return hMassError;
}

TH1F* Analysis::GetH1MassAlt() {
  return hMassAlt;
}

TH1F* Analysis::GetH1MassAltError() {
  return hMassAltError;
}

TH1F* Analysis::GetH1MassSigma() {
  return hMassSigma;
}

TH1F* Analysis::GetH1JES() {
  return hJES;
}

TH1F* Analysis::GetH1JESError() {
  return hJESError;
}

TH1F* Analysis::GetH1MassCalibrated() {
  return hMassCalibrated;
}

TH1F* Analysis::GetH1MassErrorCalibrated() {
  return hMassErrorCalibrated;
}

TString Analysis::GetIdentifier() {
  return fIdentifier;
}
