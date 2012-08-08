#include "Analysis.h"

Analysis::Analysis(po::variables_map vm) {
  fMethod = vm["method"].as<std::string>();
  fLumi   = vm["lumi"].as<double>();
  fBins   = vm["bins"].as<int>();
  fWeight = vm["weight"].as<std::string>();
  fSig    = vm["fsig"].as<double>();
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
  
  //fChainB->Add("/scratch/hh/current/cms/user/mseidel/Summer11_sT_muon/analyzeTop.root");
  //fChain->Add("/scratch/hh/current/cms/user/mseidel/Summer11_WJets_1.00_1.00_2b/analyzeTop.root");
  //fChain->Add("root/VQQJets.root");
  
  CreateHistos();
}

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
      fIdentifier(identifier), fMethod(method), fBins(bins), fLumi(lumi)
{
  fFileMuon = "/scratch/hh/lustre/cms/user/mseidel/";
  fFileMuon += fIdentifier;
  fFileMuon += "/analyzeTop.root";
    
  fChain = new TChain("analyzeHitFit/eventTree");
  fChain->Add(fFileMuon);
  //fChain->Add("/scratch/hh/current/cms/user/mseidel/STop.root");
  //fChain->Add("/scratch/hh/current/cms/user/mseidel/Summer11_WJets_1.00_1.00_2b/analyzeTop.root");
  //fChain->Add("root/VQQJets.root");
  
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
   
  Helper* helper = new Helper(fBins);
  helper->SetTDRStyle();
  
  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();
  
  double smearBins = 1.;
  double rangeX = 3.;
  double rangeY = 3.;

  double smearX = smearBins/fBins*rangeX;
  double smearY = smearBins/fBins*rangeY;
  
  double minEntries = 25;
  
  /*if (!strcmp(fMethod, "Ideogram")) {
    minEntries = 1500;
  }*/

  TString observableX = "deltaThetaHadWHadB";
  TString observableY = "deltaThetaHadQHadQBar";
  double bdisc = vm["bdisc"].as<double>();
  
  for(int i = 0; i < fBins; i++) {
    for(int j = 0; j < fBins; j++) {
      // calculate cuts
      TString cuts;
      std::stringstream stream;
      stream << (rangeX/fBins)*(i)-smearX << "<" << observableX << "&"
             << observableX << "<" << (rangeX/fBins)*(i+1)+smearX << " & "
             << rangeY/fBins*(j)-smearY << "<" << observableY << "&" 
             << observableY << "<" << rangeY/fBins*(j+1)+smearY;
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
        cuts += " & hadQBCSV<"; cuts += bdisc;
        cuts += " & hadQBarBCSV<"; cuts += bdisc;
        cuts += " & hadBBCSV>"; cuts += bdisc;
        cuts += " & lepBBCSV>"; cuts += bdisc;   // b-TAG
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

      hEntries->SetCellContent(i+1, j+1, entries);
      
      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j);

        hMass     ->SetCellContent(i+1, j+1, fAnalyzer->GetMass());
        hMassError->SetCellContent(i+1, j+1, fAnalyzer->GetMassError());
        hMassSigma->SetCellContent(i+1, j+1, fAnalyzer->GetMassSigma());
        hMassAlt  ->SetCellContent(i+1, j+1, fAnalyzer->GetMassAlt());
        hMassAltError->SetCellContent(i+1, j+1, fAnalyzer->GetMassAltError());
        hJES      ->SetCellContent(i+1, j+1, fAnalyzer->GetJES());
        hJESError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESError());
        
        std::cout << "Measured mass: " << fAnalyzer->GetMass() << " +/- "
                  << fAnalyzer->GetMassError() << " GeV" << std::endl;
        std::cout << "Measured JES: " << fAnalyzer->GetJES() << " +/- "
                  << fAnalyzer->GetJESError() << " " << std::endl;
        std::cout << "Measured mass (alt): " << fAnalyzer->GetMassAlt() << " +/- "
                  << fAnalyzer->GetMassAltError() << " GeV" << std::endl;
      }
    }
  }
  
  canvas->Clear();
  canvas->Divide(2,2);
  
  canvas->cd(1);
  hEntries->Draw("COLZ");

  canvas->cd(2);
  hMass->Draw("COLZ,TEXT");
  hMass->SetAxisRange(hMass->GetMinimum(0.05), hMass->GetMaximum(), "Z");

  canvas->cd(3);
  hMassError->Draw("COLZ,TEXT");
  hMassError->SetAxisRange(0.05, 5, "Z");
  
  canvas->cd(4);
  hMassSigma->Draw("COLZ,TEXT");
  hMassSigma->SetAxisRange(hMassSigma->GetMinimum(0.05), hMassSigma->GetMaximum(), "Z");
  
  TString path("plot/"); path += fMethod; path += "_"; path += fIdentifier; path += ".eps";
  canvas->Print(path);
  
  delete canvas;
  delete fAnalyzer;
  //delete tempFile;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(fBins);

  hEntries = helper->GetH2("Entries");
  hMass = helper->GetH2("Mass");
  hMassError = helper->GetH2("MassError");
  hMassSigma = helper->GetH2("MassSigma");
  hMassAlt = helper->GetH2("MassAlt");
  hMassAltError = helper->GetH2("MassAltError");
  
  hJES = helper->GetH2("JES");
  hJESError = helper->GetH2("JESError");
  
  hMassCalibrated = helper->GetH2("Mass (Calibrated)");
  hMassErrorCalibrated = helper->GetH2("MassError (Calibrated)");
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
      DrawEvents(fTreeWmu,  eventsPEMuon*(1.-fSig)*3./3.);
      DrawEvents(fTreeSTmu, eventsPEMuon*(1.-fSig)*0./3.);
    }
    
    if (!strcmp(fLepton, "electron") || !strcmp(fLepton, "all")) {
      DrawEvents(fTreeTTe, eventsPEElectron*fSig);
      DrawEvents(fTreeWe,  eventsPEElectron*(1.-fSig)*3./3.);
      DrawEvents(fTreeSTe, eventsPEElectron*(1.-fSig)*0./3.);
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
  
  boost::progress_display progress(nEventsPE, std::cout);
  
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
  
  std::cout << file << ": " << chain->GetEntries("leptonPt>30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2 & combi==0") << " events" << std::endl;
  
  return chain->CopyTree("leptonPt>30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2");
}

TH2F* Analysis::GetH2Mass() {
  return hMass;
}

TH2F* Analysis::GetH2MassError() {
  return hMassError;
}

TH2F* Analysis::GetH2MassAlt() {
  return hMassAlt;
}

TH2F* Analysis::GetH2MassAltError() {
  return hMassAltError;
}

TH2F* Analysis::GetH2MassSigma() {
  return hMassSigma;
}

TH2F* Analysis::GetH2JES() {
  return hJES;
}

TH2F* Analysis::GetH2JESError() {
  return hJESError;
}

TH2F* Analysis::GetH2MassCalibrated() {
  return hMassCalibrated;
}

TH2F* Analysis::GetH2MassErrorCalibrated() {
  return hMassErrorCalibrated;
}

TString Analysis::GetIdentifier() {
  return fIdentifier;
}
