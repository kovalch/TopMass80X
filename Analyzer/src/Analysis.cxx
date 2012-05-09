#include "Analysis.h"

Analysis::Analysis(po::variables_map vm) {
  fMethod = vm["method"].as<std::string>();
  fLumi   = vm["lumi"].as<double>();
  fBins   = vm["bins"].as<int>();
  fWeight = vm["weight"].as<std::string>();
  fSig    = vm["fsig"].as<double>();
  fIdentifier = vm["input"].as<std::string>();
  fLepton = vm["lepton"].as<std::string>();
  
  fFileMuon = "/scratch/hh/lustre/cms/user/mseidel/";
  fFileMuon += fIdentifier;
  fFileMuon += "_muon/analyzeTop.root";
  
  fFileElectron = "/scratch/hh/lustre/cms/user/mseidel/";
  fFileElectron += fIdentifier;
  fFileElectron += "_electron/analyzeTop.root";
  
  fChain = new TChain("analyzeHitFit/eventTree");
  if (!strcmp(fLepton, "muon") || !strcmp(fLepton, "all")) fChain->Add(fFileMuon);
  if (!strcmp(fLepton, "electron") || !strcmp(fLepton, "all")) fChain->Add(fFileElectron);
    
  fChainB = new TChain("analyzeHitFit/eventTree");
  fChainB->Add("/scratch/hh/current/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root");
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
        //cuts += " & jetMultiplicity == 4";
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
        hJES      ->SetCellContent(i+1, j+1, fAnalyzer->GetJES());
        hJESError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESError());
        
        std::cout << "Measured mass: " << fAnalyzer->GetMass() << " +/- "
                  << fAnalyzer->GetMassError() << " GeV" << std::endl;
        std::cout << "Measured JES: " << fAnalyzer->GetJES() << " +/- "
                  << fAnalyzer->GetJESError() << " " << std::endl;
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
  delete tempFile;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(fBins);

  hEntries = helper->GetH2("Entries");
  hMass = helper->GetH2("Mass");
  hMassError = helper->GetH2("MassError");
  hMassSigma = helper->GetH2("MassSigma");
  
  hJES = helper->GetH2("JES");
  hJESError = helper->GetH2("JESError");
  
  hMassCalibrated = helper->GetH2("Mass (Calibrated)");
  hMassErrorCalibrated = helper->GetH2("MassError (Calibrated)");
}

void Analysis::CreateRandomSubset() {
  std::cout << "Create random subset..." << std::endl;
  
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("target", 1);
  fChain->SetBranchStatus("hadTopMass", 1);
  fChain->SetBranchStatus("hadWRawMass", 1);
  fChain->SetBranchStatus("leptonPt", 1);
  fChain->SetBranchStatus("leptonC", 1);
  fChain->SetBranchStatus("leptonId", 1);
  fChain->SetBranchStatus("hitFitProb", 1);
  fChain->SetBranchStatus("run", 1);
  fChain->SetBranchStatus("luminosityBlock", 1);
  fChain->SetBranchStatus("event", 1);
  fChain->SetBranchStatus("combi", 1);
  fChain->SetBranchStatus("deltaThetaHadWHadB", 1);
  fChain->SetBranchStatus("deltaThetaHadQHadQBar", 1);
  fChain->SetBranchStatus("PUWeight", 1);
  fChain->SetBranchStatus("PUWeightUp", 1);
  fChain->SetBranchStatus("PUWeightDown", 1);
  fChain->SetBranchStatus("muWeight", 1);
  fChain->SetBranchStatus("bWeight", 1);
  fChain->SetBranchStatus("bWeight_bTagSFUp", 1);
  fChain->SetBranchStatus("bWeight_bTagSFDown", 1);
  fChain->SetBranchStatus("bWeight_misTagSFUp", 1);
  fChain->SetBranchStatus("bWeight_misTagSFDown", 1);
  fChain->SetBranchStatus("MCWeight", 1);
  fChain->SetBranchStatus("nVertex", 1);
  fChain->SetBranchStatus("hadQBCSV", 1);
  fChain->SetBranchStatus("hadQBarBCSV", 1);
  fChain->SetBranchStatus("hadBBCSV", 1);
  fChain->SetBranchStatus("lepBBCSV", 1);
  
  fChainB->SetBranchStatus("*",0);
  fChainB->SetBranchStatus("target", 1);
  fChainB->SetBranchStatus("hadTopMass", 1);
  fChainB->SetBranchStatus("hadWRawMass", 1);
  fChainB->SetBranchStatus("leptonPt", 1);
  fChainB->SetBranchStatus("leptonC", 1);
  fChainB->SetBranchStatus("leptonId", 1);
  fChainB->SetBranchStatus("hitFitProb", 1);
  fChainB->SetBranchStatus("run", 1);
  fChainB->SetBranchStatus("luminosityBlock", 1);
  fChainB->SetBranchStatus("event", 1);
  fChainB->SetBranchStatus("combi", 1);
  fChainB->SetBranchStatus("deltaThetaHadWHadB", 1);
  fChainB->SetBranchStatus("deltaThetaHadQHadQBar", 1);
  fChainB->SetBranchStatus("PUWeight", 1);
  fChainB->SetBranchStatus("PUWeightUp", 1);
  fChainB->SetBranchStatus("PUWeightDown", 1);
  fChainB->SetBranchStatus("muWeight", 1);
  fChainB->SetBranchStatus("bWeight", 1);
  fChainB->SetBranchStatus("bWeight_bTagSFUp", 1);
  fChainB->SetBranchStatus("bWeight_bTagSFDown", 1);
  fChainB->SetBranchStatus("bWeight_misTagSFUp", 1);
  fChainB->SetBranchStatus("bWeight_misTagSFDown", 1);
  fChainB->SetBranchStatus("MCWeight", 1);
  fChainB->SetBranchStatus("nVertex", 1);
  fChainB->SetBranchStatus("hadQBCSV", 1);
  fChainB->SetBranchStatus("hadQBarBCSV", 1);
  fChainB->SetBranchStatus("hadBBCSV", 1);
  fChainB->SetBranchStatus("lepBBCSV", 1);

  
  if (fLumi>0) {
    TTree* tempTree = fChain->CopyTree("leptonPt>30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2");
    TTree* tempTreeB = fChainB->CopyTree("leptonPt>30 & hadQBCSV<0.679 & hadQBarBCSV<0.679 & hadBBCSV>0.679 & lepBBCSV>0.679 & hitFitProb>0.2");
    //TTree* tempTree = fChain->CopyTree("leptonPt>30"); // CMSSW 4.1.4, no b-jets in selection code
    
    tempFile = new TFile("tempTree.root", "RECREATE");
    
    TRandom3* random = new TRandom3(0);
    
    std::vector<std::string> vWeight;
    boost::split( vWeight, fWeight, boost::is_any_of("-*"));
    
    double maxMCWeight = 1.; // calculate upper bound for combined MCWeight
    double maxMCWeightB = 1.;
    
    for (int i = 0; i < vWeight.size(); ++i) {
      std::cout << vWeight[i] << std::endl;
      maxMCWeight  *= tempTree ->GetMaximum((TString) vWeight[i]);
      maxMCWeightB *= tempTreeB->GetMaximum((TString) vWeight[i]);
    }
    
    std::cout << "maxMCWeightS(" << fWeight << "):" << maxMCWeight  << std::endl;
    std::cout << "maxMCWeightB(" << fWeight << "):" << maxMCWeightB << std::endl;
    
    if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
    if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }
    
    // DATA
    // eventTree->GetEntries("leptonPt>30 & bottomCSVJetMultiplicity > 1 & hitFitProb>0.2 & combi==0")
    // (Long64_t)5144
    double nEventsDataMuon      = 2906.;
    double nEventsDataElectron  = 2268.;
    
    int eventsPEMuon      = random->Poisson(nEventsDataMuon/5000.*fLumi);
    int eventsPEElectron  = random->Poisson(nEventsDataElectron/5000.*fLumi);
    
    if (!strcmp(fLepton, "muon"))     eventsPEElectron  = 0;
    if (!strcmp(fLepton, "electron")) eventsPEMuon      = 0;
    
    int permsMC  = tempTree ->GetEntries("");
    int permsMCB = tempTreeB->GetEntries("");
    
    fTree = tempTree->CloneTree(0);
    
    int target, run, luminosityBlock, event, combi, nVertex, leptonId;
    
    tempTree->SetBranchAddress("target", &target);
    tempTree->SetBranchAddress("run", &run);
    tempTree->SetBranchAddress("luminosityBlock", &luminosityBlock);
    tempTree->SetBranchAddress("event", &event);
    tempTree->SetBranchAddress("combi", &combi);
    tempTree->SetBranchAddress("nVertex", &nVertex);
    tempTree->SetBranchAddress("leptonId", &leptonId);
    
    tempTreeB->SetBranchAddress("target", &target);
    tempTreeB->SetBranchAddress("run", &run);
    tempTreeB->SetBranchAddress("luminosityBlock", &luminosityBlock);
    tempTreeB->SetBranchAddress("event", &event);
    tempTreeB->SetBranchAddress("combi", &combi);
    tempTreeB->SetBranchAddress("nVertex", &nVertex);
    tempTreeB->SetBranchAddress("leptonId", &leptonId);
    
    fTree->SetBranchAddress("target", &target);
    fTree->SetBranchAddress("run", &run);
    fTree->SetBranchAddress("luminosityBlock", &luminosityBlock);
    fTree->SetBranchAddress("event", &event);
    fTree->SetBranchAddress("combi", &combi);
    fTree->SetBranchAddress("nVertex", &nVertex);
    fTree->SetBranchAddress("leptonId", &leptonId);
    
    double hadTopMass, hadWRawMass, leptonPt, leptonC, hitFitProb, deltaThetaHadWHadB, deltaThetaHadQHadQBar, PUWeight, PUWeightUp, PUWeightDown, muWeight, bWeight, bWeight_bTagSFUp, bWeight_bTagSFDown, bWeight_misTagSFUp, bWeight_misTagSFDown, MCWeight, hadQBCSV, hadQBarBCSV, hadBBCSV, lepBBCSV;
    
    tempTree->SetBranchAddress("hadTopMass", &hadTopMass);
    tempTree->SetBranchAddress("hadWRawMass", &hadWRawMass);
    tempTree->SetBranchAddress("leptonPt", &leptonPt);
    tempTree->SetBranchAddress("leptonC", &leptonC);
    tempTree->SetBranchAddress("hitFitProb", &hitFitProb);
    tempTree->SetBranchAddress("deltaThetaHadWHadB", &deltaThetaHadWHadB);
    tempTree->SetBranchAddress("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar);
    tempTree->SetBranchAddress("PUWeight", &PUWeight);
    tempTree->SetBranchAddress("PUWeightUp", &PUWeightUp);
    tempTree->SetBranchAddress("PUWeightDown", &PUWeightDown);
    tempTree->SetBranchAddress("muWeight", &muWeight);
    tempTree->SetBranchAddress("bWeight", &bWeight);
    tempTree->SetBranchAddress("bWeight_bTagSFUp", &bWeight_bTagSFUp);
    tempTree->SetBranchAddress("bWeight_bTagSFDown", &bWeight_bTagSFDown);
    tempTree->SetBranchAddress("bWeight_misTagSFUp", &bWeight_misTagSFUp);
    tempTree->SetBranchAddress("bWeight_misTagSFDown", &bWeight_misTagSFDown);
    tempTree->SetBranchAddress("MCWeight", &MCWeight);
    tempTree->SetBranchAddress("hadQBCSV", &hadQBCSV);
    tempTree->SetBranchAddress("hadQBarBCSV", &hadQBarBCSV);
    tempTree->SetBranchAddress("hadBBCSV", &hadBBCSV);
    tempTree->SetBranchAddress("lepBBCSV", &lepBBCSV);
    
    tempTreeB->SetBranchAddress("hadTopMass", &hadTopMass);
    tempTreeB->SetBranchAddress("hadWRawMass", &hadWRawMass);
    tempTreeB->SetBranchAddress("leptonPt", &leptonPt);
    tempTreeB->SetBranchAddress("leptonC", &leptonC);
    tempTreeB->SetBranchAddress("hitFitProb", &hitFitProb);
    tempTreeB->SetBranchAddress("deltaThetaHadWHadB", &deltaThetaHadWHadB);
    tempTreeB->SetBranchAddress("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar);
    tempTreeB->SetBranchAddress("PUWeight", &PUWeight);
    tempTreeB->SetBranchAddress("PUWeightUp", &PUWeightUp);
    tempTreeB->SetBranchAddress("PUWeightDown", &PUWeightDown);
    tempTreeB->SetBranchAddress("muWeight", &muWeight);
    tempTreeB->SetBranchAddress("bWeight", &bWeight);
    tempTreeB->SetBranchAddress("bWeight_bTagSFUp", &bWeight_bTagSFUp);
    tempTreeB->SetBranchAddress("bWeight_bTagSFDown", &bWeight_bTagSFDown);
    tempTreeB->SetBranchAddress("bWeight_misTagSFUp", &bWeight_misTagSFUp);
    tempTreeB->SetBranchAddress("bWeight_misTagSFDown", &bWeight_misTagSFDown);
    tempTreeB->SetBranchAddress("MCWeight", &MCWeight);
    tempTreeB->SetBranchAddress("hadQBCSV", &hadQBCSV);
    tempTreeB->SetBranchAddress("hadQBarBCSV", &hadQBarBCSV);
    tempTreeB->SetBranchAddress("hadBBCSV", &hadBBCSV);
    tempTreeB->SetBranchAddress("lepBBCSV", &lepBBCSV);
    
    fTree->SetBranchAddress("hadTopMass", &hadTopMass);
    fTree->SetBranchAddress("hadWRawMass", &hadWRawMass);
    fTree->SetBranchAddress("leptonPt", &leptonPt);
    fTree->SetBranchAddress("leptonC", &leptonC);
    fTree->SetBranchAddress("hitFitProb", &hitFitProb);
    fTree->SetBranchAddress("deltaThetaHadWHadB", &deltaThetaHadWHadB);
    fTree->SetBranchAddress("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar);
    fTree->SetBranchAddress("PUWeight", &PUWeight);
    fTree->SetBranchAddress("PUWeightUp", &PUWeightUp);
    fTree->SetBranchAddress("PUWeightDown", &PUWeightDown);
    fTree->SetBranchAddress("muWeight", &muWeight);
    fTree->SetBranchAddress("bWeight", &bWeight);
    fTree->SetBranchAddress("bWeight_bTagSFUp", &bWeight_bTagSFUp);
    fTree->SetBranchAddress("bWeight_bTagSFDown", &bWeight_bTagSFDown);
    fTree->SetBranchAddress("bWeight_misTagSFUp", &bWeight_misTagSFUp);
    fTree->SetBranchAddress("bWeight_misTagSFDown", &bWeight_misTagSFDown);
    fTree->SetBranchAddress("MCWeight", &MCWeight);
    fTree->SetBranchAddress("hadQBCSV", &hadQBCSV);
    fTree->SetBranchAddress("hadQBarBCSV", &hadQBarBCSV);
    fTree->SetBranchAddress("hadBBCSV", &hadBBCSV);
    fTree->SetBranchAddress("lepBBCSV", &lepBBCSV);
    
    std::cout << "while (eventsDrawn < eventsPE)..." << std::endl;
    
    int eventsDrawnMuon     = 0;
    int eventsDrawnElectron = 0;
    
    int nAttempts = 0;
    int nAttemptsS = 0;
    int nAttemptsB = 0;
    
    while (eventsDrawnMuon < eventsPEMuon*fSig || eventsDrawnElectron < eventsPEElectron*fSig) {
      int drawn = random->Integer(permsMC);
      tempTree->GetEntry(drawn);
      ++nAttemptsS;
      ++nAttempts;
      //if (combi != 0) tempTree->GetEntry(drawn-combi);
      if (combi == 0 && ((leptonId==11 && eventsDrawnElectron < eventsPEElectron*fSig) || (leptonId==13 && eventsDrawnMuon < eventsPEMuon*fSig)) ) {
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
        }
        
        if (eventWeight > random->Uniform(0, maxMCWeight)) {
          for (int iComb = 0; iComb < 4; iComb++) {
	          tempTree->GetEntry(drawn + iComb);
	          
            if ((iComb != 0 && combi == 0)) {
              break;
            }
            
            fTree->Fill();
          }
          (leptonId==11) ? ++eventsDrawnElectron : ++eventsDrawnMuon;
        }
      }
    }  
    while (eventsDrawnMuon < eventsPEMuon) {
      int drawn = random->Integer(permsMCB);
      tempTreeB->GetEntry(drawn);
      ++nAttemptsB;
      ++nAttempts;
      //if (combi != 0) tempTree->GetEntry(drawn-combi);
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
        }
        
        if (eventWeight > random->Uniform(0, maxMCWeightB)) {
          for (int iComb = 0; iComb < 4; iComb++) {
	          tempTreeB->GetEntry(drawn + iComb);
	          
            if ((iComb != 0 && combi == 0)) {
              break;
            }
            
            fTree->Fill();
          }
          ++eventsDrawnMuon;
        }
      }
    }  
    std::cout << eventsDrawnMuon << " muon events drawn in " << nAttempts << " attempts." << std::endl;
    std::cout << eventsDrawnElectron << " electron events drawn in " << nAttempts << " attempts." << std::endl;
    std::cout << nAttemptsS << " - " << nAttemptsB << std::endl;
    tempFile->Write();
    //tempFile = new TFile("tempTree.root", "RECREATE");
    delete tempTree;
  }
  else {
    fTree = fChain;
  }
  std::cout << "Created random subset." << std::endl;
  //delete file;
}

TH2F* Analysis::GetH2Mass() {
  return hMass;
}

TH2F* Analysis::GetH2MassError() {
  return hMassError;
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
