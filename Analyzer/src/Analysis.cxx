#include "Analysis.h"

Analysis::Analysis(po::variables_map vm) {
  fMethod = vm["method"].as<std::string>();
  fLumi   = vm["lumi"].as<double>();
  fBins   = vm["bins"].as<int>();
  fIdentifier = vm["input"].as<std::string>();

  fFile = "/scratch/hh/lustre/cms/user/mseidel/";
  fFile += fIdentifier;
  fFile += "/analyzeTop.root";
    
  fChain = new TChain("analyzeHitFit/eventTree");
  fChain->Add(fFile);
  //fChain->Add("/scratch/hh/current/cms/user/mseidel/STop.root");
  //fChain->Add("/scratch/hh/current/cms/user/mseidel/Summer11_WJets_1.00_1.00_2b/analyzeTop.root");
  //fChain->Add("root/VQQJets.root");
  
  CreateHistos();
}

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
      fIdentifier(identifier), fMethod(method), fBins(bins), fLumi(lumi)
{
  fFile = "/scratch/hh/lustre/cms/user/mseidel/";
  fFile += fIdentifier;
  fFile += "/analyzeTop.root";
    
  fChain = new TChain("analyzeHitFit/eventTree");
  fChain->Add(fFile);
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
    fAnalyzer = new GenMatchAnalyzer(fIdentifier, fTree);
  }
  else if (!strcmp(fMethod, "MVA")) {
    fAnalyzer = new MVAAnalyzer(fIdentifier, fTree);
  }
  else if (!strcmp(fMethod, "Ideogram")) {
    fAnalyzer = new IdeogramAnalyzer(fIdentifier, fTree);
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
        cuts += " & hadQBSSV<1.74 & hadQBarBSSV<1.74 & hadBBSSV>1.74 & lepBBSSV>1.74";   // b-TAG
        //cuts += " & abs(hadQF)!=5 & abs(hadQBarF)!=5 & abs(hadBF)==5 & abs(lepBF)==5";   // b-TRUTH
        //cuts += " & run < 168000";
        //cuts += " & nlJetPt/(hadQPt+hadQBarPt)>0.3";
        cuts += " & leptonPt > 30";
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
  //fChain->SetBranchStatus("run", 1);
  fChain->SetBranchStatus("hadTopMass", 1);
  /*
	fChain->SetBranchStatus("hadTopPt", 1);
  fChain->SetBranchStatus("lepTopPt", 1);
  fChain->SetBranchStatus("sumBPt", 1);
  fChain->SetBranchStatus("TTBarPt", 1);
  fChain->SetBranchStatus("hadWPt", 1);
  fChain->SetBranchStatus("hadWE", 1);
  fChain->SetBranchStatus("lepWPt", 1);
  fChain->SetBranchStatus("hadBPt", 1);
  fChain->SetBranchStatus("lepBPt", 1);
	*/
  fChain->SetBranchStatus("hadWRawMass", 1);
  //fChain->SetBranchStatus("nlJetPt", 1);
  //fChain->SetBranchStatus("hadQPt", 1);
  //fChain->SetBranchStatus("hadQBarPt", 1);
  fChain->SetBranchStatus("leptonPt", 1);
  //fChain->SetBranchStatus("hitFitChi2", 1);
  fChain->SetBranchStatus("hitFitProb", 1);
  //fChain->SetBranchStatus("bProbSSV", 1);
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
  fChain->SetBranchStatus("MCWeight", 1);
  //fChain->SetBranchStatus("bWeight_bTagSFUp", 1);
  //fChain->SetBranchStatus("bWeight_bTagSFDown", 1);
  //fChain->SetBranchStatus("bWeight_misTagSFUp", 1);
  //fChain->SetBranchStatus("bWeight_misTagSFDown", 1);
  fChain->SetBranchStatus("jetMultiplicity", 1);
  fChain->SetBranchStatus("bottomSSVJetMultiplicity", 1);
  fChain->SetBranchStatus("nVertex", 1);
  fChain->SetBranchStatus("noPtEtaJetPt", 1);
  
  //fChain->SetBranchStatus("pdfWeights", 1);
  
  fChain->SetBranchStatus("hadQBSSV", 1);
  fChain->SetBranchStatus("hadQBarBSSV", 1);
  fChain->SetBranchStatus("hadBBSSV", 1);
  fChain->SetBranchStatus("lepBBSSV", 1);
  
  fChain->SetBranchStatus("hadQF", 1);
  fChain->SetBranchStatus("hadQBarF", 1);
  fChain->SetBranchStatus("hadBF", 1);
  fChain->SetBranchStatus("lepBF", 1);

  if (fLumi>0) {
    TTree* tempTree = fChain->CopyTree("leptonPt>30 & bottomSSVJetMultiplicity > 1 & hitFitProb>0.2");
    //TTree* tempTree = fChain->CopyTree("leptonPt>30"); // CMSSW 4.1.4, no b-jets in selection code
    
    tempFile = new TFile("tempTree.root", "RECREATE");
    
    TRandom3* random = new TRandom3(0);
    
    double maxMCWeight = tempTree->GetMaximum("MCWeight");
    
    if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }
    
    // DATA
    // eventTree->GetEntries("leptonPt>30 & bottomSSVJetMultiplicity > 1 & hitFitProb>0.2 & combi==0")
    // (Long64_t)5144

    int eventsPE = random->Poisson(5144./4700.*fLumi); // add poisson
    int permsMC = tempTree->GetEntries("");
    int eventsDrawn = 0;
    
    fTree = tempTree->CloneTree(0);
    
    int combi;
    tempTree->SetBranchAddress("combi", &combi);
    fTree->SetBranchAddress("combi", &combi);
    
    double MCWeight;
    tempTree->SetBranchAddress("MCWeight", &MCWeight);
    fTree->SetBranchAddress("MCWeight", &MCWeight);
    
    while (eventsDrawn < eventsPE) {
      int drawn = random->Integer(permsMC);
      tempTree->GetEntry(drawn);
      if (combi == 0) {
        if (MCWeight > random->Uniform(0, maxMCWeight)) {
          for (int iComb = 0; iComb < 24; iComb++) {
	          tempTree->GetEntry(drawn + iComb);
	          
            if ((iComb != 0 && combi == 0)) {
              break;
            }
            
            fTree->Fill();
          }
          eventsDrawn++;
          //std::cout << eventsDrawn << std::endl;
        }
      }
    }  
    tempFile->Write();
    tempFile = new TFile("tempTree.root", "RECREATE");
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
