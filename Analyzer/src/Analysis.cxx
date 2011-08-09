#include "Analysis.h"

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
      fIdentifier(identifier), fFile(file), fMethod(method), fBins(bins), fLumi(lumi), fAnalyzed(false)
{
  fChain = new TChain("analyzeHitFit/eventTree");
  fChain->Add(fFile);
  //fChain->Add("root/WJets.root");
  //fChain->Add("root/VQQJets.root");
}

void Analysis::Analyze(bool reanalyze) {

  if (fAnalyzed && !reanalyze) {
    std::cout << "Analysis " << fIdentifier << " has already been done, no reanalyze forced" << std::endl;
    return;
  }

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
  
  if (!fAnalyzed) CreateHistos();
  
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
        cuts += " & sumBPt > 50 & TTBarPt < 50";
        cuts += " & (bProbSSV * hitFitProb) > 0.05";
      }
      
      int entries = fTree->GetEntries(cuts);
      std::cout << entries << std::endl;

      hEntries->SetCellContent(i+1, j+1, entries);
      
      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j);

        hMass     ->SetCellContent(i+1, j+1, fAnalyzer->GetMass());
        hMassError->SetCellContent(i+1, j+1, fAnalyzer->GetMassError());
        hMassSigma->SetCellContent(i+1, j+1, fAnalyzer->GetMassSigma());
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

  fAnalyzed = true;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(fBins);

  hEntries = helper->GetH2("Entries");
  hMass = helper->GetH2("Mass");
  hMassError = helper->GetH2("MassError");
  hMassSigma = helper->GetH2("MassSigma");
  
  hMassCalibrated = helper->GetH2("Mass (Calibrated)");
  hMassErrorCalibrated = helper->GetH2("MassError (Calibrated)");
}

void Analysis::CreateRandomSubset() {
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("target", 1);
  fChain->SetBranchStatus("hadTopMass", 1);
  fChain->SetBranchStatus("hadTopPt", 1);
  fChain->SetBranchStatus("lepTopPt", 1);
  fChain->SetBranchStatus("sumBPt", 1);
  fChain->SetBranchStatus("TTBarPt", 1);
  fChain->SetBranchStatus("hadWRawMass", 1);
  fChain->SetBranchStatus("hitFitChi2", 1);
  fChain->SetBranchStatus("hitFitProb", 1);
  fChain->SetBranchStatus("bProbSSV", 1);
  fChain->SetBranchStatus("event", 1);
  fChain->SetBranchStatus("combi", 1);
  fChain->SetBranchStatus("deltaThetaHadWHadB", 1);
  fChain->SetBranchStatus("deltaThetaHadQHadQBar", 1);
  fChain->SetBranchStatus("PUWeight", 1);

  if (fLumi>0) {
    TRandom3* random = new TRandom3(0);
    double events = 208./35.*fLumi;
    double fullEvents = fChain->GetEntries("combi==0");

    fTree = fChain->CloneTree(0);
    
    int combi;
    fChain->SetBranchAddress("combi", &combi);
    fTree->SetBranchAddress("combi", &combi);
    
    for (int iEntry = 0; iEntry < fChain->GetEntries(); iEntry++) {
      fChain->GetEntry(iEntry);
      if (combi!=0) continue;
      if (random->Rndm() < events/fullEvents) {
        for (int iComb = 0; iComb < 24; iComb++) {
	        fChain->GetEntry(iEntry + iComb);
	        
          if ((iComb != 0 && combi == 0)) {
            iEntry = iEntry + iComb - 1;
            break;
          }
          
          fTree->Fill();
        }
      }
    }
  }
  else {
    fTree = fChain->CloneTree();
  }
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

TH2F* Analysis::GetH2MassCalibrated() {
  return hMassCalibrated;
}

TH2F* Analysis::GetH2MassErrorCalibrated() {
  return hMassErrorCalibrated;
}

TString Analysis::GetIdentifier() {
  return fIdentifier;
}
