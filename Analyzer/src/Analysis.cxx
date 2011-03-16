#include "Analysis.h"

void Analysis::Analyze(bool reanalyze) {
  if (fAnalyzed) {
    std::cout << "Analysis " << fIdentifier << " has already been done, no reanalyze forced" << std::endl;
    return;
  }

  std::cout << "Analyze " << fIdentifier << " with method " << fMethod << std::endl;
  
  if (!strcmp(fMethod, "GenMatch")) {
    fChain = new TChain("analyzeKinFit/eventTree");
    fAnalyzer = new GenMatchAnalyzer(fIdentifier, fChain);
  }
  else if (!strcmp(fMethod, "MVA")) {
    fChain = new TChain("analyzeMVADisc/eventTree");
    fAnalyzer = new MVAAnalyzer(fIdentifier, fChain);
  }
  else {
    return;
  }
  
  fChain->Add(fFile);
  
  gROOT ->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  gStyle->SetPaintTextFormat(".2f");
  
  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();
  
  double smearBins = 1.;
  double rangeX = 4.;
  double rangeY = 4.;

  double smearX = smearBins/fBins*rangeX;
  double smearY = smearBins/fBins*rangeY;
  
  double minEntries = 500;

  TString observableX = "deltaRHadWHadB";
  TString observableY = "deltaRHadQHadQBar";
  
  CreateHistos();
  
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
        cuts += " & genMatchDr > 0";
      }
      else if (!strcmp(fMethod, "MVA")) {
        cuts += " & mvaDisc > 0";
      }
      
      int entries = fChain->GetEntries(cuts);

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
  hMass->SetAxisRange(hMass->GetMinimum(0), hMass->GetMaximum(), "Z");

  canvas->cd(3);
  hMassError->Draw("COLZ,TEXT");
  hMassError->SetAxisRange(hMassError->GetMinimum(0), hMassError->GetMaximum(), "Z");
  
  canvas->cd(4);
  hMassSigma->Draw("COLZ,TEXT");
  hMassSigma->SetAxisRange(hMassSigma->GetMinimum(0), hMassSigma->GetMaximum(), "Z");
  
  TString path("plot/"); path += fMethod; path += "_"; path += fIdentifier; path += ".eps";
  canvas->Print(path);
  
  fAnalyzed = true;
}

void Analysis::CreateHistos() {
  double rangeY = 4.;
  double rangeXstart = 0.;
  double rangeXend = 4.;
  TString observableX = "deltaRHadWHadB";
  TString observableY = "deltaRHadQHadQBar";

  hEntries = new TH2F();
  hEntries->SetBins(fBins, 0, rangeY, fBins, rangeXstart, rangeXend);
  hEntries->SetStats(false);
  hEntries->SetTitle("Entries");
  hEntries->SetXTitle(observableX);
  hEntries->SetYTitle(observableY);
  
  hMass = new TH2F();
  hMass->SetBins(fBins, 0, rangeY, fBins, rangeXstart, rangeXend);
  hMass->SetStats(false);
  hMass->SetTitle("Mass");
  hMass->SetXTitle(observableX);
  hMass->SetYTitle(observableY);
  
  hMassError = new TH2F();
  hMassError->SetBins(fBins, 0, rangeY, fBins, rangeXstart, rangeXend);
  hMassError->SetStats(false);
  hMassError->SetTitle("MassError");
  hMassError->SetXTitle(observableX);
  hMassError->SetYTitle(observableY);
  
  hMassSigma = new TH2F();
  hMassSigma->SetBins(fBins, 0, rangeY, fBins, rangeXstart, rangeXend);
  hMassSigma->SetStats(false);
  hMassSigma->SetTitle("MassSigma");
  hMassSigma->SetXTitle(observableX);
  hMassSigma->SetYTitle(observableY);
  
  hMassCalibrated = new TH2F();
  hMassCalibrated->SetBins(fBins, 0, rangeY, fBins, rangeXstart, rangeXend);
  hMassCalibrated->SetStats(false);
  hMassCalibrated->SetTitle("Mass (Calibrated)");
  hMassCalibrated->SetXTitle(observableX);
  hMassCalibrated->SetYTitle(observableY);
  
  hMassErrorCalibrated = new TH2F();
  hMassErrorCalibrated->SetBins(fBins, 0, rangeY, fBins, rangeXstart, rangeXend);
  hMassErrorCalibrated->SetStats(false);
  hMassErrorCalibrated->SetTitle("MassError (Calibrated)");
  hMassErrorCalibrated->SetXTitle(observableX);
  hMassErrorCalibrated->SetYTitle(observableY);
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
