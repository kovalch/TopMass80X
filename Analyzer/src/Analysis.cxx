// Analysis.C:
#include "Analysis.h"

void Analysis::Analyze() {
  std::cout << "Analyze..." << std::endl;
  
  if (!strcmp(fMethod, "GenMatch")) {
    fChain = new TChain("analyzeGenMatch/eventTree");
    fAnalyzer = new GenMatchAnalyzer(fChain);
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
  double rangeY = 4.;
  double rangeXstart = 0.;
  double rangeXend = 4.;
  double rangeX = rangeXend-rangeXstart;

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
      stream << rangeXstart+(rangeX/fBins)*i-smearX << "<" << observableY << "&"
             << observableY << "<" << rangeXstart+(rangeX/fBins)*(i+1)+smearX << " & "
             << rangeY/fBins*j-smearY << "<" << observableX << "&" 
             << observableX << "<" << rangeY/fBins*(j+1)+smearY
             << " & genMatchDr < 1";
      cuts = stream.str();
      
      int entries = fChain->GetEntries(cuts);

      hEntries->Fill(rangeY/fBins*(j+0.5), rangeXstart+rangeX/fBins*(i+0.5), entries);
      
      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j);
        
        hMass     ->Fill(rangeY/fBins*(j+0.5), rangeXstart+rangeX/fBins*(i+0.5), fAnalyzer->GetMass());
        hMassError->Fill(rangeY/fBins*(j+0.5), rangeXstart+rangeX/fBins*(i+0.5), fAnalyzer->GetMassError());
        hMassSigma->Fill(rangeY/fBins*(j+0.5), rangeXstart+rangeX/fBins*(i+0.5), fAnalyzer->GetMassSigma());
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
  
  canvas->Print("plot/analyzeTop.ps");
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

int main(int argc, char** argv)
{
  Analysis* a = new Analysis("root/analyzeTop_1725.root", "GenMatch", 8);
  a->Analyze();
}
