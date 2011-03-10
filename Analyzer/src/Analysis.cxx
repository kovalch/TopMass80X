// Analysis.C:
#include "Analysis.h"

void Analysis::Print() const {
  std::cout << "fBins = " << fBins << std::endl;
}

void Analysis::Analyze() {
  std::cout << "Analyze..." << std::endl;
  
  fChain = new TChain("analyzeGenMatch/eventTree");
  fChain->Add("root/analyzeTop_1725.root");
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();
  
  double smearbins = 1.;
  double rangey = 4.;
  double rangexstart = 0.;
  double rangexend = 4.;
  double rangex = rangexend-rangexstart;

  double smearx = smearbins/fBins*rangex;
  double smeary = smearbins/fBins*rangey;
  
  double minEntries = 250;

  TString observable1 = "deltaRHadQHadQBar";
  TString observable2 = "deltaRHadWHadB";
  
  CreateHistos();
  
  MassAnalyzer* analyzer = new GenMatchAnalyzer(fChain);

  for(int i = 0; i < fBins; i++) {
    for(int j = 0; j < fBins; j++) {
      // calculate cuts
      TString cuts1, cuts2;
      std::stringstream stream1, stream2;
      stream1 << rangexstart+(rangex/fBins)*i << "<" << observable1 << "&" << observable1 << "<" << rangexstart+(rangex/fBins)*(i+1) << " & "
             << rangey/fBins*j << "<" << observable2 << "&" << observable2 << "<" << rangey/fBins*(j+1)
             << " & genMatchDr < 1";
      cuts1 = stream1.str();

      stream2 << rangexstart+(rangex/fBins)*i-smearx << "<" << observable1 << "&" << observable1 << "<" << rangexstart+(rangex/fBins)*(i+1)+smearx << " & "
             << rangey/fBins*j-smeary << "<" << observable2 << "&" << observable2 << "<" << rangey/fBins*(j+1)+smeary
             << " & genMatchDr < 1";
      cuts2 = stream2.str();
      
      
      int entries1 = fChain->GetEntries(cuts1);
      int entries2 = fChain->GetEntries(cuts2);

      hEntries->Fill(rangey/fBins*(j+0.5), rangexstart+rangex/fBins*(i+0.5), entries2);
      
      if (entries1 && entries2 > minEntries) {
        analyzer->Analyze(cuts2);
        
        std::cout << analyzer->GetMass() << std::endl;
        
        hMassMean     ->Fill(rangey/fBins*(j+0.5), rangexstart+rangex/fBins*(i+0.5), analyzer->GetMass());
        hMassMeanError->Fill(rangey/fBins*(j+0.5), rangexstart+rangex/fBins*(i+0.5), analyzer->GetMassError());
        hMassSigma    ->Fill(rangey/fBins*(j+0.5), rangexstart+rangex/fBins*(i+0.5), analyzer->GetMassSigma());
      }
    }
  }
  
  canvas->Clear();
  canvas->Divide(2,2);
  
  canvas->cd(1);
  hEntries->Draw("COLZ");

  canvas->cd(2);
  hMassMean->Draw("COLZ");
  hMassMean->SetAxisRange(hMassMean->GetMinimum(0), hMassMean->GetMaximum(), "Z");

  canvas->cd(3);
  hMassMeanError->Draw("COLZ");
  hMassMeanError->SetAxisRange(hMassMeanError->GetMinimum(0), hMassMeanError->GetMaximum(), "Z");
  
  canvas->cd(4);
  hMassSigma->Draw("COLZ");
  hMassSigma->SetAxisRange(hMassSigma->GetMinimum(0), hMassSigma->GetMaximum(), "Z");
  
  canvas->Print("plot/analyzeTop.ps");
}

void Analysis::CreateHistos() {
  double rangey = 4.;
  double rangexstart = 0.;
  double rangexend = 4.;
  TString observable1 = "deltaRHadQHadQBar";
  TString observable2 = "deltaRHadWHadB";

  hEntries = new TH2F();
  hEntries->SetBins(fBins, 0, rangey, fBins, rangexstart, rangexend);
  hEntries->SetStats(false);
  hEntries->SetTitle("Entries");
  hEntries->SetXTitle(observable2);
  hEntries->SetYTitle(observable1);
  
  hMassMean = new TH2F();
  hMassMean->SetBins(fBins, 0, rangey, fBins, rangexstart, rangexend);
  hMassMean->SetStats(false);
  hMassMean->SetTitle("MassMean");
  hMassMean->SetXTitle(observable2);
  hMassMean->SetYTitle(observable1);
  
  hMassMeanError = new TH2F();
  hMassMeanError->SetBins(fBins, 0, rangey, fBins, rangexstart, rangexend);
  hMassMeanError->SetStats(false);
  hMassMeanError->SetTitle("MassMeanError");
  hMassMeanError->SetXTitle(observable2);
  hMassMeanError->SetYTitle(observable1);
  
  hMassSigma = new TH2F();
  hMassSigma->SetBins(fBins, 0, rangey, fBins, rangexstart, rangexend);
  hMassSigma->SetStats(false);
  hMassSigma->SetTitle("MassSigma");
  hMassSigma->SetXTitle(observable2);
  hMassSigma->SetYTitle(observable1);
}

int main(int argc, char** argv)
{
  Analysis* a = new Analysis();
  a->Analyze();
}
