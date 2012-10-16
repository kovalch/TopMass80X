#include <vector>
#include <iostream>
#include <iomanip>

#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom3.h"

void toySystematics() {
  TRandom3* random = new TRandom3(0);
  
  double syst[] = {-1, 1.};
  double stat[] = {1.0, 1.0};
  
  // Initialize histogram
  TH1F* hMean = new TH1F("hMean", "(abs(down) + abs(up))/2.", 100, 0, 5);
  TH1F* hMax  = new TH1F("hMax", "Max(abs(down), abs(up))", 100, 0, 5);
  
  std::cout << std::setiosflags(std::ios::left)
              << std::setw(10) << "true down"
              << std::setw(10) << "true up"
              << std::setw(10) << "toy down"
              << std::setw(10) << "toy up"
              << std::setw(10) << "mean"
              << std::setw(10) << "max"
              << std::endl;
  
  // Loop over all permutations
  for (int i = 0; i < 1000000; i++) {
    
    double down = random->Gaus(syst[0], stat[0]);
    double up   = random->Gaus(syst[1], stat[1]);
    
    double mean = (abs(down) + abs(up))/2.;
    double max  = TMath::Max(abs(down), abs(up));
    
    // Print
    if (i == 0) {
      std::cout << std::setiosflags(std::ios::left)
                << std::setw(10) << syst[0]
                << std::setw(10) << syst[1]
                << std::setw(10) << down
                << std::setw(10) << up
                << std::setw(10) << mean
                << std::setw(10) << max
                << std::endl;
    }
    
    // Fill histogram
    hMean->Fill(mean);
    hMax ->Fill(max);
  }
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 500);
  canvas->Divide(2, 1);
  
  canvas->cd(1);
  hMean->Draw("");
  
  canvas->cd(2);
  hMax->Draw("");
}
