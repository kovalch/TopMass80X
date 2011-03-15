#include "TopMass.h"

void TopMass::Calibrate(TString method) {
 
  std::vector<TH2F*> h2Mass;
  std::vector<TH2F*> h2MassError;
  
  Analysis* a1665 = new Analysis("1665", "root/analyzeTop_1665.root", method, fBins);
  Analysis* a1725 = new Analysis("1725", "root/analyzeTop_1725.root", method, fBins);
  Analysis* a1785 = new Analysis("1785", "root/analyzeTop_1785.root", method, fBins);
  
  analyses.push_back(a1665);
  analyses.push_back(a1725);
  analyses.push_back(a1785);
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);
  
  double genMass[] = {166.5, 172.5, 178.5};
  double genMassError[] = {0.0001, 0.0001, 0.0001};
  double hadTopMass[3];
  double hadTopMassError[3];
  
  for(int i = 0; i < 3; i++){
    analyses.at(i)->Analyze();
    h2Mass.push_back(analyses.at(i)->GetH2Mass());
    h2MassError.push_back(analyses.at(i)->GetH2MassError());
  }
    
  canvasFit->cd();

  TGraphErrors* ghadTopMass;
  
  TF1* linearFit = new TF1("linearFit", "172.5+[0]+(x-172.5)*[1]");
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      if (h2Mass.at(0)->GetCellContent(i+1, j+1) > 0 && h2Mass.at(2)->GetCellContent(i+1, j+1) > 0) {
        for (int k = 0; k < 3; k++) {
          hadTopMass[k] = h2Mass.at(k)->GetCellContent(i+1, j+1);
          hadTopMassError[k] = h2MassError.at(k)->GetCellContent(i+1, j+1);
        }
        ghadTopMass = new TGraphErrors(3, hadTopMass, genMass, hadTopMassError, genMassError);
        ghadTopMass->Draw("A*");
        
        ghadTopMass->Fit("linearFit");
        
        TString path("plot/"); path += method; path += "/"; path += "fit_"; path += i; path += "_"; path += j; path += ".png";
        canvasFit->Print(path);
        
        for (int l = 0; l < 2; l++) {
          fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        }
      }
    }
  }
}

void TopMass::Measure(TString method) {
  Analysis* a = new Analysis("sim", "root/analyzeTop_1725.root", method, 8);
  a->Analyze();
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);
  
  TH2F* h2Mass = a->GetH2Mass();
  TH2F* h2MassError = a->GetH2MassError();

  TH2F* hMassCalibrated = new TH2F();
  hMassCalibrated->SetBins(fBins, 0, 4, fBins, 0, 4);
  hMassCalibrated->SetStats(false);
  hMassCalibrated->SetTitle("Mass");
  hMassCalibrated->SetXTitle("deltaRHadWHadB");
  hMassCalibrated->SetYTitle("deltaRHadQHadQBar");
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      if (fCalibFitParameter[i][j][0] && fCalibFitParameter[i][j][1]) {
        double calib = 172.5 + fCalibFitParameter[i][j][0] + fCalibFitParameter[i][j][1] * (h2Mass->GetCellContent(i+1, j+1) - 172.5);
        double caliberror = sqrt(pow(fCalibFitParError[i][j][0], 2) + pow((h2Mass->GetCellContent(i+1, j+1) - 172.5)*fCalibFitParError[i][j][1], 2) + pow(fCalibFitParameter[i][j][1]*h2MassError->GetCellContent(i+1, j+1), 2));
        
        std::cout << "Measured TopMass: " << calib << " +/- " << caliberror << " GeV" << std::endl;
        
        hMassCalibrated->SetCellContent(i+1, j+1, caliberror);
      }
    }
  }
  
  hMassCalibrated->Draw("COLZ,TEXT");
  hMassCalibrated->SetAxisRange(hMassCalibrated->GetMinimum(0), hMassCalibrated->GetMaximum(), "Z");
  
  canvasFit->Print("test.eps");
}

int main(int argc, char** argv)
{
  TopMass* top = new TopMass(8);
  top->Calibrate("GenMatch");
  top->Measure("GenMatch");
}
