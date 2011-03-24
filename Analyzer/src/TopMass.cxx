#include "TopMass.h"

TopMass::TopMass(TString method, int bins, double lumi) : fMethod(method), fBins(bins), fLumi(lumi) {


/*  a1665_jes_up = new Analysis("1665_jes_up", "root/analyzeTop_1665_jes_up.root", fMethod, fBins, fLumi);
  a1725_jes_up = new Analysis("1725_jes_up", "root/analyzeTop_1725_jes_up.root", fMethod, fBins, fLumi);
  a1785_jes_up = new Analysis("1785_jes_up", "root/analyzeTop_1785_jes_up.root", fMethod, fBins, fLumi);
  
  a1665_jes_down = new Analysis("1665_jes_down", "root/analyzeTop_1665_jes_down.root", fMethod, fBins, fLumi);
  a1725_jes_down = new Analysis("1725_jes_down", "root/analyzeTop_1725_jes_down.root", fMethod, fBins, fLumi);
  a1785_jes_down = new Analysis("1785_jes_down", "root/analyzeTop_1785_jes_down.root", fMethod, fBins, fLumi);*/
  
  //aSim = new Analysis("sim", "root/analyzeTop_1725.root", fMethod, fBins);
  
  //Calibrate();
  //Measure(aSim);
  //Systematics();
  
  WriteEnsembleTestTree();
  //EvalEnsembleTest();
}


void TopMass::WriteEnsembleTestTree() {
  massPoint m1665(166.5, "1665");
  massPoint m1725(172.5, "1725");
  massPoint m1785(178.5, "1785");
   
  massPoints.push_back(m1665);
  massPoints.push_back(m1725);
  massPoints.push_back(m1785);
  
  int nEnsembles = 100;
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->analysis = new Analysis(iMassPoint->identifier, iMassPoint->fileName, fMethod, fBins, fLumi);
    iMassPoint->h3Mass = new TH3F("h3Mass_" + iMassPoint->identifier, "h3Mass_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, 150, 200);
    iMassPoint->h3MassPull = new TH3F("h3MassPull_" + iMassPoint->identifier, "h3MassPull_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, -5, 5);
    for (int n = 0; n < nEnsembles; n++) {
      iMassPoint->analysis->Analyze(true);
      for (int i = 0; i < fBins; i++) {
        for (int j = 0; j < fBins; j++) {
          iMassPoint->h3Mass->Fill(3./fBins*i, 3./fBins*j, iMassPoint->analysis->GetH2Mass()->GetCellContent(i+1, j+1));
          iMassPoint->h3MassPull->Fill(3./fBins*i, 3./fBins*j, (iMassPoint->analysis->GetH2Mass()->GetCellContent(i+1, j+1)-iMassPoint->genMass)/iMassPoint->analysis->GetH2MassError()->GetCellContent(i+1, j+1));
        }
      }
    }

  }
  
  TFile* ensembleFile = new TFile("ensemble.root","recreate");
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->h3Mass->Write();
    iMassPoint->h3MassPull->Write();
  }
}

void TopMass::EvalEnsembleTest() {

  gROOT ->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  gStyle->SetPaintTextFormat(".2f");

  TFile* ensembleFile = new TFile("ensemble.root");
  TTree* tree = (TTree*) ensembleFile->Get("tree");
  
  double genMass;
  tree->Branch("genMass", &genMass, "genMass/D");
  
  TH2F* h2Mass = 0;
  tree->Branch("h2Mass", &h2Mass);
  
  TCanvas* canvas = new TCanvas("canvas", "hadronic top hMass", 500, 500);
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      TF1* gaus = new TF1("gaus", "gaus");
      TString getter("h2Mass.GetCellContent("); getter += i; getter += ","; getter += j; getter += ")";
      
      tree->Fit("gaus", getter, "genMass==172.5");
      
      if (gaus->GetParameter(1) > 150) {
        TString path("plot/"); path += fMethod; path += "/"; path += "ensemble_1725_"; path += i; path += "_"; path += j; path += ".png";
        canvas->Print(path);
      }
    }
  }
}


void TopMass::Calibrate() {
 
  std::vector<TH2F*> hMass;
  std::vector<TH2F*> hMassError;
  
  calibrationAnalyses.push_back(a1665);
  calibrationAnalyses.push_back(a1725);
  calibrationAnalyses.push_back(a1785);
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);
  
  double genMass[] = {166.5, 172.5, 178.5};
  double genMassError[] = {0.0001, 0.0001, 0.0001};
  double hadTopMass[3];
  double hadTopMassError[3];
  
  for(int i = 0; i < 3; i++){
    calibrationAnalyses.at(i)->Analyze();
    hMass.push_back(calibrationAnalyses.at(i)->GetH2Mass());
    hMassError.push_back(calibrationAnalyses.at(i)->GetH2MassError());
  }
    
  canvasFit->cd();

  TGraphErrors* ghadTopMass;
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      if (hMass.at(0)->GetCellContent(i+1, j+1) > 0 && hMass.at(2)->GetCellContent(i+1, j+1) > 0
          && hMassError.at(0)->GetCellContent(i+1, j+1) > 0 && hMassError.at(2)->GetCellContent(i+1, j+1) > 0) {
        for (int k = 0; k < 3; k++) {
          hadTopMass[k] = hMass.at(k)->GetCellContent(i+1, j+1);
          hadTopMassError[k] = hMassError.at(k)->GetCellContent(i+1, j+1);
        }
        ghadTopMass = new TGraphErrors(3, hadTopMass, genMass, hadTopMassError, genMassError);
        ghadTopMass->Draw("A*");
        
        TF1* linearFit = new TF1("linearFit", "172.5+[0]+(x-172.5)*[1]");        
        ghadTopMass->Fit("linearFit");
        
        TString path("plot/"); path += fMethod; path += "/"; path += "fit_"; path += i; path += "_"; path += j; path += ".png";
        canvasFit->Print(path);
        
        for (int l = 0; l < 2; l++) {
          fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        }
      }
    }
  }
}



TH2F* TopMass::Measure(Analysis* a) {
  a->Analyze();
  
  TCanvas* canvas = new TCanvas("canvas", "Hadronic top mass", 1000, 500);
  
  TH2F* hMass = a->GetH2Mass();
  TH2F* hMassError = a->GetH2MassError();

  TH2F* hMassCalibrated = a->GetH2MassCalibrated();
  TH2F* hMassErrorCalibrated = a->GetH2MassErrorCalibrated();
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      if (fCalibFitParameter[i][j][0] && fCalibFitParameter[i][j][1]) {
        double calib = 172.5 + fCalibFitParameter[i][j][0] + fCalibFitParameter[i][j][1] * (hMass->GetCellContent(i+1, j+1) - 172.5);
        double caliberror = sqrt(pow(fCalibFitParError[i][j][0], 2) + pow((hMass->GetCellContent(i+1, j+1) - 172.5)*fCalibFitParError[i][j][1], 2) + pow(fCalibFitParameter[i][j][1]*hMassError->GetCellContent(i+1, j+1), 2));
        
        std::cout << "Measured TopMass: " << calib << " +/- " << caliberror << " GeV" << std::endl;
        
        hMassCalibrated->SetCellContent(i+1, j+1, calib);
        hMassErrorCalibrated->SetCellContent(i+1, j+1, caliberror);
      }
    }
  }
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hMassCalibrated->Draw("COLZ,TEXT");
  hMassCalibrated->SetAxisRange(hMassCalibrated->GetMinimum(150), hMassCalibrated->GetMaximum(200), "Z");
  
  canvas->cd(2);
  hMassErrorCalibrated->Draw("COLZ,TEXT");
  hMassErrorCalibrated->SetAxisRange(0.05, 5, "Z");
  
  TString path("plot/"); path += fMethod; path += "_"; path += a->GetIdentifier(); path +="_calibrated.eps";
  canvas->Print(path);
  
  return hMassCalibrated;
}


void TopMass::Systematics() {
  TH2F* hMassJESdown = Measure(a1725_jes_down);
  TH2F* hMassJESnorm = Measure(a1725);
  TH2F* hMassJESup   = Measure(a1725_jes_up);
  
  hMassJESup->Add(hMassJESnorm, -1.000001);
  hMassJESnorm->Add(hMassJESdown, -1.000001);
  
  hMassJESup->SetTitle("JES up MassError");
  hMassJESnorm->SetTitle("JES down MassError");
  
  TCanvas* canvas = new TCanvas("canvas", "Hadronic top mass JES error", 1000, 500);
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hMassJESnorm->Draw("COLZ,TEXT");
  hMassJESnorm->SetAxisRange(0.05, 5, "Z");
  
  canvas->cd(2);
  hMassJESup->Draw("COLZ,TEXT");
  hMassJESup->SetAxisRange(0.05, 5, "Z");
  
  TString path("plot/"); path += fMethod; path += "_jeserror.eps";
  canvas->Print(path);
}

int main(int argc, char** argv)
{
  if (argc > 1) {
    TopMass* top = new TopMass(argv[1], 6, 500);
  }
  else {
    TopMass* top = new TopMass("GenMatch", 6, 500);
  }
}
