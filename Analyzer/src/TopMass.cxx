#include "TopMass.h"

TopMass::TopMass(TString method, int bins, double lumi) : fMethod(method), fBins(bins), fLumi(lumi) {


/*  a1665_jes_up = new Analysis("1665_jes_up", "root/analyzeTop_1665_jes_up.root", fMethod, fBins, fLumi);
  a1725_jes_up = new Analysis("1725_jes_up", "root/analyzeTop_1725_jes_up.root", fMethod, fBins, fLumi);
  a1785_jes_up = new Analysis("1785_jes_up", "root/analyzeTop_1785_jes_up.root", fMethod, fBins, fLumi);
  
  a1665_jes_down = new Analysis("1665_jes_down", "root/analyzeTop_1665_jes_down.root", fMethod, fBins, fLumi);
  a1725_jes_down = new Analysis("1725_jes_down", "root/analyzeTop_1725_jes_down.root", fMethod, fBins, fLumi);
  a1785_jes_down = new Analysis("1785_jes_down", "root/analyzeTop_1785_jes_down.root", fMethod, fBins, fLumi);*/
  
  aSim = new Analysis("sim", "root/analyzeTop_1725.root", fMethod, fBins, fLumi);
  
  //Calibrate();
  //Systematics();
  
  //WriteEnsembleTestTree();
  EvalEnsembleTest();
  Measure(aSim);
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
  int nEnsemble = 10000;

  gROOT ->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  gStyle->SetPaintTextFormat(".2f");
  
  massPoint m1665(166.5, "1665");
  massPoint m1725(172.5, "1725");
  massPoint m1785(178.5, "1785");
  
  m1665.genLumi = 2250;
  m1725.genLumi = 6100;
  m1785.genLumi = 1900;
   
  massPoints.push_back(m1665);
  massPoints.push_back(m1725);
  massPoints.push_back(m1785);
  
  TFile* ensembleFile = new TFile("root/ensemble10k_6.root");
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->h3Mass = (TH3F*) ensembleFile->Get("h3Mass_" + iMassPoint->identifier);
    iMassPoint->h3MassPull = (TH3F*) ensembleFile->Get("h3MassPull_" + iMassPoint->identifier);
  }
    
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      TCanvas* canvas = new TCanvas("canvas", "hadronic top hMass", 2000, 1000);
      canvas->Divide(4,2);
      
      TVectorD genMass(massPoints.size());
      TVectorD genMassError(massPoints.size());
      TVectorD hadTopMassBias(massPoints.size());
      TVectorD hadTopMassError(massPoints.size());
      TVectorD hadTopMassChi2NDF(massPoints.size());
      
      TVectorD hadTopMassPullWidth(massPoints.size());
      TVectorD hadTopMassPullWidthError(massPoints.size());
      
      for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
        int k = iMassPoint - massPoints.begin();
        
        iMassPoint->h3Mass = (TH3F*) ensembleFile->Get("h3Mass_" + iMassPoint->identifier);
        iMassPoint->h3MassPull = (TH3F*) ensembleFile->Get("h3MassPull_" + iMassPoint->identifier);
        
        TH1* hMass = iMassPoint->h3Mass->ProjectionZ("hMass", i+1, i+1, j+1, j+1);
        TH1* hMassPull = iMassPoint->h3MassPull->ProjectionZ("hMassPull", i+1, i+1, j+1, j+1);
        
        TF1* gaus = new TF1("gaus", "gaus");
        TF1* gausPull = new TF1("gausPull", "gaus");
        
        gaus->SetLineWidth(1);
        gaus->SetLineColor(kRed);
        
        gausPull->SetLineWidth(1);
        gausPull->SetLineColor(kRed);
        
        canvas->cd(1+k);
        hMass->Fit("gaus");
        
        canvas->cd(5+k);
        hMassPull->Fit("gausPull");
        
        genMass[k] = gaus->GetParameter(1);
        genMassError[k] = gaus->GetParError(1)*TMath::Sqrt(1+nEnsemble*fLumi/iMassPoint->genLumi);
        hadTopMassBias[k] = gaus->GetParameter(1)-iMassPoint->genMass;
        hadTopMassError[k] = gaus->GetParError(1)*TMath::Sqrt(1+nEnsemble*fLumi/iMassPoint->genLumi);
        hadTopMassChi2NDF[k] = gaus->GetChisquare()/gaus->GetNDF();
        
        hadTopMassPullWidth[k] = gausPull->GetParameter(2);
        hadTopMassPullWidthError[k] = gausPull->GetParError(2)*TMath::Sqrt(1+nEnsemble*TMath::Power(fLumi/iMassPoint->genLumi, 2));        
      }
          
      if (hadTopMassError > 0 && hadTopMassChi2NDF < 10) {
        canvas->cd(4);
        
        TGraphErrors* gBias = new TGraphErrors(genMass, hadTopMassBias, genMassError, hadTopMassError);
        if (hadTopMassBias.Min() > 0) {
          gBias->GetYaxis()->SetRangeUser(0.5, hadTopMassBias.Max()+hadTopMassError.Max()+0.5);
        }
        else if (hadTopMassBias.Max() < 0) {
          gBias->GetYaxis()->SetRangeUser(hadTopMassBias.Min()-hadTopMassError.Max()-0.5, 0.5);
        }
        gBias->SetMarkerStyle(2);
        gBias->SetMarkerSize(0.4);
        gBias->Draw("AP");
        
        TF1* zero_line = new TF1("zero_line", "0", 150, 200);
        zero_line->SetLineWidth(1);
        zero_line->SetLineStyle(3);
        zero_line->SetLineColor(kBlack);
        zero_line->Draw("SAME");
        
        TF1* linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");      
        linearFit->SetLineWidth(1);
        linearFit->SetLineColor(kRed);
        
        gBias->Fit("linearFit");
        
        canvas->cd(8);
        
        TGraphErrors* gPull = new TGraphErrors(genMass, hadTopMassPullWidth, genMassError, hadTopMassPullWidthError);
        gPull->GetYaxis()->SetRangeUser(0, 2);
        gPull->SetMarkerStyle(2);
        gPull->SetMarkerSize(0.4);
        gPull->Draw("AP");
        
        TF1* unity_line = new TF1("unity_line", "1", 150, 200);
        unity_line->SetLineWidth(1);
        unity_line->SetLineStyle(3);
        unity_line->SetLineColor(kBlack);
        unity_line->Draw("SAME");
        
        TF1* constFit = new TF1("constFit", "[0]");
        constFit->SetLineWidth(1);
        constFit->SetLineColor(kRed);
        
        gPull->Fit("constFit");
        
        TString path("plot/"); path += fMethod; path += "/"; path += "ensembletest_"; path += i; path += "_"; path += j; path += "_uncalibrated.png";
        canvas->Print(path);
        
        for (int l = 0; l < 2; l++) {
          fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        }
      }

      delete canvas;
    }
  }
  
  ensembleFile->Close();
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
        double calib = hMass->GetCellContent(i+1, j+1) - fCalibFitParameter[i][j][0] + (fCalibFitParameter[i][j][1]*(172.5-hMass->GetCellContent(i+1, j+1)));
        double caliberror = sqrt(pow((1-fCalibFitParameter[i][j][1])*hMassError->GetCellContent(i+1, j+1), 2) + pow(fCalibFitParError[i][j][0], 2) + pow(fCalibFitParError[i][j][1]*(172.5-hMass->GetCellContent(i+1, j+1)), 2));
        
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
