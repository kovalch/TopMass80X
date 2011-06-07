#include "TopMass.h"

TopMass::TopMass(TString method, int bins, double lumi) : fMethod(method), fBins(bins), fLumi(lumi) {
  
  //aSim = new Analysis("sim", "root/analyzeTop_1725.root", fMethod, fBins, fLumi);
  
  //QuickCalibration();
  //LoadXML();
  //QuickSystematics();
  
  WriteEnsembleTest(false);
  //EvalEnsembleTest(true);
  //Measure(aSim);
  
  /*
  Analysis* a1725 = new Analysis("1725", "root/analyzeTop_1725.root", fMethod, fBins, 21);
  Measure(a1725);
  //*/
  
  /*
  Analysis* run2010B = new Analysis("run2010B", "root/analyzeTop_Run2010B.root", fMethod, fBins, 36);
  Measure(run2010B);
  //*/
  
}


void TopMass::WriteEnsembleTest(bool readCalibration) {
  if (readCalibration) LoadXML();

  massPoint m1665(166.5, "1665");
  massPoint m1725(172.5, "1725");
  massPoint m1785(178.5, "1785");
   
  massPoints.push_back(m1665);
  massPoints.push_back(m1725);
  massPoints.push_back(m1785);
  
  int nEnsembles = 10;
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->analysis = new Analysis(iMassPoint->identifier, iMassPoint->fileName, fMethod, fBins, fLumi);
    iMassPoint->h3Mass = new TH3F("h3Mass_" + iMassPoint->identifier, "h3Mass_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, 150, 200);
    iMassPoint->h3MassError = new TH3F("h3MassError_" + iMassPoint->identifier, "h3MassError_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, 0, 10);
    iMassPoint->h3MassPull = new TH3F("h3MassPull_" + iMassPoint->identifier, "h3MassPull_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, -5, 5);
    for (int n = 0; n < nEnsembles; n++) {
      iMassPoint->analysis->Analyze(true);
      for (int i = 0; i < fBins; i++) {
        for (int j = 0; j < fBins; j++) {
          double mass = iMassPoint->analysis->GetH2Mass()->GetCellContent(i+1, j+1);
          double massError = iMassPoint->analysis->GetH2MassError()->GetCellContent(i+1, j+1);
          double massPull = (mass-iMassPoint->genMass)/massError;
          
          if (readCalibration) {
            if (fCalibFitParameter[i][j][0] && fCalibFitParameter[i][j][1]) {
              massError = sqrt(pow((1-fCalibFitParameter[i][j][1])*massError, 2) + pow(fCalibFitParError[i][j][0], 2) + pow(fCalibFitParError[i][j][1]*(172.5-mass), 2));
              mass = mass - fCalibFitParameter[i][j][0] + (fCalibFitParameter[i][j][1]*(172.5-mass));
              massPull = (mass-iMassPoint->genMass)/massError;
            }
            else {
              mass = 0;
              massError = 0;
              massPull = 999;
            }
          }
          
          
          iMassPoint->h3Mass->Fill(3./fBins*i, 3./fBins*j, mass);
          iMassPoint->h3MassError->Fill(3./fBins*i, 3./fBins*j, massError);
          iMassPoint->h3MassPull->Fill(3./fBins*i, 3./fBins*j, massPull);
        }
      }
    }

  }
  
  TFile* ensembleFile = new TFile("ensemble.root","recreate");
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->h3Mass->Write();
    iMassPoint->h3MassError->Write();
    iMassPoint->h3MassPull->Write();
  }
}

void TopMass::EvalEnsembleTest(bool writeCalibration) {
  Helper* helper = new Helper(fBins);
  helper->SetTDRStyle();
  gStyle->SetOptFit(0);

  TiXmlDocument doc;
  TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
  doc.LinkEndChild( decl );
  
  TiXmlElement* calibration = new TiXmlElement( "calibration" );
  doc.LinkEndChild( calibration );
  
  massPoint m1665(166.5, "1665");
  massPoint m1725(172.5, "1725");
  massPoint m1785(178.5, "1785");
  
  int nEnsemble = 1000;
  
  m1665.genLumi = 2250;
  m1725.genLumi = 6100;
  m1785.genLumi = 1900;
   
  massPoints.push_back(m1665);
  massPoints.push_back(m1725);
  massPoints.push_back(m1785);
  
  TFile* ensembleFile = new TFile("root/ensemble.root");
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->h2Mass = helper->GetH2("Mass");
    iMassPoint->h2MassError = helper->GetH2("MassError");
    iMassPoint->h3Mass = (TH3F*) ensembleFile->Get("h3Mass_" + iMassPoint->identifier);
    iMassPoint->h3MassPull = (TH3F*) ensembleFile->Get("h3MassPull_" + iMassPoint->identifier);
  }
    
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      TCanvas* canvas = new TCanvas("canvas", "hadronic top hMass", 1000, 1000);
      canvas->Divide(2,2);
      
      TVectorD hadTopMass(massPoints.size());
      TVectorD hadTopMassMeanError(massPoints.size());
      TVectorD hadTopMassError(massPoints.size());
      TVectorD hadTopMassBias(massPoints.size());
      TVectorD hadTopMassChi2NDF(massPoints.size());
      
      TVectorD hadTopMassPullWidth(massPoints.size());
      TVectorD hadTopMassPullWidthError(massPoints.size());
      
      for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
        TCanvas* canvasTemp = new TCanvas("canvasTemp", "hadronic top hMass");
        canvasTemp->cd();
        
        int k = iMassPoint - massPoints.begin();
        
        iMassPoint->hMass = iMassPoint->h3Mass->ProjectionZ("hMass_" + iMassPoint->identifier, i+1, i+1, j+1, j+1);
        iMassPoint->hMass->Rebin(2);
        iMassPoint->hMass->GetXaxis()->SetTitle("m_{t}");
        iMassPoint->hMass->GetYaxis()->SetTitle("Pseudo-experiments");
        iMassPoint->hMass->SetFillColor(kRed - 11 + massPoints.size() - k);
        
        iMassPoint->hMassPull = iMassPoint->h3MassPull->ProjectionZ("hMassPull_" + iMassPoint->identifier, i+1, i+1, j+1, j+1);
        iMassPoint->hMassPull->Rebin(2);
        iMassPoint->hMassPull->GetXaxis()->SetTitle("m_{t} pull");
        iMassPoint->hMassPull->GetYaxis()->SetTitle("Pseudo-experiments");
        iMassPoint->hMassPull->SetFillColor(kRed - 11 + massPoints.size() - k);
        
        TF1* gaus = new TF1("gaus", "gaus");
        TF1* gausPull = new TF1("gausPull", "gaus");
        
        gaus->SetLineWidth(1);
        gaus->SetLineColor(kBlack);
        gaus->SetLineWidth(2);
        
        gausPull->SetLineWidth(1);
        gausPull->SetLineColor(kBlack);
        gausPull->SetLineWidth(2);
        
        iMassPoint->hMass->Fit("gaus");
        
        iMassPoint->hMassPull->Fit("gausPull");
        
        hadTopMass[k] = gaus->GetParameter(1);
        hadTopMassMeanError[k] = gaus->GetParError(1)*TMath::Sqrt(1+nEnsemble*fLumi/iMassPoint->genLumi);
        hadTopMassError[k] = gaus->GetParameter(2);
        hadTopMassBias[k] = gaus->GetParameter(1)-iMassPoint->genMass;
        hadTopMassChi2NDF[k] = gaus->GetChisquare()/gaus->GetNDF();
        
        hadTopMassPullWidth[k] = gausPull->GetParameter(2);
        hadTopMassPullWidthError[k] = gausPull->GetParError(2)*TMath::Sqrt(1+nEnsemble*TMath::Power(fLumi/iMassPoint->genLumi, 2));
        
        iMassPoint->h2Mass->SetCellContent(i+1, j+1, hadTopMass[k]);
        iMassPoint->h2MassError->SetCellContent(i+1, j+1, hadTopMassError[k]);
      }
          
      if (hadTopMassMeanError > 0 && hadTopMassChi2NDF < 10) {
        canvas->cd(1);
        
        TLegend *leg0 = new TLegend(0.65, 0.7, 0.95, 0.9);
        leg0->SetFillStyle(0);
        leg0->SetBorderSize(0);
        
        for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
          if (iMassPoint == massPoints.begin()) {
            iMassPoint->hMass->Draw();
          }
          else {
            iMassPoint->hMass->Draw("SAME");
          }
          TString legend("m_{t,gen} = "); legend += iMassPoint->genMass; legend += " GeV";
          leg0->AddEntry( iMassPoint->hMass, legend, "F");
        }
        
        leg0->Draw();
        
        canvas->cd(3);
        for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
          if (iMassPoint == massPoints.begin()) {
            iMassPoint->hMassPull->Draw();
          }
          else {
            iMassPoint->hMassPull->Draw("SAME");
          }
        }
        
        leg0->Draw();
        
        gStyle->SetOptFit(1);
        
        canvas->cd(2);
        
        TGraphErrors* gBias = new TGraphErrors(hadTopMass, hadTopMassBias, hadTopMassMeanError, hadTopMassMeanError);
        if (hadTopMassBias.Min() > 0) {
          gBias->GetYaxis()->SetRangeUser(0.5, hadTopMassBias.Max()+hadTopMassMeanError.Max()+0.5);
        }
        else if (hadTopMassBias.Max() < 0) {
          gBias->GetYaxis()->SetRangeUser(hadTopMassBias.Min()-hadTopMassMeanError.Max()-0.5, 0.5);
        }
        gBias->GetXaxis()->SetTitle("m_{t}");
        gBias->GetYaxis()->SetRangeUser(-6, 6);
        gBias->GetYaxis()->SetTitle("bias");
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
        
        canvas->cd(4);
        
        TGraphErrors* gPull = new TGraphErrors(hadTopMass, hadTopMassPullWidth, hadTopMassMeanError, hadTopMassPullWidthError);
        gPull->GetXaxis()->SetTitle("m_{t}");
        gPull->GetYaxis()->SetRangeUser(0.75, 1.25);
        gPull->GetYaxis()->SetTitle("pull width");
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
        
        TString path("plot/"); path += fMethod; path += "/"; path += "ensembletest_"; path += i; path += "_"; path += j; path += ".png";
        canvas->Print(path);
        
        for (int l = 0; l < 2; l++) {
          fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        }
        
        if (writeCalibration) {
          TiXmlElement* bin = new TiXmlElement( "bin" );
          calibration->LinkEndChild( bin );
          bin->SetAttribute("binx", i);
          bin->SetAttribute("biny", j);
          bin->SetDoubleAttribute("p0", linearFit->GetParameter(0));
          bin->SetDoubleAttribute("p0error", linearFit->GetParError(0));
          bin->SetDoubleAttribute("p1", linearFit->GetParameter(1));
          bin->SetDoubleAttribute("p1error", linearFit->GetParError(1));
        }
      }

      delete canvas;
    }
  }
  
  gStyle->SetPadRightMargin(0.15);
  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 450);
  canvas->Divide(2,1);
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    canvas->cd(1);
    iMassPoint->h2Mass->Draw("COLZ, TEXT");
    iMassPoint->h2Mass->SetAxisRange(160, 185, "Z");
    
    canvas->cd(2);
    iMassPoint->h2MassError->Draw("COLZ, TEXT");
    iMassPoint->h2MassError->SetAxisRange(0.05, 5, "Z");
    
    TString path("plot/"); path += fMethod; path += "_"; path += "ensembletest_"; path += iMassPoint->identifier; path += ".eps";
    canvas->Print(path);
  }
  
  ensembleFile->Close();
  
  if (writeCalibration) doc.SaveFile( "calibration.xml" );
}


void TopMass::QuickCalibration() {
 
  std::vector<TH2F*> hMass;
  std::vector<TH2F*> hMassError;
  
  Analysis* a1665 = new Analysis("1665", "root/analyzeTop_1665.root", fMethod, fBins, 100000);
  Analysis* a1725 = new Analysis("1725", "root/analyzeTop_1725.root", fMethod, fBins, 100000);
  Analysis* a1785 = new Analysis("1785", "root/analyzeTop_1785.root", fMethod, fBins, 100000);
  
  calibrationAnalyses.push_back(a1665);
  calibrationAnalyses.push_back(a1725);
  calibrationAnalyses.push_back(a1785);
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);
  
  double genMass[] = {166.5, 172.5, 178.5};
  double genMassError[] = {0.0001, 0.0001, 0.0001};
  double hadTopMass[3];
  double hadTopMassBias[3];
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
          hadTopMassBias[k] = hMass.at(k)->GetCellContent(i+1, j+1) - genMass[k];
          hadTopMassError[k] = hMassError.at(k)->GetCellContent(i+1, j+1);
        }
        ghadTopMass = new TGraphErrors(3, hadTopMass, hadTopMassBias, hadTopMassError, hadTopMassError);
        ghadTopMass->Draw("A*");
        
        ghadTopMass->GetYaxis()->SetRangeUser(-6, 6);
        
        TF1* linearFit = new TF1("linearFit", "[0]+(x-172.5)*[1]");    
        ghadTopMass->Fit("linearFit");
        
        for (int l = 0; l < 2; l++) {
          fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        }
        
        TString path("plot/"); path += fMethod; path += "/"; path += "fit_"; path += i; path += "_"; path += j; path += ".eps";
        canvasFit->Print(path);
      }
    }
  }
  
  Measure(a1665);
  Measure(a1725);
  Measure(a1785);
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

      
      
      double mass = hMass->GetCellContent(i+1, j+1);
      double massError = hMassError->GetCellContent(i+1, j+1);
      
      std::cout << "Measured TopMass: " << mass << " +/- " << massError << " GeV" << std::endl;
    
      if (fCalibFitParameter[i][j][0] && fCalibFitParameter[i][j][1] && mass > 0) {
        massError = sqrt(pow((1-fCalibFitParameter[i][j][1])*massError, 2) + pow(fCalibFitParError[i][j][0], 2) + pow(fCalibFitParError[i][j][1]*(172.5-mass), 2));
        mass = mass - fCalibFitParameter[i][j][0] + (fCalibFitParameter[i][j][1]*(172.5-mass));
        
        std::cout << "Calibrated TopMass: " << mass << " +/- " << massError << " GeV" << std::endl;
      }
      else {
        mass = 0;
        massError = 0;
        std::cout << "No Calibration data or mass available" << std::endl;
      }
        
      hMassCalibrated->SetCellContent(i+1, j+1, mass);
      hMassErrorCalibrated->SetCellContent(i+1, j+1, massError);
      
    }
  }
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hMassCalibrated->Draw("COLZ,TEXT");
  hMassCalibrated->SetAxisRange(160, 185, "Z");
  
  canvas->cd(2);
  hMassErrorCalibrated->Draw("COLZ,TEXT");
  hMassErrorCalibrated->SetAxisRange(0.05, 5, "Z");
  
  TString path("plot/"); path += fMethod; path += "_"; path += a->GetIdentifier(); path +="_calibrated.eps";
  canvas->Print(path);
  
  return hMassCalibrated;
}


void TopMass::QuickSystematics() {

  Analysis* a1665 = new Analysis("1665", "root/analyzeTop_1665.root", fMethod, fBins, 10000);
  Analysis* a1725 = new Analysis("1725", "root/analyzeTop_1725.root", fMethod, fBins, 10000);
  Analysis* a1785 = new Analysis("1785", "root/analyzeTop_1785.root", fMethod, fBins, 10000);

  Analysis* a1665_jes_up = new Analysis("1665_jes_up", "root/analyzeTop_1665_jes_up.root", fMethod, fBins, 10000);
  Analysis* a1725_jes_up = new Analysis("1725_jes_up", "root/analyzeTop_1725_jes_up.root", fMethod, fBins, 10000);
  Analysis* a1785_jes_up = new Analysis("1785_jes_up", "root/analyzeTop_1785_jes_up.root", fMethod, fBins, 10000);
  
  Analysis* a1665_jes_down = new Analysis("1665_jes_down", "root/analyzeTop_1665_jes_down.root", fMethod, fBins, 10000);
  Analysis* a1725_jes_down = new Analysis("1725_jes_down", "root/analyzeTop_1725_jes_down.root", fMethod, fBins, 10000);
  Analysis* a1785_jes_down = new Analysis("1785_jes_down", "root/analyzeTop_1785_jes_down.root", fMethod, fBins, 10000);
  
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

void TopMass::LoadXML() {
  TString xmlFileName;
  if (fexists("/scratch/hh/lustre/cms/user/mseidel/calibration.xml")) {
    xmlFileName = "/scratch/hh/lustre/cms/user/mseidel/calibration.xml";
  }
  else xmlFileName = "calibration.xml";

  TiXmlDocument doc(xmlFileName);
  bool loadOkay = doc.LoadFile();
  
  TiXmlElement *pRoot, *pParm;
  
  pRoot = doc.FirstChildElement( "calibration" );
  
  pParm = pRoot->FirstChildElement("bin");
  while ( pParm ) {
    int i, j;
    double p0, p0error, p1, p1error;
    pParm->Attribute("binx", &i);
    pParm->Attribute("biny", &j);
    pParm->Attribute("p0", &p0);
    pParm->Attribute("p0error", &p0error);
    pParm->Attribute("p1", &p1);
    pParm->Attribute("p1error", &p1error);
    
    fCalibFitParameter[i][j][0] = p0;
    fCalibFitParameter[i][j][1] = p1;
    fCalibFitParError[i][j][0]  = p0error;
    fCalibFitParError[i][j][1]  = p1error;
    
    pParm = pParm->NextSiblingElement( "bin" );
  }
}

int main(int argc, char** argv)
{
  if (argc > 1) {
    TopMass* top = new TopMass(argv[1], 1, 21);
  }
  else {
    TopMass* top = new TopMass("GenMatch", 6, 500);
  }
}
