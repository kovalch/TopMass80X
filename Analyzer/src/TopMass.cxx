#include "TopMass.h"

TopMass::TopMass(TString method, int bins, double lumi) : fMethod(method), fBins(bins), fLumi(lumi) {
  
  //QuickCalibration();
  //LoadXML();
  //QuickSystematics();
  
  WriteEnsembleTest(false);
  //EvalEnsembleTest(true);
  //Measure(aSim);
  
  //analyses.push_back(new Analysis("Summer11_1725_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root", fMethod, fBins, 0));
  
  /*  Systematic samples
  analyses.push_back(new Analysis("1725_flavordown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_flavor:down/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_flavorup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_flavor:up/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_jesdown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_jes:down/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_jesup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_jes:up/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_resdown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_res:down/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_resup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_res:up/analyzeTop.root", fMethod, fBins, 0));
  //*/
  
  /*  Uncorrelated systematic samples
  analyses.push_back(new Analysis("1725_scaledown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_scaledown/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_scaleup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_scaleup/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_matchingdown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_matchingdown/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("1725_matchingup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_matchingup/analyzeTop.root", fMethod, fBins, 0));
  //*/
  
  /*  Spring11
  analyses.push_back(new Analysis("Spring11_D6T", "/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_D6T/analyzeTop.root", fMethod, fBins, 0));
  analyses.push_back(new Analysis("Spring11_Z2", "/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_Z2/analyzeTop.root", fMethod, fBins, 0));
  //*/
  
  for (int i = 0; i < analyses.size(); i++) {
    Measure(analyses[i]);
  }

  /* DATA, full 2011
  Analysis* aRun2011 = new Analysis("Run2011", "/scratch/hh/current/cms/user/mseidel/Run2011/analyzeTop.root", fMethod, fBins, 0);
  Measure(aRun2011);
  //*/
}


void TopMass::WriteEnsembleTest(bool readCalibration) {
  if (readCalibration) LoadXML();
  
  massPoint m1665j096(166.5, 0.96, "1665_0.96"); massPoints.push_back(m1665j096);
  massPoint m1665j100(166.5, 1.00, "1665_1.00"); massPoints.push_back(m1665j100);
  massPoint m1665j104(166.5, 1.04, "1665_1.04"); massPoints.push_back(m1665j104);
  
  massPoint m1725j096(172.5, 0.96, "1725_0.96"); massPoints.push_back(m1725j096);
  massPoint m1725j100(172.5, 1.00, "1725_1.00"); massPoints.push_back(m1725j100);
  massPoint m1725j104(172.5, 1.04, "1725_1.04"); massPoints.push_back(m1725j104);
  
  massPoint m1785j096(178.5, 0.96, "1785_0.96"); massPoints.push_back(m1785j096);
  massPoint m1785j100(178.5, 1.00, "1785_1.00"); massPoints.push_back(m1785j100);
  massPoint m1785j104(178.5, 1.04, "1785_1.04"); massPoints.push_back(m1785j104);
  
  int nEnsembles = 4;
  
  double genMass, mass, massError, massPull, genJES, JES, JESError, JESPull;
  TTree* tree = new TTree("tree", "tree");
  tree->Branch("genMass", &genMass, "genMass/D");
  tree->Branch("mass", &mass, "mass/D");
  tree->Branch("massError", &massError, "massError/D");
  tree->Branch("massPull", &massPull, "massPull/D");
  tree->Branch("genJES", &genJES, "genJES/D");
  tree->Branch("JES", &JES, "JES/D");
  tree->Branch("JESError", &JESError, "JESError/D");
  tree->Branch("JESPull", &JESPull, "JESPull/D");
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->analysis = new Analysis(iMassPoint->identifier, iMassPoint->fileName, fMethod, fBins, fLumi);
    
    iMassPoint->h3Mass = new TH3F("h3Mass_" + iMassPoint->identifier, "h3Mass_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, 150, 200);
    iMassPoint->h3MassError = new TH3F("h3MassError_" + iMassPoint->identifier, "h3MassError_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 200, 0, 2);
    iMassPoint->h3MassPull = new TH3F("h3MassPull_" + iMassPoint->identifier, "h3MassPull_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, -5, 5);
    
    iMassPoint->h3JES = new TH3F("h3JES_" + iMassPoint->identifier, "h3JES_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, 0.9, 1.1);
    iMassPoint->h3JESError = new TH3F("h3JESError_" + iMassPoint->identifier, "h3JESError_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 200, 0, 0.02);
    iMassPoint->h3JESPull = new TH3F("h3JESPull_" + iMassPoint->identifier, "h3JESPull_" + iMassPoint->identifier, fBins, 0, 3, fBins, 0, 3, 100, -5, 5);
    
    for (int n = 0; n < nEnsembles; n++) {
      iMassPoint->analysis->Analyze(true);
      for (int i = 0; i < fBins; i++) {
        for (int j = 0; j < fBins; j++) {
          genMass   = iMassPoint->genMass;
          mass      = iMassPoint->analysis->GetH2Mass()->GetCellContent(i+1, j+1);
          massError = iMassPoint->analysis->GetH2MassError()->GetCellContent(i+1, j+1);
          massPull  = (mass - genMass)/massError;
          
          genJES    = iMassPoint->genJES;
          JES       = iMassPoint->analysis->GetH2JES()->GetCellContent(i+1, j+1);
          JESError  = iMassPoint->analysis->GetH2JESError()->GetCellContent(i+1, j+1);
          JESPull   = (JES - genJES)/JESError;
          
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
          
          iMassPoint->h3JES->Fill(3./fBins*i, 3./fBins*j, JES);
          iMassPoint->h3JESError->Fill(3./fBins*i, 3./fBins*j, JESError);
          iMassPoint->h3JESPull->Fill(3./fBins*i, 3./fBins*j, JESPull);
          
          tree->Fill();
        }
      }
    }

  }
  
  TFile* ensembleFile = new TFile("ensemble.root","recreate");
  
  for (iMassPoint = massPoints.begin(); iMassPoint != massPoints.end(); ++iMassPoint) {
    iMassPoint->h3Mass->Write();
    iMassPoint->h3MassError->Write();
    iMassPoint->h3MassPull->Write();
    
    iMassPoint->h3JES->Write();
    iMassPoint->h3JESError->Write();
    iMassPoint->h3JESPull->Write();
  }
  
  tree->Write();
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
  
  massPoint m1665(166.5, 1., "1665");
  massPoint m1725(172.5, 1., "1725");
  massPoint m1785(178.5, 1., "1785");
  
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
	enum styles          { kDown, kNominal, kUp};
	int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
	int marker_  [ 3 ] = { 23, 20, 22};
	
  Helper* helper = new Helper(fBins);
  helper->SetTDRStyle();  
  
  std::vector< std::vector<TH2F*> > hMass;
  std::vector< std::vector<TH2F*> > hMassError;
  
  std::vector< std::vector<TH2F*> > hJES;
  std::vector< std::vector<TH2F*> > hJESError;
  
  for (int i = 0; i < 3; i++) {
    calibrationAnalyses.push_back(std::vector<Analysis*>());
    hMass              .push_back(std::vector<TH2F*>());
    hMassError         .push_back(std::vector<TH2F*>());
    hJES               .push_back(std::vector<TH2F*>());
    hJESError          .push_back(std::vector<TH2F*>());
  }
  
	
  Analysis* a1615_096 = new Analysis("1615_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1635_096 = new Analysis("1635_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1665_096 = new Analysis("1665_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1695_096 = new Analysis("1695_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_0.96/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1725_096 = new Analysis("1725_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1755_096 = new Analysis("1755_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_0.96/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1785_096 = new Analysis("1785_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1815_096 = new Analysis("1815_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1845_096 = new Analysis("1845_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_0.96/analyzeTop.root", fMethod, fBins, 0);
  
  Analysis* a1615_100 = new Analysis("1615_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1635_100 = new Analysis("1635_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1665_100 = new Analysis("1665_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1695_100 = new Analysis("1695_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_1.00/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1725_100 = new Analysis("1725_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1755_100 = new Analysis("1755_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_1.00/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1785_100 = new Analysis("1785_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1815_100 = new Analysis("1815_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1845_100 = new Analysis("1845_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_1.00/analyzeTop.root", fMethod, fBins, 0);
  
  Analysis* a1615_104 = new Analysis("1615_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1635_104 = new Analysis("1635_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1665_104 = new Analysis("1665_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1695_104 = new Analysis("1695_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_1.04/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1725_104 = new Analysis("1725_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1755_104 = new Analysis("1755_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_1.04/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1785_104 = new Analysis("1785_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1815_104 = new Analysis("1815_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1845_104 = new Analysis("1845_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_1.04/analyzeTop.root", fMethod, fBins, 0);
	  
  calibrationAnalyses[0].push_back(a1615_096);
  calibrationAnalyses[0].push_back(a1635_096);
  calibrationAnalyses[0].push_back(a1665_096);
  calibrationAnalyses[0].push_back(a1695_096);
  calibrationAnalyses[0].push_back(a1725_096);
  calibrationAnalyses[0].push_back(a1755_096);
  calibrationAnalyses[0].push_back(a1785_096);
  calibrationAnalyses[0].push_back(a1815_096);
  calibrationAnalyses[0].push_back(a1845_096);
  
  calibrationAnalyses[1].push_back(a1615_100);
  calibrationAnalyses[1].push_back(a1635_100);
  calibrationAnalyses[1].push_back(a1665_100);
  calibrationAnalyses[1].push_back(a1695_100);
  calibrationAnalyses[1].push_back(a1725_100);
  calibrationAnalyses[1].push_back(a1755_100);
  calibrationAnalyses[1].push_back(a1785_100);
  calibrationAnalyses[1].push_back(a1815_100);
  calibrationAnalyses[1].push_back(a1845_100);
  
  calibrationAnalyses[2].push_back(a1615_104);
  calibrationAnalyses[2].push_back(a1635_104);
  calibrationAnalyses[2].push_back(a1665_104);
  calibrationAnalyses[2].push_back(a1695_104);
  calibrationAnalyses[2].push_back(a1725_104);
  calibrationAnalyses[2].push_back(a1755_104);
  calibrationAnalyses[2].push_back(a1785_104);
  calibrationAnalyses[2].push_back(a1815_104);
  calibrationAnalyses[2].push_back(a1845_104);
  
  std::cout << "vectors filled" << std::endl;
  
  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);
  
  double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
  double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
  
  double genJES[]       = {0.96, 1.00, 1.04};
  double genJESError[]  = {1e-6, 1e-6, 1e-6};
  
  double hadTopMass[3][9];
  double hadTopMassBias[3][9];
  double hadTopMassError[3][9];
  
  double JES[3][9];
  double JESBias[3][9];
  double JESError[3][9];
  
  //*
  for(int iJES = 0; iJES < 3; iJES++) {
    for(int iMass = 0; iMass < 9; iMass++) {
      calibrationAnalyses[iJES][iMass]->Analyze();
			
      hMass[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2Mass());
      hMassError[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2MassError());
			
			hJES[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2JES());
      hJESError[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2JESError());
    }
  }
  //*/
  
  /* TEST
  calibrationAnalyses[1][2]->Analyze();
  for(int iJES = 0; iJES < 3; iJES++) {
    for(int iMass = 0; iMass < 9; iMass++) {
      hMass[iJES].push_back(calibrationAnalyses[1][2]->GetH2Mass());
      hMassError[iJES].push_back(calibrationAnalyses[1][2]->GetH2MassError());
			
			hJES[iJES].push_back(calibrationAnalyses[1][2]->GetH2JES());
      hJESError[iJES].push_back(calibrationAnalyses[1][2]->GetH2JESError());
    }
  }
  //*/
  
  canvasFit->cd();

  std::vector<TGraphErrors*> gMass;
	std::vector<TGraphErrors*> gJES;
	
  gStyle->SetOptFit(1);
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      //*
      if (   hMass[0][0]->GetCellContent(i+1, j+1) > 0      && hMass[2][2]->GetCellContent(i+1, j+1) > 0
          && hMassError[0][0]->GetCellContent(i+1, j+1) > 0 && hMassError[2][2]->GetCellContent(i+1, j+1) > 0) {
        for (int iJES = 0; iJES < 3; iJES++) {
          for (int iMass = 0; iMass < 9; iMass++) {
            hadTopMass[iJES][iMass]      = genMass[iMass];
            hadTopMassBias[iJES][iMass]  = hMass[iJES][iMass]->GetCellContent(i+1, j+1) - genMass[iMass];
            hadTopMassError[iJES][iMass] = hMassError[iJES][iMass]->GetCellContent(i+1, j+1);
            
            //JES[iJES][iMass]      = genJES[iJES];
            JESBias[iJES][iMass]  = hJES[iJES][iMass]->GetCellContent(i+1, j+1) - genJES[iJES];
            JESError[iJES][iMass] = hJESError[iJES][iMass]->GetCellContent(i+1, j+1);
          }
        }
        //*/
      
				TMultiGraph *mgJES = new TMultiGraph();
				mgJES->SetTitle(";m_{t,gen};JES_{meas}-JES_{gen}");
				
				for (int iJES = 0; iJES < 3; iJES++) {
					double genMassMod[9];
					for (int i = 0; i < 9; i++) {
						genMassMod[i] = genMass[i] + 0.2*(iJES-1);
					}
					
        	gJES.push_back(new TGraphErrors(9, genMassMod, JESBias[iJES], genMassError, JESError[iJES]));
					
					gJES[iJES]->SetMarkerStyle(marker_[iJES]);
					gJES[iJES]->SetMarkerColor(color_ [iJES]);
					gJES[iJES]->SetLineColor  (color_ [iJES]);
										
					TString sFit("[0]+(x-172.5-"); sFit += 0.2*(iJES-1); sFit += ")*[1]";
					
        	TF1* linearFit = new TF1("linearFit", sFit);
        	linearFit->SetParLimits(1, -1, 1);
					linearFit->SetParNames("offset", "slope");
					linearFit->SetLineColor(color_[iJES]);
        	gJES[iJES]->Fit("linearFit", "EM");
					
					mgJES->Add(gJES[iJES]);

					/*
        	for (int l = 0; l < 2; l++) {
          	fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          	fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        	}
					*/
        }
        
				mgJES->SetMinimum(-0.05);
				mgJES->SetMaximum( 0.05);
				
				mgJES->Draw("AP");
				
        canvasFit->Update();
        
        TPaveStats* stats0 = (TPaveStats*) gJES[0]->GetListOfFunctions()->FindObject("stats");
        stats0->SetTextColor(color_[0]);
        stats0->SetX1NDC(0.16);
        stats0->SetY1NDC(0.7);
        stats0->SetX2NDC(0.4233);
        stats0->SetY2NDC(0.825);
        
        TPaveStats* stats1 = (TPaveStats*) gJES[1]->GetListOfFunctions()->FindObject("stats");
        stats1->SetTextColor(color_[1]);
        stats1->SetX1NDC(0.4233);
        stats1->SetY1NDC(0.7);
        stats1->SetX2NDC(0.6867);
        stats1->SetY2NDC(0.825);
        
        TPaveStats* statsGlobal = (TPaveStats*) gJES[2]->GetListOfFunctions()->FindObject("stats");
        
        TPaveStats* stats2 = (TPaveStats*) statsGlobal->Clone();
        stats2->SetTextColor(color_[2]);
        stats2->SetX1NDC(0.6867);
        stats2->SetY1NDC(0.7);
        stats2->SetX2NDC(0.95);
        stats2->SetY2NDC(0.825);
        stats2->Draw();
        
        canvasFit->Modified();
        
				TLegend *leg = new TLegend(0.16, 0.825, 0.555, 0.95);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->AddEntry( gJES[0], "JES=0.96", "LEP");
        leg->AddEntry( gJES[1], "JES=1.00", "LEP");
				leg->AddEntry( gJES[2], "JES=1.04", "LEP");
        leg->Draw();
        
				TF1* constFit = new TF1("constFit", "[0]");
				constFit->SetParNames("offset");
				constFit->SetLineColor(kBlack);
				constFit->SetLineWidth(3);
				constFit->SetLineStyle(7);
				mgJES->Fit("constFit", "EM");
				
        statsGlobal->SetX1NDC(0.555);
        statsGlobal->SetY1NDC(0.825);
        statsGlobal->SetX2NDC(0.95);
        statsGlobal->SetY2NDC(0.95);
				
				TString path("plot/"); path += fMethod; path += "/ensembletest/"; path += "fit_JES_"; path += i; path += "_"; path += j; path += ".eps";
        canvasFit->Print(path);
        
        canvasFit->Clear();
        
        TMultiGraph *mgMass = new TMultiGraph();
				mgMass->SetTitle("m_{t} bias;m_{t,gen};m_{t,meas}-m_{t,gen}");
								
				for (int iJES = 0; iJES < 3; iJES++) {
					double genMassMod[9];
					for (int i = 0; i < 9; i++) {
						genMassMod[i] = genMass[i] + 0.2*(iJES-1);
					}
					
        	gMass.push_back(new TGraphErrors(9, genMassMod, hadTopMassBias[iJES], genMassError, hadTopMassError[iJES]));
					
					gMass[iJES]->SetMarkerStyle(marker_[iJES]);
					gMass[iJES]->SetMarkerColor(color_ [iJES]);
					gMass[iJES]->SetLineColor  (color_ [iJES]);
										
					TString sFit("[0]+(x-172.5-"); sFit += 0.2*(iJES-1); sFit += ")*[1]";
					
        	TF1* linearFit = new TF1("linearFit", sFit);
        	linearFit->SetParLimits(1, -1, 1);
					linearFit->SetParNames("offset", "slope");
					linearFit->SetLineColor(color_[iJES]);
        	gMass[iJES]->Fit("linearFit", "EM");
					
					mgMass->Add(gMass[iJES]);

					/*
        	for (int l = 0; l < 2; l++) {
          	fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
          	fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
        	}
					*/
        }
				
				//mgJES->GetXaxis()->SetLimits( 0.95, 1.05);
				//mgJES->GetYaxis()->SetLimits(-0.05, 0.05);
				
				mgMass->SetMinimum(-5);
				mgMass->SetMaximum( 5);
				
				mgMass->Draw("AP");
				
        canvasFit->Update();
        TPaveStats* statsMass0 = (TPaveStats*) gMass[0]->GetListOfFunctions()->FindObject("stats");
        statsMass0->SetTextColor(color_[0]);
        statsMass0->SetX1NDC(0.16);
        statsMass0->SetY1NDC(0.7);
        statsMass0->SetX2NDC(0.4233);
        statsMass0->SetY2NDC(0.825);
        TPaveStats* statsMass1 = (TPaveStats*) gMass[1]->GetListOfFunctions()->FindObject("stats");
        statsMass1->SetTextColor(color_[1]);
        statsMass1->SetX1NDC(0.4233);
        statsMass1->SetY1NDC(0.7);
        statsMass1->SetX2NDC(0.6867);
        statsMass1->SetY2NDC(0.825);
        TPaveStats* statsMassGlobal = (TPaveStats*) gMass[2]->GetListOfFunctions()->FindObject("stats");
        TPaveStats* statsMass2 = (TPaveStats*) statsMassGlobal->Clone();
        statsMass2->SetTextColor(color_[2]);
        statsMass2->SetX1NDC(0.6867);
        statsMass2->SetY1NDC(0.7);
        statsMass2->SetX2NDC(0.95);
        statsMass2->SetY2NDC(0.825);
        statsMass2->Draw();
        canvasFit->Modified();
				
				TLegend *legMass = new TLegend(0.16, 0.825, 0.555, 0.95);
        legMass->SetFillColor(kWhite);
        legMass->SetBorderSize(1);
        legMass->AddEntry( gMass[0], "JES=0.96", "LEP");
        legMass->AddEntry( gMass[1], "JES=1.00", "LEP");
				legMass->AddEntry( gMass[2], "JES=1.04", "LEP");
        legMass->Draw();
				
				mgMass->Fit("constFit", "EM");
				
        statsMassGlobal->SetX1NDC(0.555);
        statsMassGlobal->SetY1NDC(0.825);
        statsMassGlobal->SetX2NDC(0.95);
        statsMassGlobal->SetY2NDC(0.95);
				
				path = "plot/"; path += fMethod; path += "/ensembletest/"; path += "fit_Mass_"; path += i; path += "_"; path += j; path += ".eps";
        canvasFit->Print(path);
      }
    }
  }
  
  /*
  Measure(a1665);
  Measure(a1725);
  Measure(a1785);
  */
}



TH2F* TopMass::Measure(Analysis* a) {
  
  a->Analyze();
  
  TCanvas* canvas = new TCanvas("canvas", "Hadronic top mass", 1000, 500);
  
  TH2F* hMass = a->GetH2Mass();
  TH2F* hMassError = a->GetH2MassError();
  
  TH2F* hJES = a->GetH2JES();
  TH2F* hJESError = a->GetH2JESError();

  TH2F* hMassCalibrated = a->GetH2MassCalibrated();
  TH2F* hMassErrorCalibrated = a->GetH2MassErrorCalibrated();
  
  std::cout << "================================================" << std::endl;
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      
      double mass = hMass->GetCellContent(i+1, j+1);
      double massError = hMassError->GetCellContent(i+1, j+1);
      
      std::cout << "Measured mass: " << mass << " +/- " << massError << " GeV" << std::endl;
      
      double JES = hJES->GetCellContent(i+1, j+1);
      double JESError = hJESError->GetCellContent(i+1, j+1);
      
      std::cout << "Measured JES: " << JES << " +/- " << JESError << " " << std::endl;
    
      if (fCalibFitParameter[i][j][0] && fCalibFitParameter[i][j][1] && mass > 0) {
        massError = sqrt(pow((1-fCalibFitParameter[i][j][1])*massError, 2) + pow(fCalibFitParError[i][j][0], 2) + pow(fCalibFitParError[i][j][1]*(172.5-mass), 2));
        mass = mass - fCalibFitParameter[i][j][0] + (fCalibFitParameter[i][j][1]*(172.5-mass));
        
        std::cout << "Calibrated TopMass: " << mass << " +/- " << massError << " GeV" << std::endl;
      }
      else {
        mass = 0;
        massError = 0;
      }
        
      hMassCalibrated->SetCellContent(i+1, j+1, mass);
      hMassErrorCalibrated->SetCellContent(i+1, j+1, massError);
      
    }
  }
  
  std::cout << "================================================" << std::endl;
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hMassCalibrated->Draw("COLZ,TEXT");
  hMassCalibrated->SetAxisRange(160, 185, "Z");
  
  canvas->cd(2);
  hMassErrorCalibrated->Draw("COLZ,TEXT");
  hMassErrorCalibrated->SetAxisRange(0.05, 5, "Z");
  
  TString path("plot/"); path += fMethod; path += "_"; path += a->GetIdentifier(); path +="_calibrated.eps";
  canvas->Print(path);
  
  //return hMassCalibrated;
  return hMass;
}


void TopMass::QuickSystematics() {

  Analysis* a1725 = new Analysis("1725", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root", fMethod, fBins, 100000);
  Analysis* a1725_jes_up = new Analysis("1725_jes_up", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725-S4_flavor:up/analyzeTop.root", fMethod, fBins, 100000);
  Analysis* a1725_jes_down = new Analysis("1725_jes_down", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725-S4_flavor:down/analyzeTop.root", fMethod, fBins, 100000);
  
  TH2F* hMassJESdown = Measure(a1725_jes_down);
  TH2F* hMassJESnorm = Measure(a1725);
  TH2F* hMassJESup   = Measure(a1725_jes_up);
  
  hMassJESup->Add(hMassJESnorm, -1.000000);
  hMassJESnorm->Add(hMassJESdown, -1.000000);
  
  hMassJESup->SetTitle("JES up MassError");
  hMassJESnorm->SetTitle("JES down MassError");
  
  TCanvas* canvas = new TCanvas("canvas", "Hadronic top mass JES error", 1000, 500);
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hMassJESnorm->Draw("COLZ,TEXT");
  hMassJESnorm->SetAxisRange(0.05, 2, "Z");
  
  //*
  canvas->cd(2);
  hMassJESup->Draw("COLZ,TEXT");
  hMassJESup->SetAxisRange(0.05, 5, "Z");
  //*/
  
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
    TopMass* top = new TopMass(argv[1], 1, 4700);
  }
  else {
    TopMass* top = new TopMass("GenMatch", 6, 500);
  }
}
