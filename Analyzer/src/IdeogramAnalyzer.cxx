#include "IdeogramAnalyzer.h"

double IdeogramAnalyzer::GetMass() {
  return fMass;
}

void IdeogramAnalyzer::Analyze(TString cuts, int i, int j) {

  bool debug = false;
  int nDebug = 1;
  int maxDebug = 1000;
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int firstbin = 150;
  int lastbin  = 200;

  int bins = 100;
  
  IdeogramCombLikelihood* fptr = new IdeogramCombLikelihood();
  TF1* combLikelihood = new TF1("combLikelihood",fptr,&IdeogramCombLikelihood::Evaluate,150,200,4);
  //TF1* combBackground = new TF1("combBackground",fptr,&IdeogramCombLikelihood::CrystalBall,150,200,1);

  TF1* fitParabola = new TF1("fitParabola", "abs([1])*(x-[0])^2+[2]");
  fitParabola->SetParLimits(0, firstbin, lastbin);
  fitParabola->SetParLimits(1, 0.001, 100);
  fitParabola->SetLineColor(kRed+1);

  TF1* null = new TF1("null", "0");
  TF1* unity = new TF1("unity", "1");
  TH1D* hUnity = new TH1D("hUnity","hUnity", bins, firstbin, lastbin);
  hUnity->Eval(unity);

  TH1D* eventLikelihood = new TH1D("eventLikelihood","eventLikelihood", bins, firstbin, lastbin);
  TH1D* logEventLikelihood = new TH1D("logEventLikelihood", "logEventLikelihood", bins, firstbin, lastbin);
  TH1D* sumLogLikelihood = new TH1D("sumLogLikelihood", "sumLogLikelihood", bins, firstbin, lastbin);
  sumLogLikelihood->Eval(null);

  double hadTopMass, fitChi2, fitProb, bProb, hadBProb, bProbSSV, weight, currentWeight;
  double hitFitProb, hitFitMT, hitFitSigMT;
  double hadQBTCHE, hadQbarBTCHE, hadBBTCHE, lepBBTCHE;
  int event, currentEvent;
  int combi, previousCombi = -1;
  int nEvents = 0;
  
  TTree* eventTree = fTree->CopyTree(cuts);
  
  eventTree->SetBranchAddress("hadTopMass", &hadTopMass);
  eventTree->SetBranchAddress("hadQBTCHE", &hadQBTCHE);
  eventTree->SetBranchAddress("hadQbarBTCHE", &hadQbarBTCHE);
  eventTree->SetBranchAddress("hadBBTCHE", &hadBBTCHE);
  eventTree->SetBranchAddress("lepBBTCHE", &lepBBTCHE);
  eventTree->SetBranchAddress("fitChi2", &fitChi2);
  eventTree->SetBranchAddress("fitProb", &fitProb);
  eventTree->SetBranchAddress("hitFitProb", &hitFitProb);
  eventTree->SetBranchAddress("hitFitMT", &hitFitMT);
  eventTree->SetBranchAddress("hitFitSigMT", &hitFitSigMT);
  eventTree->SetBranchAddress("bProb", &bProb);
  eventTree->SetBranchAddress("hadBProb", &hadBProb);
  eventTree->SetBranchAddress("bProbSSV", &bProbSSV);
  eventTree->SetBranchAddress("event", &event);
  eventTree->SetBranchAddress("combi", &combi);
  
  // Build Likelihood
  for (int iEntry = 0; iEntry < eventTree->GetEntries(); iEntry++) {
    eventTree->GetEntry(iEntry);

    if (event == currentEvent) continue;
    currentEvent = event;
    nEvents++;
    if (debug && iEntry%nDebug == 0 && iEntry < 100) std::cout << iEntry << " - " << event << std::endl;
    
    eventLikelihood->Eval(null);
    eventLikelihood->SetFillColor(0);
    weight = 0;
    currentWeight = 0;

    for (int iComb = 0; iComb < 24; iComb++) {
      if (eventTree->GetEntries() < iEntry + iComb + 1) break;
      eventTree->GetEntry(iEntry + iComb);
      
      if (event != currentEvent) break;
      
      if (debug && iEntry%nDebug == 0 && iEntry < 100) {
        std::cout << "Combi: " << combi << "\tMass: " << hadTopMass
                  << "\thitFitProb: " << hitFitProb
                  << "\tbProb: " << bProb << std::endl;
      }
      
      //if (bProb * fitProb < 1e-3) continue;
      currentWeight = bProb * hitFitProb;
      if (currentWeight > weight) weight = currentWeight;
      if (currentWeight != 0) {
        /*
        combBackground->SetParameter(0, hadTopMass);
        double bkgIntegral = combBackground->Integral(0, 10000);
        //*/
        
        combLikelihood->SetParameters(hadTopMass, hitFitSigMT, currentWeight);
        eventLikelihood->Eval(combLikelihood, "A"); // add combi pdf
      }
    }
    
    if (weight == 0) continue;
    
    logEventLikelihood->Eval(null);

    for (int i = 0; i<=bins; i++) {
  	  logEventLikelihood->SetBinContent(i, -2*TMath::Log(eventLikelihood->GetBinContent(i)));
    }

    sumLogLikelihood->Add(logEventLikelihood, weight); // add weight here
    
    if (debug && iEntry%nDebug == 0 && iEntry < 100) {
      TCanvas* eventCanvas = new TCanvas("eventCanvas", "eventCanvas", 1200, 400);
      eventCanvas->Divide(3, 1);
      
      eventCanvas->cd(1);
      eventLikelihood->Draw();
      if (weight > 1./16) eventLikelihood->SetFillColor(kGreen);
      if (weight > 1./8) eventLikelihood->SetFillColor(kYellow);
      if (weight > 1./4) eventLikelihood->SetFillColor(kOrange);
      if (weight > 1./2) eventLikelihood->SetFillColor(kRed);
      
      eventCanvas->cd(2);
      logEventLikelihood->Draw();
      
      eventCanvas->cd(3);
      sumLogLikelihood->Draw();
      
      TString eventPath("plot/Ideogram/"); eventPath += fIdentifier; eventPath += "_"; eventPath += iEntry; eventPath += "_"; eventPath += currentEvent; eventPath += ".png";
      eventCanvas->Print(eventPath);
      
      delete eventCanvas;
    }
  }
  
  ctemp->cd();
  
  sumLogLikelihood->Add(hUnity, -sumLogLikelihood->GetMinimum(0));
  //sumLogLikelihood->SetAxisRange(0, 100, "Y");
  sumLogLikelihood->SetAxisRange(sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) - 10, sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) + 10, "X");
  sumLogLikelihood->Draw();
  
  std::cout << "Minimum likelihood: " << sumLogLikelihood->GetMinimum(0) << "\tMaximum likelihood (in range): " << sumLogLikelihood->GetMaximum() << std::endl;
  std::cout << "Total number of events: " << nEvents << std::endl;
  
  fitParabola->SetParameter(2, sumLogLikelihood->GetMinimum(0));
  fitParabola->SetParameter(1, 100);
  
  fitParabola->SetRange(sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) - 3, sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) + 3);

  sumLogLikelihood->Fit("fitParabola","WEMR");
  
  if (firstbin+1 < fitParabola->GetParameter(0) && fitParabola->GetParameter(0) < lastbin-1) {
    fMass = fitParabola->GetParameter(0);
    if (TMath::Sqrt(1/fitParabola->GetParameter(1)) < 2*TMath::Sqrt(fMass)) {
      fMassError = TMath::Sqrt(1/fitParabola->GetParameter(1));
    }
    else fMassError = -1;
  }
  else {
    fMass = -1;
    fMassError = -1;
  }
  fMassSigma = -1;
  
  TString path("plot/Ideogram/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
  
  delete fitParabola;
  delete null;
  delete eventTree;
  delete ctemp;
  delete eventLikelihood;
  delete logEventLikelihood;
  delete sumLogLikelihood;
}

double IdeogramAnalyzer::QBTagProbability(double bDiscriminator) {
  if (bDiscriminator == -100) return 0.787115;
  if (bDiscriminator < 0) return 1;
  
  double p0 = 5.91566e+00;
  double p1 = 5.94611e-01;
  double p2 = 3.53592e+00;
  
  return p0 * TMath::Voigt(bDiscriminator, p1, p2);
}
