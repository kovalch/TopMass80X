#include "IdeogramAnalyzer.h"

double IdeogramAnalyzer::GetMass() {
  return fMass;
}

void IdeogramAnalyzer::Analyze(TString cuts, int i, int j) {

  bool debug = false;
  int nDebug = 1;
  int maxDebug = 200;
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  double firstBinMass = 154;
  double lastBinMass  = 190;
  double resolMass    = 1.;
  int binsMass     = (lastBinMass-firstBinMass)/resolMass;
  
  double firstBinJes = 0.9;
  double lastBinJes  = 1.1;
  double resolJes    = 0.005;
  int binsJes        = (lastBinJes-firstBinJes)/resolJes;
  
  double pullWidth   = 1.365*1.07; //*1.053; //1.37068*1.07696;//1.24956e+00; //0.715;
  
  // Pile up corrections
  double mWnVertex   = 0.; // IdeogramAnalyzer::mWnVertex() - 0.1;
  
  //*
  if (debug) {
    firstBinMass = 100;
    lastBinMass  = 350;
    binsMass     = 125;
    
    firstBinJes = 0.7;
    lastBinJes  = 1.3;
    binsJes     = 60;
  }
  //*/
  
  IdeogramCombLikelihood* fptr = new IdeogramCombLikelihood();
  TF2* combLikelihood = new TF2("combLikelihood",fptr,&IdeogramCombLikelihood::Evaluate, firstBinMass, lastBinMass, firstBinJes, lastBinJes, 4, "IdeogramCombLikelihood", "Evaluate");
  //TF1* combBackground = new TF1("combBackground",fptr,&IdeogramCombLikelihood::CrystalBall,150,200,1);

  TF2* fitParabola = new TF2("fitParabola", "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2");
  TF2* systParabola = new TF2("systParabola", "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2");
  fitParabola->SetNpx(300);
  fitParabola->SetNpy(300);
  fitParabola->SetParNames("mass", "massCurv", "jes", "jesCurv", "alpha");
  fitParabola->SetParLimits(0, firstBinMass, lastBinMass);
  fitParabola->SetParLimits(1, 0.0001, 1000);
  fitParabola->SetParLimits(2, 0.5, 1.5);
  fitParabola->SetParLimits(3, 1, 10000000);
  fitParabola->SetParLimits(4, 0, 6.3);
  fitParabola->SetParameter(4, 1.5);
  fitParabola->SetLineColor(kWhite);
  fitParabola->SetLineWidth(3);

  TF2* null = new TF2("null", "0 + 0*x + 0*y");
  TF2* unity = new TF2("unity", "1 + 0*x + 0*y");
  TH2D* hUnity = new TH2D("hUnity","hUnity", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  hUnity->Eval(unity);

  TH2D* eventLikelihood = new TH2D("eventLikelihood","eventLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* logEventLikelihood = new TH2D("logEventLikelihood", "logEventLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* sumLogLikelihood = new TH2D("sumLogLikelihood", "sumLogLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  sumLogLikelihood->Eval(null);
  
  eventLikelihood->SetTitle("L(m_{t}|event)");
  eventLikelihood->SetXTitle("m_{t}");
  eventLikelihood->SetYTitle("JES");
  logEventLikelihood->SetTitle("-2#upointln{L(m_{t}|event)}");
  logEventLikelihood->SetXTitle("m_{t}");
  logEventLikelihood->SetYTitle("JES");
  sumLogLikelihood->SetTitle("-2#upointln{L(m_{t}|sample)}");
  sumLogLikelihood->SetXTitle("m_{t}");
  sumLogLikelihood->SetYTitle("JES");

  double hadTopMass, hadTopPt, lepTopPt, hadWPt, lepWPt, hadBPt, lepBPt, hadWRawMass, topPtAsymmetry, bScaleEstimator;
  double hadWE, deltaThetaHadWHadB, sinThetaStar;
  double hitFitChi2, hitFitProb, MCWeight, PUWeight, muWeight, bWeight, bWeight_bTagSFUp, bWeight_bTagSFDown, bWeight_misTagSFUp, bWeight_misTagSFDown, weight, currentWeight;
  double pdfWeights[44];
  int event, currentEvent, nVertex;
  int combi;
  int nEvents = 0;
  
  TTree* eventTree = fTree->CopyTree(cuts);
  
  eventTree->SetBranchAddress("hadTopMass", &hadTopMass);
  eventTree->SetBranchAddress("hadWRawMass", &hadWRawMass);
  eventTree->SetBranchAddress("deltaThetaHadWHadB", &deltaThetaHadWHadB);
  eventTree->SetBranchAddress("hitFitChi2", &hitFitChi2);
  eventTree->SetBranchAddress("hitFitProb", &hitFitProb);
  eventTree->SetBranchAddress("event", &event);
  eventTree->SetBranchAddress("combi", &combi);
  eventTree->SetBranchAddress("PUWeight", &PUWeight);
  eventTree->SetBranchAddress("muWeight", &muWeight);
  eventTree->SetBranchAddress("bWeight", &bWeight);
  eventTree->SetBranchAddress("bWeight_bTagSFUp", &bWeight_bTagSFUp);
  eventTree->SetBranchAddress("bWeight_bTagSFDown", &bWeight_bTagSFDown);
  eventTree->SetBranchAddress("bWeight_misTagSFUp", &bWeight_misTagSFUp);
  eventTree->SetBranchAddress("bWeight_misTagSFDown", &bWeight_misTagSFDown);
  eventTree->SetBranchAddress("nVertex", &nVertex);
  eventTree->SetBranchAddress("pdfWeights", &pdfWeights);
  
  // Build Likelihood
  for (int iEntry = 0; iEntry < eventTree->GetEntries(); iEntry++) {
    eventTree->GetEntry(iEntry);

    if (event == currentEvent) continue;
    currentEvent = event;
    nEvents++;
    if ((debug && iEntry%nDebug == 0 && iEntry < maxDebug) || iEntry%1000 == 0) std::cout << iEntry << " - " << event << std::endl;
    
    eventLikelihood->Eval(null);
    eventLikelihood->SetFillColor(0);
    weight = 0;
    currentWeight = 0;
    
    if (debug && iEntry%nDebug == 0 && iEntry < maxDebug) {
    std::cout << std::setiosflags(std::ios::left)
              << std::setw(04) << "i"
              << std::setw(10) << "mt"
              << std::setw(10) << "mW"
              << std::setw(12) << "fitProb"
              << std::setw(11) << "bProb"
              << std::setw(11) << "weight"
              << std::endl;
    }
    
    for (int iComb = 0; iComb < 24; iComb++) {
      if (eventTree->GetEntries() < iEntry + iComb + 1) break;
      eventTree->GetEntry(iEntry + iComb);
      
      if (event != currentEvent) break;
      
      MCWeight = (PUWeight != -100.) ? PUWeight * muWeight * bWeight : 1;
      //std::cout << MCWeight << std::endl;
      currentWeight = hitFitProb * MCWeight;
      //currentWeight = hitFitProb * PUWeight * bWeight * muWeight;
      weight += currentWeight;
      
      if (debug && iEntry%nDebug == 0 && iEntry < maxDebug) {
        std::cout << std::setw(04) << combi
                  << std::setw(10) << hadTopMass
                  << std::setw(10) << hadWRawMass
                  << std::setw(12) << hitFitProb
                  << std::setw(11) << currentWeight
                  << std::endl;
      }
      
      if (currentWeight != 0) {
        /*
        topPtAsymmetry  = (hadTopPt-lepTopPt)/(hadTopPt+lepTopPt);
        bScaleEstimator = -(hadWPt-lepWPt)/(hadBPt-lepBPt);
        sinThetaStar    = sqrt( pow(sin(deltaThetaHadWHadB),2)*pow(172.5,2) / pow(1/2*(pow(172.5,2)-pow(80.4,2))/(hadWE-sqrt(pow(hadWE,2) - pow(80.4,2))) - (hadWE-sqrt(pow(hadWE,2) - pow(80.4,2))*cos(deltaThetaHadWHadB))*cos(deltaThetaHadWHadB),2) + 1/pow(((hadWE-(sqrt(pow(hadWE,2) - pow(80.4,2)))*cos(deltaThetaHadWHadB)))/172.5,2));
        //std::cout << sinThetaStar << std::endl;
        */
        bScaleEstimator = 1;
        
        // Set Likelihood parameters
        combLikelihood->SetParameters(currentWeight, hadTopMass, hadWRawMass, bScaleEstimator);
        
        // add permutation to event likelihood
        eventLikelihood->Eval(combLikelihood, "A");
      }
    }
    
    if (weight == 0) continue;
    
    logEventLikelihood->Eval(null);

    for (int i = 0; i<=binsMass; i++) {
      for (int j = 0; j<=binsJes; j++) {
    	  logEventLikelihood->SetBinContent(i, j, -2*TMath::Log(eventLikelihood->GetBinContent(i, j)));
    	}
    }

    sumLogLikelihood->Add(logEventLikelihood, weight/(pullWidth*pullWidth)); // add weight here
    
    if (debug && iEntry%nDebug == 0 && iEntry < maxDebug) {
      TCanvas* eventCanvas = new TCanvas("eventCanvas", "eventCanvas", 1200, 400);
      eventCanvas->Divide(3, 1);
      
      eventCanvas->cd(1);
      eventLikelihood->SetEntries(1);
      eventLikelihood->Draw("COLZ");
      eventLikelihood->SetFillColor(kRed+1);
      /*
      if (weight > 1./16) eventLikelihood->SetFillColor(kGreen);
      if (weight > 1./8) eventLikelihood->SetFillColor(kYellow);
      if (weight > 1./4) eventLikelihood->SetFillColor(kOrange);
      if (weight > 1./2) eventLikelihood->SetFillColor(kRed);
      */
      
      eventCanvas->cd(2);
      eventLikelihood->SetEntries(1);
      logEventLikelihood->Draw("COLZ");
      
      eventCanvas->cd(3);
      eventLikelihood->SetEntries(1);
      sumLogLikelihood->Draw("COLZ");
      
      TString eventPath("plot/Ideogram/"); eventPath += fIdentifier; eventPath += "_"; eventPath += iEntry; eventPath += "_"; eventPath += currentEvent; eventPath += ".eps";
      eventCanvas->Print(eventPath);
      
      delete eventCanvas;
    }
  }
  
  ctemp->cd();
  
  std::cout << "Finishing..." << std::endl;
  
  sumLogLikelihood->Add(hUnity, -sumLogLikelihood->GetMinimum(0) + 1e-2);
  //sumLogLikelihood->SetAxisRange(0, 100, "Y");
  
  int minBinX;
  int minBinY;
  int minBinZ;
  
  sumLogLikelihood->GetMinimumBin(minBinX, minBinY, minBinZ);
  
  double minMass = sumLogLikelihood->GetXaxis()->GetBinCenter(minBinX);
  double minJes  = sumLogLikelihood->GetYaxis()->GetBinCenter(minBinY);
  
  std::cout << "minMass: " << minMass << std::endl;
  std::cout << "minJes: " << minJes << std::endl;
  
  /*
  sumLogLikelihood->SetAxisRange(minMass - 3, minMass + 3, "X");
  sumLogLikelihood->SetAxisRange(minJes - 0.03, minJes + 0.03, "Y");
  sumLogLikelihood->SetAxisRange(0, 20, "Z");
  //*/
  
  sumLogLikelihood->SetAxisRange(0, 100, "Z");
  
  //sumLogLikelihood->SetMarkerStyle(20);
  //sumLogLikelihood->SetMarkerColor(kRed+1);
  sumLogLikelihood->Draw("COLZ");
  
  std::cout << "Minimum likelihood: " << sumLogLikelihood->GetMinimum(0) << "\tMaximum likelihood (in range): " << sumLogLikelihood->GetMaximum() << std::endl;
  std::cout << "Total number of events: " << nEvents << std::endl;
  
  fitParabola->SetParLimits(0, minMass-2*resolMass, minMass+2*resolMass);
  fitParabola->SetParameter(0, minMass);
  fitParabola->SetParameter(2, minJes);
  fitParabola->SetParameter(3, 1000000);
  
  //fitParabola->SetRange(minMass - 1, minJes - 0.01, minMass + 1, minJes + 0.01);
  fitParabola->SetRange(minMass - 4, minJes - 0.04, minMass + 4, minJes + 0.04);
  //fitParabola->SetRange(minMass - 20, minJes - 0.2, minMass + 20, minJes + 0.2);

  sumLogLikelihood->Fit("fitParabola","WEMR0");
  
  double semiMajor, semiMinor, alpha;
  
  if (firstBinMass+1 < fitParabola->GetParameter(0) && fitParabola->GetParameter(0) < lastBinMass-1) {
    fMass = fitParabola->GetParameter(0);
    fJES  = fitParabola->GetParameter(2);
    if (TMath::Sqrt(1/fitParabola->GetParameter(1)) < 2*TMath::Sqrt(fMass)) {
      semiMajor = TMath::Sqrt(1/fitParabola->GetParameter(1));
      semiMinor = TMath::Sqrt(1/fitParabola->GetParameter(3));
      alpha     = fitParabola->GetParameter(4);
      
      fMassError = sqrt(pow(semiMajor * cos(alpha), 2) + pow(semiMinor * sin(alpha), 2));
      fJESError  = sqrt(pow(semiMajor * sin(alpha), 2) + pow(semiMinor * cos(alpha), 2));
    }
    else fMassError = -1;
  }
  else {
    fMass = -1;
    fMassError = -1;
  }
  fMassSigma = -1;
  
  // Fit again with previous result as range
  double sigmaLevel = 4;
  if (3*fMassError < resolMass) sigmaLevel = 8;
	fitParabola->SetRange(fMass - sigmaLevel*fMassError, fJES - sigmaLevel*fJESError,
	                      minMass + sigmaLevel*fMassError, fJES + sigmaLevel*fJESError);

  sumLogLikelihood->Fit("fitParabola","WEMR0");
  
  double contours[3] = {1, 4, 9};
  fitParabola->SetContour(3, contours);
  fitParabola->Draw("cont3 same");
  
  if (firstBinMass+1 < fitParabola->GetParameter(0) && fitParabola->GetParameter(0) < lastBinMass-1) {
    fMass = fitParabola->GetParameter(0);
    fJES  = fitParabola->GetParameter(2);
    if (TMath::Sqrt(1/fitParabola->GetParameter(1)) < 2*TMath::Sqrt(fMass)) {
      semiMajor = TMath::Sqrt(1/fitParabola->GetParameter(1));
      semiMinor = TMath::Sqrt(1/fitParabola->GetParameter(3));
      alpha     = fitParabola->GetParameter(4);
      
      fMassError = sqrt(pow(semiMajor * cos(alpha), 2) + pow(semiMinor * sin(alpha), 2));
      fJESError  = sqrt(pow(semiMajor * sin(alpha), 2) + pow(semiMinor * cos(alpha), 2));
    }
    else fMassError = -1;
  }
  else {
    fMass = -1;
    fMassError = -1;
  }
  fMassSigma = -1;
  
  // stat+syst ellipsis  
  double mSyst = 1.07;
  double jSyst = 0.009;
  
  double sm2 = fMassError*fMassError + mSyst*mSyst;
  double sj2 = fJESError*fJESError + jSyst*jSyst;
  
  fitParabola->Copy(*systParabola);
  systParabola->SetRange(fMass - 40, fJES - 0.4, minMass + 40, fJES + 0.4);
  systParabola->SetParameter(1, 2*cos(2*alpha)/(sm2 - sj2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
  systParabola->SetParameter(3, 2*cos(2*alpha)/(sj2 - sm2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
  systParabola->SetLineColor(kBlack);
  systParabola->SetLineWidth(5);
  systParabola->Draw("cont3 same");
  fitParabola->Draw("cont3 same");
  
  // create legend
  TLegend *leg0 = new TLegend(0.2, 0.15, 0.45, 0.25);
  leg0->SetFillStyle(0);
  leg0->SetBorderSize(0);
  //leg0->AddEntry((TObject*)0, "1, 2, 3#sigma", "");
  leg0->AddEntry(fitParabola, "stat", "L");
  leg0->AddEntry(systParabola, "stat + syst", "L");
  leg0->Draw();
  
  std::cout << "fMassError: " << fMassError << std::endl;
  
  TString path("plot/Ideogram/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".eps";
  ctemp->Print(path);
  
  delete fitParabola;
  delete null;
  delete eventTree;
  delete ctemp;
  delete eventLikelihood;
  delete logEventLikelihood;
  delete sumLogLikelihood;
  
  std::cout << "IdeogramAnalyzer done" << std::endl;
}

double IdeogramAnalyzer::QBTagProbability(double bDiscriminator) {
  if (bDiscriminator == -100) return 0.787115;
  if (bDiscriminator < 0) return 1;
  
  double p0 = 5.91566e+00;
  double p1 = 5.94611e-01;
  double p2 = 3.53592e+00;
  
  return p0 * TMath::Voigt(bDiscriminator, p1, p2);
}

double IdeogramAnalyzer::mWnVertex() {
  TF1* linearFit = new TF1("linearFit", "[0]+(x-0.)*[1]");
  
  fTree->Draw("hadWRawMass:nVertex >> h2(15, 0, 15, 100, 50, 150)", "(hitFitProb>0.2)");
  TH2D* h2 = (TH2D*)gDirectory->Get("h2");
  h2->FitSlicesY();
  TH1D *h2_1 = (TH1D*)gDirectory->Get("h2_1");
  h2_1->Fit("linearFit");

  return linearFit->GetParameter(1);
}


