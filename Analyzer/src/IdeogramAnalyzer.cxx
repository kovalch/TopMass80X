#include "IdeogramAnalyzer.h"

double IdeogramAnalyzer::GetMass() {
  return fMass;
}

void IdeogramAnalyzer::Analyze(TString cuts, int i, int j) {
  Scan(cuts, i, j, 154, 190, 2, 0.9, 1.1, 0.02);
  Scan(cuts, i, j, fMass-2, fMass+2, 0.25, fJES-0.015, fJES+0.015, 0.0015);
  //Scan(cuts, i, j, fMass-2, fMass+2, 0.1, fJES-0.015, fJES+0.015, 0.00075);
  double epsilon = 1e-6;
  //Scan(cuts, i, j, 154, 190, 2, 0.99, 1.01, 0.01);
  Scan(cuts, i, j, fMass-5, fMass+5, 0.25, 1.-epsilon, 1.+epsilon, epsilon, false);
  
  //Scan(cuts, i, j, fMass-2, fMass+2, 0.1, fJES-0.015, fJES+0.015, 0.0005);
  //Scan(cuts, i, j, 170, 176, 0.1, 0.999, 1.001, 0.0001);
}

void IdeogramAnalyzer::Scan(TString cuts, int i, int j, double firstBinMass, double lastBinMass,
              double resolMass, double firstBinJes, double lastBinJes, double resolJes, bool fit2D) {
  //*
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //*/
  
  /*
  Double_t r[]    = {0., 0.0, 1.0, 1.0, 1.0};
  Double_t g[]    = {0., 0.0, 0.0, 1.0, 1.0};
  Double_t b[]    = {0., 1.0, 0.0, 0.0, 1.0};
  Double_t stop[] = {0., .25, .50, .75, 1.0};
  TColor::CreateGradientColorTable(5, stop, r, g, b, 100);
  //*/
  
  bool blackWhite = false;
  bool syst       = false;
  bool cmsPrel    = true;
  bool useWeight  = false;
  
  bool debug = false;
  int nDebug = 1;
  int maxDebug = 200;
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int binsMass     = (lastBinMass-firstBinMass)/resolMass;
  
  int binsJes        = (lastBinJes-firstBinJes)/resolJes;
  
  double pullWidth   = 1.;//1.02;//06; //1.46;
  
  /*
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
  TF2* gausJESConstraint = new TF2("gausJESConstraint", "x*0 + ((y-2+0.995811)/0.016)**2", firstBinMass, lastBinMass, firstBinJes, lastBinJes);
  //TF1* combBackground = new TF1("combBackground",fptr,&IdeogramCombLikelihood::CrystalBall,150,200,1);

  TF1* fitParabola = new TF1("fitParabola", "abs([1])*(x-[0])^2+[2]");
  fitParabola->SetParLimits(0, firstBinMass, lastBinMass);
  fitParabola->SetParLimits(1, 0.0001, 1000);
  fitParabola->SetLineColor(kRed+1);
  
  TF2* fitParaboloid = new TF2("fitParaboloid", "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2 + [5]");
  TF2* systParaboloid = new TF2("systParaboloid", "abs([1])*((x-[0])*cos([4])-(y-[2])*sin([4]))^2 + abs([3])*((x-[0])*sin([4])+(y-[2])*cos([4]))^2 + [5]");
  fitParaboloid->SetNpx(300);
  fitParaboloid->SetNpy(300);
  fitParaboloid->SetParNames("mass", "massCurv", "jes", "jesCurv", "alpha", "minL");
  fitParaboloid->SetParLimits(0, firstBinMass, lastBinMass);
  fitParaboloid->SetParLimits(1, 0.0001, 1000);
  fitParaboloid->SetParLimits(2, 0.5, 1.5);
  fitParaboloid->SetParLimits(3, 1, 10000000);
  fitParaboloid->SetParLimits(4, 0, 6.3);
  fitParaboloid->SetParameter(4, 1.5);
  if (blackWhite) fitParaboloid->SetLineColor(kRed+1);
  else fitParaboloid->SetLineColor(kWhite);
  fitParaboloid->SetLineWidth(3);

  TF2* null = new TF2("null", "0 + 0*x + 0*y");
  TF2* unity = new TF2("unity", "1 + 0*x + 0*y");
  TH2D* hUnity = new TH2D("hUnity","hUnity", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  hUnity->Eval(unity);

  TH2D* eventLikelihood = new TH2D("eventLikelihood","eventLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* logEventLikelihood = new TH2D("logEventLikelihood", "logEventLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  TH2D* sumLogLikelihood = new TH2D("sumLogLikelihood", "sumLogLikelihood", binsMass, firstBinMass, lastBinMass, binsJes, firstBinJes, lastBinJes);
  sumLogLikelihood->Eval(null);
  
  eventLikelihood->SetTitle("L(m_{t}|event)");
  eventLikelihood->SetXTitle("m_{t} [GeV]");
  eventLikelihood->SetYTitle("JES");
  logEventLikelihood->SetTitle("-2#upointln{L(m_{t}|event)}");
  logEventLikelihood->SetXTitle("m_{t} [GeV]");
  logEventLikelihood->SetYTitle("JES");
  sumLogLikelihood->SetTitle("-2#upointln{L(m_{t}|sample)}");
  sumLogLikelihood->SetXTitle("m_{t} [GeV]");
  sumLogLikelihood->SetYTitle("JES");

  double hadTopMass, hadTopPt, lepTopPt, hadWPt, lepWPt, hadBPt, lepBPt, hadWRawMass, topPtAsymmetry, bScaleEstimator;
  double hadWE, deltaThetaHadWHadB, sinThetaStar;
  double hitFitChi2, hitFitProb, MCWeight, PUWeight, muWeight, bWeight, bWeight_bTagSFUp, sumMCWeight, meanMCWeight, bWeight_bTagSFDown, bWeight_misTagSFUp, bWeight_misTagSFDown, weight, currentWeight, fitWeight;
  double eventWeight;
  double pdfWeights[44];
  int event, currentEvent, nVertex, run, luminosityBlock, leptonId;
  int combi;
  int nEvents = 0;
  double productWeights = 1.;
  double sumWeights = 0.;
  double mcWeight = 1;
  
  //TFile* file = new TFile("tree.root", "UPDATE");
  TTree* eventTree = fTree->CopyTree(cuts);
  //file->Write();
  
  eventTree->SetBranchAddress("hadTopMass", &hadTopMass);
  eventTree->SetBranchAddress("hadWRawMass", &hadWRawMass);
  //eventTree->SetBranchAddress("deltaThetaHadWHadB", &deltaThetaHadWHadB);
  //eventTree->SetBranchAddress("hitFitChi2", &hitFitChi2);
  eventTree->SetBranchAddress("hitFitProb", &hitFitProb);
  
  eventTree->SetBranchAddress("run", &run);
  eventTree->SetBranchAddress("luminosityBlock", &luminosityBlock);
  eventTree->SetBranchAddress("event", &event);
  eventTree->SetBranchAddress("combi", &combi);
  eventTree->SetBranchAddress("PUWeight", &PUWeight);
  eventTree->SetBranchAddress("mcWeight", &mcWeight);
  eventTree->SetBranchAddress("muWeight", &muWeight);
  eventTree->SetBranchAddress("bWeight", &bWeight);
  eventTree->SetBranchAddress("leptonId", &leptonId);
  /*
  eventTree->SetBranchAddress("bWeight_bTagSFUp", &bWeight_bTagSFUp);
  eventTree->SetBranchAddress("bWeight_bTagSFDown", &bWeight_bTagSFDown);
  eventTree->SetBranchAddress("bWeight_misTagSFUp", &bWeight_misTagSFUp);
  eventTree->SetBranchAddress("bWeight_misTagSFDown", &bWeight_misTagSFDown);
  eventTree->SetBranchAddress("nVertex", &nVertex);
  eventTree->SetBranchAddress("pdfWeights", &pdfWeights);
  //*/
  
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
    fitWeight = 0;
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
      fitWeight += hitFitProb;
      
      if (debug && iEntry%nDebug == 0 && iEntry < maxDebug) {
        std::cout << std::setw(04) << combi
                  << std::setw(10) << hadTopMass
                  << std::setw(10) << hadWRawMass
                  << std::setw(12) << hitFitProb
                  << std::setw(11) << currentWeight
                  << std::endl;
      }
      
      if (hitFitProb != 0) {
        // Set Likelihood parameters
        combLikelihood->SetParameters(hitFitProb, hadTopMass, hadWRawMass, leptonId);
        
        // add permutation to event likelihood
        eventLikelihood->Eval(combLikelihood, "A");
      }
    }
    
    eventLikelihood->Scale(1./fitWeight);
    sumWeights += fitWeight;
    //if (weight == 0) continue;
    
    logEventLikelihood->Eval(null);

    for (int i = 0; i<=binsMass; i++) {
      for (int j = 0; j<=binsJes; j++) {
    	  logEventLikelihood->SetBinContent(i, j, -2*TMath::Log(eventLikelihood->GetBinContent(i, j)));
    	}
    }
    
    TString sEvent("(run=="); sEvent += run; sEvent += " & luminosityBlock=="; 
    sEvent += luminosityBlock; sEvent += " & event=="; sEvent += event; sEvent += ")";
    
    TString sEventWeighted = sEvent; sEventWeighted += "*("; sEventWeighted += "1"; sEventWeighted += ")";
    
    //double eventWeight = eventTree->GetEntries(sEventWeighted)/eventTree->GetEntries(sEvent);
    
    useWeight ? MCWeight = PUWeight*muWeight*bWeight : MCWeight = 1;
    
    //*
    switch(leptonId) {
      case 11:
        pullWidth = 1.04847e+00;
        break;
      case 13:
        pullWidth = 1.00543e+00;
        break;
      default:
        pullWidth = 1.;
    }
    //*/
    
    //std::cout << "mcWeight: " << mcWeight << std::endl;
    
    sumLogLikelihood->Add(logEventLikelihood, fitWeight*MCWeight/(pullWidth*pullWidth) * mcWeight/fabs(mcWeight)); // add weight here
    
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
  
  std::cout << "Sum of weights: " << sumWeights << std::endl;
  std::cout << "Total number of events: " << nEvents << std::endl;
  
  sumLogLikelihood->Scale(nEvents/sumWeights);
  
  /* JES constraint
  eventLikelihood->Eval(gausJESConstraint);
  for (int i = 0; i<=binsMass; i++) {
    for (int j = 0; j<=binsJes; j++) {
      sumLogLikelihood->SetBinContent(i, j, sumLogLikelihood->GetBinContent(i, j) + eventLikelihood->GetBinContent(i, j));
    }
  }
  //*/
  
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
  
  //sumLogLikelihood->SetMarkerStyle(20);
  //sumLogLikelihood->SetMarkerColor(kRed+1);
  
  std::cout << "Minimum likelihood: " << sumLogLikelihood->GetMinimum(0) << "\tMaximum likelihood (in range): " << sumLogLikelihood->GetMaximum() << std::endl;
  
  Helper* helper = new Helper(1);
  
  if (fit2D) {  
  
    fitParaboloid->SetParLimits(0, minMass-2*resolMass, minMass+2*resolMass);
    fitParaboloid->SetParameter(0, minMass);
    fitParaboloid->SetParLimits(2, minJes-2*resolJes, minJes+2*resolJes);
    fitParaboloid->SetParameter(2, minJes);
    fitParaboloid->SetParameter(3, 1000000);
    fitParaboloid->SetParLimits(5, sumLogLikelihood->GetMinimum(0)-1., sumLogLikelihood->GetMinimum(0)+1.);
    fitParaboloid->SetParameter(5, sumLogLikelihood->GetMinimum(0));
    
    //fitParaboloid->SetRange(minMass - 1, minJes - 0.01, minMass + 1, minJes + 0.01);
    fitParaboloid->SetRange(minMass - 4, minJes - 0.04, minMass + 4, minJes + 0.04);
    //fitParaboloid->SetRange(minMass - 20, minJes - 0.2, minMass + 20, minJes + 0.2);

    sumLogLikelihood->Fit("fitParaboloid","EMR0");
    
    double semiMajor, semiMinor, alpha;
  
    if (firstBinMass+1 < fitParaboloid->GetParameter(0) && fitParaboloid->GetParameter(0) < lastBinMass-1) {
      fMass = fitParaboloid->GetParameter(0);
      fJES  = fitParaboloid->GetParameter(2);
      if (TMath::Sqrt(1/fitParaboloid->GetParameter(1)) < 2*TMath::Sqrt(fMass)) {
        semiMajor = TMath::Sqrt(1/fitParaboloid->GetParameter(1));
        semiMinor = TMath::Sqrt(1/fitParaboloid->GetParameter(3));
        alpha     = fitParaboloid->GetParameter(4);
        
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
	  fitParaboloid->SetRange(fMass - sigmaLevel*fMassError, fJES - sigmaLevel*fJESError,
	                        minMass + sigmaLevel*fMassError, fJES + sigmaLevel*fJESError);

    sumLogLikelihood->Fit("fitParaboloid","EMR0");
    
    double contours[3] = {1, 4, 9};
    fitParaboloid->SetContour(3, contours);
    
    if (firstBinMass+1 < fitParaboloid->GetParameter(0) && fitParaboloid->GetParameter(0) < lastBinMass-1) {
      fMass = fitParaboloid->GetParameter(0);
      fJES  = fitParaboloid->GetParameter(2);
      if (TMath::Sqrt(1/fitParaboloid->GetParameter(1)) < 2*TMath::Sqrt(fMass)) {
        semiMajor = TMath::Sqrt(1/fitParaboloid->GetParameter(1));
        semiMinor = TMath::Sqrt(1/fitParaboloid->GetParameter(3));
        alpha     = fitParaboloid->GetParameter(4);
        
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
    
    fitParaboloid->SetParameter(5, 0);
    
    // stat+syst ellipsis  
    double mSyst = 1.18;
    double jSyst = 0.012;
    
    double sm2 = fMassError*fMassError + mSyst*mSyst;
    double sj2 = fJESError*fJESError + jSyst*jSyst;
    
    fitParaboloid->Copy(*systParaboloid);
    systParaboloid->SetRange(fMass - 40, fJES - 0.4, minMass + 40, fJES + 0.4);
    systParaboloid->SetParameter(1, 2*cos(2*alpha)/(sm2 - sj2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
    systParaboloid->SetParameter(3, 2*cos(2*alpha)/(sj2 - sm2 + sj2*cos(2*alpha) + sm2*cos(2*alpha)));
    systParaboloid->SetLineColor(kBlack);
    systParaboloid->SetLineStyle(7);
    //systParaboloid->SetLineWidth(5);
    
    //* Set minL to 0
    sumLogLikelihood->Add(hUnity, -sumLogLikelihood->GetMinimum(0) + 1e-2);
    sumLogLikelihood->SetAxisRange(0, 25, "Z");
    if (blackWhite) sumLogLikelihood->Draw("AXIG");
    else {
      sumLogLikelihood->Draw("COLZ");
      //sumLogLikelihood->Draw("CONT3, SAME");
    }
    if (syst) systParaboloid->Draw("cont3 same");
    fitParaboloid->Draw("cont3 same");
    //*/
    
    // create legend
    TLegend *leg0 = new TLegend(0.2, 0.15, 0.45, 0.25);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    //leg0->AddEntry((TObject*)0, "1, 2, 3#sigma", "");
    leg0->AddEntry(fitParaboloid, "stat", "L");
    if (syst)    leg0->AddEntry(systParaboloid, "stat + syst", "L");
    //leg0->Draw();
    
    helper->DrawCMSPrel();
    
    TString path("plot/Ideogram/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".eps";
    ctemp->Print(path);
  }
  
  // 1D top mass
  
  if (!fit2D) {
  
    TH1D* sumLogLikelihood1D = sumLogLikelihood->ProjectionX("sumLogLikelihood1D", sumLogLikelihood->GetYaxis()->FindBin(1.), sumLogLikelihood->GetYaxis()->FindBin(1.));
    sumLogLikelihood1D->Draw("E");
    
    fitParabola->SetParameter(2, sumLogLikelihood1D->GetMinimum(0));
    fitParabola->SetParameter(1, 100);
    fitParabola->SetRange(minMass - 2, minMass + 2);

    sumLogLikelihood1D->Fit("fitParabola","EMR");
    
    if (firstBinMass+1 < fitParabola->GetParameter(0) && fitParabola->GetParameter(0) < lastBinMass-1) {
      fMassAlt      = fitParabola->GetParameter(0);
      fMassAltError = TMath::Sqrt(1/fitParabola->GetParameter(1));
    }
    else {
      fMassAlt      = -1;
      fMassAltError = -1;
    }
    
    helper->DrawCMSPrel();
    
    TString path1D("plot/Ideogram/"); path1D+= fIdentifier; path1D += "_"; path1D += i; path1D += "_"; path1D += j; path1D += "_1D.eps";
    ctemp->Print(path1D);
    
    sumLogLikelihood1D->Delete();
  }
  
  delete fitParaboloid;
  delete null;
  delete eventTree;
  delete ctemp;
  delete eventLikelihood;
  delete logEventLikelihood;
  delete sumLogLikelihood;
  
  std::cout << "IdeogramAnalyzer done" << std::endl;
}


