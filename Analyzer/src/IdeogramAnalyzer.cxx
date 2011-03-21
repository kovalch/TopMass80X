#include "IdeogramAnalyzer.h"

double IdeogramAnalyzer::GetMass() {
  return fMass;
}

void IdeogramAnalyzer::Analyze(TString cuts, int i, int j) {
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int firstbin = 150;
  int lastbin  = 200;

  int bins = 50;
  
  IdeogramCombLikelihood* fptr = new IdeogramCombLikelihood();
  TF1* combLikelihood = new TF1("combLikelihood",fptr,&IdeogramCombLikelihood::Evaluate,150,200,3);

  TF1* fitParabola = new TF1("fitParabola", "abs([1])*(x-[0])^2+[2]");
  fitParabola->SetParLimits(0, firstbin, lastbin);
  fitParabola->SetParameter(0, (lastbin+firstbin)/2);
  fitParabola->SetParLimits(1, 0.01, 10^6);

  TF1* null = new TF1("null", "0");

  TH1F* eventLikelihood = new TH1F("eventLikelihood","eventLikelihood", bins, firstbin, lastbin);
  TH1F* logEventLikelihood = new TH1F("logEventLikelihood", "logEventLikelihood", bins, firstbin, lastbin);
  TH1F* sumLogLikelihood = new TH1F("sumLogLikelihood", "sumLogLikelihood", bins, firstbin, lastbin);
  sumLogLikelihood->Eval(null);

  double hadTopMass, fitChi2, fitProb, weight, eventLikelihoodWeight;
  int combi;
  
  TTree* eventTree = fTree->CopyTree(cuts);
  
  eventTree->SetBranchAddress("hadTopMass", &hadTopMass);
  eventTree->SetBranchAddress("fitChi2", &fitChi2);
  eventTree->SetBranchAddress("fitProb", &fitProb);
  eventTree->SetBranchAddress("combi", &combi);
  
//  fTree->Print();
//  eventTree->Print();
  
  for (int iEntry = 0; iEntry < eventTree->GetEntries(); iEntry++) {
    eventTree->GetEntry(iEntry);
    
    if (combi!=0) continue;
    
    if (fitProb > 0.05) { // skip bad events
      eventLikelihoodWeight = 1;
      eventLikelihood->Eval(null);
 
      for (int iComb = 0; iComb < 12; iComb++) {
//        std::cout << combi << std::endl;
	      eventTree->GetEntry(iEntry + iComb);
	      
	      if (iComb!=0 && combi==0 || fitProb < 0.05) break;
	      
        weight = fitProb;
	      combLikelihood->SetParameters(hadTopMass,12,weight);
	      eventLikelihood->Eval(combLikelihood, "A"); // add combi pdf
      }

      logEventLikelihood->Eval(null);

      for (int i = 0; i<=bins; i++) {
    	  logEventLikelihood->SetBinContent(i, -2*TMath::Log(eventLikelihood->GetBinContent(i)));
      }

      sumLogLikelihood->Add(logEventLikelihood);
    }
  }

  sumLogLikelihood->SetAxisRange(sumLogLikelihood->GetMinimum(0), sumLogLikelihood->GetMaximum(), "Y");
  sumLogLikelihood->Draw();
  
  fitParabola->SetParameter(2, sumLogLikelihood->GetMinimum(0));
  fitParabola->SetParameter(1, 1000);
  
  fitParabola->SetRange(sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) - 2, sumLogLikelihood->GetBinCenter(sumLogLikelihood->GetMinimumBin()) + 2);

  sumLogLikelihood->Fit("fitParabola","QBWR");
  
  if (firstbin+1 < fitParabola->GetParameter(0) && fitParabola->GetParameter(0) < lastbin-1) {
    fMass = fitParabola->GetParameter(0);
    if (TMath::Sqrt(1/fitParabola->GetParameter(1)) < TMath::Sqrt(fMass)) {
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
  
  delete ctemp;
  delete eventLikelihood;
  delete logEventLikelihood;
  delete sumLogLikelihood;
}
