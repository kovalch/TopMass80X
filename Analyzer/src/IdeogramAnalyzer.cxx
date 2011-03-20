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

  TF1* combLikelihood = new TF1("combLikelihood",
    "[2] * 0.7 * (TMath::Gaus([0],x,[1],1) + 0.3 * 1/(3.5+1) * (3.5 * TMath::Gaus([0], (6.329e+01 + [0]*0.6317), 25, 1) + TMath::Gaus([0], 225, 46, 1)))");
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
  fChain->SetBranchAddress("hadTopMass", &hadTopMass);
  fChain->SetBranchAddress("fitChi2", &fitChi2);
  fChain->SetBranchAddress("fitProb", &fitProb);
  
  int combi;
  fChain->SetBranchAddress("combi", &combi);
  
  TTree* eventTree = fChain->CopyTree(cuts);

  for (int ievent = 0; ievent < eventTree->GetEntries(); ievent++) {
    if (ievent%1000==0) std::cout << "Processing ideogram combi " << ievent << std::endl;
    eventTree->GetEntry(ievent);
    if (combi!=0) continue;
    if (fitProb > 0.05) { // skip bad events
      eventLikelihoodWeight = 1;
      eventLikelihood->Eval(null);
 
      for (int icombi = 0; icombi < 12; icombi++) {
	      eventTree->GetEntry(ievent + icombi);
	      
	      if (icombi!=0 && combi==0 || fitProb < 0.05) break;
	      
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

  sumLogLikelihood->Fit("fitParabola","BWR");
  
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
}
