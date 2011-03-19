#include "IdeogramAnalyzer.h"

double IdeogramAnalyzer::GetMass() {
  return fMass;
}

void IdeogramAnalyzer::Analyze(TString cuts, int i, int j) {
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int firstbin = 100;
  int lastbin  = 250;

  int bins = 50;

  TF1* model = new TF1("model",
    "[2] * 0.7 * (TMath::Gaus([0],x,[1],1) + 0.3 * 1/(3.5+1) * (3.5 * TMath::Gaus([0], (6.329e+01 + [0]*0.6317), 25, 1) + TMath::Gaus([0], 225, 46, 1)))");
  // TODO? differential background
  TF1* para = new TF1("para", "abs([1])*(x-[0])^2+[2]");
  para->SetParLimits(0, firstbin, lastbin);
  para->SetParameter(0, (lastbin+firstbin)/2);
  para->SetParLimits(1, 0.01, 10^6);

  TF1* null = new TF1("null", "0");

  TH1F* sum = new TH1F("sum","sum", bins, firstbin, lastbin);
  TH1F* logevent = new TH1F("logevent", "logevent", bins, firstbin, lastbin);
  TH1F* logsum = new TH1F("logsum", "logsum", bins, firstbin, lastbin);
  logsum->Eval(null);

  double hadTopMass, fitChi2, fitProb, weight, sumWeight;
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
      sumWeight = 1;
      sum->Eval(null);
 
      for (int icombi = 0; icombi < 12; icombi++) {
	      eventTree->GetEntry(ievent + icombi);
	      
	      if (icombi!=0 && combi==0) break;
	      
	      if (fitProb > 0.05) {
	        weight = fitProb;
	        model->SetParameters(hadTopMass,12,weight);
	        sum->Eval(model, "A"); // add combi pdf
	      }
      }

      logevent->Eval(null);

      for (int i = 0; i<=bins; i++) {
    	  logevent->SetBinContent(i, -2*TMath::Log(sum->GetBinContent(i)));
      }

      logsum->Add(logevent);
    }

  }

  logsum->SetAxisRange(logsum->GetMinimum(0), logsum->GetMaximum(), "Y");
  logsum->Draw();
  
  para->SetParameter(2, logsum->GetMinimum(0));
  para->SetParameter(1, 1000);
  
  para->SetRange(firstbin+(lastbin-firstbin)/bins*(logsum->GetMinimumBin()-6), firstbin+(lastbin-firstbin)/bins*(logsum->GetMinimumBin()+5));

  logsum->Fit("para","BWR");
  
  fMass      = para->GetParameter(0);
  if (TMath::Sqrt(1/para->GetParameter(1)) < TMath::Sqrt(fMass)) {
    fMassError = TMath::Sqrt(1/para->GetParameter(1));
  }
  fMassSigma = 0;
  
  TString path("plot/Ideogram/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
}
