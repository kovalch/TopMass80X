#include "IdeogramAnalyzer.h"

double IdeogramAnalyzer::GetMass() {
  return fMass;
}

void IdeogramAnalyzer::Analyze(TString cuts, int i, int j) {
/*  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();

  gaus = new TF1("gaus", "gaus");

  fChain->Draw("hadTopMass", cuts);

  fChain->Fit("gaus", "hadTopMass", cuts);
  
  TString path("plot/Ideogram/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
        
  fMass      = gaus->GetParameter(1);
  fMassError = gaus->GetParError(1);
  fMassSigma = gaus->GetParameter(2);
  
  ==================================================*/
  
  TCanvas* ctemp = new TCanvas("ctemp", "Top mass", 500, 500);
  ctemp->cd();
  
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  int firstbin = 100;
  int lastbin  = 250;

  int bins = 50;

  TF1* model = new TF1("model",
    "[2] * (TMath::Gaus([0],x,[1],1) + TMath::GammaDist([0], 2, -2451.53 + 29.2008*[0] - 0.0825*[0]^2, 1153.57 - 12.8479*[0] + 0.0370833*[0]^2))");
  // TMath::Landau([0],190,28/190,1)
  // TMath::GammaDist([0], 2, 1.3775 + 0.738333*[0], 51.0038 - 0.0541667*[0])
  // TMath::GammaDist([0], 2, -2451.53 + 29.2008*[0] - 0.0825*[0]^2, 1153.57 - 12.8479*[0] + 0.0370833*[0]^2)
  // TODO signal fraction
  TF1* para = new TF1("para", "abs([1])*(x-[0])^2+[2]");
  para->SetParLimits(0, firstbin, lastbin);
  para->SetParameter(0, (lastbin+firstbin)/2);
  para->SetParLimits(1, 0.05, 10^6);

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
	        model->SetParameters(hadTopMass,18,weight);
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
  
  para->SetRange(firstbin+(lastbin-firstbin)/bins*(logsum->GetMinimumBin()-5), firstbin+(lastbin-firstbin)/bins*(logsum->GetMinimumBin()+5));

  logsum->Fit("para","BWR");
  
  fMass      = para->GetParameter(0);
  if (TMath::Sqrt(1/para->GetParameter(1)) < TMath::Sqrt(fMass)) {
    fMassError = TMath::Sqrt(1/para->GetParameter(1));
  }
  fMassSigma = 0;
  
  TString path("plot/Ideogram/"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += ".png";
  ctemp->Print(path);
}
