#include "basicFunctions.h"

double getSoB(TH1F* sig, TH1F* bkg, double probCut, TString opt);
void drawArrow(const double xmin, const double y, const double xmax, const unsigned int color, const double lineWidth, const unsigned int lineStyle, double stretchFactor);

void optimizeProbCut(double forceProb = 0.02, unsigned int ndof = 2, int verbose=1, TString LepDecayChannel = "combined", TString inputFolderName=AnalysisFolder, bool save=true)
{
  // configure style
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  // set optimization parameter
  // optimize the product signal efficiency times signal purity:
  TString optimize = "#frac{sig}{#sqrt{sig+bkg}}";
  // optimize the product signal efficiency times signal to background ratio:
  // TString optimize = "#frac{sig}{#sqrt{bkg}}"; 
  // optimize signal purity (don't use this one!):
  // TString optimize = "#frac{sig}{sig+bkg}";

  // read files
  TString path = groupSpace+inputFolderName+"/";
  std::vector<TFile*> files_;
  if(LepDecayChannel=="combined" || LepDecayChannel=="muon"    )  files_.push_back(new TFile(path+TopFilename(kSig, sysNo,"muon"    )));
  if(LepDecayChannel=="combined" || LepDecayChannel=="muon"    )  files_.push_back(new TFile(path+TopFilename(kBkg, sysNo,"muon"    )));
  if(LepDecayChannel=="combined" || LepDecayChannel=="electron")  files_.push_back(new TFile(path+TopFilename(kSig, sysNo,"electron")));
  if(LepDecayChannel=="combined" || LepDecayChannel=="electron")  files_.push_back(new TFile(path+TopFilename(kBkg, sysNo,"electron")));

  // get trees
  vector<TTree*> trees_;
  float decayChannel;
  float chi2;
  float qAssignment;
  float weight;
  for(UInt_t i = 0; i < files_.size(); i++){
    trees_.push_back((TTree*)files_[i]->Get("analyzeTopRecoKinematicsKinFit/tree"));
    trees_[i]->SetBranchAddress("decayChannel",&decayChannel);
    trees_[i]->SetBranchAddress("chi2",&chi2);
    trees_[i]->SetBranchAddress("qAssignment",&qAssignment);
    trees_[i]->SetBranchAddress("weight",&weight);
  }

  // initialize histograms
  // only right permutaions
  TH1F* sigHisto = new TH1F("probHistSig","signal probability",1000000,0.,1.);
  // wrong permutations and ttbar bkg
  TH1F* bkgHisto = new TH1F("probHistBkg","background probability",1000000,0.,1.);
  // wrong permutations
  TH1F* bkgPermuHisto = new TH1F("probHistBkgPermu","non correct permutation probability",1000000,0.,1.);
  // all ttbar signal (right+wrong permutations)
  TH1F* ttSigHisto = new TH1F("histTtSig","ttSig probability",1000000,0.,1.);
  // all ttbar background
  TH1F* ttBkgHisto = new TH1F("histTtBkg","ttBkg probability",1000000,0.,1.);
  // kinFit permutations
  TH1F* permutation     = new TH1F("permutation","permutation", 10, 0., 10. );
  permutation->GetXaxis()->SetBinLabel(1, "#splitline{all jets correct}{assigned}");
  permutation->GetXaxis()->SetBinLabel(2, "b^{lep}#leftrightarrowb^{had}");
  permutation->GetXaxis()->SetBinLabel(3, "b^{lep}#leftrightarrowq"   );
  permutation->GetXaxis()->SetBinLabel(4, "b^{had}#leftrightarrowq"   );
  //permutation->GetXaxis()->SetBinLabel(5, "#splitline{b^{lep}#rightarrowb^{had},q#rightarrowb^{lep}}{& b^{had}#rightarrowq}");
  //permutation->GetXaxis()->SetBinLabel(6, "#splitline{b^{had}#rightarrowb^{lep},q#rightarrowb^{had}}{& b^{lep}#rightarrowq}"  );

  permutation->GetXaxis()->SetBinLabel(5, "#splitline{    gen: b^{lep}qb^{had}}{#rightarrow rec: b^{had}b^{lep}q}");
  permutation->GetXaxis()->SetBinLabel(6, "#splitline{    gen: b^{had}qb^{lep}}{#rightarrow rec: b^{lep}b^{had}q}");
  permutation->GetXaxis()->SetBinLabel(7, "#splitline{all jets}{swapped}"       );
  permutation->GetXaxis()->SetBinLabel(8, "#splitline{#geq1 jet not}{in lead 5}" );
  permutation->GetXaxis()->SetBinLabel(9, "#splitline{#geq1 wrong}{jet chosen}");
  permutation->GetXaxis()->SetBinLabel(10,"#splitline{no jet-parton}{matching}"   );
  // ttbar decaymode
  TH1F* decaychannel = new TH1F("decaychannel", "decaychannel",  25,   -4.5,   20.5);
  decaychannel->GetXaxis()->SetBinLabel(1 , "??"      );
  decaychannel->GetXaxis()->SetBinLabel(2 , "ntt"     );
  decaychannel->GetXaxis()->SetBinLabel(3 , "ll?"     );
  decaychannel->GetXaxis()->SetBinLabel(4 , "l?j"     );
  decaychannel->GetXaxis()->SetBinLabel(6 , "#splitline{e+j}{BG}"  );
  decaychannel->GetXaxis()->SetBinLabel(7 , "#splitline{#mu+j}{BG}");
  decaychannel->GetXaxis()->SetBinLabel(8 , "#tauj"   );
  decaychannel->GetXaxis()->SetBinLabel(16, "ee"      );
  decaychannel->GetXaxis()->SetBinLabel(17, "e#mu"    );
  decaychannel->GetXaxis()->SetBinLabel(18, "e#tau"   );
  decaychannel->GetXaxis()->SetBinLabel(19, "#mu#mu"  );
  decaychannel->GetXaxis()->SetBinLabel(20, "#mu#tau" );
  decaychannel->GetXaxis()->SetBinLabel(21, "#tau#tau");
  decaychannel->GetXaxis()->SetBinLabel(25, "jj"      );
  TH1F* decaychannelProb =(TH1F*)decaychannel->Clone("decaychannelProb");
  
  // kinFit permutations for a optimal probability selection
  TH1F* permutationProb = (TH1F*)(permutation->Clone("permutationProb"));
  
  // fill histograms
  Long64_t nevent = 0;
  for(UInt_t i = 0; i < files_.size(); i++)
    {
      nevent = (Long64_t)trees_[i]->GetEntries();
      for(Long64_t ientry = 0; ientry < nevent; ++ientry){
	trees_[i]->GetEntry(ientry);
	double prob = TMath::Prob(chi2,ndof);
	// ttbar SG
	//TString fileName=files_[i]->GetName();
	//if( (fileName.Contains("muon")&&fileName.Contains("Sig")&&(LepDecayChannel=="muon"||LepDecayChannel=="combined")) || (fileName.Contains("elec")&&fileName.Contains("Sig")&&(LepDecayChannel=="electron"||LepDecayChannel=="combined")) ){
	if((i>1 && decayChannel==1) || (i<=1 && decayChannel==2 && !(LepDecayChannel=="electron")) || (LepDecayChannel=="electron" && decayChannel==1)){
	  // ttbar SG all permutations
	  ttSigHisto->Fill(prob,weight);
	  // ttbar SG correct permutation
	  if(qAssignment==0) sigHisto->Fill(prob,weight);
	  else{
	    // ttbar SG wrong permutations
	    bkgPermuHisto->Fill(prob,weight);
	    // ttbar BG + ttbar SG wrong permutations
	    bkgHisto->Fill(prob,weight);
	  }
	}
	// ttbar BG
	else{
	  // ttbar BG
	  ttBkgHisto->Fill(prob,weight);
	  // ttbar BG + ttbar SG wrong permutations
	  bkgHisto->Fill(prob,weight);
	}
	if(verbose>1){
	  std::cout << files_[i]->GetName() << std::endl;
	  std::cout << "lepton: " << LepDecayChannel << std::endl;
	  std::cout << "decayChannel: " << decayChannel << std::endl;
	  std::cout << "qAssignment: "  << qAssignment  << std::endl;
	}	  
      }
    }

  // find optimal probability cut
  double goldSec     = (3.-sqrt(5.))/2.;
  double optimalProb = -1.;
  double optimalSoB  = -1.;
  double optimalEff  = -1.;
  double SoBLowProb  = 0.;
  double SoBHighProb = 0.;
  double lowProb     = 0.;
  double highProb    = 1.;
  double newLowProb  = goldSec;
  double newHighProb = 1 - goldSec;
  double sigEvents   = sigHisto->Integral();
  double ttSigEvents = ttSigHisto->Integral();
  double ttBkgEvents = ttBkgHisto->Integral();
  vector<double> probVec;
  vector<double> SoBVec;
  vector<double> SeffVec;
  vector<double> BGeffVec;
  vector<double> SGeffVec;

  if(verbose>1){
    std::cout << std::endl << "optimizing prob cut for " << optimize << std::endl;  
    std::cout << " (sig=ttbar SG correct permutations)" << std::endl;  
    std::cout << " (bkg=ttbar BG + ttbar SG wrong permutations" << std::endl << std::endl;
  }

  // calculate some points for the plots
  for ( int i = 0 ; i < 100 ; i++ ) {
    probVec.push_back(i/100.);
    SoBVec.push_back(getSoB(sigHisto,bkgHisto,i/100.,optimize));
    SeffVec.push_back(sigHisto->Integral(sigHisto->FindBin(i/100.),sigHisto->GetNbinsX()+1)/sigEvents);
    SGeffVec.push_back(ttSigHisto->Integral(ttSigHisto->FindBin(i/100.),ttSigHisto->GetNbinsX()+1)/ttSigEvents);
    BGeffVec.push_back(ttBkgHisto->Integral(ttBkgHisto->FindBin(i/100.),ttBkgHisto->GetNbinsX()+1)/ttBkgEvents);
  }

  if(forceProb>=0 && forceProb<=1){
    optimalProb=forceProb;
    optimalSoB=getSoB(sigHisto,bkgHisto,forceProb,optimize);
    optimalEff=sigHisto->Integral(sigHisto->FindBin(forceProb),sigHisto->GetNbinsX()+1)/sigEvents;
    // calculate some points in the region of the forced cut value
    for ( int i = 0 ; i < 400 ; i++ ) {
      probVec.push_back(optimalProb-0.01+i/20000.);
      SoBVec.push_back(getSoB(sigHisto,bkgHisto,optimalProb-0.01+i/20000.,optimize));
      SeffVec.push_back(sigHisto->Integral(sigHisto->FindBin(optimalProb-0.01+i/20000.),sigHisto->GetNbinsX()+1)/sigEvents);
      SGeffVec.push_back(ttSigHisto->Integral(ttSigHisto->FindBin(optimalProb-0.01+i/20000.),ttSigHisto->GetNbinsX()+1)/ttSigEvents);
      BGeffVec.push_back(ttBkgHisto->Integral(ttBkgHisto->FindBin(optimalProb-0.01+i/20000.),ttBkgHisto->GetNbinsX()+1)/ttBkgEvents);
    }
  }
  else{
    // golden section search for optimal cut value
    for ( int i = 0 ; i < 100 ; i++ ) {
      if( optimalProb!=newLowProb ){
	// calculate s/sqrt(s+b) or s/sqrt(b) for lowProb
	SoBLowProb = getSoB(sigHisto,bkgHisto,newLowProb,optimize);
	probVec.push_back(newLowProb);
	SoBVec.push_back(SoBLowProb);
	SeffVec .push_back(sigHisto  ->Integral(sigHisto  ->FindBin(newLowProb),sigHisto  ->GetNbinsX()+1)/sigEvents);
	SGeffVec.push_back(ttSigHisto->Integral(ttSigHisto->FindBin(newLowProb),ttSigHisto->GetNbinsX()+1)/ttSigEvents);
	BGeffVec.push_back(ttBkgHisto->Integral(ttBkgHisto->FindBin(newLowProb),ttBkgHisto->GetNbinsX()+1)/ttBkgEvents);
	cout << "newLowProb: " << newLowProb << " -> "+optimize+": " << SoBLowProb << endl;
      }
      if(optimalProb!=newHighProb){
	// calculate s/sqrt(s+b) or s/sqrt(b) for highProb
	SoBHighProb = getSoB(sigHisto,bkgHisto,newHighProb,optimize);
	probVec.push_back(newHighProb);
	SoBVec.push_back(SoBHighProb);
	SeffVec .push_back(sigHisto  ->Integral(sigHisto  ->FindBin(newHighProb),sigHisto  ->GetNbinsX()+1)/sigEvents);
	SGeffVec.push_back(ttSigHisto->Integral(ttSigHisto->FindBin(newHighProb),ttSigHisto->GetNbinsX()+1)/ttSigEvents);
	BGeffVec.push_back(ttBkgHisto->Integral(ttBkgHisto->FindBin(newHighProb),ttBkgHisto->GetNbinsX()+1)/ttBkgEvents);
	cout << "newHighProb: " << newHighProb << " -> "+optimize+": " << SoBHighProb << endl;
      }
      if(SoBLowProb>SoBHighProb){
	highProb=newHighProb;
	optimalProb=newLowProb;
	optimalSoB=SoBLowProb;
	newHighProb=newLowProb;
	optimalEff=sigHisto->Integral(sigHisto->FindBin(newLowProb),sigHisto->GetNbinsX()+1)/sigEvents;
	if(TMath::Abs(SoBLowProb-SoBHighProb)/(SoBLowProb+SoBHighProb)<0.0000001)break;
	SoBHighProb=SoBLowProb;
	newLowProb=lowProb+goldSec*(highProb-lowProb);
      } else {
	lowProb=newLowProb;
	optimalProb=newHighProb;
	optimalSoB=SoBHighProb;
	newLowProb=newHighProb;
	optimalEff=sigHisto->Integral(sigHisto->FindBin(newHighProb),sigHisto->GetNbinsX()+1)/sigEvents;
	if(TMath::Abs(SoBLowProb-SoBHighProb)/(SoBLowProb+SoBHighProb)<0.0000001)break;
	SoBLowProb=SoBHighProb;
	newHighProb=highProb-goldSec*(highProb-lowProb);
      }
    }
  }
  double ttSigEff = ttSigHisto->Integral(ttSigHisto->FindBin(optimalProb),ttSigHisto->GetNbinsX()+1)/ttSigEvents;
  double ttBkgEff = ttBkgHisto->Integral(ttBkgHisto->FindBin(optimalProb),ttBkgHisto->GetNbinsX()+1)/ttBkgEvents;
  double ttEff    = (ttSigHisto->Integral(ttSigHisto->FindBin(optimalProb),ttSigHisto->GetNbinsX()+1)+ttBkgHisto->Integral(ttBkgHisto->FindBin(optimalProb),ttBkgHisto->GetNbinsX()+1))/(ttSigEvents+ttBkgEvents);
  double wrongEff=bkgPermuHisto->Integral(bkgPermuHisto->FindBin(optimalProb),bkgPermuHisto->GetNbinsX()+1)/bkgPermuHisto->Integral(0,bkgPermuHisto->GetNbinsX()+1);
  if(verbose>0){
    cout << endl << "Optimal probability cut: prob > " << optimalProb << endl;
    cout << " -> "<< endl;
    cout << " eff(all ttbar          ): " << ttEff << endl;
    cout << " eff(correct permutation): " << optimalEff << endl;
    cout << " eff(wrong   permutation): " << wrongEff << endl;
    cout << " eff(ttbar-signal       ): " << ttSigEff << endl;
    cout << " eff(ttbar-bkg          ): " << ttBkgEff << endl;
    cout << " "<< optimize << "        : " << getSoB(sigHisto,bkgHisto,-1,optimize);
    cout << "  (" << getSoB(sigHisto,bkgHisto,optimalProb,optimize) << "   for prob>" << optimalProb << ")" << endl;
    cout << " ttbar SG / ttbar BG               : " <<  ttSigEvents/ttBkgEvents;
    cout << "  (" <<  ttSigHisto->Integral(ttSigHisto->FindBin(optimalProb),ttSigHisto->GetNbinsX()+1)/ttBkgHisto->Integral(ttBkgHisto->FindBin(optimalProb),ttBkgHisto->GetNbinsX()+1)  << "  for prob>" << optimalProb << ")" << endl;
    cout << " correct permutation / all ttbar SG: " <<  sigHisto->Integral(0,sigHisto->GetNbinsX()+1)/ttSigEvents;
    cout << " (" <<  sigHisto->Integral(sigHisto->FindBin(optimalProb),sigHisto->GetNbinsX()+1)/ttSigHisto->Integral(ttSigHisto->FindBin(optimalProb),ttSigHisto->GetNbinsX()+1) << " for prob>" << optimalProb << ")" << endl;
  }

  // fill permutation histos
  for(UInt_t i = 0; i < files_.size(); i++)
    {
      nevent = (Long64_t)trees_[i]->GetEntries();
      for(Long64_t ientry = 0; ientry < nevent; ++ientry){
	trees_[i]->GetEntry(ientry);
	double prob = TMath::Prob(chi2,ndof);
	if((i>1 && decayChannel==1) || (i<=1 && decayChannel==2 && !(LepDecayChannel=="electron")) || (LepDecayChannel=="electron" && decayChannel==1)){
	  permutation->Fill(qAssignment,weight);
	  if(prob>optimalProb) permutationProb->Fill(qAssignment,weight);
	}
	else{
	  decaychannel->Fill(decayChannel,weight);
	  if(prob>optimalProb) decaychannelProb->Fill(decayChannel,weight);
	}
      }
    }
  // get efficiency per permutation
  TH1F* effPermu = (TH1F*)(permutationProb->Clone("effPermu"));
  effPermu->Divide((TH1F*)(permutation->Clone()));
  // get efficiency per decay channel
  TH1F* effDecay = (TH1F*)(decaychannelProb->Clone("effDecay"));
  effDecay->Divide((TH1F*)(decaychannel->Clone()));
  // normalize permutation histos
  permutation    ->Scale(1./permutation    ->Integral(0,permutation    ->GetNbinsX()+1));
  permutationProb->Scale(1./permutationProb->Integral(0,permutationProb->GetNbinsX()+1));

  // normalize decay channel histos
  decaychannel    ->Scale(1./decaychannel    ->Integral(0,decaychannel    ->GetNbinsX()+1));
  decaychannelProb->Scale(1./decaychannelProb->Integral(0,decaychannelProb->GetNbinsX()+1));

  // set up canvas
  TCanvas *canv = new TCanvas("canv","probabilityCutOptimisation",10,10,1200,600);
  canv->Divide(2,1);
  canv->cd(1)->SetLeftMargin  (0.15);
  canv->cd(1)->SetRightMargin (0.03);
  canv->cd(1)->SetBottomMargin(0.10);
  canv->cd(1)->SetTopMargin   (0.05);
  canv->cd(2)->SetLeftMargin  (0.15);
  canv->cd(2)->SetRightMargin (0.03);
  canv->cd(2)->SetBottomMargin(0.10);
  canv->cd(2)->SetTopMargin   (0.05);
  TCanvas *canv2 = new TCanvas("canv2","probability",10,10,900,600);
  canv2->cd()->SetLeftMargin  (0.08);
  canv2->cd()->SetRightMargin (0.03);
  canv2->cd()->SetBottomMargin(0.10);
  canv2->cd()->SetTopMargin   (0.05);
  canv2->cd()->SetLogy();
  TCanvas *canv3 = new TCanvas("canv3","probabilityFraction",10,10,900,600);
  canv3->cd()->SetLeftMargin  (0.08);
  canv3->cd()->SetRightMargin (0.03);
  canv3->cd()->SetBottomMargin(0.10);
  canv3->cd()->SetTopMargin   (0.05);
  canv3->cd()->SetLogy();
  TCanvas *canv4 = new TCanvas("canv4","permutation",10,10,900,600);
  canv4->cd()->SetLeftMargin  (0.20);
  canv4->cd()->SetRightMargin (0.03);
  canv4->cd()->SetBottomMargin(0.10);
  canv4->cd()->SetTopMargin   (0.05);
  TCanvas *canv5 = new TCanvas("canv5","effPerPermutation",10,10,900,600);
  canv5->cd()->SetLeftMargin  (0.20);
  canv5->cd()->SetRightMargin (0.03);
  canv5->cd()->SetBottomMargin(0.10);
  canv5->cd()->SetTopMargin   (0.05);
  TCanvas *canv6 = new TCanvas("canv6","ttbarBGcomposition",10,10,900,600);
  canv6->cd()->SetLeftMargin  (0.08);
  canv6->cd()->SetRightMargin (0.03);
  canv6->cd()->SetBottomMargin(0.10);
  canv6->cd()->SetTopMargin   (0.05);
  TCanvas *canv7 = new TCanvas("canv7","ttbarBGeff",10,10,900,600);
  canv7->cd()->SetLeftMargin  (0.08);
  canv7->cd()->SetRightMargin (0.03);
  canv7->cd()->SetBottomMargin(0.10);
  canv7->cd()->SetTopMargin   (0.05);

  // draw graph s/sqrt(s+b) or s/sqrt(b) vs probability
  canv->cd(1);
  TGraph* optSoB = new TGraph((int)probVec.size(),&probVec.front(),&SoBVec.front());
  optSoB->SetMarkerStyle(20);
  sigHisto->SetMaximum(1.1*optimalSoB);
  sigHisto->GetXaxis()->SetTitle("probability");
  sigHisto->GetYaxis()->SetTitle(optimize);
  sigHisto->GetYaxis()->SetTitleOffset(1.65);
  sigHisto->DrawClone("axis");
  optSoB->DrawClone("p same");
  TPad *rPad = new TPad("rPad","",0.2,0.12,0.7,0.62);
  rPad->SetFillStyle(0);
  rPad->SetFillColor(0);
  rPad->SetBorderSize(0);
  rPad->SetBorderMode(0);
  rPad->SetLogy(0);
  rPad->SetLogx(0);
  rPad->SetTicky(1);
  rPad->Draw("");
  rPad->cd();
  sigHisto->GetXaxis()->SetRangeUser(optimalProb-0.01,optimalProb+0.01);
  sigHisto->GetYaxis()->SetRangeUser(optimalSoB-0.4, optimalSoB+0.1);
  sigHisto->GetYaxis()->SetTitle("");
  sigHisto->DrawClone("axis");
  drawLine(optimalProb, optimalSoB-0.4, optimalProb, optimalSoB+0.05, kRed, 3, 1);
  optSoB->DrawClone("p same");  
  // draw graph signal efficiency vs probability
  canv->cd(2);
  TGraph* effCorrect = new TGraph((int)probVec.size(),&probVec.front(),&SeffVec.front() );
  TGraph* effSG      = new TGraph((int)probVec.size(),&probVec.front(),&SGeffVec.front());
  TGraph* effBG      = new TGraph((int)probVec.size(),&probVec.front(),&BGeffVec.front());
  effCorrect->SetMarkerStyle(20);
  effSG->SetMarkerStyle(21);
  effBG->SetMarkerStyle(22);
  effSG->SetMarkerColor(kRed);
  effBG->SetMarkerColor(kBlue);
  sigHisto->SetMaximum(0.6);
  sigHisto->GetYaxis()->SetTitle("efficiency");
  sigHisto->GetXaxis()->SetRangeUser(0.,1.);
  sigHisto->DrawClone("axis");
  effCorrect->DrawClone("p same");
  effSG->DrawClone("p same");
  effBG->DrawClone("p same");
  TLegend *legEff  = new TLegend(); 
  legEff->SetX1NDC(0.3);
  legEff->SetY1NDC(0.65);
  legEff->SetX2NDC(0.9);
  legEff->SetY2NDC(0.9);
  legEff ->SetFillStyle(1001);
  legEff ->SetFillColor(10);
  legEff ->SetBorderSize(0);
  legEff ->SetTextAlign(12);
  legEff ->SetHeader(" t#bar{t} Simulation");
  legEff ->AddEntry(effCorrect, "SG correct permutation"  , "P");
  legEff ->AddEntry(effSG, "SG all permutations"  , "P");
  legEff ->AddEntry(effBG, "BG"  , "P");
  legEff->DrawClone("same");
  // draw prob distribution
  canv2->cd();
  sigHisto->SetTitle("");
  sigHisto->GetXaxis()->SetTitle("Probability (best hypothesis)");
  sigHisto->GetYaxis()->SetTitle("Events");
  sigHisto->GetYaxis()->SetTitleOffset(0.85);
  TH1F* sigFracHisto=(TH1F*)sigHisto->Clone("sigFracHisto");
  sigFracHisto->SetLineColor(kRed+1);
  sigFracHisto->SetLineWidth(3);
  sigFracHisto->SetMinimum(0.001);
  sigFracHisto->SetMaximum(1.00);
  sigFracHisto->GetYaxis()->SetTitle("rel. Events");
  sigHisto->SetFillColor(kRed+1);
  sigHisto->Add(bkgHisto);
  sigHisto->Rebin(10000);
  sigHisto->SetMinimum(10);
  sigFracHisto->Rebin(10000);
  sigHisto->DrawClone("h");
  TH1F* bkgFracHisto=(TH1F*)bkgPermuHisto->Clone("bkgFracHisto");
  bkgFracHisto->SetLineWidth(3);
  bkgFracHisto->SetLineStyle(9);
  bkgFracHisto->SetLineColor(kRed);
  bkgHisto->SetFillStyle(1001);
  bkgHisto->SetFillColor(10);
  bkgFracHisto->Rebin(10000);
  bkgHisto->Rebin(10000);
  bkgHisto->DrawClone("h same");
  bkgHisto->SetFillStyle(3001);
  bkgHisto->SetFillColor(kRed);
  bkgHisto->DrawClone("h same");
  TH1F* ttBkgFracHisto=(TH1F*)ttBkgHisto->Clone("ttBkgFracHisto");
  int ttBkgColor=kBlue;//kRed-7
  ttBkgFracHisto->SetLineColor(ttBkgColor);
  ttBkgFracHisto->SetLineWidth(3);
  ttBkgFracHisto->SetLineStyle(2);
  ttBkgHisto->SetFillColor(10);
  ttBkgHisto->SetFillStyle(1001);
  ttBkgHisto->Rebin(10000);
  ttBkgFracHisto->Rebin(10000);
  ttBkgHisto->DrawClone("h same");
  ttBkgHisto->SetFillColor(ttBkgColor);
  ttBkgHisto->SetFillStyle(3144);
  ttBkgHisto->DrawClone("h same");
  // draw legend
  TLegend *leg  = new TLegend(); 
  leg->SetX1NDC(0.3);
  leg->SetY1NDC(0.6);
  leg->SetX2NDC(0.9);
  leg->SetY2NDC(0.9);
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(10);
  leg ->SetBorderSize(0);
  //leg ->SetTextSize(0.1);
  leg ->SetTextAlign(12);
  leg ->SetHeader(" t#bar{t} Simulation");
  TLegend *leg2=(TLegend *)leg->Clone();
  leg ->AddEntry(sigHisto, "SG correct permutation"  , "F");
  leg ->AddEntry(bkgHisto, "SG other permutation"    , "F");
  leg ->AddEntry(ttBkgHisto, "BG (non prompt l+jets signal)", "F");
  leg ->AddEntry(effPermu, "optimal cut: prob>"+getTStringFromDouble(optimalProb, 2, false), "L");
  leg2 ->AddEntry(sigFracHisto, "SG correct permutation"  , "L");
  leg2 ->AddEntry(bkgFracHisto, "SG other permutation"    , "L");
  leg2 ->AddEntry(ttBkgFracHisto, "BG (non prompt l+jets signal)", "L");
  leg2 ->AddEntry(effPermu, "optimal cut: prob>"+getTStringFromDouble(optimalProb, 2, false), "L");
  leg ->Draw("same");
  drawLine(optimalProb, sigHisto->GetMinimum(), optimalProb, 0.6*sigHisto->GetMaximum(), kBlack, 2, 1);
  drawArrow(optimalProb, sigHisto->GetMaximum()/10, optimalProb+0.05, kBlack, 2, 1, 1.2);
  sigHisto->DrawClone("axis same");

  // draw normalized prob distribution
  canv3->cd();
  sigFracHisto->DrawClone("axis");
  sigFracHisto->DrawNormalized("same");
  bkgFracHisto->DrawNormalized("same");
  ttBkgFracHisto->DrawNormalized("same");
  drawLine(optimalProb, sigFracHisto->GetMinimum(), optimalProb, 0.6*sigFracHisto->GetMaximum(), kBlack, 2, 1);
  drawArrow(optimalProb, sigFracHisto->GetMaximum()/10, optimalProb+0.05, kBlack, 2, 1, 1.2);
  sigFracHisto->DrawClone("axis same");
  leg2->Draw("same");
  // draw permutation distribution
  canv4->cd();
  double maxPerm     = permutation->GetMaximum();
  double maxPermProb = permutationProb->GetMaximum();
  permutation->SetTitle("");
  permutation->GetXaxis()->SetTitle("jet permutation wrt. machted parton truth");
  permutation->GetYaxis()->SetTitle("rel. Events");
  permutation->GetXaxis()->SetTitleOffset(2.5);
  permutation->GetXaxis()->SetLabelOffset(0.008);
  permutation->SetLineWidth(3);
  permutation->SetLineColor(kBlue);
  permutation->SetFillColor(kBlue);
  permutation->SetLineStyle(1);
  permutation->SetFillStyle(1001);
  permutation->SetBarWidth(0.4);
  permutation->SetBarOffset(0.55);
  permutation->SetMaximum(1.05*(maxPerm<maxPermProb ? maxPermProb : maxPerm));
  permutation->DrawClone("hbar");
  permutationProb->SetLineWidth(3);
  permutationProb->SetLineStyle(1);
  permutationProb->SetLineColor(kRed);
  permutationProb->SetFillColor(kRed);
  permutationProb->SetFillStyle(1001);
  permutationProb->SetBarWidth(0.45);
  permutationProb->SetBarOffset(0.1);
  permutationProb->Draw("hbar same");
  permutation    ->DrawClone("hbar same");
  TLegend *leg3  = new TLegend(); 
  leg3->SetX1NDC(0.5);
  leg3->SetY1NDC(0.4);
  leg3->SetX2NDC(0.9);
  leg3->SetY2NDC(0.7);
  leg3 ->SetFillStyle(1001);
  leg3 ->SetFillColor(10);
  leg3 ->SetBorderSize(0);
  leg3 ->SetTextAlign(12);
  leg3->AddEntry(permutation, "t#bar{t}#rightarrowl+jets prompt", "F");
  leg3->AddEntry(permutationProb, "  + prob>"+getTStringFromDouble(optimalProb, 2, false), "F");
  canv4->SetGrid(1,1);
  leg3->Draw("same");
  // draw efficiency per permutation
  canv5->cd();
  effPermu->SetTitle("");
  effPermu->GetXaxis()->SetTitle("jet permutation wrt. machted parton truth");
  effPermu->GetYaxis()->SetTitle("eff(prob>"+getTStringFromDouble(optimalProb, 2, false)+")");
  effPermu->GetXaxis()->SetTitleOffset(2.5);
  effPermu->GetXaxis()->SetLabelOffset(0.008);
  effPermu->SetBarWidth(0.5);
  effPermu->SetBarOffset(0.3);
  effPermu->SetLineWidth(3);
  effPermu->SetLineColor(kBlack);
  effPermu->SetFillColor(kBlue);
  effPermu->SetFillStyle(1001);
  effPermu->SetLineStyle(1);
  effPermu->Draw("hbar");
  canv5->SetGrid(1,1);
  // draw ttbarBG composition
  canv6->cd();
  decaychannel->SetTitle("");
  decaychannel->GetXaxis()->SetTitle("t#bar{t} BG decay channel");
  decaychannel->GetYaxis()->SetTitle("rel. Events");
  decaychannel->SetLineWidth(3);
  decaychannel->SetLineColor(kBlue);
  decaychannel->SetLineStyle(1);
  decaychannel->GetXaxis()->SetRange(5,25);
  decaychannel->SetMaximum(1.2*decaychannelProb->GetMaximum());
  decaychannel    ->Draw();
  decaychannelProb->SetLineWidth(3);
  decaychannelProb->SetLineStyle(2);
  decaychannelProb->SetLineColor(kRed);
  decaychannelProb->Draw("same");
  TLegend *leg4  = new TLegend(); 
  leg4->SetX1NDC(0.35);
  leg4->SetY1NDC(0.6);
  leg4->SetX2NDC(0.7);
  leg4->SetY2NDC(0.9);
  leg4 ->SetFillStyle(1001);
  leg4 ->SetFillColor(10);
  leg4 ->SetBorderSize(0);
  leg4 ->SetTextAlign(12);
  leg4->AddEntry(decaychannel, "t#bar{t} BG (non prompt l+jets signal)", "L");
  leg4->AddEntry(decaychannelProb, "  + prob>"+getTStringFromDouble(optimalProb, 2, false), "L");
  canv6->SetGrid(1,1);
  leg4->Draw("same");
  // draw efficiency per ttbarBG composition
  canv7->cd();
  effDecay->SetTitle("");
  effDecay->GetXaxis()->SetTitle("t#bar{t} BG decay channel");
  effDecay->GetYaxis()->SetTitle("eff(prob>"+getTStringFromDouble(optimalProb, 2, false)+")");
  effDecay->GetXaxis()->SetRange(5,25);
  effDecay    ->SetLineWidth(3);
  effDecay    ->SetLineColor(kBlack);
  effDecay    ->SetLineStyle(1);
  effDecay    ->Draw();
  canv7->SetGrid(1,1);

  if(save){
    if(verbose>0) std::cout << "will save all plots as png/eps/pdf" << std::endl;
    if(verbose<=1) gErrorIgnoreLevel=kWarning;
    TString outputFolder=TString("./diffXSecFromSignal/plots/")+LepDecayChannel+"/2012/kinFitPerformance/";
    // as eps
    canv ->Print(outputFolder+(TString)(canv->GetTitle())+".eps");
    canv2->Print(outputFolder+(TString)(canv2->GetTitle())+".eps");
    canv3->Print(outputFolder+(TString)(canv3->GetTitle())+".eps");
    canv4->Print(outputFolder+(TString)(canv4->GetTitle())+".eps");
    canv5->Print(outputFolder+(TString)(canv5->GetTitle())+".eps");
    canv6->Print(outputFolder+(TString)(canv6->GetTitle())+".eps");
    canv7->Print(outputFolder+(TString)(canv7->GetTitle())+".eps");
    // as png
    canv ->Print(outputFolder+(TString)(canv->GetTitle())+".png");
    canv2->Print(outputFolder+(TString)(canv2->GetTitle())+".png");
    canv3->Print(outputFolder+(TString)(canv3->GetTitle())+".png");
    canv4->Print(outputFolder+(TString)(canv4->GetTitle())+".png");
    canv5->Print(outputFolder+(TString)(canv5->GetTitle())+".png");
    canv6->Print(outputFolder+(TString)(canv6->GetTitle())+".png");
    canv7->Print(outputFolder+(TString)(canv7->GetTitle())+".png");
    // as pdf
    canv ->Print(outputFolder+"optimalProbCut.pdf(");
    canv2->Print(outputFolder+"optimalProbCut.pdf" );
    canv3->Print(outputFolder+"optimalProbCut.pdf" );
    canv4->Print(outputFolder+"optimalProbCut.pdf" );
    canv5->Print(outputFolder+"optimalProbCut.pdf" );
    canv6->Print(outputFolder+"optimalProbCut.pdf" );
    canv7->Print(outputFolder+"optimalProbCut.pdf)");
  }
}

// calculate s/sqrt(s+b) or s/sqrt(+b)
double getSoB(TH1F* sig, TH1F* bkg, double probCut, TString opt)
{
  double signal = sig->Integral(sig->FindBin(probCut),sig->GetNbinsX()+1);
  double backgr = bkg->Integral(sig->FindBin(probCut),sig->GetNbinsX()+1);
  double result = 0.;
  if(opt=="#frac{sig}{#sqrt{bkg}}")result = signal/TMath::Sqrt(backgr);
  if(opt=="#frac{sig}{#sqrt{sig+bkg}}")result = signal/TMath::Sqrt(signal+backgr);
  return result;
}

void drawArrow(const double xmin, const double y, const double xmax, const unsigned int color, const double lineWidth, const unsigned int lineStyle, double stretchFactor)
{
  // this function draws an arrow line with the chosen coordinates,
  drawLine(xmin, y, xmax, y, color, lineWidth, lineStyle);
  drawLine(0.8*xmax, y*stretchFactor, xmax, y, color, lineWidth, lineStyle);
  drawLine(0.8*xmax, y/stretchFactor, xmax, y, color, lineWidth, lineStyle);
}
