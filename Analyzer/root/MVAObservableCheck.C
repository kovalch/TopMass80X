#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <Math/Boost.h>
#include <../../../TopAnalysis/TopAnalyzer/interface/TopAngles.h>

#include "tdrstyle.C"

TopAngles::TopAngles(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& inVector1,
		     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& inVector2,
		     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& inVector3,
		     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& inVector4,
		     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& inVector5,
		     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& inVector6):
  inputVector1(inVector1),
  inputVector2(inVector2),
  inputVector3(inVector3),
  inputVector4(inVector4),
  inputVector5(inVector5),
  inputVector6(inVector6),
  W1(inVector2+inVector3),
  W2(inVector5+inVector6),
  top1(inVector1+inVector2+inVector3),
  top2(inVector4+inVector5+inVector6),
  ttbar(inVector1+inVector2+inVector3+inVector4+inVector5+inVector6),
  CoMBoostTtbar(ttbar.BoostToCM()),
  CoMBoostTop1 (top1.BoostToCM()),
  CoMBoostTop2 (top2.BoostToCM()),
  CoMBoostW1   (W1.BoostToCM()),
  CoMBoostW2   (W2.BoostToCM())
{
}

  /// define getter functions (different for all decay channels)
    /// define common getter functions for all decay channels
    double TopAngles::ttDetFrame()   const{ return ROOT::Math::VectorUtil::Angle(top1,top2);}
    double TopAngles::bbDetFrame()   const{ return ROOT::Math::VectorUtil::Angle(inputVector1,inputVector4);}
    double TopAngles::bbTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector1), CoMBoostTtbar(inputVector4));}
    double TopAngles::WWTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(W1),CoMBoostTtbar(W2));}
    double TopAngles::tBBranch1TtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(top1),CoMBoostTtbar(inputVector1));}
    double TopAngles::tBBranch2TtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(top2),CoMBoostTtbar(inputVector4));}
    double TopAngles::bWBranch1TtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(W1)  ,CoMBoostTtbar(inputVector1));}
    double TopAngles::bWBranch2TtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(W2)  ,CoMBoostTtbar(inputVector4));}
    double TopAngles::tWBranch1TopInTtbarFrameWInTopFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(top1),CoMBoostTop1(W1));}
    double TopAngles::tWBranch2TopInTtbarFrameWInTopFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(top2),CoMBoostTop2(W2));}
    /// define getter functions for the semileptonic decay channel
    double TopAngles::qQbarTopFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTop1(inputVector2),CoMBoostTop1(inputVector3));}
    double TopAngles::qQbarDetFrame() const{ return ROOT::Math::VectorUtil::Angle(inputVector2,inputVector3);}
    double TopAngles::blepQTtbarFrame()    const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector4),CoMBoostTtbar(inputVector2));}
    double TopAngles::blepQbarTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector4),CoMBoostTtbar(inputVector3));}
    double TopAngles::bhadQTopFrame()    const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTop1(inputVector1),CoMBoostTop1(inputVector2));}
    double TopAngles::bhadQbarTopFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTop1(inputVector1),CoMBoostTop1(inputVector3));}
    double TopAngles::lepBlepTopFrame()  const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTop2(inputVector5),CoMBoostTop2(inputVector4));}
    double TopAngles::lepBlepTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector5),CoMBoostTtbar(inputVector4));}
    double TopAngles::lepBhadTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector5),CoMBoostTtbar(inputVector1));}
    double TopAngles::lepQTtbarFrame()    const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector5),CoMBoostTtbar(inputVector2));}
    double TopAngles::lepQbarTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector5),CoMBoostTtbar(inputVector3));}
    double TopAngles::lepNuTopFrame()    const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTop2(inputVector5),CoMBoostTop2(inputVector6));}
    double TopAngles::nuBlepTopFrame()   const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTop2(inputVector6),CoMBoostTop2(inputVector4));}
    double TopAngles::nuBhadTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector6),CoMBoostTtbar(inputVector1));}
    double TopAngles::nuQTtbarFrame()    const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector6),CoMBoostTtbar(inputVector2));}
    double TopAngles::nuQBarTtbarFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostTtbar(inputVector6),CoMBoostTtbar(inputVector3));}
    //double TopAngles::lepWlepLepInWFrameWInDetFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostW2(inputVector5), W2);}
    double TopAngles::lepWlepLepInWFrameWInDetFrame() const{ return ROOT::Math::VectorUtil::Angle(CoMBoostW1(inputVector1), CoMBoostW1(inputVector2));}
    double TopAngles::nuWlepNuInWFrameWInDetFrame()   const{ return ROOT::Math::VectorUtil::Angle(CoMBoostW2(inputVector6), W2);}
    double TopAngles::qWhadQInWFrameWInDetFrame()     const{ return ROOT::Math::VectorUtil::Angle(CoMBoostW1(inputVector2), W1);}
    double TopAngles::qbarWhadQInWFrameWInDetFrame()  const{ return ROOT::Math::VectorUtil::Angle(CoMBoostW1(inputVector3), W1);}

void MVAObservableCheck() {
  //*
  TFile* file = new TFile("/scratch/hh/current/cms/user/mseidel/Fall11_TTJets1845_1.00_muon/analyzeTop.root", "READ");
  TTree* fileTree = (TTree*) file->Get("analyzeHitFit/eventTree");
  TFile *f = new TFile("temp1845.root", "RECREATE");
  TTree* eventTree = fileTree->CopyTree("leptonPt>30 & bottomSSVJetMultiplicity > 1 & hitFitProb>0.1");
  //*/
  
  /*
  TFile *f = new TFile("temp1725.root", "READ");
  TTree* eventTree = (TTree*) f->Get("eventTree");
  //*/
  
  int event;
  eventTree->SetBranchAddress("event", &event);
  
  int combi;
  eventTree->SetBranchAddress("combi", &combi);
  
  double leptonPt;
  eventTree->SetBranchAddress("leptonPt", &leptonPt);
  
  double hadQBSSV;
  eventTree->SetBranchAddress("hadQBSSV", &hadQBSSV);
  
  double hadQBarBSSV;
  eventTree->SetBranchAddress("hadQBarBSSV", &hadQBarBSSV);
  
  double hadBBSSV;
  eventTree->SetBranchAddress("hadBBSSV", &hadBBSSV);
  
  double lepBBSSV;
  eventTree->SetBranchAddress("lepBBSSV", &lepBBSSV);
  
  double hitFitProb;
  eventTree->SetBranchAddress("hitFitProb", &hitFitProb);
  
  int target;
  eventTree->SetBranchAddress("target", &target);
  
  double hadTopMass;
  eventTree->SetBranchAddress("hadTopMass", &hadTopMass);
  
  TLorentzVector* hadTop  = new TLorentzVector(); 
  TLorentzVector* hadB    = new TLorentzVector(); 
  TLorentzVector* hadW    = new TLorentzVector(); 
  TLorentzVector* hadQ    = new TLorentzVector(); 
  TLorentzVector* hadQBar = new TLorentzVector(); 
  TLorentzVector* lepTop  = new TLorentzVector(); 
  TLorentzVector* lepB    = new TLorentzVector(); 
  TLorentzVector* lepW    = new TLorentzVector(); 
  TLorentzVector* lepton  = new TLorentzVector(); 
  TLorentzVector* nu      = new TLorentzVector();
  eventTree->SetBranchAddress("hadTop", &hadTop);
  eventTree->SetBranchAddress("hadB.", &hadB);
  eventTree->SetBranchAddress("hadW.", &hadW);
  eventTree->SetBranchAddress("hadQ.", &hadQ);
  eventTree->SetBranchAddress("hadQBar.", &hadQBar);
  eventTree->SetBranchAddress("lepTop.", &lepTop);
  eventTree->SetBranchAddress("lepB.", &lepB);
  eventTree->SetBranchAddress("lepW.", &lepW);
  eventTree->SetBranchAddress("lepton.", &lepton);
  eventTree->SetBranchAddress("nu.", &nu);
  
  int entries = eventTree->GetEntries();
  TH1F* hSig = new TH1F("hSig", "hSig", 80, 0, 4);
  hSig->SetLineColor(kRed+1);
  TH1F* hBkg = new TH1F("hBkg", "hBkg", 80, 0, 4);
  hBkg->SetLineColor(kBlue+1);
  
  TH1F* hMass = new TH1F("hMass", "hMass", 70, 50, 400);
  TH2F* hCorSig = new TH2F("hCorSig", "hCorSig", 400, 0, 4, 300, 100, 250);
  hCorSig->SetMarkerColor(kRed+1);
  TH2F* hCorBkg = new TH2F("hCorBkg", "hCorBkg", 400, 0, 4, 300, 100, 250);
  hCorBkg->SetMarkerColor(kBlue+1);
  
  for (int i = 0; i < entries; i++) {
    if (i%10000 == 0) std::cout << entries - i << "... " << std::endl;
    
    eventTree->GetEntry(i);
    
    if (leptonPt > 30
        && hadQBSSV < 1.74 && hadQBarBSSV < 1.74 && hadBBSSV > 1.74 && lepBBSSV > 1.74 
        && hitFitProb > 0.1) {
      //std::cout << lepton->Angle(nu->Vect()) << std::endl;
      
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > hadTop_(hadTop->Px(), hadTop->Py(), hadTop->Pz(), hadTop->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lepTop_(lepTop->Px(), lepTop->Py(), lepTop->Pz(), lepTop->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > hadW_(hadW->Px(), hadW->Py(), hadW->Pz(), hadW->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lepW_(lepW->Px(), lepW->Py(), lepW->Pz(), lepW->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > hadB_(hadB->Px(), hadB->Py(), hadB->Pz(), hadB->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > hadQ_(hadQ->Px(), hadQ->Py(), hadQ->Pz(), hadQ->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > hadQBar_(hadQBar->Px(), hadQBar->Py(), hadQBar->Pz(), hadQBar->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lepB_(lepB->Px(), lepB->Py(), lepB->Pz(), lepB->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lepton_(lepton->Px(), lepton->Py(), lepton->Pz(), lepton->E());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > nu_(nu->Px(), nu->Py(), nu->Pz(), nu->E());
      
      //*
      TopAngles* angles = new TopAngles::TopAngles(hadB_, hadQ_, hadQBar_, lepB_, lepton_, nu_);
      double angle = angles->bWBranch1TtbarFrame();
      double angle2 = angles->bWBranch2TtbarFrame();
      //*/
      
      /*
      const ROOT::Math::Boost CoMBoostHadTop(hadTop_.BoostToCM());
      const ROOT::Math::Boost CoMBoostLepTop(lepTop_.BoostToCM());
      
      double amwtHadQ    = TMath::Max(0., 4*hadTop_.M()
                  *CoMBoostHadTop(hadQ_).E()*(hadTop_.M2()-hadB_.M2()
                    -2*hadTop_.M()*CoMBoostHadTop(hadQ_).E())
                  /(TMath::Power(hadTop_.M2()-hadB_.M2(),2)+hadW_.M2()*(hadTop_.M2()
                    -hadB_.M2())-TMath::Power(hadW_.M2(),2)));
      double amwtHadQBar = TMath::Max(0., 4*hadTop_.M()
                  *CoMBoostHadTop(hadQBar_).E()*(hadTop_.M2()-hadB_.M2()
                    -2*hadTop_.M()*CoMBoostHadTop(hadQBar_).E())
                  /(TMath::Power(hadTop_.M2()-hadB_.M2(),2)+hadW_.M2()*(hadTop_.M2()
                    -hadB_.M2())-TMath::Power(hadW_.M2(),2)));
      double amwtLepton  = TMath::Max(0., 4*lepTop_.M()
                  *CoMBoostLepTop(lepton_).E()*(lepTop_.M2()-lepB_.M2()
                    -2*lepTop_.M()*CoMBoostLepTop(lepton_).E())
                  /(TMath::Power(lepTop_.M2()-lepB_.M2(),2)+lepW_.M2()*(lepTop_.M2()
                    -lepB_.M2())-TMath::Power(lepW_.M2(),2)));
      */
      
      /*
      ttDetFrame, tBBranch1TtbarFrame, tBBranch2TtbarFrame
      bWBranch2TtbarFrame, bWBranch1TtbarFrame o unmatched korreliert
      qQbarTopFrame, lepNuTopFrame + korreliert
      Angle(CoMBoostW2(inputVector5), CoMBoostTop2(inputVector5))
      Angle(CoMBoostW1(inputVector2), CoMBoostTop1(inputVector2))
      */
      
      if (target == 1) {
        hSig->Fill(angle);
        hCorSig->Fill(angle, hadTopMass);
      }
      else {
        hBkg->Fill(angle);
        hCorBkg->Fill(angle, hadTopMass);
      }
      
      //if (angle < 2 & angle2 < 2)
      hMass->Fill(hadTopMass);
    }
  }
  
  /*
  file->Close();
  f->Write();
  //*/
  
  std::cout << hSig->Integral() << std::endl;
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
  canvas->cd();
  
  hSig->Scale(1./hSig->Integral());
  hBkg->Scale(1./hBkg->Integral());
  
  hSig->Draw("");
  hBkg->Draw("same");
  
  TCanvas* canvasMass = new TCanvas("canvasMass", "canvasMass", 500, 500);
  canvasMass->cd();
  hMass->Fit("gaus");
  
  TCanvas* canvasCor = new TCanvas("canvasCor", "canvasCor", 500, 500);
  canvasCor->cd();
  hCorBkg->Draw();
  hCorSig->Draw("SAME");
}
