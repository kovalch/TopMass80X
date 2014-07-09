#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include <iostream>
#include <string>

void optimizeSelectionAllJets()
{
  TFile* fileD = TFile::Open("/nfs/dust/cms/user/eschliec/TopMass/2012/Skim_04/MJP12_PTRESI_v1_data.root");
  TFile* fileM = TFile::Open("/nfs/dust/cms/user/eschliec/TopMass/2012/Skim_04/Z2_S12_ABS_JES_100_172_5_MadSpin_sig.root");

  TTree* treeD = (TTree*)fileD->Get("analyzeKinFit/eventTree");
  TTree* treeM = (TTree*)fileM->Get("analyzeKinFit/eventTree");

  double startI = 0.1;
  double stepI  = 0.1;
  double nI     = 9;

  double startJ = 2.0;
  double stepJ  = 0.2;
  double nJ     = 10;

  std::string selBase = "&& max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2) < 6 && jet.alternativeJet[3].Pt() > 60 && jet.jet[3].Pt() > 60)";

  TH2F* mapfSig = new TH2F("mapfSig",";P(#chi^{2});#Delta R_{bb}",nI,startI,startI+stepI*nI,nJ,startJ,startJ+stepJ*nJ);
  TH2F* mapfCP  = (TH2F*)mapfSig->Clone("mapfCP");
  TH2F* mapNcp  = (TH2F*)mapfSig->Clone("mapNcp");

  for(int i=0; i<nI; ++i){
    for(int j=0; j<nJ; ++j){

      double valI = startI + i*stepI;
      char bufferI [19];
      sprintf(bufferI, "(top.fitProb > %4.2f", valI);

      double valJ = startJ + j*stepJ;
      char bufferJ [118];
      sprintf(bufferJ, "&& (pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2)) > %4.2f*%4.2f", valJ, valJ);
      /*
      std::cout << std::string(bufferI,19) << "|" << std::endl;
      std::cout << std::string(bufferJ,118) << "|" << std::endl;
      std::cout << selBase << std::endl;

      std::cout << (std::string(bufferI,19)+std::string(bufferJ,118)+selBase).c_str() << std::endl;
      std::cout << (std::string("weight.combinedWeight*")+std::string(bufferI,19)+std::string(bufferJ,118)+selBase).c_str() << std::endl;
      */
      treeD->Draw("top.combinationType>>hComboD(20,-10,10)",(std::string(bufferI,19)+std::string(bufferJ,118)+selBase).c_str(),"goff");
      treeM->Draw("top.combinationType>>hComboM(20,-10,10)",(std::string("weight.combinedWeight*( (top.combinationType!=0) ? sqrt(exp(0.156-0.00137*top.genpartonTop1.Pt())*exp(0.156-0.00137*top.genpartonTop2.Pt())) : 1.0)*")+std::string(bufferI,19)+std::string(bufferJ,118)+selBase).c_str(),"goff");

      TH1F* hD = (TH1F*)gDirectory->Get("hComboD");
      TH1F* hM = (TH1F*)gDirectory->Get("hComboM");

      mapfSig->SetBinContent(i+1,j+1,hM->Integral(0,21)/hD->GetBinContent(11));
      mapfCP ->SetBinContent(i+1,j+1,hM->GetBinContent(12)/hM->Integral(0,21));
      mapNcp ->SetBinContent(i+1,j+1,1./sqrt(hM->GetBinContent(12)));
    }
  }

  TCanvas* canv1 = new TCanvas("canv1","fSig",600,600);
  canv1->cd();
  mapfSig->Draw("colz text");
  TCanvas* canv2 = new TCanvas("canv2","fCP",600,600);
  canv2->cd();
  mapfCP->Draw("colz text");
  TCanvas* canv3 = new TCanvas("canv3","Ncp",600,600);
  canv3->cd();
  mapNcp->Draw("colz text");
}
