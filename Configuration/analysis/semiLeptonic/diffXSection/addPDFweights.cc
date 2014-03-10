// please set following paths
// LD_LIBRARY_PATH=/nfs/dust/cms/group/tophh/LHAPDF-605/lib:$LD_LIBRARY_PATH
// export LHAPATH=/nfs/dust/cms/group/tophh/LHAPDF-605/share/LHAPDF/
// execute macro by:
// make clean
// make AddPDFweights
// ./AddPDFweights

#include "addPDFweights.h"
#include "basicFunctions.h"

#include "TFile.h"
#include "TTree.h"

#include "iostream"

#include <boost/algorithm/string/replace.hpp>

#include "LHAPDF/LHAPDF.h"

AddPDFweights::AddPDFweights()
{
  // CT10nlo  CT10nnlo  HERAPDF15NNLO_VAR  MSTW2008nnlo68cl  NNPDF23_nnlo_as_0118  cteq6l1
  // Maximum  5 concurrent sets
  LHAPDF::initPDFSet(1, "cteq6l1"); // 44
  LHAPDF::initPDFSet(2, "CT10nlo"); // 52
  //LHAPDF::initPDFSet(2, "CT10nnlo"); // 50
  LHAPDF::initPDFSet(3, "MSTW2008nlo68cl"); // 40
  //LHAPDF::initPDFSet(3, "MSTW2008nnlo68cl"); // 40
  LHAPDF::initPDFSet(4, "NNPDF23_nlo_as_0118"); // 100
  //LHAPDF::initPDFSet(4, "NNPDF23_nnlo_as_0118"); // 100
  LHAPDF::initPDFSet(5, "HERAPDF15NLO_VAR"); // 10
  //LHAPDF::initPDFSet(5, "HERAPDF15NNLO_VAR"); // 10
  addPDFweights();
}

void AddPDFweights::addPDFweights()
{
  std::string samplePath (groupSpace);
  std::vector<std::string> inputSample;
  inputSample.push_back("combinedDiffXSecSigSummer12PF");
  //inputSample.push_back("muonDiffXSecSigSummer12PF");
  //inputSample.push_back("elecDiffXSecSigSummer12PF");

  for(unsigned iFile=0; iFile <inputSample.size(); ++iFile){

    TFile* fromFile = new TFile((samplePath+inputSample[iFile]+std::string(".root")).c_str(), "READ");
    TFile* toFile = new TFile((samplePath+std::string("CommonFiles/combinedDiffXSecPDFEnvelopeStudy.root")).c_str(), "RECREATE");

    std::vector<TString> folders_;
    folders_.push_back("analyzeTopPartonLevelKinematics");
    folders_.push_back("compositedHadronGenPhaseSpace");
    folders_.push_back("analyzeTopHadronLevelKinematicsLeptonPhaseSpace");
    folders_.push_back("analyzeTopHadronLevelKinematicsBjetsPhaseSpace");

    // container for xSec variables (histogram names)
    std::map<TString, std::vector<TString>> variables_;
    variables_[folders_[0]].insert(variables_[folders_[0]].end(), xSecVariablesKinFit, xSecVariablesKinFit + sizeof(xSecVariablesKinFit)/sizeof(TString));
    variables_[folders_[1]].push_back("Ngenjets");
    variables_[folders_[1]].push_back("rhosGen");
    variables_[folders_[2]].push_back("lepPtGen");
    variables_[folders_[2]].push_back("lepEtaGen");
    variables_[folders_[3]].push_back("bqPtGen");
    variables_[folders_[3]].push_back("bqEtaGen");
    variables_[folders_[3]].push_back("lbMassGen");
    variables_[folders_[3]].push_back("bbbarPtGen");
    variables_[folders_[3]].push_back("bbbarMassGen");

    // container for xSec variables (names in tree)
    std::map<TString, std::vector<TString>> treeVariables_;
    treeVariables_["topPtLead"].push_back("topPtLead");
    treeVariables_["topPtSubLead"].push_back("topPtSubLead");
    treeVariables_["topPtTtbarSys"].push_back("topPtTtbarSys");
    treeVariables_["topPt"].push_back("topPtLep");
    treeVariables_["topPt"].push_back("topPtHad");
    treeVariables_["topY"].push_back("topYLep");
    treeVariables_["topY"].push_back("topYHad");
    treeVariables_["ttbarPt"].push_back("ttbarPt");
    treeVariables_["ttbarY"].push_back("ttbarY");
    treeVariables_["ttbarMass"].push_back("ttbarMass");
    treeVariables_["ttbarDelPhi"].push_back("ttbarDelPhi");
    treeVariables_["ttbarPhiStar"].push_back("ttbarPhiStar");
    treeVariables_["Ngenjets"].push_back("NjetsTrue");
    treeVariables_["rhosGen"].push_back("rhosTrue");
    treeVariables_["lepPtGen"].push_back("lepPtGen");
    treeVariables_["lepEtaGen"].push_back("lepEtaGen");
    treeVariables_["bqPtGen"].push_back("bqPtGen");
    treeVariables_["bqPtGen"].push_back("bbarqPtGen");
    treeVariables_["bqEtaGen"].push_back("bqEtaGen");
    treeVariables_["bqEtaGen"].push_back("bbarqEtaGen");
    treeVariables_["lbMassGen"].push_back("lbMassGen");
    treeVariables_["bbbarPtGen"].push_back("bbbarPtGen");
    treeVariables_["bbbarMassGen"].push_back("bbbarMassGen");

    // container for pur variable names
    std::map<TString, TString> pureVariables_;
    pureVariables_["topPtLead"] = "topPtLead";
    pureVariables_["topPtSubLead"] = "topPtSubLead";
    pureVariables_["topPtTtbarSys"] = "topPtTtbarSys";
    pureVariables_["topPt"] = "topPt";
    pureVariables_["topY"] = "topY";
    pureVariables_["ttbarPt"] = "ttbarPt";
    pureVariables_["ttbarY"] = "ttbarY";
    pureVariables_["ttbarMass"] = "ttbarMass";
    pureVariables_["ttbarDelPhi"] = "ttbarDelPhi";
    pureVariables_["ttbarPhiStar"] = "ttbarPhiStar";
    pureVariables_["Ngenjets"] = "Njets";
    pureVariables_["rhosGen"] = "rhos";
    pureVariables_["lepPtGen"] = "lepPt";
    pureVariables_["lepEtaGen"] = "lepEta";
    pureVariables_["bqPtGen"] = "bqPt";
    pureVariables_["bqEtaGen"] = "bqEta";
    pureVariables_["lbMassGen"] = "lbMass";
    pureVariables_["bbbarPtGen"] = "bbbarPt";
    pureVariables_["bbbarMassGen"] = "bbbarMass";
    
    // container for pur variable names
    std::map<TString,std::vector<double>> binning_;
    binning_=makeVariableBinning();
    
    for(unsigned iFolder=0; iFolder < folders_.size(); ++iFolder){

      TTree* fromTree = (TTree*) fromFile->Get(folders_[iFolder]+"/tree");

      toFile->mkdir(folders_[iFolder]);
      toFile->cd(folders_[iFolder]);

      float Q;
      int id1;
      int id2;
      double x1;
      double x2;
      float weight;
      // container for values read from tree
      std::map< TString, float > value_;
      // container for histograms
      std::map<TString,std::vector<TH1F*>> hist_PDF;
      std::map<TString,TH1F*> hist_PDFdown;
      std::map<TString,TH1F*> hist_PDFup;
      std::map<TString,TCanvas*> canv_;

      fromTree->SetBranchAddress("Q",&Q);
      fromTree->SetBranchAddress("x1",&x1);
      fromTree->SetBranchAddress("x2",&x2);
      fromTree->SetBranchAddress("id1",&id1);
      fromTree->SetBranchAddress("id2",&id2);
      fromTree->SetBranchAddress("weight",&weight);
      for(unsigned iVar=0; iVar < variables_[folders_[iFolder]].size(); iVar++){
	for(unsigned iTreeVar=0; iTreeVar < treeVariables_[variables_[folders_[iFolder]][iVar]].size(); iTreeVar++){
	  value_[treeVariables_[variables_[folders_[iFolder]][iVar]][iTreeVar]];
	  fromTree->SetBranchAddress(treeVariables_[variables_[folders_[iFolder]][iVar]][iTreeVar],&value_[treeVariables_[variables_[folders_[iFolder]][iVar]][iTreeVar]]);
	}
	hist_PDF[variables_[folders_[iFolder]][iVar]].push_back( (TH1F*)fromFile->Get(folders_[iFolder]+"/"+variables_[folders_[iFolder]][iVar])->Clone() );
	for(unsigned i=1; i <=52; ++i) {
	  hist_PDF[variables_[folders_[iFolder]][iVar]].push_back( (TH1F*)fromFile->Get(folders_[iFolder]+"/"+variables_[folders_[iFolder]][iVar])->Clone() );
	  hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Reset();
	  TString name = variables_[folders_[iFolder]][iVar]+"_CT10nlo";
	  name+=i;
	  hist_PDF[variables_[folders_[iFolder]][iVar]][i]->SetName(name);
	}
	hist_PDFdown[variables_[folders_[iFolder]][iVar]] = (TH1F*)fromFile->Get(folders_[iFolder]+"/"+variables_[folders_[iFolder]][iVar])->Clone();
	  TString nameDown = variables_[folders_[iFolder]][iVar]+"PDFDown";
	  hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetName(nameDown);
	hist_PDFup[variables_[folders_[iFolder]][iVar]] = (TH1F*)fromFile->Get(folders_[iFolder]+"/"+variables_[folders_[iFolder]][iVar])->Clone();
	  TString nameUp = variables_[folders_[iFolder]][iVar]+"PDFUp";
	  hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetName(nameUp);
	  
	  canv_[variables_[folders_[iFolder]][iVar]] = new TCanvas(variables_[folders_[iFolder]][iVar]+"_canv",variables_[folders_[iFolder]][iVar]+"_canv",600,600);
      }

      //   for(unsigned i=1; i <=40; ++i) {
      //     hist_PDF[variables_[folders_[iFolder]][iVar]].push_back( (TH1F*)fromFile->Get("analyzeTopHadronLevelKinematicsBjetsPhaseSpace/bqPtGen")->Clone() );
      //     hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Reset();
      //     TString name = "bqPtGen_MSTW2008nlo68cl";
      //     name+=i;
      //     hist_PDF[variables_[folders_[iFolder]][iVar]][i]->SetName(name);
      //   }
      //   for(unsigned i=1; i <=100; ++i) {
      //     hist_PDF[variables_[folders_[iFolder]][iVar]].push_back( (TH1F*)fromFile->Get("analyzeTopHadronLevelKinematicsBjetsPhaseSpace/bqPtGen")->Clone() );
      //     hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Reset();
      //     TString name = "bqPtGen_NNPDF23_nlo_as_0118";
      //     name+=i;
      //     hist_PDF[variables_[folders_[iFolder]][iVar]][i]->SetName(name);
      //   }
      //   for(unsigned i=1; i <=10; ++i) {
      //     hist_PDF[variables_[folders_[iFolder]][iVar]].push_back( (TH1F*)fromFile->Get("analyzeTopHadronLevelKinematicsBjetsPhaseSpace/bqPtGen")->Clone() );
      //     hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Reset();
      //     TString name = "bqPtGen_HERAPDF15NLO_VAR";
      //     name+=i;
      //     hist_PDF[variables_[folders_[iFolder]][iVar]][i]->SetName(name);
      //   }

      for (int iEntry = 0, nEntries = fromTree->GetEntries(); iEntry < nEntries; ++iEntry) {
	fromTree->GetEntry(iEntry);
	//double xpdf1Central = LHAPDF::xfx(1, weightEvent->x1, weightEvent->Q, weightEvent->id1);
	//double xpdf2Central = LHAPDF::xfx(1, weightEvent->x2, weightEvent->Q, weightEvent->id2);
	//double w0Central = 1./(xpdf1Central * xpdf2Central);
	double w0 = 0.;
	for(unsigned i=0; i <=52; ++i) {
	  LHAPDF::usePDFMember(2,i);
	  double xpdf1 = LHAPDF::xfx(2, x1, Q, id1);
	  double xpdf2 = LHAPDF::xfx(2, x2, Q, id2);
	  if(i<1){
	    w0 = xpdf1 * xpdf2;
	  }
	  else{
	    double pdfWeight = ( ( ( (xpdf1 * xpdf2 / w0) - 1.0 ) / 1.645 ) + 1.0 ); // CT has 90% CL uncertainties, therefore they are scaled down to 68% CL
	    for(unsigned iVar=0; iVar < variables_[folders_[iFolder]].size(); iVar++){
	      //std::cout<<variables_[folders_[iFolder]][iVar]<<": "<<std::endl;
	      for(unsigned iTreeVar=0; iTreeVar < treeVariables_[variables_[folders_[iFolder]][iVar]].size(); iTreeVar++){
		//if(iEntry==0)std::cout<<treeVariables_[variables_[folders_[iFolder]][iVar]][iTreeVar]<<": "<<value_[treeVariables_[variables_[folders_[iFolder]][iVar]][iTreeVar]]<<", weight: "<<weight<<std::endl;
		hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Fill(value_[treeVariables_[variables_[folders_[iFolder]][iVar]][iTreeVar]],weight*pdfWeight);
	      }
	    }
	  }
	}
	//     for(unsigned i=1; i <=40; ++i) {
	//       LHAPDF::usePDFMember(3,i);
	//       double xpdf1 = LHAPDF::xfx(3, x1, Q, id1);
	//       double xpdf2 = LHAPDF::xfx(3, x2, Q, id2);
	//       //if(i<1)
	//       //  w0 = xpdf1 * xpdf2;
	//       //else
	//         pdfWeight.push_back(xpdf1 * xpdf2 / w0);
	//     }
	//     for(unsigned i=0; i <=100; ++i) {
	//       LHAPDF::usePDFMember(4,i);
	//       double xpdf1 = LHAPDF::xfx(4, x1, Q, id1);
	//       double xpdf2 = LHAPDF::xfx(4, x2, Q, id2);
	//       //if(i<1)
	//       //  w0 = xpdf1 * xpdf2;
	//       //else
	//         pdfWeight.push_back(xpdf1 * xpdf2 / w0);
	//     }
	//     for(unsigned i=1; i <=10; ++i) {
	//       LHAPDF::usePDFMember(5,i);
	//       double xpdf1 = LHAPDF::xfx(5, x1, Q, id1);
	//       double xpdf2 = LHAPDF::xfx(5, x2, Q, id2);
	//       //if(i<1)
	//       //  w0 = xpdf1 * xpdf2;
	//       //else
	//         pdfWeight.push_back(xpdf1 * xpdf2 / w0);
	//     }

	//std::cout << x1 << " " << id1 << " " << Q << " " << x2 << " " << id2 << std::endl;
      }
      
      for(unsigned iVar=0; iVar < variables_[folders_[iFolder]].size(); iVar++){
	for(unsigned i=1; i <=52; ++i) {
	  hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Scale(hist_PDF[variables_[folders_[iFolder]][iVar]][0]->Integral()/hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Integral());
	}
	unsigned nBins = hist_PDF[variables_[folders_[iFolder]][iVar]][0]->GetNbinsX();
	for(unsigned iBin=1; iBin <= nBins; ++iBin) {
	  float min = hist_PDF[variables_[folders_[iFolder]][iVar]][0]->GetBinContent(iBin);
	  float max = hist_PDF[variables_[folders_[iFolder]][iVar]][0]->GetBinContent(iBin);
	  for(unsigned i=1; i <=52; ++i) {
	    if(hist_PDF[variables_[folders_[iFolder]][iVar]][i]->GetBinContent(iBin)<min)min=hist_PDF[variables_[folders_[iFolder]][iVar]][i]->GetBinContent(iBin);
	    if(hist_PDF[variables_[folders_[iFolder]][iVar]][i]->GetBinContent(iBin)>max)max=hist_PDF[variables_[folders_[iFolder]][iVar]][i]->GetBinContent(iBin);
	  }
	  unsigned first = hist_PDF[variables_[folders_[iFolder]][iVar]][0]->FindBin(binning_[pureVariables_[variables_[folders_[iFolder]][iVar]]].front());
	  unsigned last = hist_PDF[variables_[folders_[iFolder]][iVar]][0]->FindBin(binning_[pureVariables_[variables_[folders_[iFolder]][iVar]]].back());
	  std::cout<<first<<", "<<last<<std::endl;
	  if(iBin<first){
	    hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetBinContent(iBin,min);
	    hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetBinContent(iBin,max);
	  }
	  else if(iBin>last){
	    hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetBinContent(iBin,max);
	    hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetBinContent(iBin,min);
	  }
	  else{
	    hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetBinContent(iBin,min+(max-min)*(iBin-first)/(last-first));
	    hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetBinContent(iBin,max-(max-min)*(iBin-first)/(last-first));
	  }
	}
	hist_PDFdown[variables_[folders_[iFolder]][iVar]]->Scale(hist_PDF[variables_[folders_[iFolder]][iVar]][0]->Integral()/hist_PDFdown[variables_[folders_[iFolder]][iVar]]->Integral());
	hist_PDFup[variables_[folders_[iFolder]][iVar]]->Scale(hist_PDF[variables_[folders_[iFolder]][iVar]][0]->Integral()/hist_PDFup[variables_[folders_[iFolder]][iVar]]->Integral());
	canv_[variables_[folders_[iFolder]][iVar]]->cd();
	hist_PDF[variables_[folders_[iFolder]][iVar]][0]->Draw();
	for(unsigned i=1; i <=52; ++i) {
	  hist_PDF[variables_[folders_[iFolder]][iVar]][i]->SetLineColor(kBlack);
	  hist_PDF[variables_[folders_[iFolder]][iVar]][i]->Draw("same");
	}
	hist_PDF[variables_[folders_[iFolder]][iVar]][0]->SetLineWidth(3);
	hist_PDF[variables_[folders_[iFolder]][iVar]][0]->SetLineStyle(2);
	hist_PDF[variables_[folders_[iFolder]][iVar]][0]->SetLineColor(kBlue);
	hist_PDF[variables_[folders_[iFolder]][iVar]][0]->Draw("same");
	hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetLineWidth(3);
	hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetLineStyle(2);
	hist_PDFdown[variables_[folders_[iFolder]][iVar]]->SetLineColor(kRed);
	hist_PDFdown[variables_[folders_[iFolder]][iVar]]->Draw("same");
	hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetLineWidth(3);
	hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetLineStyle(2);
	hist_PDFup[variables_[folders_[iFolder]][iVar]]->SetLineColor(kGreen-3);
	hist_PDFup[variables_[folders_[iFolder]][iVar]]->Draw("same");
	canv_[variables_[folders_[iFolder]][iVar]]->Write();
      }
    }
    toFile->Write();
    toFile->Close();
    fromFile->Close();
  }
}


