#include "addPDFweights.h"

#include "TFile.h"
#include "TTree.h"

#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "ProgramOptionsReader.h"

#include "iostream"

#include <boost/algorithm/string/replace.hpp>

#include "LHAPDF/LHAPDF.h"

typedef ProgramOptionsReader po;

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
  std::string samplePath (po::GetOption<std::string>("analysisConfig.samplePath"));
  std::string inputSample(po::GetOption<std::string>("input"));

  TFile* fromFile = new TFile((samplePath+inputSample+std::string(".root")).c_str(), "READ");
  TTree* fromTree = (TTree*) fromFile->Get("analyzeKinFit/eventTree");
  
  fromTree->SetBranchStatus("*", 0);
  fromTree->SetBranchStatus("jet.*"   , 1);
  fromTree->SetBranchStatus("top.*"   , 1);
  fromTree->SetBranchStatus("weight.*", 1);

  JetEvent    *jetEvent    = new JetEvent();
  TopEvent    *topEvent    = new TopEvent();
  WeightEvent *weightEvent = new WeightEvent();

  fromTree->SetBranchAddress("jet."   , &jetEvent);
  fromTree->SetBranchAddress("top."   , &topEvent);
  fromTree->SetBranchAddress("weight.", &weightEvent);
  
  TFile* toFile = new TFile((samplePath+inputSample+std::string("_pdf.root")).c_str(), "RECREATE");
  toFile->mkdir("analyzeKinFit");
  toFile->cd("analyzeKinFit");
  
  TTree* toTree = fromTree->CloneTree(0);
  toTree->SetBranchAddress("jet."   , &jetEvent);
  toTree->SetBranchAddress("top."   , &topEvent);
  toTree->SetBranchAddress("weight.", &weightEvent);
  
  for (int i = 0, nEntries = fromTree->GetEntries(); i < nEntries; ++i) {
    fromTree->GetEntry(i);

    //double xpdf1Central = LHAPDF::xfx(1, weightEvent->x1, weightEvent->Q, weightEvent->id1);
    //double xpdf2Central = LHAPDF::xfx(1, weightEvent->x2, weightEvent->Q, weightEvent->id2);
    //double w0Central = 1./(xpdf1Central * xpdf2Central);
    double w0 = 0.;
    for(unsigned i=0; i <=52; ++i) {
      LHAPDF::usePDFMember(2,i);
      double xpdf1 = LHAPDF::xfx(2, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(2, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      if(i<1)
        w0 = xpdf1 * xpdf2;
      else
        weightEvent->pdfWeight.push_back( ( ( ( (xpdf1 * xpdf2 / w0) - 1.0 ) / 1.645 ) + 1.0 ) ); // CT has 90% CL uncertainties, therefore they are scaled down to 68% CL
    }
    for(unsigned i=1; i <=40; ++i) {
      LHAPDF::usePDFMember(3,i);
      double xpdf1 = LHAPDF::xfx(3, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(3, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      //if(i<1)
      //  w0 = xpdf1 * xpdf2;
      //else
        weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    for(unsigned i=0; i <=100; ++i) {
      LHAPDF::usePDFMember(4,i);
      double xpdf1 = LHAPDF::xfx(4, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(4, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      //if(i<1)
      //  w0 = xpdf1 * xpdf2;
      //else
        weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    for(unsigned i=1; i <=10; ++i) {
      LHAPDF::usePDFMember(5,i);
      double xpdf1 = LHAPDF::xfx(5, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(5, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      //if(i<1)
      //  w0 = xpdf1 * xpdf2;
      //else
        weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }

    //std::cout << weightEvent->x1 << " " << weightEvent->id1 << " " << weightEvent->Q << " " << weightEvent->x2 << " " << weightEvent->id2 << std::endl;
    toTree->Fill();
  }
  
  toFile->Write();
  toFile->Close();
  fromFile->Close();
}


