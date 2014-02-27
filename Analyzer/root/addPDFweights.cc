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
  // Maximum  3 concurrent set(s) ...
  LHAPDF::initPDFSet(1, "cteq66.LHgrid"); // 44
  LHAPDF::initPDFSet(2, "MSTW2008lo68cl.LHgrid"); // 40
  LHAPDF::initPDFSet(3, "NNPDF23_nnlo_as_0118.LHgrid"); // 100
  //LHAPDF::initPDFSet(4, "HERAPDF15NNLO_VAR.LHgrid"); // 10
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

    double w0 = 0.;
    for(unsigned i=0; i <=44; ++i) {
      LHAPDF::usePDFMember(1,i);
      double xpdf1 = LHAPDF::xfx(1, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(1, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      if(i<1)
        w0 = xpdf1 * xpdf2;
      else
        weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    for(unsigned i=0; i <=40; ++i) {
      LHAPDF::usePDFMember(2,i);
      double xpdf1 = LHAPDF::xfx(2, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(2, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      if(i<1)
        w0 = xpdf1 * xpdf2;
      else
        weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    for(unsigned i=0; i <=100; ++i) {
      LHAPDF::usePDFMember(3,i);
      double xpdf1 = LHAPDF::xfx(3, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(3, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      if(i<1)
        w0 = xpdf1 * xpdf2;
      else
        weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }

    //std::cout << weightEvent->x1 << " " << weightEvent->id1 << " " << weightEvent->Q << " " << weightEvent->x2 << " " << weightEvent->id2 << std::endl;
    toTree->Fill();
  }
  
  toFile->Write();
  toFile->Close();
  fromFile->Close();
}


