#include "addPDFweights.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"

#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "ProgramOptionsReader.h"
#include "Helper.h"

#include "iostream"

#include <boost/algorithm/string/replace.hpp>

#include "LHAPDF/LHAPDF.h"

typedef ProgramOptionsReader po;

AddPDFweights::AddPDFweights()
{
  // CT10nlo  CT10nnlo  HERAPDF15NNLO_VAR  MSTW2008nnlo68cl  NNPDF23_nnlo_as_0118  cteq6l1
  // Maximum  5 concurrent sets
  //LHAPDF::initPDFSet(1, "cteq6l1"); // 44
  LHAPDF::initPDFSet( 1, "CT10nlo"); // 1+52
  LHAPDF::initPDFSet( 2, "CT10nlo_as_0117"); // 1
  LHAPDF::initPDFSet( 3, "CT10nlo_as_0119"); // 1
  //LHAPDF::initPDFSet(2, "CT10nnlo"); // 50
  LHAPDF::initPDFSet( 4, "MSTW2008nlo68cl"); // 1+40
  LHAPDF::initPDFSet( 5, "MSTW2008nlo68cl_asmz-68cl"); // 1(+40)
  LHAPDF::initPDFSet( 6, "MSTW2008nlo68cl_asmz+68cl"); // 1(+40)
  //LHAPDF::initPDFSet(3, "MSTW2008nnlo68cl"); // 40
  LHAPDF::initPDFSet( 7, "NNPDF23_nlo_as_0116"); // 1
  LHAPDF::initPDFSet( 8, "NNPDF23_nlo_as_0117"); // 4
  LHAPDF::initPDFSet( 9, "NNPDF23_nlo_as_0118"); // 12
  LHAPDF::initPDFSet(10, "NNPDF23_nlo_as_0119"); // 16
  LHAPDF::initPDFSet(11, "NNPDF23_nlo_as_0120"); // 12
  LHAPDF::initPDFSet(12, "NNPDF23_nlo_as_0121"); // 4
  LHAPDF::initPDFSet(13, "NNPDF23_nlo_as_0122"); // 1
  //LHAPDF::initPDFSet(4, "NNPDF23_nnlo_as_0118"); // 100
  //LHAPDF::initPDFSet(14, "HERAPDF15NLO_VAR"); // 10
  //LHAPDF::initPDFSet(5, "HERAPDF15NNLO_VAR"); // 10
  
  addPDFweights();
}

void AddPDFweights::addPDFweights()
{
  Helper* helper = new Helper();
  int channelID = Helper::channelID();
  
  std::string samplePath (po::GetOption<std::string>("analysisConfig.samplePath"));
  std::string inputSample(po::GetOption<std::string>("input"));
  std::string channel    (po::GetOption<std::string>("channel"));
  std::string treeDir;
  
  TChain* fromTree; int nFiles = 0;
  if (channelID == Helper::kAllJets) {
    fromTree = new TChain("analyzeKinFit/eventTree");
    nFiles = fromTree->Add((samplePath+inputSample+std::string(".root")).c_str());
    treeDir = "analyzeKinFit";
  }
  else {
	  fromTree = new TChain("analyzeHitFit/eventTree");
	  if (channelID != Helper::kElectronJets) nFiles += fromTree->Add((samplePath+inputSample+std::string("_muon/job_*.root")).c_str());
	  if (channelID != Helper::kMuonJets    ) nFiles += fromTree->Add((samplePath+inputSample+std::string("_electron/job_*.root")).c_str());
	  treeDir = "analyzeHitFit";
  }
  std::cout << "Adding " << nFiles << " files for " << inputSample << std::endl;
  
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
  
  TFile* toFile;
  if (channelID == Helper::kAllJets) {
    toFile = new TFile((samplePath+inputSample+std::string("_pdf.root")).c_str(), "RECREATE");
  }
  else {
    gSystem->mkdir((samplePath+inputSample+std::string("_pdf_")+channel).c_str());
    toFile = new TFile((samplePath+inputSample+std::string("_pdf_")+channel+std::string("/job_all.root")).c_str(), "RECREATE");
  }
  toFile->mkdir(treeDir.c_str());
  toFile->cd(treeDir.c_str());
  
  TTree* toTree = fromTree->CloneTree(0);
  toTree->SetBranchAddress("jet."   , &jetEvent);
  toTree->SetBranchAddress("top."   , &topEvent);
  toTree->SetBranchAddress("weight.", &weightEvent);
  
  for (int i = 0, nEntries = fromTree->GetEntries(); i < nEntries; ++i) {
    fromTree->GetEntry(i);
    double w0 = 0.;
    
    // CT10nlo 0--52
    for(unsigned j=0; j<=52; ++j) {
      LHAPDF::usePDFMember(1,j);
      double xpdf1 = LHAPDF::xfx(1, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(1, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      if(j<1)
        w0 = xpdf1 * xpdf2;
      else
        weightEvent->pdfWeight.push_back( ( ( ( (xpdf1 * xpdf2 / w0) - 1.0 ) / 1.645 ) + 1.0 ) ); // CT has 90% CL uncertainties, therefore they are scaled down to 68% CL
    }
    
    // CT10nlo_as_0117 53
    {
      LHAPDF::usePDFMember(2,0);
      double xpdf1 = LHAPDF::xfx(2, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(2, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back( ( ( ( (xpdf1 * xpdf2 / w0) - 1.0 ) / (5./6.) ) + 1.0 ) ); // alpha_s uncertainty scaled to 68% CL
    }
    
    // CT10nlo_as_0119 54
    {
      LHAPDF::usePDFMember(3,0);
      double xpdf1 = LHAPDF::xfx(3, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(3, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back( ( ( ( (xpdf1 * xpdf2 / w0) - 1.0 ) / (5./6.) ) + 1.0 ) ); // alpha_s uncertainty scaled to 68% CL
    }
    
    // MSTW2008nlo68cl 55--95
    for(unsigned j=0; j<=40; ++j) {
      LHAPDF::usePDFMember(4,j);
      double xpdf1 = LHAPDF::xfx(4, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(4, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // MSTW2008nlo68cl_asmz-68cl 96
    {
      LHAPDF::usePDFMember(5,0);
      double xpdf1 = LHAPDF::xfx(5, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(5, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back( ( ( ( (xpdf1 * xpdf2 / w0) - 1.0 ) / (5./4.) ) + 1.0 ) ); // alpha_s down uncertainty scaled to 68% CL
    }
    
    // MSTW2008nlo68cl_asmz+68cl 97
    {
      LHAPDF::usePDFMember(6,0);
      double xpdf1 = LHAPDF::xfx(6, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(6, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0116 98
    for(unsigned j=0; j<1; ++j) {
      LHAPDF::usePDFMember(7,j);
      double xpdf1 = LHAPDF::xfx(7, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(7, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0117 99--102
    for(unsigned j=0; j<4; ++j) {
      LHAPDF::usePDFMember(8,j);
      double xpdf1 = LHAPDF::xfx(8, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(8, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0118 103--114
    for(unsigned j=0; j<12; ++j) {
      LHAPDF::usePDFMember(9,j);
      double xpdf1 = LHAPDF::xfx(9, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(9, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0119 115--130
    for(unsigned j=0; j<16; ++j) {
      LHAPDF::usePDFMember(10,j);
      double xpdf1 = LHAPDF::xfx(10, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(10, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0120 131--142
    for(unsigned j=0; j<12; ++j) {
      LHAPDF::usePDFMember(11,j);
      double xpdf1 = LHAPDF::xfx(11, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(11, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0121 143--146
    for(unsigned j=0; j<4; ++j) {
      LHAPDF::usePDFMember(12,j);
      double xpdf1 = LHAPDF::xfx(12, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(12, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    // NNPDF23_nlo_as_0122 147
    for(unsigned j=0; j<1; ++j) {
      LHAPDF::usePDFMember(13,j);
      double xpdf1 = LHAPDF::xfx(13, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(13, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    
    
    /* HERAPDF
    for(unsigned j=0; j<=10; ++j) {
      LHAPDF::usePDFMember(5,j);
      double xpdf1 = LHAPDF::xfx(5, weightEvent->x1, weightEvent->Q, weightEvent->id1);
      double xpdf2 = LHAPDF::xfx(5, weightEvent->x2, weightEvent->Q, weightEvent->id2);
      weightEvent->pdfWeight.push_back(xpdf1 * xpdf2 / w0);
    }
    */

    //std::cout << weightEvent->x1 << " " << weightEvent->id1 << " " << weightEvent->Q << " " << weightEvent->x2 << " " << weightEvent->id2 << std::endl;
    toTree->Fill();
  }
  
  toFile->Write();
  toFile->Close();
}


