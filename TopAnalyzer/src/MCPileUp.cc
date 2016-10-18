#ifndef MCPileUp_cc_
#define MCPileUp_cc_

#include "TopAnalysis/TopAnalyzer/interface/MCPileUp.h"

// =============================
//  constructors and destructor
// =============================

MCPileUp::MCPileUp(const edm::ParameterSet& iConfig)
{
  inTag_PUSource = iConfig.getParameter<edm::InputTag>("PUSource");  
  mayConsume<edm::View<PileupSummaryInfo>>(inTag_PUSource);
}

MCPileUp::~MCPileUp(){}

// =============================

void MCPileUp::beginJob()
{
  edm::Service<TFileService> fs;

  histoNPUEvents = fs->make<TH1F>("pileup","Number of PU Events",71,-0.5,70.5);
  histoNPUEvents -> SetTitle("");
  histoNPUEvents -> GetXaxis() -> SetTitle("Number of PU Events");
  histoNPUEvents -> GetYaxis() -> SetTitle("Frequency");
  histoNPUEvents -> SetFillColor(2);

  histoNPUEvents_true = fs->make<TH1F>("pileup_true","Number of PU Events (true)",71,-0.5,70.5);
  histoNPUEvents_true -> SetTitle("");
  histoNPUEvents_true -> GetXaxis() -> SetTitle("Number of PU Events (true)");
  histoNPUEvents_true -> GetYaxis() -> SetTitle("Frequency");
  histoNPUEvents_true -> SetFillColor(2);

  histoNPUEvents_LessBins = fs->make<TH1F>("pileup_LessBins","Number of PU Events",51,-0.5,50.5);
  histoNPUEvents_LessBins -> SetTitle("");
  histoNPUEvents_LessBins -> GetXaxis() -> SetTitle("Number of PU Events");
  histoNPUEvents_LessBins -> GetYaxis() -> SetTitle("Frequency");
  histoNPUEvents_LessBins -> SetFillColor(2);

  histoNPUEvents_LessBins_true = fs->make<TH1F>("pileup_LessBins_true","Number of PU Events (true)",51,-0.5,50.5);
  histoNPUEvents_LessBins_true -> SetTitle("");
  histoNPUEvents_LessBins_true -> GetXaxis() -> SetTitle("Number of PU Events (true)");
  histoNPUEvents_LessBins_true -> GetYaxis() -> SetTitle("Frequency");
  histoNPUEvents_LessBins_true -> SetFillColor(2);
}

// =============================

void MCPileUp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<PileupSummaryInfo> > pPUInfo;
   iEvent.getByLabel(inTag_PUSource, pPUInfo);

   edm::View<PileupSummaryInfo>::const_iterator iterPU;

   for(iterPU = pPUInfo->begin(); iterPU != pPUInfo->end(); ++iterPU)
   { 
     if (iterPU->getBunchCrossing()==0) // -1: previous BX, 0: current BX,  1: next BX
     {
       histoNPUEvents          -> Fill(iterPU->getPU_NumInteractions());
       histoNPUEvents_LessBins -> Fill(iterPU->getPU_NumInteractions());

       histoNPUEvents_true          -> Fill(iterPU->getTrueNumInteractions());
       histoNPUEvents_LessBins_true -> Fill(iterPU->getTrueNumInteractions());
     }
   }
}

// =============================

void MCPileUp::endJob(){}

#endif
