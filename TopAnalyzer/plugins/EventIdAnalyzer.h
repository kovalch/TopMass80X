#ifndef EventIdAnalyzer_h  
#define EventIdAnalyzer_h

#include "TH1F.h"
#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>
#include <fstream>
#include <string>

/// This analyzer writes the event ID (i.e. run, lumi and event numbers) into a txt file

class EventIdAnalyzer : public edm::EDAnalyzer {

 public:
  
  explicit EventIdAnalyzer(const edm::ParameterSet&);
  ~EventIdAnalyzer();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
           	
  edm::InputTag input_;
  bool verbose_;
  std::string textFileName_;
  edm::InputTag weight_;
  
  /// output histogram
  std::map<std::string, TH1F*> hists_;


  
  std::ofstream myfile;
  
  TTree * tree;
  unsigned int run;
  unsigned int lumi;
  unsigned int event;

};  

#endif  
