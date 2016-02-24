#ifndef JetEventAnalyzer_h
#define JetEventAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TopMass/TopEventTree/interface/JetEvent.h"

//#include <utility>

class JetEventAnalyzer : public edm::EDAnalyzer {

  public:

  explicit JetEventAnalyzer(const edm::ParameterSet&);
  ~JetEventAnalyzer();

 private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  std::pair<TVector2, TVector2> getPullVector( std::vector<pat::Jet>::const_iterator patJet );
  std::pair<TVector2, TVector2> getGenPullVector( std::vector<pat::Jet>::const_iterator patJet );

  edm::Service<TreeRegistryService> trs;

  edm::EDGetTokenT< std::vector<pat::Jet > > jets_;
  edm::InputTag alternativeJets_;
  edm::EDGetTokenT<std::vector<pat::MET> > met_;
  std::string gluonTagName_;

  // max possible number of jets in events
  const int kJetMAX_;

  // THE JetEvent to store the information
  JetEvent* jet;

  // check only once per module run if the needed collections are available
  bool checkedIsPFJet, checkedJERSF, checkedJESSF, checkedTotalSF, checkedQGTag, checkedBReg;
  bool        isPFJet,     hasJERSF,     hasJESSF,     hasTotalSF,     hasQGTag,     hasBReg;
  bool alternativeJetsAvailable;

};

#endif
