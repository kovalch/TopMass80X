#ifndef JetEnergyResolutionBiasAnalyzer_h
#define JetEnergyResolutionBiasAnalyzer_h

#include <map>

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TLorentzVector.h"
#include "TRandom3.h"

/**
   \class   JetEnergyResolutionBiasAnalyzer JetEnergyResolutionBiasAnalyzer.h "TopAnalysis/TopAnalyzer/interface/JetEnergyResolutionBiasAnalyzer.h"
*/

class TH1;
class TH3;

class JetEnergyResolutionBiasAnalyzer : public edm::EDAnalyzer {

public:
  /// default constructor
  explicit JetEnergyResolutionBiasAnalyzer(const edm::ParameterSet&);
  /// default destructor
  ~JetEnergyResolutionBiasAnalyzer(){};

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void smearEnergy(TLorentzVector&);
  
protected:
  /// histogram container
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TH3*> hists3d_;
  /// random number generator object
  TRandom3 *randNumGen_;
  /// parameters read from config file
  double etaCut_ ;
  double resFact_;
  double resParS_;
  double resParN_;
  double resParC_;
};

#endif
