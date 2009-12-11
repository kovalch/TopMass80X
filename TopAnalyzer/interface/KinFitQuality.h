#ifndef KinFitQuality_h
#define KinFitQuality_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TopAnalysis/TopAnalyzer/interface/DoubleObject.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadronicEvent.h"

/**
   \class   KinFitQuality KinFitQuality.h "TopAnalysis/TopAnalyzer/interface/KinFitQuality.h"

   \brief   Derived class to analyze the kinematic fit hypothesis of fullhadronic ttbar events

   The structure keeps control histograms to analyze the relation between the fitted four vectors 
   of the kinematic fit hypothesis of fullhadronic ttbar events and the original input vectors 
   of the selected pat::Jets. These histograms are to be filled from the TtFullHadronicEvent and 
   from edm::View<pat::Jet>'s. The class is derived from the DoubleObject<CollectionA, CollectionB> 
   interface, which makes it usable in fwfull or fwlite.  This class is expected to be run with 
   the same jet collection from which the kinematic fit hypothesis was derived, otherwise it will 
   lead to wrong results!
*/

class KinFitQuality : public DoubleObject<TtFullHadronicEvent, const edm::View<pat::Jet> > {

 public:
  /// default constructor for fw lite
  explicit KinFitQuality();
  /// default constructor for fw full
  explicit KinFitQuality(const edm::ParameterSet& configFile);
  /// default destructor
  ~KinFitQuality(){};

  /**
     The following functions have to be implemented for any class
     derived from DoubleObject<CollectionA, CollectionB>
  **/
  /// histogramm booking for fwlite 
  void book();
  /// histogramm booking for fw
  void book(edm::Service<TFileService>& fileService);
  /// histogram filling interface for reconstruction level for access with fw or fwlite
  void fill(const TtFullHadronicEvent& tops, const edm::View<pat::Jet>& jets, const double& weight=1.);
  /// everything which needs to be done after the event loop
  void process(){};

 private:
  /// nothing to be done here for the moment
};

#endif
