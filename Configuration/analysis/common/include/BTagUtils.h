#ifndef BTagUtils_h
#define BTagUtils_h

#include <string>
#include <vector>


class TString;
class TSelectorList;

#include "classesFwd.h"
#include "TopAnalysis/ZTopUtils/interface/bTagBase.h"





class BTagSFGeneric : public ztop::bTagBase{

public:

  /// Constructor
  BTagSFGeneric(const char* btagEfficiencyInputDir,
                const char* btagEfficiencyOutputDir,
                const std::vector<std::string>& channels,
                TString systematic);

  /// Destructor
  ~BTagSFGeneric(){}



  /// Whether to produce b-tag efficiencies during Analysis
  /// If efficiencies can be found, they are used in the Analysis, but not produced
  /// If efficiencies cannot be found, they are not used in the Analysis, but produced
  bool makeEfficiencies();



  /// Fill histograms needed for b-tag efficiencies
  void fillBtagHistograms(const std::vector<int>& jetIndices,
                          const std::vector<double>& bTagDiscriminant,
                          const VLV& jets,
                          const std::vector<int>& jetPartonFlavours,
                          const double weight);

  /// Produce b-tag efficiencies
  void produceBtagEfficiencies(const std::string& channel);



  /// Get b-tag per-event scale factor
  double calculateBtagSF(const std::vector<int>& jetIndices,
                         const VLV& jets,
                         const std::vector<int>& jetPartonFlavours);

  /// Method takes the indices of b-tagged jets,
  /// and overwrites them with the b-tagged jets after randomised tag flipping
  /// Method explained in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFUtil
  /// and in: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
  void indexOfBtags(std::vector<int>& bjetIndices,
                    const std::vector<int>& jetIndices,
                    const VLV& jets,
                    const std::vector<int>& jetPartonFlavours,
                    const std::vector<double>& btagDiscriminants)const;

  /// Prepare b-tagging scale factors (efficiency histograms)
  void prepareBTags(TSelectorList* output, const std::string& channel);



private:

  /// Input directory for files holding histograms of b-tagging efficiencies
  std::string inputDirName_;

  /// Output directory for files holding histograms, in case b-tagging efficiencies are produced in Analysis
  std::string outputDirName_;

  /// Name of the file with btag histograms
  std::string fileName_;



  /// Whether to produce btag efficiencies or whether to read in already produced ones
  bool makeEfficiencies_;

  /// Pointer for bookkeeping of histograms
  TSelectorList* selectorList_;
  
  /// Map of the file names for each channel
  std::map<std::string, std::string> channelFileNames_;
  
  /// Map of the sample names for each channel
  std::map<std::string, std::string> channelSampleNames_;

};




#endif