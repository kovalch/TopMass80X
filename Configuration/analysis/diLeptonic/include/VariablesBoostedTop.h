#ifndef VariablesBoostedTop_h
#define VariablesBoostedTop_h

#include <map>
#include <vector>
#include <string>

#include "VariablesBase.h"
#include "analysisStructsFwd.h"
#include "../../common/include/classesFwd.h"

class EventMetadata;
class RecoObjects;
class TopGenObjects;
class CommonGenObjects;
class KinematicReconstructionSolutions;
namespace ttbar{
    class RecoObjectIndices;
    class GenObjectIndices;
    class GenLevelWeights;
    class RecoLevelWeights;
}









/// Class holding the variables for Boosted Top analysis
class VariablesBoostedTop : public VariablesBase{
    
public:
    
    /// Empty constructor
    VariablesBoostedTop();
    
    /// Constructor setting up input variables from physics objects
    VariablesBoostedTop(const EventMetadata& eventMetadata, 
                   const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinRecoObjects,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight);
    
    /// Destructor
    ~VariablesBoostedTop(){}
    
    /// Fill the input structs for all jet combinations of one event
    static std::vector<VariablesBase*> fillVariables(const EventMetadata& eventMetadata,
                                                     const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                                          const TopGenObjects& topGenObjects,
                                                          const KinematicReconstructionSolutions& kinRecoObjects,
                                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                                          const double& weight);
    
    VariableFloat lep_pt_;
    VariableFloat anti_lep_pt_;
    
    VariableFloat top_pt_;
    VariableFloat topbar_pt_;
    VariableFloat ttbar_delta_phi_;
    VariableFloat ttbar_pt_;
    VariableFloat top_rapidity_;
    VariableFloat ttbar_delta_eta_;
    VariableFloat ttbar_rapidity_;
    VariableFloat ttbar_mass_;
    VariableInt jet_multiplicity_;
    VariableFloat x1_;
    VariableFloat x2_;
    VariableFloat mlblbmet_;
    
    VariableFloat gen_top_pt_;
    VariableFloat gen_topbar_pt_;
    VariableFloat gen_ttbar_delta_phi_;
    VariableFloat gen_ttbar_pt_;
    VariableFloat gen_top_rapidity_;
    VariableFloat gen_ttbar_delta_eta_;
    VariableFloat gen_ttbar_rapidity_;
    VariableFloat gen_ttbar_mass_;
    VariableInt   gen_jet_multiplicity_;
    VariableFloat gen_x1_;
    VariableFloat gen_x2_;
    VariableFloat gen_mlblbmet_;

    VariableInt entry_;
    VariableInt isTopGen_;
    VariableInt isKinReco_;
    VariableFloat trueLevelWeight_;
    
    
    
private:
    
    // The names associated to the variables
    
    static constexpr const char* name_lep_pt_ = "lep_pt";
    static constexpr const char* name_anti_lep_pt_ = "anti_lep_pt";
    
    static constexpr const char* name_top_pt_ = "top_pt";
    static constexpr const char* name_topbar_pt_ = "topbar_pt";
    static constexpr const char* name_ttbar_delta_phi_ = "ttbar_delta_phi";
    static constexpr const char* name_ttbar_pt_ = "ttbar_pt";
    static constexpr const char* name_top_rapidity_ = "top_rapidity";
    static constexpr const char* name_ttbar_delta_eta_ = "ttbar_delta_eta";
    static constexpr const char* name_ttbar_rapidity_ = "ttbar_rapidity";
    static constexpr const char* name_ttbar_mass_ = "ttbar_mass";
    static constexpr const char* name_jet_multiplicity_ = "jet_multiplicity";
    static constexpr const char* name_x1_ = "x1";
    static constexpr const char* name_x2_ = "x2";
    static constexpr const char* name_mlblbmet_ = "mlblbmet";
    
    static constexpr const char* name_gen_top_pt_ = "gen_top_pt";
    static constexpr const char* name_gen_topbar_pt_ = "gen_topbar_pt";
    static constexpr const char* name_gen_ttbar_delta_phi_ = "gen_ttbar_delta_phi";
    static constexpr const char* name_gen_ttbar_pt_ = "gen_ttbar_pt";
    static constexpr const char* name_gen_top_rapidity_ = "gen_top_rapidity";
    static constexpr const char* name_gen_ttbar_delta_eta_ = "gen_ttbar_delta_eta";
    static constexpr const char* name_gen_ttbar_rapidity_ = "gen_ttbar_rapidity";
    static constexpr const char* name_gen_ttbar_mass_ = "gen_ttbar_mass";
    static constexpr const char* name_gen_jet_multiplicity_ = "gen_jet_multiplicity";
    static constexpr const char* name_gen_x1_ = "gen_x1";
    static constexpr const char* name_gen_x2_ = "gen_x2";
    static constexpr const char* name_gen_mlblbmet_ = "gen_mlblbmet";
    
    static constexpr const char* name_entry_ = "entry";
    static constexpr const char* name_isTopGen_ = "isTopGen";
    static constexpr const char* name_isKinReco_ = "isKinReco";
    static constexpr const char* name_trueLevelWeight_ = "trueLevelWeight";
    
    
};






#endif







