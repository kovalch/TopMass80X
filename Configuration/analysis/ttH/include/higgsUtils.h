#ifndef tth_utils_h
#define tth_utils_h

#include <string>
#include <vector>
#include <utility>

class TString;





namespace tth{
    
    /// Return the path where relevant input data is stored
    const std::string DATA_PATH_TTH();
    
    
    
    /// Get vector of pairs with < stepName, objectName from ROOT-file (e.g. histo or tree) >
    std::vector<std::pair<TString, TString> > nameStepPairs(const char* filename,
                                                            const char* objectNameFragment,
                                                            const std::vector<TString>& selectedSteps =std::vector<TString>());
    
    
    
    /// Assign a category name for one category, negative value means no category and returns empty string
    TString categoryName(const int& category);
    
    /// Assign a category name for a set of categories
    TString categoryName(const std::vector<int>& v_category);
    
    /// Assign a step name for a given short name of the step, and potentially for a specific category
    TString stepName(const TString& stepShort, const int& category =-1);
    
    /// Assign a step name for a given short name of the step, and potentially for a set of categories
    TString stepName(const TString& stepShort, const std::vector<int>& v_category);
    
    
    
    /// Get from a TString the selection step of pattern "step*"
    TString extractSelectionStep(const TString& name);
    
    /// Get from a TString the jet category of pattern "cate*"
    TString extractJetCategory(const TString& name);
    
    /// Get from a TString the selection step of pattern "step*", combined with 0, 1 or several categories of pattern "cate*"
    TString extractSelectionStepAndJetCategory(const TString& name);
    
    /// Get from a TString the number of categories of pattern "cate*"
    int numberOfCategories(const TString& name);
    
    /// Helper functions only needed by functions defined in this file
    namespace helper{
        
        /// Helper function to get the fragment containing a searchPattern,
        /// fragments are separated by given token,
        /// Flag allowMultiple to switch whether pattern can exist several times or only once
        TString searchFragmentByToken(const TString& name, const TString& searchPattern,
                                      const TString& token, const bool allowMultiple =false);
    }
}





#endif




