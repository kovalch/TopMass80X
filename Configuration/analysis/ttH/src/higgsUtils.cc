#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include <TObjArray.h>
#include <TString.h>
#include <Rtypes.h>

#include "higgsUtils.h"
#include "../../common/include/utils.h"
#include "../../common/include/RootFileReader.h"





const std::string tth::DATA_PATH_TTH()
{
    std::string result(common::CMSSW_BASE());
    result.append("/src/TopAnalysis/Configuration/analysis/ttH/data");
    return result;
}



std::vector<std::pair<TString, TString> > tth::nameStepPairs(const char* filename,
                                                             const char* objectNameFragment,
                                                             const std::vector<TString>& selectedSteps)
{
    std::vector<std::pair<TString, TString> > result;
    std::vector<TString> v_step;
    std::stringstream ss_step;
    std::vector<TString> v_unselectedStep;
    std::stringstream ss_unselected;
    
    // Search for all steps
    RootFileReader* fileReader(RootFileReader::getInstance());
    const std::vector<TString>& v_objectName = fileReader->findObjects(filename, objectNameFragment, false);
    for(std::vector<TString>::const_iterator i_objectName = v_objectName.begin(); i_objectName != v_objectName.end(); ++i_objectName){
        const TString& step = tth::extractSelectionStepAndJetCategory(*i_objectName);
        
        // Reject steps in case only selected steps should be chosen
        if(selectedSteps.size() && std::find(selectedSteps.begin(), selectedSteps.end(), step) == selectedSteps.end()){
            if(std::find(v_unselectedStep.begin(), v_unselectedStep.end(), step) == v_unselectedStep.end()){
                ss_unselected<<step<<", ";
                v_unselectedStep.push_back(step);
            }
            continue;
        }
        
        // Add selection steps
        if(std::find(v_step.begin(), v_step.end(), step) == v_step.end()){
            ss_step<<step<<", ";
            v_step.push_back(step);
        }
        result.push_back(std::make_pair(*i_objectName, step));
    }
    
    std::cout<<"Found chosen selection steps:\n"<<ss_step.str()<<std::endl;
    if(selectedSteps.size() && ss_unselected.str()!="") std::cout<<"Found rejected selection steps:\n"<<ss_unselected.str()<<std::endl;
    return result;
}



TString tth::categoryName(const int& category)
{
    std::ostringstream result;
    if(category >= 0){
        result<<"_cate"<<category;
    }
    return result.str().c_str();
}



TString tth::categoryName(const std::vector<int>& v_category)
{
    std::ostringstream result;
    std::vector<int> v_clone = v_category;
    std::sort(v_clone.begin(), v_clone.end());
    for(std::vector<int>::iterator i_category = v_clone.begin(); i_category != v_clone.end(); ++i_category){
        const int category = *i_category;
        if(std::find(i_category+1, v_clone.end(), category) != v_clone.end()){
            std::cerr<<"Error in tth::categoryName()! The following category exists twice: "
                     <<category<<"\n...break\n"<<std::endl;
            exit(841);
        }
        if(category < 0){
            std::cerr<<"Error in tth::categoryName()! Category needs to be >= 0, but is: "
                     <<category<<"\n...break\n"<<std::endl;
            exit(841);
        }
        result<<categoryName(category);
    }
    return result.str().c_str();
}



TString tth::stepName(const TString& stepShort, const int& category)
{
    std::ostringstream result;
    result<<"_step"<<stepShort;
    result<<categoryName(category);
    return result.str().c_str();
}



TString tth::stepName(const TString& stepShort, const std::vector<int>& v_category)
{
    std::ostringstream result;
    result<<"_step"<<stepShort;
    result<<categoryName(v_category);
    return result.str().c_str();
}



TString tth::extractSelectionStep(const TString& name)
{
    TString result(name);
    // Remove e.g. things like .root ??
    result = helper::searchFragmentByToken(result, "step", ".");
    // Further separation by "_" to find string containing step
    result = helper::searchFragmentByToken(result, "step", "_");
    //std::cout<<"The extracted selection step is (step/histogram name): "<<result<<" / "<<name<<std::endl;
    return result;
}



TString tth::extractJetCategory(const TString& name)
{
    TString result(name);
    // Remove e.g. things like .root ??
    result = helper::searchFragmentByToken(result, "cate", ".", true);
    // Further separation by "_" to find string containing step
    result = helper::searchFragmentByToken(result, "cate", "_", true);
    //std::cout<<"The extracted category bin is (bin/histogram name): "<<result<<" / "<<name<<std::endl;
    return result;
}



TString tth::extractSelectionStepAndJetCategory(const TString& name)
{
    const TString step(extractSelectionStep(name));
    const TString cate(extractJetCategory(name));
    const TString result = step+cate;
    return result;
}



int tth::numberOfCategories(const TString& name)
{
    TString category = extractJetCategory(name);
    TObjArray* a_nameFragment = TString(category).Tokenize("_");
    int result = a_nameFragment->GetEntriesFast();
    a_nameFragment->Delete();
    return result;
}



TString tth::helper::searchFragmentByToken(const TString& name, const TString& searchPattern,
                                           const TString& token, const bool allowMultiple)
{
    TString result;
    TObjArray* a_nameFragment = TString(name).Tokenize(token);
    bool alreadyFound(false);
    for(Int_t iFragment= 0; iFragment < a_nameFragment->GetEntriesFast(); ++iFragment){
        const TString& fragment = a_nameFragment->At(iFragment)->GetName();
        if(fragment.Contains(searchPattern)){
            if(!allowMultiple && alreadyFound){
                std::cerr<<"ERROR in searchFragmentByToken()! Ambiguous string, name contains search pattern \""
                         <<searchPattern<<"\" in more than one fragment: "<<name<<"\n...break\n"<<std::endl;
                exit(33);
            }
            alreadyFound = true;
            result.Append("_").Append(fragment);
        }
    }
    a_nameFragment->Delete();
    return result;
}







