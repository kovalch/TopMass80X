#ifndef analysisStructs_h
#define analysisStructs_h

#include <string>
#include <vector>

#include "analysisStructsFwd.h"





namespace ttbar{
    
    struct GenLevelWeights{
        GenLevelWeights(const double& weightMadgraphCorrection, const double& weightPileup, const double& weightGenerator,
                        const double& trueLevelWeightNoPileup, const double& trueLevelWeight);
        ~GenLevelWeights(){}
        
        #ifndef __CINT__
        const double& weightMadgraphCorrection_;
        const double& weightPileup_;
        const double& weightGenerator_;
        
        const double& trueLevelWeightNoPileup_;
        const double& trueLevelWeight_;
        #endif
    };
    
    
    
    struct RecoLevelWeights{
        RecoLevelWeights(const double& weightLeptonSF, const double& weightTriggerSF, const double& weightBtagSF,
                         const double& weightNoPileup,const double& weight);
        ~RecoLevelWeights(){}
        
        #ifndef __CINT__
        const double& weightLeptonSF_;
        const double& weightTriggerSF_;
        const double& weightBtagSF_;
        
        const double& weightNoPileup_;
        const double& weight_;
        #endif
    };
    
    
    
    struct GenObjectIndices{
        GenObjectIndices(const int& genBjetFromTopIndex, const int& genAntiBjetFromTopIndex,
                         const int& recoBjetFromTopIndex, const int& recoAntiBjetFromTopIndex,
                         const int& genBjetFromHiggsIndex, const int& genAntiBjetFromHiggsIndex,
                         const int& recoBjetFromHiggsIndex, const int& recoAntiBjetFromHiggsIndex,
                         const std::vector<int>& genVisJetIndices
                        );
        ~GenObjectIndices(){}
        
        bool uniqueGenTopMatching()const;
        bool uniqueRecoTopMatching()const;
        
        bool uniqueGenMatching()const;
        
        bool uniqueGenHiggsMatching()const;
        bool uniqueRecoHiggsMatching()const;
        
        bool uniqueRecoMatching()const;
        
        bool isCorrectPairFromTop(const int& bIndex, const int& antiBIndex)const;
        bool isSwappedPairFromTop(const int& bIndex, const int& antiBIndex)const;
        bool isPairFromTop(const int& bIndex, const int& antiBIndex)const;
        
        bool isCorrectPairFromHiggs(const int& bIndex, const int& antiBIndex)const;
        bool isSwappedPairFromHiggs(const int& bIndex, const int& antiBIndex)const;
        bool isPairFromHiggs(const int& bIndex, const int& antiBIndex)const;
        
        #ifndef __CINT__
        const int& genBjetFromTopIndex_;
        const int& genAntiBjetFromTopIndex_;
        const int& recoBjetFromTopIndex_;
        const int& recoAntiBjetFromTopIndex_;
        
        
        const int& genBjetFromHiggsIndex_;
        const int& genAntiBjetFromHiggsIndex_;
        const int& recoBjetFromHiggsIndex_;
        const int& recoAntiBjetFromHiggsIndex_;
        
        const std::vector<int>& genVisJetIndices_;
        #endif
    };
    
    
    
    struct RecoObjectIndices{
        RecoObjectIndices(const std::vector<int>& allLeptonIndices,
                          const std::vector<int>& leptonIndices, const std::vector<int>& antiLeptonIndices,
                          const int& leptonIndex, const int& antiLeptonIndex,
                          const int& leadingLeptonIndex, const int& nLeadingLeptonIndex,
                          const int& leptonXIndex, const int& leptonYIndex,
                          const std::vector<int>& jetIndices,
                          const std::vector<int>& bjetIndices);
        ~RecoObjectIndices(){}
        
        #ifndef __CINT__
        const std::vector<int>& allLeptonIndices_;
        const std::vector<int>& leptonIndices_;
        const std::vector<int>& antiLeptonIndices_;
        const int& leptonIndex_;
        const int& antiLeptonIndex_;
        const int& leadingLeptonIndex_;
        const int& nLeadingLeptonIndex_;
        const int& leptonXIndex_;
        const int& leptonYIndex_;
        
        const std::vector<int>& jetIndices_;
        
        const std::vector<int>& bjetIndices_;
        
        #endif
    };
    
}










#endif








