#include "analysisStructs.h"





tth::GenLevelWeights::GenLevelWeights(const double& weightMadgraphCorrection, const double& weightPileup,
                                      const double& weightGenerator, const double& weightTopPt,
                                      const double& trueLevelWeightNoPileup, const double& trueLevelWeight):
weightMadgraphCorrection_(weightMadgraphCorrection),
weightPileup_(weightPileup),
weightGenerator_(weightGenerator),
weightTopPt_(weightTopPt),
trueLevelWeightNoPileup_(trueLevelWeightNoPileup),
trueLevelWeight_(trueLevelWeight)
{}



tth::RecoLevelWeights::RecoLevelWeights(const double& weightLeptonSF, const double& weightTriggerSF,
                                        const double& weightBtagSF, const double& weightKinReco,
                                        const double& weightNoPileup, const double& weight):
weightLeptonSF_(weightLeptonSF),
weightTriggerSF_(weightTriggerSF),
weightBtagSF_(weightBtagSF),
weightKinReco_(weightKinReco),
weightNoPileup_(weightNoPileup),
weight_(weight)
{}



tth::GenObjectIndices::GenObjectIndices(const std::vector<int>& genJetIndices,
                                        const std::vector<std::vector<int> >& genJetBhadronIndices,
                                        const std::vector<int>& allGenBjetIndices,
                                        const std::vector<int>& genBjetIndices,
                                        const std::vector<std::vector<int> >& genJetChadronIndices,
                                        const std::vector<int>& allGenCjetIndices,
                                        const std::vector<int>& genCjetIndices,
                                        const int& genBjetFromTopIndex, const int& genAntiBjetFromTopIndex,
                                        const int& genBjetFromHiggsIndex, const int& genAntiBjetFromHiggsIndex):
genJetIndices_(genJetIndices),
genJetBhadronIndices_(genJetBhadronIndices),
allGenBjetIndices_(allGenBjetIndices),
genBjetIndices_(genBjetIndices),
genJetMatchedRecoBjetIndices_(std::vector<int>()),
genJetChadronIndices_(genJetChadronIndices),
allGenCjetIndices_(allGenCjetIndices),
genCjetIndices_(genCjetIndices),
genJetMatchedRecoCjetIndices_(std::vector<int>()),
genBjetFromTopIndex_(genBjetFromTopIndex),
genAntiBjetFromTopIndex_(genAntiBjetFromTopIndex),
recoBjetFromTopIndex_(-1),
recoAntiBjetFromTopIndex_(-1),
genBjetFromHiggsIndex_(genBjetFromHiggsIndex),
genAntiBjetFromHiggsIndex_(genAntiBjetFromHiggsIndex),
recoBjetFromHiggsIndex_(-1),
recoAntiBjetFromHiggsIndex_(-1)
{}



tth::GenObjectIndices::GenObjectIndices(const GenObjectIndices& genObjectNoRecoMatchIndices,
                                        const std::vector<int>& genJetMatchedRecoBjetIndices,
                                        const std::vector<int>& genJetMatchedRecoCjetIndices,
                                        const int& recoBjetFromTopIndex, const int& recoAntiBjetFromTopIndex,
                                        const int& recoBjetFromHiggsIndex, const int& recoAntiBjetFromHiggsIndex):
genJetIndices_(genObjectNoRecoMatchIndices.genJetIndices_),
genJetBhadronIndices_(genObjectNoRecoMatchIndices.genJetBhadronIndices_),
allGenBjetIndices_(genObjectNoRecoMatchIndices.allGenBjetIndices_),
genBjetIndices_(genObjectNoRecoMatchIndices.genBjetIndices_),
genJetMatchedRecoBjetIndices_(genJetMatchedRecoBjetIndices),
genJetChadronIndices_(genObjectNoRecoMatchIndices.genJetChadronIndices_),
allGenCjetIndices_(genObjectNoRecoMatchIndices.allGenCjetIndices_),
genCjetIndices_(genObjectNoRecoMatchIndices.genCjetIndices_),
genJetMatchedRecoCjetIndices_(genJetMatchedRecoCjetIndices),
genBjetFromTopIndex_(genObjectNoRecoMatchIndices.genBjetFromTopIndex_),
genAntiBjetFromTopIndex_(genObjectNoRecoMatchIndices.genAntiBjetFromTopIndex_),
recoBjetFromTopIndex_(recoBjetFromTopIndex),
recoAntiBjetFromTopIndex_(recoAntiBjetFromTopIndex),
genBjetFromHiggsIndex_(genObjectNoRecoMatchIndices.genBjetFromHiggsIndex_),
genAntiBjetFromHiggsIndex_(genObjectNoRecoMatchIndices.genAntiBjetFromHiggsIndex_),
recoBjetFromHiggsIndex_(recoBjetFromHiggsIndex),
recoAntiBjetFromHiggsIndex_(recoAntiBjetFromHiggsIndex)
{}



bool tth::GenObjectIndices::uniqueGenTopMatching()const
{
    return genBjetFromTopIndex_>=0 && genAntiBjetFromTopIndex_>=0 && genBjetFromTopIndex_!=genAntiBjetFromTopIndex_;
}



bool tth::GenObjectIndices::uniqueRecoTopMatching()const
{
    if(!this->uniqueGenTopMatching()) return false;
    
    return recoBjetFromTopIndex_>=0 && recoAntiBjetFromTopIndex_>=0 && recoBjetFromTopIndex_!=recoAntiBjetFromTopIndex_;
}



bool tth::GenObjectIndices::uniqueGenHiggsMatching()const
{
    return genBjetFromHiggsIndex_>=0 && genAntiBjetFromHiggsIndex_>=0 && genBjetFromHiggsIndex_!=genAntiBjetFromHiggsIndex_;
}



bool tth::GenObjectIndices::uniqueRecoHiggsMatching()const
{
    if(!this->uniqueGenHiggsMatching()) return false;
    
    return recoBjetFromHiggsIndex_>=0 && recoAntiBjetFromHiggsIndex_>=0 && recoBjetFromHiggsIndex_!=recoAntiBjetFromHiggsIndex_;
}



bool tth::GenObjectIndices::uniqueGenMatching()const
{
    if(!this->uniqueGenTopMatching() || !this->uniqueGenHiggsMatching()) return false;
    
    return !(genBjetFromTopIndex_==genBjetFromHiggsIndex_ || genBjetFromTopIndex_==genAntiBjetFromHiggsIndex_ ||
             genAntiBjetFromTopIndex_==genBjetFromHiggsIndex_ || genAntiBjetFromTopIndex_==genAntiBjetFromHiggsIndex_);
}



bool tth::GenObjectIndices::uniqueRecoMatching()const
{
    if(!this->uniqueRecoTopMatching() || !this->uniqueRecoHiggsMatching() || !this->uniqueGenMatching()) return false;
    
    return !(recoBjetFromTopIndex_==recoBjetFromHiggsIndex_ || recoBjetFromTopIndex_==recoAntiBjetFromHiggsIndex_ ||
             recoAntiBjetFromTopIndex_==recoBjetFromHiggsIndex_ || recoAntiBjetFromTopIndex_==recoAntiBjetFromHiggsIndex_);
}



bool tth::GenObjectIndices::isCorrectPairFromTop(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoBjetFromTopIndex_ && antiBIndex==recoAntiBjetFromTopIndex_;
}



bool tth::GenObjectIndices::isSwappedPairFromTop(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoAntiBjetFromTopIndex_ && antiBIndex==recoBjetFromTopIndex_;
}



bool tth::GenObjectIndices::isPairFromTop(const int& bIndex, const int& antiBIndex)const
{
    return this->isCorrectPairFromTop(bIndex, antiBIndex) || this->isSwappedPairFromTop(bIndex, antiBIndex);
}



bool tth::GenObjectIndices::isCorrectPairFromHiggs(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoBjetFromHiggsIndex_ && antiBIndex==recoAntiBjetFromHiggsIndex_;
}



bool tth::GenObjectIndices::isSwappedPairFromHiggs(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoAntiBjetFromHiggsIndex_ && antiBIndex==recoBjetFromHiggsIndex_;
}



bool tth::GenObjectIndices::isPairFromHiggs(const int& bIndex, const int& antiBIndex)const
{
    return this->isCorrectPairFromHiggs(bIndex, antiBIndex) || this->isSwappedPairFromHiggs(bIndex, antiBIndex);
}



tth::RecoObjectIndices::RecoObjectIndices(const std::vector<int>& allLeptonIndices,
                                  const std::vector<int>& leptonIndices, const std::vector<int>& antiLeptonIndices,
                                  const int& leptonIndex, const int& antiLeptonIndex,
                                  const int& leadingLeptonIndex, const int& nLeadingLeptonIndex,
                                  const int& leptonXIndex, const int& leptonYIndex,
                                  const std::vector<int>& jetIndices, const IndexPairs& jetIndexPairs,
                                  const std::vector<int>& bjetIndices):
allLeptonIndices_(allLeptonIndices),
leptonIndices_(leptonIndices),
antiLeptonIndices_(antiLeptonIndices),
leptonIndex_(leptonIndex),
antiLeptonIndex_(antiLeptonIndex),
leadingLeptonIndex_(leadingLeptonIndex),
nLeadingLeptonIndex_(nLeadingLeptonIndex),
leptonXIndex_(leptonXIndex),
leptonYIndex_(leptonYIndex),
jetIndices_(jetIndices),
jetIndexPairs_(jetIndexPairs),
bjetIndices_(bjetIndices)
{}





