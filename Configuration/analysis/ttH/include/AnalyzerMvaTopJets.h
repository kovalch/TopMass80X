#ifndef AnalyzerMvaTopJets_h
#define AnalyzerMvaTopJets_h

#include <vector>
#include <map>

class TString;
class TH1;
class TH2D;

#include "AnalyzerBaseClass.h"
#include "MvaVariablesTopJets.h"

class JetCategories;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class HiggsGenObjects;
class KinRecoObjects;
class MvaVariablesTopJets;
class MvaReaderBase;
namespace tth{
    class RecoLevelWeights;
    class GenLevelWeights;
    class RecoObjectIndices;
    class GenObjectIndices;
}









/// Class for basic histograms that are filled simultaneously for any step
class AnalyzerMvaTopJets : public AnalyzerBaseClass{

public:

    /// Constructor
    AnalyzerMvaTopJets(const char* mva2dWeightsFile,
                       const std::vector<TString>& selectionStepsNoCategories,
                       const std::vector<TString>& stepsForCategories = std::vector<TString>(),
                       const JetCategories* jetCategories =0);

    /// Destructor
    ~AnalyzerMvaTopJets(){}



private:

    /// Book all histograms for given selection step
    virtual void bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram);

    /// Fill all histograms for given selection step
    virtual void fillHistos(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                            const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                            const KinRecoObjects& kinRecoObjects,
                            const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                            const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                            const double& weight, const TString& step,
                            std::map<TString, TH1*>& m_histogram);


    /// Book 1-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    void bookHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                            const TString& name, const TString& title,
                            const int& nBinX, const double& xMin, const double& xMax);

    /// Book 2-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    void bookHistosInclExcl2D(std::map<TString, TH1*>& m_histogram, const TString& prefix, const TString& step,
                              const TString& name, const TString& title,
                              const int& nBinX, const double& xMin, const double& xMax,
                              const int& nBinY, const double& yMin, const double& yMax);

    /// Fill 1-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    void fillHistosInclExcl(std::map<TString, TH1*>& m_histogram, const TString& name,
                            const float& variable,
                            const MvaVariablesTopJets* mvaTopJetsVariables, const double& weight =1.);

    /// Fill 2-D histograms exclusively for correct, swapped and wrong combinations, and inclusively
    void fillHistosInclExcl2D(std::map<TString, TH1*>& m_histogram, const TString& name,
                              const float& variable1, const float& variable2,
                              const MvaVariablesTopJets* mvaTopJetsVariables, const double& weight =1.);



    struct MvaWeightsStruct{
    public:

        MvaWeightsStruct(const std::string& stepName, const std::vector<std::string>& v_nameCorrect,
                         const std::vector<std::string>& v_nameSwapped, const char* mva2dWeightsFile);
        ~MvaWeightsStruct(){}

        std::string stepName()const{return stepName_;}
        const std::map<std::string, MvaReaderBase*>& correctWeights()const{return m_correct_;}
        const std::map<std::string, MvaReaderBase*>& swappedWeights()const{return m_swapped_;}
        const std::map<std::string, std::map<std::string, TH2D*> >& combinedWeights()const{return m_combined_;}


    private:

        std::string stepName_;

        std::map<std::string, MvaReaderBase*> m_correct_;
        std::map<std::string, MvaReaderBase*> m_swapped_;
        std::map<std::string, std::map<std::string, TH2D*> > m_combined_;
    };


    std::vector<MvaWeightsStruct> v_mvaWeightsStruct_;

    void bookHistosPerSet(const TString& step, std::map<TString, TH1*>& m_histogram, const MvaWeightsStruct& mvaWeightsStruct);

    void bookMvaSpecificHistos(const TString& step, std::map<TString, TH1*>& m_histogram,
                               const TString& mvaConfigName, const TString& mvaType);

    void fillHistosPerSet(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                          const TopGenObjects& topGenObjects, const HiggsGenObjects& higgsGenObjects,
                          const KinRecoObjects& kinRecoObjects,
                          const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                          const tth::GenLevelWeights& genLevelWeights, const tth::RecoLevelWeights& recoLevelWeights,
                          const double& weight, const TString& step,
                          std::map<TString, TH1*>& m_histogram,
                          const MvaWeightsStruct& mvaWeightsStruct
                         );

    void fillWeightHistos(const MvaVariablesTopJetsPerEvent& mvaTopJetsVariablesPerEvent,
                          const std::vector<float>& v_mvaWeight, const size_t maxWeightIndex,
                          const double& weight, const tth::RecoObjectIndices& recoObjectIndices,
                          const tth::GenObjectIndices& genObjectIndices, std::map<TString, TH1*>& m_histogram,
                          const TString& mvaType, const std::string& mvaConfigName1, const std::string& mvaConfigName2 ="");

    void fillBestWeightHistos(const std::vector<float>& v_mvaWeights,
                              const RecoObjects& recoObjects,
                              const tth::RecoObjectIndices& recoObjectIndices, const tth::GenObjectIndices& genObjectIndices,
                              const double& weight,
                              std::map<TString, TH1*>& m_histogram,
                              const std::string& mvaType, const std::string& mvaConfigName);

};




#endif







