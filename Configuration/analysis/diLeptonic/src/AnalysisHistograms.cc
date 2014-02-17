#define AnalysisHistograms_cxx

#include <iostream>
#include <algorithm>
#include <cstdlib>

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TSelectorList.h>

#include "AnalysisHistograms.h"
#include "analysisStructs.h"
#include "ttbarUtils.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/analysisUtils.h"
#include "../../common/include/classes.h"





// --------------------------- Methods for AnalysisHistogramBase ---------------------------------------------



AnalysisHistogramsBase::AnalysisHistogramsBase(const TString& prefix,
                                               const std::vector<TString>& selectionStepsNoCategories):
prefix_(prefix),
selectorList_(0),
selectionSteps_(selectionStepsNoCategories)
{}



void AnalysisHistogramsBase::book(TSelectorList* output)
{
    // Set pointer to output, so that histograms are owned by it
    selectorList_ = output;
    
    // Book histograms for steps not separated in JetCategories
    for(const auto& stepShort : selectionSteps_){
        const TString step = ttbar::stepName(stepShort);
        this->addStep(step);
    }
    
    // Cleanup
    selectorList_ = 0;
}



void AnalysisHistogramsBase::addStep(const TString& step)
{
    // Check whether step already exists
    if(this->checkExistence(step)){
        std::cout<<"Warning in addStep()! Selection step already contained: "<<step
                 <<"\n...skip this one\n";
        return;
    }

    // Book the histograms of the specific analyser
    std::map<TString, TH1*>& m_histogram = m_stepHistograms_[step].m_histogram_;
    this->bookHistos(step, m_histogram);
}



bool AnalysisHistogramsBase::checkExistence(const TString& step)const
{
    return m_stepHistograms_.find(step) != m_stepHistograms_.end();
}



void AnalysisHistogramsBase::bookHistos(const TString&, std::map<TString, TH1*>&)
{
    // WARNING: this is empty template method, overwrite for inherited histogram class
    
    std::cerr<<"ERROR! Dummy method bookHistos() in AnalysisHistogramBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(567);
}



void AnalysisHistogramsBase::fill(const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                  const TopGenObjects& topGenObjects,
                                  const KinRecoObjects& kinRecoObjects,
                                  const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                  const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                  const double& weight, const TString& stepShort)
{
    // Set up step name and check if step exists
    const TString step = ttbar::stepName(stepShort);
    const bool stepExists(this->checkExistence(step));
    if(!stepExists) return;
    
    // Fill the histograms of the specific analyser
    std::map<TString, TH1*>& m_histogram = m_stepHistograms_[step].m_histogram_;
    this->fillHistos(recoObjects, commonGenObjects,
                     topGenObjects,
                     kinRecoObjects,
                     recoObjectIndices, genObjectIndices,
                     genLevelWeights, recoLevelWeights,
                     weight, step,
                     m_histogram);
}



void AnalysisHistogramsBase::fillHistos(const RecoObjects&, const CommonGenObjects&,
                                        const TopGenObjects&,
                                        const KinRecoObjects&,
                                        const ttbar::RecoObjectIndices&, const ttbar::GenObjectIndices&,
                                        const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
                                        const double&, const TString&,
                                        std::map<TString, TH1*>&)
{
    // WARNING: this is empty template method, overwrite for inherited histogram class
    
    std::cerr<<"ERROR! Dummy method fillHistos() in AnalysisHistogramBase is called, but overridden one should be used\n"
             <<"...break\n"<<std::endl;
    exit(568);
}



void AnalysisHistogramsBase::clear()
{
    for(auto stepHistograms : m_stepHistograms_){
        stepHistograms.second.m_histogram_.clear();
    }
    m_stepHistograms_.clear();
    selectorList_ = 0;
}







// --------------------------- Methods for EventYieldHistograms ---------------------------------------------



EventYieldHistograms::EventYieldHistograms(const std::vector<TString>& selectionStepsNoCategories):
AnalysisHistogramsBase("events_", selectionStepsNoCategories)
{
    std::cout<<"--- Beginning setting up event yield histograms\n";
    std::cout<<"=== Finishing setting up event yield histograms\n\n";
}



void EventYieldHistograms::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    TString name;
    name = "weighted";
    m_histogram[name] = this->store(new TH1D(prefix_+name+step,"event yield;;# events",8,0,8));
    m_histogram[name]->Sumw2();
}



void EventYieldHistograms::fillHistos(const RecoObjects&, const CommonGenObjects&,
                                      const TopGenObjects&,
                                      const KinRecoObjects&,
                                      const ttbar::RecoObjectIndices&, const ttbar::GenObjectIndices&,
                                      const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
                                      const double& weight, const TString&,
                                      std::map<TString, TH1*>& m_histogram)
{
    TString name;
    
    name = "weighted";
    m_histogram[name]->Fill(1., weight);
}









// --------------------------- Methods for DyScalingHistograms ---------------------------------------------



DyScalingHistograms::DyScalingHistograms(const std::vector<TString>& selectionSteps, const TString& looseStep):
AnalysisHistogramsBase("dyScaling_", selectionSteps),
looseStep_(looseStep)
{
    std::cout<<"--- Beginning setting up Drell-Yan scaling histograms\n";
    if(std::find(selectionSteps_.begin(), selectionSteps_.end(), looseStep_) == selectionSteps_.end()){
        std::cerr<<"ERROR in constructor of DyScalingHistograms!"
                 <<"Could not find in selection steps the step specified for loose histogram: "<<looseStep_
                 <<"\n...break\n"<<std::endl;
        exit(987);
    }
    std::cout<<"=== Finishing setting up Drell-Yan scaling histograms\n\n";
}



void DyScalingHistograms::bookHistos(const TString& step, std::map<TString, TH1*>& m_histogram)
{
    const TString looseStep = ttbar::stepName(looseStep_) + "zWindow";
    if(step == looseStep){
        m_histogram["h_loose"] = this->bookHisto(m_histogram["h_loose"], prefix_+"Looseh1");
    }
    
    if(step.Contains("zWindow")){
        m_histogram["h_zWindow"] = this->bookHisto(m_histogram["h_zWindow"], prefix_+"Zh1"+step);
    }
    else{
        m_histogram["h_all"] = this->bookHisto(m_histogram["h_all"], prefix_+"Allh1"+step);
        m_histogram["h_zVeto"] = this->bookHisto(m_histogram["h_zVeto"], prefix_+"TTh1"+step);
        this->bookHistos(step + "zWindow", m_stepHistograms_[step + "zWindow"].m_histogram_);
    }
}



TH1* DyScalingHistograms::bookHisto(TH1* histo, const TString& name)
{
    histo = this->store(new TH1D(name,"Dilepton Mass; m_{ll}; events",40,0,400));
    return histo;
}



void DyScalingHistograms::fillHistos(const RecoObjects& recoObjects, const CommonGenObjects&,
                                     const TopGenObjects&,
                                     const KinRecoObjects&,
                                     const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices&,
                                     const ttbar::GenLevelWeights&, const ttbar::RecoLevelWeights&,
                                     const double& weight, const TString& step,
                                     std::map<TString, TH1*>& m_histogram)
{
    const int leadingLeptonIndex = recoObjectIndices.leadingLeptonIndex_;
    const int nLeadingLeptonIndex = recoObjectIndices.nLeadingLeptonIndex_;
    
    // FIXME: should use one common function in TopAnalysis/HiggsAnalysis and here
    //bool hasLeptonPair(false);
    if(leadingLeptonIndex!=-1 && nLeadingLeptonIndex!=-1){
        ;//hasLeptonPair = true;
    }
    else return;
    
    const LV dilepton = recoObjects.allLeptons_->at(leadingLeptonIndex) + recoObjects.allLeptons_->at(nLeadingLeptonIndex);
    const double dileptonMass = dilepton.M();
    const bool isZregion = dileptonMass>76. && dileptonMass<106.;
    
    
    // Fill histograms
    const TString looseStep = ttbar::stepName(looseStep_) + "zWindow";
    if(step == looseStep) m_stepHistograms_[looseStep].m_histogram_["h_loose"]->Fill(dileptonMass, weight);
    
    if(step.Contains("zWindow")){
        m_histogram["h_zWindow"]->Fill(dileptonMass, weight);
        TString stepNoZ(step);
        stepNoZ.ReplaceAll("zWindow", "");
        m_stepHistograms_[stepNoZ].m_histogram_["h_all"]->Fill(dileptonMass, weight);
    }
    else if(!isZregion){
        m_histogram["h_zVeto"]->Fill(dileptonMass, weight);
        m_histogram["h_all"]->Fill(dileptonMass, weight);
    }
}








