#include <map>
#include <utility>
#include <iostream>
#include <iomanip>

#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TObjArray.h>
#include <TFractionFitter.h>
#include <TCanvas.h>

#include "HfFracScaleFactors.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"





HfFracScaleFactors::HfFracScaleFactors(const Samples& samples, RootFileReader* const rootFileReader):
rootFileReader_(rootFileReader)
{
    std::cout<<"--- Beginning production of Heavy-Flavour fraction scale factors\n\n";
    
    this->produceScaleFactors(samples);
    
    std::cout<<"\n=== Finishing production of Heavy-Flavour fraction scale factors\n\n";
}



void HfFracScaleFactors::produceScaleFactors(const Samples& samples)
{
    // Extract steps for Heavy-Flavour fraction scaling from first file in map
    const SystematicChannelSamples& m_systematicChannelSamples = samples.getSystematicChannelSamples();
    const TString& filename = m_systematicChannelSamples.begin()->second.begin()->second.begin()->inputFile();
    const std::vector<std::pair<TString, TString> > v_nameStepPair = tth::nameStepPairs(filename, "hfFracScaling_bTag_mult_step");
    
    // Produce scale factors
    for(const auto& nameStepPair : v_nameStepPair) this->produceScaleFactors(nameStepPair.second, samples);
    
    // Print table
    std::cout<<"Step   \t\tSystematic\tChannel\t\tScale factor (ttbb | ttb | ttother)\n";
    std::cout<<"-------\t\t----------\t-------\t\t-----------------------------------\n";
    for(auto hfFracScaleFactorsPerStep : m_hfFracScaleFactors_){
        const TString& step(hfFracScaleFactorsPerStep.first);
        for(auto hfFracScaleFactorsPerSystematic : hfFracScaleFactorsPerStep.second){
           const Systematic::Systematic& systematic(hfFracScaleFactorsPerSystematic.first);
           for(auto hfFracScaleFactorsPerChannel : hfFracScaleFactorsPerSystematic.second){
               const Channel::Channel& channel(hfFracScaleFactorsPerChannel.first);
               if(channel != Channel::emu) continue;
               std::cout<<step<<"\t\t"<<Systematic::convertSystematic(systematic)<<"\t\t"
                        <<Channel::convertChannel(channel)<<"\t\t   "
                        <<std::fixed<<std::setprecision(3)
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttbb).val<<" |    "
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttb).val<<" |    "
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttother).val<<" | \n";
                        
//                 std::cout<<"\t\t\t\t\t\t+- "
//                         <<std::fixed<<std::setprecision(3)
//                         <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttbb).err<<" | +- "
//                         <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttb).err<<" | +- "
//                         <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttother).err<<" | \n";
           }
        }
    }
}



void HfFracScaleFactors::produceScaleFactors(const TString& step, const Samples& samples)
{
    TH1::AddDirectory(false);
    // Get the scale factors from the samples
    const SystematicChannelFactors globalWeights = samples.globalWeights(step).first;
    
    for(auto systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        const auto& channelSamples(systematicChannelSamples.second);
        
//         const std::vector<Channel::Channel> v_channel = {Channel::ee, Channel::emu, Channel::mumu, Channel::combined};
        const std::vector<Channel::Channel> v_channel = {Channel::emu};
        
        TH1 *h_data_comb=0, *h_bkg_comb=0, *h_ttbb_comb=0, *h_ttb_comb=0, *h_tto_comb=0;
        
        for(Channel::Channel channel : v_channel){
            TH1 *h_data=0, *h_bkg=0, *h_ttbb=0, *h_ttb=0, *h_tto=0;
            const std::vector<Sample>& v_sample(channelSamples.at(channel));
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                const Sample& sample = v_sample.at(iSample);
                const Sample::SampleType& sampleType = sample.sampleType();
                
                TH1* h = rootFileReader_->GetClone<TH1D>(sample.inputFile(), TString("hfFracScaling_bTag_mult").Append(step));
                
                if(sampleType != Sample::data){
                    const double& weight = globalWeights.at(systematic).at(channel).at(iSample);
                    h->Scale(weight);
                }
                
                switch(sampleType) {
                    case Sample::data : 
                        if(h_data) h_data->Add(h); else h_data = h;
                        if(h_data_comb) h_data_comb->Add(h); else h_data_comb = (TH1*)h->Clone("h_data_comb");
                        break;
                    case Sample::ttbb : 
                        if(h_ttbb) h_ttbb->Add(h); else h_ttbb = h;
                        if(h_ttbb_comb) h_ttbb_comb->Add(h); else h_ttbb_comb = (TH1*)h->Clone("h_ttbb_comb");
                        break;
                    case Sample::ttb : 
                        if(h_ttb) h_ttb->Add(h); else h_ttb = h;
                        if(h_ttb_comb) h_ttb_comb->Add(h); else h_ttb_comb = (TH1*)h->Clone("h_ttb_comb");
                        break;
                    case Sample::ttother : 
                        if(h_tto) h_tto->Add(h); else h_tto = h;
                        if(h_tto_comb) h_tto_comb->Add(h); else h_tto_comb = (TH1*)h->Clone("h_tto_comb");
                        break;
                    
                    default:
                        if(h_bkg) h_bkg->Add(h); else h_bkg = h;
                        if(h_bkg_comb) h_bkg_comb->Add(h); else h_bkg_comb = (TH1*)h->Clone("h_bkg_comb");
                        break;
                }
                
            }
            
            SampleTypeValueMap sampleSFs = getScaleFactorsFromHistos(h_data, h_ttbb, h_ttb, h_tto, h_bkg, step, channel);
            
            m_hfFracScaleFactors_[step][systematic][channel][Sample::SampleType::ttbb] = sampleSFs.at(Sample::SampleType::ttbb);
            m_hfFracScaleFactors_[step][systematic][channel][Sample::SampleType::ttb] = sampleSFs.at(Sample::SampleType::ttb);
            m_hfFracScaleFactors_[step][systematic][channel][Sample::SampleType::ttother] = sampleSFs.at(Sample::SampleType::ttother);
            
            delete h_data;
            delete h_ttbb;
            delete h_ttb;
            delete h_tto;
            delete h_bkg;
            
        }
        
        SampleTypeValueMap sampleSFs = getScaleFactorsFromHistos(h_data_comb, h_ttbb_comb, h_ttb_comb, h_tto_comb, h_bkg_comb, step, Channel::combined);
        
        m_hfFracScaleFactors_[step][systematic][Channel::combined][Sample::SampleType::ttbb] = sampleSFs.at(Sample::SampleType::ttbb);
        m_hfFracScaleFactors_[step][systematic][Channel::combined][Sample::SampleType::ttb] = sampleSFs.at(Sample::SampleType::ttb);
        m_hfFracScaleFactors_[step][systematic][Channel::combined][Sample::SampleType::ttother] = sampleSFs.at(Sample::SampleType::ttother);
        
        delete h_data_comb;
        delete h_ttbb_comb;
        delete h_ttb_comb;
        delete h_tto_comb;
        delete h_bkg_comb;
    }
}


const HfFracScaleFactors::SampleTypeValueMap HfFracScaleFactors::getScaleFactorsFromHistos(TH1* h_data, TH1* h_ttbb, TH1* h_ttb, 
                                                                                           TH1* h_tto, TH1* h_bkg, 
                                                                                           const TString& step, 
                                                                                           const Channel::Channel channel)const
{
    std::vector<Sample::SampleType> sampleTypes;
    sampleTypes.push_back(Sample::SampleType::ttbb);
    sampleTypes.push_back(Sample::SampleType::ttb);
    sampleTypes.push_back(Sample::SampleType::ttother);
    // Initialisation of identity scale factors
    SampleTypeValueMap map;
    for(auto sample : sampleTypes) map[sample] = {1., 1.};
    
    // Subtracting background from data
    h_data->Add(h_bkg, -1);
    // Putting all histograms to the vector
    std::vector<TH1*> v_hist;
    v_hist.push_back(h_data);
    v_hist.push_back(h_ttbb);
    v_hist.push_back(h_ttb);
    v_hist.push_back(h_tto);
    
    // Setting the integrals of the initial histograms
    double hInt_data = h_data->Integral();
    std::vector<double> hInts;
    
    if(hInt_data<=0) {
        std::cout << "    WARNING: Background larger than data. Skipping.." << std::endl;
        return map;
    }
    
    TObjArray* mc = new TObjArray(3);
    for(size_t hId=0; hId<v_hist.size(); ++hId) {
        TH1* histo = v_hist.at(hId);
        histo->SetLineColor(hId+1);
//         normalize(histo);
//         histo->Scale(histo->GetEntries());
        if(hId==0) continue;
        hInts.push_back(histo->Integral());
        mc->Add(histo);
    }
    
    // Fitting the histograms
    TFractionFitter* fitter = new TFractionFitter(h_data, mc, "Q");
    
    int sampleId = 2;
    double sampleFrac = hInts.at(sampleId)/hInt_data;
    fitter->Constrain(sampleId+1, 0.8*sampleFrac, 1.2*sampleFrac);     // tt+other
    sampleId = 1;
    sampleFrac = hInts.at(sampleId)/hInt_data;
    fitter->Constrain(sampleId+1, 0.8*sampleFrac, 6*sampleFrac);     // tt+b
    sampleId = 0;
    sampleFrac = hInts.at(sampleId)/hInt_data;
    fitter->Constrain(sampleId+1, 0.8*sampleFrac, 6*sampleFrac);     // tt+bb
    
//     fitter->SetRangeX(1,7);
    Int_t status = fitter->Fit();
    fitter->ErrorAnalysis(1.);
    if(status != 0) {
        std::cerr << "ERROR!!! Fit failed with status: " << status << " [step: " << step << "; channel: " << channel << "]" << std::endl;
        return map;
//         exit(12);
    }
    TH1* h_fit = (TH1*)fitter->GetPlot();
    h_fit->SetLineColor(6);
    h_fit->SetLineStyle(7);
    // Extracting fit result for each sample
    for(size_t sampleId = 0; sampleId<3; ++sampleId) {
        Double_t v(0.), e(0.);
        fitter->GetResult(sampleId, v, e);
        
        map[sampleTypes.at(sampleId)].val = v / (hInts.at(sampleId) / hInt_data);
        map[sampleTypes.at(sampleId)].err = e / (hInts.at(sampleId) / hInt_data);
    }
    
    return map;
    
    // Plotting the fit result (for testing)
    TCanvas* c = new TCanvas("c", "", 800, 600);
    
    TH1* h_sum = (TH1*)h_ttbb->Clone("h_sum");
    h_sum->SetLineColor(5);
    h_sum->Scale(map[Sample::SampleType::ttbb].val);
    h_sum->Add(h_ttb, map[Sample::SampleType::ttb].val);
    h_sum->Add(h_tto, map[Sample::SampleType::ttother].val);
    
    
    h_data->Draw();
    for(size_t hId=1; hId<v_hist.size(); ++hId) {
        v_hist.at(hId)->Draw("same");
    }
    h_sum->Draw("same");
    h_fit->Draw("same");
    
    char fileName[150];
    sprintf(fileName, "plots/%s_%d.eps", h_data->GetName(), channel);
    c->SaveAs(fileName);
    delete c;
    
    return map;
}


const double& HfFracScaleFactors::hfFracScaleFactor(const TString& step,
                                                    const Systematic::Systematic& systematic,
                                                    const Channel::Channel& channel,
                                                    const Sample::SampleType& sampleType)const
{
    if(m_hfFracScaleFactors_.at(step).find(systematic) == m_hfFracScaleFactors_.at(step).end()){
        std::cerr<<"Heavy-Flavour fraction scale factor requested, but not existent for Systematic: "<<Systematic::convertSystematic(systematic)
                 <<"\n...break\n"<<std::endl;
        exit(16);
    }
    if(m_hfFracScaleFactors_.at(step).at(systematic).find(channel) == m_hfFracScaleFactors_.at(step).at(systematic).end()) {
        std::cerr<<"Heavy-Flavour fraction scale factor requested, but not existent for Channel: "<<Channel::convertChannel(channel)
                 <<"\n...break\n"<<std::endl;
        exit(16);
    }
    
    return m_hfFracScaleFactors_.at(step).at(systematic).at(channel).at(sampleType).val;
}



int HfFracScaleFactors::applyScaleFactor(double& weight,
                                         const TString& fullStepname,
                                         const Sample& sample,
                                         const Systematic::Systematic& systematic)const
{
    // Check if step without category can be extracted from full step name
    const TString step = tth::extractSelectionStep(fullStepname);
    if(step == ""){
        std::cerr<<"Heavy-Flavour fraction scale factor requested, but step could not be extracted from full step name: "<<fullStepname
                 <<"\n...break\n"<<std::endl;
        exit(14);
    }
    
    // Check whether the sample is a Heavy-Flavour fraction sample and whether it should be scaled
    const bool isTt = sample.sampleType()==Sample::ttbb || sample.sampleType()==Sample::ttb || sample.sampleType()==Sample::ttother;
    if(!isTt) return 0;
    
    // Check whether Heavy-Flavour fraction scale factor exists for extracted step
    if(m_hfFracScaleFactors_.find(step) == m_hfFracScaleFactors_.end()) return -1;
    
    // Access the Heavy-Flavour fraction scale factor
//     weight *= this->hfFracScaleFactor(step, systematic, sample.finalState(), sample.sampleType());
    weight *= this->hfFracScaleFactor(step, systematic, Channel::emu, sample.sampleType());
    return 1;
}

void HfFracScaleFactors::normalize ( TH1* histo )const
{
    double integral = histo->Integral();
    histo->Scale ( 1.0/integral );
}