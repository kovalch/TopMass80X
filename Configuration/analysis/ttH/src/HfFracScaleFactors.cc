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
#include <THStack.h>
#include <TFile.h>
#include <TError.h>

#include "HfFracScaleFactors.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"





HfFracScaleFactors::HfFracScaleFactors(const Samples& samples, RootFileReader* const rootFileReader):
rootFileReader_(rootFileReader)
{
    std::cout<<"--- Beginning production of Heavy-Flavour fraction scale factors\n\n";
    
    // Setting name of the histogram used for template fit
    histoTemplateName_ = "hfFracScaling_btag_multiplicity";
    
    // Setting id for each sample type as it will appear in the list of histograms for the fit
    // Data MUST go first
    sampleTypeIds_[Sample::data] = 0;
    // Combine samples by assigning the same id
    sampleTypeIds_[Sample::ttbb] = 1;
    sampleTypeIds_[Sample::ttb] = 2;
    sampleTypeIds_[Sample::tt2b] = 2;
    sampleTypeIds_[Sample::ttcc] = 3;
    sampleTypeIds_[Sample::ttother] = 3;
    // Backgrounds MUST go last
    sampleTypeIds_[Sample::dummy] = 4;
    
    // FIXME: Check that there are no gaps in the list of ids
    
    this->produceScaleFactors(samples);
    
    std::cout<<"\n=== Finishing production of Heavy-Flavour fraction scale factors\n\n";
}



void HfFracScaleFactors::produceScaleFactors(const Samples& samples)
{
    // Extract steps for Heavy-Flavour fraction scaling from first file in map
    const SystematicChannelSamples& m_systematicChannelSamples = samples.getSystematicChannelSamples();
    const TString& filename = m_systematicChannelSamples.begin()->second.begin()->second.begin()->inputFile();
    const std::vector<std::pair<TString, TString> > v_nameStepPair = tth::nameStepPairs(filename, histoTemplateName_);
    
    // Produce scale factors
    for(const auto& nameStepPair : v_nameStepPair) {
        if(nameStepPair.second.Contains("_cate")) continue;
        this->produceScaleFactors(nameStepPair.second, samples);
    }
    
    // Print table
    std::cout<<"Step   \t\tSystematic\tChannel\t\tScale factor (ttbb | ttb | ttother)\n";
    std::cout<<"-------\t\t----------\t-------\t\t-----------------------------------\n";
    for(auto hfFracScaleFactorsPerStep : m_hfFracScaleFactors_){
        const TString& step(hfFracScaleFactorsPerStep.first);
        for(auto hfFracScaleFactorsPerSystematic : hfFracScaleFactorsPerStep.second){
           const Systematic::Systematic& systematic(hfFracScaleFactorsPerSystematic.first);
           for(auto hfFracScaleFactorsPerChannel : hfFracScaleFactorsPerSystematic.second){
               const Channel::Channel& channel(hfFracScaleFactorsPerChannel.first);
//                if(channel != Channel::emu) continue;
               std::cout<<step<<"\t\t"<<systematic.name()<<"\t\t"
                        <<Channel::convert(channel)<<"      \t   "
                        <<std::fixed<<std::setprecision(3)
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttbb).val<<" |    "
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttb).val<<" |    "
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttother).val<<" | \n";
                        
                std::cout<<"\t\t\t\t\t\t+- "
                        <<std::fixed<<std::setprecision(3)
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttbb).err<<" | +- "
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttb).err<<" | +- "
                        <<hfFracScaleFactorsPerChannel.second.at(Sample::SampleType::ttother).err<<" | \n";
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
        const std::vector<Channel::Channel> v_channel = {Channel::combined};

        std::vector<TH1*> histos_comb;
        
        
        for(Channel::Channel channel : v_channel){
            std::vector<TH1*> histos;
            const std::vector<Sample>& v_sample(channelSamples.at(channel));
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                const Sample& sample = v_sample.at(iSample);
                const Sample::SampleType& sampleType = sample.sampleType();
                
                TH1* h = rootFileReader_->GetClone<TH1D>(sample.inputFile(), TString(histoTemplateName_).Append(step));
                h->Sumw2();
                
                if(sampleType != Sample::data){
                    const double& weight = globalWeights.at(systematic).at(channel).at(iSample);
                    h->Scale(weight);
                }

                int sampleId = -1;
                
                // Treating the sample as background if it doesn't have id (will be substracted from data)
                if (sampleTypeIds_.count(sampleType) == 0) sampleId = sampleTypeIds_[Sample::dummy];
                else sampleId = sampleTypeIds_.at(sampleType);
                
                // Adding histogram for the appropriate sample according to its id
                char tempName[20];
                if (sampleId < (int)histos.size()) {
                    histos.at(sampleId)->Add(h); 
                } else {
                    sprintf(tempName, "h_%d", sampleId);
                    histos.push_back((TH1*)h->Clone(tempName));
                }
                
                if (sampleId < (int)histos_comb.size()) {
                    histos_comb.at(sampleId)->Add(h); 
                } else {
                    sprintf(tempName, "h_%d_comb", sampleId);
                    histos_comb.push_back((TH1*)h->Clone(tempName));
                }
                
            }
            
            const std::vector<HfFracScaleFactors::ValErr> sampleSFs = getScaleFactorsFromHistos(histos, step, channel);
            
            for(auto sampleTypeId : sampleTypeIds_) {
                if(sampleTypeId.first == Sample::dummy) break;
                m_hfFracScaleFactors_[step][systematic][channel][sampleTypeId.first] = sampleSFs.at(sampleTypeId.second);
            }
            
            histos.clear();
            
        }
        
        const std::vector<HfFracScaleFactors::ValErr> sampleSFs = getScaleFactorsFromHistos(histos_comb, step, Channel::combined);
        
        for(auto sampleTypeId : sampleTypeIds_) {
            if(sampleTypeId.first == Sample::dummy) break;
            m_hfFracScaleFactors_[step][systematic][Channel::combined][sampleTypeId.first] = sampleSFs.at(sampleTypeId.second);
        }
        
        histos_comb.clear();
    }
}


const std::vector<HfFracScaleFactors::ValErr> HfFracScaleFactors::getScaleFactorsFromHistos(const std::vector<TH1*> histos, const TString& step, 
                                                                                           const Channel::Channel channel)const
{
    if(histos.size()<2) {
        std::cerr<<"Error in fitting function getScaleFactorsFromHistos! Too few histograms provided: "<< histos.size() <<"\n...break\n"<<std::endl;
        exit(17);
    }
    // Initialisation of identity scale factors
    const size_t nSamples = histos.size()-1;
    std::vector<HfFracScaleFactors::ValErr> result(nSamples, HfFracScaleFactors::ValErr(1., 1.));
    
    // Setting sets of histograms for the TFractionFitter
    TH1* h_data = histos.at(0);
    TObjArray* h_mc = new TObjArray(nSamples-1);
    
    // Subtracting background from data
    histos.at(0)->Add(histos.at(nSamples), -1);
    
    // Setting the integrals of the initial histograms
    std::vector<double> hInts;
    for(size_t iH = 0; iH < nSamples; ++iH) {
        TH1* histo = (TH1*)histos.at(iH)->Clone();
        double hInt = histo->Integral();
        hInts.push_back(hInt);
        histo->SetLineColor(iH);
        // Scaling to the number of entries to have poisson errors closer to reality
        double scale = poissonErrorScale(histo);
        histo->Scale(scale);
        
        if(iH==0) continue;
        h_mc->Add(histo);
    }
    
    if(hInts.at(0)<=0) {
        std::cout << "    WARNING: Background larger than data. Skipping.." << std::endl;
        return result;
    }
    
    // Setting up the TFractionFitter
    TFractionFitter* fitter = new TFractionFitter(h_data, h_mc, "Q");
    
    // Constraining the components of the fit
    for(size_t sampleId = 1; sampleId<nSamples; ++sampleId) {
       fitter->Constrain(sampleId, 0., 1.);
    }
    // Setting the bin range used in the fit
    fitter->SetRangeX(2,5);

    Int_t status = fitter->Fit();
// //     fitter->ErrorAnalysis(1.);
    if(status != 0) {
        std::cerr << "WARNING!!! Fit failed with status: " << status << " [step: " << step << "; channel: " << channel << "]" << std::endl;
        return result;
    }
    
    // Extracting fit result for each sample
    for(size_t sampleId = 1; sampleId<nSamples; ++sampleId) {
        Double_t v(0.), e(0.);
        fitter->GetResult(sampleId-1, v, e);
        result.at(sampleId).val = v / (hInts.at(sampleId) / hInts.at(0));
        result.at(sampleId).err = e / (hInts.at(sampleId) / hInts.at(0));
    }
    
    return result;
    
    // Storing the plot showing the result of the fit
    
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
    
    TFile* out_root = new TFile(histoTemplateName_+"_source.root", "RECREATE");
    TCanvas* canvas = new TCanvas("","");
    canvas->SetLogy();
    
    THStack* stack = new THStack("def", "def");
    for(size_t iHisto = 1; iHisto < nSamples; ++iHisto) {
        histos.at(iHisto)->Scale(result.at(iHisto).val);
        histos.at(iHisto)->SetLineColor(iHisto+1);
        histos.at(iHisto)->SetFillColor(iHisto+1);
        stack->Add(histos.at(iHisto));
    }
    histos.at(0)->SetMarkerStyle(20);
    histos.at(0)->SetLineWidth(2);
    histos.at(0)->Draw();
    stack->Draw("same HIST");
    histos.at(0)->Draw("same");
    canvas->Print(histoTemplateName_+step+"_stack.eps");
    canvas->Write(step+"_stack");
    
    canvas->Clear();
    canvas->SetLogy(0);
    
    for(size_t iHisto = 0; iHisto < nSamples; ++iHisto) {
        normalize(histos.at(iHisto));
        histos.at(iHisto)->SetLineWidth(2);
        histos.at(iHisto)->SetLineColor(iHisto+1);
        histos.at(iHisto)->SetMarkerStyle(iHisto+19);
        histos.at(iHisto)->SetMarkerColor(iHisto+1);
        if(iHisto == 0) histos.at(iHisto)->Draw("LP");
        else histos.at(iHisto)->Draw("sameLP");
    }
    
    canvas->Write(step+"_shapes");
    canvas->Print(histoTemplateName_+step+"_shapes.eps");
    out_root->Close();
    
    delete out_root;
    delete canvas;
    
    return result;
}



const double& HfFracScaleFactors::hfFracScaleFactor(const TString& step,
                                                    const Systematic::Systematic& systematic,
                                                    const Channel::Channel& channel,
                                                    const Sample::SampleType& sampleType)const
{
    if(m_hfFracScaleFactors_.at(step).find(systematic) == m_hfFracScaleFactors_.at(step).end()){
        std::cerr<<"Heavy-Flavour fraction scale factor requested, but not existent for Systematic: "<<systematic.name()
                 <<"\n...break\n"<<std::endl;
        exit(16);
    }
    if(m_hfFracScaleFactors_.at(step).at(systematic).find(channel) == m_hfFracScaleFactors_.at(step).at(systematic).end()) {
        std::cerr<<"Heavy-Flavour fraction scale factor requested, but not existent for Channel: "<<Channel::convert(channel)
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
    
    if(m_hfFracScaleFactors_.find(step) == m_hfFracScaleFactors_.end()) return -1;
    
    // Access the Heavy-Flavour fraction scale factor
//     weight *= this->hfFracScaleFactor(step, systematic, sample.finalState(), sample.sampleType());
    weight *= this->hfFracScaleFactor(step, systematic, Channel::combined, sample.sampleType());
    return 1;
}


void HfFracScaleFactors::normalize ( TH1* histo )const
{
    double integral = histo->Integral();
    histo->Scale ( 1.0/integral );
}

double HfFracScaleFactors::poissonErrorScale(const TH1* histo)const
{
    double scale = 0.;
    const int nBins = histo->GetNbinsX();
    // Finding the largest scale to have no bins with poisson uncertainty larger than the weighted one
    for(int iBin = 1; iBin <= nBins; ++iBin) {
        if(histo->GetBinContent(iBin) == 0) continue;
        double error_real = histo->GetBinError(iBin);
        double error_poisson = sqrt(histo->GetBinContent(iBin));
        double scale_new = error_poisson/error_real;
        if(scale_new < scale) continue;
        scale = scale_new;
    }
    
    return scale;
}