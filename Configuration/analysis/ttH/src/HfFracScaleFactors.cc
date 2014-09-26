#include <map>
#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>

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
rootFileReader_(rootFileReader),
workingDirectory_("HfFracScaleFactors")
{
    std::cout<<"--- Beginning production of Heavy-Flavour fraction scale factors\n\n";
    
    // Setting name of the histogram used for template fit
    histoTemplateName_ = "hfFracScaling_btag_multiplicity";
//     histoTemplateName_ = "basic_bjet_multiplicity";
    
    // Setting id for each sample type as it will appear in the list of histograms for the fit
    // Data MUST go first
    sampleTypeIds_[Sample::data] = 0;
    // Combine samples by assigning the same id
    sampleTypeIds_[Sample::ttbb] = 1;
    sampleTypeIds_[Sample::ttb] = 1;
    sampleTypeIds_[Sample::tt2b] = 2;
    sampleTypeIds_[Sample::ttcc] = 3;
    sampleTypeIds_[Sample::ttother] = 3;
    // Backgrounds MUST go last
    sampleTypeIds_[Sample::dummy] = 4;
    
    // Setting names to each template
    templateNames_.push_back("data_obs");
    templateNames_.push_back("ttbb");
    templateNames_.push_back("tt2b");
    templateNames_.push_back("ttOther");
    templateNames_.push_back("bkg");

    // Setting variation limit of for each template
    templateScaleLimits_.push_back(1.);	  // Doesn't mean anything
    templateScaleLimits_.push_back(3.);
    templateScaleLimits_.push_back(3.);
    templateScaleLimits_.push_back(3.);
    templateScaleLimits_.push_back(1.01);
    
    // FIXME: Check that there are no gaps in the list of ids
    
    this->produceScaleFactors(samples);
    
    // Hiding all standard ROOT warnings
    gErrorIgnoreLevel = kWarning;
    
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
    std::cout<<"Step   \t\tSystematic\tChannel\t\tScale factor | ";
    for(size_t i = 1; i<templateNames_.size(); ++i) std::cout << templateNames_.at(i) << " | ";
    std::cout << std::endl;
    std::cout<<"-------\t\t----------\t-------\t\t-------------------------------------------------------\n";
    for(auto hfFracScaleFactorsPerStep : m_hfFracScaleFactors_){
        const TString& step(hfFracScaleFactorsPerStep.first);
        for(auto hfFracScaleFactorsPerSystematic : hfFracScaleFactorsPerStep.second){
            const Systematic::Systematic& systematic(hfFracScaleFactorsPerSystematic.first);
            for(auto hfFracScaleFactorsPerChannel : hfFracScaleFactorsPerSystematic.second){
               
                const Channel::Channel& channel(hfFracScaleFactorsPerChannel.first);
                std::cout<<step<<"\t\t"<<systematic.name()<<"\t\t"
                        <<Channel::convert(channel)<<"      \t   "
                        <<std::fixed<<std::setprecision(3);
                for(int i = 1; i<(int)templateNames_.size(); ++i) {
                    std::cout << hfFracScaleFactorsPerChannel.second.at(sampleTypeForId(i)).val<<" |    ";
                }
                std::cout << std::endl;
                        
                std::cout<<"\t\t\t\t\t\t+- " << std::fixed<<std::setprecision(3);
                for(int i = 1; i<(int)templateNames_.size(); ++i) {
                    std::cout << hfFracScaleFactorsPerChannel.second.at(sampleTypeForId(i)).err<<" |    ";
                }
                std::cout << std::endl;
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

        for(Channel::Channel channel : v_channel){
            std::vector<TH1*> histos(sampleTypeIds_.at(Sample::dummy) + 1);
            const std::vector<Sample>& v_sample(channelSamples.at(channel));
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                const Sample& sample = v_sample.at(iSample);
                const Sample::SampleType& sampleType = sample.sampleType();
                
                TH1* h = rootFileReader_->GetClone<TH1D>(sample.inputFile(), TString(histoTemplateName_).Append(step));
                h->Sumw2();
                // Removing the bin 6
                h->SetBinContent(6, 0.0);
                h->SetBinError(6, 0.0);
                
                
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
                if (histos.at(sampleId)) {
                    histos.at(sampleId)->Add(h); 
                } else {
                    sprintf(tempName, "h_%d", sampleId);
                    histos.at(sampleId) = (TH1*)h->Clone(tempName);
                }
                
            }
            
            const std::vector<HfFracScaleFactors::ValErr> sampleSFs = getScaleFactorsFromHistos(histos, step, channel, systematic);
            
            for(int type = Sample::data; type <= Sample::dummy; ++type) {
                Sample::SampleType sampleType = (sampleTypeIds_.count((Sample::SampleType)type) == 0) ? Sample::dummy : (Sample::SampleType)type;
                m_hfFracScaleFactors_[step][systematic][channel][sampleType] = sampleSFs.at(sampleTypeIds_.at(sampleType));
            }
            
            histos.clear();
            
        }
    }
}


const std::vector<HfFracScaleFactors::ValErr> HfFracScaleFactors::getScaleFactorsFromHistos(const std::vector<TH1*> histos, 
                                                                                            const TString& step, const Channel::Channel channel,
                                                                                            const Systematic::Systematic& systematic)
{
    if(histos.size()<2) {
        std::cerr<<"Error in fitting function getScaleFactorsFromHistos! Too few histograms provided: "<< histos.size() <<"\n...break\n"<<std::endl;
        exit(17);
    }
    if(histos.size() < templateNames_.size()) {
        std::cerr<<"Not all histograms to be stored for the fit have names specified. Stopping..."<<std::endl;
        exit(17);
    }
    // Initialisation of standard scale factor values
    std::vector<HfFracScaleFactors::ValErr> result(histos.size(), HfFracScaleFactors::ValErr(1., 1.));
    
    // Creating the folder structure where templates should be stored/accessed
    TString rootFileFolder = common::accessFolder(workingDirectory_, channel, systematic, true);
    if(rootFileFolder == "") rootFileFolder = common::assignFolder(workingDirectory_, channel, systematic);
    TString rootFileName = histoTemplateName_+step+".root";
    TString rootFilePath = rootFileFolder+rootFileName;
    // Checking whether file contains all the proper histograms
    std::vector<TH1*> histosInFile;
    bool sameHistogramsInFile = true;
    for(TString histoName : templateNames_) {
        TH1* histo = rootFileReader_->GetClone<TH1F>(rootFilePath, histoName, true);
        if(!histo) break;
        histosInFile.push_back(histo);
    }
    // Checking whether each histogram has the same number of entries/bins and integral
    if(histosInFile.size() == histos.size()) {
        for(size_t histoId = 0; histoId < histos.size(); ++histoId) {
            if(histogramsAreIdentical(histos.at(histoId), histosInFile.at(histoId))) continue;
            sameHistogramsInFile = false;
            break;
        }
    } else sameHistogramsInFile = false;
    // If file already contains histograms which are different: STOP
    if(!sameHistogramsInFile && histosInFile.size()>0) {
        std::cerr << "Templates used for the fit don't match to the templates in the analysis.\n";
        std::cerr << "Remove the folder: " << workingDirectory_ << " and rerun the tool to create proper templates.\n";
        exit(20);
    }
    // Storing input templates to the root file if not there already
    if(!sameHistogramsInFile) {
        TFile* out_root = new TFile(rootFilePath, "UPDATE");
        for(size_t iHisto=0; iHisto<histos.size(); ++iHisto) {
            TH1F* hF = new TH1F();
            ((TH1D*)histos.at(iHisto))->Copy(*hF);
            if(iHisto == 0) hF->Write("dummy");
            hF->Write(templateNames_.at(iHisto));
            // Adding names of original templates
            TString dataCardEntryName("FRAC_");
            dataCardEntryName+=templateNames_.at(iHisto);
            if(iHisto == 0) continue;   // Data statistical uncertainties taken directly from the original histogram
            templateVariationNameId_[dataCardEntryName] = (int)iHisto;
            if(systematic.type() != Systematic::Type::nominal) continue;
            // Storing separate versions of the template with statistical variation of each bin (for Nominal only)
            storeStatisticalTemplateVariations(hF, templateNames_.at(iHisto), iHisto);
        }
        out_root->Close();
        delete out_root;
        
        // Writing the datacard for the specified histograms
        TString datacardName(rootFilePath);
        datacardName.ReplaceAll(".root", ".txt");
        writeDatacardWithHistos(histos, datacardName, rootFileName);
    } else std::cout << "All templates already exist as input for the fit. Nothing overwritten.\n";
    
    // Opening the root file with the fit results
    TString fitFileName = rootFileFolder+histoTemplateName_+step+"/mlfittest.root";
    // Reading all fitted templates if available
    std::vector<TH1*> histosInFit(1, histos.at(0));
    for(size_t i = 1; i<templateNames_.size(); ++i) {
        TH1* histo = rootFileReader_->GetClone<TH1F>(fitFileName, "shapes_fit_b/HF/"+templateNames_.at(i), true);
        if(!histo) break;
        histosInFit.push_back(histo);
    }
    bool fitResultsInFile = histosInFit.size() == histos.size();
    if(!fitResultsInFile) {
        std::cout << "\n### Fit results don't exist in file: " << fitFileName << std::endl;
        std::cout << "### Run the fit first:\n";
        std::cout << "    source scripts/fitHfFrac.sh\n";
        std::cout << "### Then run this tool again\n\n";
        return result;
    }
//     TH1* histoSum = rootFileReader_->GetClone<TH1F>(fitFileName, "shapes_fit_b/HF/total_background", true);
    
    // Extracting fit result for each sampleType
    for(size_t sampleId = 1; sampleId<templateNames_.size(); ++sampleId) {
        double v(0.), e(0.);
        const int binId = histosInFit.at(sampleId)->GetMaximumBin();
        v = histosInFit.at(sampleId)->GetBinContent(binId);
        e = histosInFit.at(sampleId)->GetBinError(binId);
        v /= histos.at(sampleId)->GetBinContent(binId);
        e /= histos.at(sampleId)->GetBinContent(binId);
        result.at(sampleId).val = v;
        result.at(sampleId).err = e;
    }
    
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


void HfFracScaleFactors::writeDatacardWithHistos(const std::vector<TH1*> histos, const TString& fileName, const TString& rootFileName)const
{
    if(histos.size() < templateNames_.size()) {
        std::cerr<<"Not all histograms to be stored for the fit have names specified. Stopping..."<<std::endl;
        exit(17);
    }
    
    const int nHistos = histos.size();
    
    ofstream file (fileName);
    
    file << "imax 1   number of channels\n";
    file << "jmax *   number of backgrounds\n";
    file << "kmax *   number of nuisance parameters\n";
    file << "----------------------------------------\n";
    file << "observation " << histos.at(0)->Integral() << "\n";
    file << "---------------------------------------------------------------\n";
    file << "shapes * * " << rootFileName << " $PROCESS    $PROCESS_$SYSTEMATIC\n";
    file << "---------------------------------------------------------------\n";
    file << "bin      \tHF";
    for(int i = 1; i<nHistos; ++i) file << "\tHF";
    file << "\n";
    file << "process  \tdummy";
    for(int i = 1; i<nHistos; ++i) file << "\t" << templateNames_.at(i);
    file << "\n";
    file << "process  \t0";
    for(int i = 1; i<nHistos; ++i) file << "\t" << i;
    file << "\n";
    file << "rate     \t" << histos.at(0)->Integral();
    for(int i = 1; i<nHistos; ++i) file << "\t" << histos.at(i)->Integral();
    file << "\n";
    file << "---------------------------------------------------------------\n";
    // Adding nuisance parameters (fraction of each temlpate) with variation range
//     for(int i = 1; i<nHistos; ++i) {
    for(auto nameId : templateVariationNameId_) {
        // Normal distribution for the fixed background. Uniform for fit parameters
        std::string dependence;
        float factor;
        if(nameId.first.BeginsWith("FRAC_")) {
            dependence = nameId.second < nHistos - 1 ? "lnU" : "lnN";
            factor = templateScaleLimits_.at(nameId.second);
        }
        else if(nameId.first.BeginsWith("STAT_")) {
            dependence = "shape";
            factor = 0.5;
        }
        file << nameId.first << "\t" << dependence << " -";
        // Limits of scale variation of the template
        for(int j = 1; j < nHistos; ++j) {
            file << "\t";
            if(nameId.second == j || dependence == "shape") file << factor;
            else file << "-";
        }
        file << "\n";
    }
    file << "---------------------------------------------------------------\n";
}


bool HfFracScaleFactors::histogramsAreIdentical(TH1* histo1, TH1* histo2)const
{
    if(histo1->GetEntries() != histo2->GetEntries()) return false;
    if(histo1->GetNbinsX() != histo2->GetNbinsX()) return false;
    if(std::fabs((histo1->Integral() - histo2->Integral())/histo1->Integral()) > 0.00001) return false;
    
    return true;
}


Sample::SampleType HfFracScaleFactors::sampleTypeForId(const int id)const
{
    for (std::map<Sample::SampleType, int>::const_iterator it = sampleTypeIds_.begin(); it != sampleTypeIds_.end(); ++it )
    if (it->second == id) return it->first;
    return Sample::dummy;
}

int HfFracScaleFactors::storeStatisticalTemplateVariations(const TH1* histo, const TString& name, const int templateId)
{
    int nHistosStored = 0;
    const int nBins = histo->GetNbinsX();
    for(int dirId = 0; dirId < 2; ++dirId) {
        const TString dirStr = dirId == 0 ? "Up" : "Down";
        const float dirFactor = dirId == 0 ? 1.0 : -1.0;
        for(int binId = 1; binId <= nBins; ++binId) {
            const float content = histo->GetBinContent(binId);
            const float error = histo->GetBinError(binId);
            char histoName[150];
            char dataCardName[150];
            sprintf(histoName, "%s_STAT_bin%d%s", name.Data(), binId, dirStr.Data());
            sprintf(dataCardName, "STAT_bin%d", binId);
            TH1* histoVar = (TH1*)histo->Clone();
            histoVar->SetBinContent(binId, content + dirFactor*error);
            histoVar->Write(histoName);
            nHistosStored++;
            if(templateVariationNameId_.count(dataCardName) > 0) continue;
            templateVariationNameId_[dataCardName] = templateId;
        }
    }
    
    return nHistosStored;
}
