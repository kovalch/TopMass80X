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
#include <TKey.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>


#include "HfFracScaleFactors.h"
#include "higgsUtils.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"





HfFracScaleFactors::HfFracScaleFactors(const Samples& samples, RootFileReader* const rootFileReader):
rootFileReader_(rootFileReader),
workingDirectory_("HfFracScaleFactors"),
scaleFactorsUsable_(true)
{
    std::cout<<"--- Beginning production of Heavy-Flavour fraction scale factors\n\n";
    
    // Setting name of the histogram used for template fit
    histoTemplateName_ = "hfFracScaling_btag_multiplicity";
//     histoTemplateName_ = "hfFracScaling_jetCategories";
    
    // Setting id for each sample type as it will appear in the list of histograms for the fit
    // Data MUST go first
    sampleTypeIds_[Sample::data] = 0;
    // Combine samples by assigning the same id
    sampleTypeIds_[Sample::ttbb] = 1;
    sampleTypeIds_[Sample::ttb] = 1;
    sampleTypeIds_[Sample::tt2b] = 1;
    sampleTypeIds_[Sample::ttcc] = 2;
    sampleTypeIds_[Sample::ttother] = 2;
    // Backgrounds MUST go last
    sampleTypeIds_[Sample::dummy] = 3;
    
    // Setting names to each template
    templateNames_.push_back("data_obs");
    templateNames_.push_back("tthf");
//     templateNames_.push_back("tt2b");
//     templateNames_.push_back("ttcc");
    templateNames_.push_back("ttOther");
    templateNames_.push_back("bkg");

    // Setting variation limit of for each template
    templateScaleLimits_.push_back(1.);	  // Doesn't mean anything
    templateScaleLimits_.push_back(100.);
//     templateScaleLimits_.push_back(100.);
//     templateScaleLimits_.push_back(100.);
    templateScaleLimits_.push_back(100.);
    templateScaleLimits_.push_back(1.1);
    
    // Initial scale factors for each template (to check dependence on starting parameters)
    templateInitialScaleFactors_.push_back(1.);
    templateInitialScaleFactors_.push_back(1.);
//     templateInitialScaleFactors_.push_back(1.);
//     templateInitialScaleFactors_.push_back(1.);
    templateInitialScaleFactors_.push_back(1.);
    templateInitialScaleFactors_.push_back(1.);
    
    // Prescale of the tt2b and ttcc templates to be varied as a systematic
    sampleTypePrescale_[Sample::tt2b] = 1.;
    sampleTypePrescale_[Sample::ttcc] = 1.;
    
    // Setting up systematics to be used for each template
    templateSystematics_["tthf"] = std::vector<Systematic::Type>(0);
    templateSystematics_.at("tthf").push_back(Systematic::jes);
    templateSystematics_.at("tthf").push_back(Systematic::btagDiscrPurity);
    templateSystematics_.at("tthf").push_back(Systematic::btagDiscrBstat1);
    templateSystematics_.at("tthf").push_back(Systematic::btagDiscrBstat2);
    
    templateSystematics_["ttOther"] = std::vector<Systematic::Type>(0);
    templateSystematics_.at("ttOther").push_back(Systematic::jes);
    templateSystematics_.at("ttOther").push_back(Systematic::btagDiscrPurity);
    templateSystematics_.at("ttOther").push_back(Systematic::btagDiscrBstat1);
    templateSystematics_.at("ttOther").push_back(Systematic::btagDiscrBstat2);
    templateSystematics_.at("ttOther").push_back(Systematic::btagDiscrLstat1);
    templateSystematics_.at("ttOther").push_back(Systematic::btagDiscrLstat2);
    
    // FIXME: Check that there are no gaps in the list of ids
    
    this->produceScaleFactors(samples);
    
    if(!scaleFactorsUsable_) {
        std::cerr << "\n### Not all input has been provided to produce usable scale factors. Stopping..." << std::endl;
        exit(1);
    }
    
    // Hiding all standard ROOT warnings
//     gErrorIgnoreLevel = kWarning;
    
    // Suppress default info that canvas is printed
    gErrorIgnoreLevel = 1001;
    
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
        if(nameStepPair.first != histoTemplateName_+nameStepPair.second) continue;
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
    const SystematicChannelFactors globalWeights = samples.globalWeights(step, false, false).first;
    
    for(auto systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        const auto& channelSamples(systematicChannelSamples.second);
        
//         const std::vector<Channel::Channel> v_channel = {Channel::ee, Channel::emu, Channel::mumu, Channel::combined};
        const std::vector<Channel::Channel> v_channel = {Channel::combined};

        for(Channel::Channel channel : v_channel){
            std::vector<TH1*> histos(sampleTypeIds_.at(Sample::dummy) + 1);
            if(channelSamples.count(channel) < 1) {
                std::cerr << "\nERROR! tt+HF scale factors have to be determined at all channels combined.\n"; 
                std::cerr << "PLEASE, include COMBINED channel in the list of processed channels.\n";
                exit(3);
            }
            const std::vector<Sample>& v_sample(channelSamples.at(channel));
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                const Sample& sample = v_sample.at(iSample);
                const Sample::SampleType& sampleType = sample.sampleType();
                
                if(sampleType == Sample::pseudodata) {
                    std::cerr << "\nERROR! Currently pseudodata can't be used to determine tt+HF scale factor. Stopping...\n";
                    exit(3);
                }
                
                TH1* h = rootFileReader_->GetClone<TH1D>(sample.inputFile(), TString(histoTemplateName_).Append(step));
                if(h->GetSumw2N() < 1) h->Sumw2();
                
                // Prescaling the sample if needed
                if(sampleTypePrescale_.count(sampleType)) h->Scale(sampleTypePrescale_.at(sampleType));
                
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


const std::vector<HfFracScaleFactors::ValErr> HfFracScaleFactors::getScaleFactorsFromHistos(std::vector<TH1*> histos, 
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
    
    // Rebinning the histograms to have no empty bins
    Double_t bins[6] = {1.0, 2.0, 3.0, 4., 5., 6.};      // btag multiplicity
    for(size_t iHisto = 0; iHisto < histos.size(); ++iHisto) {
        TH1* histo = histos.at(iHisto);
        if(!histo) continue;
        // Applying the initial scale factor for the fit
        histo->Scale(templateInitialScaleFactors_.at(iHisto));
        histos.at(iHisto) = histo->Rebin(5, TString(histo->GetName())+="_rebin", bins);
    }
    
    // Creating the folder structure where templates should be stored/accessed
    TString rootFileFolder = common::accessFolder(workingDirectory_, channel, Systematic::convertType(Systematic::Type::nominal), true);
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
    // Checking whether each histogram has the same number of entries/bins and integral for Nominal templates
    if(histosInFile.size() == histos.size()) {
        if(systematic.type() == Systematic::Type::nominal) {
            for(size_t histoId = 0; histoId < histos.size(); ++histoId) {
                if(histogramsAreIdentical(histos.at(histoId), histosInFile.at(histoId))) continue;
                sameHistogramsInFile = false;
                break;
            }
        }
    } else sameHistogramsInFile = false;
    // If file already contains histograms which are different: STOP
    if(!sameHistogramsInFile && histosInFile.size()>0) {
        std::cerr << "Templates used for the fit don't match to the templates in the analysis.\n";
        std::cerr << "Remove the folder: " << workingDirectory_ << " and rerun the tool to create proper templates.\n";
        exit(20);
    }
    // Storing input templates to the root file if not there already
    if(!sameHistogramsInFile || systematic.type() != Systematic::Type::nominal) {
        TFile* out_root = new TFile(rootFilePath, "UPDATE");
        for(size_t iHisto=0; iHisto<histos.size(); ++iHisto) {
            // Converting histogram to TH1F
            TH1F* hF = new TH1F();
            if(histos.at(iHisto)) ((TH1D*)histos.at(iHisto))->Copy(*hF);
            if(systematic.type() == Systematic::Type::nominal) {
                // Storing Nominal templates with corresponding fit parameters
                if(iHisto == histos.size()-1) hF->Write("dummy");
                hF->Write(templateNames_.at(iHisto));
                // Adding names of original templates
                TString dataCardEntryName("FRAC_");
                dataCardEntryName+=templateNames_.at(iHisto);
                if(iHisto == 0) continue;   // Data is not a fitting parameter
                templateVariationNameId_[dataCardEntryName] = (int)iHisto;
                // Storing separate versions of the template with statistical variation of each bin
                storeStatisticalTemplateVariations(hF, templateNames_.at(iHisto), iHisto);
            } else {
                // Skipping sample if the current systematic is not enabled for it
                if(!templateNameHasSystematic(templateNames_.at(iHisto), systematic.type())) continue;
                // Storing systematic templates for each bin variation and adding as shape nuisance parameters
                storeSystematicTemplateVariations(hF, systematic, templateNames_.at(iHisto), iHisto);
            }
            delete hF;
        }
        out_root->Close();
        delete out_root;
        
        // Writing the datacard for the specified histograms
        TString datacardName(rootFilePath);
        datacardName.ReplaceAll(".root", ".txt");
        writeDatacardWithHistos(histos, datacardName, rootFileName, systematic.type());
        plotInputTemplates(rootFilePath);
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
        scaleFactorsUsable_ = false;
        return result;
    }
//     TH1* histoSum = rootFileReader_->GetClone<TH1F>(fitFileName, "shapes_fit_b/HF/total_background", true);
    
    // Extracting fit result for each sampleType
    for(size_t sampleId = 1; sampleId<templateNames_.size(); ++sampleId) {
        double v(0.), e(0.);
        TH1* histoInFit = histosInFit.at(sampleId);
        v = histoInFit->Integral();
        // Error is a linear sum of errors in each bin (as fully correlated)
        for(int iBin = 1; iBin <= histoInFit->GetNbinsX(); ++iBin) e+=histoInFit->GetBinError(iBin);
        v /= histos.at(sampleId)->Integral()/templateInitialScaleFactors_.at(sampleId);
        e /= histos.at(sampleId)->Integral()/templateInitialScaleFactors_.at(sampleId);
        result.at(sampleId).val = v;
        result.at(sampleId).err = e;
    }
    
    return result;
}



const HfFracScaleFactors::ValErr& HfFracScaleFactors::hfFracScaleFactor(const TString& step,
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
    
    return m_hfFracScaleFactors_.at(step).at(systematic).at(channel).at(sampleType);
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
    
    Sample::SampleType sampleType = sample.sampleType();
    if(sampleTypeIds_.count(sampleType) < 1) sampleType = Sample::SampleType::dummy;
    
    // Check whether Heavy-Flavour fraction scale factor exists for extracted step
    if(m_hfFracScaleFactors_.find(step) == m_hfFracScaleFactors_.end()) return -1;
    
    if(m_hfFracScaleFactors_.find(step) == m_hfFracScaleFactors_.end()) return -1;
    
    // Determining whether uncertainty on the scale factor should be included for the current systematic
    double uncertaintyFactor = 0.;
    if( (systematic.type() == Systematic::frac_tthf && templateNames_.at(sampleTypeIds_.at(sampleType)) == "tthf") || 
        (systematic.type() == Systematic::frac_ttother && templateNames_.at(sampleTypeIds_.at(sampleType)) == "ttOther")
    ) uncertaintyFactor = 1.;
    // Determining in which direction central value should be varied
    if(systematic.variation() == Systematic::up) uncertaintyFactor *= 1.;
    else if(systematic.variation() == Systematic::down) uncertaintyFactor *= -1.;
    else uncertaintyFactor *= 0.;
    
    // Applying the Heavy-Flavour fraction scale factor from the combined channel
    HfFracScaleFactors::ValErr scaleFactor = this->hfFracScaleFactor(step, systematic, Channel::combined, sampleType);
    weight *= scaleFactor.val + scaleFactor.err * uncertaintyFactor;
    return 1;
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


void HfFracScaleFactors::writeDatacardWithHistos(const std::vector<TH1*> histos, const TString& fileName, 
                                                 const TString& rootFileName, const Systematic::Type systematicType)const
{
    if(histos.size() < templateNames_.size()) {
        std::cerr<<"Not all histograms to be stored for the fit have names specified. Stopping..."<<std::endl;
        exit(17);
    }
    
    const int nHistos = histos.size();
    
    std::ofstream file;
    if(systematicType == Systematic::nominal) {
        file.open(fileName, std::ofstream::out);
        
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
        file << "rate     \t" << histos.at(nHistos-1)->Integral();
        for(int i = 1; i<nHistos; ++i) file << "\t" << histos.at(i)->Integral();
        file << "\n";
        file << "---------------------------------------------------------------\n";
    } else {
        file.open(fileName, std::ofstream::app);
    }
    // Adding nuisance parameters (fraction of each temlpate) with variation range
    for(auto nameId : templateVariationNameId_) {
        // Normal distribution for the fixed background. Uniform for fit parameters
        std::string dependence;
        float factor;
        TString name = nameId.first;
        const int id = nameId.second;
        bool factorForAll = false;
        if(name.BeginsWith("FRAC_")) {
            dependence = id < nHistos - 1 ? "lnU" : "lnN";
            factor = templateScaleLimits_.at(id);
        } else if(name.BeginsWith("STAT_")) {
            factorForAll = true;
            dependence = "shape";
            factor = 1.0;       // corresponds to 1 sigma variation of statistical uncertainty
        } else {
            dependence = "shape";
            factor = 1.0;       // corresponds to 1 sigma variation of systematic uncertainty
        }
        // Determining whether this systematic for this sample should be enabled
        file << nameId.first << "\t" << dependence << " -";
        // Limits of scale variation of the template
        for(int j = 1; j < nHistos; ++j) {
            const bool hasSystematic = templateNameHasSystematic(templateNames_.at(j), systematicType);
            file << "\t";
            if(id == j || (hasSystematic && dependence == "shape") || factorForAll ) file << factor;
            else file << "-";
        }
        // No systematic variation for the last template (background)
//         if(dependence == "shape" && !name.BeginsWith("STAT_")) file << "\t-";
        file << "\n";
    }
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

bool HfFracScaleFactors::templateNameHasSystematic(const TString& templateName, const Systematic::Type systematicType)const
{
    if(templateSystematics_.count(templateName) < 1) return false;
    std::vector<Systematic::Type> activeSystematics = templateSystematics_.at(templateName);
    if(std::find(activeSystematics.begin(), activeSystematics.end(), systematicType) == activeSystematics.end()) return false;
    
    return true;
}


int HfFracScaleFactors::storeStatisticalTemplateVariations(const TH1* histo, const TString& name, const int templateId)
{
    int nHistosStored = 0;
    const int nBins = histo->GetNbinsX();
    for(int dirId = 0; dirId < 2; ++dirId) {
        const TString dirStr = dirId == 0 ? "Up" : "Down";
        const float dirFactor = dirId == 0 ? 1.0 : -1.0;
        // Storing a copy of the initial histogram with a signle bin moved up/down according to its error
        for(int binId = 1; binId <= nBins; ++binId) {
            const float error = histo->GetBinError(binId);
            const float content = std::max(histo->GetBinContent(binId), 1e-30);
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


int HfFracScaleFactors::storeSystematicTemplateVariations(const TH1* histo_systematic, const Systematic::Systematic& systematic, 
                                                          const TString& name, const int templateId)
{
    int nHistosStored = 0;
    if(systematic.type() == Systematic::Type::nominal) return nHistosStored;
    
    TString systematicName = Systematic::convertType(systematic.type());
    // Determining the proper histogram name for the shape nuisance parameter in the datacard
    TString dirStr;
    if(systematic.variation() == Systematic::Variation::up) {
        dirStr = "Up";
    } else if(systematic.variation() == Systematic::Variation::down) {
        dirStr = "Down";
    } else {
        std::cout << "Specified systematic variation is neither UP nor DOWN: " << systematicName << std::endl;
        exit(12);
    }
    // Removing UP/DOWN parts of the name
    systematicName.ReplaceAll(Systematic::convertVariation(systematic.variation()), "");
    char histoName[150];
    char dataCardName[150];
    sprintf(histoName, "%s_%s%s", name.Data(), systematicName.Data(), dirStr.Data());
    sprintf(dataCardName, "%s", systematicName.Data());
    histo_systematic->Write(histoName, TObject::kOverwrite);
    nHistosStored++;
    // Adding the proper name to the datacard
    if(templateVariationNameId_.count(dataCardName) < 1 && systematic.variation() == Systematic::Variation::up) {
        templateVariationNameId_[dataCardName] = templateId;
    }
    
//     const int nBins = histo_original->GetNbinsX();
//     for(int binId = 1; binId <= nBins; ++binId) {
//         const float contentSystematic = std::max(histo_systematic->GetBinContent(binId), 1e-30);
//         char histoName[150];
//         char dataCardName[150];
//         sprintf(histoName, "%s_%s_bin%d%s", name.Data(), systematicName.Data(), binId, dirStr.Data());
//         sprintf(dataCardName, "%s_bin%d", systematicName.Data(), binId);
//         TH1* histoVar = (TH1*)histo_original->Clone();
//         histoVar->SetBinContent(binId, contentSystematic);
//         histoVar->Write(histoName, TObject::kOverwrite);
//         nHistosStored++;
//         if(templateVariationNameId_.count(dataCardName) > 0) continue;
//         // Adding to the datacard only once (for UP variations)
//         if(systematic.variation() != Systematic::Variation::up) continue;
//         templateVariationNameId_[dataCardName] = templateId;
//     }
    
    return nHistosStored;
}


void HfFracScaleFactors::plotInputTemplates(const TString rootFileName)const
{
    // Preparing nice colours
    std::vector<Color_t> colors;
    colors.push_back(1);
    colors.push_back(2);
    colors.push_back(kAzure+2);
    colors.push_back(kTeal+4);
    colors.push_back(kOrange);
    colors.push_back(3);
    
    TString outputFileName(rootFileName);
    outputFileName.ReplaceAll(".root", "");
    
    std::map<TString, TH1*> m_nameHisto;
    std::map<TString, TGraphAsymmErrors*> m_nameGraph_stat;
    std::map<TString, TGraphAsymmErrors*> m_nameGraph_statSyst;
    THStack* stack = new THStack("def", "def");
    TH1* stacksum(0);
    
    TFile* rootFile = new TFile(rootFileName, "READONLY");
    // Getting all objects from the ROOT file
    TList* list = rootFile->GetListOfKeys();
//     TKey* key;
//     TIter nextKey(rootFile->GetListOfKeys());
//     int histoId = 0;
    for(int keyId = 0; keyId < list->GetEntries(); ++keyId) {
//     while(key = (TKey*)nextKey()) {
        TKey* key = (TKey*)list->At(keyId);
        TString name(key->GetName());
        // Skipping histogram for binwise statistical variations (statistical errors taken direclty from histograms)
        if(name.Contains("_STAT_bin")) continue;
        if(name.Contains("dummy")) continue;
        TH1* histo = (TH1*)key->ReadObj();
        m_nameHisto[name] = histo;
        if(name==templateNames_.at(0)) {
        } else if(std::find(templateNames_.begin(), templateNames_.end(), name) != templateNames_.end()) {
            stack->Add(histo);
            if(stacksum) stacksum->Add(histo);
            else stacksum = (TH1*)histo->Clone("stacksum");
        }
        
//         histoId++;
    }
    
    float yMin = 1e10;
    float yMax = 1e-10;
    // Creating the graph with asymmetric errors for each sample
    int templateId = 0;
    for(TString name : templateNames_) {
        TH1* histoOriginal = m_nameHisto.at(name);
        TGraphAsymmErrors* graph_stat = new TGraphAsymmErrors(histoOriginal);
        TGraphAsymmErrors* graph_statSyst = new TGraphAsymmErrors(histoOriginal);
        const int nBins = graph_stat->GetN();
        std::vector<float> deviationsSquaredUp(nBins, 0.);
        std::vector<float> deviationsSquaredDown(nBins, 0.);
        // Adding in quadrature all systematic variations for each bin of this histogram
        for(auto nameHisto : m_nameHisto) {
            if(!nameHisto.first.Contains(name)) continue;
            for(int iBin = 0; iBin<nBins; ++iBin) {
                TH1* systHisto = (TH1*)nameHisto.second->Clone();
                // Normalizing histogram to nominal (for shape uncertainties)
                common::normalize(systHisto, histoOriginal->Integral());
                float deviation = systHisto->GetBinContent(iBin+1) - histoOriginal->GetBinContent(iBin+1);
                if(deviation > 0) {
                    deviationsSquaredUp.at(iBin) = deviationsSquaredUp.at(iBin) + deviation*deviation;
                } else if (deviation < 0) {
                    deviationsSquaredDown.at(iBin) = deviationsSquaredDown.at(iBin) + deviation*deviation;
                }
            }
        }
        // Updating errors of the graph with systematic uncertainties
        for(int iBin = 0; iBin<nBins; ++iBin) {
            float value = histoOriginal->GetBinContent(iBin+1);
            float e_stat = histoOriginal->GetBinError(iBin+1);
            float e_systSq_up = deviationsSquaredUp.at(iBin);
            float e_systSq_down = deviationsSquaredDown.at(iBin);
            float e_total_up = sqrt(e_stat*e_stat + e_systSq_up);
            float e_total_down = sqrt(e_stat*e_stat + e_systSq_down);
            graph_statSyst->SetPointEYhigh(iBin, e_total_up);
            graph_statSyst->SetPointEYlow(iBin, e_total_down);
            // Updating Y axis limits
            if(value+e_total_up > yMax) yMax = value+e_total_up;
            if(value-e_total_down < yMin) yMin = value-e_total_down;
        }
        graph_stat->SetLineColor(colors.at(templateId));
        graph_stat->SetLineStyle(1);
        graph_stat->SetMarkerColor(colors.at(templateId));
        graph_stat->SetMarkerStyle(20+templateId);
        m_nameGraph_stat[name] = graph_stat;
        
        graph_statSyst->SetLineColor(colors.at(templateId));
        graph_statSyst->SetLineStyle(2);
        graph_statSyst->SetMarkerColor(colors.at(templateId));
        graph_statSyst->SetMarkerStyle(20+templateId);
        m_nameGraph_statSyst[name] = graph_statSyst;
        
        templateId++;
    }
    
    // Setting up the canvas
    TCanvas* canvas = new TCanvas("canvas","", 800, 800);
    TLegend* legend = new TLegend(0.70,0.7,0.92,0.85);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    legend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - legend->GetNRows()*0.04);
    legend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    legend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
    legend->Clear();
    canvas->Clear();
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    canvas->SetName("");
    canvas->SetTitle("");

    // Drawing data first
    TGraph* g_data = m_nameGraph_stat.at(templateNames_.at(0));
    common::normalize(g_data);
    g_data->Draw("AP");
    legend->AddEntry(g_data, templateNames_.at(0), "pl");
    double normScale_min = 1e30;
    for(TString name : templateNames_) {
        if(name == templateNames_.at(0)) continue;
        TGraph* graph_stat = m_nameGraph_stat.at(name);
        TGraph* graph_statSyst = m_nameGraph_statSyst.at(name);
        common::normalize(graph_statSyst);
        double normScale = common::normalize(graph_stat);
        if(normScale < normScale_min) normScale_min = normScale;
        graph_statSyst->Draw("P");
        graph_stat->Draw("P");
        legend->AddEntry(graph_stat, name, "pl");
    }
    legend->Draw();
    
    yMax *= normScale_min;
    yMin *= normScale_min;
    
    canvas->Print(outputFileName+"_templates.eps");
    canvas->SetLogy();
    canvas->Print(outputFileName+"_templates_log.eps");
    
    canvas->Clear();
    
    delete canvas;
    delete stack;
    
}
