#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>
#include <TColorWheel.h>
#include <TH1.h>

#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"





/// Input base for the file lists containing the samples to be processed
constexpr const char* FileListBASE = "FileLists_plot/HistoFileList_";





Sample::Sample():
legendEntry_(""),
color_(0),
crossSection_(0),
sampleType_(dummy),
finalState_(Channel::undefined),
systematic_(Systematic::undefined),
inputFileName_(""),
luminosityWeightPerInversePb_(-999.)
{}



Sample::Sample(const TString legendEntry, const int color, const double crossSection, const SampleType sampleType):
legendEntry_(legendEntry),
color_(color),
crossSection_(crossSection),
sampleType_(sampleType),
finalState_(Channel::undefined),
systematic_(Systematic::undefined),
inputFileName_(""),
luminosityWeightPerInversePb_(-999.)
{}



std::vector<std::pair<TString, Sample> > Samples::setSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Define all samples as differential as they are needed
    // Samples with same legend will appear as one single sample (requires also same colour)
    Sample data("Data", kBlack, 1., Sample::data);
    Sample ttbarsignalPlusBbbar("t#bar{t}b#bar{b}", 18, 234.0, Sample::ttbb);
    Sample ttbarsignalPlusB("t#bar{t}b", 12, 234.0, Sample::ttb);
    Sample ttbarsignalPlusOther("t#bar{t}Other", 23, 234.0, Sample::ttother);
    Sample ttbarbkg("t#bar{t}Other", 23, 234.0, Sample::ttother);
    Sample singletop("Single Top", kViolet-3, 11.1);
    Sample ww("Diboson", 10, 54.838);
    Sample wz("Diboson", 10, 33.21);
    Sample zz("Diboson", 10, 17.654);
    Sample dyee1050("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 860.5, Sample::dyee);
    Sample dyee50inf("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 3532.8, Sample::dyee);
    Sample dymumu1050("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 860.5, Sample::dymumu);
    Sample dymumu50inf("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 3532.8, Sample::dymumu);
    Sample dytautau1050("Z / #gamma* #rightarrow #tau#tau", kAzure+10, 860.5, Sample::dytautau);
    Sample dytautau50inf("Z / #gamma* #rightarrow #tau#tau", kAzure+10, 3532.8, Sample::dytautau);
    Sample wlnu("W+Jets", kSpring+2, 36257.2);
    Sample qcdmu15("QCD Multijet", kOrange-2, 3.640E8*3.7E-4);
    Sample qcdmu2030("QCD Multijet", kOrange-2, 2.870E8*6.500E-3);
    Sample qcdmu3050("QCD Multijet", kOrange-2, 6.609E7*12.20E-3);
    Sample qcdmu5080("QCD Multijet", kOrange-2, 8.802E6*21.80E-3);
    Sample qcdmu80120("QCD Multijet", kOrange-2, 1.024E6*39.50E-3);
    Sample qcdmu120170("QCD Multijet", kOrange-2, 1.578E5*47.30E-3);
    Sample qcdem2030("QCD Multijet", kOrange-2, 2.886E8*10.10E-3);
    Sample qcdem3080("QCD Multijet", kOrange-2, 7.433E7*62.10E-3);
    Sample qcdem80170("QCD Multijet", kOrange-2, 1.191E6*153.9E-3);
    Sample qcdbcem2030("QCD Multijet", kOrange-2, 2.886E8*5.800E-4);
    Sample qcdbcem3080("QCD Multijet", kOrange-2, 7.424E7*2.250E-3);
    Sample qcdbcem80170("QCD Multijet", kOrange-2, 1.191E6*10.90E-3);
    Sample ttbarW("t#bar{t}W", kViolet-4, 0.232);
    Sample ttbarZ("t#bar{t}Z", kTeal+1, 0.2057);
    Sample ttbarH125inclusiveOther("t#bar{t}H Other", kTeal+3, 0.1293, Sample::ttHother);
    Sample ttbarH125inclusiveBbbar("t#bar{t}H (b#bar{b} via incl.)", kSpring+9, 0.1293, Sample::ttHbb);
    Sample ttbarH125tobbbar("t#bar{t}H (b#bar{b})", 2, 0.1293*0.577, Sample::ttHbb);
    Sample ttbarH110inclusiveOther("t#bar{t}H110 Other", kTeal+4, 0.1871, Sample::ttHother);
    Sample ttbarH110inclusiveBbbar("t#bar{t}H110 (b#bar{b} via incl.)", kSpring+10, 0.1871, Sample::ttHbb);
    Sample ttbarH110tobbbar("t#bar{t}H110 (b#bar{b})", 3, 0.1871*0.744, Sample::ttHbb);
    Sample ttbarH115inclusiveOther("t#bar{t}H115 Other", kTeal+5, 0.1651, Sample::ttHother);
    Sample ttbarH115inclusiveBbbar("t#bar{t}H115 (b#bar{b} via incl.)", kSpring+11, 0.1651, Sample::ttHbb);
    Sample ttbarH115tobbbar("t#bar{t}H115 (b#bar{b})", 4, 0.1651*0.703, Sample::ttHbb);
    Sample ttbarH120inclusiveOther("t#bar{t}H120 Other", kTeal+6, 0.1459, Sample::ttHother);
    Sample ttbarH120inclusiveBbbar("t#bar{t}H120 (b#bar{b} via incl.)", kSpring+12, 0.1459, Sample::ttHbb);
    Sample ttbarH120tobbbar("t#bar{t}H120 (b#bar{b})", 5, 0.1459*0.648, Sample::ttHbb);
    Sample ttbarH1225inclusiveOther("t#bar{t}H122.5 Other", kTeal+7, 0.1373, Sample::ttHother);
    Sample ttbarH1225inclusiveBbbar("t#bar{t}H122.5 (b#bar{b} via incl.)", kSpring+13, 0.1373, Sample::ttHbb);
    Sample ttbarH1275inclusiveOther("t#bar{t}H127.5 Other", kTeal+8, 0.1218, Sample::ttHother);
    Sample ttbarH1275inclusiveBbbar("t#bar{t}H127.5 (b#bar{b} via incl.)", kSpring+14, 0.1218, Sample::ttHbb);
    Sample ttbarH130inclusiveOther("t#bar{t}H130 Other", kTeal+9, 0.1149, Sample::ttHother);
    Sample ttbarH130inclusiveBbbar("t#bar{t}H130 (b#bar{b} via incl.)", kSpring+15, 0.1149, Sample::ttHbb);
    Sample ttbarH130tobbbar("t#bar{t}H130 (b#bar{b})", 6, 0.1149*0.494, Sample::ttHbb);
    Sample ttbarH135inclusiveOther("t#bar{t}H135 Other", kTeal+10, 0.1024, Sample::ttHother);
    Sample ttbarH135inclusiveBbbar("t#bar{t}H135 (b#bar{b} via incl.)", kSpring+16, 0.1024, Sample::ttHbb);
    Sample ttbarH135tobbbar("t#bar{t}H135 (b#bar{b})", 7, 0.1024*0.404, Sample::ttHbb);
    Sample ttbarH140inclusiveOther("t#bar{t}H140 Other", kTeal+11, 0.09150, Sample::ttHother);
    Sample ttbarH140inclusiveBbbar("t#bar{t}H140 (b#bar{b} via incl.)", kSpring+17, 0.09150, Sample::ttHbb);
    
    
    // Access fileList containing list of input root files
    const TString histoListName(FileListBASE + Systematic::convertSystematic(systematic) + "_" + Channel::convertChannel(channel) + ".txt");
    std::cout << "Reading file: " << histoListName << std::endl;
    ifstream fileList(histoListName);
    if (fileList.fail()) {
        std::cerr << "Error reading file: " << histoListName << std::endl;
        exit(1);
    }
    
    // Read in fileList to a vector
    std::vector<TString> v_filename;
    while(!fileList.eof()){
        TString filename;
        fileList>>filename;
        if(filename == ""){continue;} // Skip empty lines
        if(filename.BeginsWith("#")){continue;} // Comment lines in FileList with '#'
        v_filename.push_back(filename);
    }
    
    // Fill vector with samples defined above, based on their filename
    // Samples will appear in the order defined here (however samples with same legend will be combined in blocks later)
    std::vector<std::pair<Sample, std::vector<TString> > > v_sampleNamepatternsPair;
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(data, {"run"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdmu15, {"qcdmu15"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdmu2030, {"qcdmu2030"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdmu3050, {"qcdmu3050"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdmu5080, {"qcdmu5080"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdmu80120, {"qcdmu80120"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdmu120170, {"qcdmu120170"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdem2030, {"qcdem2030"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdem3080, {"qcdem3080"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdem80170, {"qcdem80170"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdbcem2030, {"qcdbcem2030"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdbcem3080, {"qcdbcem3080"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(qcdbcem80170, {"qcdbcem80170"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(wlnu, {"wtolnu"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(dyee1050, {"dyee", "1050"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(dyee50inf, {"dyee", "50inf"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(dymumu1050, {"dymumu", "1050"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(dymumu50inf, {"dymumu", "50inf"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(dytautau1050, {"dytautau", "1050"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(dytautau50inf, {"dytautau", "50inf"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(singletop, {"single"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ww, {"ww"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(wz, {"wz"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(zz, {"zz"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarbkg, {"ttbarbg"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarsignalPlusOther, {"ttbarsignal", "PlusOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarsignalPlusB, {"ttbarsignal", "PlusB."}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarsignalPlusBbbar, {"ttbarsignal", "PlusBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarW, {"ttbarW"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarZ, {"ttbarZ"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH125inclusiveOther, {"ttbarH125inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH125inclusiveBbbar, {"ttbarH125inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH125tobbbar, {"ttbarH125tobbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH110inclusiveOther, {"ttbarH110inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH110inclusiveBbbar, {"ttbarH110inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH110tobbbar, {"ttbarH110tobbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH115inclusiveOther, {"ttbarH115inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH115inclusiveBbbar, {"ttbarH115inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH115tobbbar, {"ttbarH115tobbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH120inclusiveOther, {"ttbarH120inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH120inclusiveBbbar, {"ttbarH120inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH120tobbbar, {"ttbarH120tobbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH1225inclusiveOther, {"ttbarH1225inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH1225inclusiveBbbar, {"ttbarH1225inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH1275inclusiveOther, {"ttbarH1275inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH1275inclusiveBbbar, {"ttbarH1275inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH130inclusiveOther, {"ttbarH130inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH130inclusiveBbbar, {"ttbarH130inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH130tobbbar, {"ttbarH130tobbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH135inclusiveOther, {"ttbarH135inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH135inclusiveBbbar, {"ttbarH135inclusiveBbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH135tobbbar, {"ttbarH135tobbbar"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH140inclusiveOther, {"ttbarH140inclusiveOther"}));
    v_sampleNamepatternsPair.push_back(std::pair<Sample, std::vector<TString> >(ttbarH140inclusiveBbbar, {"ttbarH140inclusiveBbbar"}));
    
    // Return the input files and corresponding samples, ordered by the name patterns
    return this->samplesByNamePatterns(v_filename, v_sampleNamepatternsPair);
}








// ------------------------------------ Further methods of class Sample ---------------------------------------



TString Sample::legendEntry()const{return legendEntry_;}



int Sample::color()const{return color_;}



double Sample::crossSection()const{return crossSection_;}



Sample::SampleType Sample::sampleType()const{return sampleType_;}



void Sample::setFinalState(const Channel::Channel& finalState){finalState_ = finalState;}



Channel::Channel Sample::finalState()const{return finalState_;}



void Sample::setSystematic(const Systematic::Systematic& systematic){systematic_ = systematic;}



Systematic::Systematic Sample::systematic()const{return systematic_;}



void Sample::setInputFile(const TString& inputFileName)
{
    // Check if input file really exists
    ifstream ifile(inputFileName);
    if(ifile.fail()){
        std::cerr<<"Cannot find requested file: "<<inputFileName<<"\n...breaking\n";
        exit(17);
    }
    inputFileName_ = inputFileName;
}



TString Sample::inputFile()const{return inputFileName_;}



double Sample::luminosityWeight(const double& luminosityInInversePb)const
{
    if(sampleType_ == data) return 1.;
    if(luminosityWeightPerInversePb_ < 0.){
        std::cerr<<"ERROR in Sample::luminosityWeight()! Value is negative (probably not calculated): "
                 <<luminosityWeightPerInversePb_<<"\n...break\n"<<std::endl;
        exit(886);
    }
    return luminosityInInversePb*luminosityWeightPerInversePb_;
}



void Sample::calculateLuminosityWeight()
{
    if(sampleType_ == SampleType::data) return;
    if(crossSection_ <= 0.){
        std::cerr<<"ERROR: Sample XSection is <0. Can't calculate luminosity weight!!\n...break\n"<<std::endl;
        exit(887);
    }
    
    RootFileReader* fileReader = RootFileReader::getInstance();
    const TH1* const h_weightedEvents = fileReader->Get<TH1>(inputFileName_, "weightedEvents");
    const double weightedEvents(h_weightedEvents->GetBinContent(1));
    
    luminosityWeightPerInversePb_ = crossSection_/weightedEvents;
    //std::cout<<"Input file: "<<inputFileName_<<std::endl;
    //std::cout<<"Xsection, weighted events, lumi weight: "
    //         <<crossSection<<" , "<<weightedEvents<<" , "<<luminosityWeightPerInversePb_<<std::endl;
}







// ------------------------------------ Further methods of class Samples ---------------------------------------



Samples::Samples():
globalScaleFactors_(0)
{}



Samples::Samples(const std::vector<Channel::Channel>& v_channel,
                 const std::vector<Systematic::Systematic>& v_systematic,
                 const GlobalScaleFactors* globalScaleFactors):
globalScaleFactors_(globalScaleFactors)
{
    std::cout<<"--- Beginning to set up the samples\n\n";
    
    for(const auto& systematic : v_systematic){
        for(const auto& channel : v_channel){
            this->addSamples(channel, systematic);
        }
    }
    
    std::cout<<"\n=== Finishing to set up the samples\n\n";
}



void Samples::setGlobalWeights(const GlobalScaleFactors* globalScaleFactors)
{
    globalScaleFactors_ = globalScaleFactors;
}



std::vector<std::pair<TString, Sample> > Samples::samplesByNamePatterns(const std::vector<TString>& v_filename,
                                                                        const std::vector<std::pair<Sample, std::vector<TString> > >& v_sampleNamepatternsPair)
{
    std::vector<std::pair<TString, Sample> > result;
    
    std::vector<size_t> v_selectedFileIndex;
    for(const auto& sampleNamepatternsPair : v_sampleNamepatternsPair){
        const Sample& sample = sampleNamepatternsPair.first;
        const std::vector<TString>& v_namepattern = sampleNamepatternsPair.second;
        if(!v_namepattern.size()){
            std::cerr<<"ERROR in Samples::sampleByNamePatterns()! No name patterns are specified\n...break\n"<<std::endl;
            exit(512);
        }
        
        for(size_t iFile = 0; iFile < v_filename.size(); ++iFile){
            const TString& filename = v_filename.at(iFile);
            
            bool allPatternsContained(true);
            for(const TString& namePattern : v_namepattern){
                if(!filename.Contains(namePattern)){
                    allPatternsContained = false;
                    break;
                }
            }
            
            if(allPatternsContained){
                if(std::find(v_selectedFileIndex.begin(), v_selectedFileIndex.end(), iFile) != v_selectedFileIndex.end()){
                    std::cerr<<"ERROR in Samples::sampleByNamePatterns()! Input file selected twice for different samples: "
                             <<filename<<"\n...break\n"<<std::endl;
                    exit(513);
                }
                v_selectedFileIndex.push_back(iFile);
                result.push_back(std::make_pair(filename, sample));
            }
        }
    }
    
    if(result.size() != v_filename.size()){
        std::cout<<"WARNING: Not all samples of input file list will be used (# input, # usage): "<<v_filename.size()<<" , "<<result.size()<<"\n";
    }
    
    return result;
}



void Samples::addSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Add all samples as they are defined in setSamples()
    std::vector<std::pair<TString, Sample> > v_filenameSamplePair(this->setSamples(channel, systematic));
    
    // Set sample options via filename
    std::vector<Sample> v_sample(this->setSampleOptions(systematic, v_filenameSamplePair));
    
    // Order files by legendEntry
    this->orderByLegend(v_sample);
    
    // Create map of maps, containing Sample per channel per systematic
    m_systematicChannelSample_[systematic][channel] = v_sample;
    
    //for(auto systematicChannelSample : m_systematicChannelSample_){
    //    for(auto channelSample : systematicChannelSample.second){
    //        for(auto sample : channelSample.second){
    //            std::cout<<"We have samples: "<<Systematic::convertSystematic(systematic)
    //                <<" , "<<Channel::convertChannel(channel)
    //                <<" , "<<sample.inputFile()
    //                <<" , "<<Channel::convertChannel(sample.finalState())<<"\n";
    //        }
    //    }
    //}
}



std::vector<Sample> Samples::setSampleOptions(const Systematic::Systematic& systematic,
                                              const std::vector<std::pair<TString, Sample> >& v_filenameSamplePair)
{
    std::vector<Sample> v_sample;

    for(auto filenameSamplePair : v_filenameSamplePair){
        TString filename(filenameSamplePair.first);
        Sample sample(filenameSamplePair.second);
        
        // Assign dilepton final state to each sample
        sample.setFinalState(this->assignFinalState(filename));
        
        // Assign specific systematic to each sample and adjust filename accordingly
        sample.setSystematic(this->assignSystematic(filename, systematic));
        
        // Check if input file really exists and set it
        sample.setInputFile(filename);
        
        // Calculate the luminosity weight of the sample
        sample.calculateLuminosityWeight();
        
        v_sample.push_back(sample);
    }

    return v_sample;
}



void Samples::orderByLegend(std::vector<Sample>& v_sample)
{
    // Associate vector constituents to legend entry
    // and store all unequal legend entries
    std::vector<std::pair<TString, Sample> > v_legendSamplePair;
    std::vector<TString> v_legendEntry;
    for(auto sample : v_sample){
        const TString& legendEntry(sample.legendEntry());
        //std::cout<<"Legends before: "<<legendEntry<<std::endl;
        v_legendSamplePair.push_back(std::pair<TString, Sample>(legendEntry, sample));
        if(std::find(v_legendEntry.begin(), v_legendEntry.end(), legendEntry) != v_legendEntry.end())continue;
        else v_legendEntry.push_back(legendEntry);
    }

    // Clear vector and fill it again in correct order
    v_sample.clear();
    for(auto legendEntry : v_legendEntry){
        int testColor(-999);
        for(auto legendSamplePair : v_legendSamplePair){
             if(legendSamplePair.first != legendEntry)continue;

             //std::cout<<"Legends after: "<<legendEntry<<std::endl;
             v_sample.push_back(legendSamplePair.second);

             // Check if all with same legend do have same colour
             int color = legendSamplePair.second.color();
             if(testColor == -999)testColor = color;
             else if(testColor != color){
                 std::cerr<<"ERROR! Samples have different colors but same legends: "<<legendEntry
                          <<"\n...break\n";
                 exit(999);
             }
         }
    }
}



Channel::Channel Samples::assignFinalState(const TString& filename)
{
    std::vector<Channel::Channel> v_channel {Channel::ee, Channel::emu, Channel::mumu};
    for(auto channel : v_channel){
        TString finalState(Channel::convertChannel(channel));
        finalState.Prepend("/");
        finalState.Append("/");
        if(filename.Contains(finalState)){
            return channel;
        }
    }
    return Channel::undefined;
}



Systematic::Systematic Samples::assignSystematic(TString&, const Systematic::Systematic&)
//Systematic::Systematic Samples::assignSystematic(TString& filename, const Systematic::Systematic& systematic)
{
    // FIXME: adjust filename corresponding to specific systematic

    return Systematic::undefined;
}



const SystematicChannelSamples& Samples::getSystematicChannelSamples()const
{
    return m_systematicChannelSample_;
}



const std::vector<Sample>& Samples::getSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)const
{
    if(m_systematicChannelSample_.find(systematic) == m_systematicChannelSample_.end()){
        std::cerr<<"ERROR in getSamples! No samples available for requested systematic: "<<Systematic::convertSystematic(systematic)
                 <<"\n...break\n"<<std::endl;
        exit(321);
    }
    if(m_systematicChannelSample_.at(systematic).find(channel) == m_systematicChannelSample_.at(systematic).end()){
        std::cerr<<"ERROR in getSamples! No samples available for requested channel: "<<Channel::convertChannel(channel)
                 <<"\n...break\n"<<std::endl;
        exit(322);
    }
    return m_systematicChannelSample_.at(systematic).at(channel);
}



std::pair<SystematicChannelFactors, bool> Samples::globalWeights(const TString& objectname,
                                                                 const bool dyCorrection,
                                                                 const bool ttbbCorrection)const
{
    const TString fullStepname = tth::extractSelectionStepAndJetCategory(objectname);
    
    return globalScaleFactors_->scaleFactors(*this, fullStepname, dyCorrection, ttbbCorrection);
}




std::pair<SystematicChannelFactors, bool> Samples::globalWeights(const TString& objectname)const
{
    const TString fullStepname = tth::extractSelectionStepAndJetCategory(objectname);
    
    return globalScaleFactors_->scaleFactors(*this, fullStepname);
}







