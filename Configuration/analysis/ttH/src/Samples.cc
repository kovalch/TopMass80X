#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>
#include <TColorWheel.h>
#include <TH1.h>

#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"





/// Input base for the file lists containing the samples to be processed
constexpr const char* FileListBASE = "FileLists_plot/HistoFileList_";





Sample::Sample():
legendEntry_(""), color_(0), crossSection_(0),
sampleType_(dummy), finalState_(Channel::undefined), systematic_(Systematic::undefined),
inputFileName_("")
{}



Sample::Sample(const TString legendEntry, const int color, const double crossSection, const SampleType sampleType):
legendEntry_(legendEntry), color_(color), crossSection_(crossSection),
sampleType_(sampleType), finalState_(Channel::undefined), systematic_(Systematic::undefined),
inputFileName_("")
{}



std::vector<std::pair<TString, Sample> > Samples::setSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Define all samples as differential as they are needed
    Sample data("Data", kBlack, 1., Sample::data);
    Sample ttbarsignalPlusBbbar("t#bar{t}b#bar{b}", 18, 234.0);
    Sample ttbarsignalPlusB("t#bar{t}b", 12, 234.0);
    Sample ttbarsignalPlusOther("t#bar{t}Other", 23, 234.0);
    Sample ttbarbkg("t#bar{t} Bkg", kPink-9, 234);
    Sample singletop("Single Top", kViolet-3, 11.1);
    Sample ww("Diboson", 10, 54.838);
    Sample wz("Diboson", 10, 33.21);
    Sample zz("Diboson", 10, 17.654);
    Sample dyee1050("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 860.5, Sample::dyee);
    Sample dyee50inf("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 3532.8, Sample::dyee);
    Sample dymumu1050("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 860.5, Sample::dymumu);
    Sample dymumu50inf("Z / #gamma* #rightarrow ee/#mu#mu", kAzure+2, 3532.8, Sample::dymumu);
    Sample dytautau1050("Z / #gamma* #rightarrow #tau#tau", kAzure+8, 860.5, Sample::dytautau);
    Sample dytautau50inf("Z / #gamma* #rightarrow #tau#tau", kAzure+8, 3532.8, Sample::dytautau);
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
    Sample ttbarW("t#bar{t}W", kViolet+7, 0.232);
    Sample ttbarZ("t#bar{t}Z", kTeal+1, 0.2057);
    Sample ttbarHinclusiveOther("t#bar{t}H Other", kSpring-6, 0.1302, Sample::higgssignal);
    // Sample ttbarHinclusiveBbbar("t#bar{t}H (b#bar{b} via incl.)", kSpring+9, 0.1302, Sample::higgssignal);
    Sample ttbarHtobbbar("t#bar{t}H (b#bar{b})", 2, 0.1302*0.577, Sample::higgssignal);


    // Access FileList containing list of input root files
    const TString histoListName(FileListBASE + Systematic::convertSystematic(systematic) + "_" + Channel::convertChannel(channel) + ".txt");
    std::cout << "Reading file: " << histoListName << std::endl;
    ifstream fileList(histoListName);
    if (fileList.fail()) {
        std::cerr << "Error reading file: " << histoListName << std::endl;
        exit(1);
    }


    // Fill vector with samples defined above, based on their filename
    std::vector<std::pair<TString, Sample> > v_filenameSamplePair;
    while(!fileList.eof()){
        TString filename;
        fileList>>filename;
        if(filename==""){continue;} // Skip empty lines
        if(filename.BeginsWith("#")){continue;} // Comment lines in FileList with '#'
        if(filename.Contains("run"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, data));
        else if(filename.Contains("ttbarsignal") && filename.Contains("PlusBbbar"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarsignalPlusBbbar));
        else if(filename.Contains("ttbarsignal") && filename.Contains("PlusB"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarsignalPlusB));
        else if(filename.Contains("ttbarsignal") && filename.Contains("PlusOther"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarsignalPlusOther));
        else if(filename.Contains("ttbarbg"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarbkg));
        else if(filename.Contains("single"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, singletop));
        else if(filename.Contains("ww"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ww));
        else if(filename.Contains("wz"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, wz));
        else if(filename.Contains("zz"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, zz));
        else if(filename.Contains("dyee") && filename.Contains("1050"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, dyee1050));
        else if(filename.Contains("dyee") && filename.Contains("50inf"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, dyee50inf));
        else if(filename.Contains("dymumu") && filename.Contains("1050"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, dymumu1050));
        else if(filename.Contains("dymumu") && filename.Contains("50inf"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, dymumu50inf));
        else if(filename.Contains("dytautau") && filename.Contains("1050"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, dytautau1050));
        else if(filename.Contains("dytautau") && filename.Contains("50inf"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, dytautau50inf));
        else if(filename.Contains("wtolnu"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, wlnu));

        else if(filename.Contains("qcdmu15"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdmu15));
        else if(filename.Contains("qcdmu2030"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdmu2030));
        else if(filename.Contains("qcdmu3050"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdmu3050));
        else if(filename.Contains("qcdmu5080"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdmu5080));
        else if(filename.Contains("qcdmu80120"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdmu80120));
        else if(filename.Contains("qcdmu120170"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdmu120170));
        else if(filename.Contains("qcdem2030"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdem2030));
        else if(filename.Contains("qcdem3080"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdem3080));
        else if(filename.Contains("qcdem80170"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdem80170));
        else if(filename.Contains("qcdbcem2030"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdbcem2030));
        else if(filename.Contains("qcdbcem3080"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdbcem3080));
        else if(filename.Contains("qcdbcem80170"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, qcdbcem80170));

        else if(filename.Contains("ttbarW"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarW));
        else if(filename.Contains("ttbarZ"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarZ));
        else if(filename.Contains("ttbarH") && filename.Contains("inclusiveOther"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarHinclusiveOther));
        // else if(filename.Contains("ttbarH") && filename.Contains("inclusiveBbbar"))
            // v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarHinclusiveBbbar));
        else if(filename.Contains("ttbarH") && filename.Contains("tobbbar"))
            v_filenameSamplePair.push_back(std::pair<TString, Sample>(filename, ttbarHtobbbar));
        else{
            std::cout<<"Sample of FileList ("<<filename<<") is not defined in default samples, skip it...\n";
            continue;
        }
    }

    return v_filenameSamplePair;
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









// ------------------------------------ Further methods of class Samples ---------------------------------------



Samples::Samples(){}



Samples::Samples(const std::vector< Channel::Channel >& v_channel, const std::vector< Systematic::Systematic >& v_systematic)
{
    std::cout<<"--- Beginning to set up the samples\n\n";

    for (auto systematic : v_systematic) {
        for (auto channel : v_channel) {
            this->addSamples(channel, systematic);
        }
    }

    std::cout<<"\n=== Finishing to set up the samples\n\n";
}



void Samples::addSamples(const Channel::Channel& channel, const Systematic::Systematic& systematic)
{
    // Add all samples as they are defined in setSamples()
    std::vector< std::pair< TString, Sample > > v_filenameSamplePair(this->setSamples(channel, systematic));

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



std::vector<Sample> Samples::setSampleOptions(const Systematic::Systematic& systematic, const std::vector< std::pair<TString, Sample> >& v_filenameSamplePair)
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








// ------------------------------------ Functions defined in namespace Tools ---------------------------------------



double Tools::luminosityWeight(const Sample& sample, const double luminosity, RootFileReader* fileReader)
{
    double luminosityWeight=0;
    if(sample.sampleType() == Sample::SampleType::data)return 1;
    const TString& inputFile(sample.inputFile());
    const double crossSection(sample.crossSection());
    if(crossSection<=0.){
        std::cout<<"Sample XSection is <0. Can't calculate luminosity weight!! returning\n";
        return 0;
    }
    const TH1* h_weightedEvents = fileReader->Get<TH1>(inputFile, "weightedEvents");
    const double weightedEvents(h_weightedEvents->GetBinContent(1));
    luminosityWeight = luminosity*crossSection/weightedEvents;
    //std::cout<<"Input file: "<<inputFile<<std::endl;
    //std::cout<<"Luminosity, Xsection, weighted events, lumi weight: "
    //         <<luminosity<<" , "<<crossSection<<" , "<<weightedEvents<<" , "<<luminosityWeight<<std::endl;

    return luminosityWeight;
}











