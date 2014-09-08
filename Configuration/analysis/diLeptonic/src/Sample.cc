#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>

#include "Sample.h"
// #include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"





Sample::Sample():
legendEntry_(""),
color_(0),
crossSection_(0),
sampleType_(dummy),
finalState_(Channel::undefined),
systematic_(),
inputFileName_("")
{}



Sample::Sample(const TString& legendEntry,
               const int color,
               const double& crossSection,
               const std::vector<TString>& v_filename,
               const SampleType& sampleType):
legendEntry_(legendEntry),
color_(color),
crossSection_(crossSection),
sampleType_(sampleType),
finalState_(Channel::undefined),
systematic_(),
inputFileName_(""),
v_filename_(v_filename)
{}



TString Sample::legendEntry()const{return legendEntry_;}



int Sample::color()const{return color_;}



double Sample::crossSection()const{return crossSection_;}



Sample::SampleType Sample::sampleType()const{return sampleType_;}



bool Sample::checkFilename(const TString& filename)const
{
    // Access the real channel in which the sample was processed
    const Channel::Channel& channel = common::finalState(filename);
    
    // Access the basic filename by stripping off the folders
    TString filenameBase(filename);
    if(filenameBase.Contains('/')){
        Ssiz_t last = filenameBase.Last('/');
        filenameBase = filenameBase.Data() + last + 1;
    }
    
    // Strip off the corresponding channel prefix of format "channel_"
    const TString channelPrefix = Channel::convert(channel) + "_";
    filenameBase.ReplaceAll(channelPrefix, "");
    
    // Return whether the filename is associated to this sample
    return std::find(v_filename_.begin(), v_filename_.end(), filenameBase) != v_filename_.end();
}



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







