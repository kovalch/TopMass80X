#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>

#include "Sample.h"
#include "higgsUtils.h"
#include "../../common/include/sampleHelpers.h"





Sample::Sample():
legendEntry_(""),
color_(0),
crossSection_(0.),
crossSectionRelativeUp_(0.),
crossSectionRelativeDown_(0.),
sampleType_(dummy),
finalState_(Channel::undefined),
systematic_(),
inputFileName_("")
{}



Sample::Sample(const TString& legendEntry,
               const int color,
               const double& crossSection,
               const double& crossSectionRelativeUp,
               const double& crossSectionRelativeDown,
               const std::vector<TString>& v_filename,
               const SampleType& sampleType):
legendEntry_(legendEntry),
color_(color),
crossSection_(crossSection),
crossSectionRelativeUp_(crossSectionRelativeUp),
crossSectionRelativeDown_(crossSectionRelativeDown),
sampleType_(sampleType),
finalState_(Channel::undefined),
systematic_(),
inputFileName_(""),
v_filename_(v_filename)
{
    if(crossSectionRelativeDown < -0.9) crossSectionRelativeDown_ = crossSectionRelativeUp;
}


void Sample::setLegendEntry(const TString& legendEntry) {legendEntry_ = legendEntry;}

TString Sample::legendEntry()const{return legendEntry_;}



void Sample::setColor(const int color) {color_ = color;}

int Sample::color()const{return color_;}



double Sample::crossSection(const Systematic::Systematic& systematic)const{
    // Check if it is cross section systematic
    const Systematic::Type type = systematic.type();
    if(std::find(Systematic::crossSectionTypes.begin(), Systematic::crossSectionTypes.end(), type) == Systematic::crossSectionTypes.end()) return crossSection_;
    
    // Apply systematic variation for specific samples
    double factor(1.);
    if((type==Systematic::xsec_tt2b && sampleType_==tt2b) ||
       (type==Systematic::xsec_ttcc && sampleType_==ttcc) ||
       (type==Systematic::xsec_ttother && sampleType_==ttother) ||
       (type==Systematic::xsec_ttH && (sampleType_==ttHbb || sampleType_==ttHother) ) ||
       (type==Systematic::xsec_ttZ && sampleType_==ttZ) ){
        if(systematic.variation() == Systematic::up) factor += crossSectionRelativeUp_;
        else if(systematic.variation() == Systematic::down) factor -= crossSectionRelativeDown_;
        else{
            std::cerr<<"ERROR in Sample::crossSection()! Systematic is cross section type, but variation is neither up nor down\n...break\n"<<std::endl;
            exit(813);
        }
    }
    return crossSection_*factor;
}



void Sample::setSampleType(SampleType sampleType) {sampleType_ = sampleType;}

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
        std::cerr<<"ERROR in Sample::setInputFile()! Cannot find requested file: "<<inputFileName<<"\n...break\n"<<std::endl;
        exit(17);
    }
    ifile.close();
    inputFileName_ = inputFileName;
}



TString Sample::inputFile()const{return inputFileName_;}



bool Sample::containsFilenamesOfSample(const Sample& sample, const bool checkReweightedFiles)const
{
    const std::vector<TString>& filenamesToCheck = sample.getFileNames();
    
    for(TString filename : filenamesToCheck) {
        if(!checkReweightedFiles && std::find(v_filename_.begin(), v_filename_.end(), filename) != v_filename_.end()) return true;

        // Contrsucting a filename template representing reweighted sample (should match the naming scheme in load_Analysis.C)
        TString reweightedTemplate(filename);
        reweightedTemplate.ReplaceAll(".root", "_reweighted");
        // Checking whether the sample has any file matching the reweighted and original versions of the filename
        for(TString theFilename : v_filename_) {
            if(theFilename == filename) return true;
            if(theFilename.BeginsWith(reweightedTemplate)) return true;
        }
    }
    
    return false;
}



const std::vector<TString>& Sample::getFileNames()const
{
    return v_filename_;
}







