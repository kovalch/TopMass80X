#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

#include <TString.h>

#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "EventYields.h"
#include "plotterHelpers.h"
#include "PlotterDiffXSSystematic.h"
#include "HistoListReader.h"
#include "higgsUtils.h"
#include "../../common/include/utils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"



void HistoSystematic(const std::vector<std::string>& v_plot, 
                     const std::vector<Channel::Channel>& v_channel,
                     const std::vector<Systematic::Systematic>& v_systematic)
{
    // Set up systematic variations
    std::vector<Systematic::Variation> v_variation;
    v_variation.push_back(Systematic::Variation::up);
    v_variation.push_back(Systematic::Variation::down);
    
    // Setting the list of names of the input root files
    const TString fileList_base = "FileLists_plot_systematic";
    std::map<Channel::Channel, std::map<Systematic::Systematic, std::map<TString, std::pair<TString, TString> > > > m_inputRootFileNames;
    
    for(Channel::Channel channel : v_channel) {
        const TString name_channel = Channel::convert(channel);
        for(Systematic::Systematic systematic : v_systematic) {
            const TString name_systematic = Systematic::convertType(systematic.type());
            for(Systematic::Variation variation : v_variation) {
                TString name_variation = Systematic::convertVariation(variation);
                // Changing the name for systematics that don't have variations
                if(std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), systematic.type()) == Systematic::upDownTypes.end()) {
                    name_variation = "";
                }
                const TString inputFileListName = fileList_base+"/"+"HistoFileList_"+name_systematic+name_variation+"_"+name_channel+".txt";
                std::ifstream file(inputFileListName.Data());
                if(!file) {
                    std::cerr << "### File list not found: " << inputFileListName << " Breaking...\n\n";
                    exit(1);
                }
                // Reading each line of the file corresponding to a separate histogram
                std::string line_;
                while(std::getline(file, line_)) {
                    TString fileName(line_);
                    // Extracting the histogram name
                    TString histoName(fileName);
                    histoName.Replace(0, histoName.Last('/')+1, "");
                    histoName.Resize(histoName.Last('_'));
                    
                    if(variation == Systematic::Variation::up) m_inputRootFileNames[channel][systematic][histoName].first = fileName;
                    else if(variation == Systematic::Variation::down) m_inputRootFileNames[channel][systematic][histoName].second = fileName;
                }
            }
        }
    }
    
    // Create Plotter
    PlotterDiffXSSystematic plotter("Plots", m_inputRootFileNames);
    
    // Access the histoList specifying printing parameters of histograms
    const std::string histoListFile(tth::DATA_PATH_TTH() + "/" + "HistoList_DiffXS_systematic");
    const HistoListReader histoList(histoListFile.data());
    if(histoList.isZombie()){
        std::cerr<<"Error in Histo! Cannot find HistoList with name: "<<histoListFile<<"\n...break\n"<<std::endl;
        exit(12);
    }
    
    // Loop over all histograms in histoList and print them
    std::cout<<"--- Beginning with the plotting\n\n";
    for(auto it = histoList.begin(); it != histoList.end(); ++it){
        
        // Access plot properties from histoList and check whether histogram name contains name pattern
        const PlotProperties& plotProperties = it->second;
        std::cout << "\nchecking " << plotProperties.name << std::endl;
        bool found = false;
        for(const auto& plot : v_plot){
            if(plot.size() && plot[0] == '+'){
                if(plotProperties.name.CompareTo(&plot[1], TString::kIgnoreCase) == 0){
                    found = true;
                    break;
                }
            }
            else if(plotProperties.name.Contains(plot, TString::kIgnoreCase)){
                found = true;
                break;
            }
        }
        if(!found){
            std::cout<<"... no histograms found, continue with next\n";
            continue;
        }
        
        // Set plot properties
        plotter.setOptions(plotProperties.name, plotProperties.specialComment, plotProperties.ytitle, plotProperties.xtitle, 
                           plotProperties.rebin, plotProperties.do_dyscale, plotProperties.logX, plotProperties.logY, 
                           plotProperties.ymin, plotProperties.ymax, plotProperties.xmin, plotProperties.xmax);
        
        // Loop over all systematics and all channels and write histograms
        plotter.producePlots();
    }
    std::cout<<"\n=== Finishing with the plotting\n\n";
}



/// All systematics allowed for plotting
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal, all, allAvailable,
        pu, lept, trig,
        jer, jes,
        btag, 
        btagPt, btagEta,
        btagLjet, 
        btagLjetPt, btagLjetEta,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrPurity,
        kin,
        
        lumi,
        frac_tthf, frac_ttother,
        xsec_tt2b, xsec_ttcc,
        
        topPt,
        pdf
    };
}



int main(int argc, char** argv){
    
    // Get and check configuration parameters
    CLParameter<std::string> opt_plot("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly", false, 1, 100);
    CLParameter<std::string> opt_systematic("s", "Systematic variation - default is allAvailabel. Use 'all' for all, use 'Nominal' for no systematics", false, 1, 100,
        common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up plots
    std::vector<std::string> v_plot {""};
    if(opt_plot.isSet()){
        v_plot = opt_plot.getArguments();
        std::cout<< "Processing only histograms containing in name: ";
        for(auto plot : v_plot)std::cout<< plot << " ";
        std::cout << "\n\n";
    }
    
    // Set up channels
    std::vector<Channel::Channel> v_channel;
    v_channel.push_back(Channel::combined);
    std::cout << "Processing channels: ";
    for(auto channel : v_channel) std::cout << Channel::convert(channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic = Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics);
    if(!opt_systematic.isSet() || opt_systematic[0] == Systematic::convertType(Systematic::allAvailable)){
        v_systematic = common::findSystematicsFromFilelists("FileLists_plot_systematic", v_channel, v_systematic);
        // Adding systematics that do not require specific root files
        for(Systematic::Type type : Systematic::fileIndependentTypes) {
            if(std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), type) != Systematic::upDownTypes.end()) {
                v_systematic.push_back(Systematic::Systematic(type, Systematic::up));
                v_systematic.push_back(Systematic::Systematic(type, Systematic::down));
            } else if(std::find(Systematic::centralTypes.begin(), Systematic::centralTypes.end(), type) != Systematic::centralTypes.end()) {
                v_systematic.push_back(Systematic::Systematic(type, Systematic::central));
            }
        }
    } else{
        if(opt_systematic[0] == Systematic::convertType(Systematic::all)) ; // do nothing
        else v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
    }
    std::cout << "Processing systematics: "; 
    for(auto systematic : v_systematic) std::cout << systematic.name() << " ";
    std::cout << "\n\n";
    
    
    // Start analysis
    HistoSystematic(v_plot, v_channel, v_systematic);
}



