#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>

#include <vector>
#include <iostream>
#include <set>

#include <future>

#include "plotterclass.h"
#include "HistoListReader.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/sampleHelpers.h"

using namespace std;

///Please put the variation of each systematics one after each other, satarting from the UP variation.
///    NOT valid example: MATCH_UP, MASS_DOWN, MASS_UP
const std::vector<const char*> VectorOfValidSystematics 
    {"Nominal",
    "JER_UP", "JER_DOWN", "JES_UP", "JES_DOWN",
    "PU_UP", "PU_DOWN", "TRIG_UP", "TRIG_DOWN", "LEPT_UP", "LEPT_DOWN",
    "DY_UP", "DY_DOWN", "BG_UP", "BG_DOWN", 
    "KIN_UP", "KIN_DOWN",
    "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
    "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
    "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
//     "BTAG_BEFF_UP", "BTAG_BEFF_DOWN", "BTAG_CEFF_UP", "BTAG_CEFF_DOWN", "BTAG_LEFF_UP", "BTAG_LEFF_DOWN",
    "MASS_UP", "MASS_DOWN", "MATCH_UP", "MATCH_DOWN", "SCALE_UP", "SCALE_DOWN", 
    "POWHEG", "POWHEGHERWIG", "MCATNLO",// "PERUGIA11", // "SPINCORR", 
    "all"};
    
void Histo(bool doControlPlots, bool doUnfold, bool doDiffXSPlotOnly,
           std::vector<std::string> plots, 
           std::vector<std::string> systematics, 
           std::vector<std::string> channels,
           const bool drawUncBand)
{
    //to stay compatible with old code
    std::set<TString> SetOfValidSystematics;
    for (auto s: VectorOfValidSystematics) SetOfValidSystematics.insert(s);

    HistoListReader histoList(doControlPlots ? "HistoList_control" : "HistoList");
    if (histoList.IsZombie()) exit(12);
    for (auto it = histoList.begin(); it != histoList.end(); ++it) {
        const PlotProperties& p = it->second;
        std::cout << "checking " << p.name << std::endl;
        bool found = false;
        for (auto plot : plots) {
            if (plot.size() && plot[0] == '+') {
                if (p.name.CompareTo(&plot[1], TString::kIgnoreCase) == 0) {
                    found = true;
                    break;
                }
            } else if (p.name.Contains(plot, TString::kIgnoreCase)) {
                found = true;
                break;
            }
        }
        if (!found) continue;

        // Create Plotter 
        Plotter h_generalPlot;
        h_generalPlot.ListOfSystematics(SetOfValidSystematics);
        
        /////////////////////////////////////////////////////
        /////////   UNFOLDING OPTIONS     ///////////////////
        /////////////////////////////////////////////////////

        TString outpath = "";
        h_generalPlot.UnfoldingOptions(doUnfold);
        h_generalPlot.SetOutpath("");

        /// Do fit in the ratio plot, if ratio plot done
        h_generalPlot.DoFitInRatio(0);

        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        ///////////////////////////////////////////////////// 
        h_generalPlot.setOptions(p.name,p.specialComment,p.ytitle,p.xtitle, 
                                    p.rebin, p.do_dyscale, p.logX, p.logY, 
                                    p.ymin, p.ymax, p.xmin, p.xmax, p.bins, p.xbinbounds, p.bincenters);
        h_generalPlot.DYScaleFactor(p.specialComment);
        //need preunfolding for ALL channels before unfolding!!
        for (auto channel : channels) {
            for (auto systematic : systematics) {
                h_generalPlot.setDrawUncBand(drawUncBand);
                h_generalPlot.preunfolding(channel, systematic);
            }
        }
        if (doUnfold) {
            for (auto systematic:systematics){
                //unfold all channels except combined in parallel
                std::vector<std::future<void>> unfoldJobs;
                for (auto channel:channels) {
/// This part of the code was only necessary for the statistical combination of the 3 channels
//                     if (channel != "combined") {
                        TString ch(channel.c_str());
                        TString sys(systematic.c_str());
                        unfoldJobs.push_back(std::async(std::launch::async, [ch, sys](Plotter p) -> void { 
                            p.unfolding(ch,sys);
                        }, h_generalPlot));
                        //FIXME:
                        //it seems unfolding is not thread safe! 
                        //fix that and remove the next line
                        unfoldJobs.at(unfoldJobs.size()-1).wait();
//                     }
                }
                //wait for the 3 channels to finish
                for (auto &i : unfoldJobs) {
                    i.get();
                }
/// This part of the code was only necessary for the statistical combination of the 3 channels
//                 //only now do combined (if requested)
//                 for (auto channel:channels) {
//                     if (channel == "combined") h_generalPlot.unfolding(channel,systematic);
//                 }
            }
        std::cout << "Done with the unfolding\n";
        }
        if (doDiffXSPlotOnly) {
            for (auto channel:channels){
                std::vector<TString> syst (systematics.begin(), systematics.end());
                std::vector<TString> valid_sys;
                for (size_t sys=0; sys<syst.size(); ++sys){
                    ///Ugly method to use the systematic convention used up to now
                    if(syst.at(sys) == TString("Nominal") || syst.at(sys).Contains("_DOWN") ||
                      syst.at(sys) == "MCATNLO" || syst.at(sys).Contains("HERWIG") || syst.at(sys).Contains("PERUGIA11") || syst.at(sys).Contains("SPINCORR") ){continue;}

                    if (!syst.at(sys).Contains("POWHEG")) {valid_sys.push_back(syst.at(sys).Remove(syst.at(sys).Length()-2, 2));}
                    else {valid_sys.push_back(syst.at(sys));};
                };
                h_generalPlot.PlotDiffXSec(channel, valid_sys);
            }
        }
    }
}

/**
 * Helper function to create a function which checks if a string found is in the
 * passed vector of string.
 * 
 * @param allowed a vector of allowed strings (char*s)
 * @return a function taking a std::string and returning a bool
 */
std::function<bool(const std::string &s)> makeStringChecker(const std::vector<const char *> allowed) {
    return [allowed](const std::string &test) {
        return std::find(begin(allowed), end(allowed), test) != end(allowed);
    };
}

int main(int argc, char** argv) {
    CLParameter<std::string> opt_type("t", "cp|unfold|plot - required, cp=contol plots, unfold, or only plot diffXS", true, 1, 1,
        makeStringChecker({"cp", "unfold", "plot"}));
    CLParameter<std::string> opt_plots("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly", false, 1, 100);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
        makeStringChecker({"ee", "emu", "mumu", "combined"}));
    CLParameter<std::string> opt_sys("s", "Systematic variation - default is Nominal, use 'all' for all", false, 1, 100,
        makeStringChecker(VectorOfValidSystematics));
    CLParameter<bool> opt_band("b", "If existing, draw uncertainty band if '-t cp' ", false, 0, 0);
    CLAnalyser::interpretGlobal(argc, argv);
    
    std::vector<std::string> channels { "emu", "ee", "mumu", "combined" };
    if (opt_channel.isSet()) channels = opt_channel.getArguments();
    std::cout << "Processing channels: "; 
    for (auto ch: channels) std::cout << ch << " "; std::cout << "\n";
        
    std::vector<std::string> systematics { "Nominal" };
    if (opt_sys.isSet()) {
        systematics = opt_sys.getArguments();
        if (systematics[0] == "all") {
            systematics.clear();
            for (string syst: VectorOfValidSystematics) {
                if (syst != "all" && syst != "POWHEGHERWIG") systematics.push_back(syst);
            }
        }
    }
    std::cout << "Processing systematics (use >>-s all<< to process all knwon systematics): "; 
    for (auto sys: systematics) std::cout << sys << " "; std::cout << "\n";

    std::vector<std::string> plots { "" };
    if (opt_plots.isSet()) plots = opt_plots.getArguments();

    bool doControlPlots = opt_type[0] == "cp";
    bool doDiffXSPlotOnly = opt_type[0] == "plot";
    bool doUnfold = opt_type[0] == "unfold";
    bool drawUncBand = (opt_band.isSet() && doControlPlots);
    Histo(doControlPlots, doUnfold, doDiffXSPlotOnly, plots, systematics, channels, drawUncBand);
}
