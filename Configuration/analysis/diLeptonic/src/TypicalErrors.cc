#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <TString.h>

#include "../../common/include/CommandLineParameters.h"

using namespace std;


std::vector<TString> Channels {"ee", "emu", "mumu", "combined"};

std::vector<TString> Variables_FullPS {"ToppT", "TopRapidity","ToppTTTRestFrame", "ToppTLead", "ToppTNLead",
                                       "TTBarpT", "TTBarRapidity", "TTBarMass", "TTBarDeltaPhi"
                                      };

std::vector<TString> Variables_VisiblePS {"BJetpT","BJetEta","LeptonpT","LeptonEta",
                                          "LLBarpT", "LLBarMass","LeptonBjetMass", "BBBarMass", "BBBarpT"
                                          };

std::vector<TString> Systematics {"TRIG_", "LEPT_","BG_", "DY_", "JES_", "JER_", "PU_",
                                  "BTAG_", "BTAG_LJET_", "BTAG_PT_", "BTAG_LJET_PT_",
                                  "KIN_",
                                  "HAD_", "MASS_", "SCALE_", "MATCH_", "PDF_"
                                  };


std::vector<TString> Files(TString channel = "", TString variable = "", TString phaseSpace = "")
{
    std::vector<TString> WhichVariable;
    std::vector<TString> FileVector;

    if (variable!=""){WhichVariable.push_back(variable);}
    else if (phaseSpace == "full")    {WhichVariable = Variables_FullPS;}
    else if (phaseSpace == "visible") {WhichVariable = Variables_VisiblePS;}
    else{
        WhichVariable = Variables_FullPS;
        WhichVariable.reserve(WhichVariable.size()+Variables_VisiblePS.size());
        WhichVariable.insert(WhichVariable.end(), Variables_VisiblePS.begin(), Variables_VisiblePS.end());
    }

    for (size_t j=0; j<WhichVariable.size(); j++){
        FileVector.push_back(TString("Plots/FinalResults/").Append(channel).Append("/").Append(WhichVariable.at(j)).Append("_SystematicsLaTeX.txt"));
    }

    return FileVector;
}


std::vector<string> SplitLine(string Line)
{
    /// Returns a std::vector which its elements will be the content of 'Line' separated by blank space ' '

    std::vector<string> output_vector;
    istringstream iss(Line);
    copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter<std::vector<string> >(output_vector));
    return output_vector;
}




std::vector<double> ReadLineFromFile (TString Filename, TString Systematic)
{
    /// Returns typical error for systematic 'Systematic' in file 'Filename'

    if (Filename == "" || Systematic == ""){
        std::cout<<"\n\n******** ERROR ******** ERROR ******** ERROR ******** ERROR ********"<<std::endl;
        std::cout<<"You didn't provide a file path neither a systematic. Exiting!"<<std::endl;
        std::cout<<"\n\n******** ERROR ******** ERROR ******** ERROR ******** ERROR ********"<<std::endl;
        exit(9);
    }

    string LineToSplit;
    ifstream infile;
    std::vector<string> SplittedLine;
    std::vector<double> errors;

    infile.open(Filename, ios_base::in);
    if (infile.fail()) {
        std::cout<<"\n\n******** WARNING ******** WARNING ******** WARNING ******** WARNING ********"<<std::endl;
        std::cout<<"The file "<<Filename<<" you requested doesn't exist."<<std::endl;
        std::cout<<"Tis file will be skiped in the calculation of 'typical error'"<<std::endl;
        std::cout<<"******** WARNING ******** WARNING ******** WARNING ******** WARNING ********\n\n"<<std::endl;
        exit(23);
    }
    while (!infile.eof()) {
        LineToSplit.clear();
        getline(infile, LineToSplit);
        SplittedLine.clear();
        SplittedLine = SplitLine(LineToSplit);
        if(SplittedLine.size()<=0){break;}
        if (SplittedLine.at(0) != Systematic){continue;}
        for(size_t i=1; i<SplittedLine.size(); i++){
            if(SplittedLine.at(i) == "Lin.Avg.(%)="){break;}
            errors.push_back(atof(SplittedLine.at(i).c_str()));
        }
    }
    infile.close();
    return errors;
}



void SanityCheck( TString channel = "", TString systematic = "", TString variable = "")
{
    std::vector<TString> ValidVariable = Variables_FullPS;
    ValidVariable.reserve(ValidVariable.size()+Variables_VisiblePS.size());
    ValidVariable.insert(ValidVariable.end(), Variables_VisiblePS.begin(), Variables_VisiblePS.end());
    
    if (find(Channels.begin(), Channels.end(), channel) == Channels.end() && channel != ""){
        std::cout<<"\n\nThe proposed channel '"<<channel<<"' is not valid. Exiting!\n"<<std::endl;
        exit(2);
    }

    if (find(Systematics.begin(), Systematics.end(), systematic) == Systematics.end() && systematic != ""){
        std::cout<<"\n\nThe proposed systematic '"<<systematic<<"' is not valid (or is not implemented yet). Exiting!\n"<<std::endl;
        exit(22);
    }

    if (find(ValidVariable.begin(), ValidVariable.end(), variable) == ValidVariable.end() && variable != ""){
        std::cout<<"\n\nThe proposed variable '"<<variable<<"' is not valid (or is not implemented yet). Exiting!\n"<<std::endl;
        exit(222);
    }

}

void TypicalError( TString channel = "", TString systematic = "", TString variable = "", TString phaseSpace = "")
{
    SanityCheck(channel, systematic, variable);
    
    std::vector<TString> Channel, Systematic;

    if ( channel != ""){Channel.push_back(channel);}
    else { Channel = Channels; }
    
    if ( systematic != "") {Systematic.push_back(systematic);}
    else { Systematic = Systematics;}
    
    for (size_t l=0; l<Channel.size(); l++){
        std::vector<TString> FileList = Files(Channel.at(l), variable, phaseSpace);
        std::vector<double> error;
        std::cout<<"----------------------------------------"<<std::endl;

        for (size_t j=0; j<Systematic.size(); j++){
            error.clear();
            for (size_t i=0; i<FileList.size(); i++){
                std::vector<double> Typ_error = ReadLineFromFile(FileList.at(i), Systematic.at(j));
                if(Typ_error.size()) error.insert(error.end(), Typ_error.begin(), Typ_error.end());
            }
            std::sort(error.begin(), error.end());

            int extra = (error.size()%2) ? 0 : 1;
            int meanPoint =error.size()/2;
            printf("Total typical error for systematic %*s in channel %s is: %2.2f %%\n", 15, Systematic.at(j).Data(), Channel.at(l).Data(), 0.5*(error.at(meanPoint-extra) + error.at(meanPoint)));
        }
    }
}


int main(int argc, char** argv) {
    CLParameter<std::string> opt_v("v", "Return the typical error for certain variable, e.g. 'ToppTLead', 'LLBarMass', ...", false, 1, 1);
    CLParameter<std::string> opt_s("s", "Return the typical systematic uncertainty for a certain systematic variation, e.g. 'PU_', 'TRIG_', 'BTAG_LJET_ETA_', 'BTAG_PT_', ...", false, 1, 1);
    CLParameter<std::string> opt_c("c", "Return the typical systematic uncertainty for an specific channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            [](const std::string &ch){return ch == "" || ch == "ee" || ch == "emu" || ch == "mumu" || ch == "combined";});
    CLParameter<std::string> opt_ps("ps", "Specify set of variables according to the phase space in which they are measured. Valid options: visible, full", false, 1, 1,
            [](const std::string &ps){return ps == "full" || ps == "visible";});
    CLAnalyser::interpretGlobal(argc, argv);
    
    TString ValidSystematics = opt_s.isSet() ? opt_s[0] : "";
    TString ValidVariable    = opt_v.isSet() ? opt_v[0] : "";
    TString ValidChannel     = opt_c.isSet() ? opt_c[0] : "";
    TString ValidPhaseSpace  = (!opt_v.isSet() && opt_ps.isSet()) ? opt_ps[0] : "";
        
    TypicalError(ValidChannel, ValidSystematics, ValidVariable, ValidPhaseSpace);
    
    return 0;
}
