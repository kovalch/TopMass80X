/**
*   Macro to obtain the theory predictions for all parton level distributions from the dileptonic analysis
*   
*   Input:
*       $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic/selectionRoot/<theory>/<channel>/<channel>_ttbarsignalplustau.root
*
*   within this ROOT file the histograms "VisGen*" must exist.
*   Usually these files are produced after analysing the NTuples: install/bin/load_Analysis
*
*   Output:
*       $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic/theoryHistograms/<channel>/VisGen<variable>.root
*
*   Execute via:
*       cd $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic/
*       root -l -b -q macros/obtainTheoryCurves.C
*/



#include <iostream>
#include <fstream>
#include <cstdio>

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"

std::string basedir = "./analysis_postCWR/selectionRoot/";
std::string outputDir = "./theoryHistograms/";

const size_t nmodels = 4, nparticles = 9, nchannels = 3;
std::string models [4] = {"Nominal", "MCATNLO", "POWHEG", "POWHEGHERWIG"};
std::string particles [9] = {"VisGenToppT", "VisGenToppTTTRestFrame", "VisGenToppTLead", "VisGenToppTNLead", "VisGenTopRapidity",
                            "VisGenTTBarDeltaPhi", "VisGenTTBarpT", "VisGenTTBarRapidity", "VisGenTTBarMass"}; 
std::string channels [3] = {"ee", "emu", "mumu"};



/// Check existence file. Return 1 is exists, 0 if not. 
bool checkFileExistence(std::string file)
{
    std::ifstream f(file.c_str());
    if (f.good())
    {
        f.close();
        return true;
    } else {
        std::cout<<"File '"<<file<<"'\nDoes not exist.\nEXIT!"<<std::endl;
        f.close();
        return false;
    }
}


void getCurves(std::string particle)
{
    std::cout<<"Getting theory curves for '"<<particle<<"'"<<std::endl;
    
    /// Iterate in all channels
    for (size_t channel = 0; channel < nchannels; channel++)
    {
        gSystem->mkdir((outputDir+channels[channel]).c_str(), kTRUE);
        TFile* foutput = new TFile((outputDir+channels[channel]+"/"+particle+".root").c_str(), "RECREATE");
        foutput->Close();

        /// Iterate in all theory models
        for (size_t model = 0; model<nmodels; model++)
        {
            std::string mmodel = "";
            if(models[model] == "POWHEG")       mmodel = "_powheg";
            if(models[model] == "POWHEGHERWIG") mmodel = "_powhegHerwig";
            if(models[model] == "MCATNLO")      mmodel = "_mcatnlo";
        
            std::string inputFile = basedir+models[model]+"/"+channels[channel]+"/"+channels[channel]+"_ttbarsignalplustau"+mmodel+".root";
            if(!checkFileExistence(inputFile)) exit(10);
    
            TFile *f_input = new TFile(inputFile.c_str());
            if(!f_input)
            {
                std::cout<<"File '"<<inputFile<<"'\nCannot be opened.\nEXIT!"<<std::endl;
                exit(11);
            }
            
            TH1D *histogram = (TH1D *) f_input->Get(particle.c_str())->Clone("histogram");
            if(!histogram)
            {
                std::cout<<"Histogram '"<<particle<<"'\ncannot be obtained from file '"<<inputFile<<"'\nEXIT!"<<std::endl;
                exit(12);
            }
            if(particle.find("Top") != std::string::npos && particle.find("Lead") == std::string::npos)
            {
                std::string antiparticle = particle;
                antiparticle = antiparticle.replace(antiparticle.begin(), antiparticle.begin()+9, "VisGenAntiTop");
                histogram->Add((TH1D *) f_input->Get(antiparticle.c_str()));
            }

            std::string modelName = (models[model] == "Nominal"       ? "Madgraph" :
                                    (models[model] == "POWHEG"        ? "Powheg" :
                                    (models[model] == "POWHEGHERWIG"  ? "PowhegHerwig" :
                                    (models[model] == "MCATNLO"       ? "Mcatnlo" : "unknown")
                                    )));
            foutput = new TFile((outputDir+channels[channel]+"/"+particle+".root").c_str(), "UPDATE");
            histogram->Write(modelName.c_str(), TObject::kOverwrite);
            foutput->Close();

            delete histogram;
            f_input->Close();
            delete f_input;
        }
        
    }
}


void getCurvesCombined(std::string particle)
{

    gSystem->mkdir((outputDir+"combined").c_str(), kTRUE);
    TFile *foutput = new TFile((outputDir+"combined/"+particle+".root").c_str(), "RECREATE");
    foutput->Close();
    std::string mmodel = "";

    /// Iterate in all channels
    for (size_t model = 0; model<nmodels; model++)
    {
        if(models[model] == "Nominal")      mmodel = "Madgraph";
        if(models[model] == "POWHEG")       mmodel = "Powheg";
        if(models[model] == "POWHEGHERWIG") mmodel = "PowhegHerwig";
        if(models[model] == "MCATNLO")      mmodel = "Mcatnlo";

        std::string ee_filename = outputDir+"ee/"+particle+".root";     if(!checkFileExistence(ee_filename)) exit(13);
        std::string emu_filename = outputDir+"emu/"+particle+".root";   if(!checkFileExistence(emu_filename)) exit(14);
        std::string mumu_filename = outputDir+"mumu/"+particle+".root"; if(!checkFileExistence(mumu_filename)) exit(15);

        TFile *ee_file   = TFile::Open(ee_filename.c_str());
        TFile *emu_file  = TFile::Open(emu_filename.c_str());
        TFile *mumu_file = TFile::Open(mumu_filename.c_str());

        TH1D *ee_histo   = (TH1D *) ee_file->Get(mmodel.c_str())->Clone("ee");
        TH1D *emu_histo  = (TH1D *) emu_file->Get(mmodel.c_str())->Clone("emu");
        TH1D *mumu_histo = (TH1D *) mumu_file->Get(mmodel.c_str())->Clone("mumu");

        TH1D *histogram = (TH1D *)ee_histo->Clone("histogram");
        histogram->Add(emu_histo);
        histogram->Add(mumu_histo);

        foutput = new TFile((outputDir+"combined/"+particle+".root").c_str(), "UPDATE");
        if(!foutput) exit(16);
        histogram->Write(mmodel.c_str(), TObject::kOverwrite);
        foutput->Close();
    }
}




void obtainTheoryCurves()
{
    std::cout<<"\n\n\n\033[1;31m******************************************************************\n"
             <<"******************************************************************\033[1;m\n\n"
             <<"\033[1;34m Macro to obtain the theory predictions for all parton level distributions from the dileptonic analysis\033[1;m\n"
             <<"\033[1;34m Input:\033[1;m\n"
             <<"     $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic/selectionRoot/<theory>/<channel>/<channel>_ttbarsignalplustau.root\n"
             <<" within this ROOT file the histograms 'VisGen*' must exist.\n"
             <<" Usually these files are produced after analysing the NTuples: install/bin/load_Analysis\n"
             <<"\033[1;34m Output\033[1;m :\n"
             <<"    $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic/theoryHistograms/<channel>/VisGen<variable>.root\n"
             <<" Execute via:\n"
             <<"     cd $CMSSW_BASE/src/TopAnalysis/Configuration/analysis/diLeptonic/\n"
             <<"     root -l -b -q macros/obtainTheoryCurves.C++g\n\n"
             <<"\033[1;31m******************************************************************\n"
             <<"******************************************************************\033[1;m\n\n\n\n\n"<<std::endl;

    for (size_t particle = 0; particle < nparticles; particle++)
    {
        getCurves(particles[particle]);
        getCurvesCombined(particles[particle]);
    }

}