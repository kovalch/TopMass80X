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
#include <map>
#include <TH1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TSystem.h>

#include "../../common/include/CommandLineParameters.h"

const std::vector<std::string> valid_systs {
    "TRIG_UP", "TRIG_DOWN", 
    "LEPT_UP", "LEPT_DOWN", 
    "BG_UP", "BG_DOWN", 
    "DY_UP", "DY_DOWN", 
    "JES_UP", "JES_DOWN", 
    "JER_UP", "JER_DOWN", 
    "PU_UP","PU_DOWN",
    "BTAG_UP", "BTAG_DOWN", 
    "BTAG_LJET_UP","BTAG_LJET_DOWN",
    "BTAG_PT_UP", "BTAG_PT_DOWN", 
    "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", 
    "KIN_UP", "KIN_DOWN", 
    "POWHEG", "MCATNLO", 
    "MASS_UP", "MASS_DOWN", 
    "SCALE_UP", "SCALE_DOWN", 
    "MATCH_UP","MATCH_DOWN"
};




/// Calculate the relative difference, with respect to the nominal, of a given systematic uncertainty
std::vector<double> CalcRelDiff(std::string nominal, std::string variation)
{
    std::ifstream infile1, infile2;
    std::vector<double> result(10, 1);

    infile1.open(nominal, std::ios_base::in);
    if (infile1.fail()) {
        std::cout<<"\n\n******** WARNING ******** WARNING ******** WARNING ******** WARNING ********"<<std::endl;
        std::cout<<"The file "<<nominal<<" you requested doesn't exist."<<std::endl;
        std::cout<<"Tis file will be skiped in the calculation of 'typical error'"<<std::endl;
        std::cout<<"******** WARNING ******** WARNING ******** WARNING ******** WARNING ********\n\n"<<std::endl;
        return result;
    }
    infile2.open(variation, std::ios_base::in);
    if (infile2.fail()) {
        infile1.close();
        std::cout<<"\n\n******** WARNING ******** WARNING ******** WARNING ******** WARNING ********"<<std::endl;
        std::cout<<"The file "<<variation<<" you requested doesn't exist."<<std::endl;
        std::cout<<"Tis file will be skiped in the calculation of 'typical error'"<<std::endl;
        std::cout<<"******** WARNING ******** WARNING ******** WARNING ******** WARNING ********\n\n"<<std::endl;
        return result;
    }

    std::string dummy;
    double tmplowBin = 0, tmphighBin = 0, tmpbinCtr = 0, tmpXsec = 0, tmpStat = 0;
    std::vector<double> nomxsec, nomstatError, nomgenvalue, nomlowBin, nomhighBin, nombinCenter;
    while (!infile1.eof()) {
        infile1>>dummy>>tmplowBin>>dummy>>tmphighBin>>dummy>>tmpbinCtr>>dummy>>tmpXsec>>dummy>>tmpStat>>dummy>>tmpbinCtr;
        nomxsec.push_back(tmpXsec);
        nomlowBin.push_back(tmplowBin);
        nomhighBin.push_back(tmphighBin);
        nomstatError.push_back(tmpStat);
    }
    infile1.close();

    std::vector<double> varxsec, varstatError, vargenvalue, varlowBin, varhighBin, varbinCenter;
    while (!infile2.eof()) {
        infile2>>dummy>>tmplowBin>>dummy>>tmphighBin>>dummy>>tmpbinCtr>>dummy>>tmpXsec>>dummy>>tmpStat>>dummy>>tmpbinCtr;
        varxsec.push_back(tmpXsec);
        varlowBin.push_back(tmplowBin);
        varhighBin.push_back(tmphighBin);
        varstatError.push_back(tmpStat);
    }
    infile2.close();

    result.resize(nomxsec.size());
    for(size_t iter = 0; iter<result.size(); iter++)
    {
        double diff = varxsec.at(iter) - nomxsec.at(iter);
        result.at(iter) = 100*diff/nomxsec.at(iter);
    }

    return result;
}

/// Get the bin boundaries for a given distribution
std::vector<double> getBinBoundaries(std::string file)
{
    std::ifstream infile1;
    infile1.open(file, std::ios_base::in);
    std::vector<double> nomlowBin, nomhighBin;
    if (infile1.fail()) {
        std::cout<<"\n\n******** WARNING ******** WARNING ******** WARNING ******** WARNING ********"<<std::endl;
        std::cout<<"The file "<<file<<" you requested doesn't exist."<<std::endl;
        std::cout<<"******** WARNING ******** WARNING ******** WARNING ******** WARNING ********\n\n"<<std::endl;
        return nomlowBin;
    }

    std::string dummy;
    double tmplowBin,tmphighBin, dummyDouble;
    while (!infile1.eof()) {
        infile1>>dummy>>dummyDouble>>dummy>>tmplowBin>>dummy>>tmphighBin>>dummy>>dummyDouble>>dummy>>dummyDouble>>dummy>>dummyDouble;
        nomlowBin.push_back(tmplowBin);
        nomhighBin.push_back(tmphighBin);
    }
    infile1.close();
    nomlowBin.push_back(nomhighBin.back());

    return nomlowBin;
}

/// Plot the relative difference of different systematic variations
void PlotDifferenceToCentralResult(std::string channel, std::string variable)
{
    const std::string outDir = "SystComparison/";
    const std::string appendix = "/"+channel+"/Hyp"+variable+"Results.txt";
    std::string filenominal = "UnfoldingResults/Nominal/"+appendix;

    std::vector<std::string> systematic = valid_systs;
    const int nSyst = systematic.size();
    std::vector<double> binRanges = getBinBoundaries(filenominal);
    const int nbins = binRanges.size()-1;
    double bins[nbins+1];
    std::copy(binRanges.begin(), binRanges.end(), bins);

    std::string xtitle = "";

    if(!variable.find("Top") )        {xtitle = xtitle+"^{t}";}
    else if(!variable.find("TTBar"))  {xtitle = xtitle+"^{t#bar{t}}";}
    else if(!variable.find("Lepton")) {xtitle = xtitle+"^{l}";}
    else if(!variable.find("LLBar"))  {xtitle = xtitle+"^{l#bar{l}}";}
    else if(!variable.find("BJet"))   {xtitle = xtitle+"^{b}";}
    else if(!variable.find("BBBar"))  {xtitle = xtitle+"^{b#bar{b}}";}

    if(variable.find("pT")  && variable.find("pT") != std::string::npos)                  {xtitle = "p_{T}"+xtitle;}
    else if(variable.find("Eta")  && variable.find("Eta") != std::string::npos)           {xtitle = "#eta"+xtitle;}
    else if(variable.find("Rapidity")  && variable.find("Rapidity") != std::string::npos) {xtitle = "y"+xtitle;}
    else if(variable.find("Mass")  && variable.find("Mass") != std::string::npos)         {xtitle = "m"+xtitle;}

    TLegend *leg = new TLegend(0.1, 0.7, 0.9, 0.9);
    leg->SetNColumns(8);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    TH1F *histos[nSyst];
    int color = 0;
    for (int iter = 0; iter<nSyst; iter++)
    {
        std::string syst = systematic.at(iter);
        std::string fileup = "UnfoldingResults/"+syst+appendix;

        std::vector<double> results;
        results.clear();
        results = CalcRelDiff(filenominal, fileup);
        histos[iter] = new TH1F(syst.c_str(), syst.c_str(), nbins, bins);
        histos[iter]->SetTitle("Relative uncertainty (in %)");
        histos[iter]->GetYaxis()->SetRangeUser(-10,10);
        histos[iter]->GetYaxis()->SetTitle("Unc. (%)");
        histos[iter]->GetXaxis()->SetTitle(xtitle.c_str());
        histos[iter]->GetXaxis()->SetTitleOffset(1.1);
        if(iter==10) color +=10;
        if(iter%2==0)color++;
        histos[iter]->SetLineColor(color);
        if(syst.find("UP") && syst.find("UP") != std::string::npos){
            histos[iter]->SetLineStyle(1);
        } else if(syst.find("DOWN")){
            histos[iter]->SetLineStyle(3);
        }
        leg->AddEntry(histos[iter], syst.c_str(), "l");
        for(int binIter=0; binIter<nbins; binIter++){
            histos[iter]->SetBinContent(binIter+1, results.at(binIter));
        }
    }

    TCanvas *c1 = new TCanvas("c", "c", 1000, 1000);
    gStyle->SetOptStat(0);
    histos[0]->Draw();
    for (int iter = 1; iter< nSyst; iter++){histos[iter]->Draw("same");}
    leg->Draw("same");
    c1->Update();
    TLine *line = new TLine(c1->GetUxmin(), 0., c1->GetUxmax(), 0.);
    line->SetLineColor(kBlack);
    line->Draw("same");

    gSystem->mkdir(outDir.c_str());
    c1->Print((outDir+"systematicComparison_"+channel+"_"+variable+".eps").c_str());

    c1->Close(); delete c1;
    leg->Clear(); delete leg;
    for (int iter = 0 ; iter<nSyst; iter++){
        delete histos[iter];
    }
}




int main(int argc, char** argv) {
    CLParameter<std::string> opt_var("v", "Return the typical error for certain variable, e.g. 'ToppTLead', 'LLBarMass', ...", true, 1, 100);
    CLParameter<std::string> opt_chan("c", "Return the typical systematic uncertainty for an specific channel: ee, emu, mumu, combined", true, 1, 4,
            [](const std::string &ch){return ch == "ee" || ch == "emu" || ch == "mumu" || ch == "combined";});
    CLAnalyser::interpretGlobal(argc, argv);

    std::vector<std::string> channels, variables;
    if (opt_var.isSet())  {variables = opt_var.getArguments();}
    if (opt_chan.isSet()) {channels = opt_chan.getArguments();}

    for (auto chan : channels){
        for (auto var : variables){
            PlotDifferenceToCentralResult(chan, var);
        }
    }
    return 0;
}
