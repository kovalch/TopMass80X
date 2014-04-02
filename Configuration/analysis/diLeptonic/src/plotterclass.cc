#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <memory>

#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TExec.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TList.h>
#include <TGraphAsymmErrors.h>
#include <TError.h>

#include "plotterclass.h"
#include "ttbarUtils.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/plotterUtils.h"
#include "../../common/include/ScaleFactors.h"


// DAVID
#include "DilepSVDFunctions.h"

using namespace std;

const bool Plotter::doClosureTest = 0; //Signal is MC - so add BG to it and dont do DY scaling

void Plotter::setLumi(double newLumi, double xSec)
{
    this->lumi = newLumi;
    this->topxsec = xSec;
}

void Plotter::setDrawUncBand(bool drawUncBand)
{
    drawUncBand_ = drawUncBand;
}

// DAVID
void Plotter::UnfoldingOptions(bool doSVD)
{
  doUnfolding = doSVD;
  drawNLOCurves = true; // boolean to draw/not-draw extra theory curves in the Diff.XSection plots


  drawPlotRatio    = true;
  drawSmoothMadgraph = false;
  drawMCATNLO      = true;
  drawKidonakis    = true;
  drawAhrens       = true;
  drawPOWHEG       = true;;
  drawPOWHEGHERWIG = false;
  drawPERUGIA11 = false;
  drawMadScaleMatching = false;
  drawMadMass     = false;
}

void Plotter::DoFitInRatio(bool doFit)
{
    doFit_ = doFit;
}

// DAVID
void Plotter::SetOutpath(TString path)
{
  outpath = path;
}


void Plotter::unfolding(TString channel, TString systematic)
{

    std::cout << "\n\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Starting Calculation of Diff. X Section for '" << name << "' in Channel '" << channel << "' and Systematic '"<< systematic <<"':" << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;

    if(channel == "combined")
    {
        ifstream ResultsEE ("UnfoldingResults/"+systematic+"/ee/"+name+"Results.txt");
        ifstream ResultsEMu("UnfoldingResults/"+systematic+"/emu/"+name+"Results.txt");
        ifstream ResultsMuMu("UnfoldingResults/"+systematic+"/mumu/"+name+"Results.txt");

        if (!ResultsEE.is_open()){
            Plotter::preunfolding("ee", systematic);
            CalcDiffXSec("ee",systematic);
        } else {
            ResultsEE.close();
        }
        if (!ResultsEMu.is_open()){
            Plotter::preunfolding("emu", systematic);
            CalcDiffXSec("emu",systematic);
        } else {
            ResultsEMu.close();
        }
        if (!ResultsMuMu.is_open()){
            Plotter::preunfolding("mumu", systematic);
            CalcDiffXSec("mumu",systematic);
        } else {
            ResultsMuMu.close();
        }
    }

    CalcDiffXSec(channel,systematic);

    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Finished Calculation of Diff. X Section for '" << name << "' in Channel '" << channel << "' and Systematic '"<< systematic <<"':" << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;

}

void Plotter::preunfolding(TString Channel, TString Systematic)
{
    write(Channel, Systematic);
}


void Plotter::DYScaleFactor(TString SpecialComment){

    DYScale = {1,1,1,1};

    if(!doDYScale || doClosureTest) return; //need to make a switch for control plots that don't want DYScale

    TString nameAppendix = "";
    if ( !SpecialComment.BeginsWith("_post") &&  SpecialComment != "Standard" ){
        std::cout<<"\n\n*******************************************************************"<<std::endl;
        std::cout<<"ERROR: When calculating the DY Scale factor you must specify in which step you want to calculate the DY SF:"<<std::endl;
        std::cout<<" '_postZcut', '_post2jets', '_postMET', '_post1btag', '_postKinReco' or 'Standard' = _postKinReco"<<std::endl;
        std::cout<<"*******************************************************************\n\n"<<std::endl;
        exit(444);
    }
    if (SpecialComment.BeginsWith("_post")){
        if(SpecialComment.EqualTo("_postZcut"))nameAppendix = "_step4";
        if(SpecialComment.EqualTo("_post2jets"))nameAppendix = "_step5";
        if(SpecialComment.EqualTo("_postMET"))nameAppendix = "_step6";
        if(SpecialComment.EqualTo("_post1btag"))nameAppendix = "_step7";
        if(SpecialComment.EqualTo("_postKinReco"))nameAppendix = "_step8";
        
    } else if ( SpecialComment == "Standard") {
        nameAppendix = "_step8";
    }

    std::cout<<"\n\nBegin DYSCALE FACTOR calculation at selection step "<<nameAppendix<<std::endl;
    
    std::vector<TString> Vec_Files = InputFileList("combined", "Nominal");//Read the hardcoded list of files
    if(Vec_Files.size()<1) {std::cout<<"WARNING(in DYScaleFactor)!!! No datasets available to calculate DY SF. EXITING!!"<<std::endl; return;}
    
    double NoutEEDYMC=0, NinEEDYMC=0, NoutMuMuDYMC=0, NinMuMuDYMC=0;//Number of events in/out of z-veto region for the DY MC
    double NinEE=0, NinMuMu=0, NinEMu=0;//Number of events in z-veto region for data
    double NinEEloose=0, NinMuMuloose=0;//Number of data events in Z-Veto region with MET cut
    double NinEEMC=0, NinMuMuMC=0;//All other MC events

    for(size_t i=0; i < Vec_Files.size(); i++){
        double LumiWeight = CalcLumiWeight(Vec_Files.at(i));
        double allWeights=LumiWeight;//calculate here all the flat-weights we apply: Lumi*others*...
        if(Vec_Files.at(i).Contains("ee_") || Vec_Files.at(i).Contains("mumu_")){
            if(Vec_Files.at(i).Contains("run")){
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), "dyScaling_Looseh1");
                ApplyFlatWeights(htemp, allWeights);
                ApplyFlatWeights(htemp1, allWeights);
                if(Vec_Files.at(i).Contains("ee_")){
                    NinEE+=htemp->Integral();
                    NinEEloose+=htemp1->Integral();
                }
                if(Vec_Files.at(i).Contains("mumu_")){
                    NinMuMu+=htemp->Integral();
                    NinMuMuloose+=htemp1->Integral();
                }
                delete htemp; delete htemp1;
            }
            else if(Vec_Files.at(i).Contains("dy")){
                if(Vec_Files.at(i).Contains("50inf")){
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                    TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    ApplyFlatWeights(htemp1, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee_")){
                        NinEEDYMC+=htemp->Integral();
                        NoutEEDYMC+=htemp1->Integral();
                    }
                    if(Vec_Files.at(i).Contains("mumu_")){
                        NinMuMuDYMC+=htemp->Integral();
                        NoutMuMuDYMC+=htemp1->Integral();
                    }
                    delete htemp; delete htemp1;
                }
                else{
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee_")){   NoutEEDYMC+=htemp->Integral();}
                    if(Vec_Files.at(i).Contains("mumu_")){ NoutMuMuDYMC+=htemp->Integral();}
                    delete htemp;
                }
            }
            else{
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                ApplyFlatWeights(htemp, LumiWeight);
                if(Vec_Files.at(i).Contains("ee_")){   NinEEMC+=htemp->Integral();   }
                if(Vec_Files.at(i).Contains("mumu_")){ NinMuMuMC+=htemp->Integral(); }
                delete htemp;
            }
        }
        
        if(Vec_Files.at(i).Contains("emu_") && Vec_Files.at(i).Contains("run")){
            TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
            ApplyFlatWeights(htemp, LumiWeight);
            NinEMu+=htemp->Integral();
            delete htemp;
        }

    }
    
    
    double kee = sqrt(NinEEloose/NinMuMuloose);
//     printf("kee = sqrt(%.2f/%.2f) = %.2f\n", NinEEloose, NinMuMuloose, kee);
    
    double kmumu = sqrt(NinMuMuloose/NinEEloose);
//     printf("kmumu = sqrt(%.2f/%.2f) = %.2f\n", NinMuMuloose, NinEEloose, kmumu);
    
    double RoutinEE = NoutEEDYMC/NinEEDYMC;
//     printf("RoutinEE = %.2f/%.2f = %.2f\n", NoutEEDYMC, NinEEDYMC, RoutinEE);
    
    double RoutinMuMu = NoutMuMuDYMC/NinMuMuDYMC;
//     printf("RoutinMuMu = %.2f/%.2f = %.2f\n", NoutMuMuDYMC, NinMuMuDYMC, RoutinMuMu);
    
    double NoutMCEE = RoutinEE*(NinEE - 0.5*NinEMu*kee);
    double NoutMCMuMu = RoutinMuMu*(NinMuMu - 0.5*NinEMu*kmumu);

    double DYSFEE = NoutMCEE/NoutEEDYMC;
    double DYSFMuMu = NoutMCMuMu/NoutMuMuDYMC;
    double DYSFEMu = std::sqrt(DYSFEE * DYSFMuMu);

    std::cout << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Calculation of DY Scale Factors for '" << name << "'  at selection step "<<nameAppendix << std::endl;

    std::cout<<"DYSFEE:                 "<<DYSFEE<<std::endl;
    std::cout<<"DYSFMuMu:               "<<DYSFMuMu<<std::endl;
    std::cout<<"DYSFEMu:                "<<DYSFEMu<<std::endl;

    std::cout<<"NinEEloose:             "<<NinEEloose<<std::endl;
    std::cout<<"NinMMloose:             "<<NinMuMuloose<<std::endl;

    std::cout<<"kee:                    "<<kee<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<std::endl;
    std::cout<<"kmumu:                  "<<kmumu<<" +- "<<0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<std::endl;

    std::cout<<"Rout/Rin ee:            "<<RoutinEE<<std::endl;
    std::cout<<"Rout/Rin Mumu:          "<<RoutinMuMu<<std::endl;

    std::cout<<"Est. From Data(ee):     "<<NoutMCEE<<std::endl;
    std::cout<<"Est. From Data(mumu):   "<<NoutMCMuMu<<std::endl;

    std::cout<<"Est. From MC(ee):       "<<NoutEEDYMC<<std::endl;
    std::cout<<"Est. From MC(mumu):     "<<NoutMuMuDYMC<<std::endl;

    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << std::endl;

    if (DYSFEE < 0.2 || DYSFMuMu < 0.2) {
        std::cout << "The DY SF is too low (below 0.2). Something is probably wrong.\n";
        std::exit(1);
    }
    
    DYScale.at(0)=DYSFEE;
    DYScale.at(1)=DYSFMuMu;
    DYScale.at(2)=DYSFEMu;
    DYScale.at(3)=(DYSFEE+DYSFMuMu+DYSFEMu)/3;//not correct, but close, fix later

    std::cout<<"End DYSCALE FACTOR calculation\n"<<std::endl;

}



void Plotter::CalcDiffSystematics(TString Channel, TString Systematic, TString SystematicUp, TString SystematicDown, double flat_Syst){

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "    Preparing to Calculate " << Systematic << "-Uncertainty ... " << std::endl;

    ofstream ResultsFile;

    string ResultsFilestring = outpathResults.Data();
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(Systematic);
    ResultsFilestring.append("/");
    ResultsFilestring.append(Channel);
    gSystem->mkdir(ResultsFilestring.c_str(), true);
    ResultsFilestring.append("/");
    ResultsFilestring.append(name);
    ResultsFilestring.append("Results.txt");
    std::cout<<ResultsFilestring<<std::endl;
    ResultsFile.open(ResultsFilestring.c_str());
        
    double Xbins[XAxisbins.size()];
    for(size_t i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];} 

    if(flat_Syst > 0.0){// flat systematics will be: flat_Syst != 0.0
        for(size_t bin = 0; bin < XAxisbinCenters.size(); bin++){
            ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" SystematicError: "<<flat_Syst<<std::endl;
        }
    }

    double Sys_Error;//, Sum_Errors;

    TString newname = name;
    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp",3,"",0);
    }
    
    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}
    
    // DAVID. Guckst du hier! 
    if ( doUnfolding == true && flat_Syst == 0.0 && channelType != 3) {//Only do it in non combined channel
    
        // SVD Helper Class
        DilepSVDFunctions mySVDFunctions; 
        mySVDFunctions.SetOutputPath(outpath);
                    
        // Variables for the needed histograms
        TH1D* theDataHist = nullptr;
        TH1D* theBgrHist = nullptr;
        TH1D* theBgrHistUp = nullptr;
        TH1D* theBgrHistDown = nullptr;
        TH1D* theTtBgrHist = nullptr;
        TH1D* theTtBgrHistUp = nullptr;
        TH1D* theTtBgrHistDown = nullptr;
        TH1D* theRecHist = nullptr;
        TH1D* theRecHistUp = nullptr;
        TH1D* theRecHistDown = nullptr;
        TH1D* theGenHist = nullptr;
        TH1D* theGenHistUp = nullptr;
        TH1D* theGenHistDown = nullptr;
        TH2D* theRespHist = nullptr;
        TH2D* theRespHistUp = nullptr;
        TH2D* theRespHistDown = nullptr; 

        // DAVID:
        // Data, Signal and Background
        // can be obtained from vectors that Tyler fills.
        // These are the vectors
        // "hists", "systhistsUp" amd "systhistsDown"
        // Notice, that the DY Background will be scaled with
        // the DYScale.

        TString filename = "preunfolded/Nominal/"+Channel+"/"+name+"_UnfoldingHistos.root";
        theDataHist = fileReader->GetClone<TH1D>(filename, "aDataHist");
        theBgrHist =  fileReader->GetClone<TH1D>(filename, "aBgrHist");
        theTtBgrHist =  fileReader->GetClone<TH1D>(filename, "aTtBgrHist");
        theRecHist =  fileReader->GetClone<TH1D>(filename, "aRecHist");
        theGenHist =  fileReader->GetClone<TH1D>(filename, "aGenHist");
        theRespHist = fileReader->GetClone<TH2D>(filename, "aRespHist");

        TString filenameUp = "preunfolded/"+SystematicUp+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
        theBgrHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aBgrHist");
        theTtBgrHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aTtBgrHist");
        theRecHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aRecHist");
        theGenHistUp =  fileReader->GetClone<TH1D>(filenameUp, "aGenHist");
        theRespHistUp =  fileReader->GetClone<TH2D>(filenameUp, "aRespHist");

        TString filenameDown = "preunfolded/"+SystematicDown+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
        theBgrHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aBgrHist");
        theTtBgrHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aTtBgrHist");
        theRecHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aRecHist");
        theGenHistDown =  fileReader->GetClone<TH1D>(filenameDown, "aGenHist");
        theRespHistDown =  fileReader->GetClone<TH2D>(filenameDown, "aRespHist");

        // Get the binning 
        double* theBins = Xbins;
        int numberBins = bins;

        // Names and Labels
        TString channelLabelStr(channelLabel[channelType]);
        TString theChannelName = Channel;
        //if ( channelLabelStr.Contains("#mu#mu")  ) theChannelName = "mumu";
        //if ( channelLabelStr.Contains("e#mu")    ) theChannelName = "emu";
        //if ( channelLabelStr.Contains("ee")      ) theChannelName = "ee";
        //if ( channelLabelStr.Contains("Dilepton Combined")    ) theChannelName = "combined";
        TString theParticleName = "";
        if ( name.Contains("Lepton")  ) theParticleName = "Leptons";
        if ( name.Contains("LLBar")   ) theParticleName = "LepPair";
        if ( name.Contains("Top")     ) theParticleName = "TopQuarks";
        if ( name.Contains("TTBar")   ) theParticleName = "TtBar";
        if ( name.Contains("BJet")    ) theParticleName = "BJets";
        if ( name.Contains("BBBar")   ) theParticleName = "BBbar";
        if ( name.Contains("JetMult") ) theParticleName = "Jets";
        if ( name.Contains("ExtraJet" ) ) theParticleName = "ExtraJets";
        if ( name.Contains("DeltaRJet12") ) theParticleName = "ExtraJets";
        TString theQuantityName = "";
        if ( name.Contains("pT")      ) theQuantityName = "Pt";
        if ( name.Contains("Eta")     ) theQuantityName = "Eta";
        if ( name.Contains("Rapidity")) theQuantityName = "Rapidity";
        if ( name.Contains("Mass")    ) theQuantityName = "Mass";
        if ( name.Contains("JetMult")    ) theQuantityName = "Mult";
        if ( name.Contains("ExtraJetEta") ) theQuantityName = "Eta";
        if ( name.Contains("ExtraJetpT") ) theQuantityName = "Pt";
        if ( name.Contains("DeltaRJet12") ) theQuantityName = "DeltaR";
        TString theSpecialPostfix = "";
        theSpecialPostfix = name;
        if ( specialComment.CompareTo("Standard") != 0 ) {
            //theSpecialPostfix = specialComment;
        } 

        if ( name.Contains("ExtraJetEta2") ) theSpecialPostfix = "JExtra2";
        if ( name.Contains("ExtraJetpT2") ) theSpecialPostfix = "JExtra2";
        if ( name.Contains("ExtraJetEta3") ) theSpecialPostfix = "JExtra3";
        if ( name.Contains("ExtraJetpT3") ) theSpecialPostfix = "JExtra3";

        TString theSystematicName = Systematic; 

        std::cout << std::endl;
        std::cout << std::endl;
    
        // Get the integrals for the normalization
        double totalDataEventsNom  = TopSVDFunctions::SVD_Integral1D(theDataHist, 0, false);
        double totalBgrEventsNom   = TopSVDFunctions::SVD_Integral1D(theBgrHist, 0, false);
        double totalBgrEventsUp    = TopSVDFunctions::SVD_Integral1D(theBgrHistUp, 0, false);
        double totalBgrEventsDown  = TopSVDFunctions::SVD_Integral1D(theBgrHistDown, 0, false);
        double totalTtBgrEventsNom = TopSVDFunctions::SVD_Integral1D(theTtBgrHist, 0, false);
        double totalTtBgrEventsUp  = TopSVDFunctions::SVD_Integral1D(theTtBgrHistUp, 0, false);
        double totalTtBgrEventsDown= TopSVDFunctions::SVD_Integral1D(theTtBgrHistDown, 0, false);
        double totalRecEventsNom   = TopSVDFunctions::SVD_Integral1D(theRecHist, 0, false);
        double totalRecEventsUp    = TopSVDFunctions::SVD_Integral1D(theRecHistUp, 0, false);
        double totalRecEventsDown  = TopSVDFunctions::SVD_Integral1D(theRecHistDown, 0, false);
        double totalGenEventsNom   = TopSVDFunctions::SVD_Integral1D(theGenHist, 0, false);
        double totalGenEventsUp    = TopSVDFunctions::SVD_Integral1D(theGenHistUp, 0, false);
        double totalGenEventsDown  = TopSVDFunctions::SVD_Integral1D(theGenHistDown, 0, false);

        // UNFOLDING OF SYSTEMATICS
        // Retrieve histograms with the unfolded quantities.
        // Note: The unfolded histograms have additional side bins!
        // Keep this in mind when accessing bin content via indices 
        TH1D* symmSysErrors = NULL;

        mySVDFunctions.SVD_DoUnfoldSys(
                                    theDataHist,
                                    theBgrHist, theBgrHistUp, theBgrHistDown, 
                                    theTtBgrHist, theTtBgrHistUp, theTtBgrHistDown, 
                                    theGenHist, theGenHistUp, theGenHistDown, 
                                    theRecHist, theRecHistUp, theRecHistDown, 
                                    theRespHist, theRespHistUp, theRespHistDown, 
                    totalDataEventsNom, 
                    totalBgrEventsNom,  totalBgrEventsUp,  totalBgrEventsDown, 
                    totalTtBgrEventsNom,  totalTtBgrEventsUp,  totalTtBgrEventsDown, 
                    totalRecEventsNom,  totalRecEventsUp,  totalRecEventsDown, 
                    totalGenEventsNom,  totalGenEventsUp,  totalGenEventsDown,  
                                    theBins, numberBins,
                                    symmSysErrors,  
                                    theChannelName, theParticleName, theQuantityName, theSpecialPostfix, theSystematicName
                                    ); 
        
        
        
        
        //Symetrize Eta and Rapidity distributions
        if (theQuantityName == "Eta" || theQuantityName == "Rapidity" ){
            for(int j=0; j<(int) symmSysErrors->GetNbinsX(); ++j){
                std::cout<<"In bin "<<j<<" binCenter "<<symmSysErrors->GetBinCenter(j+1)<<" Content "<<symmSysErrors->GetBinContent(j+1)<<std::endl;
            }

            int Nbins = theDataHist->GetNbinsX();
            std::cout<<"Nbins in "<<symmSysErrors->GetName()<<" = "<<symmSysErrors->GetNbinsX()<<std::endl;
            //There are 2 extra bins coming from the unfolding ==>  skip the underflow+1 bin from left and and overflow+1 bin from right
            for(int i=0; i<Nbins; ++i){
                std::cout<<"(2nd loop) In bin "<<i<<" binCenter "<<symmSysErrors->GetBinCenter(i+2)<<" Content "<<symmSysErrors->GetBinContent(i+2)<<std::endl;
                std::cout<<"                     binCenter "<<symmSysErrors->GetBinCenter(Nbins-i+1)<<" Content "<<symmSysErrors->GetBinContent(Nbins-i+1)<<std::endl;
                Sys_Error = 0.5*(symmSysErrors->GetBinContent(i+2)+symmSysErrors->GetBinContent(Nbins+1-i));
                std::cout<<"Symetrized error "<<Sys_Error<<std::endl;

                //std::cout<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[i]<<" bin: "<<Xbins[i]<<" to "<<Xbins[i+1]<<" SystematicError: "<<Sys_Error<<std::endl;
                ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[i]<<" bin: "<<Xbins[i]<<" to "<<Xbins[i+1]<<" SystematicRelError: "<<Sys_Error<<std::endl;
            }
        }
        else{
            // Save the shifts in Tyler's triple-matrix ...
            for(Int_t bin = 0; bin < theDataHist->GetNbinsX(); ++bin) {
                Sys_Error = symmSysErrors->GetBinContent(bin+2); // Keep in mind the extra layer of OF bins
                // Save it
                //std::cout<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" SystematicError: "<<Sys_Error<<std::endl;
                ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" SystematicRelError: "<<Sys_Error<<std::endl;
            }
        }
        
        //delete objects
        delete symmSysErrors;
        delete theDataHist;
        delete theBgrHist;
        delete theBgrHistUp;
        delete theBgrHistDown;
        delete theTtBgrHist;
        delete theTtBgrHistUp;
        delete theTtBgrHistDown;
        delete theRecHist;
        delete theRecHistUp;
        delete theRecHistDown;
        delete theGenHist;
        delete theGenHistUp;
        delete theGenHistDown;
        delete theRespHist;
        delete theRespHistUp;
        delete theRespHistDown; 
    }
    
    
    
    if (channelType==3){//for 'combined' channel: do an statistical combination of the the 3 independent channels using the central value's statistical error as weight
        
        TString eefilename="UnfoldingResults/Nominal/ee/"+name+"Results.txt";
        TString emufilename="UnfoldingResults/Nominal/emu/"+name+"Results.txt";
        TString mumufilename="UnfoldingResults/Nominal/mumu/"+name+"Results.txt";
        TString combfilename="UnfoldingResults/Nominal/combined/"+name+"Results.txt";
        TString eeErrfilename="UnfoldingResults/"+Systematic+"/ee/"+name+"Results.txt";
        TString emuErrfilename="UnfoldingResults/"+Systematic+"/emu/"+name+"Results.txt";
        TString mumuErrfilename="UnfoldingResults/"+Systematic+"/mumu/"+name+"Results.txt";
        
        //check the existence of the file
        if(gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename) ||
           gSystem->AccessPathName(eeErrfilename) || gSystem->AccessPathName(emuErrfilename) || gSystem->AccessPathName(mumuErrfilename)){
            std::cout<<"WARNING (in CalcDiffSystematics)!!"<<std::endl;
            std::cout<<"One of the input files you use for the combined XSection measurement doesn't exist!! Exiting!!"<<std::endl;
            exit(22);
        }
        
        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);
        ifstream ResultsComb(combfilename);
        ifstream SystErrEE(eeErrfilename);
        ifstream SystErrEMu(emuErrfilename);
        ifstream SystErrMuMu(mumuErrfilename);
        
        double perChannelDiffXSecPlot[4][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecSystError[3][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecStatError[4][bins]; //perChannelDiffXSecStatError[channel][bin]
        double perChannelGenDiffXSec[4][bins];       //perChannelGenDiffXSec[channel][bin]
        TString Dummy="";
        
        for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
            ResultsEE>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[0][bin]>>Dummy>>perChannelDiffXSecStatError[0][bin]>>Dummy>>perChannelGenDiffXSec[0][bin];
            ResultsEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[1][bin]>>Dummy>>perChannelDiffXSecStatError[1][bin]>>Dummy>>perChannelGenDiffXSec[1][bin];
            ResultsMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[2][bin]>>Dummy>>perChannelDiffXSecStatError[2][bin]>>Dummy>>perChannelGenDiffXSec[2][bin];
            ResultsComb>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[3][bin]>>Dummy>>perChannelDiffXSecStatError[3][bin]>>Dummy>>perChannelGenDiffXSec[3][bin];
            SystErrEE>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecSystError[0][bin];
            SystErrEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecSystError[1][bin];
            SystErrMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecSystError[2][bin];
        }
        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

//         double binWidth[bins];
        double DiffXSecVecVaried[bins];
        //do the actual combined Diff.XSection calculation
        for(int i=0; i<bins; i++){
//             binWidth[i] = Xbins[i+1]-Xbins[i];
            for(int j=0; j<3; j++){//check if any stat error is 0, in this case set their contribution to 0!!
                if(perChannelDiffXSecStatError[j][i] == 0){
                    perChannelDiffXSecStatError[j][i] = 1e100;
                    perChannelDiffXSecPlot[j][i] = 0;
                }
            }
            DiffXSecVecVaried[i] =((( 1 + perChannelDiffXSecSystError[0][i]) * perChannelDiffXSecPlot[0][i]/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))
                                  +(( 1 + perChannelDiffXSecSystError[1][i]) * perChannelDiffXSecPlot[1][i]/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))
                                  +(( 1 + perChannelDiffXSecSystError[2][i]) * perChannelDiffXSecPlot[2][i]/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])))/
                                  ((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                   (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                   (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])));

            double Sys_Error = std::fabs(DiffXSecVecVaried[i]-perChannelDiffXSecPlot[3][i])/perChannelDiffXSecPlot[3][i];
            ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[i]<<" bin: "<<Xbins[i]<<" to "<<Xbins[i+1]<<" SystematicRelError: "<<Sys_Error<<std::endl;
        }
    }
    
    ResultsFile.close();
}

Plotter::Plotter()
{
    gErrorIgnoreLevel = 1001;
    name="defaultName";
    specialComment="Standard";
    rangemin=0;
    rangemax=3;
    YAxis="N_{events}";
    initialized=false;
    datafiles = 0;

    channelLabel.insert(channelLabel.begin(),4, "");

    // DAVID
    outpath = "";
    outpathPlots = "Plots";
    outpathPlots = "UnfoldingResults";
    subfolderChannel = "";
    subfolderSpecial = "";

    fileReader = RootFileReader::getInstance();
}

void Plotter::setOptions(TString name_, TString specialComment_, TString YAxis_, TString XAxis_, int rebin_, bool doDYScale_, bool logX_, bool logY_, double ymin_, double ymax_, double rangemin_, double rangemax_, int bins_, std::vector<double> XAxisbins_, std::vector<double> XAxisbinCenters_)
{
    name=name_; //Variable name you want to plot
    specialComment=specialComment_; 
    YAxis=YAxis_; //Y-axis title
    XAxis=XAxis_; //X-axis title
    rebin=rebin_; //Nr. of bins to be merged together
    doDYScale = doDYScale_; //Apply DY scale factor?
    logX = logX_; //Draw X-axis in Log scale
    logY = logY_; //Draw Y-axis in Log scale
    ymin = ymin_; //Min. value in Y-axis
    ymax = ymax_; //Max. value in Y-axis
    rangemin=rangemin_; //Min. value in X-axis
    rangemax=rangemax_; //Max. value in X-axis
    bins=bins_; //Number of bins to plot
    XAxisbins.clear(); 
    XAxisbins = XAxisbins_; // Bins edges=bins+1
    XAxisbinCenters.clear();
    XAxisbinCenters = XAxisbinCenters_; //Central point for BinCenterCorrection=bins

    //Modify the X/Y-axis labels
    if(XAxis.Contains("band#bar{b}")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("band#bar{b}",11,"b and #bar{b}",13);
    }
    if(XAxis.Contains("tand#bar{t}")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("tand#bar{t}",11,"t and #bar{t}",13);
    }
    if(XAxis.Contains("l^{+}andl^{-}")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("l^{+}andl^{-}",13,"l^{+} and l^{-}",15);
    }
    if(XAxis.Contains("p_{T}^{t}(t#bar{t}c.o.m.)")){//Histogram naming convention has to be smarter
        XAxis.ReplaceAll("p_{T}^{t}(t#bar{t}c.o.m.)",25,"p_{T}^{t}(t#bar{t} c.o.m.)",26);
    }
    if(YAxis.Contains("Toppairs")){
        YAxis.ReplaceAll("Toppairs",8,"Top-quark pairs",15);
    }
    if(YAxis.Contains("Topquarks")){
        YAxis.ReplaceAll("Topquarks",9, "Top quarks",10);
    }
    if(YAxis.Contains("Numberof")){
        YAxis.ReplaceAll("Numberof", 8, "Number of ",10);
    }

    DYScale.insert(DYScale.begin(), 4, 1.);//Initialize the DY scale-factor to (1., 1., 1., 1.)
}


void Plotter::setDataSet(std::vector<TString> dataset_, std::vector<double> scales_, std::vector<TString> legends_, std::vector<int> colors_, TString DYEntry_)
{
    dataset.clear();
    scales.clear();
    legends.clear();
    legendsSyst.clear();
    colors.clear();
    dataset=dataset_;
    scales=scales_;
    legends=legends_;
    colors=colors_;
    DYEntry=DYEntry_;
}

void Plotter::setDataSet(TString mode, TString Systematic)
{
    initialized=false;
    legendsSyst.clear();

    if(channelLabel.size()<4){channelLabel.insert(channelLabel.begin(), 4, "");}

    if(mode =="ee"){channelType=0;channelLabel.at(0)="ee";}
    if(mode =="mumu"){channelType=1;channelLabel.at(1)="#mu#mu";}
    if(mode =="emu"){channelType=2;channelLabel.at(2)="e#mu";}
    if(mode =="combined"){channelType=3;channelLabel.at(3)="Dilepton Combined";}

    // Set dataset specific subfolders
    outpathPlots = "./Plots";
    outpathResults = "./UnfoldingResults";
    subfolderChannel = mode;
    subfolderChannel.Prepend("/");
    subfolderSpecial = "";
    if ( specialComment.CompareTo("Standard") != 0 ) {
        //subfolderSpecial = specialComment.Prepend("/");
    }

    DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";

    if(Systematic.Contains("DY_") || Systematic.Contains("BG_")){Systematic = "Nominal";}//We just need to vary the nominal DY and BG systematics

    TString histoListName = "FileLists/HistoFileList_"+Systematic+"_"+mode+".txt";
    std::cout << "reading " << histoListName << std::endl;
    ifstream FileList(histoListName);
    if (FileList.fail()) {
        std::cerr << "Error reading " << histoListName << std::endl;
        exit(1);
    }
    TString filename;
    datafiles=0;

    dataset.clear();
    legends.clear();
    colors.clear();

    while(!FileList.eof()){
        FileList>>filename;
        if(filename==""){continue;}//Skip empty lines
        dataset.push_back(filename);
        if(filename.Contains("run")){legends.push_back("Data"); colors.push_back(kBlack);datafiles++;}
        else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} Signal"); colors.push_back(kRed+1);}
        else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} Other"); colors.push_back(kRed-7);}
        else if(filename.Contains("single")){legends.push_back("Single Top"); colors.push_back(kMagenta);}
        else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); colors.push_back(10);}
        else if(filename.Contains("dytautau")){legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
        else if(filename.Contains("dymumu")||filename.Contains("dyee")){legends.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colors.push_back(kAzure-2);}
        else if(filename.Contains("wtolnu")){legends.push_back("W+Jets"); colors.push_back(kGreen-3);}
        else if(filename.Contains("qcd")){legends.push_back("QCD Multijet"); colors.push_back(kYellow);}
        else if(filename.Contains("ttbarZ") ||filename.Contains("ttbarW") || filename.Contains("ttgjets")){legends.push_back("t#bar{t}+V"); colors.push_back(kOrange-2);}
    }
    FileList.close();
}


void Plotter::ListOfSystematics(set<TString> listofsystematics){
    
    //List of systematics valid for running. Please check (or add any new systematic) the function 'ListOfSystematics' in Histo.C
    
    ListOfSyst = listofsystematics;
    
}


std::vector<TString> Plotter::InputFileList(TString mode, TString Systematic)
{
    //Hard code the input file list. This functions handles correctly the data files. So no symlinks are anymore needed
    //The file structure must be: selectionRoot/SYSTEMATIC/channel/XXX.root
    //This function will take the data (aka: run2012A, run2012B, run2012C, ...) from the Nominal systematic for ALL the systematic
    //This function will take every sample, except the ttbar, from Nominal for any signal systematic: POWHEG, MCATNLO, MATCH, SCALE, MASS
    
    std::vector<TString> FileVector;
    FileVector.clear();
    
    if(ListOfSyst.find(Systematic) == ListOfSyst.end()){
        std::cout<<"WARNING(in InputFileList)!!! Using a non valid systematic: "<<Systematic<<std::endl;
        std::cout<<"Please use one among the supported ones:"<<std::endl;
        for(set<TString>::iterator iter = ListOfSyst.begin(); iter!=ListOfSyst.end(); ++iter){ std::cout<<*iter<<std::endl;}
        return FileVector;
    }
    
    if( mode.CompareTo("combined") && mode.CompareTo("ee") && mode.CompareTo("emu") && mode.CompareTo("mumu")){
        std::cout<<"The decay channel you provided is not supported."<<std::endl;
        std::cout<<"Please use: ee, emu, mumu, combined"<<std::endl;
        return FileVector;
    }
    
    if(!mode.CompareTo("combined")){
        std::vector<TString> eemode   = Plotter::InputFileList(TString("ee"), Systematic);
        std::vector<TString> emumode  = Plotter::InputFileList(TString("emu"), Systematic);
        std::vector<TString> mumumode = Plotter::InputFileList(TString("mumu"), Systematic);
        FileVector.insert(FileVector.end(), eemode.begin(), eemode.end());
        FileVector.insert(FileVector.end(), emumode.begin(), emumode.end());
        FileVector.insert(FileVector.end(), mumumode.begin(), mumumode.end());
        //for (auto s: FileVector) std::cout << "FileVector = " << s << "\n";
        //shouldn't this be sorted? i.e. first runs of all channels, then ...
        //need to ask Ivan !!!FIXME
        return FileVector;
    }

    //data is only stored in the Nominal directory
    TString nominalPath = TString("selectionRoot/Nominal/") + mode + "/" + mode;
    if (!doClosureTest) {
        FileVector.push_back(nominalPath + "_run2012A.root");
        FileVector.push_back(nominalPath + "_run2012B.root");
        FileVector.push_back(nominalPath + "_run2012C.root");
        FileVector.push_back(nominalPath + "_run2012D.root");
    } else {
        FileVector.push_back(nominalPath + "_ttbarsignalplustau_fakerun_nominal.root");
    }
    
    //MC depends on the specific Systematic: Signal systematics only use different signal samples
    TString tempName;
    if(!Systematic.CompareTo("JER_UP") || !Systematic.CompareTo("JER_DOWN") || !Systematic.CompareTo("JES_UP") || !Systematic.CompareTo("JES_DOWN") ||
       !Systematic.CompareTo("PU_UP") || !Systematic.CompareTo("PU_DOWN") || !Systematic.CompareTo("BTAG_UP") || !Systematic.CompareTo("BTAG_DOWN") ||
       !Systematic.CompareTo("BTAG_PT_UP") || !Systematic.CompareTo("BTAG_PT_DOWN") || !Systematic.CompareTo("BTAG_ETA_UP") || !Systematic.CompareTo("BTAG_ETA_DOWN")
    ){
        tempName = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode;
    }
    else{
        tempName = TString("selectionRoot/Nominal/") + mode + "/" + mode;
    }
    
    FileVector.push_back(tempName + "_dyee1050.root");
    FileVector.push_back(tempName + "_dyee50inf.root");
    FileVector.push_back(tempName + "_dymumu1050.root");
    FileVector.push_back(tempName + "_dymumu50inf.root");
    FileVector.push_back(tempName + "_singleantitop_tw.root");
    FileVector.push_back(tempName + "_singletop_tw.root");
    FileVector.push_back(tempName + "_wtolnu.root");
    FileVector.push_back(tempName + "_wwtoall.root");
    FileVector.push_back(tempName + "_wztoall.root");
    FileVector.push_back(tempName + "_zztoall.root");
    FileVector.push_back(tempName + "_dytautau1050.root");
    FileVector.push_back(tempName + "_dytautau50inf.root");
    FileVector.push_back(tempName + "_qcdbcem2030.root");
    FileVector.push_back(tempName + "_qcdbcem3080.root");
    FileVector.push_back(tempName + "_qcdbcem80170.root");
    FileVector.push_back(tempName + "_qcdem2030.root");
    FileVector.push_back(tempName + "_qcdem3080.root");
    FileVector.push_back(tempName + "_qcdem80170.root");
    FileVector.push_back(tempName + "_qcdmu15.root");
//     FileVector.push_back(tempName + "_ttbarW.root");
//     FileVector.push_back(tempName + "_ttbarZ.root");
//     FileVector.push_back(tempName + "_ttgjets.root");
    
    
    TString ttbarsignalplustau = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode + "_ttbarsignalplustau.root";
    TString ttbarbgviatau      = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode + "_ttbarbgviatau.root";
    TString ttbarbg            = TString("selectionRoot/") + Systematic + "/" + mode + "/" + mode + "_ttbarbg.root";

    //Add extra filename for signal systematic filenames
    if(!Systematic.CompareTo("POWHEG") || !Systematic.CompareTo("MCATNLO") || 
       !Systematic.CompareTo("MASS_UP") || !Systematic.CompareTo("MASS_DOWN") ||
       !Systematic.CompareTo("MATCHINGUP") || !Systematic.CompareTo("MATCHINGDOWN") ||
       !Systematic.CompareTo("MATCH_UP") || !Systematic.CompareTo("MATCH_DOWN") ||
       !Systematic.CompareTo("SCALE_UP") || !Systematic.CompareTo("SCALE_DOWN")){
        
        if(!Systematic.CompareTo("MATCHINGUP") || !Systematic.CompareTo("MATCHINGDOWN")){ Systematic.Remove(5, 3);}
        TString tmpSyst = Systematic;
        tmpSyst.Prepend("_").ToLower();
        
        ttbarsignalplustau.ReplaceAll(".root", tmpSyst + ".root");
        ttbarbg.ReplaceAll(".root", tmpSyst + ".root");
        ttbarbgviatau.ReplaceAll(".root", tmpSyst + ".root");
    }
    FileVector.push_back(ttbarsignalplustau);
    FileVector.push_back(ttbarbg);
    FileVector.push_back(ttbarbgviatau);
    
    return FileVector;
    
}


bool Plotter::fillHisto()
{   
    TH1::AddDirectory(kFALSE);
    if (initialized) { return true; }
    hists.clear();
    for(unsigned int i=0; i<dataset.size(); i++){
//         std::cout << i << ": " << dataset.at(i) << std::endl;
        TH1D *hist = fileReader->GetClone<TH1D>(dataset.at(i), name, true);
        if (!hist) return false;
        if (!name.Contains("_step") && !name.Contains("bcp_") && !name.Contains("Lead") && !name.EndsWith("bkr") && !name.EndsWith("akr")
            && (name.Contains("Lepton") || name.Contains("BJet") || name.Contains("Top")))
        {
            TString stemp = name;
            if(name.Contains("Lepton"))    {stemp.ReplaceAll("Lepton",6,"AntiLepton",10);}
            else if(name.Contains("BJet")) {stemp.ReplaceAll("BJet",4,"AntiBJet",8);}
            else if(name.Contains("Top"))  {stemp.ReplaceAll("Top",3,"AntiTop",7);}
            const TH1D* other = fileReader->Get<TH1D>(dataset.at(i), stemp);
            if (other) hist->Add(other); 
            else std::cerr << "Cannot find corresponding anti-quantity histogram: "<<stemp<<std::endl;
        }

        //Rescaling to the data luminosity
        double LumiWeight = CalcLumiWeight(dataset.at(i));
        //std::cout << "File " << dataset.at(i) << " has weight " << LumiWeight << "\n";

        ApplyFlatWeights(hist, LumiWeight);

        common::setHHStyle(*gStyle);
        hists.push_back(*hist);
        delete hist;
    }
    if (doClosureTest) {
        for (size_t i = 1; i<dataset.size(); ++i) {
            if (!dataset.at(i).Contains("ttbarsignalplustau")){
                hists[0].Add(&hists[i]);
            }
        }
    }
    initialized=true;
    return true;
}

void Plotter::addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis)
{
    if (!addToThis) addToThis = addThis;
    else {
        addToThis->Add(addThis);
        delete addThis;
    }
}


void Plotter::write(TString Channel, TString Systematic) // do scaling, stacking, legending, and write in file 
{
    setDataSet(Channel,Systematic);
    if (!fillHisto()) return;
    if (hists.size() == 0) { 
        std::cerr << "***ERROR! No histograms available! " << Channel << "/" << Systematic << std::endl; 
        exit(11);
    }
        
    auto c = make_shared<TCanvas>("","");
    auto stack = make_shared<THStack>("def", "def");
    TLegend *leg  = new TLegend();

    std::vector<TH1*> drawhists;
    drawhists.resize(hists.size());
    
    std::stringstream ss;
    ss << DYScale[channelType];
    TString scale;
    scale=(TString)ss.str();
    int legchange=0;
    leg->Clear();
    c->Clear();
    c->SetName("");

    c->SetTitle("");
    TH1* aDataHist = NULL;
    TH1* aBgrHist = NULL;
    TH1* aTtBgrHist = NULL;
    TH1* aRecHist = NULL;
    TH1* aGenHist = NULL; 
    TH1* aRespHist = NULL;

    double Xbins[XAxisbins.size()];
    TString newname = name;

    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
        newname.ReplaceAll("Hyp",3,"",0);
    }

    bool init=false;

    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}

    std::unique_ptr<TH1> sumttbar;
    std::unique_ptr<TH1> allttbar;

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        drawhists[i]=(TH1D*) hists[i].Clone();//rebin and scale the histograms
        if(rebin>1) drawhists[i]->Rebin(rebin);
        if(XAxisbins.size()>1) drawhists[i] = drawhists[i]->Rebin(bins,"tmp",Xbins);
        setStyle(drawhists[i], i, true);
        if (legends.at(i) == "t#bar{t} Signal") {
            if (sumttbar.get()) sumttbar->Add(drawhists[i]);
            else sumttbar = std::unique_ptr<TH1>{static_cast<TH1*>(drawhists[i]->Clone())};
        }
        if (legends.at(i).Contains("t#bar{t}")){
            if (allttbar.get()) allttbar->Add(drawhists[i]);
            else allttbar = std::unique_ptr<TH1>{static_cast<TH1*>(drawhists[i]->Clone())};
        }
    }
    for(unsigned int i=0; i<hists.size() ; ++i){ // prepare histos and leg
//         std::cout << "Legend ["<<i<<"] = " << legends.at(i) << std::endl;
    if(legends.at(i) != "Data"){
        //      drawhists[i]->Scale(12.1/5.1);

        if(XAxisbins.size()>1){//only distributions we want to unfold will have a binning vector
            if(legends.at(i) == "t#bar{t} Signal" && doUnfolding){
                TString ftemp = dataset.at(i);
                double LumiWeight = CalcLumiWeight(dataset.at(i));
                if (!init) {
                    aRespHist = fileReader->GetClone<TH2>(ftemp, "GenReco"+newname);
                    aGenHist = fileReader->GetClone<TH1D>(ftemp, "VisGen"+newname);
                    //Rebin(bins,"aTtBgrHist",Xbins);
                    if(!newname.Contains("Lead") && (newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet"))){
                        aRespHist->Add(fileReader->Get<TH2>(ftemp, "GenRecoAnti"+newname));
                        aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGenAnti"+newname));
                    }
                    init =true;
                } else {//account for more than one signal histogram
                    aRespHist->Add(fileReader->Get<TH2>(ftemp, "GenReco"+newname));
                    aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGen"+newname));
                    if(!newname.Contains("Lead") && (newname.Contains("Lepton")||newname.Contains("Top")||newname.Contains("BJet"))){
                        aGenHist->Add(fileReader->Get<TH1D>(ftemp, "VisGenAnti"+newname));
                        aRespHist->Add(fileReader->Get<TH2>(ftemp, "GenRecoAnti"+newname));
                    }
                }
                aRespHist->Scale(LumiWeight);
                aGenHist->Scale(LumiWeight);
            }
            //std::cout<<"Legend: "<<legends.at(i)<<std::endl;
            if(legends.at(i) == "t#bar{t} Signal"){
                signalHist = i; //What is happening for the combined channel?? only 1 integer value??
                aRecHist = drawhists[i]->Rebin(bins,"aRecHist",Xbins); //What happens in the combined channel? Are the 3 channels summed up?? NO!!
                //std::cout<<"Added "<<legends.at(i)<<" to aRecHist"<<std::endl;
            }else if(legends.at(i) == "t#bar{t} Other"){//IMPORTANT: TTbar Other are added to the ttbarbackground histogram AND the Background Hist gram
                addAndDelete_or_Assign(aTtBgrHist, drawhists[i]->Rebin(bins,"aTtBgrHist",Xbins));
                addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins,"aTtBgrHist",Xbins));
            }else if((legends.at(i) == DYEntry)){
                drawhists[i]->Scale(DYScale.at(channelType));

                //Here we take into account the systematic shifts needed for DY systematic because it only modifies the nominal dataset
                if(Systematic == "DY_UP"){
                    drawhists[i]->Scale(1.3);
                }
                if(Systematic == "DY_DOWN"){
                    drawhists[i]->Scale(0.7);
                }
                addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
            }else{
                //Here we take into account the systematic shifts needed for BG systematic because it only modifies the nominal dataset
                if(Systematic == "BG_UP"){
                    drawhists[i]->Scale(1.3);
                }
                if(Systematic == "BG_DOWN"){
                    drawhists[i]->Scale(0.7);
                }
                addAndDelete_or_Assign(aBgrHist, drawhists[i]->Rebin(bins,"aBgrHist",Xbins));
            }
        }

        if(i > 1){
            if(legends.at(i) != legends.at(i-1)){
                legchange = i;
                if((legends.at(i) == DYEntry)&& DYScale.at(channelType) != 1) leg->AddEntry(drawhists[i], legends.at(i),"f");
                else leg->AddEntry(drawhists[i], legends.at(i) ,"f");
            }else{
                drawhists[legchange]->Add(drawhists[i]);
            }
        }
        if(i!=(hists.size()-1)){
            if(legends.at(i)!=legends.at(i+1)){
                drawhists[i]->SetLineColor(1);
            }
        }else{
            drawhists[i]->SetLineColor(1);
        }
        if(legends.at(i) != legends.at(i-1) ){
            drawhists[i]->SetLineColor(1);
            if(!legends.at(i).Contains("QCD") || addQCDToControlPlot()){stack->Add(drawhists[i]);};

        }
    }
    else{
        if(i==0) leg->AddEntry(drawhists[i], legends.at(i) ,"pe");
        if(i>0){
            if(legends.at(i) != legends.at(i-1)){
                leg->AddEntry(drawhists[i], legends.at(i) ,"pe");
            }
            if(legends.at(i) == legends.at(0)){
                drawhists[0]->Add(drawhists[i]);
            }
        }
    }
    }

    if(XAxisbins.size()>1 && doUnfolding){//only distributions we want to unfold will have a binning vector
        aDataHist = drawhists[0]->Rebin(bins,"aDataHist",Xbins);

        TString outdir = common::assignFolder("preunfolded", Channel, Systematic);
        TFile *f15 = new TFile(outdir.Copy()+name+"_UnfoldingHistos.root","RECREATE");
        aDataHist->Write("aDataHist"); delete aDataHist;
        aTtBgrHist->Write("aTtBgrHist"); delete aTtBgrHist;
        aBgrHist->Write("aBgrHist"); delete aBgrHist;
        aGenHist->Write("aGenHist"); delete aGenHist;
        aRespHist->Write("aRespHist"); delete aRespHist;
        aRecHist->Write("aRecHist"); delete aRecHist;

        f15->Close();
        delete f15;
    }
    if(doUnfolding) return;
    TLegend *leg1 = (TLegend*)leg->Clone("leg1");
    TLegend *leg2 = (TLegend*)leg->Clone("leg2");
    setControlPlotLegendStyle(drawhists, legends, leg, leg1, leg2);

    if(name.Contains("HypjetMultiXSec")){

        double InclusiveXsectionWrite[4], InclusiveXsectionStatErrorWrite[4];
        CalcXSec(dataset, InclusiveXsectionWrite, InclusiveXsectionStatErrorWrite, Systematic,"");

        Plotter::MakeTable(Channel, Systematic);
        if(channelType==3) Plotter::PlotXSec(Channel);
    }

    drawhists[0]->SetMinimum(ymin);

    if(rangemin!=0 || rangemax!=0) {drawhists[0]->SetAxisRange(rangemin, rangemax, "X");}

    if(logY)c->SetLogy();
    if(ymax==0){
        if(logY){drawhists[0]->SetMaximum(18  * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));}
        else    {drawhists[0]->SetMaximum(1.5 * drawhists[0]->GetBinContent(drawhists[0]->GetMaximumBin()));}
    }
    else{drawhists[0]->SetMaximum(ymax);}

    if(name.Contains("HypTopRapidity") || name.Contains("HypTTBarRapidity")) drawhists[0]->GetXaxis()->SetNdivisions(511);

    drawhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    TGaxis::SetMaxDigits(4);


    //Removal of extra ticks in JetMult plots
    if(name.Contains("HypJetMultpt")) {
        drawhists[0]->GetXaxis()->SetNdivisions(drawhists[0]->GetNbinsX(),0,0, 1);

        TString TitBin = "";
        for(int bin = 1; bin <= drawhists[0]->GetNbinsX(); bin++) {
            if( bin == drawhists[0]->GetNbinsX()) {TitBin += "#geq"; TitBin += drawhists[0]->GetBinCenter(bin); drawhists[0]->GetXaxis()->SetBinLabel(bin,TitBin);}
            else{TitBin += drawhists[0]->GetBinCenter(bin);
            drawhists[0]->GetXaxis()->SetBinLabel(bin,TitBin);
            }
            TitBin  = "";
        }
    }

    //Add the binwidth to the yaxis in yield plots

    TString ytitle = drawhists[0]->GetYaxis()->GetTitle();
    double binwidth = drawhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binwidth;

    if(name.Contains("Rapidity") || name.Contains("Eta") || name.Contains("Phi") || (name.Contains("Mass") && name.Contains("1st"))){ytitle.Append(" / ").Append(width.str());}
    else if(name.Contains("pT") || (name.Contains("Mass") && ! name.Contains("1st")) || name.Contains("mass") || name.Contains("MET") || name.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
    drawhists[0]->GetYaxis()->SetTitle(ytitle);
    drawhists[0]->Draw("e1");
    gStyle->SetEndErrorSize(0);

    stack->Draw("same HIST");

    TH1* stacksum = common::summedStackHisto(stack.get());
    TH1* uncBand = nullptr, *uncBandPlot = nullptr;
    if(drawUncBand_){
        uncBand = dynamic_cast<TH1*>(allttbar->Clone("uncBand"));//common::summedStackHisto(stack.get());
        getSignalUncertaintyBand(uncBand, Channel);
        uncBand->SetFillStyle(3354);
        uncBand->SetFillColor(kBlack);
        uncBand->SetLineColor(kBlack);
        gStyle->SetHatchesLineWidth(1);
        gStyle->SetHatchesSpacing(0.8);
        uncBand->SetMarkerStyle(0);

        uncBandPlot = dynamic_cast<TH1*> (uncBand->Clone("uncBandPlot"));
        for (int i=0; i<=stacksum->GetNbinsX(); i++){
            uncBandPlot->SetBinContent(i, stacksum->GetBinContent(i));
        }
        uncBandPlot->Draw("same,e2");
        leg->AddEntry(uncBand, "Uncertainty", "f");
        if(leg2)leg2->AddEntry(uncBand, "Uncertainty", "f");
    }

    gPad->RedrawAxis();
    TExec *setex1 { new TExec("setex1","gStyle->SetErrorX(0.5)") };//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    if(drawUncBand_)uncBandPlot->Draw("same,e2");
    TExec *setex2 { new TExec("setex2","gStyle->SetErrorX(0.)") };
    setex2->Draw();
    drawhists[0]->Draw("same,e1");

    DrawCMSLabels(0, 8);
    DrawDecayChLabel(channelLabel[channelType]);
    if(name.Contains("JetMult")) {
        TString legtit = "";
        if(name.Contains("pt60")) legtit += "p_{T}^{jet}> 60 GeV";
        else if(name.Contains("pt100")) legtit += "p_{T}^{jet}> 100 GeV";
        else legtit += "p_{T}^{jet}> 30 GeV";
        leg->SetHeader(legtit);
    }
//    leg->Draw("SAME");
    if(leg1) leg1->Draw("SAME");
    if(leg2) leg2->Draw("SAME");

    if(drawPlotRatio){
        double yminCP_ = 0.49, ymaxCP_ = 1.51;
        yRangeControlPlotRatio(yminCP_, ymaxCP_);
        common::drawRatio(drawhists[0], stacksum, uncBand, yminCP_, ymaxCP_, doFit_);
    }

    // Create Directory for Output Plots 
    TString outdir = common::assignFolder(outpathPlots, Channel, Systematic);
    c->Print(outdir.Copy()+name+".eps");

    // Get the ratio plot from the canvas
    TPad* tmpPad = dynamic_cast<TPad*>(c->GetPrimitive("rPad"));
    TH1* ratio = nullptr;
    if(tmpPad) ratio = dynamic_cast<TH1*>(tmpPad->GetPrimitive("ratio"));

    //save Canvas AND sources in a root file
    TFile out_root(outdir.Copy()+name+"_source.root", "RECREATE");
    drawhists[0]->Write(name+"_data");
    sumttbar->Write(name+"_signalmc");
    allttbar->Write(name+"_allttbar");
    stacksum->SetName(name);
    stacksum->Write(name+"_allmc");
    if(ratio && ratio->GetEntries())ratio->Write("ratio");
    c->Write(name + "_canvas");
    out_root.Close();

    c->Clear();

    for (TH1* h : drawhists) delete h;
}

void Plotter::setStyle(TH1 *hist, unsigned int i, bool isControlPlot)
{
    hist->SetFillColor(colors[i]);
    hist->SetLineColor(colors[i]);
    hist->SetLineWidth(1);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.08);
    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetXaxis()->SetLabelOffset(0.007);
    hist->GetYaxis()->SetLabelOffset(0.007);

    if(legends.at(i) == "Data"){
        hist->SetFillColor(0);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.2);
        hist->SetLineWidth(2);
        if ((name.Contains("pT") || name.Contains("Mass")) && (!name.Contains("1st") && !name.Contains("Rapidity") && !name.Contains("Phi"))) {
            hist->GetXaxis()->SetTitle(XAxis+" #left[GeV#right]");
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}"+" #left[GeV^{-1}#right]"); 
        } else {
            hist->GetXaxis()->SetTitle(XAxis);
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+XAxis+"}");
        }
        if (isControlPlot) hist->GetYaxis()->SetTitle(YAxis);
    }
}


void Plotter::PlotXSec(TString Channel){

    TH1::AddDirectory(kFALSE);

    std::vector<TString> vec_systematic {"PDF_", "HAD_", "MATCH_", "MASS_", "SCALE_", "BTAG_", "BTAG_LJET_", "KIN_", "LEPT_", "TRIG_", "BG_", "DY_", "PU_", "JER_", "JES_"};//For the time being uintil all systematics are finalished
    std::vector<TString> vec_channel {"ee","mumu","emu","combined"};

    double BR_Error = 0.015;
    double Lumi_Error = 0.026;

    double InclusiveXsectionPlot[4] = {0.}, InclusiveXsectionStatErrorPlot[4] = {0.}, InclusiveXsectionSysErrorPlot[4] = {0.}, InclusiveXsectionTotalErrorPlot[4] = {0.};
    for (int j=0; j<(int)vec_channel.size(); j++){
        TString outdir = common::assignFolder(outpathPlots, vec_channel.at(j), TString("FinalResults"));
        ifstream SysResultsList("Plots/Nominal/"+vec_channel.at(j)+"/InclusiveXSec.txt");
        TString DUMMY;
        SysResultsList>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>InclusiveXsectionPlot[j]>>DUMMY>>InclusiveXsectionStatErrorPlot[j];
        SysResultsList.close();

        std::ofstream OutputFile(outdir.Copy()+"InclusiveXSecResultLateX.txt", std::ofstream::trunc);
        OutputFile<<"Inclusive XSection Numerical Results for channel "<<vec_channel.at(j)<<std::endl;

        double syst_square_for_channel=0.0;
        for (int i=0; i<(int) vec_systematic.size(); i++){

            ifstream SysUP, SysDOWN;

            if(vec_systematic.at(i) != "HAD_"){
                SysUP.open("Plots/"+vec_systematic.at(i)+"UP/"+vec_channel.at(j)+"/InclusiveXSec.txt");
                SysDOWN.open("Plots/"+vec_systematic.at(i)+"DOWN/"+vec_channel.at(j)+"/InclusiveXSec.txt");
                if(!SysUP.is_open() || !SysDOWN.is_open()) continue;
            } else{
                SysUP.open("Plots/MCATNLO/"+vec_channel.at(j)+"/InclusiveXSec.txt");
                SysDOWN.open("Plots/POWHEG/"+vec_channel.at(j)+"/InclusiveXSec.txt");
            }
            double VarUp = 0, VarDown = 0, StatErrUp = 0, StatErrDown = 0;

            SysUP>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>VarUp>>DUMMY>>StatErrUp;
            SysDOWN>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>DUMMY>>VarDown>>DUMMY>>StatErrDown;
            SysUP.close();
            SysDOWN.close();

            //systematic error in %
            double sys_err=(TMath::Abs(InclusiveXsectionPlot[j]-VarUp)+TMath::Abs(InclusiveXsectionPlot[j]-VarDown))*0.5/InclusiveXsectionPlot[j];
            syst_square_for_channel+=sys_err*sys_err;
            OutputFile<<vec_systematic.at(i)<<" (%): "<<setprecision(3)<<sys_err*100<<std::endl;
        }
        OutputFile<<"BranchingRatio (%): "<<setprecision(3)<<BR_Error*100<<std::endl;
        OutputFile<<"Luminosity (%): "<<setprecision(3)<<Lumi_Error*100<<std::endl;

        InclusiveXsectionSysErrorPlot[j]=TMath::Sqrt(syst_square_for_channel + BR_Error * BR_Error + Lumi_Error * Lumi_Error);

        InclusiveXsectionTotalErrorPlot[j] = sqrt(InclusiveXsectionStatErrorPlot[j]*InclusiveXsectionStatErrorPlot[j] +
                                                  InclusiveXsectionPlot[j]*InclusiveXsectionSysErrorPlot[j]*InclusiveXsectionPlot[j]*InclusiveXsectionSysErrorPlot[j]
                                                 );

        OutputFile<<"\n\n*******************************************************************************\n\n";
        OutputFile<<" InclXsec[pb]     Stat.[pb]    Syst.[pb]   Total[pb]"<<std::endl;
        OutputFile<<setprecision(6)<<InclusiveXsectionPlot[j]<<" +- "<<setprecision(3)<<InclusiveXsectionStatErrorPlot[j]<<" +- "<<setprecision(4)<<InclusiveXsectionSysErrorPlot[j]*InclusiveXsectionPlot[j]<<" +- "<<setprecision(4)<<InclusiveXsectionTotalErrorPlot[j]<<std::endl;
        OutputFile.close();
    }

    // measured results with statistical error
    Double_t mx[]   = {      0.50,       1.50,       2.50,       3.50};
    Double_t mexl[] = {      0.00,       0.00,       0.00,       0.00};
    Double_t mexh[] = {      0.00,       0.00,       0.00,       0.00};

    TGraphAsymmErrors *mplot = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh,InclusiveXsectionStatErrorPlot, InclusiveXsectionStatErrorPlot);
    mplot->SetMarkerStyle(20);
    mplot->GetYaxis()->SetNoExponent(kTRUE);
    mplot->SetMarkerColor(kBlack);
    mplot->SetMarkerSize(1.5);
    mplot->SetLineColor(kBlack);

    TGraphAsymmErrors *mplotwithsys = new TGraphAsymmErrors(4, mx, InclusiveXsectionPlot, mexl, mexh,InclusiveXsectionTotalErrorPlot, InclusiveXsectionTotalErrorPlot);
    mplotwithsys->SetMarkerStyle(20);
    mplotwithsys->SetMarkerColor(kBlack);
    mplotwithsys->SetMarkerSize(1.5);
    mplotwithsys->SetLineColor(kBlack);

    // kidonakis
    Double_t kidonmean = 234;
    Double_t kidonx[]   = {    -0.5,     0.5,   1.5,     2.5,     3.5,     4.5};
    Double_t kidony[]   = {kidonmean,kidonmean,kidonmean,kidonmean,kidonmean,kidonmean};
    Double_t kidonexl[] = {      .4,    .4,      .5,      .5,      .5,      .5};
    Double_t kidonexh[] = {      .5,    .5,      .5,      .5,      .4,      .4};
    Double_t kidoneyl[] = {    15.6,    15.6,    15.6,  15.6,    15.6,    15.6};
    Double_t kidoneyh[] = {    13.9,    13.9,    13.9,  13.9,    13.9,    13.9};

    TGraphAsymmErrors *kidonplot = new TGraphAsymmErrors(6, kidonx, kidony, kidonexl, kidonexh, kidoneyl, kidoneyh);
    kidonplot->SetLineColor(kGreen+1);
    kidonplot->SetLineWidth(4);
    kidonplot->SetFillColor(kGreen+1);
    kidonplot->SetFillStyle(3004);

    // Full NNLO
    Double_t nnlomean = 245.794;
    Double_t errDown = 10.656;
    Double_t errUp   = 8.652;
    Double_t nnlox[]   = {    -0.5,     0.5,   1.5,     2.5,     3.5,     4.5};
    Double_t nnloy[]   = {nnlomean,nnlomean,nnlomean,nnlomean,nnlomean,nnlomean};
    Double_t nnloexl[] = {      .4,    .4,      .5,      .5,      .5,      .5};
    Double_t nnloexh[] = {      .5,    .5,      .5,      .5,      .4,      .4};
    Double_t nnloeyl[] = { errDown, errDown, errDown, errDown, errDown, errDown};
    Double_t nnloeyh[] = {   errUp,   errUp,   errUp,   errUp,   errUp,   errUp};

    TGraphAsymmErrors *nnloplot = new TGraphAsymmErrors(6, nnlox, nnloy, nnloexl, nnloexh, nnloeyl, nnloeyh);
    nnloplot->SetLineColor(kGreen+1);
    nnloplot->SetLineWidth(4);
    nnloplot->SetFillColor(kGreen+1);
    nnloplot->SetFillStyle(3004);

    // TopLHC working group, prescription for m_top = 172.5 GeV
    //   https://indico.cern.ch/getFile.py/access?contribId=4&sessionId=1&resId=0&materialId=slides&confId=280522

    Double_t toplhcwgmean = 252.89;
    Double_t toplhcwgDown = 15.313;
    Double_t toplhcwgUp   = 16.266;
    Double_t toplhcwgx[]   = {    -0.5,     0.5,   1.5,     2.5,     3.5,     4.5};
    Double_t toplhcwgy[]   = {toplhcwgmean,toplhcwgmean,toplhcwgmean,toplhcwgmean,toplhcwgmean,toplhcwgmean};
    Double_t toplhcwgexl[] = {      .4,    .4,      .5,      .5,      .5,      .5};
    Double_t toplhcwgexh[] = {      .5,    .5,      .5,      .5,      .4,      .4};
    Double_t toplhcwgeyl[] = { toplhcwgDown, toplhcwgDown, toplhcwgDown, toplhcwgDown, toplhcwgDown, toplhcwgDown};
    Double_t toplhcwgeyh[] = {   toplhcwgUp,   toplhcwgUp,   toplhcwgUp,   toplhcwgUp,   toplhcwgUp,   toplhcwgUp};

    TGraphAsymmErrors *toplhcwgplot = new TGraphAsymmErrors(6, toplhcwgx, toplhcwgy, toplhcwgexl, toplhcwgexh, toplhcwgeyl, toplhcwgeyh);
    toplhcwgplot->SetLineColor(kGreen+1);
    toplhcwgplot->SetLineWidth(4);
    toplhcwgplot->SetFillColor(kGreen+1);
    toplhcwgplot->SetFillStyle(3004);

    // mcfm
    Double_t mcfmmean = 225.197;
    Double_t mcfmx[]   = {      -0.5,     0.5,     1.5,     2.5,     3.5,     4.5};
    Double_t mcfmy[]   = {mcfmmean,mcfmmean,mcfmmean,mcfmmean,mcfmmean,mcfmmean};
    Double_t mcfmexl[] = {        .4,      .4,      .5,      .5,      .5,      .5};
    Double_t mcfmexh[] = {        .5,      .5,      .5,      .5,      .4,      .4};
    Double_t mcfmeyl[] = {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
    Double_t mcfmeyh[] = {   0.0, 0.0,  0.0,  0.0,  0.0,  0.0};

    TGraphAsymmErrors *mcfmplot = new TGraphAsymmErrors(6, mcfmx, mcfmy, mcfmexl, mcfmexh, mcfmeyl, mcfmeyh);
    mcfmplot->SetLineColor(kBlue+1);
    mcfmplot->SetLineWidth(4);
    mcfmplot->SetFillColor(kBlue+1);
    mcfmplot->SetFillStyle(3005);

    TH1F* framehist = new TH1F("framehist","",4,0.,4.);
    framehist->SetMinimum(100);
    framehist->SetMaximum(380);
    framehist->GetXaxis()->SetTickLength(0);
    framehist->GetXaxis()->SetBinLabel(1,"");
    framehist->GetXaxis()->SetBinLabel(2,"");
    framehist->GetXaxis()->SetBinLabel(3,"");
    framehist->GetYaxis()->SetTitle("#sigma [pb]");
    framehist->GetYaxis()->CenterTitle(kTRUE);
    framehist->GetYaxis()->SetNoExponent(kTRUE);

    TPaveText* box1 = new TPaveText(0.25,0.20,0.33,0.30,"NDC");
    box1->SetFillColor(10);
    box1->SetTextSize(0.04);
    box1->AddText("ee");

    TPaveText* box2 = new TPaveText(0.44,0.20,0.52,0.30,"NDC");
    box2->SetFillColor(10);
    box2->SetTextSize(0.04);
    box2->AddText("#mu#mu");

    TPaveText* box3 = new TPaveText(0.62,0.20,0.72,0.30,"NDC");
    box3->SetFillColor(10);
    box3->SetTextSize(0.04);
    box3->AddText("e#mu");

    TPaveText* box4 = new TPaveText(0.82,0.20,0.90,0.30,"NDC");
    box4->SetFillColor(10);
    box4->SetTextSize(0.04);
    box4->AddText("combined");

    TLegend* leg =  new TLegend( 0.56, 0.65, 0.89, 0.85);
    leg->SetBorderSize( 0 );
    leg->SetFillColor( 0 );
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->AddEntry( mplot,       "Measurements",            "p"  );
    leg->AddEntry( mcfmplot, "MCFM #otimes CTQE66M", "lf" );
//    leg->AddEntry( kidonplot,    "Kidonakis #otimes MSTW2008 NNLO",     "lf" );
//    leg->AddEntry( nnloplot,    "NNLO #otimes MSTW2008 NNLO",     "lf" );
    leg->AddEntry( toplhcwgplot,    "TOP LHC WG",     "lf" );

    TCanvas* c = new TCanvas("plot", "plot", 1200, 800);
    framehist->Draw();
    mcfmplot->Draw("C,2,SAME");
    //kidonplot->Draw("C,2,SAME");
    //nnloplot->Draw("C,2,SAME");
    toplhcwgplot->Draw("C,2,SAME");
    gStyle->SetEndErrorSize(8);
    mplot->Draw("p,SAME");
    mplotwithsys->Draw("p,SAME,Z");
    leg ->Draw("SAME");

    box1->Draw("SAME");
    box2->Draw("SAME");
    box3->Draw("SAME");
    box4->Draw("SAME");

    TString outdir = common::assignFolder(outpathPlots, Channel, TString("FinalResults"));
    c->Print(outdir.Copy()+"InclusiveXSec.eps");
    c->Print(outdir.Copy()+"InclusiveXSec.C");
    c->Clear();
    delete c;

}


void Plotter::MakeTable(TString Channel, TString Systematic){

    TH1D *numhists5[hists.size()];
    TH1D *numhists6[hists.size()];
    TH1D *numhists7[hists.size()];
    TH1D *numhists8[hists.size()];
    TH1D *numhists9[hists.size()];

    for(unsigned int i=0; i<dataset.size(); i++){

        TH1D *temp_hist5 = fileReader->GetClone<TH1D>(dataset[i], "events_weighted_step4");//TH1D *temp_hist5 = fileReader->GetClone<TH1D>(dataset[i], "step5");
        TH1D *temp_hist6 = fileReader->GetClone<TH1D>(dataset[i], "events_weighted_step5");//TH1D *temp_hist6 = fileReader->GetClone<TH1D>(dataset[i], "step6");
        TH1D *temp_hist7 = fileReader->GetClone<TH1D>(dataset[i], "events_weighted_step6");//TH1D *temp_hist7 = fileReader->GetClone<TH1D>(dataset[i], "step7");
        TH1D *temp_hist8 = fileReader->GetClone<TH1D>(dataset[i], "events_weighted_step7");//TH1D *temp_hist8 = fileReader->GetClone<TH1D>(dataset[i], "step8");
        TH1D *temp_hist9 = fileReader->GetClone<TH1D>(dataset[i], "events_weighted_step8");//TH1D *temp_hist9 = fileReader->GetClone<TH1D>(dataset[i], "step9");

        double LumiWeight = CalcLumiWeight(dataset.at(i));
        ApplyFlatWeights(temp_hist5, LumiWeight);
        ApplyFlatWeights(temp_hist6, LumiWeight);
        ApplyFlatWeights(temp_hist7, LumiWeight);
        ApplyFlatWeights(temp_hist8, LumiWeight);
        ApplyFlatWeights(temp_hist9, LumiWeight);

        numhists5[i]=temp_hist5;
        numhists6[i]=temp_hist6;
        numhists7[i]=temp_hist7;
        numhists8[i]=temp_hist8;
        numhists9[i]=temp_hist9;

    }

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        if(legends.at(i) == DYEntry){
            //numhists5[i]->Scale(DYScale[channelType]);//DYscale not applied in step5 and 6?
            //numhists6[i]->Scale(DYScale[channelType]);
            numhists7[i]->Scale(DYScale.at(channelType));
            numhists8[i]->Scale(DYScale.at(channelType));
            numhists9[i]->Scale(DYScale.at(channelType));
        }
    }

    ////////////////////////////Make output for tables
    double tmp_num5 = 0;
    double tmp_num6 = 0;
    double tmp_num7 = 0;
    double tmp_num8 = 0;
    double tmp_num9 = 0;

    TString outdir = common::assignFolder(outpathPlots, Channel, Systematic);
    ofstream EventFile5; EventFile5.open(outdir.Copy()+"Events5.txt");
    ofstream EventFile6; EventFile6.open(outdir.Copy()+"Events6.txt");
    ofstream EventFile7; EventFile7.open(outdir.Copy()+"Events7.txt");
    ofstream EventFile8; EventFile8.open(outdir.Copy()+"Events8.txt");
    ofstream EventFile9; EventFile9.open(outdir.Copy()+"Events9.txt");
    
    double bg_num5 = 0;
    double bg_num6 = 0;
    double bg_num7 = 0;
    double bg_num8 = 0;
    double bg_num9 = 0;

    for(unsigned int i=0; i<hists.size() ; i++){
        tmp_num5+=numhists5[i]->Integral();
        tmp_num6+=numhists6[i]->Integral();
        tmp_num7+=numhists7[i]->Integral();
        tmp_num8+=numhists8[i]->Integral();
        tmp_num9+=numhists9[i]->Integral();

        if(i==(hists.size()-1)){
            EventFile5<<legends.at(i)<<": "<<tmp_num5<<std::endl;
            EventFile6<<legends.at(i)<<": "<<tmp_num6<<std::endl;
            EventFile7<<legends.at(i)<<": "<<tmp_num7<<std::endl;
            EventFile8<<legends.at(i)<<": "<<tmp_num8<<std::endl;
            EventFile9<<legends.at(i)<<": "<<tmp_num9<<std::endl;
            bg_num5+=tmp_num5;
            bg_num6+=tmp_num6;
            bg_num7+=tmp_num7;
            bg_num8+=tmp_num8;
            bg_num9+=tmp_num9;
            tmp_num5=0;
            tmp_num6=0;
            tmp_num7=0;
            tmp_num8=0;
            tmp_num9=0;
        }
        else if(legends.at(i)!=legends.at(i+1)){
            EventFile5<<legends.at(i)<<": "<<tmp_num5<<std::endl;
            EventFile6<<legends.at(i)<<": "<<tmp_num6<<std::endl;
            EventFile7<<legends.at(i)<<": "<<tmp_num7<<std::endl;
            EventFile8<<legends.at(i)<<": "<<tmp_num8<<std::endl;
            EventFile9<<legends.at(i)<<": "<<tmp_num9<<std::endl;
            if(legends.at(i)!="Data"){
                bg_num5+=tmp_num5;
                bg_num6+=tmp_num6;
                bg_num7+=tmp_num7;
                bg_num8+=tmp_num8;
                bg_num9+=tmp_num9;
            }
            tmp_num5=0;
            tmp_num6=0;
            tmp_num7=0;
            tmp_num8=0;
            tmp_num9=0;
        }
    }
    EventFile5<<"Total background: "<<bg_num5<<std::endl;
    EventFile5.close();
    EventFile6<<"Total background: "<<bg_num6<<std::endl;
    EventFile6.close();
    EventFile7<<"Total background: "<<bg_num7<<std::endl;
    EventFile7.close();
    EventFile8<<"Total background: "<<bg_num8<<std::endl;
    EventFile8.close();
    EventFile9<<"Total background: "<<bg_num9<<std::endl;
    EventFile9.close();
    std::cout<<"\nEvent yields saved in "<<outdir<<std::endl;
}

double Plotter::CalcXSec(std::vector<TString> datasetVec, double InclusiveXsectionVec[4],double InclusiveXsectionStatErrorVec[4], TString Systematic, TString Shift){

//    double BranchingFraction[4]={0.01582, 0.01573, 0.03155, 0.06310};//[ee, mumu, emu, combined] including tau into electron/muon
    double BranchingFraction[4]={0.01166, 0.01166, 0.02332, 0.04666};//[ee, mumu, emu, combined] not including tau

    double NrOfEvts_VisGen_afterSelection_noweight = 0, NrOfEvts_VisGen_afterSelection = 0;
    double NrOfEvts_afterSelection_noweight = 0, NrOfEvts_afterSelection = 0;
    double NrOfEvts_Gen_afterRecoSelection_noweight = 0, NrOfEvts_Gen_afterRecoSelection = 0;
    double NrOfEvts = 0;

    TH1D *numhists[hists.size()];
    double numbers[5]={0., 0., 0., 0., 0.};//[0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]background(non-ttbar)
    double error_numbers[5]={0., 0., 0., 0., 0.};//Square of error: [0]=data, [1]=Signal, [2]Signal(only lumi & PU weights), [3]ttbar background, [4]background(non-ttbar)
//     double TTbarBGnum =0;

    if (Systematic.Contains("UP"))  { Shift="Up";}
    if (Systematic.Contains("DOWN")){ Shift="Down";}

    for(unsigned int i=0; i<datasetVec.size(); i++){
        TH1D *hist = fileReader->GetClone<TH1D>(datasetVec[i], "events_weighted_step8");
        ApplyFlatWeights(hist, CalcLumiWeight(datasetVec.at(i)));
        numhists[i]=hist;
    }

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg 
        if(legends.at(i) == "Data"){
            numbers[0]+=numhists[i]->Integral();
            error_numbers[0]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
        }
        else if(legends.at(i) == "t#bar{t} Signal"){
            numbers[1]+=numhists[i]->Integral();
            error_numbers[1]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C

            TH1D *GenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll");
            TH1D *GenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_noweight");
            TH1D *VisGenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "VisGenAll");
            TH1D *VisGenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "VisGenAll_noweight");
            TH1D *RecoGenPlot = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_RecoCuts");
            TH1D *RecoGenPlot_noweight = fileReader->GetClone<TH1D>(datasetVec.at(i), "GenAll_RecoCuts_noweight");
            TH1 *h_NrOfEvts = fileReader->GetClone<TH1>(datasetVec.at(i), "weightedEvents");

            double LumiWeight = CalcLumiWeight(datasetVec.at(i));
            ApplyFlatWeights(GenPlot, LumiWeight);
            ApplyFlatWeights(GenPlot_noweight, LumiWeight);
            ApplyFlatWeights(VisGenPlot, LumiWeight);
            ApplyFlatWeights(VisGenPlot_noweight, LumiWeight);
            ApplyFlatWeights(RecoGenPlot, LumiWeight);
            ApplyFlatWeights(RecoGenPlot_noweight, LumiWeight);
            ApplyFlatWeights(h_NrOfEvts, LumiWeight);

            NrOfEvts += h_NrOfEvts->GetBinContent(1);
            NrOfEvts_afterSelection += GenPlot->Integral();
            NrOfEvts_afterSelection_noweight += GenPlot_noweight->Integral();
            NrOfEvts_VisGen_afterSelection += VisGenPlot->Integral();
            NrOfEvts_Gen_afterRecoSelection += RecoGenPlot->Integral();
            NrOfEvts_Gen_afterRecoSelection_noweight += RecoGenPlot_noweight->Integral();
            NrOfEvts_VisGen_afterSelection_noweight += VisGenPlot_noweight->Integral();

            numbers[2]+=GenPlot->Integral();
            error_numbers[2]+=GenPlot->GetBinError(18) * GenPlot->GetBinError(18); //This bin selection is hardcoded please change it if changes when filling in Analysis.C

        }
        else if(legends.at(i) == "t#bar{t} Other"){
            numbers[3]+=numhists[i]->Integral();
            error_numbers[3]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
        }
        else {
            if((legends.at(i) == DYEntry)){
                numhists[i]->Scale(DYScale.at(channelType));
            }
            if((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Up"){
                numhists[i]->Scale(1.3);
            }
            if((legends.at(i) == DYEntry) && Systematic.Contains("DY_") && Shift == "Down"){
                numhists[i]->Scale(0.7);
            }
            if(Systematic.Contains("BG_") && Shift=="Up" && legends.at(i)!= "t#bar{t} Other" && legends.at(i) != DYEntry){
                numhists[i]->Scale(1.3);
            }
            if(Systematic.Contains("BG_") && Shift=="Down" && legends.at(i)!= "t#bar{t} Other" && legends.at(i) != DYEntry){
                numhists[i]->Scale(0.7);
            }
            numbers[4]+=numhists[i]->Integral();
            error_numbers[4]+=numhists[i]->GetBinError(2) * numhists[i]->GetBinError(2); //This bin selection is hardcoded please change it if changes when filling in Analysis.C
        }
    }

    ////////////////////////////Make output for tables

    double tmp_num = 0;

    ofstream EventFile, XSecFile;
    TString outdir = common::assignFolder(outpathPlots, subfolderChannel.Copy().Remove(0,1), Systematic);
    EventFile.open(outdir.Copy()+"Events.txt");
    XSecFile.open(outdir.Copy()+"InclusiveXSec.txt");

    double bg_num = 0;
    for(unsigned int i=0; i<hists.size() ; i++){
        tmp_num+=numhists[i]->Integral();

        if(i==(hists.size()-1)){
            EventFile<<legends.at(i)<<": "<<tmp_num<<std::endl;
            bg_num+=tmp_num;
            tmp_num=0;
        }
        else if(legends.at(i)!=legends.at(i+1)){
            EventFile<<legends.at(i)<<": "<<tmp_num<<std::endl;
            if(legends.at(i)!="Data")bg_num+=tmp_num;
            tmp_num=0;
        }
    }
    EventFile<<"Total MCs: "<<bg_num<<std::endl;
    EventFile<<"\nDataEvents= "<<numbers[0]<<std::endl;
    EventFile<<"SignalReco= "<<numbers[1]<<std::endl;
    EventFile<<"SignalGen = "<<numbers[2]<<std::endl;
    EventFile<<"ttbar bags= "<<numbers[3]<<std::endl;
    EventFile<<"All Backgd= "<<numbers[4]<<std::endl;
    EventFile<<"Efficiency= "<<(numbers[1]/numbers[2])<<std::endl;
    EventFile<<"BrancRatio= "<<BranchingFraction[channelType]<<std::endl;
    EventFile<<"Total Gen Events (no weights)= "<<NrOfEvts<<std::endl;
    EventFile<<"Gen Events after Selection (no weights)= "<<NrOfEvts_afterSelection_noweight<<std::endl;
    EventFile<<"Visible Gen Events after Selection (no weights)= "<<NrOfEvts_VisGen_afterSelection_noweight<<std::endl;
    EventFile<<"Acceptance= "<<NrOfEvts_afterSelection_noweight/NrOfEvts<<std::endl;
    EventFile<<"Visible Acceptance= "<<NrOfEvts_VisGen_afterSelection_noweight/NrOfEvts<<std::endl;
    EventFile<<"------------------------------------------------------------------------------"<<std::endl;
    EventFile<<"Efficiency and acceptance definitions as proposed by the TopXSection conveners\n"<<std::endl;
    EventFile<<"N_rec = "<<numbers[1]<<std::endl;
    EventFile<<"N_gen (with cuts at parton level) = "<<NrOfEvts_VisGen_afterSelection<<std::endl;
    EventFile<<"N_gen (with cuts at parton level, no weights) = "<<NrOfEvts_VisGen_afterSelection_noweight<<std::endl;
    EventFile<<"N_gen (with cuts at reco level) = "<<NrOfEvts_Gen_afterRecoSelection<<std::endl;
    EventFile<<"N_gen (with cuts at reco level, no weights) = "<<NrOfEvts_Gen_afterRecoSelection_noweight<<std::endl;
    EventFile<<"N_gen = "<<numbers[2]<<std::endl;
    EventFile<<"\nEfficiency = N_rec / N_gen (with cuts at parton level) = "<<numbers[1]/NrOfEvts_VisGen_afterSelection<<std::endl;
    EventFile<<"Efficiency = N_rec / N_gen (with cuts at parton level && noweights) = "<<numbers[1]/NrOfEvts_VisGen_afterSelection_noweight<<std::endl;
    EventFile<<"\nEfficiency' = N_rec / N_gen (with cuts at reco level) = "<<numbers[1]/NrOfEvts_Gen_afterRecoSelection<<std::endl;
    EventFile<<"Efficiency' = N_rec / N_gen (with cuts at reco level && noweights) = "<<numbers[1]/NrOfEvts_Gen_afterRecoSelection_noweight<<std::endl;
    EventFile<<"\nAcceptance = N_gen (with cuts at parton level) / N_gen = "<<NrOfEvts_VisGen_afterSelection/numbers[2]<<std::endl;
    EventFile<<"Acceptance = N_gen (with cuts at parton level && noweights) / N_gen = "<<NrOfEvts_VisGen_afterSelection_noweight/numbers[2]<<std::endl;
    EventFile<<"Eff * Acc = "<<numbers[1] / numbers[2]<<std::endl;
    EventFile<<"------------------------------------------------------------------------------"<<std::endl;


    double xsec = ( (numbers[0]-numbers[4]) * (numbers[1]/(numbers[1]+numbers[3])) ) / ( (numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);
    double xsecstaterror = TMath::Sqrt(error_numbers[0]) * (numbers[1]/(numbers[1]+numbers[3])) / ( (numbers[1]/numbers[2])*BranchingFraction[channelType]*lumi);

    if(channelType!=3){
        InclusiveXsectionVec[channelType] = xsec;
        InclusiveXsectionStatErrorVec[channelType] = xsecstaterror;
    }
    else{
        TString eefilename="Plots/"+Systematic+"/ee/InclusiveXSec.txt";
        TString mumufilename="Plots/"+Systematic+"/mumu/InclusiveXSec.txt";
        TString emufilename="Plots/"+Systematic+"/emu/InclusiveXSec.txt";

        //check the existence of the file
        if( gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename)){
            std::cout<<"WARNING (in CalcXSec)!!"<<std::endl;
            std::cout<<"One of the input files you use for the combined XSection measurement doesn't exist!!\nExiting!!"<<std::endl;
            exit(888);
        }

        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);

        TString Dummy="";

        ResultsEE>>Dummy>>Dummy>>Dummy>>Dummy>>Dummy>>InclusiveXsectionVec[0]>>Dummy>>InclusiveXsectionStatErrorVec[0];
        ResultsMuMu>>Dummy>>Dummy>>Dummy>>Dummy>>Dummy>>InclusiveXsectionVec[1]>>Dummy>>InclusiveXsectionStatErrorVec[1];
        ResultsEMu>>Dummy>>Dummy>>Dummy>>Dummy>>Dummy>>InclusiveXsectionVec[2]>>Dummy>>InclusiveXsectionStatErrorVec[2];

        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        InclusiveXsectionVec[channelType] =( InclusiveXsectionVec[0]/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0])
                                            +InclusiveXsectionVec[1]/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1])
                                            +InclusiveXsectionVec[2]/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2]) 
                                            )/
                                            ( 1/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0])
                                            + 1/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1])
                                            + 1/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2])
                                            );

        InclusiveXsectionStatErrorVec[channelType] =1/(TMath::Sqrt(
                                            (1/(InclusiveXsectionStatErrorVec[0]*InclusiveXsectionStatErrorVec[0]))
                                            +(1/(InclusiveXsectionStatErrorVec[1]*InclusiveXsectionStatErrorVec[1]))
                                            +(1/(InclusiveXsectionStatErrorVec[2]*InclusiveXsectionStatErrorVec[2]))
                                                                    ));
    }

    EventFile<<"XSection  = "<<InclusiveXsectionVec[channelType]<<std::endl;
    EventFile<<"XSecStaErr= "<<InclusiveXsectionStatErrorVec[channelType]<<std::endl;
    EventFile.close();
    XSecFile<<"Systematic: "<<Systematic<<" Channel: "<<subfolderChannel<<" InclXSection: "<<InclusiveXsectionVec[channelType]<<" AbsStatError: "<<InclusiveXsectionStatErrorVec[channelType]<<std::endl;
    XSecFile.close();
    std::cout<<"\nInclusive XSection information saved in: "<<outdir<<std::endl;
    return xsec;
}

int Plotter::CalcDiffXSec(TString Channel, TString Systematic){

    double Xbins[XAxisbins.size()];
    double binWidth[XAxisbinCenters.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
    double GenSignalSum[XAxisbinCenters.size()];

    TString ftemp = "preunfolded/"+Systematic+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
    if(gSystem->AccessPathName(ftemp)){
        std::cout<<"WARNING (in CalcDiffXSec)!!"<<std::endl;
        std::cout<<"File: "<<ftemp<<" doesn't exist"<<std::endl;
        std::cout<<"One of the input files you use for the combined XSection measurement doesn't exist!!"<<std::endl;
        return -1;
    }
    TH1D* theDataHist =  fileReader->GetClone<TH1D>(ftemp, "aDataHist");
    TH1D* theBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aBgrHist");
    TH1D* theTtBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aTtBgrHist");
    TH1D* theRecHist =  fileReader->GetClone<TH1D>(ftemp, "aRecHist");
    TH1D* theGenHist =  fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH2D* theRespHist =  fileReader->GetClone<TH2D>(ftemp, "aRespHist");

    std::unique_ptr<TH1> theGenHistRebinned { theGenHist->Rebin(bins,"aDataHist",Xbins) };
    for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
        GenSignalSum[bin] = theGenHistRebinned->GetBinContent(bin+1);
    }

    double GenDiffXSecVec[4][bins];
    double DiffXSecVec[4][bins];
    double DiffXSecStatErrorVec[4][bins];

    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}

    // DAVID 
    if ( doUnfolding == true && channelType !=3)//do the unfolding only in the individual channels: ee, emu, mumu
    {
        // SVD Helper Class
        DilepSVDFunctions mySVDFunctions;
        mySVDFunctions.SetOutputPath(outpath);

        // Binning
        double* theBins = Xbins;
        int numberBins = bins;

        // Names and Labels
        TString channelLabelStr(channelLabel[channelType]);
        TString theChannelName = Channel;
        TString theParticleName = "";
        if ( name.Contains("Lepton")  ) theParticleName = "Leptons";
        if ( name.Contains("LLBar")   ) theParticleName = "LepPair";
        if ( name.Contains("Top")     ) theParticleName = "TopQuarks";
        if ( name.Contains("TTBar")   ) theParticleName = "TtBar";
        if ( name.Contains("BBBar")   ) theParticleName = "BBbar";
        if ( name.Contains("BJet")    ) theParticleName = "BJets"; 
        if ( name.Contains("JetMult")    ) theParticleName = "Jets";
        if ( name.Contains("ExtraJet") ) theParticleName = "ExtraJets";
        if ( name.Contains("DeltaRJet12") ) theParticleName = "DeltaR";
        TString theQuantityName = "";
        if ( name.Contains("pT")      ) theQuantityName = "Pt";
        if ( name.Contains("Eta")     ) theQuantityName = "Eta";
        if ( name.Contains("Rapidity")) theQuantityName = "Rapidity";
        if ( name.Contains("Mass")    ) theQuantityName = "Mass"; 
        if ( name.Contains("JetMult")    ) theQuantityName = "Mult";
        if ( name.Contains("ExtraJetEta") ) theQuantityName = "Eta";
        if ( name.Contains("ExtraJetpT") ) theQuantityName = "Pt";
        if ( name.Contains("DeltaRJet12") ) theQuantityName = "DeltaR";
        TString theSpecialPostfix = "";
        //if (name.Contains("Lead")) theSpecialPostfix = name;
        theSpecialPostfix = name;
        if ( specialComment.CompareTo("Standard") != 0 ) {
            //theSpecialPostfix = specialComment;
        }

        double totalDataEventsNom[1]  = {TopSVDFunctions::SVD_Integral1D(theDataHist, 0, false)}; 
        double totalBgrEventsNom[1]   = {TopSVDFunctions::SVD_Integral1D(theBgrHist, 0, false)};
        double totalTtBgrEventsNom[1]   = {TopSVDFunctions::SVD_Integral1D(theTtBgrHist, 0, false)};
        double totalRecEventsNom[1]   = {TopSVDFunctions::SVD_Integral1D(theRecHist, 0, false)};
        double totalGenEventsNom[1]  = {TopSVDFunctions::SVD_Integral1D(theGenHist, 0, true)}; 

        // UNFOLDING 
        // Retrieve a histogram with the unfolded quantities.
        // Note: The unfolded histogram has additional side bins!
        // Keep this in mind when accessing bin content via indices
        TH1D* unfoldedDistribution = NULL;
        TH1D* unfoldedDistributionNormalized = NULL;
        int numSystematics = 0;
        mySVDFunctions.SVD_DoUnfold(
                                    theDataHist, 
                                    theBgrHist, 
                                    theTtBgrHist, 
                                    theGenHist, 
                                    theRecHist, 
                                    theRespHist, 
                    totalDataEventsNom,
                    totalBgrEventsNom,
                    totalTtBgrEventsNom,
                    totalRecEventsNom,
                    totalGenEventsNom,
                                    theBins, numberBins,  
                                    unfoldedDistribution, 
                                    unfoldedDistributionNormalized,
                                    numSystematics,
                                    theChannelName, theParticleName, theQuantityName, theSpecialPostfix, "");


        // Make a vector from the result
        double UnfoldingResult[XAxisbinCenters.size()];
        double UnfoldingError[XAxisbinCenters.size()];
        for ( size_t i = 0; i < XAxisbinCenters.size() ; i++ ) {
            UnfoldingResult[i] = unfoldedDistributionNormalized->GetBinContent(i+2);//account for extra row in SVD unfolding
            UnfoldingError[i] = unfoldedDistributionNormalized->GetBinError(i+2); //absolute bin error
            //UnfoldingResult[i] = unfoldedDistribution->GetBinContent(i+2);//account for extra row in SVD unfolding
            //UnfoldingError[i] = unfoldedDistribution->GetBinError(i+2);
        }

        SignalEventswithWeight=0;
        // CROSS SECTION CALCULATION
        for (Int_t i=0; i<bins; ++i) {
            SignalEventswithWeight+=GenSignalSum[i];
        }

        for (Int_t i=0; i<bins; ++i) {
            //      if(Channel!="combined"){
            binWidth[i] = Xbins[i+1]-Xbins[i];
            DiffXSecVec[channelType][i] = UnfoldingResult[i]/(binWidth[i]);
            DiffXSecStatErrorVec[channelType][i] = UnfoldingError[i]/(binWidth[i]); // absolute statistical error
            GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)

            if(name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")){
                GenDiffXSecVec[channelType][i] = GenDiffXSecVec[channelType][i]/2.;
            }
        }
    }

    if (doUnfolding && channelType==3){//for 'combined' channel: do an statistical combination of the the 3 independent channels

        TString eefilename="UnfoldingResults/"+Systematic+"/ee/"+name+"Results.txt";
        TString emufilename="UnfoldingResults/"+Systematic+"/emu/"+name+"Results.txt";
        TString mumufilename="UnfoldingResults/"+Systematic+"/mumu/"+name+"Results.txt";

        //check the existence of the file
        if(gSystem->AccessPathName(eefilename) || gSystem->AccessPathName(emufilename) || gSystem->AccessPathName(mumufilename) ){
            std::cout<<"WARNING (in CalcDiffXSec)!!"<<std::endl;
            std::cout<<"One of the input files you use for the combined XSection measurement doesn't exist!!"<<std::endl;
            return -1;
        }

        ifstream ResultsEE(eefilename);
        ifstream ResultsEMu(emufilename);
        ifstream ResultsMuMu(mumufilename);
        double perChannelDiffXSecPlot[3][bins];      //perChannelDiffXSecPlot[channel][bin]
        double perChannelDiffXSecStatError[3][bins]; //perChannelDiffXSecStatError[channel][bin]
        double perChannelGenDiffXSec[3][bins];       //perChannelGenDiffXSec[channel][bin]
        TString Dummy="";
        for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
            ResultsEE>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[0][bin]>>Dummy>>perChannelDiffXSecStatError[0][bin]>>Dummy>>perChannelGenDiffXSec[0][bin];
            ResultsMuMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[1][bin]>>Dummy>>perChannelDiffXSecStatError[1][bin]>>Dummy>>perChannelGenDiffXSec[1][bin];
            ResultsEMu>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>perChannelDiffXSecPlot[2][bin]>>Dummy>>perChannelDiffXSecStatError[2][bin]>>Dummy>>perChannelGenDiffXSec[2][bin];
        }
        ResultsEE.close(); ResultsEMu.close(); ResultsMuMu.close();

        //for gen level distribution
        for (Int_t i=0; i<bins; ++i) {
            SignalEventswithWeight+=GenSignalSum[i];
        }

        //do the actual combined Diff.XSection calculation
        for(int i=0; i<bins; i++){
            binWidth[i] = Xbins[i+1]-Xbins[i];
            for(int j=0; j<3; j++){//check if any stat error is 0, in this case set their contribution to 0!!
                if(perChannelDiffXSecStatError[j][i] == 0){
                    perChannelDiffXSecStatError[j][i] = 1e100;
                    perChannelDiffXSecPlot[j][i] = 0;
                }
            }
            DiffXSecVec[channelType][i] =(perChannelDiffXSecPlot[0][i]/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i])
                                         +perChannelDiffXSecPlot[1][i]/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i])
                                         +perChannelDiffXSecPlot[2][i]/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i]))/
                                         ((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                          (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                          (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i])));
            DiffXSecStatErrorVec[channelType][i]=1/(TMath::Sqrt((1/(perChannelDiffXSecStatError[0][i]*perChannelDiffXSecStatError[0][i]))+
                                                                (1/(perChannelDiffXSecStatError[1][i]*perChannelDiffXSecStatError[1][i]))+
                                                                (1/(perChannelDiffXSecStatError[2][i]*perChannelDiffXSecStatError[2][i]))));
            GenDiffXSecVec[channelType][i] = (GenSignalSum[i]*topxsec)/(SignalEventswithWeight*binWidth[i]);//DIRTY (signal*topxsec)/(total events*binwidth)
        }
    }

    ofstream ResultsFile, ResultsLateX;

    gSystem->mkdir("UnfoldingResults/"+Systematic+"/"+Channel, true);

    string ResultsFilestring = outpathResults.Data();
    ResultsFilestring.append(subfolderSpecial.Data());
    ResultsFilestring.append("/");
    ResultsFilestring.append(Systematic);
    ResultsFilestring.append("/");
    ResultsFilestring.append(Channel);
    ResultsFilestring.append("/");
    ResultsFilestring.append(name);
    ResultsFilestring.append("Results.txt");
    ResultsFile.open(ResultsFilestring.c_str());


    string ResultsFilestringLatex = outpathPlots.Data();
    ResultsFilestringLatex.append(subfolderChannel.Data());
    ResultsFilestringLatex.append(subfolderSpecial.Data());
    ResultsFilestringLatex.append("/");
    ResultsFilestringLatex.append(name);
    ResultsFilestringLatex.append("ResultsLaTeX.txt");
    ResultsLateX.open(ResultsFilestringLatex.c_str());
    ResultsLateX<<"Bin Center & Bin & 1/#sigma d#sigma/dX & stat(%) & syst(%) & total(%)"<<std::endl;
    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
        ResultsFile<<"XAxisbinCenters[bin]: "<<XAxisbinCenters[bin]<<" bin: "<<Xbins[bin]<<" to "<<Xbins[bin+1]<<" DiffXsec: "<<DiffXSecVec[channelType][bin]<<" StatError: "<<DiffXSecStatErrorVec[channelType][bin]<<" GenDiffXsec: "<<GenDiffXSecVec[channelType][bin]<<std::endl;
    }
    ResultsFile.close();
    ResultsLateX.close();

    //clean up
    delete theDataHist;
    delete theBgrHist;
    delete theTtBgrHist;
    delete theRecHist;
    delete theGenHist;
    delete theRespHist;
    
    Plotter::PlotSingleDiffXSec(Channel, Systematic);

    return 1;
}

void Plotter::PlotDiffXSec(TString Channel, std::vector<TString>vec_systematic){

    setDataSet(Channel,"Nominal");
    if (!fillHisto()) return;

    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}

    TH1::AddDirectory(kFALSE);
    TGaxis::SetMaxDigits(2);

    double Xbins[XAxisbins.size()];
    for(unsigned int i = 0; i<XAxisbins.size();i++){Xbins[i]=XAxisbins[i];}
    double binCenters[XAxisbinCenters.size()];
    if ( drawSmoothMadgraph ) {
        for ( unsigned int i = 0; i<XAxisbinCenters.size();i++ ) {
            binCenters[i] = XAxisbinCenters[i];
        }
    } else {
        for ( unsigned int i = 0; i<XAxisbinCenters.size();i++ ) {
            binCenters[i] = XAxisbins[i] + (XAxisbins[i+1]-XAxisbins[i])/2;
        }
    }

    double DataSum[XAxisbinCenters.size()];
    double GenSignalSum[XAxisbinCenters.size()];
    double BGSum[XAxisbinCenters.size()];
    bool init = false;
    TH1 *varhists[hists.size()];
    TString newname = name;
    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
      newname.ReplaceAll("Hyp",3,"",0);
    }

    TString ftemp = "preunfolded/Nominal/"+Channel+"/"+name+"_UnfoldingHistos.root";

//     theDataHist =  fileReader->GetClone<TH1D>(ftemp, "aDataHist");
//     theBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aBgrHist");
//     theTtBgrHist =  fileReader->GetClone<TH1D>(ftemp, "aTtBgrHist");
//     RecoPlot =  fileReader->GetClone<TH1D>(ftemp, "aRecHist");
    TH1 *GenPlotTheory =  fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH2 *genReco2d =  fileReader->GetClone<TH2D>(ftemp, "aRespHist");


    for (unsigned int i =0; i<hists.size(); i++){
      varhists[i]=hists[i].Rebin(bins,"varhists",Xbins);
      setStyle(varhists[i], i);
    }

    std::unique_ptr<TH1> GenPlot { GenPlotTheory->Rebin(bins,"genplot",Xbins) };

    THStack * stack=  new THStack("def", "def");
    TLegend *leg = new TLegend();
    int legchange = 0;
    std::vector<TH1 *> varhistsPlotting;
    varhistsPlotting.resize(hists.size());

    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        setStyle(varhists[i], i);
        varhistsPlotting[i]=(TH1*)varhists[i]->Clone();
        if(legends.at(i) != "Data"){
            if((legends.at(i) == DYEntry)){
                varhists[i]->Scale(DYScale.at(channelType));
                varhistsPlotting[i]->Scale(DYScale.at(channelType));
            }

            if(i!=(hists.size()-1)){
                if(legends.at(i)!=legends.at(i+1)){
                //std::cout<<legends.at(i)<<std::endl;
                varhistsPlotting[i]->SetLineColor(1);
                }
            }else{
                varhistsPlotting[i]->SetLineColor(1);
            }

            if(legends.at(i) != legends.at(i-1)){
                varhistsPlotting[i]->SetLineColor(1);
                stack->Add(varhistsPlotting[i]);
            }
            if(i > 1){
                if(legends.at(i) != legends.at(i-1)){
                    legchange = i;
                    if( (legends.at(i) == DYEntry) && DYScale.at(channelType)!= 1){
                        leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                    }
                    else leg->AddEntry(varhistsPlotting[i], legends.at(i) ,"f");
                }
                else{varhistsPlotting[legchange]->Add(varhistsPlotting[i]);}
            }
        }
        else{ if(i==0) leg->AddEntry(varhistsPlotting[i], legends.at(i) ,"pe");}
    }

    ///////////////////////////////////
    //purity and stability plots as taken from CombinedCrossSection...

    TH1* genHist = (TH1*)GenPlot->Clone();
    TH1* genRecHist = new TH1D("","",bins,Xbins);
    int intbinsX[XAxisbins.size()];
    int intbinsY[XAxisbins.size()];

    // fill the elements of the main diagonal of the 2d hist into binned 1D histogram
    for (unsigned int i=0; i<XAxisbins.size(); ++i) {
        intbinsX[i] = genReco2d->GetXaxis()->FindBin(Xbins[i]+0.001);
        intbinsY[i] = genReco2d->GetYaxis()->FindBin(Xbins[i]+0.001);

        if (i>0) {
            genRecHist->SetBinContent(i,((TH2D*)genReco2d)->Integral( intbinsX[i-1],intbinsX[i]-1,intbinsY[i-1],intbinsY[i]-1));
        }
    }

    TH1* genPseHist = ((TH2D*)genReco2d)->ProjectionY();
    TH1* recPseHist = ((TH2D*)genReco2d)->ProjectionX();
    
    TH1* genBinHist    = genPseHist->Rebin(bins,"genBinHist", Xbins);
    TH1* recBinHist    = recPseHist->Rebin(bins,"recBinHist", Xbins);

    genRecHist->SetBinContent(0,      0);
    genRecHist->SetBinContent(bins+1,0);
    genBinHist->SetBinContent(0,      0);
    genBinHist->SetBinContent(bins+1,0);
    recBinHist->SetBinContent(0,      0);
    recBinHist->SetBinContent(bins+1,0);
    genHist   ->SetBinContent(0,      0);
    genHist   ->SetBinContent(bins+1,0);

    // this is realy ugly but necessary:
    // As it seems, somewhere a double is tranformed into a float so that
    // efficiencies can be larger than 1.
    for(Int_t i=1; i<=genRecHist->GetNbinsX(); ++i){
      if(genRecHist->GetBinContent(i) > recBinHist->GetBinContent(i)){
        genRecHist->SetBinContent(i,recBinHist->GetBinContent(i));
        std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin" << i
             << " = " << genRecHist->GetBinContent(i) << " is larger than number of reconstructed events in that bin"
             << " = " << recBinHist->GetBinContent(i) << std::endl;
      }
      if(genRecHist->GetBinContent(i) > genBinHist->GetBinContent(i)){
        genRecHist->SetBinContent(i,genBinHist->GetBinContent(i));
            std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin " << i
                 << " is larger than number of generated events in that bin" << std::endl;
      }
    }
    // efficiency, purity, stability
    TGraphAsymmErrors* grE; // for efficiency
    TGraphAsymmErrors* grP; // for purity
    TGraphAsymmErrors* grS; // for stability

    // efficiency
    grE = new TGraphAsymmErrors(recBinHist, genHist);
    grE->SetMinimum(0);
    grE->SetMaximum(1);
    grE->SetLineColor(8);
    grE->SetLineWidth(2);
    grE->SetMarkerSize(2);
    grE->SetMarkerStyle(21);
    grE->SetMarkerColor(8);

    // purity
    grP = new TGraphAsymmErrors(genRecHist, recBinHist);
    grP->SetLineColor(4);
    grP->SetLineWidth(2);
    grP->SetMarkerSize(2);
    grP->SetMarkerStyle(23);
    grP->SetMarkerColor(4);

    // stability
    grS = new TGraphAsymmErrors(genRecHist, genBinHist);
    grS->SetLineColor(2);
    grS->SetLineWidth(2);
    grS->SetMarkerSize(2);
    grS->SetMarkerStyle(22);
    grS->SetMarkerColor(2);


    grE->GetXaxis()->SetTitle(XAxis);
    TCanvas * cESP = new TCanvas("ESP","ESP");

    // this is a dummy to get the x axis range corrct

    recBinHist->Reset();
    recBinHist->Draw();
    recBinHist->SetMaximum(1.);
    recBinHist->GetXaxis()->SetTitle(TString("Reconstructed ").Copy().Append(XAxis));
    recBinHist->GetXaxis()->SetNoExponent(kTRUE);
    grE->GetXaxis()->SetNoExponent(kTRUE);
    grE->Draw("P,SAME");
    grP->Draw("P,SAME");
    grS->Draw("P,SAME");

    TLegend* leg3 = new TLegend();
    leg3->AddEntry(grE, "Efficiency", "p" );
    leg3->AddEntry(grP, "Purity",    "p" );
    leg3->AddEntry(grS, "Stability", "p" );
    setResultLegendStyle(leg3, 0);
    leg3->Draw("SAME");

    TString outdir = common::assignFolder(outpathPlots, Channel, TString("FinalResults"));
    cESP->Print(outdir.Copy()+"ESP_"+name+".eps");
    cESP->Clear();
    delete cESP;

    init = false;
    for (unsigned int hist =0; hist<hists.size(); hist++){
        if(legends.at(hist) == "Data"){
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                DataSum[bin]+=varhists[hist]->GetBinContent(bin+1);
            }
        }
        else if((legends.at(hist) == "t#bar{t} Signal")&&init==false){
            signalHist=hist;
            init=true;
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                GenSignalSum[bin] += GenPlot->GetBinContent(bin+1);
            }
        }
        else{
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                BGSum[bin]+=varhists[hist]->GetBinContent(bin+1);
            }
        }
    }
    double totalDataSum = 0;
    double GenDiffXSecPlot[XAxisbinCenters.size()];
    for (Int_t bin=0; bin<bins; ++bin) {
        totalDataSum+=DataSum[bin];
    }
    TH1 *h_DiffXSec    = (TH1D*)varhists[0]->Clone(); h_DiffXSec->Reset();
    TH1 *h_GenDiffXSec = (TH1D*)varhists[0]->Clone(); h_GenDiffXSec->Reset();

    //The systematic array is filled in the order in which the Stack is filled
    double DiffXSecPlot[XAxisbinCenters.size()];
    double DiffXSecStatErrorPlot[XAxisbinCenters.size()];
    double DiffXSecSysErrorbySysPlot[XAxisbinCenters.size()][(int)vec_systematic.size()];
    double DiffXSecSysErrorPlot[XAxisbinCenters.size()];
    double DiffXSecTotalErrorPlot[XAxisbinCenters.size()];

    double ModelSysPlot[XAxisbinCenters.size()];
    double ExpSysPlot[XAxisbinCenters.size()];

    //Read absolute systematic uncertainty of each bin from file
    for(int Syst=0; Syst<(int)vec_systematic.size(); Syst++){
        if(vec_systematic.at(Syst).Contains("HERWIG") || vec_systematic.at(Syst).Contains("MCATNLO")) continue;

        TString sysup = vec_systematic.at(Syst)+"UP";
        TString sysdown = vec_systematic.at(Syst)+"DOWN";
        if(vec_systematic.at(Syst) == "POWHEG")
        {
            sysup = "POWHEG";
            sysdown = "MCATNLO";
            vec_systematic.at(Syst) = "HAD_";
        };
        ifstream SysResultsList("UnfoldingResults/"+vec_systematic.at(Syst)+"/"+Channel+"/"+name+"Results.txt");
        if(!SysResultsList.is_open()){
            /// Apply symmetrization of errors only if the file doesn't exist!!
            /// If file exists means that the error is set by hand!!
            Plotter::CalcUpDownDifference(Channel, sysup, sysdown, name);
        }
        SysResultsList.close();
//         if(vec_systematic.at(Syst) == "POWHEG") {vec_systematic.at(Syst) = "HAD_";}
        SysResultsList.open("UnfoldingResults/"+vec_systematic.at(Syst)+"/"+Channel+"/"+name+"Results.txt");
        for (Int_t bin=0; bin<bins; ++bin) {
            TString DUMMY;
            SysResultsList>>DUMMY>>XAxisbinCenters[bin]>>DUMMY>>Xbins[bin]>>DUMMY>>Xbins[bin+1]>>DUMMY>>DiffXSecSysErrorbySysPlot[bin][Syst];
        }
        SysResultsList.close();
    }
    TString Dummy;
    //Read central and abolute statistical uncertainty values from Nominal
    ifstream ResultsList("UnfoldingResults/Nominal/"+Channel+"/"+name+"Results.txt");
    for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
        ResultsList>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>DiffXSecPlot[bin]>>Dummy>>DiffXSecStatErrorPlot[bin]>>Dummy>>GenDiffXSecPlot[bin];
        h_DiffXSec->SetBinContent(bin+1,DiffXSecPlot[bin]);
        h_DiffXSec->SetBinError(bin+1,DiffXSecStatErrorPlot[bin]);
        h_GenDiffXSec->SetBinContent(bin+1,GenDiffXSecPlot[bin]);
    }

    double TotalVisXSection = 1.; //this can currently be set to 1. because David's code takes care of the normalization, but just incase we need it
    //double TotalVisXSection = h_DiffXSec->Integral("width");

    h_DiffXSec->Scale(1/TotalVisXSection);

    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
        double syst_square = 0;
        ExpSysPlot[bin]=0.;
        ModelSysPlot[bin]=0.;
        for(int syst =0; syst<(int)vec_systematic.size(); syst++){
            syst_square += DiffXSecSysErrorbySysPlot[bin][syst]*DiffXSecSysErrorbySysPlot[bin][syst];
            if(vec_systematic.at(syst).Contains("JER_") || vec_systematic.at(syst).Contains("JES_") || vec_systematic.at(syst).Contains("PU_") ||
                vec_systematic.at(syst).Contains("DY_") || vec_systematic.at(syst).Contains("BG_") || vec_systematic.at(syst).Contains("TRIG_") ||
                vec_systematic.at(syst).Contains("LEPT_") || vec_systematic.at(syst).Contains("BTAG_") || vec_systematic.at(syst).Contains("KIN_")){
                ExpSysPlot[bin]+=DiffXSecSysErrorbySysPlot[bin][syst]*DiffXSecSysErrorbySysPlot[bin][syst];
            }
            else{
                ModelSysPlot[bin]+=DiffXSecSysErrorbySysPlot[bin][syst]*DiffXSecSysErrorbySysPlot[bin][syst];
            }
        }
        DiffXSecSysErrorPlot[bin]+=syst_square;
        ExpSysPlot[bin]=sqrt(ExpSysPlot[bin]);
        ModelSysPlot[bin]=sqrt(ModelSysPlot[bin]);
        DiffXSecStatErrorPlot[bin] = DiffXSecStatErrorPlot[bin]/TotalVisXSection;
        DiffXSecPlot[bin]=DiffXSecPlot[bin]/TotalVisXSection;
        DiffXSecSysErrorPlot[bin]=sqrt(DiffXSecSysErrorPlot[bin])*DiffXSecPlot[bin]; //absolute systematic error in bin 'bin'
        DiffXSecTotalErrorPlot[bin]=sqrt(DiffXSecSysErrorPlot[bin]*DiffXSecSysErrorPlot[bin] + DiffXSecStatErrorPlot[bin]*DiffXSecStatErrorPlot[bin]);//absolute total error
    }

    //The Markus plots
    TCanvas * c10 = new TCanvas("Markus","Markus");
    THStack* SystHists = new THStack("MSTACK","MSTACK");
    TLegend * leg10 =  new TLegend(0.20,0.65,0.45,0.90);

    FILE *systfile;
    systfile = fopen(outdir.Copy()+newname+"_SystematicsLaTeX.txt", "w");
    for(int systs =0; systs<(int)vec_systematic.size(); systs++){
        if (vec_systematic.at(systs) == "BTAG_ETA_" || vec_systematic.at(systs) == "BTAG_LJET_ETA_") {continue;}//Skip the BTAG_ETA systematic because it's added in quadrature to BTAG_PT
        TH1D* systtemp = (TH1D*)varhists[0]->Clone();
        systtemp->Reset();
        double TotalSyst=0.0, TotalSqSyst=0.0;
        double AvgSyst= 0.0, SqAvgSys=0.0;

        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            if(vec_systematic.at(systs) == "BTAG_PT_" || vec_systematic.at(systs) == "BTAG_LJET_PT_"){
                DiffXSecSysErrorbySysPlot[bin][systs]= TMath::Sqrt(
                    (DiffXSecSysErrorbySysPlot[bin][systs]*DiffXSecSysErrorbySysPlot[bin][systs])
                    +(DiffXSecSysErrorbySysPlot[bin][systs+1]*DiffXSecSysErrorbySysPlot[bin][systs+1])
                );
            }

            systtemp->SetBinContent(bin+1,(DiffXSecSysErrorbySysPlot[bin][systs]*DiffXSecSysErrorbySysPlot[bin][systs]));
            if(bin==0){
                fprintf(systfile, "%s", (vec_systematic.at(systs)+" ").Data());
            }
            fprintf(systfile, "%2.5f ", TMath::Sqrt(systtemp->GetBinContent(bin+1))*100);
            if(bin>0 && bin<bins-1){//Exclude the 2 side bins
                TotalSyst=TotalSyst+TMath::Sqrt(systtemp->GetBinContent(bin+1));
                TotalSqSyst=TotalSqSyst+systtemp->GetBinContent(bin+1);
            }
        }
        AvgSyst=TotalSyst/(bins-2);
        SqAvgSys=TMath::Sqrt(TotalSqSyst/(bins-2));
        fprintf(systfile, "Lin.Avg.(%%)= %.5f  Quad.Avg.(%%)=%.5f\n", 100*AvgSyst, 100*SqAvgSys);
        systtemp->SetFillColor((int)vec_systematic.size()-systs);
        SystHists->Add(systtemp);
        leg10->AddEntry(systtemp, vec_systematic.at(systs), "f");
    }
    SystHists->Draw();
    fclose(systfile);

    if(vec_systematic.size()>0){
        if(name.Contains("pT") ||name.Contains("Mass") ){
            SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis.Copy().Append(" #left[GeV#right]"));
            if(name.Contains("Rapidity")) SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis);
        }
        else  SystHists->GetHistogram()->GetXaxis()->SetTitle(XAxis);
        SystHists->GetHistogram()->GetYaxis()->SetTitle("#sum #left( #frac{#Delta #sigma}{#sigma} #right)^{2}");
        SystHists->GetXaxis()->SetNoExponent(kTRUE);

        leg10->SetFillColor(0);
        leg10->Draw("SAME");
        c10->Print(outdir.Copy()+"MSP_"+name+".eps");
        delete leg10;
    }
    delete c10;

    //The Experimental/Model/Statistical plot
    TCanvas * c11 = new TCanvas("EMS","EMS");
    TH1D* ExpHist = (TH1D*)varhists[0]->Clone();    ExpHist->Reset();
    TH1D* ModelHist = (TH1D*)varhists[0]->Clone();  ModelHist->Reset();
    TH1D* StatHist = (TH1D*)varhists[0]->Clone();   StatHist->Reset();
    TH1D* TotalHist = (TH1D*)varhists[0]->Clone();  TotalHist->Reset();
    for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
        ExpHist->SetBinContent(bin+1,100*ExpSysPlot[bin]);
        ModelHist->SetBinContent(bin+1,100*ModelSysPlot[bin]);
        StatHist->SetBinContent(bin+1,100*DiffXSecStatErrorPlot[bin]/DiffXSecPlot[bin]);
        TotalHist->SetBinContent(bin+1,100*DiffXSecTotalErrorPlot[bin]/DiffXSecPlot[bin]);
    }
    TotalHist->SetMinimum(0.);
    TotalHist->GetYaxis()->SetTitle("#frac{#Delta#sigma}{#sigma} [%]");
    TotalHist->SetLineColor(1);
    ExpHist->SetLineColor(kRed);
    StatHist->SetLineColor(kGreen);
    ModelHist->SetLineColor(kBlue);
    TLegend * leg11 =  new TLegend(0.65,0.60,0.90,0.85);
    leg11->SetFillColor(0);
    leg11->AddEntry(ExpHist->Clone(), "Experimental Uncertainty", "l");
    leg11->AddEntry(StatHist->Clone(), "Statistical Uncertainty", "l");
    leg11->AddEntry(ModelHist->Clone(), "Model Uncertainty", "l");
    leg11->AddEntry(TotalHist->Clone(), "Total Uncertainty", "l");
    TotalHist->Draw();ModelHist->Draw("SAME");ExpHist->Draw("SAME");StatHist->Draw("SAME");
    leg11->Draw("SAME");
    TotalHist->GetXaxis()->SetNoExponent(kTRUE);
    c11->Print(outdir.Copy()+"SEM_"+name+".eps");
    c11->Clear();
    delete c11;

    Double_t mexl[XAxisbinCenters.size()];
    Double_t mexh[XAxisbinCenters.size()];
    for (unsigned int j=0; j<XAxisbinCenters.size();j++){mexl[j]=0;mexh[j]=0;}
    TGraphAsymmErrors *tga_DiffXSecPlot = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecStatErrorPlot, DiffXSecStatErrorPlot);
    tga_DiffXSecPlot->SetMarkerStyle(1);
    tga_DiffXSecPlot->SetMarkerColor(kBlack);
    tga_DiffXSecPlot->SetMarkerSize(0.8);
    tga_DiffXSecPlot->SetLineWidth(2);
    tga_DiffXSecPlot->SetLineColor(kBlack);
   
    TGraphAsymmErrors *tga_DiffXSecPlotwithSys = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecTotalErrorPlot, DiffXSecTotalErrorPlot);
    tga_DiffXSecPlotwithSys->SetMarkerStyle(20);
    tga_DiffXSecPlotwithSys->SetMarkerColor(kBlack);
    tga_DiffXSecPlotwithSys->SetMarkerSize(0.8);
    tga_DiffXSecPlotwithSys->SetLineWidth(2);
    tga_DiffXSecPlotwithSys->SetLineColor(kBlack);

    //GenPlotTheory->Scale(topxsec/(SignalEventswithWeight*GenPlotTheory->GetBinWidth(1)));
    GenPlot->Scale(topxsec/(SignalEventswithWeight*GenPlot->GetBinWidth(1)));
    if( (name.Contains("Lepton")||name.Contains("Top")||name.Contains("BJet")) && !name.Contains("Lead")){
      GenPlotTheory->Scale(1./2.);
    }
//     GenPlotTheory->Rebin(2);
    GenPlotTheory->Scale(GenPlotTheory->Integral("width"));
    h_GenDiffXSec->Scale(h_GenDiffXSec->Integral("width"));

    bool binned_theory=true; //############

    TH1* madgraphhist = 0, *mcnlohist=0, *mcnlohistup=0, *mcnlohistdown=0, *powheghist=0, *powhegHerwighist=0, *perugia11hist = 0;
    TH1* mcnlohistnorm=0;
    TGraph *mcatnloBand=0;

    TH1* mcnlohistnormBinned = 0, *mcnlohistupBinned = 0;
    TH1* mcnlohistdownBinned = 0, *mcnlohistBinned = 0;
    TH1* madgraphhistBinned = 0, *powheghistBinned = 0, *powhegHerwighistBinned = 0;
    TH1* perugia11histBinned = 0;

    TH1 *Kidoth1_Binned = 0, *Ahrensth1_Binned = 0;

    TH1* madup=0, *maddown=0, *matchup=0, *matchdown=0, *match2up = 0, *match2down = 0;
    TH1* madupBinned = 0, *maddownBinned = 0, *matchupBinned = 0, *matchdownBinned = 0, *match2upBinned = 0, *match2downBinned = 0;

    madgraphhist = GetNloCurve(newname, "Nominal");
    madgraphhist->Scale(1./madgraphhist->Integral("width"));
    madgraphhistBinned = madgraphhist->Rebin(bins,"madgraphplot",Xbins);
    for (Int_t bin=0; bin<bins; bin++){
        madgraphhistBinned->SetBinContent(bin+1,madgraphhistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/madgraphhist->GetBinWidth(1)));
    }
    madgraphhistBinned->Scale(1./madgraphhistBinned->Integral("width"));

    bool canDrawMCATNLO = true;
    if (drawNLOCurves && drawMCATNLO) {
        mcnlohist = GetNloCurve(newname,"MCATNLO");
        double mcnloscale = 1./mcnlohist->Integral("width");
        if (binned_theory==false) mcnlohist->Rebin(2);mcnlohist->Scale(0.5); //#####
        mcnlohist->Scale(mcnloscale);

        mcnlohistBinned = mcnlohist->Rebin(bins,"mcnloplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            mcnlohistBinned->SetBinContent(bin+1,mcnlohistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohist->GetBinWidth(1)));
            mcnlohistBinned->SetBinError(bin+1, 1e-10);
        }
        mcnlohistBinned->Scale(1./mcnlohistBinned->Integral("width"));

        if(name.Contains("LeptonpT")){mcnlohistnorm = GetNloCurve("Leptons","Pt","MCatNLO");}//temprorary until I change the naming convention in the root file
        else if(name.Contains("LeptonEta")){mcnlohistnorm = GetNloCurve("Leptons","Eta","MCatNLO");}
        else if(name.Contains("LLBarpT")){mcnlohistnorm = GetNloCurve("LepPair","Pt","MCatNLO");}
        else if(name.Contains("LLBarMass")){mcnlohistnorm = GetNloCurve("LepPair","Mass","MCatNLO");}
        else if(name.Contains("ToppT")){mcnlohistnorm = GetNloCurve("TopQuarks","Pt","MCatNLO");}
        else if(name.Contains("TopRapidity")){mcnlohistnorm = GetNloCurve("TopQuarks","Rapidity","MCatNLO");}
        else if(name.Contains("TTBarpT")){mcnlohistnorm = GetNloCurve("TtBar","Pt","MCatNLO");}
        else if(name.Contains("TTBarRapidity")){mcnlohistnorm = GetNloCurve("TtBar","Rapidity","MCatNLO");}
        else if(name.Contains("TTBarMass")){mcnlohistnorm = GetNloCurve("TtBar","Mass","MCatNLO");}
        else if(name.Contains("BJetpT")){mcnlohistnorm = GetNloCurve("Jets","Pt","MCatNLO");}
        else if(name.Contains("BJetEta")){mcnlohistnorm = GetNloCurve("Jets","Eta","MCatNLO");}
        //else if(name.Contains("LeptonBJetMass")){mcnlohistnorm = GetNloCurve("Jets","Eta","MCatNLO");}

        //    if (binned_theory==false) mcnlohistnorm->Rebin(5);mcnlohistnorm->Scale(0.2);
        if (mcnlohistnorm) {
            mcnlohistnormBinned = mcnlohistnorm->Rebin(bins,"genBinHistNorm", Xbins);

            if(name.Contains("LeptonpT")){mcnlohistup = GetNloCurve("Leptons","Pt","MCNLOup");}//temprorary until I change the naming convention in the root file
            else if(name.Contains("LeptonEta")){mcnlohistup = GetNloCurve("Leptons","Eta","MCNLOup");}
            else if(name.Contains("LLBarpT")){mcnlohistup = GetNloCurve("LepPair","Pt","MCNLOup");}
            else if(name.Contains("LLBarMass")){mcnlohistup = GetNloCurve("LepPair","Mass","MCNLOup");}
            else if(name.Contains("ToppT")){mcnlohistup = GetNloCurve("TopQuarks","Pt","MCNLOup");}
            else if(name.Contains("TopRapidity")){mcnlohistup = GetNloCurve("TopQuarks","Rapidity","MCNLOup");}
            else if(name.Contains("TTBarpT")){mcnlohistup = GetNloCurve("TtBar","Pt","MCNLOup");}
            else if(name.Contains("TTBarRapidity")){mcnlohistup = GetNloCurve("TtBar","Rapidity","MCNLOup");}
            else if(name.Contains("TTBarMass")){mcnlohistup = GetNloCurve("TtBar","Mass","MCNLOup");}
            else if(name.Contains("BJetpT")){mcnlohistup = GetNloCurve("Jets","Pt","MCNLOup");}
            else if(name.Contains("BJetEta")){mcnlohistup = GetNloCurve("Jets","Eta","MCNLOup");}
            //    if (binned_theory==false) mcnlohistup->Rebin(5);mcnlohistup->Scale(0.2);
            mcnlohistupBinned    = mcnlohistup->Rebin(bins,"genBinHist", Xbins);


            if(name.Contains("LeptonpT")){mcnlohistdown = GetNloCurve("Leptons","Pt","MCNLOdown");}//temprorary until I change the naming convention in the root file
            else if(name.Contains("LeptonEta")){mcnlohistdown = GetNloCurve("Leptons","Eta","MCNLOdown");}
            else if(name.Contains("LLBarpT")){mcnlohistdown = GetNloCurve("LepPair","Pt","MCNLOdown");}
            else if(name.Contains("LLBarMass")){mcnlohistdown = GetNloCurve("LepPair","Mass","MCNLOdown");}
            else if(name.Contains("ToppT")){mcnlohistdown = GetNloCurve("TopQuarks","Pt","MCNLOdown");}
            else if(name.Contains("TopRapidity")){mcnlohistdown = GetNloCurve("TopQuarks","Rapidity","MCNLOdown");}
            else if(name.Contains("TTBarpT")){mcnlohistdown = GetNloCurve("TtBar","Pt","MCNLOdown");}
            else if(name.Contains("TTBarRapidity")){mcnlohistdown = GetNloCurve("TtBar","Rapidity","MCNLOdown");}
            else if(name.Contains("TTBarMass")){mcnlohistdown = GetNloCurve("TtBar","Mass","MCNLOdown");}
            else if(name.Contains("BJetpT")){mcnlohistdown = GetNloCurve("Jets","Pt","MCNLOdown");}
            else if(name.Contains("BJetEta")){mcnlohistdown = GetNloCurve("Jets","Eta","MCNLOdown");}
            //    if (binned_theory==false) mcnlohistdown->Rebin(5);mcnlohistdown->Scale(0.2);
            mcnlohistdownBinned    = mcnlohistdown->Rebin(bins,"genBinHist", Xbins);

            for (Int_t bin=0; bin<bins; bin++){
                mcnlohistupBinned->SetBinContent(bin+1,mcnlohistupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistup->GetBinWidth(1)));
                mcnlohistdownBinned->SetBinContent(bin+1,mcnlohistdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistdown->GetBinWidth(1)));
                mcnlohistnormBinned->SetBinContent(bin+1,mcnlohistnormBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistnorm->GetBinWidth(1)));
            }
            mcnlohistupBinned->Scale(1./mcnlohistnormBinned->Integral("width"));
            mcnlohistdownBinned->Scale(1./mcnlohistnormBinned->Integral("width"));
            mcnlohistnormBinned->Scale(1./mcnlohistnormBinned->Integral("width"));

            for (Int_t bin=0; bin<bins; bin++){
                mcnlohistupBinned->SetBinContent(bin+1,(mcnlohistupBinned->GetBinContent(bin+1)/mcnlohistnormBinned->GetBinContent(bin+1))*mcnlohistBinned->GetBinContent(bin+1));
                mcnlohistdownBinned->SetBinContent(bin+1,(mcnlohistdownBinned->GetBinContent(bin+1)/mcnlohistnormBinned->GetBinContent(bin+1))*mcnlohistBinned->GetBinContent(bin+1));
            }

            //Uncertainty band for MC@NLO
            Double_t x[bins];
            Double_t xband[2*bins];
            Double_t errup[bins];
            Double_t errdn[bins];
            Double_t errorband[2*bins];

            for( Int_t j = 0; j< bins; j++ ){
                x[j]=mcnlohistBinned->GetBinCenter(j+1);
                errup[j]=(mcnlohistupBinned->GetBinContent(j+1)/mcnlohistnormBinned->GetBinContent(j+1))*mcnlohistBinned->GetBinContent(j+1);
                errdn[j]=(mcnlohistdownBinned->GetBinContent(j+1)/mcnlohistnormBinned->GetBinContent(j+1))*mcnlohistBinned->GetBinContent(j+1);

                xband[j] = x[j];
                errorband[j] = errdn[j]; //lower band
                xband[2*bins-j-1] = x[j];
                errorband[2*bins-j-1] = errup[j]; //upper band
            }

            mcatnloBand = new TGraph(2*bins, xband, errorband);
            mcatnloBand->SetFillColor(kGray);
            mcatnloBand->SetFillStyle(1001);
            mcatnloBand->SetLineColor(kBlue);
            mcatnloBand->SetLineWidth(2);
            mcatnloBand->SetLineStyle(5);
            canDrawMCATNLO = false;
        } else {
            std::cout << "\n*************************\nMC@NLO Curve not available!\n**********************\n";
            canDrawMCATNLO = false;
        }
    }
    if(drawNLOCurves && drawPOWHEG){
        powheghist = GetNloCurve(newname, "POWHEG");
        double powhegscale = 1./powheghist->Integral("width");
        if (binned_theory==false) powheghist->Rebin(2);powheghist->Scale(0.5);
        powheghist->Scale(powhegscale);
        powheghistBinned = powheghist->Rebin(bins,"powhegplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            powheghistBinned->SetBinContent(bin+1,powheghistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/powheghist->GetBinWidth(1)));
            powheghistBinned->SetBinError(bin+1, 1e-10);
        }
        powheghistBinned->Scale(1./powheghistBinned->Integral("width"));
    }
    if(drawNLOCurves && drawPOWHEGHERWIG){
        powhegHerwighist = GetNloCurve(newname, "POWHEGHERWIG");
        double powhegHerwigscale = 1./powhegHerwighist->Integral("width");
        if (binned_theory==false) powhegHerwighist->Rebin(2);powhegHerwighist->Scale(0.5);
        powhegHerwighist->Scale(powhegHerwigscale);
        powhegHerwighistBinned = powhegHerwighist->Rebin(bins,"powhegHerwigplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            powhegHerwighistBinned->SetBinContent(bin+1,powhegHerwighistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/powhegHerwighist->GetBinWidth(1)));
            powhegHerwighistBinned->SetBinError(bin+1, 1e-10);
        }
        powhegHerwighistBinned->Scale(1./powhegHerwighistBinned->Integral("width"));
    }

    if(drawNLOCurves && drawPERUGIA11){
        perugia11hist = GetNloCurve(newname, "PERUGIA11");
        double perugia11histscale = 1./perugia11hist->Integral("width");
        perugia11hist->Scale(perugia11histscale);

        perugia11histBinned = perugia11hist->Rebin(bins,"perugia11plot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            perugia11histBinned->SetBinContent(bin+1,perugia11histBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/perugia11hist->GetBinWidth(1)));
            perugia11histBinned->SetBinError(bin+1, 1e-10);
        }
        perugia11histBinned->Scale(1./perugia11histBinned->Integral("width"));
    }
    if(drawNLOCurves && drawKidonakis &&
        (name== "HypToppT" || name == "HypTopRapidity") && !name.Contains("Lead") && !name.Contains("RestFrame")){
        TString kidoFile = ttbar::DATA_PATH_DILEPTONIC() + "/kidonakisNNLO_8TeV.root";
        if(name.Contains("ToppT"))             Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topPt");
        else if(name.Contains("TopRapidity"))  Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topY");
        for (int iter = 0; iter<=Kidoth1_Binned->GetNbinsX(); iter++) Kidoth1_Binned->SetBinError(iter, 1e-10);
        Kidoth1_Binned->Scale(1./Kidoth1_Binned->Integral("width"));
    }

    if(drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")){
        TString ahrensFile = ttbar::DATA_PATH_DILEPTONIC() + "/ahrensNNLL_8TeV.root";
        if(name == "HypTTBarMass")    Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarM");
        else if(name == "HypTTBarpT") Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarPt");
        for (int iter = 0; iter<=Ahrensth1_Binned->GetNbinsX(); iter++) Ahrensth1_Binned->SetBinError(iter, 1e-10);
        Ahrensth1_Binned->Scale(1./Ahrensth1_Binned->Integral("width"));
    } else {drawAhrens = 0;}

    if(drawMadScaleMatching){
//    if(drawNLOCurves && drawMadScaleMatching){
        madup = GetNloCurve(newname,"SCALE_UP");
        double madscale = 1./madup->Integral("width");
        if (binned_theory==false) madup->Rebin(2);madup->Scale(0.5);
        madup->Scale(madscale);
        madupBinned = madup->Rebin(bins,"madupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
          madupBinned->SetBinContent(bin+1,madupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/madup->GetBinWidth(1)));
        }
        madupBinned->Scale(1./madupBinned->Integral("width"));

        maddown = GetNloCurve(newname,"SCALE_DOWN");
        double maddownscale = 1./maddown->Integral("width");
        if (binned_theory==false) maddown->Rebin(2);maddown->Scale(0.5);
        maddown->Scale(maddownscale);
        maddownBinned = maddown->Rebin(bins,"maddownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            maddownBinned->SetBinContent(bin+1,maddownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/maddown->GetBinWidth(1)));
        }
        maddownBinned->Scale(1./maddownBinned->Integral("width"));

        matchup = GetNloCurve(newname,"MATCH_UP");
        double matchscale = 1./matchup->Integral("width");
        if (binned_theory==false) matchup->Rebin(2);matchup->Scale(0.5);
        matchup->Scale(matchscale);
        matchupBinned = matchup->Rebin(bins,"matchupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            matchupBinned->SetBinContent(bin+1,matchupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchup->GetBinWidth(1)));
        }
        matchupBinned->Scale(1./matchupBinned->Integral("width"));

        matchdown = GetNloCurve(newname,"MATCH_DOWN");
        double matchdownscale = 1./matchdown->Integral("width");
        if (binned_theory==false) matchdown->Rebin(2);matchdown->Scale(0.5);
        matchdown->Scale(matchdownscale);
        matchdownBinned = matchdown->Rebin(bins,"matchdownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
           matchdownBinned->SetBinContent(bin+1,matchdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchdown->GetBinWidth(1)));
        }
        matchdownBinned->Scale(1./matchdownBinned->Integral("width"));

    }

    if(drawNLOCurves && drawMadMass){
        madup = GetNloCurveMass(newname,"MASS_UP","173.5");
        madup->Scale(1./madup->Integral("width"));
        madupBinned = madup->Rebin(bins,"madupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
          madupBinned->SetBinContent(bin+1,madupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/madup->GetBinWidth(1)));
        }
        madupBinned->Scale(1./madupBinned->Integral("width"));

        maddown = GetNloCurveMass(newname,"MASS_DOWN","171.5");
        maddown->Scale(1./maddown->Integral("width"));
        maddownBinned = maddown->Rebin(bins,"maddownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            maddownBinned->SetBinContent(bin+1,maddownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/maddown->GetBinWidth(1)));
        }
        maddownBinned->Scale(1./maddownBinned->Integral("width"));

        matchup = GetNloCurveMass(newname,"MASS_UP","175.5");
        matchup->Scale(1./matchup->Integral("width"));
        matchupBinned = matchup->Rebin(bins,"matchupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            matchupBinned->SetBinContent(bin+1,matchupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchup->GetBinWidth(1)));
        }
        matchupBinned->Scale(1./matchupBinned->Integral("width"));

        matchdown = GetNloCurveMass(newname,"MASS_DOWN","169.5");
        matchdown->Scale(1./matchdown->Integral("width"));
        matchdownBinned = matchdown->Rebin(bins,"matchdownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
           matchdownBinned->SetBinContent(bin+1,matchdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchdown->GetBinWidth(1)));
        }
        matchdownBinned->Scale(1./matchdownBinned->Integral("width"));

        match2up = GetNloCurveMass(newname,"MASS_UP","178.5");
        double match2scale = 1./match2up->Integral("width");
        match2up->Scale(match2scale);
        match2upBinned = match2up->Rebin(bins,"match2upplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            match2upBinned->SetBinContent(bin+1,match2upBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/match2up->GetBinWidth(1)));
        }
        match2upBinned->Scale(1./match2upBinned->Integral("width"));
        
        match2down = GetNloCurveMass(newname,"MASS_DOWN","166.5");
        match2down->Scale(1./match2down->Integral("width"));
        match2downBinned = match2down->Rebin(bins,"match2downplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
           match2downBinned->SetBinContent(bin+1,match2downBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/match2down->GetBinWidth(1)));
        }
        match2downBinned->Scale(1./match2downBinned->Integral("width"));

    }

    TCanvas * c = new TCanvas("DiffXS","DiffXS");
    if(logY){
      c->SetLogy();
    }
    h_DiffXSec->SetMarkerStyle(tga_DiffXSecPlotwithSys->GetMarkerStyle());
    h_DiffXSec->SetMarkerSize(tga_DiffXSecPlotwithSys->GetMarkerSize());
    h_DiffXSec->SetMarkerColor(tga_DiffXSecPlotwithSys->GetMarkerColor());
    //MCFMHist->SetMarkerStyle(2);
    if(ymax!=0){
        if(logY){
            madgraphhistBinned->SetMaximum(18*madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));
        }
        else{ madgraphhistBinned->SetMaximum(1.5*madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));}
    }
    madgraphhistBinned->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("Rapidity") || name.Contains("Eta")){madgraphhistBinned->GetYaxis()->SetNoExponent(kTRUE);}
    if (name.Contains("HypJetMultpt")) {
        //madgraphhistBinned->GetXaxis()->SetNdivisions(madgraphhistBinned->GetNbinsX(),0,0, 1);
        TString TitBin = "";
        for(int bin = 1; bin <= madgraphhistBinned->GetNbinsX(); bin++) {
            if( bin == madgraphhistBinned->GetNbinsX()) {TitBin += "#geq"; TitBin += madgraphhistBinned->GetBinCenter(bin); madgraphhistBinned->GetXaxis()->SetBinLabel(bin,TitBin);}
            else{TitBin += madgraphhistBinned->GetBinCenter(bin);
            madgraphhistBinned->GetXaxis()->SetBinLabel(bin,TitBin);
            }
            TitBin  = "";
       }
    }


    if (ymax!=0) madgraphhistBinned->SetMaximum(ymax);
    if (ymin!=0) madgraphhistBinned->SetMinimum(ymin);
    gStyle->SetEndErrorSize(8);
    if (drawNLOCurves && drawMCATNLO && canDrawMCATNLO) {
    //    mcatnloBand->Draw("same, F");
        mcnlohistupBinned->SetFillColor(kGray);
        mcnlohistupBinned->SetLineColor(kGray);
        mcnlohistupBinned->Draw("same");
        mcnlohistdownBinned->SetLineColor(10);
        mcnlohistdownBinned->SetFillColor(10);
        mcnlohistdownBinned->Draw("same");
    }
    GenPlotTheory->SetLineColor(kRed+1);
    GenPlotTheory->SetLineWidth(2);
    GenPlotTheory->SetLineStyle(1);
    h_GenDiffXSec->SetLineColor(GenPlotTheory->GetLineColor());
    h_GenDiffXSec->SetLineStyle(GenPlotTheory->GetLineStyle());

    // Plot statistical and totl error of full result
    TGraphAsymmErrors *ratio_stat = 0, *ratio_tota = 0;
    if(tga_DiffXSecPlot) ratio_stat = (TGraphAsymmErrors*)tga_DiffXSecPlot->Clone("ratio_stat");
    if(tga_DiffXSecPlotwithSys) ratio_tota = (TGraphAsymmErrors*)tga_DiffXSecPlotwithSys->Clone("ratio_tota");

    if(ratio_stat){
        ratio_stat->SetFillStyle(1001);
        ratio_stat->SetFillColor(kGray+1);
        ratio_stat->SetLineColor(0);
        for (Int_t iter = 0; iter<tga_DiffXSecPlot->GetN(); iter++)
        {
            double binWidth = (XAxisbins[iter+1] - XAxisbins[iter])/2;
            double x = tga_DiffXSecPlot->GetX()[iter];
            double y_ratio_stat = tga_DiffXSecPlot->GetY()[iter] / tga_DiffXSecPlot->GetY()[iter];
            double abserr_ratio_stat = y_ratio_stat - std::abs(tga_DiffXSecPlot->GetErrorY(iter) - tga_DiffXSecPlot->GetY()[iter]) / tga_DiffXSecPlot->GetY()[iter];
            ratio_stat->SetPoint(iter, x, y_ratio_stat);
            ratio_stat->SetPointError(iter, binWidth, binWidth, abserr_ratio_stat, abserr_ratio_stat);
        }
    }
    if(ratio_tota){
        ratio_tota->SetFillStyle(1001);
        ratio_tota->SetFillColor(kOrange-4);
        ratio_tota->SetLineColor(0);
        for (Int_t iter = 0; iter<tga_DiffXSecPlotwithSys->GetN(); iter++)
        {
            double binWidth = (XAxisbins[iter+1] - XAxisbins[iter])/2;
            double x = tga_DiffXSecPlotwithSys->GetX()[iter];
            double y_ratio_tota = tga_DiffXSecPlotwithSys->GetY()[iter] / tga_DiffXSecPlotwithSys->GetY()[iter];
            double abserr_ratio_tota = y_ratio_tota - std::abs(tga_DiffXSecPlotwithSys->GetErrorY(iter) - tga_DiffXSecPlotwithSys->GetY()[iter]) / tga_DiffXSecPlotwithSys->GetY()[iter];
            ratio_tota->SetPoint(iter, x, y_ratio_tota);
            ratio_tota->SetPointError(iter, binWidth, binWidth, abserr_ratio_tota, abserr_ratio_tota);
        }
    }

    TLegend *leg2 = new TLegend();
    if(doClosureTest){
        leg2->AddEntry(h_DiffXSec, "Pseudo-Data", "p");
    } else {
        leg2->AddEntry(h_DiffXSec, "Data", "p");
    }

    setTheoryStyleAndFillLegend(h_DiffXSec, "data");
    setTheoryStyleAndFillLegend(madgraphhist, "madgraph");
    setTheoryStyleAndFillLegend(madgraphhistBinned, "madgraph", leg2);
    madgraphhistBinned->GetXaxis()->SetTitle(varhists[0]->GetXaxis()->GetTitle());
    madgraphhistBinned->GetYaxis()->SetTitle(varhists[0]->GetYaxis()->GetTitle());
    madgraphhistBinned->Draw();


    if (drawNLOCurves && drawMCATNLO) {
        setTheoryStyleAndFillLegend(mcnlohist, "mcatnloherwig");
        setTheoryStyleAndFillLegend(mcnlohistBinned, "mcatnloherwig", leg2);
        mcnlohistBinned->Draw("SAME,HISTO");
    }
    if(drawNLOCurves && drawPOWHEGHERWIG){
        setTheoryStyleAndFillLegend(powhegHerwighist, "powhegherwig");
        setTheoryStyleAndFillLegend(powhegHerwighistBinned, "powhegherwig", leg2);
        powhegHerwighistBinned->Draw("SAME,HISTO");
    }
    if(drawNLOCurves && drawPOWHEG){
        setTheoryStyleAndFillLegend(powheghist, "powhegpythia");
        setTheoryStyleAndFillLegend(powheghistBinned, "powhegpythia", leg2);
        powheghistBinned->Draw("SAME,HISTO");
    }
    if(drawNLOCurves && drawPERUGIA11){
        setTheoryStyleAndFillLegend(perugia11hist, "perugia11");
        setTheoryStyleAndFillLegend(perugia11histBinned, "perugia11", leg2);
        perugia11histBinned->Draw("SAME,HISTO");
    }
    if(drawNLOCurves && drawKidonakis &&
        (name== "HypToppT" || name == "HypTopRapidity") &&
        !name.Contains("Lead") && !name.Contains("RestFrame")){
        setTheoryStyleAndFillLegend(Kidoth1_Binned, "kidonakis", leg2);
        Kidoth1_Binned->Draw("SAME,HISTO");
    }
    if(drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT"))
    {
        setTheoryStyleAndFillLegend(Ahrensth1_Binned, "ahrens", leg2);
        Ahrensth1_Binned->Draw("SAME,HISTO");
    }
    if(drawNLOCurves && (drawMadScaleMatching || drawMadMass)){
        if(drawMadScaleMatching){
            setTheoryStyleAndFillLegend(madupBinned,    "scaleup",   leg2);
            setTheoryStyleAndFillLegend(maddownBinned,  "scaledown", leg2);
            setTheoryStyleAndFillLegend(matchupBinned,  "matchup",   leg2);
            setTheoryStyleAndFillLegend(matchdownBinned,"matchdown", leg2);
        }
        if(drawMadMass){
            setTheoryStyleAndFillLegend(madupBinned,    "mass173.5", leg2);
            setTheoryStyleAndFillLegend(maddownBinned,  "mass171.5", leg2);
            setTheoryStyleAndFillLegend(matchupBinned,  "mass175.5", leg2);
            setTheoryStyleAndFillLegend(matchdownBinned,"mass169.5", leg2);
            setTheoryStyleAndFillLegend(match2upBinned, "mass178.5", leg2);
            setTheoryStyleAndFillLegend(match2downBinned,"mass166.5", leg2);
        }
        madupBinned->Draw("SAME");
        maddownBinned->Draw("SAME");
        matchupBinned->Draw("SAME");
        matchdownBinned->Draw("SAME");
        match2upBinned->Draw("SAME");
        match2downBinned->Draw("SAME");

        ofstream OutputFileXSec(string("Plots/"+Channel+"/"+name+"DiffXsecMass.txt"));
        for(int i = 1; i < madupBinned->GetNbinsX(); i++){
        //OutputFileXSec<<"Nominal "<<"Mass 181 GeV" << " Mass 175 GeV"<< "Mass 169 GeV" << "Mass 163 GeV"<<endl;                             OutputFileXSec<<h_DiffXSec->GetBinContent(i)<< " "<<tga_DiffXSecPlot->GetErrorY(i-1)<<" "<<tga_DiffXSecPlotwithSys->GetErrorY(i-1)<< " "<<h|
        }
        OutputFileXSec.close();
    }

    if (drawNLOCurves) {
        //if (drawPOWHEG && powheghist->GetEntries())                                                     leg2->AddEntry(powheghistBinned, "tt+1jet m=178.5 GeV","l");
        //if (drawPOWHEGHERWIG && powhegHerwighist->GetEntries())                                         leg2->AddEntry(powhegHerwighistBinned, "tt+1jet m=172.5 GeV","l");
    }

    if (drawSmoothMadgraph){
        if( name.Contains("HypTTBarMass")){
            GenPlotTheory->Rebin(15);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(20, "R");
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypBJetpTNLead")){
            GenPlotTheory->Rebin(5);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if(name.Contains("HypBJetEtaLead") || name.Contains("HypBJetEtaNLead") || name.Contains("HypBJetpTLead")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonEtaNLead")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(5);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonEta") || name.Contains("HypTopRapidity") || 
                name.Contains("HypBJetEta") || name.Contains("HypTTBarRapidity")){
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if (name.Contains("HypToppT")){
            GenPlotTheory->Rebin(6);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(10);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonpTLead")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(4);
            SmoothMadgraph->Draw("SAME, L");
        }
        else if(name.Contains("HypLeptonpTNLead")){
            GenPlotTheory->Rebin(3);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if(name.Contains("HypTTBarpT")){
            GenPlotTheory->Rebin(5);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if(name.Contains("HypLLBarpT")){
            GenPlotTheory->Rebin(2);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            GenPlotTheory->Draw("same,c");
        }
        else if(name.Contains("HypLLBarMass")){
            GenPlotTheory->Rebin(4);
            GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
            TH1D *SmoothMadgraph =(TH1D*)GenPlotTheory->Clone("SmoothMadgraph");
            SmoothMadgraph->Smooth(2);
            SmoothMadgraph->Draw("same,l");
        }
        else {GenPlotTheory->Draw("SAME,C");} //### 150512 ### 
    }

//     madgraphhistBinned->Draw("SAME");
    DrawCMSLabels(0, 8);
    DrawDecayChLabel(channelLabel[channelType]);


    if(name.Contains("JetMult")) {
        TString legtit = "";
        if(name.Contains("pt60")) legtit += "p_{T}^{jet}> 60 GeV";
        else if(name.Contains("pt100")) legtit += "p_{T}^{jet}> 100 GeV";
        else legtit += "p_{T}^{jet}> 30 GeV";
        leg2->SetHeader(legtit);
    }
    setResultLegendStyle(leg2);
    leg2->Draw("same");

    if (drawNLOCurves && drawKidonakis &&  (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame")){
        DrawLabel("(arXiv:1210.7813)", leg2->GetX1NDC()+0.06, leg2->GetY1NDC()-0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }
    if (drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")) {
        DrawLabel("(arXiv:1003.5827)", leg2->GetX1NDC()+0.06, leg2->GetY1NDC()-0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }

//     madgraphhistBinned->Draw("SAME");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    tga_DiffXSecPlotwithSys->Draw("p, SAME, Z");
    gPad->RedrawAxis();

    if(drawPlotRatio) {
        TH1D *tmpKido = 0, *tmpAhrens = 0;
        if(Kidoth1_Binned){
            /// Idiot definition of temporary histogram for Kidonakis due to the larger number of bins in raw histogram
            tmpKido = (TH1D*)h_DiffXSec->Clone();
            tmpKido->SetLineColor(Kidoth1_Binned->GetLineColor());
            tmpKido->SetLineStyle(Kidoth1_Binned->GetLineStyle());
            tmpKido->SetLineWidth(Kidoth1_Binned->GetLineWidth());
            for (int i=0; i<(int)tmpKido->GetNbinsX()+2; i++){tmpKido->SetBinContent(i,Kidoth1_Binned->GetBinContent(Kidoth1_Binned->FindBin(tmpKido->GetBinCenter(i))));};
        }
        if(Ahrensth1_Binned){
            /// Idiot definition of temporary histogram for Ahrens due to the larger number of bins in raw histogram
            tmpAhrens = (TH1D*)h_DiffXSec->Clone();
            tmpAhrens->SetLineColor(Ahrensth1_Binned->GetLineColor());
            tmpAhrens->SetLineStyle(Ahrensth1_Binned->GetLineStyle());
            tmpAhrens->SetLineWidth(Ahrensth1_Binned->GetLineWidth());
            for (int i=0; i<(int)tmpAhrens->GetNbinsX()+2; i++){tmpAhrens->SetBinContent(i,Ahrensth1_Binned->GetBinContent(Ahrensth1_Binned->FindBin(tmpAhrens->GetBinCenter(i))));};
        }

        double yminRatio, ymaxRatio;
        setResultRatioRanges(yminRatio, ymaxRatio);
        common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat, ratio_tota, powheghistBinned, mcnlohistBinned, tmpKido, tmpAhrens,powhegHerwighistBinned, perugia11histBinned, yminRatio, ymaxRatio);
        c->Update();
        c->Modified();
    };

    if(drawNLOCurves && drawMadScaleMatching) common::drawRatioXSEC(h_DiffXSec,madgraphhistBinned,0,0,madupBinned,maddownBinned,matchupBinned,matchdownBinned,0,0,0.4, 1.6);
    if(drawNLOCurves && drawMadMass) common::drawRatioXSEC(madgraphhistBinned,h_DiffXSec,ratio_stat,0,madupBinned,maddownBinned,matchupBinned,matchdownBinned,match2upBinned,match2downBinned,0.4, 1.6);


    c->Print(outdir.Copy()+"DiffXS_"+name+".eps");
    //c->Print(outdir.Copy()+"DiffXS_"+name+".C");
    TFile out_source(outdir.Copy()+"DiffXS_"+name+"_source.root", "RECREATE");
    c->Write("canvas");
    tga_DiffXSecPlot->Write("data_staterror_only");
    tga_DiffXSecPlotwithSys->Write("data");
    h_GenDiffXSec->Write("mc");
    out_source.Close();
    delete c;
    gStyle->SetEndErrorSize(0);

    std::cout<<"-------------------------------------------------------------------"<<std::endl;
    std::cout<<"Starting the calculation of Chi2/ndof\n"<<std::endl;
    double chi2 = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, madgraphhistBinned);
    std::cout<<"The CHI2/ndof (vs Madgraph) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2<<std::endl;

    if(drawNLOCurves && drawPOWHEG && powheghistBinned && powheghistBinned->GetEntries()){
        double chi2Powheg = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, powheghistBinned);
        std::cout<<"The CHI2/ndof (vs Powheg+Pythia) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2Powheg<<std::endl;
    }
    if(drawNLOCurves && drawPOWHEGHERWIG && powhegHerwighistBinned && powhegHerwighistBinned->GetEntries()){
        double chi2PowhegHerwig = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, powhegHerwighistBinned);
        std::cout<<"The CHI2/ndof (vs Powheg+Herwig) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2PowhegHerwig<<std::endl;
    }
    if(drawNLOCurves && drawMCATNLO && mcnlohistBinned && mcnlohistBinned->GetEntries()){
        double chi2McAtNlo = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, mcnlohistBinned);
        std::cout<<"The CHI2/ndof (vs MC@NLO+Herwig) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2McAtNlo<<std::endl;
    }
    if(drawNLOCurves && drawKidonakis && Kidoth1_Binned && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame")){
        double chi2Kidonakis = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, Kidoth1_Binned);
        std::cout<<"The CHI2/ndof (vs Kidonakis) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2Kidonakis<<std::endl;
    }
    if(drawNLOCurves && drawMadScaleMatching){
        double chi2 = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, madupBinned);
        std::cout<<"The CHI2/ndof (vs SCALE UP) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2<<endl;
        chi2 = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, maddownBinned);
        std::cout<<"The CHI2/ndof (vs SCALE DOWN) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2<<endl;
        chi2 = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, matchupBinned);
        std::cout<<"The CHI2/ndof (vs MATCH UP) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2<<endl;
        chi2 = Plotter::GetChi2 (tga_DiffXSecPlotwithSys, matchdownBinned);
        std::cout<<"The CHI2/ndof (vs MATCH DOWN) value for '"<<name<<"' in channel '"<<subfolderChannel.Copy().Remove(0, 1)<<"' is "<<chi2<<endl;
    }
    std::cout<<"-------------------------------------------------------------------"<<std::endl;

    PrintResultTotxtFile(Channel, binCenters, tga_DiffXSecPlot, tga_DiffXSecPlotwithSys);

    TCanvas * c1 = new TCanvas("DiffXS","DiffXS");
    TH1* stacksum = common::summedStackHisto(stack);

    for(unsigned int i=1; i<hists.size() ; i++){ // sum all data plots to first histogram
        if(legends.at(i) == legends.at(0)){
            varhists[0]->Add(varhists[i]);
        }
    }
    TH1D* syshist =0;
    syshist = (TH1D*)stacksum->Clone();
    double lumierr = 0.026;
    //stat uncertainty::make a function 
    for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
        Double_t binc = 0;
        binc += stacksum->GetBinContent(i);
        syshist->SetBinContent(i, binc);
        // calculate uncertainty: lumi uncertainty
        Double_t binerr2 = binc*binc*lumierr*lumierr;
        Double_t topunc = 0; // uncertainty on top xsec

        double topxsecErr2 = 2.2*2.2 + 11.6*11.6;

        double topRelUnc =  TMath::Sqrt(topxsecErr2)/topxsec;
        //Functionality for multiple signal histograms
        topunc += varhists[signalHist]->GetBinContent(i)*topRelUnc;
        binerr2 += (topunc*topunc);
        syshist->SetLineColor(1);
        syshist->SetBinError(i, TMath::Sqrt(binerr2));
    }

    setControlPlotLegendStyle(varhistsPlotting, legends, leg);
    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);
    //leg->AddEntry( syshist, "Uncertainty", "f" );


    varhists[0]->SetMaximum(1.5*varhists[0]->GetBinContent(varhists[0]->GetMaximumBin()));

    varhists[0]->SetMinimum(0);
    varhists[0]->GetYaxis()->SetTitle("events");
    varhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    varhists[0]->Draw("e");

    //Add the binwidth to the yaxis in yield plots
    TString ytitle = TString(varhists[0]->GetYaxis()->GetTitle()).Copy();
    double binwidth = varhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binwidth;
    if(name.Contains("Rapidity") || name.Contains("Eta")){ytitle.Append(" / ").Append(width.str());}
    else if(name.Contains("pT") || name.Contains("Mass") || name.Contains("mass") || name.Contains("MET") || name.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
    varhists[0]->GetYaxis()->SetTitle(ytitle);

    stack->Draw("same HIST");

    //Only necessary if we want error bands

    /*    TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0.5)");//this is frustrating and stupid but apparently necessary...
    setex1->Draw();
    syshist->Draw("same,E2");
    TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.)");
    setex2->Draw();*/
    varhists[0]->Draw("same, e1"); //############
    //varhists[0]->Draw("same, e"); 
    DrawCMSLabels(0, 8);
    DrawDecayChLabel(channelLabel[channelType]);
    leg->Draw("SAME");
    gPad->RedrawAxis();

    c1->Print(outdir.Copy()+"preunfolded_"+name+".eps");
    TFile out_root(outdir.Copy()+"preunfolded_"+name+"_source.root", "RECREATE");

    varhists[0]->Write(name+"_data");
    stacksum->Write(name+"_allmc");
    c1->Write(name + "_canvas");
    out_root.Close();
    c1->Clear();
    delete c1;
    delete stacksum;
    for (unsigned int i =0; i<hists.size(); i++){
        delete varhists[i];
//         delete varhistsPlotting.at(i);
    }
}


void Plotter::PlotSingleDiffXSec(TString Channel, TString Systematic){

    setDataSet(Channel,Systematic);
    if (!fillHisto()) return;

    if(Channel=="ee"){channelType=0;}
    if(Channel=="mumu"){channelType=1;}
    if(Channel=="emu"){channelType=2;}
    if(Channel=="combined"){channelType=3;}

    TH1::AddDirectory(kFALSE);
    TGaxis::SetMaxDigits(2);

    double Xbins[XAxisbins.size()];
    double binCenters[XAxisbinCenters.size()];
    
    for(unsigned int i = 0; i<XAxisbins.size();i++)        {Xbins[i]=XAxisbins[i];}
    for(unsigned int i = 0; i<XAxisbinCenters.size();i++ ) {binCenters[i] = XAxisbins[i] + (XAxisbins[i+1]-XAxisbins[i])/2;}

    double DataSum[XAxisbinCenters.size()];
    double GenSignalSum[XAxisbinCenters.size()];
    double BGSum[XAxisbinCenters.size()];

    TH1 *varhists[hists.size()];
    TString newname = name;
    if(name.Contains("Hyp")){//Histogram naming convention has to be smarter
      newname.ReplaceAll("Hyp",3,"",0);
    }

    TString ftemp = "preunfolded/"+Systematic+"/"+Channel+"/"+name+"_UnfoldingHistos.root";
    TH1 *GenPlotTheory =  fileReader->GetClone<TH1D>(ftemp, "aGenHist");
    TH2 *genReco2d =  fileReader->GetClone<TH2D>(ftemp, "aRespHist");

    for (unsigned int i =0; i<hists.size(); i++){
        varhists[i]=hists[i].Rebin(bins,"varhists",Xbins);
        setStyle(varhists[i], i);
    }

    std::unique_ptr<TH1> GenPlot { GenPlotTheory->Rebin(bins,"genplot",Xbins) };

    THStack * stack=  new THStack("def", "def");
    TLegend *leg = new TLegend();
    int legchange = 0;
    std::vector<TH1 *> varhistsPlotting;
    varhistsPlotting.resize(hists.size());


    for(unsigned int i=0; i<hists.size() ; i++){ // prepare histos and leg
        setStyle(varhists[i], i);
        varhistsPlotting[i]=(TH1*)varhists[i]->Clone();
        if(legends.at(i) != "Data"){
            if((legends.at(i) == DYEntry)){
                varhists[i]->Scale(DYScale.at(channelType));
                varhistsPlotting[i]->Scale(DYScale.at(channelType));
            }

            if(i!=(hists.size()-1)){
                if(legends.at(i)!=legends.at(i+1)){varhistsPlotting[i]->SetLineColor(1);}
            }else{
                varhistsPlotting[i]->SetLineColor(1);
            }

            if(legends.at(i) != legends.at(i-1)){
                varhistsPlotting[i]->SetLineColor(1);
                stack->Add(varhistsPlotting[i]);
            }
            if(i > 1){
                if(legends.at(i) != legends.at(i-1)){
                    legchange = i;
                    if( (legends.at(i) == DYEntry) && DYScale.at(channelType)!= 1){
                        leg->AddEntry(varhistsPlotting[i], legends.at(i), "f");
                    } else {
                        leg->AddEntry(varhistsPlotting[i], legends.at(i) ,"f");
                    }
                } else{
                    varhistsPlotting[legchange]->Add(varhistsPlotting[i]);
                }
            }
        } else{
            if(i==0) leg->AddEntry(varhistsPlotting[i], legends.at(i) ,"pe");
        }
    }

    ///////////////////////////////////
    //purity and stability plots as taken from CombinedCrossSection...

    TH1* genHist = (TH1*)GenPlot->Clone();
    TH1* genRecHist = new TH1D("","",bins,Xbins);
    int intbinsX[XAxisbins.size()];
    int intbinsY[XAxisbins.size()];

    // fill the elements of the main diagonal of the 2d hist into binned 1D histogram
    for (unsigned int i=0; i<XAxisbins.size(); ++i) {
        intbinsX[i] = genReco2d->GetXaxis()->FindBin(Xbins[i]+0.001);
        intbinsY[i] = genReco2d->GetYaxis()->FindBin(Xbins[i]+0.001);

        if (i>0) genRecHist->SetBinContent(i,
                                          ((TH2D*)genReco2d)->Integral( intbinsX[i-1],intbinsX[i]-1,
                                                                        intbinsY[i-1],intbinsY[i]-1)
                                          );

    }

    TH1* genPseHist = ((TH2D*)genReco2d)->ProjectionY();
    TH1* recPseHist = ((TH2D*)genReco2d)->ProjectionX();
    
    TH1* genBinHist    = genPseHist->Rebin(bins,"genBinHist", Xbins);
    TH1* recBinHist    = recPseHist->Rebin(bins,"recBinHist", Xbins);

    genRecHist->SetBinContent(0,      0);
    genRecHist->SetBinContent(bins+1,0);
    genBinHist->SetBinContent(0,      0);
    genBinHist->SetBinContent(bins+1,0);
    recBinHist->SetBinContent(0,      0);
    recBinHist->SetBinContent(bins+1,0);
    genHist   ->SetBinContent(0,      0);
    genHist   ->SetBinContent(bins+1,0);

    // this is realy ugly but necessary:
    // As it seems, somewhere a double is tranformed into a float so that
    // efficiencies can be larger than 1.
    for(Int_t i=1; i<=genRecHist->GetNbinsX(); ++i){
      if(genRecHist->GetBinContent(i) > recBinHist->GetBinContent(i)){
        genRecHist->SetBinContent(i,recBinHist->GetBinContent(i));
        std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin" << i
                  << " = " << genRecHist->GetBinContent(i) << " is larger than number of reconstructed events in that bin"
                  << " = " << recBinHist->GetBinContent(i) << std::endl;
        }
      if(genRecHist->GetBinContent(i) > genBinHist->GetBinContent(i)){
        genRecHist->SetBinContent(i,genBinHist->GetBinContent(i));
        std::cout << "WARNING in PlotDifferentialCrossSections: number of events generated and reconstructed in bin " << i
                  << " is larger than number of generated events in that bin" << std::endl;
      }
    }

    // efficiency, purity, stability
    TGraphAsymmErrors* grE; // for efficiency
    TGraphAsymmErrors* grP; // for purity
    TGraphAsymmErrors* grS; // for stability

    // efficiency
    grE = new TGraphAsymmErrors(recBinHist, genHist);
    grE->SetMinimum(0);
    grE->SetMaximum(1);
    grE->SetLineColor(8);
    grE->SetLineWidth(2);
    grE->SetMarkerSize(2);
    grE->SetMarkerStyle(21);
    grE->SetMarkerColor(8);

    // purity
    grP = new TGraphAsymmErrors(genRecHist, recBinHist);
    grP->SetLineColor(4);
    grP->SetLineWidth(2);
    grP->SetMarkerSize(2);
    grP->SetMarkerStyle(23);
    grP->SetMarkerColor(4);

    // stability
    grS = new TGraphAsymmErrors(genRecHist, genBinHist);
    grS->SetLineColor(2);
    grS->SetLineWidth(2);
    grS->SetMarkerSize(2);
    grS->SetMarkerStyle(22);
    grS->SetMarkerColor(2);


    grE->GetXaxis()->SetTitle(XAxis);
    TCanvas * cESP = new TCanvas("ESP","ESP");

    // this is a dummy to get the x axis range corrct

    recBinHist->Reset();
    recBinHist->Draw();
    recBinHist->SetMaximum(1.);
    recBinHist->GetXaxis()->SetTitle(TString("Reconstructed ").Copy().Append(XAxis));
    recBinHist->GetXaxis()->SetNoExponent(kTRUE);
    grE->GetXaxis()->SetNoExponent(kTRUE);
    grE->Draw("P,SAME");
    grP->Draw("P,SAME");
    grS->Draw("P,SAME");
    TLegend* leg3 = new TLegend();
    leg3->AddEntry(grE, "Efficiency", "p" );
    leg3->AddEntry(grP, "Purity",    "p" );
    leg3->AddEntry(grS, "Stability", "p" );
    setResultLegendStyle(leg3, 0);
    leg3->Draw("SAME");

    TString outdir = common::assignFolder(outpathPlots, Channel, Systematic);
    cESP->Print(outdir.Copy()+"ESP_"+name+".eps");
    cESP->Clear();
    delete cESP;

    bool init = false;
    for (unsigned int hist =0; hist<hists.size(); hist++){
        if(legends.at(hist) == "Data"){
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                DataSum[bin]+=varhists[hist]->GetBinContent(bin+1);
            }
        }
        else if((legends.at(hist) == "t#bar{t} Signal")&&init==false){
            signalHist=hist;
            init=true;
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                GenSignalSum[bin] += GenPlot->GetBinContent(bin+1);
            }
        }
        else{
            for (Int_t bin=0; bin<bins; ++bin) {//poor for loop placement, but needed because genplot is the sum of all signal histograms
                BGSum[bin]+=varhists[hist]->GetBinContent(bin+1);
            }
        }
    }
    double totalDataSum = 0;
    double GenDiffXSecPlot[XAxisbinCenters.size()];
    for (Int_t bin=0; bin<bins; ++bin) {
        totalDataSum+=DataSum[bin];
    }

    TH1 *h_DiffXSec    = (TH1D*)varhists[0]->Clone(); h_DiffXSec->Reset();
    TH1 *h_GenDiffXSec = (TH1D*)varhists[0]->Clone(); h_GenDiffXSec->Reset();

    double DiffXSecPlot[XAxisbinCenters.size()];
    double DiffXSecStatErrorPlot[XAxisbinCenters.size()];

    TString Dummy;
    //Read central and absolute statistical uncertainty values from Nominal
    ifstream ResultsList("UnfoldingResults/"+Systematic+"/"+Channel+"/"+name+"Results.txt");
    if(!ResultsList.is_open())
    {
        std::cout<<"WARNING (in PlotSingleDiffXSec): File is not open.\nFix this. \nEXITING!!"<<std::endl;
        exit(123);
    }
    for (Int_t bin=0; bin<bins; bin++){//Retrieve arrays for plotting
        ResultsList>>Dummy>>XAxisbinCenters[bin]>>Dummy>>Xbins[bin]>>Dummy>>Xbins[bin+1]>>Dummy>>DiffXSecPlot[bin]>>Dummy>>DiffXSecStatErrorPlot[bin]>>Dummy>>GenDiffXSecPlot[bin];
        h_DiffXSec->SetBinContent(bin+1,DiffXSecPlot[bin]);
        h_DiffXSec->SetBinError(bin+1,DiffXSecStatErrorPlot[bin]);
        h_GenDiffXSec->SetBinContent(bin+1,GenDiffXSecPlot[bin]);
    }

    double TotalVisXSection = 1.; //this can currently be set to 1. because David's code takes care of the normalization, but just incase we need it

    h_DiffXSec->Scale(1/TotalVisXSection);

    Double_t mexl[XAxisbinCenters.size()];
    Double_t mexh[XAxisbinCenters.size()];
    for (unsigned int j=0; j<XAxisbinCenters.size();j++){mexl[j]=0;mexh[j]=0;}
    TGraphAsymmErrors *tga_DiffXSecPlot = new TGraphAsymmErrors(bins, binCenters, DiffXSecPlot, mexl, mexh, DiffXSecStatErrorPlot, DiffXSecStatErrorPlot);
    tga_DiffXSecPlot->SetMarkerStyle(1);
    //tga_DiffXSecPlot->SetMarkerStyle(20);
    tga_DiffXSecPlot->SetMarkerColor(kBlack);
    tga_DiffXSecPlot->SetMarkerSize(1);
    tga_DiffXSecPlot->SetLineWidth(2);
    tga_DiffXSecPlot->SetLineColor(kBlack);

    GenPlot->Scale(topxsec/(SignalEventswithWeight*GenPlot->GetBinWidth(1)));
    if( (name.Contains("Lepton") || name.Contains("Top") || name.Contains("BJet")) &&
        !name.Contains("Lead")) {
            GenPlotTheory->Scale(1./2.);
    }

    GenPlotTheory->Scale(1./GenPlotTheory->Integral("width"));
    h_GenDiffXSec->Scale(1./h_GenDiffXSec->Integral("width"));

    TH1* madgraphhist = 0, *mcnlohist=0, *mcnlohistup=0, *mcnlohistdown=0, *powheghist=0, *powhegHerwighist=0, *perugia11hist = 0;
    TH1* mcnlohistnorm=0;
    TGraph *mcatnloBand=0;

    TH1* madgraphhistBinned = 0, *mcnlohistnormBinned = 0, *mcnlohistupBinned = 0;
    TH1* mcnlohistdownBinned = 0, *mcnlohistBinned = 0;
    TH1* powheghistBinned = 0, *powhegHerwighistBinned = 0;
    TH1* perugia11histBinned = 0;

    TH1 *Kidoth1_Binned = 0, *Ahrensth1_Binned = 0;

    TH1* madup=0, *maddown=0, *matchup=0, *matchdown=0, *match2up = 0, *match2down = 0;
    TH1* madupBinned = 0, *maddownBinned = 0, *matchupBinned = 0, *matchdownBinned = 0, *match2upBinned = 0, *match2downBinned = 0;

    bool canDrawMCATNLO = true;
    if (drawNLOCurves && drawMCATNLO) {
        mcnlohist = GetNloCurve(newname,"MCATNLO");
        double mcnloscale = 1./mcnlohist->Integral("width");
        mcnlohist->Scale(mcnloscale);

        mcnlohistBinned = mcnlohist->Rebin(bins,"mcnloplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            mcnlohistBinned->SetBinContent(bin+1,mcnlohistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohist->GetBinWidth(1)));
            mcnlohistBinned->SetBinError(bin+1, 1e-10);
        }
        mcnlohistBinned->Scale(1./mcnlohistBinned->Integral("width"));

        if(name.Contains("LeptonpT")){mcnlohistnorm = GetNloCurve("Leptons","Pt","MCatNLO");}//temprorary until I change the naming convention in the root file
        else if(name.Contains("LeptonEta")){mcnlohistnorm = GetNloCurve("Leptons","Eta","MCatNLO");}
        else if(name.Contains("LLBarpT")){mcnlohistnorm = GetNloCurve("LepPair","Pt","MCatNLO");}
        else if(name.Contains("LLBarMass")){mcnlohistnorm = GetNloCurve("LepPair","Mass","MCatNLO");}
        else if(name.Contains("ToppT")){mcnlohistnorm = GetNloCurve("TopQuarks","Pt","MCatNLO");}
        else if(name.Contains("TopRapidity")){mcnlohistnorm = GetNloCurve("TopQuarks","Rapidity","MCatNLO");}
        else if(name.Contains("TTBarpT")){mcnlohistnorm = GetNloCurve("TtBar","Pt","MCatNLO");}
        else if(name.Contains("TTBarRapidity")){mcnlohistnorm = GetNloCurve("TtBar","Rapidity","MCatNLO");}
        else if(name.Contains("TTBarMass")){mcnlohistnorm = GetNloCurve("TtBar","Mass","MCatNLO");}
        else if(name.Contains("BJetpT")){mcnlohistnorm = GetNloCurve("Jets","Pt","MCatNLO");}
        else if(name.Contains("BJetEta")){mcnlohistnorm = GetNloCurve("Jets","Eta","MCatNLO");}

        if (mcnlohistnorm) {
            mcnlohistnormBinned = mcnlohistnorm->Rebin(bins,"genBinHistNorm", Xbins);

            if(name.Contains("LeptonpT")){mcnlohistup = GetNloCurve("Leptons","Pt","MCNLOup");}//temprorary until I change the naming convention in the root file
            else if(name.Contains("LeptonEta")){mcnlohistup = GetNloCurve("Leptons","Eta","MCNLOup");}
            else if(name.Contains("LLBarpT")){mcnlohistup = GetNloCurve("LepPair","Pt","MCNLOup");}
            else if(name.Contains("LLBarMass")){mcnlohistup = GetNloCurve("LepPair","Mass","MCNLOup");}
            else if(name.Contains("ToppT")){mcnlohistup = GetNloCurve("TopQuarks","Pt","MCNLOup");}
            else if(name.Contains("TopRapidity")){mcnlohistup = GetNloCurve("TopQuarks","Rapidity","MCNLOup");}
            else if(name.Contains("TTBarpT")){mcnlohistup = GetNloCurve("TtBar","Pt","MCNLOup");}
            else if(name.Contains("TTBarRapidity")){mcnlohistup = GetNloCurve("TtBar","Rapidity","MCNLOup");}
            else if(name.Contains("TTBarMass")){mcnlohistup = GetNloCurve("TtBar","Mass","MCNLOup");}
            else if(name.Contains("BJetpT")){mcnlohistup = GetNloCurve("Jets","Pt","MCNLOup");}
            else if(name.Contains("BJetEta")){mcnlohistup = GetNloCurve("Jets","Eta","MCNLOup");}
            mcnlohistupBinned    = mcnlohistup->Rebin(bins,"genBinHist", Xbins);


            if(name.Contains("LeptonpT")){mcnlohistdown = GetNloCurve("Leptons","Pt","MCNLOdown");}//temprorary until I change the naming convention in the root file
            else if(name.Contains("LeptonEta")){mcnlohistdown = GetNloCurve("Leptons","Eta","MCNLOdown");}
            else if(name.Contains("LLBarpT")){mcnlohistdown = GetNloCurve("LepPair","Pt","MCNLOdown");}
            else if(name.Contains("LLBarMass")){mcnlohistdown = GetNloCurve("LepPair","Mass","MCNLOdown");}
            else if(name.Contains("ToppT")){mcnlohistdown = GetNloCurve("TopQuarks","Pt","MCNLOdown");}
            else if(name.Contains("TopRapidity")){mcnlohistdown = GetNloCurve("TopQuarks","Rapidity","MCNLOdown");}
            else if(name.Contains("TTBarpT")){mcnlohistdown = GetNloCurve("TtBar","Pt","MCNLOdown");}
            else if(name.Contains("TTBarRapidity")){mcnlohistdown = GetNloCurve("TtBar","Rapidity","MCNLOdown");}
            else if(name.Contains("TTBarMass")){mcnlohistdown = GetNloCurve("TtBar","Mass","MCNLOdown");}
            else if(name.Contains("BJetpT")){mcnlohistdown = GetNloCurve("Jets","Pt","MCNLOdown");}
            else if(name.Contains("BJetEta")){mcnlohistdown = GetNloCurve("Jets","Eta","MCNLOdown");}
            mcnlohistdownBinned    = mcnlohistdown->Rebin(bins,"genBinHist", Xbins);

            for (Int_t bin=0; bin<bins; bin++){
                mcnlohistupBinned->SetBinContent(bin+1,mcnlohistupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistup->GetBinWidth(1)));
                mcnlohistdownBinned->SetBinContent(bin+1,mcnlohistdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistdown->GetBinWidth(1)));
                mcnlohistnormBinned->SetBinContent(bin+1,mcnlohistnormBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/mcnlohistnorm->GetBinWidth(1)));
            }
            mcnlohistupBinned->Scale(1./mcnlohistnormBinned->Integral("width"));
            mcnlohistdownBinned->Scale(1./mcnlohistnormBinned->Integral("width"));
            mcnlohistnormBinned->Scale(1./mcnlohistnormBinned->Integral("width"));

            for (Int_t bin=0; bin<bins; bin++){
                mcnlohistupBinned->SetBinContent(bin+1,(mcnlohistupBinned->GetBinContent(bin+1)/mcnlohistnormBinned->GetBinContent(bin+1))*mcnlohistBinned->GetBinContent(bin+1));
                mcnlohistdownBinned->SetBinContent(bin+1,(mcnlohistdownBinned->GetBinContent(bin+1)/mcnlohistnormBinned->GetBinContent(bin+1))*mcnlohistBinned->GetBinContent(bin+1));
            }

            //Uncertainty band for MC@NLO
            Double_t x[bins];
            Double_t xband[2*bins];
            Double_t errup[bins];
            Double_t errdn[bins];
            Double_t errorband[2*bins];

            for( Int_t j = 0; j< bins; j++ ){
                x[j]=mcnlohistBinned->GetBinCenter(j+1);
                errup[j]=(mcnlohistupBinned->GetBinContent(j+1)/mcnlohistnormBinned->GetBinContent(j+1))*mcnlohistBinned->GetBinContent(j+1);
                errdn[j]=(mcnlohistdownBinned->GetBinContent(j+1)/mcnlohistnormBinned->GetBinContent(j+1))*mcnlohistBinned->GetBinContent(j+1);

                xband[j] = x[j];
                errorband[j] = errdn[j]; //lower band
                xband[2*bins-j-1] = x[j];
                errorband[2*bins-j-1] = errup[j]; //upper band
            }

            mcatnloBand = new TGraph(2*bins, xband, errorband);
            mcatnloBand->SetFillColor(kGray);
            mcatnloBand->SetFillStyle(1001);
            mcatnloBand->SetLineColor(kBlue);
            mcatnloBand->SetLineWidth(2);
            mcatnloBand->SetLineStyle(5);
            canDrawMCATNLO = false;
        } else {
            std::cout << "\n*************************\nMC@NLO Curve not available!\n**********************\n";
            canDrawMCATNLO = false;
        }
    }

    madgraphhist = GetNloCurve(newname, "Nominal");
    double madgraphscale = 1./madgraphhist->Integral("width");
    madgraphhist->Scale(madgraphscale);
    madgraphhistBinned = madgraphhist->Rebin(bins,"madgraphplot",Xbins);
    for (Int_t bin=0; bin<bins; bin++){
        madgraphhistBinned->SetBinContent(bin+1,madgraphhistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/madgraphhist->GetBinWidth(1)));
    }
    madgraphhistBinned->Scale(1./madgraphhistBinned->Integral("width"));

    if(drawNLOCurves && drawPOWHEG){
        powheghist = GetNloCurve(newname, "POWHEG");
        double powhegscale = 1./powheghist->Integral("width");
        powheghist->Scale(powhegscale);
        powheghistBinned = powheghist->Rebin(bins,"powhegplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            powheghistBinned->SetBinContent(bin+1,powheghistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/powheghist->GetBinWidth(1)));
            powheghistBinned->SetBinError(bin+1, 1e-10);
        }
        powheghistBinned->Scale(1./powheghistBinned->Integral("width"));
    }
    if(drawNLOCurves && drawPOWHEGHERWIG){
        powhegHerwighist = GetNloCurve(newname, "POWHEGHERWIG");
        double powhegHerwigscale = 1./powhegHerwighist->Integral("width");
        powhegHerwighist->Scale(powhegHerwigscale);
        powhegHerwighistBinned = powhegHerwighist->Rebin(bins,"powhegHerwigplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            powhegHerwighistBinned->SetBinContent(bin+1,powhegHerwighistBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/powhegHerwighist->GetBinWidth(1)));
            powhegHerwighistBinned->SetBinError(bin+1, 1e-10);
        }
        powhegHerwighistBinned->Scale(1./powhegHerwighistBinned->Integral("width"));
    }
    if(drawNLOCurves && drawPERUGIA11){
        perugia11hist = GetNloCurve(newname, "PERUGIA11");
        double perugia11histscale = 1./perugia11hist->Integral("width");
        perugia11hist->Scale(perugia11histscale);

        perugia11histBinned = perugia11hist->Rebin(bins,"perugia11plot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            perugia11histBinned->SetBinContent(bin+1,perugia11histBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/perugia11hist->GetBinWidth(1)));
            perugia11histBinned->SetBinError(bin+1, 1e-10);
        }
        perugia11histBinned->Scale(1./perugia11histBinned->Integral("width"));
    }
    if(drawNLOCurves && drawKidonakis &&
        (name== "HypToppT" || name == "HypTopRapidity") && !name.Contains("Lead") && !name.Contains("RestFrame")){
        TString kidoFile = ttbar::DATA_PATH_DILEPTONIC() + "/kidonakisNNLO_8TeV.root";
        if(name.Contains("ToppT"))             Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topPt");
        else if(name.Contains("TopRapidity"))  Kidoth1_Binned = fileReader->GetClone<TH1>(kidoFile, "topY");
        for(int iter=0; iter <= Kidoth1_Binned->GetNbinsX(); iter++) Kidoth1_Binned->SetBinError(iter, 1e-10);
        Kidoth1_Binned->Scale(1./Kidoth1_Binned->Integral("width"));
    }
    if(drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT")){
        TString ahrensFile = ttbar::DATA_PATH_DILEPTONIC() + "/ahrensNNLL_8TeV.root";
        if(name == "HypTTBarMass")    Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarM");
        else if(name == "HypTTBarpT") Ahrensth1_Binned = fileReader->GetClone<TH1>(ahrensFile, "ttbarPt");
        for(int iter=0; iter <= Ahrensth1_Binned->GetNbinsX(); iter++) Ahrensth1_Binned->SetBinError(iter, 1e-10);
        Ahrensth1_Binned->Scale(1./Ahrensth1_Binned->Integral("width"));
    } else {drawAhrens = 0;}

    if(drawMadScaleMatching){
        madup = GetNloCurve(newname,"SCALE_UP");
        madup->Scale(1./madup->Integral("width"));
        madupBinned = madup->Rebin(bins,"madupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
          madupBinned->SetBinContent(bin+1,madupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/madup->GetBinWidth(1)));
        }
        madupBinned->Scale(1./madupBinned->Integral("width"));

        maddown = GetNloCurve(newname,"SCALE_DOWN");
        double maddownscale = 1./maddown->Integral("width");
        maddown->Scale(maddownscale);
        maddownBinned = maddown->Rebin(bins,"maddownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            maddownBinned->SetBinContent(bin+1,maddownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/maddown->GetBinWidth(1)));
        }
        maddownBinned->Scale(1./maddownBinned->Integral("width"));

        matchup = GetNloCurve(newname,"MATCH_UP");
        matchup->Scale(1./matchup->Integral("width"));
        matchupBinned = matchup->Rebin(bins,"matchupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            matchupBinned->SetBinContent(bin+1,matchupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchup->GetBinWidth(1)));
        }
        matchupBinned->Scale(1./matchupBinned->Integral("width"));

        matchdown = GetNloCurve(newname,"MATCH_DOWN");
        matchdown->Scale(1./matchdown->Integral("width"));
        matchdownBinned = matchdown->Rebin(bins,"matchdownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
           matchdownBinned->SetBinContent(bin+1,matchdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchdown->GetBinWidth(1)));
        }
        matchdownBinned->Scale(1./matchdownBinned->Integral("width"));

    }

    if(drawNLOCurves && drawMadMass){
        madup = GetNloCurveMass(newname,"MASS_UP","173.5");
        madup->Scale(1./madup->Integral("width"));
        madupBinned = madup->Rebin(bins,"madupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
          madupBinned->SetBinContent(bin+1,madupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/madup->GetBinWidth(1)));
        }
        madupBinned->Scale(1./madupBinned->Integral("width"));

        maddown = GetNloCurveMass(newname,"MASS_DOWN","171.5");
        maddown->Scale(1./maddown->Integral("width"));
        maddownBinned = maddown->Rebin(bins,"maddownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            maddownBinned->SetBinContent(bin+1,maddownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/maddown->GetBinWidth(1)));
        }
        maddownBinned->Scale(1./maddownBinned->Integral("width"));

        matchup = GetNloCurveMass(newname,"MASS_UP","175.5");
        matchup->Scale(1./matchup->Integral("width"));
        matchupBinned = matchup->Rebin(bins,"matchupplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            matchupBinned->SetBinContent(bin+1,matchupBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchup->GetBinWidth(1)));
        }
        matchupBinned->Scale(1./matchupBinned->Integral("width"));

        matchdown = GetNloCurveMass(newname,"MASS_DOWN","169.5");
        matchdown->Scale(1./matchdown->Integral("width"));
        matchdownBinned = matchdown->Rebin(bins,"matchdownplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
           matchdownBinned->SetBinContent(bin+1,matchdownBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/matchdown->GetBinWidth(1)));
        }
        matchdownBinned->Scale(1./matchdownBinned->Integral("width"));

        match2up = GetNloCurveMass(newname,"MASS_UP","178.5");
        double match2scale = 1./match2up->Integral("width");
        match2up->Scale(match2scale);
        match2upBinned = match2up->Rebin(bins,"match2upplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
            match2upBinned->SetBinContent(bin+1,match2upBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/match2up->GetBinWidth(1)));
        }
        match2upBinned->Scale(1./match2upBinned->Integral("width"));

        match2down = GetNloCurveMass(newname,"MASS_DOWN","166.5");
        double match2downscale = 1./match2down->Integral("width");
        match2down->Scale(match2downscale);
        match2downBinned = match2down->Rebin(bins,"match2downplot",Xbins);
        for (Int_t bin=0; bin<bins; bin++){//condense matrices to arrays for plotting
           match2downBinned->SetBinContent(bin+1,match2downBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/match2down->GetBinWidth(1)));
        }
        match2downBinned->Scale(1./match2downBinned->Integral("width"));

    }

    TCanvas * c = new TCanvas("DiffXS","DiffXS");
    if(logY){
      c->SetLogy();
    }
    h_DiffXSec->SetMarkerStyle(20);
    if(ymax!=0){
        if(logY){
            madgraphhistBinned->SetMaximum(18*madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));
        } else{
            madgraphhistBinned->SetMaximum(1.5*madgraphhistBinned->GetBinContent(madgraphhistBinned->GetMaximumBin()));
        }
    }
    madgraphhistBinned->GetXaxis()->SetNoExponent(kTRUE);
    if (name.Contains("Rapidity") || name.Contains("Eta")){madgraphhistBinned->GetYaxis()->SetNoExponent(kTRUE);}
    if (ymax!=0) madgraphhistBinned->SetMaximum(ymax);
    if (ymin!=0) madgraphhistBinned->SetMinimum(ymin);

    gStyle->SetEndErrorSize(8);
    if (drawNLOCurves && drawMCATNLO && canDrawMCATNLO) {
        mcnlohistupBinned->SetFillColor(kGray);
        mcnlohistupBinned->SetLineColor(kGray);
        mcnlohistupBinned->Draw("same");
        mcnlohistdownBinned->SetLineColor(10);
        mcnlohistdownBinned->SetFillColor(10);
        mcnlohistdownBinned->Draw("same");
    }
    GenPlotTheory->SetLineColor(kRed+1);
    GenPlotTheory->SetLineWidth(2);
    GenPlotTheory->SetLineStyle(1);
    h_GenDiffXSec->SetLineColor(GenPlotTheory->GetLineColor());
    h_GenDiffXSec->SetLineStyle(GenPlotTheory->GetLineStyle());

    // Plot statistical and totl error of full result
    TGraphAsymmErrors *ratio_stat = 0;
    if(tga_DiffXSecPlot) ratio_stat = (TGraphAsymmErrors*)tga_DiffXSecPlot->Clone("ratio_stat");

    if(ratio_stat){
        ratio_stat->SetFillStyle(1001);
        ratio_stat->SetFillColor(kGray+1);
        ratio_stat->SetLineColor(0);
        for (Int_t iter = 0; iter<tga_DiffXSecPlot->GetN(); iter++)
        {
            double binWidth = (XAxisbins[iter+1] - XAxisbins[iter])/2;
            double x = tga_DiffXSecPlot->GetX()[iter];
            double y_ratio_stat = tga_DiffXSecPlot->GetY()[iter] / tga_DiffXSecPlot->GetY()[iter];
            double abserr_ratio_stat = y_ratio_stat - std::abs(tga_DiffXSecPlot->GetErrorY(iter) - tga_DiffXSecPlot->GetY()[iter]) / tga_DiffXSecPlot->GetY()[iter];
            ratio_stat->SetPoint(iter, x, y_ratio_stat);
            ratio_stat->SetPointError(iter, binWidth, binWidth, abserr_ratio_stat, abserr_ratio_stat);
        }
    }

    TLegend *leg2 = new TLegend();
    leg2->SetHeader(Systematic);
    if(doClosureTest){
        leg2->AddEntry(h_DiffXSec, "Pseudo-Data", "p");
    } else {
        leg2->AddEntry(h_DiffXSec, "Data", "p");
    }

    setTheoryStyleAndFillLegend(h_DiffXSec, "data");
    setTheoryStyleAndFillLegend(madgraphhist, "madgraph");
    setTheoryStyleAndFillLegend(madgraphhistBinned, "madgraph", leg2);
    madgraphhistBinned->GetXaxis()->SetTitle(varhists[0]->GetXaxis()->GetTitle());
    madgraphhistBinned->GetYaxis()->SetTitle(varhists[0]->GetYaxis()->GetTitle());
    madgraphhistBinned->Draw();

    TH1* realTruthBinned = nullptr;
    if(doClosureTest && dataset.at(0).Contains("fakerun") && Channel != "combined")
    {
        newname.ReplaceAll("Hyp", "VisGen");
        newname.Prepend("VisGen");
        TH1* realTruth = fileReader->GetClone<TH1>(dataset.at(0), newname);
        if (!newname.Contains("Lead") && (newname.Contains("Lepton") || newname.Contains("Top") || newname.Contains("BJet")))
        {
            //loop over anti-particle histograms and add them
            TString antiName = newname.ReplaceAll("VisGen", "VisGenAnti");
            realTruth->Add(fileReader->Get<TH1>(dataset.at(0), antiName));
        }
        // Rebin histogram to final binning
        realTruth->Scale(1./realTruth->Integral("width"));
        realTruthBinned = realTruth->Rebin(bins,"realTruthBinned",Xbins);
        for (Int_t bin=0; bin<bins; bin++){
            realTruthBinned->SetBinContent(bin+1,realTruthBinned->GetBinContent(bin+1)/((Xbins[bin+1]-Xbins[bin])/realTruth->GetBinWidth(1)));
        }
        realTruthBinned->Scale(1./realTruthBinned->Integral("width"));

        // Histogram style
        realTruthBinned->SetLineColor(kRed-7);
        realTruthBinned->SetLineStyle(2);
        realTruthBinned->Draw("SAME][");
        leg2->AddEntry(realTruthBinned, "Simu. Reweighted", "l");
    }

    if (drawNLOCurves && drawMCATNLO) {
        setTheoryStyleAndFillLegend(mcnlohist, "mcatnloherwig");
        setTheoryStyleAndFillLegend(mcnlohistBinned, "mcatnloherwig", leg2);
        mcnlohistBinned->Draw("SAME,histo");
    }
    if(drawNLOCurves && drawPOWHEGHERWIG){
        setTheoryStyleAndFillLegend(powhegHerwighist, "powhegherwig");
        setTheoryStyleAndFillLegend(powhegHerwighistBinned, "powhegherwig", leg2);
        powhegHerwighistBinned->Draw("SAME,histo");
    }
    if(drawNLOCurves && drawPOWHEG){
        setTheoryStyleAndFillLegend(powheghist, "powhegpythia");
        setTheoryStyleAndFillLegend(powheghistBinned, "powhegpythia", leg2);
        powheghistBinned->Draw("SAME,histo");
    }
    if(drawNLOCurves && drawPERUGIA11){
        setTheoryStyleAndFillLegend(perugia11hist, "perugia11");
        setTheoryStyleAndFillLegend(perugia11histBinned, "perugia11", leg2);
        perugia11histBinned->Draw("SAME,histo");
    }
    if(drawNLOCurves && drawKidonakis &&
        (name== "HypToppT" || name == "HypTopRapidity") && 
        !name.Contains("Lead") && !name.Contains("RestFrame")){
        setTheoryStyleAndFillLegend(Kidoth1_Binned, "kidonakis", leg2);
        Kidoth1_Binned->Draw("SAME,histo");
    }
    if(drawNLOCurves && drawAhrens && (name == "HypTTBarMass" || name == "HypTTBarpT"))
    {
        setTheoryStyleAndFillLegend(Ahrensth1_Binned, "ahrens", leg2);
        Ahrensth1_Binned->Draw("SAME,histo");
    }
    if(drawNLOCurves && (drawMadScaleMatching || drawMadMass)){
        if(drawMadScaleMatching){
            setTheoryStyleAndFillLegend(madupBinned,    "scaleup",   leg2);
            setTheoryStyleAndFillLegend(maddownBinned,  "scaledown", leg2);
            setTheoryStyleAndFillLegend(matchupBinned,  "matchup",   leg2);
            setTheoryStyleAndFillLegend(matchdownBinned,"matchdown", leg2);
        }
        if(drawMadMass){
            setTheoryStyleAndFillLegend(madupBinned,    "mass173.5", leg2);
            setTheoryStyleAndFillLegend(maddownBinned,  "mass171.5", leg2);
            setTheoryStyleAndFillLegend(matchupBinned,  "mass175.5", leg2);
            setTheoryStyleAndFillLegend(matchdownBinned,"mass169.5", leg2);
            setTheoryStyleAndFillLegend(match2upBinned, "mass178.5", leg2);
            setTheoryStyleAndFillLegend(match2downBinned,"mass166.5", leg2);
        }
        madupBinned->Draw("SAME");
        maddownBinned->Draw("SAME");
        matchupBinned->Draw("SAME");
        matchdownBinned->Draw("SAME");
        match2upBinned->Draw("SAME");
        match2downBinned->Draw("SAME");

        ofstream OutputFileXSec(string("Plots/"+Channel+"/"+name+"DiffXsecMass.txt"));
        for(int i = 1; i < madupBinned->GetNbinsX(); i++){
        //OutputFileXSec<<"Nominal "<<"Mass 181 GeV" << " Mass 175 GeV"<< "Mass 169 GeV" << "Mass 163 GeV"<<endl;                             OutputFileXSec<<h_DiffXSec->GetBinContent(i)<< " "<<tga_DiffXSecPlot->GetErrorY(i-1)<<" "<<tga_DiffXSecPlotwithSys->GetErrorY(i-1)<< " "<<h|
        }
        OutputFileXSec.close();
    }

//     madgraphhistBinned->Draw("SAME");
    DrawCMSLabels(0, 8);
    DrawDecayChLabel(channelLabel[channelType]);

    if (drawNLOCurves) {
        //if (drawPOWHEG && powheghist->GetEntries())                                                     leg2->AddEntry(powheghistBinned, "tt+1jet m=178.5 GeV","l");
        //if (drawPOWHEGHERWIG && powhegHerwighist->GetEntries())                                         leg2->AddEntry(powhegHerwighistBinned, "tt+1jet m=172.5 GeV","l");
    }

    if(name.Contains("JetMult")) {
        TString legtit = "";
        if(name.Contains("pt60")) legtit += "p_{T}^{jet}> 60 GeV";
        else if(name.Contains("pt100")) legtit += "p_{T}^{jet}> 100 GeV";
        else legtit += "p_{T}^{jet}> 30 GeV";
        leg2->SetHeader(legtit);
    }

    setResultLegendStyle(leg2);
    leg2->Draw("same");

    if (drawNLOCurves && drawKidonakis &&  (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame")){
        DrawLabel("(arXiv:1210.7813)", leg2->GetX1NDC()+0.06, leg2->GetY1NDC()-0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }
    if (drawNLOCurves && drawAhrens &&  (name == "HypTTBarMass" || name == "HypTTBarpT")){
        DrawLabel("(arXiv:1003.5827)", leg2->GetX1NDC()+0.06, leg2->GetY1NDC()-0.025, leg2->GetX2NDC(), leg2->GetY1NDC(), 12, 0.025);
    }

//     madgraphhistBinned->Draw("same");
    gStyle->SetEndErrorSize(10);
    tga_DiffXSecPlot->Draw("p, SAME");
    ///Stupid clone to be able to see the stat error bars
    TGraphAsymmErrors *tga_DiffXSecPlotClone = (TGraphAsymmErrors*) tga_DiffXSecPlot->Clone();
    tga_DiffXSecPlotClone->SetMarkerStyle(20);
    tga_DiffXSecPlotClone->Draw("p, SAME");
    gPad->RedrawAxis();

    if(drawPlotRatio) {
        TH1D *tmpKido = 0, *tmpAhrens = 0;
        if(Kidoth1_Binned){
            /// Idiot definition of temporary histogram for Kidonakis due to the larger number of bins in raw histogram
            tmpKido = (TH1D*)h_DiffXSec->Clone();
            tmpKido->SetLineColor(Kidoth1_Binned->GetLineColor());
            tmpKido->SetLineStyle(Kidoth1_Binned->GetLineStyle());
            tmpKido->SetLineWidth(Kidoth1_Binned->GetLineWidth());
            for (int i=0; i<(int)tmpKido->GetNbinsX()+2; i++){tmpKido->SetBinContent(i,Kidoth1_Binned->GetBinContent(Kidoth1_Binned->FindBin(tmpKido->GetBinCenter(i))));};
        }
        if(Ahrensth1_Binned){
            /// Idiot definition of temporary histogram for Ahrens due to the larger number of bins in raw histogram
            tmpAhrens = (TH1D*)h_DiffXSec->Clone();
            tmpAhrens->SetLineColor(Ahrensth1_Binned->GetLineColor());
            tmpAhrens->SetLineStyle(Ahrensth1_Binned->GetLineStyle());
            tmpAhrens->SetLineWidth(Ahrensth1_Binned->GetLineWidth());
            for (int i=0; i<(int)tmpAhrens->GetNbinsX()+2; i++){tmpAhrens->SetBinContent(i,Ahrensth1_Binned->GetBinContent(Ahrensth1_Binned->FindBin(tmpAhrens->GetBinCenter(i))));};
        }

        double yminRatio, ymaxRatio;
        setResultRatioRanges(yminRatio, ymaxRatio);
        if(doClosureTest) { common::drawRatioXSEC(h_DiffXSec, realTruthBinned, ratio_stat,0,madgraphhistBinned, powheghistBinned,
                                                 mcnlohistBinned, tmpKido, tmpAhrens,powhegHerwighistBinned, 0.4, 1.6);
        } else { //common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat,0,powheghistBinned, mcnlohistBinned, tmpKido, tmpAhrens, powhegHerwighistBinned, perugia11histBinned,0.4, 1.6);
            if(drawNLOCurves && drawMadScaleMatching) common::drawRatioXSEC(h_DiffXSec,madgraphhistBinned,0,0,madupBinned,maddownBinned,matchupBinned,matchdownBinned,0,0,0.4, 1.6);
            else if(drawNLOCurves && drawMadMass) common::drawRatioXSEC(madgraphhistBinned,h_DiffXSec,0,0,madupBinned,maddownBinned,matchupBinned,matchdownBinned,match2upBinned,match2downBinned,0.4, 1.6);
            else common::drawRatioXSEC(h_DiffXSec, madgraphhistBinned, ratio_stat, 0, powheghistBinned, mcnlohistBinned, tmpKido, tmpAhrens,powhegHerwighistBinned, perugia11histBinned, yminRatio, ymaxRatio);
        };
    };
    c->Update();
    c->Modified();

    c->Print(outdir.Copy()+"DiffXS_"+name+".eps");
    TFile out_source(outdir.Copy()+"DiffXS_"+name+"_source.root", "RECREATE");
    c->Write("canvas");
    tga_DiffXSecPlot->Write("data_staterror_only");
    h_GenDiffXSec->Write("mc");
    out_source.Close();
    delete c;
    gStyle->SetEndErrorSize(0);

//    double result_Integral = Plotter::CalculateIntegral(tga_DiffXSecPlot, Xbins);
    std::cout<<"-------------------------------------------------------------------"<<std::endl;
    std::cout<<"Starting the calculation of Chi2/ndof\n"<<std::endl;
    double chi2 = Plotter::GetChi2 (tga_DiffXSecPlot, madgraphhistBinned);
    std::cout<<"The CHI2/ndof (vs Madgraph) value for '"<<name<<"' in channel '"<<Channel<<"' is "<<chi2<<std::endl;

    if(drawPOWHEG && powheghistBinned && powheghistBinned->GetEntries()){
        double chi2Powheg = Plotter::GetChi2 (tga_DiffXSecPlot, powheghistBinned);
        std::cout<<"The CHI2/ndof (vs Powheg+Pythia) value for '"<<name<<"' in channel '"<<Channel<<"' is "<<chi2Powheg<<std::endl;
    }
    if(drawPOWHEGHERWIG && powhegHerwighistBinned && powhegHerwighistBinned->GetEntries()){
        double chi2PowhegHerwig = Plotter::GetChi2 (tga_DiffXSecPlot, powhegHerwighistBinned);
        std::cout<<"The CHI2/ndof (vs Powheg+Herwig) value for '"<<name<<"' in channel '"<<Channel<<"' is "<<chi2PowhegHerwig<<std::endl;
    }
    if(drawMCATNLO && mcnlohistBinned && mcnlohistBinned->GetEntries()){
        double chi2McAtNlo = Plotter::GetChi2 (tga_DiffXSecPlot, mcnlohistBinned);
        std::cout<<"The CHI2/ndof (vs MC@NLO+Herwig) value for '"<<name<<"' in channel '"<<Channel<<"' is "<<chi2McAtNlo<<std::endl;
    }
    if(drawPERUGIA11 && perugia11histBinned && perugia11histBinned->GetEntries()){
        double chi2perugia11 = Plotter::GetChi2 (tga_DiffXSecPlot, perugia11histBinned);
        std::cout<<"The CHI2/ndof (vs MadGraph+Pythia (Perugia11 tune)) value for '"<<name<<"' in channel '"<<Channel<<"' is "<<chi2perugia11<<std::endl;
    }
    if(drawKidonakis && Kidoth1_Binned && (name.Contains("ToppT") || name.Contains("TopRapidity")) && !name.Contains("Lead") && !name.Contains("RestFrame")){
        double chi2Kidonakis = Plotter::GetChi2 (tga_DiffXSecPlot, Kidoth1_Binned);
        std::cout<<"The CHI2/ndof (vs Kidonakis) value for '"<<name<<"' in channel '"<<Channel<<"' is "<<chi2Kidonakis<<std::endl;
    }
    std::cout<<"-------------------------------------------------------------------"<<std::endl;

    TCanvas * c1 = new TCanvas("DiffXS","DiffXS");
    TH1* stacksum = common::summedStackHisto(stack);

    for(unsigned int i=1; i<hists.size() ; i++){ // sum all data plots to first histogram
        if(legends.at(i) == legends.at(0)){
            varhists[0]->Add(varhists[i]);
        }
    }
    TH1D* syshist =0;
    syshist = (TH1D*)stacksum->Clone();
    double lumierr = 0.045;
    //stat uncertainty::make a function 
    for(Int_t i=0; i<=syshist->GetNbinsX(); ++i){
        Double_t binc = 0;
        binc += stacksum->GetBinContent(i);
        syshist->SetBinContent(i, binc);
        // calculate uncertainty: lumi uncertainty
        Double_t binerr2 = binc*binc*lumierr*lumierr;
        Double_t topunc = 0; // uncertainty on top xsec

        double topxsecErr2 = 2.2*2.2 + 11.6*11.6;

        double topRelUnc =  TMath::Sqrt(topxsecErr2)/topxsec;
        //Functionality for multiple signal histograms
        topunc += varhists[signalHist]->GetBinContent(i)*topRelUnc;
        binerr2 += (topunc*topunc);
        syshist->SetLineColor(1);
        syshist->SetBinError(i, TMath::Sqrt(binerr2));
    }

    syshist->SetFillStyle(3004);
    syshist->SetFillColor(kBlack);
    //leg->AddEntry( syshist, "Uncertainty", "f" );


    varhists[0]->SetMaximum(1.5*varhists[0]->GetBinContent(varhists[0]->GetMaximumBin()));

    varhists[0]->SetMinimum(0);
    varhists[0]->GetYaxis()->SetTitle("events");
    varhists[0]->GetXaxis()->SetNoExponent(kTRUE);
    varhists[0]->Draw("e");

    //Add the binwidth to the yaxis in yield plots
    TString ytitle = TString(varhists[0]->GetYaxis()->GetTitle()).Copy();
    double binwidth = varhists[0]->GetXaxis()->GetBinWidth(1);
    std::ostringstream width;
    width<<binwidth;
    if(name.Contains("Rapidity") || name.Contains("Eta")){ytitle.Append(" / ").Append(width.str());}
    else if(name.Contains("pT") || name.Contains("Mass") || name.Contains("mass") || name.Contains("MET") || name.Contains("HT")){ytitle.Append(" / ").Append(width.str()).Append(" GeV");};
    varhists[0]->GetYaxis()->SetTitle(ytitle);

    stack->Draw("same HIST");

    //Only necessary if we want error bands

    varhists[0]->Draw("same, e1");
    DrawCMSLabels(0, 8);
    DrawDecayChLabel(channelLabel[channelType]);
    setControlPlotLegendStyle(varhistsPlotting, legends, leg);
    leg->Draw("SAME");
    gPad->RedrawAxis();

    c1->Print(outdir.Copy()+"preunfolded_"+name+".eps");
    TFile out_root(outdir.Copy()+"preunfolded_"+name+"_source.root", "RECREATE");

    varhists[0]->Write(name+"_data");
    stacksum->Write(name+"_allmc");
    c1->Write(name + "_canvas");
    out_root.Close();
    c1->Clear();
    delete c1;
    delete stacksum;
    for (unsigned int i =0; i<hists.size(); i++){
        delete varhists[i];
//         delete varhistsPlotting.at(i);
    }
}



// get generator cross section curve for NLO prediction
TH1* Plotter::GetNloCurve(const char *particle, const char *quantity, const char *generator){

    TH1::AddDirectory(kFALSE);
    TString histname;
    if (strcmp(particle, "TopQuarks") == 0 || strcmp(particle, "TtBar") == 0) {
        histname="total_";
    } else {
        histname="visible_";
    }
    histname.Append(particle);
    histname.Append(quantity);
    histname.Append("_");
    histname.Append(generator);
    
    
    TString filename;
    if(strcmp(generator, "Powheg")==0){filename = "selectionRoot/Nominal/emu/ttbarsignalplustau_powheg.root";}
    else if(strcmp(generator, "MCatNLO")==0){filename = ttbar::DATA_PATH_DILEPTONIC() + "/MCatNLO_status3_v20120729.root";}
    else if(strcmp(generator, "MCNLOup")==0){filename = ttbar::DATA_PATH_DILEPTONIC() + "/MCatNLO_Uncert_Up_status3_v20120729.root";}
    else if(strcmp(generator, "MCNLOdown")==0){filename = ttbar::DATA_PATH_DILEPTONIC() + "/MCatNLO_Uncert_Down_status3_v20120729.root";}
    
    TH1 *hist = fileReader->GetClone<TH1>(filename, histname, true);
    if (hist) {
        TH1* weight = fileReader->GetClone<TH1>(filename, TString("total_LeptonsPt_").Append(generator).Data(), true);
        if(!weight){
            std::cerr << "WARNING in GetNloCurve: histogram to extract original number of events could not be opened! No weighting applied!" << std::endl;
        }
        return hist;
    }  
    std::cerr << "WARNING in GetNloCurve: input file could not been opened! Returning dummy!" << std::endl;
    hist = new TH1D();
    return hist; //I'd prefer to return nullptr
}

TH1* Plotter::GetNloCurve(TString NewName, TString Generator){

    TString filename;
    if (Generator == "Nominal"){
        filename = "_ttbarsignalplustau.root";
    } else if (Generator == "MCATNLO") {
        filename = "_ttbarsignalplustau_mcatnlo.root";
    } else if (Generator == "POWHEG") {
        filename = "_ttbarsignalplustau_powheg.root";
    } else if (Generator == "POWHEGHERWIG") {
        filename = "_ttbarsignalplustau_powhegHerwig.root";
    } else if (Generator == "PERUGIA11") {
        filename = "_ttbarsignalplustau_Perugia11.root";
    } else if (Generator == "MATCH_UP") {
        filename = "_ttbarsignalplustau_matchingup.root";
    } else if (Generator == "MATCH_DOWN") {
        filename = "_ttbarsignalplustau_matchingdown.root";
    } else if (Generator == "SCALE_UP") {
        filename = "_ttbarsignalplustau_scaleup.root";
    } else if (Generator == "SCALE_DOWN") {
        filename = "_ttbarsignalplustau_scaledown.root";
    } else {
        std::cerr << "Unknown Generator!\n";
        std::exit(2);
    }
    
    const static std::vector<TString> channelName {"ee", "mumu", "emu"};
    std::vector<TString> files;
    assert(channelType >= 0);
    assert(channelType <= 3);
    for (int i = 0; i <= 2; ++i) {
        if (channelType == i || channelType == 3) {
            files.push_back("selectionRoot/"+Generator+"/"+channelName.at(i)+"/"+channelName.at(i) + filename);
            std::cout << "Getting NLO curve from: " << files.at(files.size()-1) << std::endl;
        }
    }
    
    TString histname("VisGen"+NewName);
    TH1* hist = fileReader->GetClone<TH1>(files.at(0), histname);
    for (size_t i = 1; i < files.size(); ++i) {
        hist->Add(fileReader->Get<TH1>(files.at(i), histname));
    }

    if (!NewName.Contains("Lead") 
        && (NewName.Contains("Lepton") || NewName.Contains("Top") || NewName.Contains("BJet")))
    {
        //loop over anti-particle histograms and add them
        TString antiName("VisGenAnti"+NewName);
        for (const auto& file : files) {
            hist->Add(fileReader->Get<TH1>(file, antiName));
        }
    }
    return hist;
}

TH1* Plotter::GetNloCurveMass(TString NewName, TString Generator, TString Mass){
    
    TString filename;
    if (Generator == "MASS_UP" && Mass == "173.5") {
        filename = "_ttbarsignalplustau_massup.root";
    } else if (Generator == "MASS_UP" && Mass == "175.5") {
        filename = "_ttbarsignalplustau_175_massup.root";
    } else if (Generator == "MASS_UP" && Mass == "178.5") {
        filename = "_ttbarsignalplustau_178_massup.root";
    } else if (Generator == "MASS_DOWN" && Mass == "171.5") {
        filename = "_ttbarsignalplustau_massdown.root";
    } else if (Generator == "MASS_DOWN" && Mass == "169.5") {
        filename = "_ttbarsignalplustau_169_massdown.root";
    } else if (Generator == "MASS_DOWN" && Mass == "166.5") {
        filename = "_ttbarsignalplustau_166_massdown.root";
    } else {
        std::cerr << "Unknown Generator!\n";
        std::exit(2);
    }

    const static std::vector<TString> channelName {"ee", "mumu", "emu"};
    std::vector<TString> files;
    assert(channelType >= 0);
    assert(channelType <= 3);
    for (int i = 0; i <= 2; ++i) {
        if (channelType == i || channelType == 3) {
            files.push_back("selectionRoot/"+Generator+"/"+channelName.at(i)+"/"+channelName.at(i) + filename);
            cout << "Getting NLO curve from: " << files.at(files.size()-1) << endl;
        }
    }

    TString histname("VisGen"+NewName);
    TH1* hist = fileReader->GetClone<TH1>(files.at(0), histname);
    for (size_t i = 1; i < files.size(); ++i) {
        hist->Add(fileReader->Get<TH1>(files.at(i), histname));
    }

    if (!NewName.Contains("Lead")
        && (NewName.Contains("Lepton") || NewName.Contains("Top") || NewName.Contains("BJet")))
    {
        //loop over anti-particle histograms and add them
        TString antiName("VisGenAnti"+NewName);
        for (const auto& file : files) {
            hist->Add(fileReader->Get<TH1>(file, antiName));
        }
    }
    return hist;
}


TH1F* Plotter::ConvertGraphToHisto(TGraphErrors *pGraph){
  // takes data from a graph, determines binning and fills data into histogram
  Int_t NPoints = pGraph->GetN();
  Double_t BinLimits[NPoints+1];
  // sort graph
  pGraph->Sort();
  // determine lower limit of histogram: half the distance to next point
  Double_t x0,x1,y;
  pGraph->GetPoint(0,x0,y);
  pGraph->GetPoint(1,x1,y);
  Double_t Distance = TMath::Abs(x0-x1);
  BinLimits[0] = x0 - Distance/2.;
  // now set upper limits for all the other points
  for (Int_t k = 0 ; k<NPoints-1;k++){
    pGraph->GetPoint(k,x0,y);
    pGraph->GetPoint(k+1,x1,y);
    Distance = TMath::Abs(x0-x1);
    BinLimits[k+1] = x0 + Distance/2.;}
  // for the last point set upper limit similar to first point:
  pGraph->GetPoint(NPoints-2,x0,y);
  pGraph->GetPoint(NPoints-1,x1,y);
  Distance = TMath::Abs(x0-x1);
  BinLimits[NPoints] = x1 + Distance/2.;
  // now we know the binning and can create the histogram:
  TString Name = "ConvertedHisto"; 
  // make name unique 
  Name+= rand();
  TH1F *ThisHist = new TH1F(Name,"Converted Histogram",NPoints,BinLimits);
  // now fill the histogram
  for (Int_t i = 0; i<pGraph->GetN();i++){
    Double_t x2,y2;
    pGraph->GetPoint(i,x2,y2);
    ThisHist->SetBinContent(i+1,y2);
  }
  return ThisHist;
}

double Plotter::GetChi2 (TGraphAsymmErrors *data, TH1 *mc){
    
    double chi2 = 0.0;
    
    for ( int i=0; i<(int)data->GetN(); ++i ) {
        if (data->GetErrorYhigh(i) == 0 || data->GetErrorYlow(i) == 0) {
            std::cout<<"When calculating the Chi2 the DATA TGraph has error 0 in bin "<<i<<std::endl;
            //exit(42);
            return 0;
        }
        double dataMinusMC = data->GetY()[i]-mc->GetBinContent(mc->FindBin(data->GetX()[i]));
        double dataError   = (data->GetErrorYhigh(i) + data->GetErrorYlow(i))/2;
        chi2 += dataMinusMC * dataMinusMC / (dataError * dataError);
    }
    return chi2/data->GetN();
}


//TH1F* Plotter::reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, const std::vector<double> &vecBinning, TString plotname, bool rescale=1){
TH1F* Plotter::reBinTH1FIrregularNewBinning(TH1F *histoOldBinning, TString plotname, bool rescale){
  //  This function rebins a histogram using a variable binning
  // 
  //  (1) It is not required to have an equidistant binning.
  //  (2) Any type of ROOT-histgramme can be used, so the template 
  //      arguments should be 
  //      (a) histoT = TH1D,   TH1F,  ....
  //      (b) varT   = double, float, ....
  //  
  //  modified quantities: none
  //  used functions:      none
  //  used enumerators:    none
  //  
  //  "histoOldBinning":   plot to be re-binned
  //  "vecBinning":        vector containing all bin edges 
  //                       from xaxis.min to xaxis.max
  //  "rescale":           rescale the rebinned histogramme
  //                       (applicable for cross-section e.g.) 
  //  std::cout << std::endl;
  //  std::cout << std::endl;
  //  std::cout << "asdfasdfasdfasdfasdf hallo david " << plotname << " " << rescale << std::endl;
  //  std::cout << "histoOldBinning = ";
  //  for ( int i = 0 ; i < histoOldBinning->GetXaxis()->GetNbins() + 1; i++ ) std::cout << " " << histoOldBinning->GetXaxis()->GetBinLowEdge(i+1);
  //  std::cout << std::endl;
  //  std::cout << std::endl;
  //  std::cout << std::endl;

  unsigned int vecIndex=0;

  // fill vector into array to use appropriate constructor of TH1-classes
  const double *binArray = XAxisbins.data();

  // create histo with new binning
  TH1F *histoNewBinning = new TH1F("histoNewBinning"+plotname,"histoNewBinning"+plotname,XAxisbins.size()-1,binArray);

  // fill contents of histoOldBinning into histoNewBinning and rescale if applicable
  for (vecIndex = 0; vecIndex < XAxisbins.size()-1; vecIndex++){

    double lowEdge      = XAxisbins[vecIndex]; 
    if (plotname=="topPt"&&vecIndex==0&&lowEdge==0.0) lowEdge+=10;  // adhoc fix to compensate for minimum top-Pt cut in NNLO curve
    double highEdge     = XAxisbins[vecIndex+1];
    double newBinWidth  = highEdge - lowEdge;
    double newBinCenter = 0.5*(highEdge+lowEdge);
    double binSum       = 0.0;  

    for (int j=1; j<histoOldBinning->GetNbinsX(); j++){

    double oldBin = histoOldBinning->GetBinCenter(j);

    if( (oldBin>=lowEdge) && (oldBin<highEdge) ){
        if (rescale) binSum+=histoOldBinning->GetBinContent(j) * histoOldBinning->GetBinWidth(j);
        else         binSum+=histoOldBinning->GetBinContent(j);
    }
    }

    if (rescale) histoNewBinning->Fill(newBinCenter,binSum/newBinWidth);
    else histoNewBinning->Fill(newBinCenter,binSum);
  }

  return (TH1F*)histoNewBinning->Clone();
}


void Plotter::setResultLegendStyle(TLegend *leg, const bool result)
{
    double x1 = 0.560, y1 = 0.655;
    double height = 0.175, width = 0.275;
    if(result){
        if(name.Contains("Eta")){
            x1 = 0.4;
            y1 = 0.39;
        }
        if(name.Contains("TopRapidity") || name.Contains("HypTTBarDeltaRapidity")){
            x1 = 0.4;
            y1 = 0.41;
            height += 0.045;
        }
        if(name.Contains("DeltaPhi")){
            x1 = 0.4;
            y1 = 0.56;
        }
        if(name.Contains("TTBarRapidity") || name.Contains("TTBarpT") || name.Contains("TTBarMass") || name.Contains("ToppT")){
            height += 0.045;
            y1 -= 0.045;
        }
    }
    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}


void Plotter::setControlPlotLegendStyle(std::vector< TH1* > drawhists, std::vector< TString > legends, TLegend* leg, TLegend *leg1, TLegend *leg2)
{
    //this appendix  to the legend aligns vertically all entries in the legend
    std::string  appendix = "                                                                                    #color[10]{#bar{MARTIN}_{GOERNER}}";

    //hardcoded ControlPlot legend
    std::vector<TString> OrderedLegends;
    OrderedLegends.push_back("Data");
    OrderedLegends.push_back("t#bar{t} Signal");
    OrderedLegends.push_back("t#bar{t} Other");
    OrderedLegends.push_back("Single Top");
    OrderedLegends.push_back("t#bar{t}+V");
    OrderedLegends.push_back("W+Jets");
    OrderedLegends.push_back("Z / #gamma* #rightarrow ee/#mu#mu");
    OrderedLegends.push_back("Z / #gamma* #rightarrow #tau#tau");
    OrderedLegends.push_back("Diboson");
    OrderedLegends.push_back("QCD Multijet");

    leg->Clear();
    if(leg1) leg1->Clear();
    if(leg2) leg2->Clear();
    for(size_t i=0; i<OrderedLegends.size(); ++i){
        for(size_t j=0; j<drawhists.size(); ++j){
            if (OrderedLegends.at(i) == legends.at(j)){
                if( OrderedLegends.at(i) == "Data"){
                    leg->AddEntry(drawhists.at(j), OrderedLegends.at(i)+appendix, "pe");
                    if(leg1)leg1->AddEntry(drawhists.at(j), OrderedLegends.at(i)+appendix, "pe");
                    break;
                }
                else{
                    leg->AddEntry(drawhists.at(j), OrderedLegends.at(i)+appendix, "f");
                    if (leg1 && i < 5) leg1->AddEntry(drawhists.at(j), OrderedLegends.at(i)+appendix, "f");
                    if (leg2 && i > 4) leg2->AddEntry(drawhists.at(j), OrderedLegends.at(i)+appendix, "f");
                    break;
                }
            }
        }
    }
    //coordinates for legend without splitting
    double x1 = 1.00 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25 + 0.005;
    double y1 = 1.00 - gStyle->GetPadTopMargin() - gStyle->GetTickLength() - 0.05 - 0.03 * leg->GetNRows();
    double x2 = 1.00 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.025;
    double y2 = 1.00 - gStyle->GetPadTopMargin() - 0.8 * gStyle->GetTickLength();

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x2);
    leg->SetY2NDC(y2);

    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextAlign(12);

    if (!leg1) return;
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextAlign(12);

    // Define shifts of legends
    double xShift=0.8*(x2-x1);
    double yShift=0.06;
    double y1mod=y1+0.45*(y2-y1);

    // coordinates for splitted legends
    leg1->SetX1NDC(x1 - xShift);
    leg1->SetY1NDC(y1mod - yShift);
    leg1->SetX2NDC(x2 - xShift);
    leg1->SetY2NDC(y2 - yShift);

    if(!leg2) return;
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextAlign(12);

    // coordinates for splitted legends
    leg2->SetX1NDC(x1 + 0.01);
    leg2->SetY1NDC(y1mod - yShift);
    leg2->SetX2NDC(x2 + 0.01);
    leg2->SetY2NDC(y2 - yShift);

}

void Plotter::DrawLabel(TString text, const double x1, const double y1, const double x2, const double y2, int centering, double textSize){
    //Function to add Kidonakis references to DiffXSection plots' legends
    TPaveText *label = new TPaveText(x1, y1, x2, y2, "br NDC");
    label->AddText(text);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(centering);
    label->Draw("same");
}



void Plotter::ApplyFlatWeights(TH1* varhists, const double weight){

    if(weight == 0) {std::cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<std::endl;}
    //if(weight >=1e3){std::cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<std::endl;}
    varhists->Scale(weight);
}

double Plotter::CalcLumiWeight(const TString& WhichSample){
    if (WhichSample.Contains("run")) return 1;
    double lumiWeight=0;
    if(WhichSample!=""){
        double XSection = SampleXSection(WhichSample);
        if(XSection <= 0.){
            std::cout<<"Sample XSection is <0. Can't calculate luminosity weight!! returning"<<std::endl;
            return 0;
        }
        
        //From 'filename' get the number of weighted (MC weights) event processed.
        const TH1 *h_NrOfEvts = fileReader->Get<TH1>(WhichSample, "weightedEvents");
        double NrOfEvts = h_NrOfEvts->GetBinContent(1);
        lumiWeight = lumi*XSection/NrOfEvts;
    }
    
    if (lumiWeight == 0) {
        std::cout << WhichSample << " has lumi weight 0\n";
    }
    return lumiWeight;
}


double Plotter::SampleXSection(const TString& filename){
    
    //MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/view/CMS/StandardModelCrossSectionsat8TeV
    //  AN-12/194    AN-12/228
    
    if(filename.Contains("run"))              {return 1;}
    else if(filename.Contains("FullLept"))    {return topxsec * 0.1049;}
    else if(filename.Contains("SemiLept"))    {return topxsec * 0.4380;}
    else if(filename.Contains("Hadronic"))    {return topxsec * 0.4570;}
    else if(filename.Contains("Perugia11") &&
        filename.Contains("signal"))          {return topxsec * 0.1049;}
    else if(filename.Contains("ttbar") && !filename.Contains("ttbarW") && 
        !filename.Contains("ttbarZ"))         {return topxsec;}
    else if(filename.Contains("single"))      {return 11.1;}
    else if(filename.Contains("ww"))          {return 54.838;}
    else if(filename.Contains("wz"))          {return 33.21;}
    else if(filename.Contains("zz"))          {return 17.654;}
    else if(filename.Contains("1050"))        {return 860.5;}
    else if(filename.Contains("50inf"))       {return 3532.8;}
    else if(filename.Contains("wtolnu"))      {return 36257.2;}
    else if(filename.Contains("qcdmu15"))     {return 3.640E8*3.7E-4;}
    else if(filename.Contains("qcdmu2030"))   {return 2.870E8*6.500E-3;}
    else if(filename.Contains("qcdmu3050"))   {return 6.609E7*12.20E-3;}
    else if(filename.Contains("qcdmu5080"))   {return 8.802E6*21.80E-3;}
    else if(filename.Contains("qcdmu80120"))  {return 1.024E6*39.50E-3;}
    else if(filename.Contains("qcdmu120170")) {return 1.578E5*47.30E-3;}
    else if(filename.Contains("qcdem2030"))   {return 2.886E8*10.10E-3;}
    else if(filename.Contains("qcdem3080"))   {return 7.433E7*62.10E-3;}
    else if(filename.Contains("qcdem80170"))  {return 1.191E6*153.9E-3;}
    else if(filename.Contains("qcdbcem2030")) {return 2.886E8*5.800E-4;}
    else if(filename.Contains("qcdbcem3080")) {return 7.424E7*2.250E-3;}
    else if(filename.Contains("qcdbcem80170")){return 1.191E6*10.90E-3;}
    else if(filename.Contains("ttbarW"))      {return 0.232;}
    else if(filename.Contains("ttbarZ"))      {return 0.2057;}
    else if(filename.Contains("ttgjets"))     {return 1.8;}
    
    return -1;
}

// Draw label for Decay Channel in upper left corner of plot
void Plotter::DrawDecayChLabel(TString decaychannel, double textSize) {

    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void Plotter::DrawCMSLabels(int cmsprelim, double energy, double textSize) {

    const char *text;
    if(cmsprelim ==2 ) {//Private work for PhDs students
        text = "Private Work, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else if (cmsprelim==1) {//CMS preliminary label
        text = "CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else {//CMS label
        text = "CMS, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    }
    
    TPaveText *label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
    label->SetY2NDC(1.0);
    label->SetTextFont(42);
    label->AddText(Form(text, lumi/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}


void Plotter::PrintResultTotxtFile (TString channel, double binCenters[], TGraphAsymmErrors *tga_DiffXSecPlot, TGraphAsymmErrors *tga_DiffXSecPlotwithSys){
    
    ofstream SavingFile;
    string filename = string("Plots/FinalResults/"+channel+"/"+name+"LaTeX.txt");
    SavingFile.open(filename.c_str(), ios_base::out);
    
    SavingFile<<"Variable: "<<name<<"   Channel: "<<channelLabel.at(channelType)<<std::endl;
    SavingFile<<"BinCenter & LowXbinEdge  &  HighXbinEdge  &   DiffXSec  &  StatError(\\%)  & SystError(\\%)  & TotalError(\\%) \\\\"<<std::endl;
    SavingFile<<"\\hline"<<std::endl;
    for(int i=0; i<(int)tga_DiffXSecPlot->GetN(); i++){
        double DiffXSec=tga_DiffXSecPlot->GetY()[i];
        double RelStatErr=0, RelSysErr=0, RelTotErr=0;
        if(DiffXSec!=0.0){
            RelStatErr = 100*(tga_DiffXSecPlot->GetErrorY(i))/DiffXSec;
            RelTotErr  = 100*(tga_DiffXSecPlotwithSys->GetErrorY(i))/DiffXSec;
            if(RelTotErr>=RelStatErr) RelSysErr = TMath::Sqrt(RelTotErr*RelTotErr - RelStatErr*RelStatErr);
        }
        SavingFile<<"$"<<setprecision(3)<<binCenters[i]<<"$ & $"<<XAxisbins.at(i)<<"$ to $"<<setprecision(3)<<XAxisbins.at(i+1)<<"$   &  ";
        SavingFile<<setprecision(5)<<DiffXSec<<"  &   "<<setprecision(3)<<RelStatErr<<" &    "<<setprecision(3)<<RelSysErr<<" &    ";
        SavingFile<<setprecision(3)<<RelTotErr<<" \\\\"<<std::endl;
    }
    SavingFile.close();
}



void Plotter::CalcUpDownDifference( TString Channel, TString Syst_Up, TString Syst_Down, TString Variable){
    
    ///Function to get the error of a certain systematic: Sqrt(0.5*(Err_Up*Err_Up + Err_Down*Err_Down))/Nominal
    /// This prescription is taken from David's SVD_DoUnfoldSys function (DilepSVDFunctions.cc)

    ifstream NominalFile, UpFile, DownFile;
    string nominalfilename = string ("UnfoldingResults/Nominal/"+Channel+"/"+Variable+"Results.txt");
    string upfilename = string ("UnfoldingResults/"+Syst_Up+"/"+Channel+"/"+Variable+"Results.txt");
    string downfilename = string ("UnfoldingResults/"+Syst_Down+"/"+Channel+"/"+Variable+"Results.txt");

    NominalFile.open(nominalfilename);
    UpFile.open(upfilename);
    DownFile.open(downfilename);

    std::vector<double> BinCenters, BinLowEdge, BinHigEdge;
    std::vector<double> RelativeError;
    TString Dummy = "";
    double dummy = 0;

    if (!NominalFile.is_open() || !UpFile.is_open() || !DownFile.is_open()){
        std::cout<<"Nominal: "<<nominalfilename<<std::endl;
        std::cout<<"Sys Up : "<<upfilename<<std::endl;
        std::cout<<"Sys Dow: "<<downfilename<<std::endl;
        std::cout<<"The input file cannot be opened. Exiting!!"<<std::endl;
        exit(433);
    }
    for (Int_t bin=0; bin<bins; bin++){

        double XAxisbinCenters_Nom = 0, XLowEdge_Nom = 0, XHigEdge_Nom = 0;
        double XAxisbinCenters_Up = 0, XLowEdge_Up = 0, XHigEdge_Up = 0;
        double XAxisbinCenters_Down = 0, XLowEdge_Down = 0, XHigEdge_Down = 0;
        double CentralValue_Nom=0, CentralValue_Up = 0, CentralValue_Down = 0;

        NominalFile>>Dummy>>XAxisbinCenters_Nom>>Dummy>>XLowEdge_Nom>>Dummy>>XHigEdge_Nom>>Dummy>>CentralValue_Nom>>Dummy>>dummy>>Dummy>>dummy;
        UpFile>>Dummy>>XAxisbinCenters_Up>>Dummy>>XLowEdge_Up>>Dummy>>XHigEdge_Up>>Dummy>>CentralValue_Up>>Dummy>>dummy>>Dummy>>dummy;
        DownFile>>Dummy>>XAxisbinCenters_Down>>Dummy>>XLowEdge_Down>>Dummy>>XHigEdge_Down>>Dummy>>CentralValue_Down>>Dummy>>dummy>>Dummy>>dummy;

        BinCenters.push_back(XAxisbinCenters_Up);
        BinLowEdge.push_back(XLowEdge_Up);
        BinHigEdge.push_back(XHigEdge_Up);

        if(CentralValue_Nom != 0){
            double up = std::fabs(CentralValue_Up - CentralValue_Nom);
            double down = std::fabs(CentralValue_Down - CentralValue_Nom);
            double rel_err = 0.5*(up + down)/CentralValue_Nom;

            RelativeError.push_back(rel_err);
        }
    }
    NominalFile.close(); UpFile.close(); DownFile.close();

    if(Syst_Up =="POWHEG" && Syst_Down == "MCATNLO")
    {
        Syst_Up = "HAD_";
    } else {
        Syst_Up.Remove(Syst_Up.Length()-2,2);
    }

    ofstream SystematicRelError (common::assignFolder("UnfoldingResults", Channel, Syst_Up)+Variable+"Results.txt");
    if(!SystematicRelError.is_open()){
        std::cout<<"The output file cannot be opened. Exiting!!"<<std::endl;
        exit(434);
    }
    for (int bin = 0; bin<(int)RelativeError.size(); bin++){
        SystematicRelError<<"XAxisbinCenters[bin]: "<<BinCenters.at(bin)<<" bin: "<<BinLowEdge.at(bin)<<" to "<<BinHigEdge.at(bin)<<" SystematicRelError: "<<RelativeError.at(bin)<<std::endl;
    }
    SystematicRelError.close();

}



double Plotter::CalculateIntegral(TGraphAsymmErrors *tga_DiffXSecPlot, double Xbins[])
{
    double tmp_integral = 0;
    for (int i =0; i<tga_DiffXSecPlot->GetN(); i++)
    {
        tmp_integral += tga_DiffXSecPlot->GetY()[i] * (Xbins[i+1]-Xbins[i]);
    }
    std::cout<<"Integral = "<<tmp_integral<<std::endl;
    return tmp_integral;
}



void Plotter::setTheoryStyleAndFillLegend(TH1* histo, TString theoryName, TLegend *leg){

    histo->GetXaxis()->SetTitleOffset(1.08);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetLabelOffset(0.007);
    histo->GetXaxis()->SetLabelSize(0.04);

    histo->GetYaxis()->SetTitleOffset(1.7);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelOffset(0.007);
    histo->GetYaxis()->SetLabelSize(0.04);

    histo->SetLineWidth(2);
    if(theoryName != "data"){
        histo->SetMarkerSize(0);
        histo->SetMarkerStyle(1);
    }

    if(theoryName == "madgraph"){
        histo->SetLineColor(kRed+1);
        histo->SetLineStyle(1);
        if(leg) leg->AddEntry(histo, "MadGraph+Pythia",  "l");
    }
    if(theoryName == "powhegpythia"){
        histo->SetLineColor(kGreen+1);
        histo->SetLineStyle(7);
        if(leg) leg->AddEntry(histo, "Powheg+Pythia",  "l");
    }
    if(theoryName == "powhegherwig"){
        histo->SetLineColor(kGreen+3);
        histo->SetLineStyle(4);
        if(leg) leg->AddEntry(histo, "Powheg+Herwig",  "l");
    }
    if(theoryName == "mcatnloherwig"){
        histo->SetLineColor(kBlue);
        histo->SetLineStyle(5);
        if(leg) leg->AddEntry(histo, "MC@NLO+Herwig",  "l");
    }
    if(theoryName == "ahrens"){
        histo->SetLineColor(kViolet-6);
        histo->SetLineStyle(6);
        if(leg) leg->AddEntry(histo, "NLO+NNLL",  "l");
    }
    if(theoryName == "kidonakis"){
        histo->SetLineColor(kViolet-6);
        histo->SetLineStyle(2);
        if(leg) leg->AddEntry(histo, "Approx. NNLO",  "l");
    }
    if(theoryName == "matchup" || theoryName == "mass175.5"){
        histo->SetLineStyle(7);
        histo->SetLineColor(kPink-7);
        if(leg && theoryName == "matchup") leg->AddEntry(histo,"Matching up",  "l");
        if(leg && theoryName == "mass175.5") leg->AddEntry(histo,"Mass = 175.5 GeV",  "l");
    }
    if(theoryName == "matchdown" || theoryName == "mass169.5"){
        histo->SetLineStyle(7);
        histo->SetLineColor(kRed-7);
        if(leg && theoryName == "matchdown") leg->AddEntry(histo,"Matching down",  "l");
        if(leg && theoryName == "mass169.5")   leg->AddEntry(histo,"Mass = 169.5 GeV",  "l");
    }
    if(theoryName == "scaleup" || theoryName == "mass173.5"){
        histo->SetLineStyle(2);
        histo->SetLineColor(kAzure+2);
        if(leg && theoryName == "scaleup") leg->AddEntry(histo,"4*Q^{2}",  "l");
        if(leg && theoryName == "mass173.5") leg->AddEntry(histo,"Mass = 173.5 GeV",  "l");
    }
    if(theoryName == "scaledown" || theoryName == "mass171.5"){
        histo->SetLineStyle(2);
        histo->SetLineColor(8);
        if(leg && theoryName == "scaledown") leg->AddEntry(histo,"Q^{2}/4",  "l");
        if(leg && theoryName == "mass171.5")   leg->AddEntry(histo,"Mass = 171.5 GeV",  "l");
    }
    if(theoryName == "mass178.5"){
        histo->SetLineStyle(3);
        histo->SetLineColor(42);
        if(leg && theoryName == "matchup") leg->AddEntry(histo,"Matching up",  "l");
        if(leg && theoryName == "mass178.5") leg->AddEntry(histo,"Mass = 178.5 GeV",  "l");
    }
    if(theoryName == "mass166.5"){
        histo->SetLineStyle(3);
        histo->SetLineColor(48);
        if(leg && theoryName == "matchdown") leg->AddEntry(histo,"Matching down",  "l");
        if(leg && theoryName == "mass166.5")   leg->AddEntry(histo,"Mass = 166.5 GeV",  "l");
    }

}

bool Plotter::addQCDToControlPlot()const
{
    if( name.Contains("_step2") || name.Contains("_step3") || name.Contains("_step4") || name.Contains("_step5") || name.Contains("_step6") || name.Contains("_step7") ||
        name.Contains("_noBTag") || name.Contains("_diLep") || name.Contains("Electron") || name.Contains("Muon") || name == "MET" ||
        (!name.Contains("Hyp") && (name.Contains("jetpT") || name.Contains("jetHT")) ))
       {
            return 1;
        }
    return 0;
}


void Plotter::getSignalUncertaintyBand(TH1* uncBand, TString channel_)
{
    if(!uncBand)  return;
    std::vector<TString> syst {"MASS_", "SCALE_", "MATCH_", "HAD_", "PDF_"};

    double norm_Events = uncBand->Integral(-1e6, 1e6);
    int nbins = uncBand->GetNbinsX();
    std::vector<double> vec_varup(nbins, 0), vec_vardown(nbins, 0);
    for (size_t iter = 0; iter<syst.size(); iter++)
    {
        TString file_up = "", file_do = "";
        if(syst.at(iter) == "HAD_"){
            file_up = TString("Plots/MCATNLO/"+channel_+"/"+name+"_source.root");
            file_do = TString("Plots/POWHEG/"+channel_+"/"+name+"_source.root");
        } else{
            file_up = TString("Plots/"+syst.at(iter)+"UP/"+channel_+"/"+name+"_source.root");
            file_do = TString("Plots/"+syst.at(iter)+"DOWN/"+channel_+"/"+name+"_source.root");
        }
        ifstream inputFileStream(file_up);
        if(inputFileStream.fail()){continue;}
        inputFileStream.close();
        ifstream inputFileStream2(file_do);
        if(inputFileStream2.fail()){continue;}
        inputFileStream2.close();

        // This lines crashes the code, some probles arises form the HistoListReader class
        TH1D *tmpUp = fileReader->GetClone<TH1D>(file_up, name+"_allttbar", 1);
        TH1D *tmpDo = fileReader->GetClone<TH1D>(file_do, name+"_allttbar", 1);

        if(!tmpUp && tmpDo){ delete tmpDo; continue;}
        if(tmpUp && !tmpDo){ delete tmpUp; continue;}
        if(!tmpUp && !tmpDo) continue;

        if (nbins != tmpUp->GetNbinsX() || nbins != tmpDo->GetNbinsX()){
            continue;
        }
        tmpDo->SetDirectory(0);
        tmpUp->SetDirectory(0);

        double upIntegral = tmpUp->Integral(-1e6, 1e6);
        double doIntegral = tmpDo->Integral(-1e6, 1e6);

        tmpUp->Scale(norm_Events/upIntegral);
        tmpDo->Scale(norm_Events/doIntegral);
        for (Int_t nbin = 0; nbin < nbins; nbin++)
        {
            double binContent = uncBand->GetBinContent(nbin+1);
            double rel_diffup = std::abs(tmpUp->GetBinContent(nbin+1) - binContent) / binContent;
            double rel_diffdo = std::abs(tmpDo->GetBinContent(nbin+1) - binContent) / binContent;
            if (binContent<1e-6){
                rel_diffdo = 0;
                rel_diffup = 0;
            }
            vec_varup.at(nbin) += rel_diffup * rel_diffup;
            vec_vardown.at(nbin) += rel_diffdo * rel_diffdo;
        }
        delete tmpDo; delete tmpUp;
    }
    for (size_t iter = 0; iter< vec_varup.size(); iter++)
    {
        double centralValue =  uncBand->GetBinContent(iter+1);
        uncBand->SetBinError(iter+1, centralValue * 0.5 * (std::sqrt(vec_vardown.at(iter)) + std::sqrt(vec_varup.at(iter))));
    }
}



void Plotter::yRangeControlPlotRatio(double &yminCP_, double &ymaxCP_)const
{
    yminCP_ = 0.51;
    ymaxCP_ = 1.49;

    if(name.Contains("HypToppT"))               { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypTopRapidity"))         { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypTTBarDeltaRapidity"))  { yminCP_ = 0.51; ymaxCP_ = 1.49;}

    if(name.Contains("HypTTBarpT"))             { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypTTBarMass"))           { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypTTBarRapidity"))       { yminCP_ = 0.51; ymaxCP_ = 1.49;}

    if(name.Contains("HypLeptonpT"))            { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypLeptonEta"))           { yminCP_ = 0.51; ymaxCP_ = 1.49;}

    if(name.Contains("HypLLBarMass"))           { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypLLBarpT"))             { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypLLBarDPhi"))           { yminCP_ = 0.51; ymaxCP_ = 1.49;}

    if(name.Contains("HypBJetpT"))              { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypBJetEta"))             { yminCP_ = 0.51; ymaxCP_ = 1.49;}

    if(name.Contains("HypBBBarMass"))           { yminCP_ = 0.51; ymaxCP_ = 1.49;}
    if(name.Contains("HypBBBarpT"))             { yminCP_ = 0.51; ymaxCP_ = 1.49;}

    if(name.Contains("HypLeptonBjetMass"))      { yminCP_ = 0.51; ymaxCP_ = 1.49;}
}



void Plotter::setResultRatioRanges(double &yminRes_, double &ymaxRes_ )const
{
    yminRes_ = 0.49;
    ymaxRes_ = 1.51;

    if(name.Contains("HypToppTTTRestFrame"))    { yminRes_ = 0.70; ymaxRes_ = 1.65; return;}
    if(name.Contains("HypToppTNLead"))          { yminRes_ = 0.70; ymaxRes_ = 1.65; return;}
    if(name.Contains("HypToppTLead"))           { yminRes_ = 0.70; ymaxRes_ = 1.65; return;}
    if(name.Contains("HypToppT"))               { yminRes_ = 0.75; ymaxRes_ = 1.65; return;}

    if(name.Contains("HypTopRapidity"))         { yminRes_ = 0.85; ymaxRes_ = 1.25; return;}

    if(name.Contains("HypTTBarDeltaRapidity"))  { yminRes_ = 0.51; ymaxRes_ = 1.49; return;}
    if(name.Contains("HypTTBarDeltaPhi"))       { yminRes_ = 0.85; ymaxRes_ = 1.25; return;}

    if(name.Contains("HypTTBarpT"))             { yminRes_ = 0.50; ymaxRes_ = 1.79; return;}
    if(name.Contains("HypTTBarMass"))           { yminRes_ = 0.70; ymaxRes_ = 1.95; return;}
    if(name.Contains("HypTTBarRapidity"))       { yminRes_ = 0.90; ymaxRes_ = 1.29; return;}

    if(name.Contains("HypLeptonpT"))            { yminRes_ = 0.85; ymaxRes_ = 1.35; return;}
    if(name.Contains("HypLeptonEta"))           { yminRes_ = 0.85; ymaxRes_ = 1.35; return;}

    if(name.Contains("HypLLBarMass"))           { yminRes_ = 0.85; ymaxRes_ = 1.25; return;}
    if(name.Contains("HypLLBarpT"))             { yminRes_ = 0.85; ymaxRes_ = 1.25; return;}
    if(name.Contains("HypLLBarDPhi"))           { yminRes_ = 0.51; ymaxRes_ = 1.49; return;}

    if(name.Contains("HypBJetpT"))              { yminRes_ = 0.70; ymaxRes_ = 1.75; return;}
    if(name.Contains("HypBJetEta"))             { yminRes_ = 0.80; ymaxRes_ = 1.25; return;}

    if(name.Contains("HypBBBarMass"))           { yminRes_ = 0.90; ymaxRes_ = 1.14; return;}
    if(name.Contains("HypBBBarpT"))             { yminRes_ = 0.60; ymaxRes_ = 1.49; return;}

    if(name.Contains("HypLeptonBjetMass"))      { yminRes_ = 0.70; ymaxRes_ = 1.45; return;}
}
