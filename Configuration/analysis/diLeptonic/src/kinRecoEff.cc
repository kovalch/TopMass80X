#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TSystem.h>

#include "EffHist.h"
#include "../../common/include/RootFileReader.h"
#include "homelessFunctions.h"

using namespace std;

double sqr(double value) { return value*value; }




int main(int argc, char *argv[])
{
  
              gSystem->Exec("mkdir Plots");
    
  RootFileReader *fileReader_ = RootFileReader::getInstance();
  
  vector<TString> systematics {"Nominal"};
  vector<TString> channels {"ee","emu","mumu"};

  TString rootFilePrePath("/data/user/korol/CMSSW_5_3_11/src/TopAnalysis/Configuration/analysis/diLeptonic/");
  
  vector<TString> effHistNames {"KinReco_nRecoEvt","KinReco_nRecoEvt","_Eff","_vs_JetMult","_vs_LepEta","_vs_JetEta","_vs_LeppT","_vs_MET","_vs_LepEta2","_vs_JetEta2","_vs_LeppT2","_vs_AntiLeptonpT","_vs_LeptonpT","_vs_LeptonEta","_vs_AntiLeptonEta","_vs_JetpT","_vs_JetpT2"};
  
  double NdataAfter = 0;
  double NdataBefore = 0;
  double dataEff = 0;
  double dataEffUnc = 0; 
  
  double NmcAfter = 0;
  double NmcBefore = 0;
  double allmcEff =  0;
  double allmcEffUnc = 0;
  
  double sfUnc = 0;
  
  char dataEffString[100];
  char allmcEffString[100];
  char sfString[100];
  
  //homelessFunctions* homelessFunc = new homelessFunctions(fileReader_,doClosureTest,doDYScale);
  homelessFunctions* homelessFunc = new homelessFunctions(fileReader_,0,1);
  
  for ( auto systematic : systematics ){
              gSystem->Exec("mkdir Plots/"+systematic);
    for ( TString channel : channels ){
        gSystem->Exec("mkdir Plots/"+systematic+"/"+channel);
        
        std::vector<TString> dataset;
        std::vector<TString> legends;
        int channelType;
        TString channelLabel;
        if(channel.Contains("ee")){channelLabel="ee";channelType=0;}
        if(channel.Contains("mumu")){channelLabel="#mu#mu";channelType=1;}
        if(channel.Contains("emu")){channelLabel="e#mu";channelType=2;}
        if(channel.Contains("combined")){channelLabel="Dilepton Combined";channelType=3;}
        
        std::vector<double> DYScale={1,1,1,1};
        homelessFunc->DYScaleFactor("Standard",DYScale,"");
        
        
        for(Size_t i=0;i<DYScale.size();++i) std::cout << "DYScale " << i << " " << DYScale.at(i) << std::endl;
        
        TString fileListName("/data/user/korol/CMSSW_5_3_11/src/TopAnalysis/Configuration/analysis/diLeptonic/FileLists/HistoFileList_");
        fileListName.Append(systematic+"_"+channel+".txt");
        ifstream FileList(fileListName);
        if (FileList.fail()) std::cerr << "Error reading " << fileListName << std::endl;
        
        for ( int j=2; j<(int)effHistNames.size();j++ ){
        
            TString filename;
            while(!FileList.eof()){
                FileList>>filename;
                if(filename==""){continue;}//Skip empty lines
                //filename.Prepend("/"+systematic+"/"+channel+"/");
                filename.Prepend(rootFilePrePath);
                dataset.push_back(filename);
                if(filename.Contains("run")){legends.push_back("Data"); }
                else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} Signal"); }
                else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} Other"); }
                else if(filename.Contains("single")){legends.push_back("Single Top"); }
                else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); }
                else if(filename.Contains("dytautau")){legends.push_back("Z / #gamma* #rightarrow #tau#tau"); }
                else if(filename.Contains("dymumu")||filename.Contains("dyee")){legends.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); }
                else if(filename.Contains("wtolnu")){legends.push_back("W+Jets"); }
                else if(filename.Contains("qcd")){legends.push_back("QCD Multijet"); }
                else if(filename.Contains("ttbarZ") ||filename.Contains("ttbarW") || filename.Contains("ttgjets")){legends.push_back("t#bar{t}+Z/W/#gamma");}
            }
            TString DYEntry = "Z / #gamma* #rightarrow ee/#mu#mu";

            TString h_NumName=effHistNames.at(0);
            TString h_DeNumName=effHistNames.at(1);
            h_NumName+=effHistNames.at(j);
            h_DeNumName+=effHistNames.at(j);
            h_NumName.Append("_step8");
            h_DeNumName.Append("_step7");
            std::vector<TH1D> histsNum;
            std::vector<TH1D> histsDeNum;
        
        for(unsigned int i=0; i<dataset.size(); i++){
            TH1D * histNum = fileReader_->GetClone<TH1D>(dataset.at(i),h_NumName);
            TH1D * histDeNum = fileReader_->GetClone<TH1D>(dataset.at(i),h_DeNumName);
            
            if (!histNum) std::cout<< "Error: Can not read histo: "<< h_NumName << " from file: " << dataset.at(i) << "\n";
            if (!histDeNum) std::cout<< "Error: Can not read histo: "<< h_DeNumName << " from file: " << dataset.at(i) << "\n";
            
            //Rescaling to the data luminosity
            double LumiWeight = homelessFunc->CalcLumiWeight(dataset.at(i));
            homelessFunc->ApplyFlatWeights(histNum, LumiWeight);
            homelessFunc->ApplyFlatWeights(histDeNum, LumiWeight);
            
            if((legends.at(i) == DYEntry)&&(channelType!=2)){
                histNum->Scale(DYScale[channelType]);
                histDeNum->Scale(DYScale[channelType]);
            }
            
            histsNum.push_back(*histNum);
            histsDeNum.push_back(*histDeNum);
        }

        if (histsNum.size() == 0)std::cerr << "***ERROR! No num-histograms available! "<< std::endl; 
        if (histsDeNum.size() == 0)std::cerr << "***ERROR! No denum-histograms available! "<< std::endl; 
    
        TH1 *mcNumHist = NULL;
        TH1 *mcDeNumHist = NULL;
        TH1 *dataNumHist = NULL;
        TH1 *dataDeNumHist = NULL;
        
        if(legends.at(0) != "Data") std::cout << "***ERROR! legends_.at(0) != Data " << std::endl;
        for(Size_t i=0; i<histsNum.size(); i++)
        {
            if(legends.at(i) == "Data"){
                if(dataNumHist==NULL){
                    dataNumHist=&histsNum.at(i);
                    dataDeNumHist=&histsDeNum.at(i);
                }
                else{ 
                    dataNumHist->Add(&histsNum.at(i));
                    dataDeNumHist->Add(&histsDeNum.at(i));
                }
            }
              else{
//           else if(legends.at(i) != "t#bar{t} Signal"){
                if(mcNumHist==NULL){
                    mcNumHist=&histsNum.at(i);
                    mcDeNumHist=&histsDeNum.at(i);
                }
                else{ 
                    mcNumHist->Add(&histsNum.at(i));
                    mcDeNumHist->Add(&histsDeNum.at(i));
                }
            }
        }
        
          EffHist effHist(effHistNames.at(j));

          TString savepath;
          savepath.Append(+"./Plots/"+systematic+"/"+channel+"/"+"KinRecoEff_"+effHistNames.at(j)+".root");
          savepath.ReplaceAll("_vs_",4,"",0);

            //TString title(channel);
            TString title("");
            if(effHistNames.at(j)=="_vs_JetMult")title.Append("  JetMult - Kin. Reco. behaviour; jet multiplicity; eff.");
            if(effHistNames.at(j)=="_vs_LepEta")title.Append("  lead-LeptonEta - Kin. Reco. behaviour; #eta(lead-Lepton); eff.");
            if(effHistNames.at(j)=="_vs_LepEta2")title.Append("  nLead-LeptonEta - Kin. Reco. behaviour; #eta(nLead-Lepton); eff.");
            if(effHistNames.at(j)=="_vs_JetEta")title.Append("  lead-BjetEta - Kin. Reco. behaviour; #eta(Ljet); eff.");
            if(effHistNames.at(j)=="_vs_JetEta2")title.Append("  nLead-BjetEta - Kin. Reco. behaviour; #eta(nLjet); eff.");
            if(effHistNames.at(j)=="_vs_JetpT")title.Append("  lead-BjetpT - Kin. Reco. behaviour; p_{T}(Ljet); eff.");
            if(effHistNames.at(j)=="_vs_JetpT2")title.Append("  nLead-BjetpT - Kin. Reco. behaviour; p_{T}(nLjet); eff.");
            if(effHistNames.at(j)=="_vs_LeppT")title.Append(" lead-LeptonpT - Kin. Reco. behaviour; p_{T}(lead-Lepton) [GeV]; eff.");
            if(effHistNames.at(j)=="_vs_LeppT2")title.Append("  nLead-LeptonpT - Kin. Reco. behaviour; p_{T}(nLead-Lepton) [GeV]; eff.");
            if(effHistNames.at(j)=="_vs_AntiLeptonpT")title.Append("  Anti-LeptonpT - Kin. Reco. behaviour; p_{T}(antiLepton) [GeV]; eff.");
            if(effHistNames.at(j)=="_vs_LeptonpT")title.Append(" LeptonpT - Kin. Reco. behaviour; p_{T}(Lepton) [GeV]; eff.");
            if(effHistNames.at(j)=="_vs_LeptonEta")title.Append("  LeptonEta - Kin. Reco. behaviour; #eta(Lepton); eff.");
            if(effHistNames.at(j)=="_vs_AntiLeptonEta")title.Append("  Anti-LeptonEta - Kin. Reco. behaviour; #eta(antiLepton); eff.");
            if(effHistNames.at(j)=="_vs_MET")title.Append("  MET - Kin. Reco. behaviour; MET; eff.");
            effHist.setTitle(title);
            std::cout << title << std::endl;
            
            
            if(j==2){
                NdataAfter=dataNumHist->GetBinContent(1);
                NdataBefore=dataDeNumHist->GetBinContent(1);
                NmcAfter=mcNumHist->GetBinContent(1);
                NmcBefore=mcDeNumHist->GetBinContent(1);
            }
            
            if(j>2){
                dataEff = NdataAfter / NdataBefore;
                dataEffUnc = sqrt(dataEff * (1-dataEff) / NdataBefore);
                
                allmcEff =  NmcAfter / NmcBefore;
                allmcEffUnc = sqrt(allmcEff * (1-allmcEff) / NmcBefore);
                
                sfUnc = sqrt(sqr(dataEffUnc/dataEff) + sqr(allmcEffUnc/allmcEff));
                
                sprintf(dataEffString, "%.2f%%", 100*dataEff);
                sprintf(allmcEffString, "%.2f%%", 100*allmcEff);
                sprintf(sfString, "%.4f +- %.4f", dataEff/allmcEff, sfUnc);
                effHist.addData((TH1D*)dataNumHist,(TH1D*)dataDeNumHist);
                effHist.addMc((TH1D*)mcNumHist,(TH1D*)mcDeNumHist);
                
                effHist.savePlotEffSF(savepath,dataEffString,allmcEffString,sfString);
            }

      }
        
    }
  }


  return 0;
}
