#include <map>
#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveText.h>


#include "UsefulTools.h"
#include "../../common/include/RootFileReader.h"

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

UsefulTools::UsefulTools( RootFileReader* rootFileReader,bool isClosureTest,bool isDYScale):

lumi(19712),

//topxsec(244.849), //again changes with normalization, must be set outside of the class
//topxsec(244.794), //Mitov, arXiv:1303.6254
//topxsec(247.998), //Measured XSection after normalization to Mitov, and statistical combination of channels
//topxsec(248.207), //Measured XSection after normalization to Mitov, and after combination of channels at yield level
topxsec(245.102), //Measured XSection using MadSpin, after normalization to Mitov and statistical combination of channels

fileReader(rootFileReader),
doClosureTest(isClosureTest),
doDYScale(isDYScale)
{
    for (auto s: VectorOfValidSystematics) ListOfSyst.insert(s);
}

void UsefulTools::fillSetListOfSystematics(std::set<TString>& set)
{

    for (auto s: VectorOfValidSystematics) set.insert(s);
}

void UsefulTools::fillVectorOfValidSystematics(std::vector<const char*>& vect)
{

    for (auto s: VectorOfValidSystematics) vect.push_back(s);
}

double UsefulTools::SampleXSection(const TString& filename){
    
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


void UsefulTools::fillLegendColorDataset(const TString& fileListName, std::vector<TString>& legends, std::vector<int>& colors, std::vector<TString>& dataset){

        std::cout << "reading " << fileListName << std::endl;
    
        ifstream FileList(fileListName);
        if (FileList.fail()){
            std::cerr << "Error reading " << fileListName << std::endl;
            exit(1);
        }
        
        TString filename;
        
        dataset.clear();
        legends.clear();
        colors.clear();
    
        while(!FileList.eof()){ 
        FileList>>filename;
        if(filename==""){continue;}//Skip empty lines
        dataset.push_back(filename);
        if(filename.Contains("run")){legends.push_back("Data"); colors.push_back(kBlack);}
        else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} Signal"); colors.push_back(kRed+1);}
        else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} Other"); colors.push_back(kRed-7);}
        else if(filename.Contains("single")){legends.push_back("Single Top"); colors.push_back(kMagenta);}
        else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); colors.push_back(10);}
        else if(filename.Contains("dytautau")){legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
        else if(filename.Contains("dymumu")||filename.Contains("dyee")){legends.push_back("Z / #gamma* #rightarrow ee/#mu#mu"); colors.push_back(kAzure-2);}
        else if(filename.Contains("wtolnu")){legends.push_back("W+Jets"); colors.push_back(kGreen-3);}
        else if(filename.Contains("qcd")){legends.push_back("QCD Multijet"); colors.push_back(kYellow);}
        else if(filename.Contains("ttbarZ") ||filename.Contains("ttbarW") || filename.Contains("ttgjets")){legends.push_back("t#bar{t}+Z/W/#gamma"); colors.push_back(kOrange-2);}
    }
    FileList.close();
}


double UsefulTools::CalcLumiWeight(const TString& WhichSample){
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

std::vector<TString> UsefulTools::InputFileList(TString mode, TString Systematic)
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
        for(std::set<TString>::iterator iter = ListOfSyst.begin(); iter!=ListOfSyst.end(); ++iter){ std::cout<<*iter<<std::endl;}
        return FileVector;
    }
    
    if( mode.CompareTo("combined") && mode.CompareTo("ee") && mode.CompareTo("emu") && mode.CompareTo("mumu")){
        std::cout<<"The decay channel you provided is not supported."<<std::endl;
        std::cout<<"Please use: ee, emu, mumu, combined"<<std::endl;
        return FileVector;
    }
    
    if(!mode.CompareTo("combined")){
        std::vector<TString> eemode   = UsefulTools::InputFileList(TString("ee"), Systematic);
        std::vector<TString> emumode  = UsefulTools::InputFileList(TString("emu"), Systematic);
        std::vector<TString> mumumode = UsefulTools::InputFileList(TString("mumu"), Systematic);
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
    FileVector.push_back(tempName + "_ttbarW.root");
    FileVector.push_back(tempName + "_ttbarZ.root");
    FileVector.push_back(tempName + "_ttgjets.root");
    
    
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

void UsefulTools::ApplyFlatWeights(TH1* varhists, const double weight){

    if(weight == 0) {std::cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<std::endl;}
    //if(weight >=1e3){std::cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<std::endl;}
    varhists->Scale(weight);
}

void UsefulTools::ApplyFlatWeights(TH2* varhists, const double weight){

    if(weight == 0) {std::cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<std::endl;}
    //if(weight >=1e3){std::cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<std::endl;}
    varhists->Scale(weight);
}

void UsefulTools::DYScaleFactor(TString SpecialComment,std::vector<double>& DYScale,TString name){

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

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void UsefulTools::DrawCMSLabels(int cmsprelim, double energy, double textSize) {

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
    //label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
    //label->SetY2NDC(1.0);
    label->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.02 );
    label->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() + 0.03 );
    
    label->SetTextFont(42);
    label->AddText(Form(text, lumi/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}


void UsefulTools::setStyle(TH1 *hist, TString Axis)
{
//     hist->SetLineColor(2);
//     hist->SetLineWidth(1);
//     hist->GetXaxis()->SetLabelFont(42);
//     hist->GetYaxis()->SetLabelFont(42);
//     hist->GetXaxis()->SetTitleFont(42);
//     hist->GetYaxis()->SetTitleFont(42);
//     hist->GetXaxis()->SetTitleSize(0.05);
//     hist->GetYaxis()->SetTitleSize(0.05);
//     hist->GetXaxis()->SetTitleOffset(1.08);
//     hist->GetYaxis()->SetTitleOffset(1.7);
//     hist->GetXaxis()->SetLabelOffset(0.007);
//     hist->GetYaxis()->SetLabelOffset(0.007);
// 
// 
//         hist->SetFillColor(0);
//     hist->SetMarkerStyle(20);
        if ((Axis.Contains("p_{T}", TString::kIgnoreCase) || Axis.Contains("m(t#bar{t})", TString::kIgnoreCase)) && 
            (!name.Contains("1st") && !name.Contains("Rapidity") && !name.Contains("Eta") && !name.Contains("Phi") && !name.Contains("JetMult"))) {
            hist->GetXaxis()->SetTitle(Axis+" #left[GeV#right]");
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+Axis+"}"+" #left[GeV^{-1}#right]"); 
            
        } else {
            hist->GetXaxis()->SetTitle(Axis);
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+Axis+"}");
        }

}


