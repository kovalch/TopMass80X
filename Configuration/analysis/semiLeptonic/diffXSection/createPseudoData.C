// ---------description----------------------------------------
// ------------------------------------------------------------
// this Makro creates pseudo data by taking the same plot from 
// various MC samples, weighting them concerning luminosity, 
// adding them and applies a poisson smearing for every point
// ------------------------------------------------------------
// adaptions are necessary within areas marked via !!! ... !!! 
// ------------------------------------------------------------

#include "basicFunctions.h"
#include <TFile.h>
#include <TRandom3.h>
#include <TKey.h>
#include <TDirectory.h>

void poisson(const std::map< TString, std::map <unsigned int, TH1F*> > histo_, const std::vector<TString> plotList_, const std::string decayChannel, TFile& outputfile, const int luminosity, const unsigned int verbose=0, bool smear=1, bool useReweightedTop=0, double avReweight=1, bool useZprime=0, double zPrimeLumiWeight=1);
void createPseudoDataFunc(double luminosity=19712., const std::string decayChannel="muon", TString specifier="NoDistort");

void createPseudoData(double luminosity=19712.){

  // collect leptons
  std::vector<std::string> leptons_;
  leptons_.push_back("muon"    );
  leptons_.push_back("electron");
  // collect relevant tests
  std::vector<TString> tests_;
  //tests_.push_back("NoDistort"    );
  //tests_.push_back("topPtUp"      );
  //tests_.push_back("topPtDown"    );
  //tests_.push_back("ttbarMassUp"  );
  //tests_.push_back("ttbarMassDown");
  //tests_.push_back("data"         );
  //tests_.push_back("1000"         );
  tests_.push_back("topMass161p5"   );
  tests_.push_back("topMass163p5"   );
  tests_.push_back("topMass166p5"   );
  tests_.push_back("topMass169p5"   );
  tests_.push_back("topMass175p5"   );
  tests_.push_back("topMass178p5"   );
  tests_.push_back("topMass181p5"   );
  tests_.push_back("topMass184p5"   );

  // create all relevant files
  for(unsigned int lepton=0; lepton<leptons_.size(); ++lepton){
    for(unsigned int test=0; test<tests_.size(); ++test){
      std::cout << "- " << leptons_[lepton] << ", " << tests_[test] << std::endl;
      createPseudoDataFunc(luminosity, leptons_[lepton], tests_[test]);
    }
  }

}

void createPseudoDataFunc(double luminosity, const std::string decayChannel, TString specifier){
  // specifier="NoDistort" for closure test; "topPtUp", "topPtDown", "ttbarMassUp", "ttbarMassDown" or "data" for shape distortions; "1000" for Zprime
  // "verbose": set detail level of output ( 0: no output, 1: std output 2: output for debugging )
  int verbose=0;
  // "smear": say if you want to do a poisson smearing for each bin or just a combination for the different samples 
  bool smear=false;
  // "dataFile": absolute path of data file, used to define plots of interest
  TString dataFile= "";
  if(decayChannel.compare("electron")==0) dataFile="/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/elecDiffXSecData2012ABCDAll.root";
  else                                    dataFile="/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/muonDiffXSecData2012ABCDAll.root";
  // "useReweightedTop": use parton level reweighted ttbar signal file in pseudo data?
  bool useReweightedTop = (specifier.Contains("NoDistort")) ? false : true;
  // "zprime": include additional Zprime in pseudo data?
  bool zprime=specifier.Contains("1000") ? true : false;
  // naming for outputfile, constructed within macro
  TString outNameExtension="";

  // check input parameter validity
  // type of closure test
  if(!(specifier.Contains("NoDistort")||specifier.Contains("1000")||specifier.Contains("data")||specifier.Contains("topPtUp")||specifier.Contains("topPtDown")||specifier.Contains("ttbarMassUp")||specifier.Contains("ttbarMassDown")||specifier.Contains("topMass161p5")||specifier.Contains("topMass163p5")||specifier.Contains("topMass166p5")||specifier.Contains("topMass169p5")||specifier.Contains("topMass175p5")||specifier.Contains("topMass178p5")||specifier.Contains("topMass181p5")||specifier.Contains("topMass184p5"))){
    std::cout << "ERROR: invalid input specifier=" << specifier << std::endl;
    std::cout << "supported are: specifier=NoDistort,1000,data,topPtUp,topPtDown,ttbarMassUp,ttbarMassDown,topMass161p5,topMass163p5,topMass166p5,topMass169p5,topMass175p5,topMass178p5,topMass181p5,topMass184p5" << std::endl;
    exit(0);
  }
  // decay channel
  if(!(decayChannel.compare("electron")==0||decayChannel.compare("muon")==0)){
    std::cout << "ERROR: invalid input decayChannel=" << decayChannel << std::endl;
    std::cout << "supported are: decayChannel=muon,electron" << std::endl;
    exit(0);
  }
  

  //  ---
  //     create container for histos and files
  //  ---
  std::map<unsigned int, TFile*> files_;
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  std::map< TString, std::map <unsigned int, TH2F*> > histo2_;

  //  ---
  //     parton level reweighted top distribution
  //  ---
  // a) name and path of rootfile
  // path
  TString nameTtbarReweighted="/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/";
  // subfolder for reweighting systematics
  if(specifier!="NoDistort"&&!specifier.Contains("p5")){
    nameTtbarReweighted+="ttbarReweight/";
    // SG reweighting test
    if(decayChannel.compare("electron")==0) nameTtbarReweighted+="elecDiffXSec";
    else                                    nameTtbarReweighted+="muonDiffXSec";
    nameTtbarReweighted+="SigSysDistort"+specifier+"Summer12PF.root";
  }
  else{
    // SG sample for corresponding top mass
    int sys=sysNo; // == "NoDistort"
    if(     specifier=="topMass161p5") sys=sysTopMassDown4;
    else if(specifier=="topMass163p5") sys=sysTopMassDown3;
    else if(specifier=="topMass166p5") sys=sysTopMassDown2;
    else if(specifier=="topMass169p5") sys=sysTopMassDown;
    else if(specifier=="topMass175p5") sys=sysTopMassUp;
    else if(specifier=="topMass178p5") sys=sysTopMassUp2;
    else if(specifier=="topMass181p5") sys=sysTopMassUp3;
    else if(specifier=="topMass184p5") sys=sysTopMassUp4;
    nameTtbarReweighted+=TopFilename(kSig, sys, decayChannel);
  }
  // BG
  TString nameTtbarBGReweighted=nameTtbarReweighted;
  nameTtbarBGReweighted.ReplaceAll("Sig","Bkg");
  if(useReweightedTop&&!specifier.Contains("p5")) outNameExtension+="Reweighted"+specifier;
  // b) get average weight for reweighted samples
  double avWeight=1;
  if(useReweightedTop && specifier!="NoDistort"&& !specifier.Contains("p5")){
    TFile* ttbarRewfile = new (TFile)(nameTtbarReweighted);
    TString weightPlot="eventWeightDileptonModelVariation/modelWeightSum";
    histo_["avWeight"][kSig] = (TH1F*)(ttbarRewfile->Get(weightPlot)->Clone());
    avWeight=histo_ ["avWeight"][kSig]->GetBinContent(2)/histo_ ["avWeight"][kSig]->GetBinContent(1);
    histo_["avWeight"].erase(kSig);
    TFile* sigfile = new (TFile)("/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/"+TopFilename(kSig, sysNo, decayChannel));
    TString partonPlot="analyzeTopPartonLevelKinematics/ttbarMass";
    gROOT->cd();
    histo_["parton"   ][kSig] = (TH1F*)(sigfile     ->Get(partonPlot)->Clone());
    histo_["partonRew"][kSig] = (TH1F*)(ttbarRewfile->Get(partonPlot)->Clone());
    double avWeight2  = histo_["partonRew"][kSig]->Integral(0, histo_["partonRew"][kSig]->GetNbinsX()+1);
    avWeight2        /= histo_["parton"   ][kSig]->Integral(0, histo_["parton"   ][kSig]->GetNbinsX()+1);
    if(verbose>1){
      std::cout << "ratio unweighted/weighted" << std::endl;
      std::cout << avWeight  << " (from " << weightPlot << ")" << std::endl;
      std::cout << avWeight2 << TString(" (from entries in ")+partonPlot+" wrt /afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/"+TopFilename(kSig, sysNo, decayChannel)+"): " << std::endl;
    }
  }

  //  ---
  //     Z prime 
  //  ---
  // xSec scaling (default value chosen for M1000W100: 5pb)
  double xSecSF=1.0;
  // path & naming
  TString zprimeMass =specifier;
  TString zprimeWidth=specifier;
  zprimeWidth.ReplaceAll("1000", "100");
  TString nameZprime="/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/zprime/";
  if(decayChannel.compare("electron")==0) nameZprime+="elec";
  else                                    nameZprime+="muon";
  nameZprime+="DiffXSecZprimeM"+zprimeMass+"W"+zprimeWidth+"Summer12PF.root";
  // event weight wrt luminosity
  double zPrimeLumiWeight=1.;
  if     (zprimeMass=="1000") zPrimeLumiWeight=(xSecSF*5*luminosity)/104043;
  // naming of oututfile 
  TString xSecSFstr="";
  if      (xSecSF==0.5 ) xSecSFstr="x0p5";
  else if (xSecSF==0.25) xSecSFstr="x0p25";
  else if (xSecSF==0.1 ) xSecSFstr="x0p1";
  else if (xSecSF==0.03) xSecSFstr="x0p03";
  else if (xSecSF>1.   ) xSecSFstr="x"+getTStringFromDouble(xSecSF,0);
  if(zprime) outNameExtension="andM"+zprimeMass+"W"+zprimeWidth+xSecSFstr+"Zprime";
  // sample iteration limit
  int kLast  =kSAToptW;
  int kZprime=kLast+1;
  if(zprime) kLast=kZprime;

  if(zprime&&zPrimeLumiWeight==1){
    std::cout << "ERROR: chosen zprime weight is exactly 1!" << std::endl;
    exit(0);
  }
  
  // -------------------------
  // !!! choose luminosity !!!
  // -------------------------
  TString lum = getTStringFromInt(roundToInt(luminosity));
  // -----------------------------------------
  // !!! add all contributing samples here !!!
  // -----------------------------------------
  TString inputFolder = "/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/";
  // save all default top analysis samples in files_
  files_ = getStdTopAnalysisFiles(inputFolder, sysNo, dataFile, decayChannel);
  // remove combined file dor single Top & DiBoson
  if(files_.count(kSTop )>0)files_.erase(kSTop );
  if(files_.count(kDiBos)>0)files_.erase(kDiBos);
  // remove combined QCD file for electron channel
  if(decayChannel.compare("electron")==0&&files_.count(kQCD)>0) files_.erase(kQCD);
  // remove all QCD files to be consistent with later treatment
  if(files_.count(kQCD)>0) files_.erase(kQCD);
  if(decayChannel.compare("electron")==0){
    for(int sample = kQCDEM1; sample<=kQCDBCE3; sample++){
      if(files_.count(sample)>0) files_.erase(sample);
    }
  }
  // add zprime
  if(zprime) files_[kZprime]=new (TFile)(nameZprime);
  // exchange default ttbar files with the reweighted ones
  if(useReweightedTop){ 
    files_[kSig]=new (TFile)(nameTtbarReweighted  );
    files_[kBkg]=new (TFile)(nameTtbarBGReweighted);
  }

  //  -----------------------------------------
  //     get list of all plots to be considered
  //  -----------------------------------------
  std::vector<TString> plotList_;
  TString currentFolder ="";
  TString currentPlot   ="";
  if(verbose>0) std::cout << std::endl << "searching for all plots in file " << files_[kData]->GetName() << std::endl;
  // go to data file for getting plot names
  files_[kData]->cd();
  // loop objects in file
  TIter fileIterator(gDirectory->GetListOfKeys());
  if(verbose>2){
    std::cout << "looping objects in in ";
    gDirectory->pwd();
  }
  TKey *fileKey;
  int count=0;
  while( (fileKey = (TKey*)fileIterator()) ){
    ++count;
    if(verbose>2) std::cout << "folderObject #" << count;
    TObject *fileObject = fileKey->ReadObj(); 
    // check if object is a directory
    if(!(fileObject->InheritsFrom("TDirectory"))&&(verbose>2)) std::cout << " is no directory" << count << std::endl;
    if(fileObject->InheritsFrom("TDirectory")){
      currentFolder = (TString)fileObject->GetName();
      if(verbose>2) std::cout << " ("+currentFolder+")" << std::endl;
      // go to directory
      ((TDirectory*)fileObject)->cd();
      if(verbose>2){
	std::cout << "looping objects in in ";
	gDirectory->pwd();
      }
      TIter folderIterator(gDirectory->GetListOfKeys());
      TKey *folderKey;
      int count2=0;
      while( (folderKey = (TKey*)folderIterator()) ) {
	++count2;
	TObject *folderObject = folderKey->ReadObj(); 
	currentPlot = (TString)folderObject->GetName();
	if(verbose>2) std::cout << "plotObject #" << count2 << " ("+currentPlot+")" << std::endl;
	// check if object is a TH1 
	if(folderObject->InheritsFrom("TH1")){
	  if(verbose>2) std::cout << "inherits from TH1";
	  // check if TH2 and neglect because a simple 
	  // re-smearing is here not possible
	  if((folderObject->InheritsFrom("TH2"))&&(verbose>2)) std::cout << " and from TH2" << std::endl;
	  else{
	    if(verbose>2) std::cout << " and NOT from TH2" << std::endl;
	    // add to list of plots
	    if(verbose>2) std::cout << "will be added to list of output plots" << std::endl;
	    plotList_.push_back(currentFolder+"/"+currentPlot);
	  }
	}
      }
    }
  }
  // close data file after getting plot names
  files_[kData]->Close();
  // remove data file from list of considered files
  if(files_.count(kData )>0)files_.erase(kData );
  // print list of files considered
  if(verbose>0){
    std::cout << std::endl << "the following files will be considered: " << std::endl;
    // loop files
    for(int sample = kSig; sample<=kLast; sample++){
      // check existence of folder in all existing files
      if(files_.count(sample)>0){
	std::cout << files_[sample]->GetName() << std::endl;
      }
    }
  }
  // print out the name of all plots
  if(verbose>0){
    std::cout << std::endl << "list of all plots to be considered: " << std::endl;
    // loop list of plots
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      std::cout << "\"" << plotList_[plot] << "\"" << std::endl;
    }
    std::cout << plotList_.size() << " plots in total" << std::endl;
  }
  // ---------------------------------------
  // !!! load all hists !!!
  // ---------------------------------------
  unsigned int N1Dplots=plotList_.size();
  int Nplots=0;
  if(verbose>0) std::cout << std::endl << "loading plots: " << std::endl;
  // a) for std analysis file (& reweighted ttbar)
  getAllPlots(files_, plotList_, histo_, histo2_, N1Dplots, Nplots, verbose-1);
  if(verbose>0){
    std::cout << Nplots << " plots loaded from " << plotList_.size();
    std::cout << " requested" << std::endl << "empty plots are not counted" << std::endl;
  }
  // b) for zprime
  if(zprime){ 
    for(unsigned int plot=0; plot<plotList_.size(); ++plot){
      histo_[plotList_[plot]][kZprime] = (TH1F*)(files_[kZprime]->Get(plotList_[plot]));
    }
  }

  // ---------------------------------------
  // !!! definition of output file(name) !!!
  // ---------------------------------------
  //TString outputfile="/afs/naf.desy.de/user/g/goerner/WGspace/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/"+(TString)decayChannel+"PseudoData"+lum+"pb"+outNameExtension+"8TeV.root";
  TString outputfile="/afs/naf.desy.de/user/g/goerner/WGspace/RecentAnalysisRun8TeV_doubleKinFit/pseudodata/"+pseudoDataFileName(specifier, decayChannel);
  TFile* out = new TFile(outputfile, "recreate");
  if(verbose>0) std::cout << std::endl << "outputfile: " << outputfile << std::endl;
  poisson(histo_, plotList_, decayChannel, *out, luminosity, verbose, smear, useReweightedTop, avWeight, zprime, zPrimeLumiWeight);

  // free memory
  for(std::map< TString, std::map <unsigned int, TH1F*> >::iterator histo=histo_.begin(); histo!=histo_.end(); ++histo){
    for(std::map <unsigned int, TH1F*>::iterator histoSub=histo->second.begin(); histoSub!=histo->second.end(); ++histoSub){
      if(histoSub->second) delete histoSub->second;
    }
    histo->second.clear();
  }
  histo_ .clear();
  for(std::map< TString, std::map <unsigned int, TH2F*> >::iterator histo2=histo2_.begin(); histo2!=histo2_.end(); ++histo2){
    for(std::map <unsigned int, TH2F*>::iterator histo2Sub=histo2->second.begin(); histo2Sub!=histo2->second.end(); ++histo2Sub){
      if(histo2Sub->second) delete histo2Sub->second;      
    }
    histo2->second.clear();
  }
  histo2_ .clear();
  for(std::map<unsigned int, TFile*>::const_iterator file=files_.begin(); file!=files_.end(); ++file){
    if(file->second){
      file->second->Close();	  
      delete file->second;
    }
  }
  files_ .clear();
  out->Close();
  delete out;

}

void poisson(std::map< TString, std::map <unsigned int, TH1F*> > histo_, std::vector<TString> plotList_, const std::string decayChannel, TFile& outputfile, int luminosity, const unsigned int verbose, bool smear, bool useReweightedTop, double avReweight, bool useZprime, double zPrimeLumiWeight){
  unsigned int kCombined=111;
  
  if(verbose>1){
    std::cout << "running function poisson" << std::endl;
    std::cout << "-----------------------------" << std::endl;
    std::cout << "luminosity: " << luminosity << std::endl;
    std::cout << "decayChannel: " << decayChannel << std::endl;
    std::cout << "smear: " << smear << std::endl;
    std::cout << "useReweightedTop: " << useReweightedTop << std::endl;
    std::cout << "avReweight: " << avReweight << std::endl;
    std::cout << "useZprime: " << useZprime << std::endl;
    std::cout << "zPrimeLumiWeight: " << zPrimeLumiWeight << std::endl;
    sleep(5);
  }
  // upper bound for looping events: in/exclude zprime
  int kLast  =kSAToptW;
  int kZprime=kLast+1;
  if(useZprime) kLast=kZprime;

  // go to output file
  outputfile.cd();
  //  ---
  //     combine plot from all files   
  //  ---
  if(verbose>0) std::cout << "start to combine plots " << std::endl;
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    if(verbose>1){
      std::cout << "plot " << plotList_[plot] << std::endl;
      std::cout << "considered samples: " << std::endl;
    }
    // indicate first plot existing in sample
    bool first=true;
    // loop samples
    for(int sample = kSig; sample<=kLast; sample++){
      // check if plot exists
      if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(sample))){
	// take care of non existing plots
	if(!(useZprime&&sample==kZprime&&plotList_[plot].Contains("PUControlDistributions"))){
	  // do lumiweighting
	  if     (useZprime       &&sample==kZprime) histo_[plotList_[plot]][sample]->Scale(zPrimeLumiWeight);
	  else histo_[plotList_[plot]][sample]->Scale(lumiweight(sample, luminosity, sysNo, decayChannel));
	  // for reweighted ttbar: weight to keep same normalization
	  if(useReweightedTop&&(sample==kSig||sample==kBkg)) histo_[plotList_[plot]][sample]->Scale(1.0/avReweight);
	  // first subsample (should be ttbar signal)
	  // -> create histogram map entry
	  if(first){ 
	    histo_[plotList_[plot]][kCombined] = (TH1F*)(histo_[plotList_[plot]][sample]->Clone());
	    first=false;
	  }
	  // add other subsample
	  else{ 
	    histo_[plotList_[plot]][kCombined]->Add((TH1F*)(histo_[plotList_[plot]][sample]->Clone()),1.0);
	  }
	  if(verbose>1){ 
	    if(useZprime&&sample==kZprime) std::cout << "z prime, weight " << zPrimeLumiWeight << std::endl;
	    else if(useReweightedTop&&(sample==kSig||sample==kBkg))  std::cout << "reweighted "+sampleLabel(sample,decayChannel)+", weight " << 1.0/avReweight << std::endl;
	    else std::cout << sampleLabel(sample,decayChannel) << ", weight " << lumiweight(sample, luminosity, sysNo, decayChannel) << std::endl;
	  }
	}
      }
    }
  }

  //  ---
  //     poisson smearing   
  //  ---
  if(verbose>0) std::cout << "start to do the poisson smearing plots " << std::endl;
  // loop plots
  for(unsigned int plot=0; plot<plotList_.size(); ++plot){
    // check if plot exists
    if((histo_.count(plotList_[plot])>0)&&(histo_[plotList_[plot]].count(kCombined))){
      if(verbose>1) std::cout << "plot " << plotList_[plot] << std::endl;
      // do the poisson smearing
      if(smear){
	TRandom3 rnd(0);
	for(int ibin=0; ibin<=histo_[plotList_[plot]][kCombined]->GetNbinsX()+1; ++ibin){
	  int evts=histo_[plotList_[plot]][kCombined]->GetBinContent(ibin);
	  evts=rnd.Poisson(histo_[plotList_[plot]][kCombined]->GetBinContent(ibin));
	  histo_[plotList_[plot]][kCombined]->SetBinContent(ibin, evts);
	  histo_[plotList_[plot]][kCombined]->SetBinError(ibin, sqrt(evts));
	}
      }
      // do the saving
      // use the same directory as in the input files
      TString directory=getStringEntry(plotList_[plot], 1);
      TString plotName=getStringEntry(plotList_[plot], 2);
      if(verbose>1) std::cout << "check existence of directory " << directory << std::endl;
      if((outputfile.GetDirectory(directory)!=0)&&(verbose>1)) std::cout << "already existing" << std::endl;
      if(outputfile.GetDirectory(directory)==0){
	if(verbose>1) std::cout << "not existing - will create it" << directory << std::endl;
	outputfile.mkdir(directory);
      }
      outputfile.cd(directory);
      if(verbose>1) std::cout << "save plot" << directory << std::endl;
      // save plot
      histo_[plotList_[plot]][kCombined]->Write(plotName);
    }
  }
}

