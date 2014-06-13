#include "basicFunctions.h"
#include "BCC.C"

void applyBCCs(TString inputFolderName=AnalysisFolder, bool save=false, TString decayChannel="combined", bool extrapolate=true, bool hadron=false, bool mergeLepHad=true, bool addCrossCheckVariables=false)
{
    // INFO: this is the macro part to process the variables needed for differenetial cross sections

    // name of file to calculate BCC from
    TString inputFile=TopFilename(kSig, 0, "muon");
    if(decayChannel=="electron") inputFile.ReplaceAll("muon", "elec");
    else if (decayChannel=="combined"){
      inputFile.ReplaceAll("muon", "combined");
      //inputFile.ReplaceAll(".root", "Large.root");
    }
    inputFile=groupSpace+inputFolderName+"/"+inputFile;
    std::cout << std::endl << "input file: " << inputFile << std::endl;
  
    // name of file to add calculted BCC x values to
    TString outputFile="diffXSecBCC.root";
    std::cout << "output file: " << outputFile << std::endl;
    TString LV="Parton";
    if(!extrapolate&&hadron) LV="Hadron";
    TString PS="";
    if(!extrapolate) PS="PhaseSpace";
      
    // list of variables of interest
    std::vector<TString> xSecVariables_;
    // a) top and ttbar quantities
    //if(!hadron){
    xSecVariables_.push_back("ttbarDelPhi");
    
    //xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFit    , xSecVariablesKinFit     + sizeof(xSecVariablesKinFit    )/sizeof(TString)-1);
    //xSecVariables_.insert(xSecVariables_.end(), xSecVariablesKinFitNorm, xSecVariablesKinFitNorm + sizeof(xSecVariablesKinFitNorm)/sizeof(TString)-1);
      //}
    // b) lepton and b-jet quantities
    //if(hadron||!extrapolate){
    //xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalState    , xSecVariablesFinalState     + sizeof(xSecVariablesFinalState    )/sizeof(TString));
    //xSecVariables_.insert(xSecVariables_.end(), xSecVariablesFinalStateNorm, xSecVariablesFinalStateNorm + sizeof(xSecVariablesFinalStateNorm)/sizeof(TString));
    //}
    // c) cross check variables (only available for parton level cross-sections)
    //if (addCrossCheckVariables && !hadron){
    //  xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVar,     xSecVariablesCCVar     + sizeof(xSecVariablesCCVar)/sizeof(TString)    );
    //  xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesCCVarNorm, xSecVariablesCCVarNorm + sizeof(xSecVariablesCCVarNorm)/sizeof(TString));
    //}
    //xSecVariables_.insert( xSecVariables_.end(),   xSecVariablesIncl,      xSecVariablesIncl      + sizeof(xSecVariablesIncl)/sizeof(TString)     );
    // collect tree folder names and branchnames of variables for BCC
    std::vector<TString> xSecVariableBranchNames_, folderNamesBCC_;
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      if(!xSecVariables_[i].Contains("Norm")&&xSecVariables_[i]!="inclusive"&&!xSecVariables_[i].Contains("Star")){
	// collect names of branch entries and names of cross section folder(s)
	if(xSecVariables_[i]=="topPt" || xSecVariables_[i]=="topY"){ 
	  xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Had");
	  xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Lep");
	  folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+PS);
	  folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+PS);
	}
	else if(xSecVariables_[i].Contains("bq")){
	  if(!extrapolate&&hadron){
	    xSecVariableBranchNames_.push_back(xSecVariables_[i]+"GenLep");
	    TString antib=xSecVariables_[i]+"GenHad";
	    antib.ReplaceAll("bq", "bbarq");
	    xSecVariableBranchNames_.push_back(antib);
	    folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+"Bjets"+PS);
	    folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+"Bjets"+PS);
	  }
	  else{
	    xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Had");
	    xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Lep");
	    folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+PS);
	    folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+PS);
	  }
	}
	else if(xSecVariables_[i].Contains("bbbar")){
	  xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Gen");
	  folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+"Bjets"+PS);
	}
	else if(xSecVariables_[i].Contains("lep")&&!extrapolate&&hadron&&!xSecVariables_[i].Contains("Plus")&&!xSecVariables_[i].Contains("Minus")){
	  xSecVariableBranchNames_.push_back(xSecVariables_[i]+"Gen");
	  folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+"Lepton"+PS);
	}
	else{
	  xSecVariableBranchNames_.push_back(xSecVariables_[i]);
  	  folderNamesBCC_.push_back("analyzeTop"+LV+"LevelKinematics"+PS); 
	}
      }
    }
    std::cout << "variables of interest: " << std::endl;
    for(unsigned int i=0; i<folderNamesBCC_.size(); ++i){
      std::cout << folderNamesBCC_[i] << "/tree/" << xSecVariableBranchNames_[i] << std::endl;
    }

    // run calculation
    BCC b(inputFile, folderNamesBCC_, xSecVariableBranchNames_, mergeLepHad, addCrossCheckVariables);
    b.runBCCCalculation();
    // container for results of BCC
    std::map<TString, std::vector<double> > correctedCenters_= b.getMapWithCorrectedCentersInX();
    std::map<TString, std::vector<double> > corrCenterErrors_= b.getMapWithCenterErrorsInX();
    for(unsigned int i=0; i<xSecVariables_.size(); ++i){
      std:: cout << xSecVariables_[i] << std::endl;
      for(int bin=0; bin<(int)correctedCenters_[i].size(); ++bin){
	std::cout << correctedCenters_[i].at(bin) << std::endl;
      }
    }
    // do saving
    if(save){
      std::cout << "saving:" << std::endl;
      // loop variables
      for(unsigned int i=0; i<xSecVariables_.size(); ++i){
	TString candidate=xSecVariables_[i].ReplaceAll("Norm","");
	std::cout << "candidate " << candidate;
	if(correctedCenters_.count(candidate)!=0){
	  std::cout << std::endl << "-> ";
	  // convert map to TArrayD
	  TVectorD* temp=new TVectorD(correctedCenters_[candidate].size(), correctedCenters_[candidate].data());
	  // cross check
// 	  for(int vec=0; vec<temp->GetNrows();++vec){	    
// 	    double wert=temp->operator[](vec);
// 	    if(vec==0) std::cout << "size: "<< correctedCenters_[candidate].size() << std::endl;
// 	    std::cout << correctedCenters_[candidate][vec] << "->" << wert << std::endl;
// 	  }
          // save to file
	  saveToRootFile(outputFile, temp, true, 1, "BCC/"+LV+PS+"/"+candidate);
	}
	else std::cout << "-> not found" << std::endl;
      }
    }
}
