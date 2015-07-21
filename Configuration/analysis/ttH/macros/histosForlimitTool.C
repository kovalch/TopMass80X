//////////////////////////////////////////////
//                                          //
//   Write by: Jasone Garay Garcia and      //
//             Christian Contreras-Campana  //
//   email: christian.contreras@desy.de     //
//   Date: 18.07.2015                       //
//                                          //
//////////////////////////////////////////////

// Description: Sowftware produce root file that is used in limit setting tool

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1.h>
#include <TObjArray.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <typeinfo>

#include <map>

#include "TSystem.h"
#include "TPRegexp.h"
#include "TClonesArray.h"
#include "TObjString.h"

#include <dirent.h>

#include <TList.h>
#include <TKey.h>
#include <TIterator.h>
#include <TObject.h>
#include <TObjString.h>

// Function handles producing root file with category, process, obserable, systematic mva information
void histosForlimitTool(TString Str = "mvaEventA_mvaWeight", TString outputFileName = "mvaEventA_mvaWeight")
{
  // Define regular expression string
  TString regexString = Str+"_(\\w+)_(\\w+)_(\\w+)_(\\w+)";

  // Types of systematic variations mapping to naming convention of ttH collaboration
  std::map<TString,TString> mapOfSystematics;

  mapOfSystematics["BTAGDISCR_BPURITY_DOWN"] = "CMS_CSVHF_Down";
  mapOfSystematics["BTAGDISCR_BPURITY_UP"]   = "CMS_CSVHF_Up";
  mapOfSystematics["BTAGDISCR_BSTAT2_DOWN"]  = "CMS_CSVHFStats2_Down";
  mapOfSystematics["BTAGDISCR_BSTAT2_UP"]    = "CMS_CSVHFStats2_Up";
  mapOfSystematics["BTAGDISCR_CERR2_DOWN"]   = "CMS_CSVCErr2_Down";
  mapOfSystematics["BTAGDISCR_CERR2_UP"]     = "CMS_CSVCErr2_Down";
  mapOfSystematics["BTAGDISCR_LSTAT1_DOWN"]  = "CMS_CSVLFStats1_Down";
  mapOfSystematics["BTAGDISCR_LSTAT1_UP"]    = "CMS_CSVLFStats1_Up";
  mapOfSystematics["BTAGDISCR_BSTAT1_DOWN"]  = "CMS_CSVHFStats1_Down";
  mapOfSystematics["BTAGDISCR_BSTAT1_UP"]    = "CMS_CSVHFStats1_Up";
  mapOfSystematics["BTAGDISCR_CERR1_DOWN"]   = "CMS_CSVCErr1_Down";
  mapOfSystematics["BTAGDISCR_CERR1_UP"]     = "CMS_CSVCErr1_Up";
  mapOfSystematics["BTAGDISCR_LPURITY_DOWN"] = "CMS_CSVLF_Down";
  mapOfSystematics["BTAGDISCR_LPURITY_UP"]   = "CMS_CSVLF_Up";
  mapOfSystematics["BTAGDISCR_LSTAT2_DOWN"]  = "CMS_CSVLFStats2_Down";
  mapOfSystematics["BTAGDISCR_LSTAT2_UP"]    = "CMS_CSVLFStats2_Up";
  mapOfSystematics["JES_DOWN"] = "scale_j_Down";
  mapOfSystematics["JES_UP"]   = "scale_j_Up";
  mapOfSystematics["Nominal"]  = "Nominal";

  // Additional systematics not included for the Baseline analysis
  mapOfSystematics["XSEC_TT2B_UP"]    = "CMS_XSEC_TT2B_Up";
  mapOfSystematics["XSEC_TT2B_DOWN"]  = "CMS_XSEC_TT2B_Down";
  mapOfSystematics["XSEC_TTCC_UP"]    = "CMS_XSEC_TTCC_Up";
  mapOfSystematics["XSEC_TTCC_DOWN"]  = "CMS_XSEC_TTCC_Down";
  mapOfSystematics["XSEC_TTH_UP"]     = "CMS_XSEC_TTH_Up";
  mapOfSystematics["XSEC_TTH_DOWN"]   = "CMS_XSEC_TTH_Down";
  mapOfSystematics["XSEC_TTOTHER_UP"] = "CMS_XSEC_TTOTHER_Up";
  mapOfSystematics["XSEC_TTOTHER_DOWN"] = "CMS_TTOTHER_Down";
  mapOfSystematics["XSEC_TTZ_UP"]       = "CMS_XSEC_TTZ_Up";
  mapOfSystematics["XSEC_TTZ_DOWN"]     = "CMS_XSEC_TTZ_Down";
  mapOfSystematics["FRAC_TTHF_UP"]      = "CMS_FRAC_TTHF_Up";
  mapOfSystematics["FRAC_TTHF_DOWN"]    = "CMS_FRAC_TTHF_Down";
  mapOfSystematics["FRAC_TTOTHER_UP"]   = "CMS_FRAC_TTOTHER_Up";
  mapOfSystematics["FRAC_TTOTHER_DOWN"] = "CMS_FRAC_TTOTHER_Down";
  mapOfSystematics["LUMI_UP"]   = "CMS_LUMI_Up";
  mapOfSystematics["LUMI_DOWN"] = "CMS_LUMI_Down";

  // Event categories based on the # of jets and # of b-tagged jets
  std::map<TString,TString> mapOfCategories;

  mapOfCategories["cate0"] = "dl_je3_te2";
  mapOfCategories["cate1"] = "dl_jge3_te3";
  mapOfCategories["cate2"] = "dl_jge4_te2";
  mapOfCategories["cate3"] = "dl_jge4_tge4";

  TFile* outputFile(0);

  // Set and change into to current working directory
  TDirectory *rootDir  = gDirectory;
  rootDir->cd();

  std::map<TString, TString>::iterator systematic;

  //std::string cwdPath = getcwd(NULL,0);
  DIR *pDIR1;
  DIR *pDIR2;
  DIR *pDIR3;

  struct dirent *entry1;
  struct dirent *entry2;
  struct dirent *entry3;

  // Begin iterating throug the Plots directory for File System hierarchy
  if((pDIR1=opendir("./Plots"))){

    while((entry1 = readdir(pDIR1))){

      if( strcmp(entry1->d_name, ".") != 0 && strcmp(entry1->d_name, "..") != 0 ) {
	//std::cout << entry1->d_name << "\n";

	if((pDIR2=opendir(std::string("./Plots/"+std::string(entry1->d_name)).c_str()))) {

	  while((entry2 = readdir(pDIR2))) {

	    if( strcmp(entry2->d_name, ".") != 0 && strcmp(entry2->d_name, "..") != 0 ) {
	      //std::cout << entry2->d_name << "\n";

	      if((pDIR3=opendir(std::string("./Plots/"+std::string(entry1->d_name)+"/"+std::string(entry2->d_name)).c_str()))) {

		while((entry3 = readdir(pDIR3))){
		  //std::cout << entry3->d_name << "\n";

		  TPRegexp subStr = TPRegexp(regexString+".root");

		  // check we access only root files of a particular type
		  if(TString(entry3->d_name).Contains(subStr)) {
		    
		    TObjArray *subStrArray = TPRegexp(regexString+".root").MatchS(TString(entry3->d_name));

		    // Extract information from filename string
		    const TString fileName  = ((TObjString *)subStrArray->At(0))->GetString();
		    const TString mvaConfig = ((TObjString *)subStrArray->At(1))->GetString();
		    const TString step = ((TObjString *)subStrArray->At(2))->GetString();
		    const TString cate = ((TObjString *)subStrArray->At(3))->GetString();

		    TString observableType = "observable";

		    // Store results in the macros directory
		    const TString outpath("./macros");                 		    

		    // Create input file handler
		    TFile* sampleFile = new TFile(TString("./Plots/")+TString(entry1->d_name)+"/"+TString(entry2->d_name)+"/"+fileName);

		    std::map<TString,TH1D*> mapOfHistograms;

		    // Check that file exists and is not corrupt
		    if(sampleFile->IsZombie()) {
		      std::cout<<"\n\tWe didn't find the "+TString(entry1->d_name)+" input!!\n";
		    }
		    else {		      
		      TString sample = (fileName.Copy()).ReplaceAll("_source.root","");
		      
		      // Update (i.e. append) to output file                                                                                           
		      outputFile = new TFile(outpath+"/"+outputFileName+".root","UPDATE");
		      
		      // Create list from list of histogram from root file
		      TList* list = sampleFile->GetListOfKeys();

		      if(false)list->Print();

		      if(!list) {
			std::cout << TString::Format("<Error> No keys found in file\n");
			exit(1);
		      }

		      // Create iterator from list
		      TIter next((TList*)list);
		      TKey* key;

		      // Iterate and extrat histograms from file
		      while ((key=(TKey*)next())){

			// String capture of process type
			TObjArray *subStrA = TPRegexp(regexString+".*").MatchS(TString(key->GetName()));
			const TString processType  = ((TObjString *)subStrA->At(4))->GetString();

			// Store only histograms not containing "canvas" label
			if (!processType.Contains("canvas"))
			  mapOfHistograms[processType]     = (TH1D*)sampleFile->Get(key->GetName());
		      }

		      // Check histogram exist
		      for (std::map<TString,TH1D*>::iterator iter = mapOfHistograms.begin(); iter != mapOfHistograms.end(); ++iter) {
			if(!mapOfHistograms[iter->first])
			std::cout << "\n\tWe didn't find the "+iter->first+" histogram!!\n";
		      }

		      // Access integral values to determine observation and rates
		      if (false) {

			std::cout << "File type: " << sample.Data() << std::endl;
			std::cout << TString::Format("The integral for data is %f\n", mapOfHistograms["data"]->Integral());
			std::cout << TString::Format("The integral for bkg is %f\n", mapOfHistograms["bkg"]->Integral());
			std::cout << TString::Format("The integral for signal is %f\n", mapOfHistograms["ttH_hbb"]->Integral());
		      }
		      
		      if(outputFile->GetDirectory(mapOfCategories[cate]+"__"+observableType))
			{
			  outputFile->cd(mapOfCategories[cate]+"__"+observableType);
			}
		      else {
			outputFile->mkdir(mapOfCategories[cate]+"__"+observableType,mapOfCategories[cate]+"__"+observableType);
			outputFile->cd(mapOfCategories[cate]+"__"+observableType);
		      }
		      
		      // Write histograms to file
		      for (std::map<TString,TH1D*>::iterator iter = mapOfHistograms.begin(); iter != mapOfHistograms.end(); ++iter) {
			
			TString label = TString(entry1->d_name).Contains("Nominal") ? iter->first : iter->first+"__"+mapOfSystematics[entry1->d_name];
			mapOfHistograms[iter->first]->Write(label,TObject::kOverwrite);
                      }

		      outputFile->Write("",TObject::kOverwrite);
		      outputFile->Close();		    
		    } // End section where histograms are written
		    
		    sampleFile->Close();
		  }
		}
	      }// Most innter directory
	      closedir(pDIR3);
	    }
	  }// Second most outer dirrectory
	  closedir(pDIR2);
	}
      }
    }// First most outter director
    closedir(pDIR1);
  }

  // Exit ROOT prompt
  std::cout << "**** Finish processing ****" << std::endl;
  gSystem->Exit(0);
}
