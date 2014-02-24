#include "basicFunctions.h"

bool isValidsample(unsigned int sample, unsigned int systematicVariation);

void treeComparison(double luminosity = 19712, bool save = true, int verbose=0, TString inputFolderName= "RecentAnalysisRun8TeV_doubleKinFit", TString dataFile= "/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/elecDiffXSecData2012ABCDAll.root:/afs/naf.desy.de/group/cms/scratch/tophh/RecentAnalysisRun8TeV_doubleKinFit/muonDiffXSecData2012ABCDAll.root", const std::string decayChannel = "combined", bool withRatioPlot=true, TString test="PV")
{
  // test= "prob" or "PV" 
  // data/MC -> MC/data
  bool invert=false;
  // linear fit in ratio?
  bool linFit=true;

  // ===================================
  // Define plotting order
  // ===================================
  std::vector<int> samples_;
  samples_.push_back(kSig);
  samples_.push_back(kBkg);
  samples_.push_back(kSTop);
  samples_.push_back(kWjets);
  samples_.push_back(kZjets);
  samples_.push_back(kDiBos);
  samples_.push_back(kQCD);
  samples_.push_back(kData);

  // ============================
  //  Set Root Style
  // ============================
		
  TStyle myStyle("HHStyle","HHStyle");
  setHHStyle(myStyle);
  myStyle.SetStripDecimals(true);
  myStyle.cd();
  gROOT->SetStyle("HHStyle");
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(2);

  // user specific configuration
  TString testQuantity=""; // name of separator in tree
  TString treePath=""; // path of tree
  // values for splitting topPt (Val0<=plot1<Val1<=plot2<Val2)
  std::vector< double > Val_;
  TString treeExt="";
  if(test=="PV"){
    testQuantity="nPV";
    treePath="compositedKinematicsKinFit/tree";
    Val_.push_back(0. );
    Val_.push_back(8. );
    Val_.push_back(13.);
    Val_.push_back(17.);
    Val_.push_back(22.);
    Val_.push_back(50.);
    treeExt="Fit";
  }
  if(test=="prob"){
    testQuantity="chi2";
    treePath="analyzeTopRecoKinematicsKinFit/tree";
    Val_.push_back(0. );
    Val_.push_back(1.386); // chi2<1.386 ~prob>0.50
    Val_.push_back(2.1);   // chi2<2.1 ~prob>0.35
    Val_.push_back(3.219); // chi2<3.219 ~prob>0.20
    Val_.push_back(7.824); // chi2<7.824 ~prob>0.02
    Val_.push_back(99999.);
    treeExt="";
  }
  // default configurations
  unsigned int systematicVariation=sysNo;
  TString ttbarMC="Madgraph";
  bool scaleTtbarToMeasured=true;
  // adjust luminosity and data files for combined control plots
  double luminosityEl=constLumiElec;
  double luminosityMu=constLumiMuon;
  if(!dataFile.Contains(":")){
    std::cout << "wrong input filenames, should be dataFileEl:dataFileMu, but is ";
    std::cout << dataFile << std::endl;
    exit(0);
  }
  TString dataFileEl=getStringEntry(dataFile,1 , ":");
  TString dataFileMu=getStringEntry(dataFile,42, ":");
  // file container
  std::map<unsigned int, TFile*> files_, filesMu_, filesEl_;
  // file vector storage
  std::vector< std::map<unsigned int, TFile*> > fileList_;
  // get analysis files
  TString inputFolder="/afs/naf.desy.de/group/cms/scratch/tophh/"+inputFolderName;
  if(verbose>0) std::cout << "loading files from " << inputFolder << std::endl;
  if(decayChannel!="combined"){
    TString dataFiletemp= decayChannel=="muon" ? dataFileMu : dataFileEl;
    files_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFiletemp, decayChannel, ttbarMC);
    fileList_.push_back(files_);
  }
  else{
    filesMu_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFileMu, "muon"    , ttbarMC);
    filesEl_ = getStdTopAnalysisFiles(inputFolder, systematicVariation, dataFileEl, "electron", ttbarMC);
    fileList_.push_back(filesMu_);
    fileList_.push_back(filesEl_);
  }
  // topPt histogram template
  if(verbose>0) std::cout << "creating temp histo" << std::endl;
  TH1F* temp= (TH1F*)(((TH1F*)(fileList_.at(0)[kSig]->Get("analyzeTopRecoKinematicsKinFit/topPt"))->Clone()));
  temp->Rebin(20);
  temp->Reset("icms");
  temp->SetTitle("");
  temp->GetXaxis()->SetTitle("p_{T}^{t} #left[GeV#right]");
  temp->GetYaxis()->SetTitle("Top quarks");
  //axesStyle(*temp, "p_{T}^{t} #left[GeV#right]", "norm. Top quarks", 0., 0.15);
  temp->GetXaxis()->SetRangeUser(0.,500.);
  //temp->GetYaxis()->SetRangeUser(0.,0.2 );
  temp->SetStats(kFALSE);
  temp->SetLineWidth(3);
  temp->SetMarkerSize(1.25);
  int binMax=temp->GetNbinsX()+1;
  // container for all histos
  std::map< TString, std::map <unsigned int, TH1F*> > histo_;
  
  // determine number of channels
  unsigned int nchannels = decayChannel=="combined" ? 2 : 1;
  
  // loop decay channels
  if(verbose>0) std::cout << "looping channels" << std::endl;
  for(unsigned int channel=0; channel<nchannels; ++channel){
    std::map<unsigned int, TFile*> tempfiles_=fileList_.at(channel);
    std::string tempChannel= decayChannel!="combined" ? decayChannel : (channel==0 ? "muon" : "electron");
    if(verbose>1) std::cout << " - " << tempChannel << std::endl;
    TString channelExt=getTStringFromInt(channel);
    // loop samples
    for(unsigned int sample=kSig; sample<=kSAToptW; ++sample){
      bool note=false;
      // check if sample is relevant
      if(isValidsample(sample, systematicVariation)){
	if(verbose>1){
	  std::cout << "  -> processing " << sampleLabel(sample, tempChannel);
	  std::cout << " (file " << tempfiles_[sample]->GetName() << ")" << std::endl;
	}
	// calculate luminosity event weight
	double lumi=decayChannel!="combined" ? luminosity : (channel==0 ? luminosityMu : luminosityEl);
	double lumiwgt=lumiweight(sample, lumi, systematicVariation, tempChannel);
	if(verbose>1) std::cout << "     (lumiweight=" << lumiwgt << ")" << std::endl;
	// get trees
	TTree* tree = (TTree*)(tempfiles_[sample]->Get(treePath));
	if(!tree){
	  std::cout << "     !ERROR: tree not found!" << std::endl;
	  exit(0);
	}
	else if(verbose>1)  std::cout << "     (tree found, contains " << tree->GetEntries() << " entries)" << std::endl;
	// container for values read from tree
	std::map< TString, float > value_;
	// initialize map entries with 0 
	value_["weight"  ]=1.;
	value_["testQuantity"    ]=0.;
	value_["topPtLep"]=0.;
	value_["topPtHad"]=0.;
	// initialize branches
	tree->SetBranchStatus ("*", 0);
	tree->SetBranchStatus ("weight"  , 1);
	tree->SetBranchStatus ("topPtLep"+treeExt, 1);
	tree->SetBranchStatus ("topPtHad"+treeExt, 1);
	tree->SetBranchStatus (testQuantity     , 1);
	tree->SetBranchAddress(testQuantity     ,(&value_["testQuantity"    ]));
	tree->SetBranchAddress("weight"  ,(&value_["weight"  ]));
	tree->SetBranchAddress("topPtLep"+treeExt,(&value_["topPtLep"]));
	tree->SetBranchAddress("topPtHad"+treeExt,(&value_["topPtHad"]));
	// initialize result plots
	histo_["topPt"+channelExt     ][sample]=(TH1F*)(temp->Clone());
	for(int plot=1; plot<(int)Val_.size(); ++plot){
	  histo_["topPtProb"+getTStringFromInt(plot)+channelExt][sample]=(TH1F*)(temp->Clone());
	}
	if(verbose>1) std::cout << "     -> looping tree" << std::endl;
	// loop all events to fill plots
	for(unsigned int event=0; event<tree->GetEntries(); ++event){
	  // get event
	  tree->GetEntry(event);
	  // check if values are reasonable
	  if(!((value_["weight"]>0&&value_["weight"]<10)||(test!="PV"&&value_["weight"]==0.))){ 
	    if(!note){ std::cout << "!!! WARNING - some weights are strange (e.g." << value_["weight"] << ") !!!"<< std::endl; note=true; }
	    value_["weight"]=1.0;
	  }
	  // get relevant quantities
	  double weight=value_["weight"]*lumiwgt;
	  double filterQuantity  =value_["testQuantity"  ];
	  double topPtLep=value_["topPtLep"];
	  double topPtHad=value_["topPtHad"];
	  if(verbose>2){
	    std::cout << "      event #" << event+1 << "/" << tree->GetEntries() << ":" << std::endl;
	    std::cout << "      weight=" << weight << ", " << "testQuantity" << "=" << filterQuantity << ", topPtLep=" << topPtLep << ", topPtHad=" << topPtHad << std::endl;
	  }
	  // fill histo for all
	  histo_["topPt"+channelExt][sample]->Fill(topPtLep, weight);
	  histo_["topPt"+channelExt][sample]->Fill(topPtHad, weight);
	  // fill histo for different ranges of the filterQuantity
	  for(int plot=1; plot<(int)Val_.size(); ++plot){
	    TString nameNr=getTStringFromInt(plot);
	    if(filterQuantity>=Val_[plot-1]&&filterQuantity<Val_[plot]){
	      histo_["topPtProb"+nameNr+channelExt][sample]->Fill(topPtLep, weight);
	      histo_["topPtProb"+nameNr+channelExt][sample]->Fill(topPtHad, weight);
	    }
	  } // end for loop separation values
	} // end for loop tree events
      } // end if is valid sample
    } // end for loop samples
  } // end for loop decay channels

  // create final plots
  // -> combine decay channels and MC samples
  unsigned int kAllMC=42;
  // loop samples
  if(verbose>0) std::cout << "combining decay channels and MC samples" << std::endl;
  for(unsigned int sample=kSig; sample<=kSAToptW; ++sample){
    // check if sample is relevant
    if(isValidsample(sample, systematicVariation)){
      // loop decay channels
      for(unsigned int channel=0; channel<nchannels; ++channel){
	std::string tempChannel= decayChannel!="combined" ? decayChannel : (channel==0 ? "muon" : "electron");
	if(verbose>1) std::cout << " -> processing " << sampleLabel(sample, tempChannel) << "(" << tempChannel << ")" << std::endl;
	TString channelExt=getTStringFromInt(channel);
	// get plots for current channels
	std::vector <TH1F*> tempHist_;
        tempHist_.push_back((TH1F*)histo_["topPt"+channelExt     ][sample]->Clone());
	for(int plot=1; plot<(int)Val_.size(); ++plot){
	  TString nameNr=getTStringFromInt(plot);
	  tempHist_.push_back((TH1F*)histo_["topPtProb"+nameNr+channelExt][sample]->Clone());
	}
	std::vector <int> nevents_;
        for(int plot=0; plot<(int)Val_.size(); ++plot){
	  nevents_.push_back(tempHist_[plot]->Integral(0,binMax));
	}
	// add all channels for final histogram
	if(channel==0){
	  for(int plot=0; plot<(int)Val_.size(); ++plot){
	    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	    histo_["topPt"+nameNr][sample]=(TH1F*)tempHist_[plot]->Clone();
	  }
	}
	else{
          for(int plot=0; plot<(int)Val_.size(); ++plot){
            TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	    histo_["topPt"+nameNr][sample]->Add((TH1F*)tempHist_[plot]->Clone());
	  }
	}
	if(verbose>2){
	  std::cout << "    (#events(" << tempChannel << ")=";
	  for(int plot=0; plot<(int)Val_.size(); ++plot){
	    std::cout <<  nevents_[plot];
	    if(plot<(int)Val_.size()-1) std::cout << " / ";
	  }
	  std::cout << ")" << std::endl;
	}
      } // end channel for loop
      std::vector <double> neventsComb_;
      for(int plot=0; plot<(int)Val_.size(); ++plot){
	TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	neventsComb_.push_back(histo_["topPt"+nameNr][sample]->Integral(0,binMax));
      }
      if(verbose>1){
	std::cout << "    (#events=";
	for(int plot=0; plot<(int)Val_.size(); ++plot){
	  std::cout <<  neventsComb_[plot];
	  if(plot<(int)Val_.size()-1) std::cout << " / ";
	}
	std::cout << ")" << std::endl;
      }
      // combine all MC samples
      if(sample!=kData){
	if(sample==kSig){
	  for(int plot=0; plot<(int)Val_.size(); ++plot){
	    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	    histo_["topPt"+nameNr][kAllMC]=(TH1F*)histo_["topPt"+nameNr][sample]->Clone();
	  }
	}
	else{
	  for(int plot=0; plot<(int)Val_.size(); ++plot){
            TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	    histo_["topPt"+nameNr][kAllMC]->Add((TH1F*)histo_["topPt"+nameNr][sample]->Clone());
	  }
	}
	// MC histogram style
	for(int plot=0; plot<(int)Val_.size(); ++plot){
	  TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	  histogramStyle(*histo_["topPt"+nameNr][sample], sample, true);
	  //histo_["topPt"+nameNr][sample]->SetLineColor(histo_["topPt"+nameNr][sample]->GetFillColor());
	  histo_["topPt"+nameNr][sample]->SetLineWidth(1);
	}
      }
      else{
	// data histogram style
        for(int plot=0; plot<(int)Val_.size(); ++plot){
          TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	  histogramStyle(*histo_["topPt"+nameNr][sample], kData, false);
	}
      }
    } // end if sample is valid
  } // end sample for loop

  // print some interesting numbers
  // chosen slices
  std::vector< double >neventsMC_;
  if(verbose>0){
    std::cout << "slices in " << testQuantity << ":  all / ";
    for(int plot=1; plot<(int)Val_.size(); ++plot){
      std::cout <<  "[" << Val_[plot-1] << ".." << Val_[plot] << "]";
      if(plot<(int)Val_.size()-1) std::cout << "/ ";
    }
    std::cout << std::endl;
  }
  // total number of MC events
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
    neventsMC_.push_back(histo_["topPt"+nameNr][kAllMC]->Integral(0,binMax));
  }
  if(verbose>0){
    std::cout << "#events( MC )=";
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      std::cout <<  neventsMC_[plot];
      if(plot<(int)Val_.size()-1) std::cout << ", ";
    }
    std::cout << std::endl;
  }
  // total number of data events
  std::vector< double >neventsData_;
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
    neventsData_.push_back(histo_["topPt"+nameNr][kData]->Integral(0,binMax));
  }
  if(verbose>0){
    std::cout << "#events(Data)=";
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      std::cout <<  neventsData_[plot];
      if(plot<(int)Val_.size()-1) std::cout << " / ";
    }
    std::cout << std::endl;
  }
  // data over MC ratio
  if(verbose>0){
    std::cout << "(data/MC ratio=";
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      std::cout <<  neventsData_[plot]/neventsMC_[plot];
      if(plot<(int)Val_.size()-1) std::cout << " / ";
    }
    std::cout << ")" << std::endl;
  }

  // scale ttbar to match total number of events
  if(scaleTtbarToMeasured){
    if(verbose>0) std::cout << "scale ttbar component to match #data events " << std::endl;
    std::vector<double>neventsTop_, SFTop_;
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
      neventsTop_.push_back(histo_["topPt"+nameNr][kSig]->Integral(0,binMax)+histo_["topPt"+nameNr][kBkg]->Integral(0,binMax));
      SFTop_.push_back((neventsTop_[plot]+(neventsData_[plot]-neventsMC_[plot]))/neventsTop_[plot]);
    }
    // scale combined and top MC plots
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
      // subtract ttbar from all MC
      histo_["topPt"+nameNr][kAllMC]->Add((TH1F*)(histo_["topPt"+nameNr][kSig]->Clone()) , -1.);
      histo_["topPt"+nameNr][kAllMC]->Add((TH1F*)(histo_["topPt"+nameNr][kBkg]->Clone()) , -1.);
      // scale ttbar 
      histo_["topPt"+nameNr][kSig]->Scale(SFTop_[plot]);
      histo_["topPt"+nameNr][kBkg]->Scale(SFTop_[plot]);
      // re-add ttbar MC
      histo_["topPt"+nameNr][kAllMC]->Add((TH1F*)(histo_["topPt"+nameNr][kSig]->Clone()) );
      histo_["topPt"+nameNr][kAllMC]->Add((TH1F*)(histo_["topPt"+nameNr][kBkg]->Clone()) );
    }
    // printout
    std::vector<double>neventsMCscaled_;
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
      neventsMCscaled_.push_back(histo_["topPt"+nameNr][kAllMC]->Integral(0,binMax));
    }
    if(verbose>1){
      std::cout << "#events(scaledMC)=";
      for(int plot=0; plot<(int)Val_.size(); ++plot){
	std::cout << neventsMCscaled_[plot];
        if(plot<(int)Val_.size()-1) std::cout << " / ";
      }
      std::cout << std::endl;
      std::cout << "(data/scaledMC ratio=";
      for(int plot=0; plot<(int)Val_.size(); ++plot){
	std::cout <<  neventsData_[plot]/neventsMCscaled_[plot];
	if(plot<(int)Val_.size()-1) std::cout << " / ";
      }
      std::cout << ")" << std::endl;
    }
  };

  // combine single top subsamples
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
    histo_["topPt"+nameNr][kSTop]=(TH1F*)histo_["topPt"+nameNr][kSAToptW]->Clone();
    for(int sample=(int)kSTops; sample<(int)kSAToptW; ++sample){
      // add to combined STop
      histo_["topPt"+nameNr][kSTop]->Add((TH1F*)histo_["topPt"+nameNr][sample]->Clone());
    } // end for loop single top subsamples
  } // end for loop plots

  // create MC histo stack plots
  int lastSample=-1;
  // loop samples
  if(verbose>0) std::cout << "creating MC stack histos" << std::endl;
  for(int sampleOri=(int)kDiBos; sampleOri>=(int)kSig; --sampleOri){
    // use previous defined order
    int sampleMod=samples_[sampleOri];
    if(verbose>1) std::cout << "processing " << sampleLabel(sampleMod, decayChannel) << "(= " << sampleMod <<", last sample=" << lastSample << ")" << std::endl;
    // exclude QCD and Diboson
    if(sampleMod!=kQCD&&sampleMod!=kDiBos){
      if(lastSample>-1){
	if(verbose>1) std::cout << "adding " << sampleLabel(lastSample, decayChannel) << " to " << sampleLabel(sampleMod, decayChannel) << std::endl;	 	
	// loop plots
	for(int plot=0; plot<(int)Val_.size(); ++plot){
	  TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
	  if(verbose>1) std::cout << "processing " << "topPt"+nameNr << std::endl;
	  // add to stack
	  if(!histo_.count("topPt"+nameNr)>0) std::cout << "WARNING: topPt"+nameNr+" does not exist in histo_!" << std::endl;
	  else if(!histo_["topPt"+nameNr].count(sampleMod )>0) std::cout << "WARNING: sample " << sampleMod << " does not exist in histo_[topPt"+nameNr+"]!" << std::endl;
	  else if(!histo_["topPt"+nameNr].count(lastSample)>0) std::cout << "WARNING: sample " << lastSample << " does not exist in histo_[topPt"+nameNr+"]!" << std::endl;
	  histo_["topPt"+nameNr][sampleMod]->Add((TH1F*)histo_["topPt"+nameNr][lastSample]->Clone());
	} // end for loop plots
	if(verbose>1) std::cout << "done" <<  std::endl;
      } // if not last sample
      if(verbose>1) std::cout << "lastSample set to " << sampleMod << std::endl;
      lastSample=sampleMod;
    } // end else if !QCD
  } // end for loop original sample ordering

  // printout
  std::vector<double>neventsMCstack_;
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
    neventsMCstack_.push_back(histo_["topPt"+nameNr][kAllMC]->Integral(0,binMax));
  }
  if(verbose>1){
    std::cout << "#events(stack MC) =";
    for(int plot=0; plot<(int)Val_.size(); ++plot){
      std::cout << neventsMCstack_[plot];
      if(plot<(int)Val_.size()-1) std::cout << " / ";
    }
    std::cout << std::endl;
  }

  // all MC histogram style
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
    histo_["topPt"+nameNr][kAllMC]->SetLineColor(kBlue);
    histo_["topPt"+nameNr][kAllMC]->SetMarkerColor(kBlue);
    histo_["topPt"+nameNr][kAllMC]->SetLineStyle(1);
  }
  
  std::vector<double> zeroerr_;
  for(int bin=0; bin<histo_["topPt"][kAllMC]->GetNbinsX(); ++bin) zeroerr_.push_back(0);

  // normalization
  // -> not done at the moment, scaled wrt Ndata for each plot separately

  //  Create canvas
  if(verbose>0)  std::cout << "creating canvas" << std::endl;
  std::vector<TCanvas*> plotCanvas_;
  for(unsigned int sample=0; sample<Val_.size(); sample++){
    addCanvas(plotCanvas_);
  }

  // create legends
  if(verbose>0)  std::cout << "creating legend" << std::endl;
  std::vector< TLegend* > leg_;
  TLegend *leg= new TLegend(0.73, 0.5, 0.91, 0.88);
  legendStyle(*leg,"");
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr   = plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);    
    TString sVallow  = plot==0 ? "" : getTStringFromDouble(Val_[plot-1], getRelevantDigits(Val_[plot-1]));
    TString sValhigh = plot==0 ? "" : getTStringFromDouble(Val_[plot  ], getRelevantDigits(Val_[plot ]));
    TString legHeader= plot==0 ? "KinFit, all" : sVallow+"#leq"+testQuantity+"<"+sValhigh;
    TLegend *templeg=(TLegend*)leg->Clone(); 
    templeg ->SetHeader(legHeader);
    templeg ->AddEntry(histo_["topPt"+nameNr][kData ], "data"   , "P");
    //templeg ->AddEntry(histo_["topPt"+nameNr][kAllMC], "all MC" , "L" );
    for(unsigned int sample=kSig; sample<kData; sample++){
      unsigned int sampleMod=samples_[sample];
      if(sampleMod!=kQCD&&sampleMod!=kDiBos) templeg ->AddEntry(histo_["topPt"+nameNr][sampleMod], sampleLabel(sampleMod, decayChannel), "F" );
    }
    leg_.push_back((TLegend*)(templeg->Clone()));
  }

  int canvasNumber=0;
  // do the plotting 
  if(verbose>0)  std::cout << "plotting " << std::endl;
  for(int plot=0; plot<(int)Val_.size(); ++plot){
    TString nameNr= plot==0 ? "" : "TestSlice"+getTStringFromInt(plot);
    TString title = plot==0 ? "topPtKinFit"+TString(decayChannel) : "topPt"+testQuantity+getTStringFromInt(plot)+TString(decayChannel);
    plotCanvas_[canvasNumber]->cd(0);
    plotCanvas_[canvasNumber]->SetTitle(title);
    // drawing
    // plots
    histo_["topPt"+nameNr][kSig]->SetMaximum(1.3*histo_["topPt"+nameNr][kData]->GetMaximum());
    histo_["topPt"+nameNr][kSig]->GetXaxis()->SetNoExponent(true);
    histo_["topPt"+nameNr][kSig]->GetYaxis()->SetNoExponent(true);
    histo_["topPt"+nameNr][kSig]->Draw("axis");
    // loop samples
    for(int sampleOri=(int)kSig; sampleOri<=(int)kDiBos; ++sampleOri){
      // use previous defined order
      int sampleMod=samples_[sampleOri];
      if(verbose>2) std::cout << "processing sample " << sampleMod << " ("+sampleLabel(sampleMod, decayChannel)+")" << std::endl; 
      // draw other MC samples, excluding QCD 
      if(sampleMod!=kQCD&&sampleMod!=kDiBos){
	if(verbose>2) std::cout << "-> drawing!" << sampleMod << std::endl;
	histo_["topPt"+nameNr][sampleMod]->Draw("hist same");	
      }
    } // end for loop ori samples
    //histo_["topPt"+nameNr][kAllMC]->Draw("hist same");
    histo_["topPt"+nameNr][kData ]->Draw("ep same");
    // legend
    leg_[plot]->Draw("same");
    // add labels for decay channel, luminosity, energy and CMS preliminary (if applicable)
    if      (decayChannel=="muon"    ) DrawDecayChLabel("#mu + Jets");
    else if (decayChannel=="electron") DrawDecayChLabel("e + Jets");
    else                               DrawDecayChLabel("e/#mu + Jets Combined");  
    DrawCMSLabels(true,luminosity);
    // draw ratio
    if(withRatioPlot){
      // labels of ratio
      TString ratioLabelNominator  ="N_{MC}";
      TString ratioLabelDenominator="N_{Data}";
      double ratMin= invert ? 0.75 : 0.30;
      double ratMax= invert ? 1.75 : 1.29;
      std::vector<double> err_;
      for(int bin=1; bin<histo_["topPt"+nameNr][kSig]->GetNbinsX(); ++bin){
	double ratio = histo_["topPt"+nameNr][kData]->GetBinContent(bin)/histo_["topPt"+nameNr][kSig]->GetBinContent(bin);
	if(invert) ratio=1./ratio;
	double val=ratio*(histo_["topPt"+nameNr][kData]->GetBinError(bin)/histo_["topPt"+nameNr][kData]->GetBinContent(bin));    
	if(val<0||val>histo_["topPt"+nameNr][kData]->GetBinContent(bin)) val=1.;
	err_.push_back(val);
      }      
      int rval1 = drawRatio(histo_["topPt"+nameNr][kData], histo_["topPt"+nameNr][kSig], ratMin, ratMax, myStyle, verbose, err_, ratioLabelNominator, ratioLabelDenominator, "p e", kBlack, true, 0.5, 505, invert, true, linFit);
      if (rval1!=0) std::cout << " Problem occured when creating ratio plot for " << nameNr << std::endl;
    }
    canvasNumber++;
  }
  
  // saving
  if(verbose>0) std::cout << "saving" << std::endl;
  if(save){
    TString outfolder="./diffXSecFromSignal/plots/combined/2012/topPtTest/";
    // eps and png
    if(verbose==0) gErrorIgnoreLevel=kWarning;
    saveCanvas(plotCanvas_, outfolder, "topPtTest"+testQuantity+TString(decayChannel), true, true, true);
  }
}

bool isValidsample(unsigned int sample, unsigned int systematicVariation){
  // exclude combined samples as they contain no tree
  // and ENDOFENUMs as they represent no sample
  if(sample==kQCD||sample==kSTop||sample==kDiBos||sample==ENDOFSAMPLEENUM||sample==ENDOFSAMPLEENUM2) return false;
  // exclude all QCD samples (statistical fluctuations)
  if(sample>=kQCDEM1&&sample<=kQCDBCE3  ) return false;
  // exclude systematic signal samples
  if(sample>=kSigPow&&sample<=kBkgPowHer) return false;
  // excude splitted single top samples for non scale variation
  if(!(systematicVariation==sysSingleTopScaleUp||systematicVariation==sysSingleTopScaleDown)&&sample>=kSToptW1&&sample<=kSAToptW3) return false;
  return true;
}
