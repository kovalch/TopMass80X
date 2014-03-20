#include "iostream"
#include "string"
#include "vector"
#include "iomanip"
#include "fstream"


#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TError.h>



/*
 *
 *  run this macro via:
 *  root -l -b -q printTreeComparison.cc++
 *
*/


std::string basedir ("/data/group/top/DiffXS/");
std::string dir1 ("2014_02_10_TAG_N004");
std::string dir2 ("2014_03_14_TAG_N005");
std::string variable_array[] = {"met.pt()", "leptons.pt()", "jets.pt()"};

int nrEntries(TTree *tree, TString variable)
{
    TH1I histo = TH1I("histo", "histo", 100, 0, 100);
    tree->Draw(variable+">>histo");
    return histo.GetEntries();
}


bool checkFile(const TString &file0)
{
    ifstream _file0(file0);
    if(!_file0.is_open()){
        std::cout<<"File:\n "<<file0<<std::endl;
        std::cout<<" doesn't exist or cannot be opened.\nContinue"<<std::endl;
        return 0;
    };
    _file0.close();
    return 1;
}

void printTreeComparison()
{

    gErrorIgnoreLevel = 1001;

    std::string samples_array[] = {
        /*
         * default data and bg simulations
         */
        "ee_run2012A","ee_run2012B","ee_run2012C","ee_run2012D",
        "emu_run2012A","emu_run2012B","emu_run2012C","emu_run2012D",
        "mumu_run2012A","mumu_run2012B","mumu_run2012C","mumu_run2012D","qcdbcem2030",
        "dy1050","dy50inf",
        "qcdbcem3080","qcdbcem80170",
        "qcdem2030","qcdem3080","qcdem80170",
        "qcdmu15",
        "wtolnu",
        "singleantitop_tw","singletop_tw",
        "wwtoall","wztoall","zztoall",
        "ttbarbg","ttbarsignalplustau",
        /*
         * ttbar signal variations
         */
        "ttbarbg_massdown","ttbarbg_massup",
        "ttbarsignalplustau_massdown","ttbarsignalplustau_massup",
        "ttbarbg_matchingdown","ttbarbg_matchingup",
        "ttbarsignalplustau_matchingdown","ttbarsignalplustau_matchingup",
        "ttbarbg_scaledown","ttbarbg_scaleup",
        "ttbarsignalplustau_scaledown","ttbarsignalplustau_scaleup",
        /*
         * other ttbar generators
         */
        "ttbarbg_mcatnlo", "ttbarsignalplustau_mcatnlo",
        "ttbarbg_powheg","ttbarsignalplustau_powheg",
        "ttbarbg_powhegHerwig","ttbarsignalplustau_powhegHerwig",
        "ttbarsignalplustau_Perugia11","ttbarsignalplustau_Perugia11NoCR",
        "ttbarbg_HadronicMadgraphWithSpinCorrelation","ttbarbg_SemiLeptMadgraphWithSpinCorrelation","ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation",
        /*
         * ttH samples
         */
        "ttbarH125inclusive","ttbarH125tobbbar",
        "ttbarW","ttbarZ","ttgjets",
    };

    std::vector<std::string> samples(samples_array, samples_array+sizeof(samples_array)/sizeof(samples_array[0]));
    std::vector<std::string> variable (variable_array, variable_array+sizeof(variable_array)/sizeof(variable_array[0]));

    std::cout<<"\n\nComparing trees from: "<<std::endl;
    std::cout<<"   "<<basedir+dir1<<std::endl;
    std::cout<<"   "<<basedir+dir2<<std::endl;
    std::cout<<"\n\n"<<std::endl;
    
    for (size_t vars=0; vars<variable.size(); vars++){
        std::cout<<"\n--------------------------------------------------------------------------------------"<<std::endl;
        std::cout<<std::setw(40)<<"Sample "<<std::setw(15)<<dir1.c_str()<<" "<<std::setw(15)<<dir2.c_str()<<std::setw(15)<<" diff(%)  variable: "<<variable.at(vars)<<std::endl;
        for (size_t iter = 0 ; iter<samples.size(); iter++){
            TString file0=basedir+dir1+"/"+samples.at(iter)+".root";
            TString file1=basedir+dir2+"/"+samples.at(iter)+".root";

            if(!checkFile(file0)) continue;
            if(!checkFile(file1)) continue;

            TFile *f0 = new TFile(file0);
            TFile *f1 = new TFile(file1);

            TTree *t0 = (TTree*)f0->Get("writeNTuple/NTuple");
            TTree *t1 = (TTree*)f1->Get("writeNTuple/NTuple");

            int nrEvents_old = nrEntries(t0, variable.at(vars));
            int nrEvents_new = nrEntries(t1, variable.at(vars));

            std::cout<<std::setw(40)<<samples.at(iter).c_str()<<std::setw(15)<<nrEvents_old<<std::setw(15)<<nrEvents_new<<std::setw(15)<<100*(1.*nrEvents_old-nrEvents_new)/nrEvents_old<<std::endl;

            delete t0; delete t1;
            f0->Close(); f1->Close();
            delete f0;   delete f1;
        }
    }

    std::cout<<"\n--------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<std::setw(40)<<"Sample "<<std::setw(15)<<dir1.c_str()<<" "<<std::setw(15)<<dir2.c_str()<<std::setw(15)<<" diff(%)  variable: unweightedEvents"<<std::endl;
    for (size_t iter = 0 ; iter<samples.size(); iter++){
        TString file0=basedir+dir1+"/"+samples.at(iter)+".root";
        TString file1=basedir+dir2+"/"+samples.at(iter)+".root";

        if(!checkFile(file0)) continue;
        if(!checkFile(file1)) continue;

        TFile *f0 = new TFile(file0);
        TFile *f1 = new TFile(file1);

        TH1I *t0 = (TH1I*)f0->Get("EventsBeforeSelection/unweightedEvents");
        TH1I *t1 = (TH1I*)f1->Get("EventsBeforeSelection/unweightedEvents");

        int nrEvents_old = t0->GetEntries();
        int nrEvents_new = t1->GetEntries();

        std::cout<<std::setw(40)<<samples.at(iter).c_str()<<std::setw(15)<<nrEvents_old<<std::setw(15)<<nrEvents_new<<std::setw(15)<<100*(1.*nrEvents_old-nrEvents_new)/nrEvents_old<<std::endl;
        delete t0; delete t1;
        f0->Close(); f1->Close();
        delete f0;   delete f1;
    }
}
