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
 *  > root -l -b printTreeComparison.C
 *  > .L printTreeComparison.C++
 *  > addVariable("variable to count")
 *  > printTreeComparison("folder with ntupleA", "folder with ntupleB")
 *
*/


int nrEntries(TTree *tree, TString variable)
{
    TH1I histo = TH1I("histo", "histo", 100, 0, 100);
    tree->Draw(variable+">>histo");
    return histo.GetEntries();
}


bool checkFile(const TString &file0, ofstream& outputFile)
{
    ifstream _file0(file0);
    if(!_file0.is_open()){
        outputFile<<"File:\n "<<file0<<std::endl;
        outputFile<<" doesn't exist or cannot be opened.\nContinue"<<std::endl;
        return 0;
    };
    _file0.close();
    return 1;
}



std::vector<std::string> v_variable(1, "runNumber");



void addVariable(const std::string& variable){
    v_variable.push_back(variable);
}



void printTreeComparison(const TString& folder1, const TString& folder2)
{

    gErrorIgnoreLevel = 1001;
    
    std::string variable_array[v_variable.size()];
    for(size_t iVariable = 0; iVariable < v_variable.size(); ++iVariable){
        variable_array[iVariable] = v_variable.at(iVariable);
    }
    
    std::string samples_array[] = {
        
        // default data and bg simulations
        "ee_run2012A","ee_run2012B","ee_run2012C","ee_run2012D",
        "emu_run2012A","emu_run2012B","emu_run2012C","emu_run2012D",
        "mumu_run2012A","mumu_run2012B","mumu_run2012C","mumu_run2012D","qcdbcem2030",
        "dy1050","dy50inf",
        "qcdbcem2030", "qcdbcem3080","qcdbcem80170",
        "qcdem2030","qcdem3080","qcdem80170",
        "qcdmu15",
        "wtolnu",
        "singleantitop_tw","singletop_tw",
        "wwtoall","wztoall","zztoall",
        "ttbarbg","ttbarsignalplustau",
        
        // ttbar signal variations
        "ttbarbg_massdown","ttbarbg_massup",
        "ttbarsignalplustau_massdown","ttbarsignalplustau_massup",
        "ttbarbg_matchingdown","ttbarbg_matchingup",
        "ttbarsignalplustau_matchingdown","ttbarsignalplustau_matchingup",
        "ttbarbg_scaledown","ttbarbg_scaleup",
        "ttbarsignalplustau_scaledown","ttbarsignalplustau_scaleup",
        
        // other ttbar generators
        "ttbarbg_mcatnlo", "ttbarsignalplustau_mcatnlo",
        "ttbarbg_powheg","ttbarsignalplustau_powheg",
        "ttbarbg_powhegHerwig","ttbarsignalplustau_powhegHerwig",
        "ttbarsignalplustau_Perugia11","ttbarsignalplustau_Perugia11NoCR",
        "ttbarbg_HadronicMadgraphWithSpinCorrelation","ttbarbg_SemiLeptMadgraphWithSpinCorrelation","ttbarsignalplustau_FullLeptMadgraphWithSpinCorrelation",
        
        // ttH samples
        "ttbarH110inclusive","ttbarH110tobbbar","ttbarH115inclusive","ttbarH115tobbbar","ttbarH120inclusive","ttbarH120tobbbar",
        "ttbarH125inclusive","ttbarH125tobbbar","ttbarH130inclusive","ttbarH130tobbbar","ttbarH135inclusive","ttbarH135tobbbar",
        "ttbarH1225tobbbar","ttbarH1275tobbbar","ttbarH140tobbbar",
        "ttbarW","ttbarZ","ttgjets",
        "ggH125zz4l", "ttww", "vbfH125zz4l", "wwgjets", "wwto2l2nu", "www", "wwz", "wzto3lnu", "zzto4l", "zzz"
    };

    std::vector<std::string> samples(samples_array, samples_array+sizeof(samples_array)/sizeof(samples_array[0]));
    std::vector<std::string> variable (variable_array, variable_array+sizeof(variable_array)/sizeof(variable_array[0]));
    
    TString outputName("ntupleValidation___");
    TString folderTmp1 = folder1;
    if(folderTmp1.EndsWith("/")) folderTmp1.Remove(folderTmp1.Last('/'));
    TString folderTmp2 = folder2;
    if(folderTmp2.EndsWith("/")) folderTmp2.Remove(folderTmp2.Last('/'));
    const TString short1 = folderTmp1.Data() + folderTmp1.Last('/') + 1;
    const TString short2 = folderTmp2.Data() + folderTmp2.Last('/') + 1;
    outputName.Append(short1).Append("___").Append(short2).Append(".txt");
    std::ofstream outputFile(outputName.Data());
    if(!outputFile.is_open()){
        std::cerr<<"ERROR! cannot open output file "<<outputName<<"\n...break\n"<<std::endl;
        exit(1);
    }
    
    outputFile<<"\n\nComparing trees from: "<<std::endl;
    outputFile<<"   ntupleA:   "<<folderTmp1<<std::endl;
    outputFile<<"   ntupleB:   "<<folderTmp2<<std::endl;
    outputFile<<"\n\n"<<std::endl;
    
    
    
    outputFile<<"\n--------------------------------------------------------------------------------------"<<std::endl;
    outputFile<<std::setw(60)<<"Sample "<<std::setw(15)<<"ntupleA"<<" "<<std::setw(15)<<"ntupleB"<<std::setw(15)<<" diff(%)"<<"  variable: unweightedEvents"<<std::endl;
    for (size_t iter = 0 ; iter<samples.size(); iter++){
        TString file0 = folderTmp1+"/"+samples.at(iter)+".root";
        TString file1 = folderTmp2+"/"+samples.at(iter)+".root";

        if(!checkFile(file0, outputFile)) continue;
        if(!checkFile(file1, outputFile)) continue;

        TFile *f0 = new TFile(file0);
        TFile *f1 = new TFile(file1);

        TH1I *t0 = (TH1I*)f0->Get("EventsBeforeSelection/unweightedEvents");
        TH1I *t1 = (TH1I*)f1->Get("EventsBeforeSelection/unweightedEvents");

        int nrEvents_old = t0->GetEntries();
        int nrEvents_new = t1->GetEntries();

        outputFile<<std::setw(40)<<samples.at(iter).c_str()<<std::setw(15)<<nrEvents_old<<std::setw(15)<<nrEvents_new<<std::setw(15)<<std::fixed<<std::setprecision(5)<<100*(1.*nrEvents_old-nrEvents_new)/nrEvents_old<<std::endl;
        delete t0; delete t1;
        f0->Close(); f1->Close();
        delete f0;   delete f1;
    }
    
    for (size_t vars=0; vars<variable.size(); vars++){
        outputFile<<"\n--------------------------------------------------------------------------------------"<<std::endl;
        outputFile<<std::setw(60)<<"Sample "<<std::setw(15)<<"ntupleA"<<" "<<std::setw(15)<<"ntupleB"<<std::setw(15)<<" diff(%)"<<"  variable: "<<variable.at(vars)<<std::endl;
        for (size_t iter = 0 ; iter<samples.size(); iter++){
            TString file0 = folderTmp1+"/"+samples.at(iter)+".root";
            TString file1 = folderTmp2+"/"+samples.at(iter)+".root";

            if(!checkFile(file0, outputFile)) continue;
            if(!checkFile(file1, outputFile)) continue;

            TFile *f0 = new TFile(file0);
            TFile *f1 = new TFile(file1);

            TTree *t0 = (TTree*)f0->Get("writeNTuple/NTuple");
            TTree *t1 = (TTree*)f1->Get("writeNTuple/NTuple");

            int nrEvents_old = nrEntries(t0, variable.at(vars));
            int nrEvents_new = nrEntries(t1, variable.at(vars));
            
            outputFile<<std::setw(40)<<samples.at(iter).c_str()<<std::setw(15)<<nrEvents_old<<std::setw(15)<<nrEvents_new<<std::setw(15)<<std::fixed<<std::setprecision(5)<<100.*(1.*nrEvents_old-nrEvents_new)/nrEvents_old<<std::endl;

            delete t0; delete t1;
            f0->Close(); f1->Close();
            delete f0;   delete f1;
        }
    }
}
