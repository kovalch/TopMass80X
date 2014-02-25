#include <string>
#include <algorithm>
#include <fstream>

#include <TString.h>
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TObject.h>

#include "RootFileReader.h"





RootFileReader::~RootFileReader()
{
//     std::cout << "Deleting RootFileReader\n";
//     for(const auto& i : fileMap_){
//         delete i.second;
//     }
//     for(const auto& i : accessed_){
//         std::cout << "accessed: " << i.first << " " << i.second << std::endl;
//     }
//     for(const auto& i : opened_){
//         std::cout << "opened: " << i.first << " " << i.second << std::endl;
//     }
}



TObject* RootFileReader::GetObj(const TString& filename, const TString& histoname, const bool allowNonexisting)
{
    ifstream inputFileStream(filename);
    if(!inputFileStream.is_open()){
        if(allowNonexisting) return nullptr;
        std::cerr << "\n\n******************* ERROR ******************* ERROR ******************* ERROR *******************\n\n"
                    << "The file " << filename << " does not exist, thus cannot get histogram " << histoname
                    << "\n\n******************* ERROR ******************* ERROR ******************* ERROR *******************\n"
                    << std::endl;
        exit(1);
    }
    inputFileStream.close();

    auto& file = fileMap_[filename];
    if(!file){
        file = TFile::Open(filename);
        ++opened_[filename];
    }
    else{
        fileOrder_.erase(std::remove_if(begin(fileOrder_), end(fileOrder_), 
                         [=](const TString &str){ return str == filename;}),
                         end(fileOrder_));
    }
    
    fileOrder_.push_back(filename);
    ++accessed_[filename];
    //at max, keep 20 files open
    if(fileMap_.size() > 20){
        //delete 0th element
        delete fileMap_[fileOrder_[0]];
        fileMap_.erase(fileOrder_[0]);
        fileOrder_.erase(fileOrder_.begin());
    }

    TObject* result = file->Get(histoname);
    TH1* as_histo = dynamic_cast<TH1*>(result);
    if(as_histo) as_histo->SetDirectory(file);
    return result;
}



std::vector<TString> RootFileReader::findObjects(const TString& filename, const TString& objectNameFragment, const bool fragmentAtBegin)
{
    ifstream inputFileStream(filename);
    if(!inputFileStream.is_open()){
        const std::string helpText(fragmentAtBegin ? "beginning with: " : "containing: ");
        std::cerr << "\n\n******************* ERROR ******************* ERROR ******************* ERROR *******************\n\n"
                  << "The file "<<filename<<" does not exist, thus cannot search for object "<<helpText<<objectNameFragment
                  << "\n\n******************* ERROR ******************* ERROR ******************* ERROR *******************\n"
                  << std::endl;
        exit(2);
    }
    inputFileStream.close();
    
    const TFile* file = TFile::Open(filename);
    std::vector<TString> result;
    const TList* objectList = file->GetListOfKeys();
    for(int i = 0; i < objectList->GetSize(); ++i){
        TObject* i_objectList = objectList->At(i);
        TString* objectName(0);
        objectName = new TString(i_objectList->GetName());
        //std::cout<<"\n\tName of object: "<<*objectName<<"\n\n";
        if(fragmentAtBegin){
            if(objectName->BeginsWith(objectNameFragment)) result.push_back(*objectName);
        }
        else{
            if(objectName->Contains(objectNameFragment)) result.push_back(*objectName);
        }
    }
    return result;
}



RootFileReader* RootFileReader::getInstance()
{
    static RootFileReader instance;
    return &instance;
}






