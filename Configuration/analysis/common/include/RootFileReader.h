#ifndef RootFileReader_h
#define RootFileReader_h

#include <map>
#include <vector>
#include <iostream>
#include <type_traits>
#include <cstdlib>

#include <TObject.h>
#include <TString.h>

class TFile;





class RootFileReader{
    
public:
    
    /// Empty Constructor
    RootFileReader(){};
    
    /// Destructor
    ~RootFileReader();
    
    
    
    /// Returns the singleton instance
    static RootFileReader* getInstance();
    
    /// Find in a given file histograms whose names contain the given fragment,
    /// either anywhere in the name, or only at the begin of the name
    std::vector<TString> findObjects(const TString& filename, const TString& objectNameFragment, const bool fragmentAtBegin =true);
    
    /// Get a histogram from the file
    template <typename T> void Get(const TString& filename, const TString& histoname, T& result, const bool allowNonexisting =false, const bool verbosity =true);

    /// Get a histogram from the file, you need to pass the type here
    /// Warning: histo will be deleted once 60 more files have been opened
    template <typename T> const T* Get(const TString& filename, const TString& histoname, const bool allowNonexisting =false, const bool verbosity =true);
    
    /// Get a clone of a histogram from the file, you need to pass the type here
    template <typename T> T* GetClone(const TString& filename, const TString& histoname, const bool allowNonexisting =false, const bool verbosity =true);
    
    
    
private:
    
    /// Get a TObject from the file
    TObject* GetObj(const TString& filename, const TString& histoname, const bool allowNonexisting);
    
    
    
    /// Input files mapped to the file name
    std::map<TString, TFile*> fileMap_;
    
    /// Ordered file names
    std::vector<TString> fileOrder_;
    
    /// Accesses counter mapped to file name
    std::map<TString, int> accessed_;
    
    /// Openings counter mapped to file name
    std::map<TString, int> opened_;
};



template <typename T> void RootFileReader::Get(const TString& filename, const TString& histoname, T& result, const bool allowNonexisting, const bool verbosity)
{
    //is_assignable seems to be missing in this gcc version
    //static_assert(std::is_assignable<TObject*, T>::value == true, "You must convert to a TObject* like type!");
    static_assert(std::is_pointer<T>::value == true, "You must convert to a TObject* like type!");
    TObject* obj = this->GetObj(filename, histoname, allowNonexisting);
    if(obj == nullptr){
        if(allowNonexisting){
            if(verbosity) std::cout << "Warning: " << histoname << " is not in " << filename << std::endl;
            result = nullptr;
            return;
        }
        else{
            std::cerr << "ERROR: " << histoname << " is not in " << filename << std::endl;
            exit(1);
        }
    }
    result = dynamic_cast<T>(obj);
    if(!result){
        std::cerr << "The histogram " << histoname << " in file " << filename 
                  << " is of incompatible type (cannot typecast)!" << std::endl;
        exit(1);
    }
}



template <typename T> const T* RootFileReader::Get(const TString& filename, const TString& histoname, const bool allowNonexisting, const bool verbosity)
{
    T* result;
    this->Get(filename, histoname, result, allowNonexisting, verbosity);
    return result;
}



template <typename T> T* RootFileReader::GetClone(const TString& filename, const TString& histoname, const bool allowNonexisting, const bool verbosity)
{
    T* result;
    this->Get(filename, histoname, result, allowNonexisting, verbosity);
    if(!result && allowNonexisting) return nullptr;
    result = static_cast<T*>(result->Clone());
    result->SetDirectory(0);
    return result;
}




#endif




