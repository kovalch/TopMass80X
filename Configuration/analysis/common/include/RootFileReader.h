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
//class TString;





class RootFileReader{
    
public:
    
    /// Empty Constructor
    RootFileReader(){};
    
    /// Destructor
    ~RootFileReader();
    
    
    
    //returns the singleton instance
    static RootFileReader* getInstance();
    
    // Find in a given file histograms whose names contain the given fragment,
    // either anywhere in the name, or only at the begin of the name
    std::vector<TString> findObjects(const char* filename, const char* objectNameFragment, const bool fragmentAtBegin =true);
    
    //get a histogram from the file
    template <typename T> void Get(const char* filename, const char* histoname, T& result, const bool allowNonexisting =false, const bool verbosity =true);

    //get a histogram from the file, you need to pass the type here
    //warning: histo will be deleted once 60 more files have been opened
    template <typename T> const T* Get(const char* filename, const char* histoname, const bool allowNonexisting =false, const bool verbosity =true);
    
    template <typename T> T* GetClone(const char* filename, const char* histoname, const bool allowNonexisting =false, const bool verbosity =true);
    
    
    
private:
    
    TObject* GetObj(const char* filename, const char* histoname, const bool allowNonexisting);
    
    std::map<TString, TFile*> fileMap_;
    std::vector<TString> fileOrder_;
    std::map<TString, int> accessed_, opened_;
};



template <typename T> void RootFileReader::Get(const char* filename, const char* histoname, T& result, const bool allowNonexisting, const bool verbosity)
{
    //is_assignable seems to be missing in this gcc version
    //static_assert(std::is_assignable<TObject*, T>::value == true, "You must convert to a TObject* like type!");
    static_assert(std::is_pointer<T>::value == true, "You must convert to a TObject* like type!");
    TObject *obj = GetObj(filename, histoname, allowNonexisting);
    if (obj == nullptr) {
        if (allowNonexisting) {
            if(verbosity) std::cout << "Warning: " << histoname << " is not in " << filename << std::endl;
            result = nullptr;
            return;
        } else {
            std::cerr << "ERROR: " << histoname << " is not in " << filename << std::endl;
            exit(1);
        }            
    }
    result = dynamic_cast<T>(obj);
    if (!result) {
        std::cerr << "The histogram " << histoname << " in file " << filename 
            << " is of incompatible type (cannot typecast)!" << std::endl;
        exit(1);
    }
}



template <typename T> const T* RootFileReader::Get(const char* filename, const char* histoname, const bool allowNonexisting, const bool verbosity)
{
    T* result;
    Get(filename, histoname, result, allowNonexisting, verbosity);
    return result;
}



template <typename T> T* RootFileReader::GetClone(const char* filename, const char* histoname, const bool allowNonexisting, const bool verbosity)
{
    T* result;
    Get(filename, histoname, result, allowNonexisting, verbosity);
    if (!result && allowNonexisting) return nullptr;
    result = static_cast<T*>(result->Clone());
    result->SetDirectory(0);
    return result;
}





#endif // RootFileReader_h




