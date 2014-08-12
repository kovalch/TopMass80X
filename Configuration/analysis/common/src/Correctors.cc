#include <string>

#include "Correctors.h"
#include "utils.h"
#include "classes.h"
#include "TopAnalysis/ZTopUtils/interface/RecoilCorrector.h"




MetRecoilCorrector::MetRecoilCorrector(const std::string& fileData, const std::string& fileMC):
recoilCorrector_(new ztop::RecoilCorrector())
{
    std::cout<<"--- Beginning preparation of MVA MET recoil correction\n";
    
    std::string inputDir(common::DATA_PATH_COMMON());
    inputDir.append("/");
    const std::string inputData = inputDir + fileData;
    const std::string inputMC = inputDir + fileMC;
    
    std::cout<<"Using file for data: "<<inputData<<"\n";
    std::cout<<"Using file for MC: "<<inputMC<<"\n";
    
    recoilCorrector_->setFiles(inputData, inputMC);
    
    std::cout<<"=== Finishing preparation of MVA MET recoil correction\n\n";
}



MetRecoilCorrector::~MetRecoilCorrector()
{
    delete recoilCorrector_;
}



double MetRecoilCorrector::applyCorrection(LV& met, const LV& genZ, const LV& dilepton, const int nJet)const
{
    float metX = met.px();
    float metY = met.py();
    const double weight = recoilCorrector_->correctMet(metX, metY, genZ.px(), genZ.py(), dilepton.px(), dilepton.py(), nJet);
    met.SetPxPyPzE(metX, metY, met.pz(), met.E());
    return weight;
}







