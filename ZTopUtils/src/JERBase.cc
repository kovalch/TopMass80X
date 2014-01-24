#include "../interface/JERBase.h"

namespace ztop {

void JERBase::setSystematics(std::string type) {
    resranges_.clear();
    resranges_ << 0 << 0.5 << 1.1 << 1.7 << 2.3 << 5;

    if (type == "def") { //standard
        resfactors_.clear();
        resfactors_ << 1.052 << 1.057 << 1.096 << 1.134 << 1.288;
        std::cout << "JER set to default" << std::endl;
    } else if (type == "down") {
        resfactors_.clear();
        resfactors_ << 1.115 << 1.114 << 1.162 << 1.228 << 1.488;
        std::cout << "JER set to syst down" << std::endl;
    } else if (type == "up") {
        resfactors_.clear();
        resfactors_ << 0.990 << 1.001 << 1.032 << 1.042 << 1.089;
        std::cout << "JER set to syst up" << std::endl;
    }

}

void JERBase::correctP4(float & recopt, float& recoeta, float & recophi, float & recom, //full lorentzvector
        const float & genpt) const{
    if (genpt < 1)
        return;


    std::vector<float>::const_iterator it=std::lower_bound(resranges_.begin(),
            resranges_.end(), recoeta);
    size_t etabin=0;
    if(recoeta==*it)
        etabin= it-resranges_.begin();
    else
        etabin= it-resranges_.begin()-1;

    double deltapt = (1. - resfactors_.at(etabin))* (recopt - genpt);
    double scale = std::max(0., recopt + deltapt) / recopt;
    recopt*=scale;
    recom*=scale;
}
} //namespace
