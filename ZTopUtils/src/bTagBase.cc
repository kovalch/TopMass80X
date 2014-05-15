#include "../interface/bTagBase.h"
#include <stdexcept>
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"

namespace ztop {
/**
 * constructor does nothing in particular
 */
bTagBase::bTagBase():showWarnings(false),debug(false),init_(false),wp_(length_wp),is2011_(false),
        syst_(nominal),makeeffs_(true),tempsamplename_(""),histp_(0),effhistp_(0),medianvecp_(0),sumStuffEff_(0.9999999),sumStuffSfEff_(0.9999999),
        nancount_(0),hf_rewfile_(0),lf_rewfile_(0){
    wpvals_.resize(workingPoints::length_wp, 0);
    minpt_.resize (workingPoints::length_wp, 0);
    maxpt_.resize (workingPoints::length_wp, 0);
    //load all variables
    initWorkingpoints();

    //setting up naming scheme:
    histoNames_.   resize(histoNames_en::length_histoNames_en, "");
    effHistoNames_.resize(effHistoNames_en::length_effHistoNames_en, "");

    // Associating histogram names
    histoNames_.at(bjets_h)="bjets2D";
    histoNames_.at(btagged_h)="bjetsTagged2D";
    histoNames_.at(cjets_h)="cjets2D";
    histoNames_.at(ctagged_h)="cjetsTagged2D";
    histoNames_.at(ljets_h)="ljets2D";
    histoNames_.at(ltagged_h)="ljetsTagged2D";
    // Adding efficiency histogram names
    effHistoNames_.at(beff_h)="beff2D";
    effHistoNames_.at(ceff_h)="ceff2D";
    effHistoNames_.at(leff_h)="leff2D";

    //the names are once defined here, NOWHERE ELSE
    //they don't need ANY additional ordering AT ALL.

    //for shape reweighting
    ///REWRITE!!!

    setupShapeReweightingBins();


}
const float& bTagBase::getWPDiscrValue()const{
    if(wp_==length_wp){ //default value
        throw std::logic_error("bTagBase::getWPDiscrValue: working point not set");
    }
    return wpvals_[wp_];
}

bTagBase::workingPoints bTagBase::getWorkingPoint() const {
    if(wp_==length_wp){ //default value
        throw std::logic_error("bTagBase::getWorkingPoint: working point not set");
    }
    return wp_;
}

bTagBase::bTagBase(const bTagBase& rhs){
    copyFrom(rhs);
}
bTagBase & bTagBase::operator = (const bTagBase& rhs){
    copyFrom(rhs);
    return *this;
}
void bTagBase::copyFrom(const bTagBase& rhs){
    if(this == &rhs) return;
    showWarnings=rhs.showWarnings;
    debug=rhs.debug;

    histos_=rhs.histos_;
    effhistos_=rhs.effhistos_;
    medianMap_=rhs.medianMap_;

    init_=rhs.init_;
    wp_=rhs.wp_;
    is2011_=rhs.is2011_;
    syst_=rhs.syst_;
    wpvals_=rhs.wpvals_;
    minpt_=rhs.minpt_;
    maxpt_=rhs.maxpt_;
    histoNames_=rhs.histoNames_;
    effHistoNames_=rhs.effHistoNames_;
    makeeffs_=rhs.makeeffs_;
    tempsamplename_=rhs.tempsamplename_;

    if(rhs.histp_){
        for(std::map<std::string, std::vector<TH2D> >::const_iterator it=rhs.histos_.begin();it!=rhs.histos_.end();++it){
            if(rhs.histp_ == &it->second){
                histp_=&histos_[it->first];
            }
        }
    }
    else{
        histp_=0;
    }
    if(rhs.effhistp_){
        for(std::map<std::string, std::vector<TH2D> >::const_iterator it=rhs.effhistos_.begin();it!=rhs.effhistos_.end();++it){
            if(rhs.effhistp_ == &it->second){
                effhistp_=&effhistos_[it->first];
            }
        }
    }
    else{
        effhistp_=0;
    }
    if(rhs.medianvecp_){
        for(std::map<std::string, std::vector<float> >::const_iterator it=rhs.medianMap_.begin();it!=rhs.medianMap_.end();++it){
            if(rhs.medianvecp_ == &it->second){
                medianvecp_=&medianMap_[it->first];
            }
        }
    }
    else{
        medianvecp_=0;
    }

    sumStuffEff_=rhs.sumStuffEff_;
    sumStuffSfEff_=rhs.sumStuffSfEff_;
    nancount_=rhs.nancount_;

    shapeRWPtBinsHF_=rhs.shapeRWPtBinsHF_;
    shapeRWPtBinsLF_=rhs.shapeRWPtBinsLF_;
    shapeRWEtaBinsLF_=rhs.shapeRWEtaBinsLF_;

    hf_rewfile_=0;
    lf_rewfile_=0;
    if(rhs.hf_rewfile_)
        hf_rewfile_=(TFile*)rhs.hf_rewfile_->Clone();
    if(rhs.lf_rewfile_)
        lf_rewfile_=(TFile*)rhs.lf_rewfile_->Clone();

    h_csv_wgt_hf_=rhs.h_csv_wgt_hf_;
    c_csv_wgt_hf_=rhs.c_csv_wgt_hf_;
    h_csv_wgt_lf_=rhs.h_csv_wgt_lf_;

}


bTagBase::~bTagBase() {
    cleanptr();

    if(hf_rewfile_) {hf_rewfile_->Close();delete hf_rewfile_;}
    if(lf_rewfile_) {lf_rewfile_->Close();delete lf_rewfile_;}
}


void bTagBase::setSystematic(systematics sys) {
    syst_ = sys;
    if(hf_rewfile_ && lf_rewfile_){ //reweight files are associated

        TString histsuff=assoHFSysHistoName(syst_);
        for(size_t i=0;i<h_csv_wgt_hf_.size();i++){
            h_csv_wgt_hf_.at(i) = tryToGet(hf_rewfile_, Form("csv_ratio_Pt%i_Eta0_%s",(int)i,histsuff.Data()));
        }
        histsuff=assoCSysHistoName(syst_);
        for(size_t i=0;i<c_csv_wgt_hf_.size();i++){
            c_csv_wgt_hf_.at(i) = tryToGet(hf_rewfile_,Form("c_csv_ratio_Pt%i_Eta0_%s",(int)i,histsuff.Data()));
        }
        histsuff=assoLFSysHistoName(syst_);
        for(size_t i=0;i<h_csv_wgt_lf_.size();i++){
            for(size_t j=0;j<h_csv_wgt_lf_.at(i).size();j++){
                h_csv_wgt_lf_.at(i).at(j) = tryToGet(lf_rewfile_, Form("csv_ratio_Pt%i_Eta%i_%s",(int)i,(int)j,histsuff.Data()));
            }
        }


    }

}

/**
 * sets the sample name. It can be used as a unique identifier when
 * reading or writing the histograms/data
 * If input is read, the name is used to identify the right set of
 * histograms (e.g. different for each sample but in the same bTagBase
 * object)
 */
int bTagBase::setSampleName(const std::string & samplename) {
    //set pointers
    tempsamplename_=samplename;
    initWorkingpoints();
    if(debug)
        std::cout << "bTagBase::setSampleName " <<  tempsamplename_<< std::endl;


    std::map<std::string, std::vector<TH2D> >::iterator sampleit = histos_.find(
            tempsamplename_);
    if (sampleit != histos_.end())
        histp_ = &(sampleit->second);
    else
        histp_ = 0;
    std::map<std::string, std::vector<TH2D> >::iterator effit = effhistos_.find(
            tempsamplename_);
    if (effit != effhistos_.end())
        effhistp_ = &(effit->second);
    else
        effhistp_ = 0;

    if (!makeeffs_ && effhistp_ && histp_) {
        std::cout << "loaded b-tag efficiency histos for " << tempsamplename_
                << std::endl;
        std::map<std::string, std::vector<float> >::iterator medianit = medianMap_.find(
                tempsamplename_);
        if ( medianit == medianMap_.end() ) {
            std::cout<<" bTagBase::setSampleName: Median map doesn't exist for "<<tempsamplename_
                    <<". exit!"<<std::endl;
            return -3;
        }
        medianvecp_ = &(medianit->second);
        return 1;
    } else if (!makeeffs_ && !effhistp_) {
        std::cout << " bTagBase::setSampleName: efficiency for " << tempsamplename_
                << " not derived, yet! exit." << std::endl;
        return -1;
    }

    // efficiencies not determined yet -> prepare efficiency histos
    // bool adddirtemp=TH1::AddDirectory;
    TH1::AddDirectory(kFALSE); //prevent some weird root behaviour
    std::cout << "preparing b-tag efficiency histos for " << tempsamplename_
            << std::endl;

    //////// put in algorithm to automatically rebin in case of low statistics !! then use a fine (equals btv binning) binning

    float effptbins[] = { 20., 50., 70., 100., 160., 210., 800. };
    unsigned int npt = 7;
    float effetabins[] = { 0.0, 0.5, 1.0, 2.5 };
    unsigned int neta = 4;

    //probably there will be less statistics for the light jets, so histos might need a coarser binning

    float l_effptbins[] = { 20., 70., 120., 800. };
    unsigned int l_npt = 4;
    float l_effetabins[] = { 0.0, 1.5, 3.0 };
    unsigned int l_neta = 3;

    std::vector<TH2D> temp;

    // Defining the histograms. Should have names as in histoNames_ and effHistoNames_
    TH2D bjets = TH2D(getJetHistoOrderedNames().at(bjets_h).c_str(), "unTagged Bjets", npt - 1,
            effptbins, neta - 1, effetabins);
    bjets.Sumw2();

    TH2D bjetstagged = TH2D(getJetHistoOrderedNames().at(btagged_h).c_str(), "Tagged Bjets",
            npt - 1, effptbins, neta - 1, effetabins);
    bjetstagged.Sumw2();

    TH2D cjets = TH2D(getJetHistoOrderedNames().at(cjets_h).c_str(), "unTagged Cjets", npt - 1,
            effptbins, neta - 1, effetabins);
    cjets.Sumw2();

    TH2D cjetstagged = TH2D(getJetHistoOrderedNames().at(ctagged_h).c_str(), "Tagged Cjets",
            npt - 1, effptbins, neta - 1, effetabins);
    cjetstagged.Sumw2();

    TH2D ljets = TH2D(getJetHistoOrderedNames().at(ljets_h).c_str(), "unTagged Ljets", l_npt - 1,
            l_effptbins, l_neta - 1, l_effetabins);
    ljets.Sumw2();

    TH2D ljetstagged = TH2D(getJetHistoOrderedNames().at(ltagged_h).c_str(), "Tagged Ljets",
            l_npt - 1, l_effptbins, l_neta - 1, l_effetabins);
    ljetstagged.Sumw2();

    temp << bjets << bjetstagged << cjets << cjetstagged << ljets
            << ljetstagged;
    histos_[tempsamplename_] = temp;
    histp_ = &(histos_.find(tempsamplename_)->second);


    TH2D beff = TH2D(getEffHistoOrderedNames().at(beff_h).c_str(), "Bjets eff", npt - 1, effptbins,
            neta - 1, effetabins);
    beff.Sumw2();

    TH2D ceff = TH2D(getEffHistoOrderedNames().at(ceff_h).c_str(), "Cjets eff", npt - 1, effptbins,
            neta - 1, effetabins);
    ceff.Sumw2();

    TH2D leff = TH2D(getEffHistoOrderedNames().at(leff_h).c_str(), "Ljets eff", l_npt - 1, l_effptbins,
            l_neta - 1, l_effetabins);
    leff.Sumw2();

    temp.clear();
    temp << beff << ceff << leff;
    effhistos_[tempsamplename_] = temp;
    effhistp_ = &(effhistos_.find(tempsamplename_)->second);

    return 0;

}

/**
 * adds an entry for a jet with p4, genpartonFlavour, bDiscrValue  and PUweight
 * to the efficiency histograms
 */
void bTagBase::fillEff(const float& pt, const float&abs_eta,
        const int &genPartonFlavour, const float &bDiscrVal,
        const float &puweight) {

    if(debug) std::cout << "bTagBase::fillEff " << std::endl;
    if (!makeeffs_)
        return;
    if(wp_==length_wp){ //default value
        throw std::logic_error("bTagBase::fillEff: working point not set");
    }
    if (!histp_) { //protection
        std::cout
        << "bTagBase::fillEff: you have to set a samplename before filling efficiency histograms!"
        << std::endl;
        throw std::runtime_error(
                "bTagBase::fillEff: you have to set a samplename before filling efficiency histograms!");
    }
    if (genPartonFlavour == 0)
        return;

    jetTypes jettype=jetType(genPartonFlavour);

    if (jettype == bjet) { // b jets
        histp_->at(0).Fill(pt, abs_eta, puweight);
        if (bDiscrVal > wpvals_[wp_])
            histp_->at(1).Fill(pt, abs_eta, puweight);
    } else if (jettype == cjet) { // c jets

        histp_->at(2).Fill(pt, abs_eta, puweight);
        if (bDiscrVal > wpvals_[wp_])
            histp_->at(3).Fill(pt, abs_eta, puweight);
    } else if (jettype == lightjet) { // light jets (including gluon jets)
        histp_->at(4).Fill(pt, abs_eta, puweight);
        if (bDiscrVal > wpvals_[wp_])
            histp_->at(5).Fill(pt, abs_eta, puweight);
    }
}
/**
 * creates efficiency histograms, to be run after all are filled
 */
void bTagBase::makeEffs() {
    if (!makeeffs_)
        return;
    if(debug) std::cout << "bTagBase::makeEffs " << std::endl;
    if (!histp_ || !effhistp_) {
        std::cout << "bTagBase::makeEffs: you have to set a samplename!"
                << std::endl;
        throw std::runtime_error(
                "bTagBase::makeEffs: you have to set a samplename!");
    }
    if(wp_==length_wp){ //default value
        throw std::logic_error("bTagBase::makeEffs working point not set");
    }

    for (unsigned int i = 0; i < effhistp_->size(); i++) {
        for (int binx = 1; binx <= effhistp_->at(i).GetNbinsX() + 1; binx++) {
            for (int biny = 1; biny <= effhistp_->at(i).GetNbinsY() + 1;
                    biny++) {
                //
                //uses histos at 2i and 2i+1
                float cont = 1;
                float err = 0.99; //to avoid zeros!
                if (histp_->at(2 * i).GetBinContent(binx, biny) > 0) {
                    cont = histp_->at(2 * i + 1).GetBinContent(binx, biny)
                        						                                                                                                                                                                         / histp_->at(2 * i).GetBinContent(binx, biny);
                    if (debug)
                        std::cout << "makeEffs: content: " << cont;
                    err = sqrt(
                            cont * (1 - cont)
                            / histp_->at(2 * i).GetBinContent(binx, biny));
                    if (debug)
                        std::cout << " error: " << err << "  " << binx << " "
                        << biny << std::endl;
                }
                if (err > 0.02 && err != 0.99) {
                    if (showWarnings)
                        std::cout
                        << "bTagBase::makeEffs: warning. error in bin ("
                        << binx << ", " << biny << ") for histogram "
                        << effhistp_->at(i).GetName()
                        << " is larger than 0.02" << std::endl;

                    //here there is space for automatic rebinning!! just bool binok, if not merge 2 bins - maybe in second step, here performance is not crucial

                }
                effhistp_->at(i).SetBinContent(binx, biny, cont);
                effhistp_->at(i).SetBinError(binx, biny, err);
            }
        }
    }
    std::vector<float> medianVec(length_median, 0);
    medianVec.at(bpt) = ztop::bTagBase::median(histp_->at(1).ProjectionX());
    medianVec.at(beta) = ztop::bTagBase::median(histp_->at(1).ProjectionY());
    medianVec.at(cpt) = ztop::bTagBase::median(histp_->at(3).ProjectionX());
    medianVec.at(ceta) = ztop::bTagBase::median(histp_->at(3).ProjectionY());
    medianVec.at(lpt) = ztop::bTagBase::median(histp_->at(5).ProjectionX());
    medianVec.at(leta) = ztop::bTagBase::median(histp_->at(5).ProjectionY());
    medianMap_[tempsamplename_] = medianVec;
}

void bTagBase::readShapeReweightingFiles(const TString& heavy,const TString& light){
    TFile * f_CSVwgt_HF = new TFile (heavy);
    if(!f_CSVwgt_HF || f_CSVwgt_HF->IsZombie()){
        throw std::runtime_error("bTagBase::readShapeReweightingFiles: heavy flavour file not ok");
    }
    TFile * f_CSVwgt_LF = new TFile (light);
    if(!f_CSVwgt_LF || f_CSVwgt_LF->IsZombie()){
        throw std::runtime_error("bTagBase::readShapeReweightingFiles: light flavour file not ok");
    }
    if(hf_rewfile_) {hf_rewfile_->Close();delete hf_rewfile_;}
    if(lf_rewfile_) {lf_rewfile_->Close();delete lf_rewfile_;}

    hf_rewfile_=f_CSVwgt_HF;
    lf_rewfile_=f_CSVwgt_LF;


    setSystematic(syst_);


}










/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// private members
///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


void bTagBase::countJet(const float & pt, const float& abs_eta,
        const int & genPartonFlavor) {

    if(makeeffs_) return;
    if(debug) std::cout << "bTagBase::countJet " << std::endl;
    if(wp_==length_wp){ //default value
        throw std::logic_error("bTagBase::countJet: working point not set");
    }
    if (genPartonFlavor == 0)
        return;

    jetTypes jettype=jetType(genPartonFlavor);
    float sf = 1;
    float eff = 0;
    sf=jetSF(pt, abs_eta,jettype);
    eff=jetEff(pt, abs_eta,jettype);


    sumStuffEff_ = sumStuffEff_ * (1 - eff);
    sumStuffSfEff_ = sumStuffSfEff_ * (1 - sf * eff);

}

bool bTagBase::jetIsTagged(const float & pt, const float& abs_eta,
        const int & genPartonFlavor, const float & tagValue, const unsigned int & seed) const {


    if(wp_==length_wp){ //default value
        throw std::logic_error("bTagBase::jetIsTagged: working point not set");
    }
    jetTypes jettype=jetType(genPartonFlavor);
    const bool isBTagged = tagValue > wpvals_[wp_];
    if(makeeffs_) return isBTagged;
    if(debug) std::cout << "bTagBase::jetIsTagged " << std::endl;

    const double Btag_eff = jetEff(pt, abs_eta, jettype);
    const double Btag_SF = jetSF(pt, abs_eta, jettype);
    const float coin = TRandom3(seed).Uniform(1.0);

    if ( Btag_SF > 1. ) {  // use this if SF>1
        if ( !isBTagged ) {
            // fraction of jets that need to be upgraded
            const float mistagPercent = (1.0 - Btag_SF) / (1.0 - (1.0/Btag_eff) );

            //upgrade to tagged
            if( coin < mistagPercent ) return true;
        }
    }
    else if ( Btag_SF < 1. ) {  // use this if SF<1
        // downgrade tagged to untagged
        if ( isBTagged && coin > Btag_SF ) return false;
    }
    else {  // no change if exactly SF==1
        return isBTagged;
    }
    return isBTagged;

}
float bTagBase::getEventSF() {
    if(makeeffs_) return 1.;
    float sf = (1 - sumStuffSfEff_) / (1 - sumStuffEff_);
    if(sf!=sf){
        sf=1;
        ++nancount_;
    }
    resetCounter();
    return sf;
}





float bTagBase::getJetDiscrShapeWeight(const float & pt, const float& abs_eta,
        const int & genPartonFlavor, const float& jetdiscr)const{

    if(genPartonFlavor==0) //covers data
        return 1;

    if(!hf_rewfile_ || !lf_rewfile_){
        throw std::logic_error("bTagBase::getJetDiscrShapeWeight: first read in weight files!");
    }

    const jetTypes & type = jetType(genPartonFlavor);

    const TH1D * histo=0;

    if(type==bjet)
        histo=getShapeHFHisto(pt, true);
    else if(type==cjet)
        histo=getShapeHFHisto(pt,false);
    else
        histo=getShapeLFHisto(pt,abs_eta);

    int useCSVBin = (jetdiscr>=0.) ? histo->FindFixBin(jetdiscr) : 1;
    float out=histo->GetBinContent(useCSVBin);
    if(out<=0) out=1;
    return out;

}





float bTagBase::jetSF(const float &pt, const float& abs_eta,const jetTypes & jettype) const{
    if(jettype==bjet) return BJetSF(pt, abs_eta);
    if(jettype==cjet) return CJetSF(pt, abs_eta);
    if(jettype==lightjet) return LJetSF(pt, abs_eta);
    else return 1.;
}

float bTagBase::jetEff(const float &pt, const float& abs_eta,const jetTypes & jettype) const{
    if(jettype == undefined)
        throw std::logic_error("bTagBase::jetEff: undefined jet type");
    if(!effhistp_)
        throw std::logic_error("bTagBase::jetEff: histogram pointer not set (bTagBase::setSampleName not called)");
    unsigned int effh = 100;
    if(jettype==bjet)  effh = 0;
    else if(jettype==cjet)  effh = 1;
    else if(jettype==lightjet)  effh = 2;

    int ptbin = 0;
    int etabin = 0;
    int bla = 0;

    effhistp_->at(effh).GetBinXYZ(effhistp_->at(effh).FindBin(pt, abs_eta),
            ptbin, etabin, bla);
    return effhistp_->at(effh).GetBinContent(ptbin, etabin);
}

float bTagBase::median(TH1* h1) const{
    int nBin = h1->GetXaxis()->GetNbins();
    std::vector<double> x(nBin);
    h1->GetXaxis()->GetCenter(&x[0]);
    TH1D* h1D = dynamic_cast<TH1D*>(h1);
    if(!h1D){
        std::cerr << "Median needs a TH1D!\n";
        throw std::logic_error("bTagBase::median: Median needs a TH1D!\n" );
    }
    const double* y = h1D->GetArray();
    // exclude underflow/overflows from bin content array y
    return TMath::Median(nBin, &x[0], &y[1]);
}





const TH1D * bTagBase::getShapeLFHisto(const float & pt, const float & fabseta)const{

    size_t ptbin = std::lower_bound(shapeRWPtBinsLF_.begin(), shapeRWPtBinsLF_.end(), pt) - shapeRWPtBinsLF_.begin();
    if ((shapeRWPtBinsLF_.at(ptbin) == pt || ptbin == shapeRWPtBinsLF_.size()) && shapeRWPtBinsLF_.at(0) != pt) //allow lowest bin to bin on boundary!
        --ptbin;
    if(ptbin>0)--ptbin; //no overflow
    size_t fetabin = std::lower_bound(shapeRWEtaBinsLF_.begin(), shapeRWEtaBinsLF_.end(), fabseta) - shapeRWEtaBinsLF_.begin();
    if ((shapeRWEtaBinsLF_.at(fetabin) == fabseta || fetabin == shapeRWEtaBinsLF_.size()) && shapeRWEtaBinsLF_.at(0) != fabseta) //allow lowest bin to bin on boundary!)
        --fetabin;
    if(fetabin>0)--fetabin; //no overflow

    try{
        return &h_csv_wgt_lf_.at(ptbin).at(fetabin);
    }catch(...){
        std::cout << "bTagBase::getShapeLFHisto: pt: " << pt << " eta: " << fabseta  << " bin: "<< ptbin << "," <<fetabin  << std::endl;
        throw std::out_of_range("bTagBase::getShapeLFHisto: pt/eta out of range");
    }

    return 0;

}
const TH1D * bTagBase::getShapeHFHisto(const float & pt, bool bjet)const {

    size_t ptbin = std::lower_bound(shapeRWPtBinsHF_.begin(), shapeRWPtBinsHF_.end(), pt) - shapeRWPtBinsHF_.begin();
    if ((shapeRWPtBinsHF_.at(ptbin) == pt || ptbin == shapeRWPtBinsHF_.size()) && shapeRWPtBinsLF_.at(0) != pt)
        --ptbin;
    if(ptbin>0)--ptbin; //no overflow
    try{
        if(bjet)
            return &h_csv_wgt_hf_.at(ptbin);
        else
            return &c_csv_wgt_hf_.at(ptbin);
    }catch(...){
        std::cout << "bTagBase::getShapeHFHisto: pt: " << pt << " bin: "<< ptbin << std::endl;
        throw std::out_of_range("bTagBase::getShapeHFHisto: pt/eta out of range");
    }
    return 0;
}


TH1D bTagBase::tryToGet(TFile * f, TString name)const{
    TH1D  * h= (TH1D  *)f->Get(name);
    if(!h || h->IsZombie()){
        std::cout << "bTagBase::tryToGet " << name << " - not found in "<< f->GetName() <<std::endl;
        throw std::runtime_error("bTagBase::tryToGet: Histogram not found " );
    }
    TH1D out=*h;
    delete h;
    return out;
}

// ///////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// BTagPOG INPUT //////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////

/*
 * For each tagger or working point, the input from the BTag POG needs to be implemented
 * as follows:
 *
 * Add the tagger name including the working point to the enum workingPoints in bTagBase.h
 * You can find the definition right on top of the class declaration. Do not touch the last entry.
 *
 * 1. Add the working point discriminator value to the function initWorkingpoints()
 *    as wpvals_.at(<your_new_enum)= discriminatorvalue;
 *    and the pt range the SF are defined in
 *    minpt_.at(<your_new_enum)=?;
 *    maxpt_.at(<your_new_enum)=?;
 * 1a add the name of the working point as a strin to getWorkingPointString()
 *
 *
 * 2. Add the corresponding scale factors for B jets to the function BJetSF
 *    - copy the full entry starting with "if(wp_== csvl_wp){"
 *    - paste it right before the marker --END OF WORKING POINTS---
 *    - edit the comments where the SF come from, give a link etc
 *    - adapt numbers and WP definition, change "if" to "else if"
 *    - don't change the parts marked as "don't change"
 *      (--->   <region> <---- indicates a region not to be changed)
 * 3. Add the scale factors to the function CJetSF/LJetSF
 *    - similar to BJetSF
 *
 */

void bTagBase::initWorkingpoints() {
    if(debug)
        std::cout << "bTagBase::initWorkingpoints: "<< wpvals_.size() << std::endl;

    // Please give some information about the SF here. E.g
    // -------
    // SF from:
    // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt  (EPS 2013 prescription)
    //  CSVL
    // -------
    wpvals_.at(csvl_wp) = 0.244;
    minpt_.at(csvl_wp) = 20;
    maxpt_.at(csvl_wp) = 800;

    wpvals_.at(csvm_wp) = 0.679;
    minpt_.at(csvm_wp) = 20;
    maxpt_.at(csvm_wp) = 800;

    wpvals_.at(csvt_wp) = 0.898;
    minpt_.at(csvt_wp) = 20;
    maxpt_.at(csvt_wp) = 800;


}

std::string bTagBase::getWorkingPointString()const{
    if     (wp_ == csvl_wp) return "csvl";
    else if (wp_ == csvm_wp) return "csvm";
    else if (wp_ == csvt_wp) return "csvt";
    else return "notDef";
}

float bTagBase::BJetSF(const float & pt, const float& abs_eta,
        float multiplier) const { //multiplier not const for a reason!

    float x = pt;

    if (is2011_) { //for 7 TeV !!!!

    } else { //this is for 8 TeV!!!
        if (wp_ == csvl_wp) { // ############### CSV Loose ################

            // Please give some information about the SF here. E.g
            // -------
            // SF from:
            // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt  (EPS 2013 prescription)
            // -------

            // please use the UPPER bin boundaries array given in the BTagPOG payload and then
            // add the lowest bin (e.g. 20) at front
            // make sure to use "static"

            static float ptbins[] = { 20, 30, 40, 50, 60, 70, 80, 100, 120, 160,
                    210, 260, 320, 400, 500, 600, 800 };

            //put the errors to the corresponding SF here (copy/paste from BTagPOG payload)
            //make sure to use static and that it has the same size as ptbins!!!
            static float SFb_error[] = { 0.033408, 0.015446, 0.0146992, 0.0183964,
                    0.0185363, 0.0145547, 0.0176743, 0.0203609,
                    0.0143342, 0.0148771, 0.0157936, 0.0176496,
                    0.0209156, 0.0278529, 0.0346877, 0.0350101 };

            //don't change ---------------------------->
            static size_t ptbinsSize = sizeof(ptbins) / sizeof(ptbins[0]);
            static size_t SFb_errorSize = sizeof(SFb_error) / sizeof(SFb_error[0]);
            if (SFb_errorSize != ptbinsSize - 1) {
                std::cout
                << "bTagBase::BJetSF: Size of SFb_error should be one less than of ptbins. throwing exception!"
                << std::endl;
                throw std::logic_error(
                        "bTagBase::BJetSF: Size of SFb_error should be one less than of ptbins.");
            }
            //force range according to BTagPOG input (dont change)
            if (x < ptbins[0]) {
                x = ptbins[0];
                multiplier *= 2;
            } else if (x >= ptbins[ptbinsSize - 1]) {
                x = ptbins[ptbinsSize - 1] - 0.0001;
                multiplier *= 2;
            }
            // <---------------------------- don't change

            //put the parametrization of the scale factor here (copy/paste from BTagPOG payload)
            //do NOT use static (should not compile anyway)
            float SF = 1.00572
                    * ((1. + (0.013676 * x)) / (1. + (0.0143279 * x)));

            //don't change:
            //takes care of systematic variations (if any)
            return calcHeavySF(ptbins, SFb_error, ptbinsSize, x, abs_eta,  SF, multiplier);
        } else if (wp_ == csvm_wp) {  // ############### CSV Medium ################

            // -------
            // SF from:
            // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt  (EPS 2013 prescription)
            // THIS IS WINTER13!!!
            // -------

            static float ptbins[] = { 20, 30, 40, 50, 60, 70, 80, 100, 120, 160,
                    210, 260, 320, 400, 500, 600, 800 };

            //put the errors to the corresponding SF here (copy/paste from BTagPOG payload)
            //make sure to use static and that it has the same size as ptbins!!!
            static float SFb_error[] = { 0.0415694, 0.023429,0.0261074, 0.0239251,
                    0.0232416, 0.0197251, 0.0217319, 0.0198108,
                    0.0193, 0.0276144, 0.0205839, 0.026915,
                    0.0312739, 0.0415054, 0.0740561, 0.0598311};

            //don't change ---------------------------->
            static size_t ptbinsSize = sizeof(ptbins) / sizeof(ptbins[0]);
            static size_t SFb_errorSize = sizeof(SFb_error) / sizeof(SFb_error[0]);
            if (SFb_errorSize != ptbinsSize - 1) {
                std::cout
                << "bTagBase::BJetSF: Size of SFb_error should be one less than of ptbins. throwing exception!"
                << std::endl;
                throw std::logic_error(
                        "bTagBase::BJetSF: Size of SFb_error should be one less than of ptbins.");
            }
            //force range according to BTagPOG input (dont change)
            if (x < ptbins[0]) {
                x = ptbins[0];
                multiplier *= 2;
            } else if (x >= ptbins[ptbinsSize - 1]) {
                x = ptbins[ptbinsSize - 1] - 0.0001;
                multiplier *= 2;
            }
            // <---------------------------- don't change

            //put the parametrization of the scale factor here (copy/paste from BTagPOG payload)
            //do NOT use static (should not compile anyway)
            float SF = (0.939158 + (0.000158694 * x))+( -2.53962e-07 * (x * x));

            //don't change:
            //takes care of systematic variations (if any)
            return calcHeavySF(ptbins, SFb_error, ptbinsSize, x, abs_eta,  SF, multiplier);
        } else if (wp_ == csvt_wp) { // ############### CSV Tight ################

            // -------
            // SF from:
            // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt  (EPS 2013 prescription)
            // -------

            static float ptbins[] = { 20, 30, 40, 50, 60, 70, 80, 100, 120, 160,
                    210, 260, 320, 400, 500, 600, 800 };

            //put the errors to the corresponding SF here (copy/paste from BTagPOG payload)
            //make sure to use static and that it has the same size as ptbins!!!
            static float SFb_error[] = { 0.0511028, 0.0306671, 0.0317498, 0.032779,
                    0.0291528, 0.0249308, 0.0301118, 0.032047,
                    0.0348072, 0.0357745, 0.0378756, 0.0412608,
                    0.0777516, 0.0860741, 0.0942209, 0.104106};

            //don't change ---------------------------->
            static size_t ptbinsSize = sizeof(ptbins) / sizeof(ptbins[0]);
            static size_t SFb_errorSize = sizeof(SFb_error) / sizeof(SFb_error[0]);
            if (SFb_errorSize != ptbinsSize - 1) {
                std::cout
                << "bTagBase::BJetSF: Size of SFb_error should be one less than of ptbins. throwing exception!"
                << std::endl;
                throw std::logic_error(
                        "bTagBase::BJetSF: Size of SFb_error should be one less than of ptbins.");
            }
            //force range according to BTagPOG input (dont change)
            if (x < ptbins[0]) {
                x = ptbins[0];
                multiplier *= 2;
            } else if (x >= ptbins[ptbinsSize - 1]) {
                x = ptbins[ptbinsSize - 1] - 0.0001;
                multiplier *= 2;
            }
            // <---------------------------- don't change

            //put the parametrization of the scale factor here (copy/paste from BTagPOG payload)
            //do NOT use static (should not compile anyway)
            float SF = (0.9203 + (-3.32421e-05 * x)) + (-7.74664e-08 * (x * x));

            //don't change:
            //takes care of systematic variations (if any)
            return calcHeavySF(ptbins, SFb_error, ptbinsSize, x, abs_eta,  SF, multiplier);
        }

        //----------------END OF WORKING POINTS-----------------//
        else { //never reached
            throw std::logic_error(
                    "bTagBase: no working point set (this should never happen so someone has screwed up the code)");
        }

    }
    return 0; //never reached
}

float bTagBase::CJetSF(const float &pt, const float &abs_eta, float multiplier) const {
    if (wp_ == csvl_wp) {
        // Please give some information about the SF here. E.g
        // -------
        // SF from:
        // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG (EPS 2013 prescription)
        // -------
        //
        // the same scale factor is supposed to be used with twice the uncertainty (multiplier*=2)
        return BJetSF(pt, abs_eta, multiplier * 2.);
    }
    else if (wp_ == csvm_wp) {
        // Please give some information about the SF here. E.g
        // -------
        // SF from:
        // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG (EPS 2013 prescription)
        // -------
        //
        // the same scale factor is supposed to be used with twice the uncertainty (multiplier*=2)
        return BJetSF(pt, abs_eta, multiplier * 2.);
    }
    else if (wp_ == csvt_wp) {
        // Please give some information about the SF here. E.g
        // -------
        // SF from:
        // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG (EPS 2013 prescription)
        // -------
        //
        // the same scale factor is supposed to be used with twice the uncertainty (multiplier*=2)
        return BJetSF(pt, abs_eta, multiplier * 2.);
    }

    //----------------END OF WORKING POINTS-----------------//
    return 0;
}

float bTagBase::LJetSF(const float& pt, const float& abs_eta, float multiplier) const {

    //don't change ---------------------------->
    float x = pt;

    systematics tempsyst = syst_;

    // For the time being use the median points of b-falvour jets!!!
    if((tempsyst == lightuppt && pt > medianvecp_->at(bpt)) ||
            (tempsyst == lightdownpt && pt < medianvecp_->at(bpt)) ||
            (tempsyst == lightupeta && abs_eta > medianvecp_->at(beta)) ||
            (tempsyst == lightdowneta && abs_eta < medianvecp_->at(beta))){
        tempsyst = lightup;
    } else if((tempsyst == lightuppt && pt < medianvecp_->at(bpt)) ||
            (tempsyst == lightdownpt && pt > medianvecp_->at(bpt)) ||
            (tempsyst == lightupeta && abs_eta < medianvecp_->at(beta)) ||
            (tempsyst == lightdowneta && abs_eta > medianvecp_->at(beta))) {
        tempsyst = lightdown;
    }
    //<---------------------------- don't change
    if (is2011_) {
        if (x >= 670)
            x = 669; //approx but should be ok for ttbar/Z etc
        //for eta > 2.4, the value for 2.4 is taken (should also be ok)
        return 1;
    }

    else {
        if (wp_ == csvl_wp) {
            // Please give some information about the SF here. E.g
            // -------
            // SF from:
            // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG (EPS 2013 prescription)
            // From https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_EPS2013.C  (EPS 2013 prescription)
            // -------
            //
            // make sure, the range is appropriate (check the prescription)

            if (x > 800 && abs_eta <= 1.5)
                x = 800;
            else if (x > 700 && abs_eta > 1.5)
                x = 700;

            // Please copy-paste the formulas here. Check the correct eta range
            // systematics are handled directly here. use "else if" where applicable

            if (tempsyst != lightdown && tempsyst != lightup) { //don't change
                if (abs_eta <= 0.5)
                    return ((1.01177 + (0.0023066 * x)) + (-4.56052e-06 * (x * x))) + (2.57917e-09 * (x * (x * x)));
                else if (abs_eta <= 1.0)
                    return ((0.975966 + (0.00196354 * x)) + (-3.83768e-06 * (x * x))) + (2.17466e-09 * (x * (x * x)));
                else if (abs_eta <= 1.5)
                    return ((0.93821 + (0.00180935 * x)) + (-3.86937e-06 * (x * x))) + (2.43222e-09 * (x * (x * x)));
                else
                    return ((1.00022 + (0.0010998 * x)) + (-3.10672e-06 * (x * x))) + (2.35006e-09 * (x * (x * x)));
            } else if (tempsyst == lightdown) { //don't change
                if (abs_eta <= 0.5)
                    return ((0.977761 + (0.00170704 * x)) + (-3.2197e-06 * (x * x))) + (1.78139e-09 * (x * (x * x)));
                else if (abs_eta <= 1.0)
                    return ((0.945135 + (0.00146006 * x)) + (-2.70048e-06 * (x * x))) + (1.4883e-09 * (x * (x * x)));
                else if (abs_eta <= 1.5)
                    return ((0.911657 + (0.00142008 * x)) + (-2.87569e-06 * (x * x))) + (1.76619e-09 * (x * (x * x)));
                else
                    return ((1.03039 + (0.0013358 * x)) + (-3.89284e-06 * (x * x))) + (3.01155e-09 * (x * (x * x)));
            } else if (tempsyst == lightup) { //don't change
                if (abs_eta <= 0.5)
                    return ((1.04582 + (0.00290226 * x)) + (-5.89124e-06 * (x * x))) + (3.37128e-09 * (x * (x * x)));
                else if (abs_eta <= 1.0)
                    return ((1.00683 + (0.00246404 * x)) + (-4.96729e-06 * (x * x))) + (2.85697e-09 * (x * (x * x)));
                else if (abs_eta <= 1.5)
                    return ((0.964787 + (0.00219574 * x)) + (-4.85552e-06 * (x * x))) + (3.09457e-09 * (x * (x * x)));
                else
                    return ((1.1388 + (0.000468418 * x)) + (-1.36341e-06 * (x * x))) + (1.19256e-09 * (x * (x * x)));
            }
        } else if (wp_ == csvm_wp) {
            // Please give some information about the SF here. E.g
            // -------
            // SF from:
            // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG (EPS 2013 prescription)
            // From https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_EPS2013.C  (EPS 2013 prescription)
            // -------
            //
            // make sure, the range is appropriate (check the prescription)

            if (x > 800 && abs_eta <= 1.5)
                x = 800;
            else if (x > 700 && abs_eta > 1.5)
                x = 700;

            // Please copy-paste the formulas here. Check the correct eta range
            // systematics are handled directly here. use "else if" where applicable

            if (tempsyst != lightdown && tempsyst != lightup) { //don't change
                if (abs_eta <= 0.8)
                    return ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
                else if (abs_eta <= 1.6)
                    return ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
                else
                    return ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
            } else if (tempsyst == lightdown) { //don't change
                if (abs_eta <= 0.8)
                    return ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
                else if (abs_eta <= 1.6)
                    return ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
                else
                    return ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
            } else if (tempsyst == lightup) { //don't change
                if (abs_eta <= 0.8)
                    return ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
                else if (abs_eta <= 1.6)
                    return ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
                else
                    return ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
            }
        } else if (wp_ == csvt_wp) {
            // Please give some information about the SF here. E.g
            // -------
            // SF from:
            // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG (EPS 2013 prescription)
            // From https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_EPS2013.C  (EPS 2013 prescription)
            // -------
            //
            // make sure, the range is appropriate (check the prescription)

            if (x > 800 && abs_eta <= 1.5)
                x = 800;
            else if (x > 700 && abs_eta > 1.5)
                x = 700;

            // Please copy-paste the formulas here. Check the correct eta range
            // systematics are handled directly here. use "else if" where applicable

            if (tempsyst != lightdown && tempsyst != lightup) { //don't change
                return ((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
            } else if (tempsyst == lightdown) { //don't change
                return ((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
            } else if (tempsyst == lightup) { //don't change
                return ((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)));
            }
        }

        //----------------END OF WORKING POINTS-----------------//

    }
    return 0; //never reached
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////CSV Shape reweighting input     //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// these functions define the IO,    ///////////////////////////////////
////////////////////////////////////////     bins and naming conventions   ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



TString bTagBase::assoHFSysHistoName(systematics sys)const {

    switch(sys){
    case systematics::jesup:         return "final_JESUp" ;
    case systematics::jesdown:        return "final_JESDown" ;
    case systematics::lightup:       return "final_LFUp" ;
    case systematics::lightdown:     return "final_LFDown" ;
    case systematics::hfstat1up:   return "final_Stats1Up" ;
    case systematics::hfstat1down: return "final_Stats1Down" ;
    case systematics::hfstat2up:   return "final_Stats2Up" ;
    case systematics::hfstat2down:  return "final_Stats2Down" ;


    default : return "final";
    }

}
TString bTagBase::assoCSysHistoName(systematics sys)const{
    switch(sys){
    case systematics::cerr1up:   return "final_cErr1Up" ;
    case systematics::cerr1down: return "final_cErr1Down" ;
    case systematics::cerr2up:   return "final_cErr2Up" ;
    case systematics::cerr2down: return "final_cErr2Down" ;


    default :  return "final";
    }
}
TString bTagBase::assoLFSysHistoName(systematics sys)const{

    switch(sys){
    case systematics::jesup:            return "final_JESUp" ;
    case systematics::jesdown:         return "final_JESDown" ;
    case systematics::heavyup:         return "final_HFUp" ;
    case systematics::heavydown:       return "final_HFDown" ;
    case systematics::lfstat1up:    return "final_Stats1Up" ;
    case systematics::lfstat1down: return "final_Stats1Down" ;
    case systematics::lfstat2up:    return "final_Stats2Up" ;
    case systematics::lfstat2down: return "final_Stats2Down" ;


    default : return "final";
    }
}




void bTagBase::setupShapeReweightingBins(){
    shapeRWPtBinsHF_  << 30 << 40 << 60 << 100 << 160 << 140000;

    shapeRWPtBinsLF_  << 30 << 40 << 60 << 140000;
    shapeRWEtaBinsLF_ << 0 << 0.8 << 1.6 << 6; //<< 2.4;


    ////the following works automativally
    TH1::AddDirectory(false);

    h_csv_wgt_hf_.resize(shapeRWPtBinsHF_.size()-1);
    c_csv_wgt_hf_.resize(shapeRWPtBinsHF_.size()-1);

    h_csv_wgt_lf_.resize(shapeRWPtBinsLF_.size()-1);
    for(size_t i=0;i<h_csv_wgt_lf_.size();i++){
        h_csv_wgt_lf_.at(i).resize(shapeRWEtaBinsLF_.size()-1);
    }


}

}
