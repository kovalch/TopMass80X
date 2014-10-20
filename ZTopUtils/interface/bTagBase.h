#ifndef BTAGBASE_H
#define BTAGBASE_H

#include <map>
#include "TH2D.h"
#include "TH1.h"
#include "TString.h"
#include "../interface/miscUtils.h"
#include <iostream>
#include <string>
class TFile;

namespace ztop {
/**
 ***** workflow to use event wise scale factors (only for at least one jet selection):
 *  -create object
 *  -set sample name              setSampleName(name)
 *  -fill efficiencies with jets  fillEff(p4,genPartonFlavour,PUweight)
 *  ---- get your event weight with your wrapper class!! ----
 *  ---functions to do that:
 *  - resetCounter()
 *  - countJet() (in jet loop)
 *  - getEventSF() (after jet loop - also resets the counter) returns 1 here
 *
 *
 *  -after event loop run makeEffs() once to create histograms
 *  - input / output should be organized by your wrapper class!
 *
 *  -load from file or use already filled object (implement in your wrapper!)
 *  -setMakeEff(false)
 * -(fillEff,makeEffs etc can remain in the loop/at the end of the loop, does nothing now)
 *  ---- get your event weight with your wrapper class!! ----
 *  ---functions to do that:
 *  - resetCounter()
 *  - countJet() (in jet loop)
 *  - getEventSF() (after jet loop - also resets the counter) now returns the proper SF!
 *
 *
 ***** workflow to use "random tagging" (works for N b-jets selection):
 *
 * -create efficiencies the same way as described above and load them
 * -instead of doing a resetCounter... combination, use bool jetIsTagged(...) function
 *  to tag/untag jets
 *
 *
 *
 *  WHATEVER you add as functions, please don't use exit() in case an error occurs.
 *  replace it with either:
 *  - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 *  - return something (-1 or another int for dubugging)
 *
 */
class bTagBase {
public:

    bTagBase();

    ~bTagBase();
    bTagBase(const bTagBase& );
    bTagBase & operator = (const bTagBase& rhs);
    void copyFrom(const bTagBase& );

    /**
     * enum for working points
     * should be named <tagger><wp>_wp
     *
     */
    enum workingPoints {
        csvt_wp,
        csvm_wp,
        csvl_wp,
        /* some other WP/taggers ..*/
        //
        /*don't touch this: needs to be last entry*/
        length_wp
    };

    /**
     * enum for systematic variations
     */
    enum systematics {
        nominal,
        
        // shape uncertainties for norm. diff measurements
        heavyup, heavydown,
        lightup, lightdown,
        heavyuppt, heavydownpt,
        heavyupeta, heavydowneta,
        lightuppt, lightdownpt,
        lightupeta, lightdowneta,

        // uncertainties on the discriminator shape reweighting
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration#Example_code
        jesup,jesdown,  //yes correlated to "standard jes"
        purityup,puritydown,
        hfstat1up,hfstat1down,
        hfstat2up,hfstat2down,
        lfstat1up,lfstat1down,
        lfstat2up,lfstat2down,

        // leave this as last entry
        length_syst
    };


    enum histoTypes {tag, eff};

    void setWorkingPoint(workingPoints wp) { wp_ = wp; }
    const float& getWPDiscrValue()const;
    std::string getWorkingPointString()const;

    workingPoints getWorkingPoint() const;
    /**
     * switches on SF for 7 TeV data
     */
    void setIs2011(bool is) {is2011_ = is;}

    void setSystematic(systematics sys);
    systematics getSystematic() const{ return syst_;}


    int setSampleName(const std::string &); //checks if effs should be made, if sample exists,..

    /**
     * enables the filling of histograms for efficiencies and disables
     * scale factor output
     */
    void setMakeEff(bool makee) {makeeffs_ = makee; }
    bool getMakeEff()const { return makeeffs_;}

    void fillEff(const float &, const float&, const int &, const float &,const float&);
    void makeEffs();

    /**
     * returns the number of nan sf
     */
    size_t getNanCount()const{return nancount_;}

    /**
     * Input as .root files with histograms. See:
     * https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration
     * can be found in ZTopUtils/data/
     * - csv_rwt_hf.root -> heavy flavour
     * - csv_rwt_lf.root -> light flavour
     */
    void readShapeReweightingFiles(const TString& heavy,const TString& light);

    /**
     * this function gets a jet weight based on CSV reweighting for MC
     * if genPartonFlavour == 0 it just returns 1 (always true for data)
     * For each jet you consider for your total discriminator distribution,
     * multiply the previous weight to obtain an event weight.
     */
    float getJetDiscrShapeWeight(const float & pt, const float& abs_eta,
            const int & genPartonFlavor, const float& jetdiscr)const;


    bool showWarnings;
    bool debug;

protected:

    /**
     * use these functions to fill the TH2D vectors in histos_ and effhistos_ in the right order
     *
     */
    const std::vector<std::string> & getJetHistoOrderedNames()const{return histoNames_;}
    const std::vector<std::string> & getEffHistoOrderedNames()const{return effHistoNames_;}

    /**
     * use these enums to associate the entries of histos_ and effhistos_ with the given names from
     * e.g. getJetHistoOrderedNames()
     * so getJetHistoOrderedNames().at(bjets_h) will be the name of the b-jet histogram!
     * And bjets_h will be the index of the histogram in histos_[<samplename>]
     */
    enum histoNames_en{bjets_h, btagged_h, cjets_h, ctagged_h, ljets_h, ltagged_h, length_histoNames_en};
    enum effHistoNames_en{beff_h, ceff_h, leff_h, length_effHistoNames_en};

    std::map<std::string, std::vector<TH2D> > histos_;      //! bjets, btagged, cjets, ctagged, ljets, ltagged
    std::map<std::string, std::vector<TH2D> > effhistos_;   //! beff, ceff, leff


    /**
     * use these enums to associate the entries of the median vector in medianMap_ with the corresponding
     * input/output
     */
    enum medians {bpt, beta, cpt, ceta, lpt, leta, length_median};
    std::map<std::string, std::vector<float> > medianMap_;

    // functions to be used in the wrapper classes inheriting from this class
    // example implementation:
    /*
     resetCounter();
     for(size_t i=0;i<jets.size();i++){
     ztop::NTJet *jet=jets.at(i);
     countJet(fabs(jet->eta()),jet->pt(),jet->genPartonFlavour());
     }
     return getEventSF();
     */

    void resetCounter() {
        sumStuffEff_ = 0.99999999;
        sumStuffSfEff_ = 0.99999999;
    }
    /**
     *
     */
    void countJet(const float&, const float&, const int & genPartonFlavor);
    float getEventSF();

    /**
     * function to be used if random tagging is applied
     */
    bool jetIsTagged(const float&, const float&, const int & genPartonFlavor,
            const float & tagValue, const unsigned int & seed) const; // to be implmented



    void cleanptr() {
        histp_ = 0;
        effhistp_ = 0;
    }

    float median(TH1 *)const ;

    std::string histoNameAtId(const int id, const histoTypes type)const {
        if(type == tag) return (int)histoNames_.size()>id ? histoNames_.at(id) : std::string("");
        if(type == eff) return (int)effHistoNames_.size()>id ? effHistoNames_.at(id) : std::string("");
        else return std::string("");
    }




private:
    enum jetTypes{bjet,cjet,lightjet,undefined};

    bool init_;

    workingPoints wp_;
    bool is2011_;
    systematics syst_;

    std::vector<float> wpvals_, minpt_, maxpt_;

    std::vector<std::string> histoNames_;
    std::vector<std::string> effHistoNames_;

    bool makeeffs_;

    std::string tempsamplename_;
    std::vector<TH2D> * histp_;
    std::vector<TH2D> * effhistp_;
    std::vector<float> * medianvecp_;

    float sumStuffEff_;
    float sumStuffSfEff_;

    //counts the fraction of nans
    size_t nancount_;




    jetTypes jetType(const int & )const;

    float calcHeavySF(float *, float *, const size_t &, const float &,
            const float &, const float &, const float &) const;

    float jetSF(const float &pt, const float& abs_eta,const jetTypes & jettype) const;
    float jetEff(const float &pt, const float& abs_eta,const jetTypes & jettype) const;

    //SFs btv input
    float BJetSF(const float &pt, const float& abs_eta,
            float multiplier = 1) const;
    float CJetSF(const float &pt, const float &abs_eta,
            float multiplier = 1) const;
    float LJetSF(const float &pt, const float &abs_eta,
            float multiplier = 1) const;

    void initWorkingpoints();

    /////////////////////////////////////////////////////////////////////////////////
    //////////////////////// private member (functions) for shape reweighting ///////
    /////////////////////////////////////////////////////////////////////////////////

    /**
     * this should be called in the constructor and sets up
     * the binning of the reweight histograms in pt (and eta)
     * bins are defined in the cc file according to the format of the input files
     */
    void setupShapeReweightingBins();

    /**
     * here the binning is stored of the shape histograms in pt and eta
     * let the last bin boundary extremely high since there is no overflow provided
     */
    std::vector<float> shapeRWPtBinsHF_,shapeRWPtBinsLF_,shapeRWEtaBinsLF_;

    /**
     * Pointers to the heavy flavour and light flavour input files
     *   are set within the read in functions to files in memory
     */
    TFile * hf_rewfile_, *lf_rewfile_;

    /**
     * vectors of the histograms. One for each pt (and eta) bin.
     * These vectors contain only the histograms for the selected systematic
     */
    std::vector<TH1D> h_csv_wgt_hf_,c_csv_wgt_hf_;
    std::vector<std::vector<TH1D> > h_csv_wgt_lf_;


    /**
     * These function return pointers to the right histogram depending on
     * pt (and eta)
     */
    const TH1D * getShapeLFHisto(const float & pt, const float & fabseta)const;
    const TH1D * getShapeHFHisto(const float & pt, bool bjet)const;

    /**
     * Helper functions to associate a systematics enum to a given histogram name
     * as provided in the input file
     */
    TString assoHFSysHistoName(systematics sys)const;
    TString assoCSysHistoName(systematics sys)const;
    TString assoLFSysHistoName(systematics sys)const;

    /**
     * small helper function to read in a histogram
     * throws an exception if histogram is not found
     */
    TH1D  tryToGet(TFile * f, TString name)const;

};

///inlined functions for performance reasons:
inline bTagBase::jetTypes bTagBase::jetType(const int & partonflavor)const{
    if(std::abs(partonflavor) == 5) return bjet;
    else if (std::abs(partonflavor) == 4) return cjet;
    else if(std::abs(partonflavor) > 0) return lightjet;
    else return undefined;
}

inline float bTagBase::calcHeavySF(float* ptbins, float * SFb_error,
        const size_t & ptbinsSize, const float & pt, const float & abseta, const float & SF,
        const float & multiplier) const {

    if (syst_ != heavyup && syst_ != heavydown && syst_ != heavyuppt && syst_ != heavydownpt && syst_ != heavyupeta && syst_ != heavydowneta)
        return SF;

    //this uses the standard histogram bin range definition from root
    size_t ptbin = std::lower_bound(ptbins, ptbins + ptbinsSize, pt)
    - &ptbins[0];
    if (ptbins[ptbin] == pt)
        --ptbin;

    if (syst_ == heavyup ||
            (syst_ == heavyuppt && pt < medianvecp_->at(bpt)) || (syst_ == heavydownpt && pt > medianvecp_->at(bpt)) ||
            (syst_ == heavyupeta &&  abseta < medianvecp_->at(beta)) || (syst_ == heavydowneta && abseta > medianvecp_->at(beta)) )
        return SF + (multiplier * SFb_error[ptbin]);
    if (syst_ == heavydown ||
            (syst_ == heavyuppt && pt > medianvecp_->at(bpt)) || (syst_ == heavydownpt && pt < medianvecp_->at(bpt)) ||
            (syst_ == heavyupeta &&  abseta > medianvecp_->at(beta)) || (syst_ == heavydowneta && abseta < medianvecp_->at(beta)) )
        return SF - (multiplier * SFb_error[ptbin]);

    return 0;    //never reaches
}

}    //namespace
////////////////////////////////////
#endif

