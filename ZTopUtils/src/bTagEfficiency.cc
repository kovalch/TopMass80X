/*
 * bTagEfficiency.cc
 *
 *  Created on: Apr 7, 2015
 *      Author: kiesej
 */


#include "../interface/bTagEfficiency.h"
#include <stdexcept>
#include "TMath.h"


namespace ztop{

bool bTagEfficiency::debug=false;

std::vector<std::string> bTagEfficiency::getJetHistoOrderedNames(){
	//setting up naming scheme:
	std::vector<std::string> out;
	out.   resize(hist_length_histoNames, "");
	// Associating histogram names
	out.at(hist_bjets_h)="bjets2D";
	out.at(hist_btagged_h)="bjetsTagged2D";
	out.at(hist_cjets_h)="cjets2D";
	out.at(hist_ctagged_h)="cjetsTagged2D";
	out.at(hist_ljets_h)="ljets2D";
	out.at(hist_ltagged_h)="ljetsTagged2D";
	return out;
}
std::vector<std::string> bTagEfficiency::getEffHistoOrderedNames(){
	std::vector<std::string>  out(effhist_length_effHistoNames);
	// Adding efficiency histogram names
	out.at(effhist_b)="beff2D";
	out.at(effhist_c)="ceff2D";
	out.at(effhist_udsg)="leff2D";
	return out;
}

bTagEfficiency::bTagEfficiency(): init_(false),makeeffs_(false){

}
bTagEfficiency::~bTagEfficiency(){

}

void bTagEfficiency::fillEff(const float& pt, const float&eta,
		const int& genpartflav, const bool &tagged,
		const float &puweight){
	if(!makeeffs_)
		return;
	if(!init_)
		initHistos();
	if(genpartflav==0)//protect against data
		return ;
	float abs_eta=fabs(eta);
	BTagEntry::JetFlavor jetflav=jetFlavor(genpartflav);
	histos_.at(assoFlavorToHist(jetflav,false)).Fill(pt,abs_eta,puweight);
	if(tagged){
		histos_.at(assoFlavorToHist(jetflav,true)).Fill(pt,abs_eta,puweight);
	}
}
void bTagEfficiency::makeEff(){
	if(!makeeffs_)
			return;
	for (unsigned int i = 0; i < effhistos_.size(); i++) {
		for (int binx = 1; binx <= effhistos_.at(i).GetNbinsX() + 1; binx++) {
			for (int biny = 1; biny <= effhistos_.at(i).GetNbinsY() + 1;
					biny++) {
				//
				//uses histos at 2i and 2i+1
				float cont = 1;
				float err = 0.99; //to avoid zeros!
				if (histos_.at(2 * i).GetBinContent(binx, biny) > 0) {
					cont = histos_.at(2 * i + 1).GetBinContent(binx, biny) /
							histos_.at(2 * i).GetBinContent(binx, biny);
					if (debug)
						std::cout << "makeEffs: content: " << cont;
					err = sqrt(
							cont * (1 - cont)
							/ histos_.at(2 * i).GetBinContent(binx, biny));
					if (debug)
						std::cout << " error: " << err << "  " << binx << " "
						<< biny << std::endl;
				}

				effhistos_.at(i).SetBinContent(binx, biny, cont);
				effhistos_.at(i).SetBinError(binx, biny, err);

			}
		}

	}


	medians_.resize(med_length_medians, 0);
	medians_.at(med_bpt) = produceMedian(histos_.at(hist_btagged_h).ProjectionX());
	medians_.at(med_beta) = produceMedian(histos_.at(hist_btagged_h).ProjectionY());
	medians_.at(med_cpt) = produceMedian(histos_.at(hist_ctagged_h).ProjectionX());
	medians_.at(med_ceta) = produceMedian(histos_.at(hist_ctagged_h).ProjectionY());
	medians_.at(med_lpt) = produceMedian(histos_.at(hist_ltagged_h).ProjectionX());
	medians_.at(med_leta) = produceMedian(histos_.at(hist_ltagged_h).ProjectionY());

}


float bTagEfficiency::getMedian(bTagEfficiency::btag_medians med)const{
	if(medians_.size()<med_length_medians)
		throw std::logic_error("bTagEfficiency::getMedian: cannot get median that was never created or read.");
	return medians_.at(med);
}
float bTagEfficiency::jetEff(const float &pt, const float& abs_eta,const BTagEntry::JetFlavor& jetflav) const{
	if(effhistos_.size() < effhist_length_effHistoNames)
		throw std::out_of_range("bTagEfficiency::jetEff: no efficiency histos!");
	int ptbin = 0;
	int etabin = 0;
	int bla = 0;
	effHistoNames hist=assoFlavorToEff(jetflav);

	effhistos_.at(hist).GetBinXYZ(effhistos_.at(hist).FindFixBin(pt, abs_eta),
			ptbin, etabin, bla);
	return effhistos_.at(hist).GetBinContent(ptbin, etabin);
}

void  bTagEfficiency::writeToTFile(TFile *f){
	if(!makeeffs_)
		return;
	if(!init_)
		initHistos();
	makeEff();
	for(size_t i=0;i<effhistos_.size();i++)
		f->WriteTObject(&effhistos_.at(i));
	for(size_t i=0;i<histos_.size();i++)
		f->WriteTObject(&histos_.at(i));
	//medians
	TH1D tmp=TH1D("medians","medians",med_length_medians+1,0,med_length_medians);
	for(size_t i=0;i<medians_.size();i++)
		tmp.SetBinContent(i+1,medians_.at(i));
	f->WriteTObject(&tmp);
}
void  bTagEfficiency::readFromTFile(TFile * f){
	initHistos();
	for(size_t i=0;i<effhistos_.size();i++)
		effhistos_.at(i)=tryToGet<TH2D>(f,effhistos_.at(i).GetName());
	for(size_t i=0;i<histos_.size();i++)
		histos_.at(i)=tryToGet<TH2D>(f,histos_.at(i).GetName());

	TH1D tmp=tryToGet<TH1D>(f,"medians");
	for(size_t i=0;i<medians_.size();i++)
		medians_.at(i)=tmp.GetBinContent(i+1);
	makeeffs_=false;
}


///CHECK!!!!!!seems to not work
float bTagEfficiency::produceMedian(TH1 * h1)const{
	int nBin = h1->GetXaxis()->GetNbins();
	std::vector<double> x(nBin);
	h1->GetXaxis()->GetCenter(&x[0]);
	TH1D* h1D = dynamic_cast<TH1D*>(h1);
	if(!h1D){
		throw std::logic_error("bTagEfficiency::produceMedian: Median needs a TH1D!\n" );
	}
	const double* y = h1D->GetArray();
	// exclude underflow/overflows from bin content array y
	float median=TMath::Median(nBin, &x[0], &y[1]);
	return median;
}


BTagEntry::JetFlavor bTagEfficiency::jetFlavor(const int & partonflavor)const{
	if(std::abs(partonflavor) == 5) return BTagEntry::FLAV_B;
	else if (std::abs(partonflavor) == 4) return BTagEntry::FLAV_C;
	else if(std::abs(partonflavor) > 0) return BTagEntry::FLAV_UDSG;
	throw std::out_of_range(" bTagSFBase::jetFlavor: undefined parton flavor");
}


void bTagEfficiency::initHistos(){

	TH1::AddDirectory(false);

	float effptbins[] = { 20., 50., 70., 100., 160., 210., 800. };
	unsigned int npt = 7;
	float effetabins[] = { 0.0, 0.5, 1.0, 2.5 };
	unsigned int neta = 4;

	//probably there will be less statistics for the light jets, so histos might need a coarser binning

	float l_effptbins[] = { 20., 70., 120., 800. };
	unsigned int l_npt = 4;
	float l_effetabins[] = { 0.0, 1.5, 3.0 };
	unsigned int l_neta = 3;

	histos_.resize(hist_length_histoNames,TH2D());

	// Defining the histograms. Should have names as in histoNames_ and effHistoNames_
	histos_.at(hist_bjets_h) = TH2D(getJetHistoOrderedNames().at(hist_bjets_h).c_str(), "unTagged Bjets", npt - 1,
			effptbins, neta - 1, effetabins);
	histos_.at(hist_bjets_h).Sumw2();

	histos_.at(hist_btagged_h) = TH2D(getJetHistoOrderedNames().at(hist_btagged_h).c_str(), "Tagged Bjets",
			npt - 1, effptbins, neta - 1, effetabins);
	histos_.at(hist_btagged_h) .Sumw2();

	histos_.at(hist_cjets_h)  = TH2D(getJetHistoOrderedNames().at(hist_cjets_h).c_str(), "unTagged Cjets", npt - 1,
			effptbins, neta - 1, effetabins);
	histos_.at(hist_cjets_h) .Sumw2();

	histos_.at(hist_ctagged_h)  = TH2D(getJetHistoOrderedNames().at(hist_ctagged_h).c_str(), "Tagged Cjets",
			npt - 1, effptbins, neta - 1, effetabins);
	histos_.at(hist_ctagged_h) .Sumw2();

	histos_.at(hist_ljets_h)  = TH2D(getJetHistoOrderedNames().at(hist_ljets_h).c_str(), "unTagged Ljets", l_npt - 1,
			l_effptbins, l_neta - 1, l_effetabins);
	histos_.at(hist_ljets_h) .Sumw2();

	histos_.at(hist_ltagged_h)  = TH2D(getJetHistoOrderedNames().at(hist_ltagged_h).c_str(), "Tagged Ljets",
			l_npt - 1, l_effptbins, l_neta - 1, l_effetabins);
	histos_.at(hist_ltagged_h) .Sumw2();

	effhistos_.resize(effhist_length_effHistoNames,TH2D());

	effhistos_.at(effhist_b) = TH2D(getEffHistoOrderedNames().at(effhist_b).c_str(), "Bjets eff", npt - 1, effptbins,
			neta - 1, effetabins);
	effhistos_.at(effhist_b).Sumw2();

	effhistos_.at(effhist_c) = TH2D(getEffHistoOrderedNames().at(effhist_c).c_str(), "Cjets eff", npt - 1, effptbins,
			neta - 1, effetabins);
	effhistos_.at(effhist_c).Sumw2();

	effhistos_.at(effhist_udsg) = TH2D(getEffHistoOrderedNames().at(effhist_udsg).c_str(), "Ljets eff", l_npt - 1, l_effptbins,
			l_neta - 1, l_effetabins);
	effhistos_.at(effhist_udsg).Sumw2();

	medians_.resize(med_length_medians,0);

	init_=true;
}


}

