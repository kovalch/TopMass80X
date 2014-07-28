#include "../interface/miscUtils.h"
#include <cmath>
#include <stdexcept>

namespace ztop {

double logPoisson(const double & k, const double& lambda){
	return k * log(lambda) - lgamma(k+1.) - lambda;
}
double poisson(const double & k, const double& lambda){
	return exp(logPoisson(k,lambda));
}

/**
 * all inputs in input coordinate system! rescaling will be done by function
 * E.g. a however scaled point has 10 entries with stat=sqrt(stat2) of 2, then stat will be 2! Nothing else
 * Try to keep shifts small here! in case, rescale before
 */
double shiftedLnPoisson(const float & centre, const float & stat, const float& evalpoint){
	// if(centre<1) return -9999999999999;
	double out=0;
	double stat2=stat*stat;
	double shift=stat2/centre;
	double realstat2=centre / shift;//centre/stat2;
	double realeval = evalpoint/shift;
	out=logPoisson(realeval, realstat2);

	//keep integral
	out-=log(shift);
	return out;
}

/**
 * all inputs in input coordinate system! rescaling will be done by function
 * no stat vlaues are squared!
 * Try to keep shifts small here! in case, rescale before
 */
double shiftedLnPoissonMCStat(const float & centre, const float & stat, const float & mcstat, const float& evalpoint){

	//convolute two shifted poissons
	//conv: shiftedLnPoisson(centre, stat2, evalpoint) otimes shiftedLnPoisson(centre, mcstat2, evalpoint)




	double out=0;
	double stat2=stat*stat;
	double mcstat2=mcstat*mcstat;
	double shift=stat2/centre;
	// go to real prediction level
	// this way get rid of any scaling
	double realcentre = centre / shift; // = npred
	// double realstat2 = realcentre = centre / shift = centre2 / stat2 = npred
	double realmcstat2 = mcstat2 / shift / shift;
	double realevalpoint = evalpoint / shift;

	// adjusted formula  from caldwell, etc arxiv 1112.2593
	// P = a / b * c / d
	//logP = parta - partb + partc - partd
	// n: initial MC statistics
	// s: scale factor with N_pred = N_gen / s
	// o: observed value in scale ~N_pred

	double s=1 / (realmcstat2/(realcentre));//*shift))  ;
	double n=realcentre * s  ;
	double o=realevalpoint;

	double parta=0,partb=0,partc=0,partd=0;

	const double log4=1.386294361119891;

	// if(fabs(fabs(s)-1) < 0.0001 ){
	//      parta=0;
	//  }
	//  else{
	parta=(n+0.5)*log(s);
	//   }
	if(parta!=parta)parta=0;


	if(o==0){
		out= (n+0.5) * log(s/(s+1));
	}
	else{
		partb=(o+n+0.5) * log(1+s);
		partc=lgamma(n+1.) + lgamma(2*(o+n) +1);
		partd = o*log4 + lgamma(o+1) + lgamma(2*n +1)+ lgamma(n+o+1);

		out=parta - partb + partc - partd;

	}
	//out-=log(shift);
	if(out > 1){
		std::cout << "\ncentre: "<< centre
				<<  "\nrealcentre: "<< realcentre
				<< "\nstat2: "<<stat2
				<< "\nmcstat2: "<<mcstat2
				<< "\nrealmcstat2: " <<realmcstat2
				<< "\nrealevalpoint: "<<o
				<< "\nshift: "<< shift
				<< "\ns: "<< s
				<< "\nn: "<< n<<std::endl;

		throw std::runtime_error("shiftedLnPoissonMCStat: you exceeded the range, where the approximations used here are valid by far");
	}
	return out;
}

float getTtbarXsec(float topmass, float energy, float* scaleerr, float * pdferr){
	/*
	 * all numbers following arxiv 1303.6254
	 *
	 */
	float mref=173.3;
	float referencexsec=0;
	float deltam=topmass-mref;


	float a1=0,a2=0;

	if(isApprox(energy,8.f,0.01)){
		a1=-1.1125;
		a2=0.070778;
		referencexsec=245.8;
		if(scaleerr)
			*scaleerr=0.034;
		if(pdferr)
			*pdferr=0.026;
	}
	else if(isApprox(energy,7.f,0.01)){
		a1=-1.24243;
		a2=0.890776;
		referencexsec=172.0;
		if(scaleerr)
			*scaleerr=0.034;
		if(pdferr)
			*pdferr=0.028;
	}

	float reldm=mref/(mref+deltam);

	float out= referencexsec* (reldm*reldm*reldm*reldm) * (1+ a1*(deltam)/mref + a2*(deltam/mref)*(deltam/mref));

	return out;
}

void addRelError(TH2D &h, double err) {
	for (int binx = 1; binx <= h.GetNbinsX() + 1; binx++) {
		for (int biny = 1; biny <= h.GetNbinsY() + 1; biny++) {
			double add = h.GetBinContent(binx, biny) * err;
			double newerr = sqrt(
					pow(h.GetBinError(binx, biny), 2) + pow(add, 2));
			h.SetBinError(binx, biny, newerr);
		}
	}
}
void addRelError(TH1D &h, double err) {
	for (int binx = 1; binx <= h.GetNbinsX() + 1; binx++) {
		double add = h.GetBinContent(binx) * err;
		double newerr = sqrt(pow(h.GetBinError(binx), 2) + pow(add, 2));
		h.SetBinError(binx, newerr);
	}
}

void displayStatusBar(Long64_t event, Long64_t nEvents, int ndiv) {

	if ((event + 1) * ndiv % nEvents < ndiv) {
		int statusbar = (event + 1) * ndiv / nEvents;
		std::cout << "\r[";
		for (int i = 0; i < statusbar * 50 / ndiv; i++) {
			std::cout << "=";
		}
		for (int i = statusbar * 50 / ndiv; i < 50; i++) {
			std::cout << " ";
		}
		std::cout << "] " << statusbar * 100 / ndiv << "%   ";
		flush(std::cout);
		statusbar++;
	}
	if (event == 0) {
		std::cout << "[                                                  ] "
				<< "0%   ";
		flush(std::cout);
	}
}

TH2D divideTH2DBinomial(TH2D &h1, TH2D &h2) { //! out = h1 / h2
	TH2D out = h1;
	if (h1.GetNbinsX() != h2.GetNbinsX() || h1.GetNbinsY() != h2.GetNbinsY()) {
		std::cout
		<< "divideTH2DBinomial: Error! histograms must have same binning!"
		<< std::endl;
		return h1;
	}
	for (int binx = 1; binx <= h1.GetNbinsX() + 1; binx++) {
		for (int biny = 1; biny <= h1.GetNbinsY() + 1; biny++) {
			double cont = 0;
			double err = 1;
			if (h2.GetBinContent(binx, biny) != 0) {
				cont = h1.GetBinContent(binx, biny)
                                                                                                                        				/ h2.GetBinContent(binx, biny);
				err = sqrt(cont * (1 - cont) / h1.GetBinContent(binx, biny));
			}
			out.SetBinContent(binx, biny, cont);
			out.SetBinError(binx, biny, err);
		}
	}
	return out;
}

TH2D divideTH2D(TH2D &h1, TH2D &h2) {
	TH2D out = h1;
	if (h1.GetNbinsX() != h2.GetNbinsX() || h1.GetNbinsY() != h2.GetNbinsY()) {
		std::cout
		<< "divideTH2DBinomial: Error! histograms must have same binning!"
		<< std::endl;
		return h1;
	}
	for (int binx = 1; binx <= h1.GetNbinsX() + 1; binx++) {
		for (int biny = 1; biny <= h1.GetNbinsY() + 1; biny++) {
			double cont = 0;
			double err = 1;
			if (h2.GetBinContent(binx, biny) != 0) {
				cont = h1.GetBinContent(binx, biny)
                                                                                                                        				/ h2.GetBinContent(binx, biny);
				err = sqrt(
						pow(
								h1.GetBinError(binx, biny)
								/ h2.GetBinContent(binx, biny), 2)
								+ pow(
										(cont / h2.GetBinContent(binx, biny)
												* h2.GetBinError(binx, biny)),
												2));
			}
			out.SetBinContent(binx, biny, cont);
			out.SetBinError(binx, biny, err);
		}
	}
	return out;
}

TH1D divideTH1DBinomial(TH1D &h1, TH1D &h2) { //! out = h1 / h2
	TH1D out = h1;
	if (h1.GetNbinsX() != h2.GetNbinsX()) {
		std::cout
		<< "divideTH1DBinomial: Error! histograms must have same binning!"
		<< std::endl;
		return h1;
	}
	for (int binx = 1; binx <= h1.GetNbinsX() + 1; binx++) {
		double cont = 0;
		double err = 1;
		if (h2.GetBinContent(binx) != 0) {
			cont = h1.GetBinContent(binx) / h2.GetBinContent(binx);
			err = sqrt(cont * (1 - cont) / h1.GetBinContent(binx));
		}
		out.SetBinContent(binx, cont);
		out.SetBinError(binx, err);
	}
	return out;
}

bool fileExists(const char * filename) {
	std::ifstream FileTest(filename);
	bool exists = FileTest;
	FileTest.close();
	return exists;
}

}
