#include "../interface/miscUtils.h"
#include <cmath>


namespace ztop {

double logPoisson(const double & k, const double& lambda){
    return k * log(lambda) - lgamma(k+1.) - lambda;
}
double poisson(const double & k, const double& lambda){
    return exp(logPoisson(k,lambda));
}
double shiftedLnPoisson(const float & centre, const float & stat2, const float& evalpoint){
    // if(centre<1) return -9999999999999;
    double out=0;
    double shift=stat2/centre;

    out=logPoisson(evalpoint*shift, stat2);
    //keep integral
    out+=log(shift);
    return out;
}

//uses gaussian approx starting at 200 for mc and data: TBI
double shiftedLnPoissonMCStat(const float & centre, const float & stat2, const float & mcstat2, const float& evalpoint){

    //convolute two shifted poissons
    //conv: shiftedLnPoisson(centre, stat2, evalpoint) otimes shiftedLnPoisson(centre, mcstat2, evalpoint)

    // adjusted formula  from caldwell, etc arxiv 1112.2593
    // P = a / b * c / d
    //logP = parta - partb + partc - partd

    double shift=stat2/centre;
    double mcshift=mcstat2/centre;

    const double log4=1.386294361119891;

    double out=0;
    if(mcstat2< 200 || stat2<200){

        double s=mcshift /shift;
        double n=mcstat2 ;
        double o=evalpoint* shift;

        double parta=0,partb=0,partc=0,partd=0;

        if(s-1 < 0.01 ){
            parta=0;
        }
        else{
            parta=(n+0.5)*log(s);
        }


        if(o==0){
            out= (n+0.5) * log(s/(s+1));
        }
        else{
            partb=(o+n+0.5) * log(1+s);
            partc=lgamma(n+1.) + lgamma(2*(o+n) +1);
            partd = o*log4 + lgamma(o+1) + lgamma(2*n +1)+ lgamma(n+o+1);

            out=parta - partb + partc - partd;

        }
        out+=log(shift);
    }
    else{ //both gaussian
        out=(log(1/sqrt(2*(3.14159265358979323846)*(stat2/(shift*shift)+mcstat2/(mcshift*mcshift))))
                + (-((long double)(evalpoint-centre))*((long double)(evalpoint-centre))/((long double)(2*(stat2/(shift*shift)+mcstat2/(mcshift*mcshift))))));
    }
    //out=logPoisson(evalpoint*shift,stat2);
    //out+=log(shift);
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
