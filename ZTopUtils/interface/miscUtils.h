#ifndef ZTOPMISCUTILS_h
#define ZTOPMISCUTILS_h

#include "TString.h"
#include <sstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "TH1D.h"
#include "TH2D.h"
#include <fstream>
#include <cmath>

/*
 WHATEVER you add here, please don't use exit() in case an error occurs.
 replace it with either:
 - throw an exception (throw std::logic_error("sometext") or std::runtime_error("");)
 - return something (-1 or another int for dubugging)

 */
#ifndef ZTOP_COUTVAR_DEF
#define ZTOP_COUTVAR_DEF
#include <iostream>
#define ZTOP_COUTVAR(x) std::cout << #x << ": " <<x << std::endl;
#endif

namespace ztop {

void coutDateTime();

inline double square(const double & in){
	return in*in;
}


double getMaxVar(bool up, const double & upvar, const double & downvar, bool& anticorr);

template<class T>
long double gammaFunc(T in, float eps=0.00001){
	if(eps && rint(in)>1 &&fabs(rint(in) - in) < eps){
		long double out=1;
		for(int i=in-1;i>0;i--)
			out*=i;
		return out;
	}
	else{
		return tgammal(in);
	}
}

/**
 *
 * @param k expectation value
 * @param lambda scanned value
 * @return the log of a continuous poisson prob for lambda
 */
double logPoisson(const double & k, const double& lambda);


/**
 *
 * @param k expectation value
 * @param lambda scanned value
 * @return  continuous poisson prob for lambda
 */
double poisson(const double & k, const double& lambda);

/**
 * returns poisson(evalpoint) with maximum at centre value, but shape according to stat
 * stat is NOT stat2! its just the stat you e.g. get from a histogram by calling getBinError()
 */
double shiftedLnPoisson(const float & centreN, const float & stat, const float& evalpoint);
/**
 * includes statistics of prediction when comparing a data point to prediction
 * stat is NOT stat2! its just the stat you e.g. get from a histogram by calling getBinError()
 * same for mcstat
 */
double shiftedLnPoissonMCStat(const float & centre, const float & stat, const float & mcstat, const float& evalpoint,bool useold=true);


/**
 * following numbers and mass dependence provided in NNLO paper arXiv:1303.6254
 * errors are NOT returned in % (so e.g. 0.026)
 */
float getTtbarXsec(float topmass, float energy=8, float* scaleerr=0, float * pdferr=0);

/**
 * slope from Hathor, central value form Kidonakis as
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
 */
float getTWXsec(float topmass, double energy=8);

template<class t>
bool isApprox(t a, t b, double eps = 0.01) {
	if (fabs(a - b) < eps)
		return true;
	else
		return false;
}
template<class t>
bool isAbsApprox(t a, t b, double eps = 0.01) {
	return (isApprox(a, b, eps) || isApprox(a, -b, eps));
}

template<class t>
t square(const t & in) {
	return in * in;
}

template<class t>
TString toTString(t in) {
	std::ostringstream s;
	s << in;
	TString out = s.str();
	return out;
}

template<class t>
std::string toString(t in) {
	std::ostringstream s;
	s << in;
	std::string out = s.str();
	return out;
}

TH2D divideTH2DBinomial(TH2D &h1, TH2D &h2);
TH2D divideTH2D(TH2D &h1, TH2D &h2);

void addRelError(TH2D &h, double err);
void addRelError(TH1D &h, double err);

TH1D divideTH1DBinomial(TH1D &h1, TH1D &h2);

template<class t>
std::vector<unsigned int> getOrderAscending(typename std::vector<t> vec) {
	typename std::vector<t> vec2 = vec;
	std::sort(vec.begin(), vec.end());
	std::vector<unsigned int> out;
	for (unsigned int i = 0; i < vec.size(); i++) {
		for (unsigned int j = 0; j < vec2.size(); j++) {
			if (vec.at(i) == vec2.at(j)) {
				out.push_back(j);
				break;
			}
		}
	}
	return out;
}


/**
 *  *Proper Delta Phi implementation instead for Delta R calculation.
 *   *Takes smallest possible Delta Phi
 *    */
double deltaPhi(const double & a, const double & b);


template<class t, class u>
double dR(t &first, u &sec) {
	return sqrt(
			(first.eta() - sec.eta()) * (first.eta() - sec.eta())
			+ deltaPhi(first.phi(),sec.phi()) * deltaPhi(first.phi(),sec.phi()));
}
template<class t, class u>
double dR(t *first, u *sec) {
	return sqrt(
			(first->eta() - sec->eta()) * (first->eta() - sec->eta())
			+ deltaPhi(first->phi(),sec->phi()) * deltaPhi(first->phi(),sec->phi()));
}

template<class t, class u>
bool noOverlap(t *first, u *sec, double deltaR) {
	bool nooverlap = true;
	if ((deltaR * deltaR)
			> square(first->eta() - sec->eta())
			+ square(deltaPhi(first->phi(), sec->phi()))) {
		nooverlap = false;
	}
	return nooverlap;
}

template<class T, class U>
bool noOverlap(T * first, typename std::vector<U*> &vecsec, double deltaR) {
	bool nooverlap = true;
	for (size_t i = 0; i < vecsec.size(); i++) {
		if (!noOverlap(first, vecsec.at(i), deltaR)) {
			nooverlap = false;
			break;
		}
	}
	return nooverlap;
}

template<class t>
int isIn(t element, typename std::vector<t> vec) {
	int IsIn = -1;
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec[i] == element) {
			IsIn = i;
			break;
		}
	}
	return IsIn;
}

template<class T, class U>
std::vector<T>& operator<<(std::vector<T>& vec, const U& x) {
	vec.push_back((T) x);
	return vec;
}

template<class T, class U>
std::vector<T>& operator<<(std::vector<T>& vec, const std::vector<U> & x) {
	vec.insert(vec.end(), x.begin(), x.end());
	return vec;
}

void displayStatusBar(Long64_t event, Long64_t nEvents, int ndiv = 100, bool force=false);

template<class T>
bool allEqual(std::vector<T> vec, T val) {
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec.at(i) != val)
			return false;
	}
	return true;
}
template<class T>
bool NoneEqual(std::vector<T> vec, T val) {
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec.at(i) == val)
			return false;
	}
	return true;
}

bool fileExists(const char * filename);

/**
 * returns index of element in vector of class U with closest deltaR to input element.
 * If dRmax is set to 0, the smallest dR value is written to this input variable
 * return -1 if no match is found.
 * the input elements are NOT changed but passing by const& leads to errors in some cases
 */
template<class T, class U>
int getClosestInDR(T* element, std::vector<U*>& coll, double & dRmax = 999,
		const double & dptrel = 200) {
	double dRmin = 9999;
	if (dRmax)
		dRmin = dRmax;
	int idx = -1;
	for (size_t i = 0; i < coll.size(); i++) {
		double dr = dR(element, coll.at(i));
		if (dr < dRmin
				&& dptrel
				> fabs(element->pt() - coll.at(i)->pt())
				/ element->pt()) {
			dRmin = dr;
			idx = i;
		}
	}
	if (!dRmax)
		dRmax = dRmin;
	return idx;
}

template<class T, class U>
int getClosestInDR(T& element, std::vector<U>& coll, double & dRmax = 999,
		const double & dptrel = 200) {
	double dRmin = 9999;
	if (dRmax)
		dRmin = dRmax;
	int idx = -1;
	for (size_t i = 0; i < coll.size(); i++) {
		double dr = dR(element, coll.at(i));
		if (dr < dRmin
				&& dptrel
				> fabs(element.pt() - coll.at(i).pt()) / element.pt()) {
			dRmin = dr;
			idx = i;
		}
	}
	if (!dRmax)
		dRmax = dRmin;
	return idx;
}
/**
 * works like std::sort but returns a vector of indecies that has the following form:
 * vector.at(unsrotedindex) == sortedindex
 * switched off in ROOT!!!
 */
#ifndef __CINT__
template<typename _RandomAccessIterator, typename _Compare>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last,
		_Compare __comp);
/**
 * works like std::sort but returns a vector of indecies that has the following form:
 * vector.at(unsrotedindex) == sortedindex
 */
template<typename _RandomAccessIterator>
inline std::vector<size_t>
retsort(_RandomAccessIterator __first, _RandomAccessIterator __last);
#endif
}

#ifndef __CINT__
#include <algorithm>
#include <iterator>

template<typename _RandomAccessIterator, typename _Compare>
inline std::vector<size_t>
ztop::retsort(_RandomAccessIterator __first, _RandomAccessIterator __last,
		_Compare __comp){
	typedef typename std::iterator_traits<_RandomAccessIterator>::value_type
			_ValueType;
	std::vector<_ValueType> copy(__first,__last); //copy
	std::vector<size_t> sortedilo;
	std::sort(copy.begin(),copy.end(),__comp);

	for(_RandomAccessIterator it=copy.begin();it!=copy.end();++it){
			//get the position in input
			size_t pos=std::find(__first,__last,*it)-__first;
			sortedilo.push_back(pos);
	}
	std::copy(copy.begin(),copy.end(),__first);
	return sortedilo;
}


template<typename _RandomAccessIterator>
inline std::vector<size_t>
ztop::retsort(_RandomAccessIterator __first, _RandomAccessIterator __last){
	typedef typename std::iterator_traits<_RandomAccessIterator>::value_type
			_ValueType;
	std::vector<_ValueType> copy(__first,__last); //copy
	std::vector<size_t> sortedilo;
	std::sort(copy.begin(),copy.end());

	for(_RandomAccessIterator it=copy.begin();it!=copy.end();++it){
		//get the position in input
		size_t pos=std::find(__first,__last,*it)-__first;
		sortedilo.push_back(pos);
	}
	std::copy(copy.begin(),copy.end(),__first);
	return sortedilo;
}

#endif

#endif
