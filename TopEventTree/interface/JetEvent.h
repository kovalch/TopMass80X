/*
 * JetEvent.h
 *
 *  Created on: Feb 6, 2013
 *      Author: eschliec
 */

#ifndef JETEVENT_H_
#define JETEVENT_H_

#include "TObject.h"
#include "TString.h"
#include "TLorentzVector.h"
#include <vector>

#include <iostream>

class JetEvent : public TObject {
public:

  JetEvent();
  void init();

  unsigned int nJet() { return jet.size(); }

  static TString deltaR(const char* vec1, const char* vec2)
  {
    return TString("sqrt(pow("+TString(vec1)+".Eta()-"+TString(vec2)+".Eta(),2) + pow(TVector2::Phi_mpi_pi("+TString(vec1)+".Phi()-"+TString(vec2)+".Phi())+TMath::Pi(),2))");
  }
  static TString deltaPhi(const char* vec1, const char* vec2)
  {
    return TString("TVector2::Phi_mpi_pi("+TString(vec1)+".Phi()-"+TString(vec2)+".Phi())+TMath::Pi()");
  }
  static TString deltaAlpha(const char* vec1, const char* vec2)
  {
    // return the angle w.r.t. another 3-vector
    TString p1    = "(("+TString(vec1)+".X()*"+TString(vec1)+".X()) + ("+TString(vec1)+".Y()*"+TString(vec1)+".Y()) + ("+TString(vec1)+".Z()*"+TString(vec1)+".Z()))";
    TString p2    = "(("+TString(vec2)+".X()*"+TString(vec2)+".X()) + ("+TString(vec2)+".Y()*"+TString(vec2)+".Y()) + ("+TString(vec2)+".Z()*"+TString(vec2)+".Z()))";
    TString ptot2 = "(TMath::Sqrt("+p1+"*"+p2+"))";
    TString dot   = "(("+TString(vec1)+".X()*"+TString(vec2)+".X()) + ("+TString(vec1)+".Y()*"+TString(vec2)+".Y()) + ("+TString(vec1)+".Z()*"+TString(vec2)+".Z()))";
    TString arg   = "TMath::ACos("+dot+"/"+ptot2+")";
    std::cout << arg << std::endl;
    return arg;
  }
  ClassDef(JetEvent,1);

  // JetEvent data

  std::vector<TLorentzVector> jet;

  std::vector<int> nConstituents;
  std::vector<int> nChargedHadrons;
  std::vector<int> nNeutralHadrons;
  std::vector<int> nElectrons;
  std::vector<int> nPhotons;
  std::vector<int> nMuons;
  std::vector<double> fChargedHadron;
  std::vector<double> fNeutralHadron;
  std::vector<double> fElectron;
  std::vector<double> fPhoton;
  std::vector<double> fMuon;

  std::vector<double> charge;
  std::vector<double> flavour;
  std::vector<double> bTagCSV;
  std::vector<double> gluonTag;

  std::vector<double> jerSF;
  std::vector<double> jesSF;
  std::vector<double> totalSF;

  std::vector<TVector2> pull;
  std::vector<TVector2> pullCharged;


private:
};

#endif /* JETEVENT_H_ */
