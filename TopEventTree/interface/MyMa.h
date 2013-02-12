#include "TString.h"

//#include <iostream>

namespace MyMa{

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
    TString p1    = "(("+TString(vec1)+".X()*"+TString(vec1)+".X()) + ("+TString(vec1)+".Y()*"+TString(vec1)+".Y()) + ("+TString(vec1)+".Z()*"+TString(vec1)+".Z()))";
    TString p2    = "(("+TString(vec2)+".X()*"+TString(vec2)+".X()) + ("+TString(vec2)+".Y()*"+TString(vec2)+".Y()) + ("+TString(vec2)+".Z()*"+TString(vec2)+".Z()))";
    TString ptot2 = "(TMath::Sqrt("+p1+"*"+p2+"))";
    TString dot   = "(("+TString(vec1)+".X()*"+TString(vec2)+".X()) + ("+TString(vec1)+".Y()*"+TString(vec2)+".Y()) + ("+TString(vec1)+".Z()*"+TString(vec2)+".Z()))";
    TString arcos = "TMath::ACos("+dot+"/"+ptot2+")";
    //std::cout << arg << std::endl;
    return arcos;
  }
}
