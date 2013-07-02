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

// invariant mass

//sqrt(pow(fitVecs[0].E () + fitVecs[1].E () + jets[fitAssigns[2]].E (), 2) - pow(fitVecs[0].Px() + fitVecs[1].Px() + jets[fitAssigns[2]].Px(), 2) - pow(fitVecs[0].Py() + fitVecs[1].Py() + jets[fitAssigns[2]].Py(), 2) - pow(fitVecs[0].Pz() + fitVecs[1].Pz() + jets[fitAssigns[2]].Pz(), 2))

//sqrt(pow(fitVecs[3].E () + fitVecs[4].E () + jets[fitAssigns[5]].E (), 2) - pow(fitVecs[3].Px() + fitVecs[4].Px() + jets[fitAssigns[5]].Px(), 2) - pow(fitVecs[3].Py() + fitVecs[4].Py() + jets[fitAssigns[5]].Py(), 2) - pow(fitVecs[3].Pz() + fitVecs[4].Pz() + jets[fitAssigns[5]].Pz(), 2))

