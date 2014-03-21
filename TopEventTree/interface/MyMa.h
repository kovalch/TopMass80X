#include "TString.h"

//#include <iostream>

namespace MyMa{

static TString deltaR(const char* vec1, const char* vec2)
{
  return TString("sqrt(pow("+TString(vec1)+".Eta()-"+TString(vec2)+".Eta(),2) + pow(TVector2::Phi_mpi_pi("+TString(vec1)+".Phi()-"+TString(vec2)+".Phi()),2))");
}
static TString deltaPhi(const char* vec1, const char* vec2)
{
  return TString("abs(TVector2::Phi_mpi_pi("+TString(vec1)+".Phi()-"+TString(vec2)+".Phi()))");
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

static TString invariantMass(const char* vec1, const char* vec2, TString scalingFactorVec2 = "1."){

	TString EnergySquared = "TMath::Power(" + TString(vec1)+".E() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".E(),2)";
	TString PtSquared = "(TMath::Power(" + TString(vec1)+".Px() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Px(),2) + TMath::Power(" + TString(vec1)+".Py() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Py(),2) + TMath::Power(" + TString(vec1)+".Pz() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Pz(),2))";
	TString InvMass = "TMath::Sqrt("+ EnergySquared + " - " + PtSquared + ")";

	return	InvMass;
}

static TString invariantMass(const char* vec1, const char* vec2, TString scalingFactorVec1 = "1.", TString scalingFactorVec2 = "1."){

	TString EnergySquared = "TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".E() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".E(),2)";
	TString PtSquared = "(TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".Px() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Px(),2) + TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".Py() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Py(),2) + TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".Pz() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Pz(),2))";
	TString InvMass = "TMath::Sqrt("+ EnergySquared + " - " + PtSquared + ")";

	return	InvMass;
}

static TString invariantMassThreeVec(const char* vec1, const char* vec2, const char* vec3, TString scalingFactorVec1 = "1.", TString scalingFactorVec2 = "1.", TString scalingFactorVec3 = "1."){

	TString EnergySquared = "TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".E() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".E() + "+ TString(scalingFactorVec3) + "*" + TString(vec3) + ".E(),2)";
	TString PtSquared = "(TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".Px() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Px() + " + TString(scalingFactorVec3) + "*" + TString(vec3)+".Px(),2) + TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".Py() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Py() + " + TString(scalingFactorVec3) + "*" + TString(vec3)+".Py(),2) + TMath::Power(" + TString(scalingFactorVec1) + "*" + TString(vec1)+".Pz() + " + TString(scalingFactorVec2) + "*" + TString(vec2)+".Pz() + " + TString(scalingFactorVec3) + "*" + TString(vec3)+".Pz(),2))";
	TString InvMass = "TMath::Sqrt("+ EnergySquared + " - " + PtSquared + ")";

	return	InvMass;
}

}

// invariant mass

//sqrt(pow(fitVecs[0].E () + fitVecs[1].E () + jets[fitAssigns[2]].E (), 2) - pow(fitVecs[0].Px() + fitVecs[1].Px() + jets[fitAssigns[2]].Px(), 2) - pow(fitVecs[0].Py() + fitVecs[1].Py() + jets[fitAssigns[2]].Py(), 2) - pow(fitVecs[0].Pz() + fitVecs[1].Pz() + jets[fitAssigns[2]].Pz(), 2))

//sqrt(pow(fitVecs[3].E () + fitVecs[4].E () + jets[fitAssigns[5]].E (), 2) - pow(fitVecs[3].Px() + fitVecs[4].Px() + jets[fitAssigns[5]].Px(), 2) - pow(fitVecs[3].Py() + fitVecs[4].Py() + jets[fitAssigns[5]].Py(), 2) - pow(fitVecs[3].Pz() + fitVecs[4].Pz() + jets[fitAssigns[5]].Pz(), 2))

