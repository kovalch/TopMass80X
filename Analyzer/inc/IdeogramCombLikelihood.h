#ifndef IDEOGRAMCOMBLIKELIHOO_H
#define IDEOGRAMCOMBLIKELIHOO_H
#include <map>
#include <utility>

class IdeogramCombLikelihood {
private:
  double PCP(double *x, double *p);
  double PWP(double *x, double *p);
  double PUN(double *x, double *p);
  
  double PJES(double *x, double *p, double *q);
  //double PCPJES(double *x, double *p, double *q);
  //double PWPJES(double *x, double *p, double *q);
  //double PUNJES(double *x, double *p, double *q);
  
  double PCPJES1(double *x, double *p);
  double PWPJES1(double *x, double *p);
  double PUNJES1(double *x, double *p);
  
  double PCPJES2(double *x, double *p);
  double PWPJES2(double *x, double *p);
  double PUNJES2(double *x, double *p);
  
  double PBKG(double *x, double *p);
  double PBKGJES(double *x, double *p, double *q);
  double PBKGJES1(double *x, double *p);
  double PBKGJES2(double *x, double *p);
  
  double PBKGintegral_;
  
  typedef std::pair<double,double> ScanPoint;
  typedef std::map<ScanPoint,double> ScanPointMap; 
  ScanPointMap PWPintegrals_;
  ScanPointMap PUNintegrals_;
  
public:
  IdeogramCombLikelihood() : PBKGintegral_(-1) {};
  
  double Evaluate(double *x, double *p);
};


#endif
