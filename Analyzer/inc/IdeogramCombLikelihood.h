#ifndef IDEOGRAMCOMBLIKELIHOOD_H
#define IDEOGRAMCOMBLIKELIHOOD_H

#include "ProgramOptionsReader.h"

class IdeogramCombLikelihood {
public:
  IdeogramCombLikelihood();
  virtual ~IdeogramCombLikelihood() {};

  virtual double Evaluate(double *x, double *p) = 0;

  std::vector<double> readParameters(const char *whichParameter);
  void SetFixedParams(double p0 = 0, double p1 = 0, double p2 = 0, double p3 = 0, double p4 = 0, double p5 = 0, double p6 = 0);
  double GetFixedParam(int index);
  void SetActive(bool active);
  bool IsActive();

protected:
  std::vector<double> parsCP_;
  std::vector<double> parsWP_;
  std::vector<double> parsUN_;
  std::vector<double> parsCPJES_;
  std::vector<double> parsWPJES_;
  std::vector<double> parsUNJES_;

  std::vector<double> massOffset_;
  std::vector<double> massSlopeMass_;
  std::vector<double> massSlopeJES_;
  std::vector<double> massSlopeMassJES_;
  std::vector<double> jesOffset_;
  std::vector<double> jesSlopeMass_;
  std::vector<double> jesSlopeJES_;
  std::vector<double> jesSlopeMassJES_;

  double fCP_, fWP_, fUN_;

  std::vector<double> fp_;
  bool useFixedParams_;
  bool active_;

};

#endif
