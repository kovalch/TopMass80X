#ifndef IDEOGRAMCOMBLIKELIHOOD_H
#define IDEOGRAMCOMBLIKELIHOOD_H

#include <vector>

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
  static std::vector<double> parsCP_;
  static std::vector<double> parsWP_;
  static std::vector<double> parsUN_;
  static std::vector<double> parsCPJES_;
  static std::vector<double> parsWPJES_;
  static std::vector<double> parsUNJES_;

  static std::vector<double> massOffset_;
  static std::vector<double> massSlopeMass_;
  static std::vector<double> massSlopeJES_;
  static std::vector<double> massSlopeMassJES_;
  static std::vector<double> jesOffset_;
  static std::vector<double> jesSlopeMass_;
  static std::vector<double> jesSlopeJES_;
  static std::vector<double> jesSlopeMassJES_;

  static double fCP_;
  static double fWP_;
  static double fUN_;

  std::vector<double> fp_;
  bool useFixedParams_;
  bool active_;

};

#endif
