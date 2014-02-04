#ifndef SKIMMER_C
#define SKIMMER_C

#include "string"

class Skimmer {
public:
  Skimmer();
  ~Skimmer() {}

private:

  void skim(std::string inputPath, std::string outputPath, std::string sample);
};

#endif /* SKIMMER_C */
