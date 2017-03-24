// Brian Goldman

// This file aggregates optimization methods
// so that they can be chosen through configuration options

#ifndef OPTIMIZATIONCOLLECTION_H_
#define OPTIMIZATIONCOLLECTION_H_

#include "Optimizer.h"
#include "LTGA.h"
#include "Pyramid.h"
#include "RandomRestartHC.h"
#include "LambdaLambda.h"
#include "HBOA.h"
#include "Popless.h"
#include "PowerLaw.h"
#include <unordered_map>

namespace optimize {
// Renaming of the function pointer used to create new optimization methods
using pointer=shared_ptr<Optimizer> (*)(Random& rand, shared_ptr<Evaluator> evaluator, Configuration& config);

// Lookup table translates strings to function pointers
static std::unordered_map<string, pointer> lookup( {
  { "LTGA", LTGA::create },
  { "Pyramid", Pyramid::create },
  { "RandomRestartHC", RandomRestartHC::create },
  { "LambdaLambda", LambdaLambda::create },
  { "HBOA", HBOA::create},
  { "Popless", Popless::create},
  { "PowerLaw", PowerLaw::create},
});
}

#endif /* OPTIMIZATIONCOLLECTION_H_ */
