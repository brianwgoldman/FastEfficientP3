// Brian Goldman

#ifndef POWERLAW_H_
#define POWERLAW_H_

#include "Optimizer.h"
#include "Util.h"

// Inherits and implements the Optimizer interface
class PowerLaw : public Optimizer {
 public:
  PowerLaw(Random& _rand, shared_ptr<Evaluator> _evaluator,
           Configuration& _config);
  virtual bool iterate() override;create_optimizer(PowerLaw);

private:
  // Uses a population size of 1
  vector<bool> solution;
  float fitness;
  vector<double> precomputed_power_law;

  double mutation_rate() const;
  size_t sample_waiting_time(double p) const;
  vector<bool> mutate(const vector<bool>& parent);
};

#endif /* POWERLAW_H_ */
