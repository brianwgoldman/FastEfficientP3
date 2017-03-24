// Brian Goldman

#include "PowerLaw.h"

// Constructs some tools used during evolution, performs initial evaluation
PowerLaw::PowerLaw(Random& _rand, shared_ptr<Evaluator> _evaluator,
                   Configuration& _config)
    : Optimizer(_rand, _evaluator, _config),
      precomputed_power_law(length / 2) {
  // create and evaluate initial solution
  solution = rand_vector(rand, length);
  fitness = evaluator->evaluate(solution);
  double beta = config.get<float>("beta");
  double running_sum = 0;
  for (size_t i = 0; i < precomputed_power_law.size(); i++) {
    precomputed_power_law[i] = pow((i + 1), -beta) + running_sum;
    running_sum = precomputed_power_law[i];
  }
  for (auto & value : precomputed_power_law) {
    value /= running_sum;
  }
}

size_t PowerLaw::sample_waiting_time(double p) const {
  double x = std::uniform_real_distribution<double>(0, 1)(rand);
  return 1 + static_cast<int>(std::log1p(-x) / std::log1p(-p));
}

double PowerLaw::mutation_rate() const {
  // sample by bisection
  double x = std::uniform_real_distribution<double>(0, 1)(rand);
  if (x < precomputed_power_law[0]) {
    return 1;
  }
  size_t low = 0;
  size_t high = precomputed_power_law.size() - 1;
  while (high - low > 1) {
    size_t mid = (low + high) / 2;
    if (x > precomputed_power_law[mid]) {
      low = mid;
    } else {
      high = mid;
    }
  }
  return static_cast<double>(high + 1) / length;
}

// Selects "flips" number of random locations to perform bit flips.
vector<bool> PowerLaw::mutate(const vector<bool>& parent) {
  vector<bool> mutant(solution);
  double p = mutation_rate();
  size_t index = sample_waiting_time(p);
  while (index < length) {
    mutant[index] = not mutant[index];
    index += sample_waiting_time(p);
  }
  return mutant;
}

// Performs a full generation of the algorithm.
bool PowerLaw::iterate() {
  const auto& mutant = mutate(solution);
  //print(mutant);
  float new_fitness = evaluator->evaluate(mutant);
  if (new_fitness >= fitness) {
    solution = mutant;
    fitness = new_fitness;
  }
  // This algorithm never reaches stagnation, so always return true
  return true;
}
