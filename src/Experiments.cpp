// Brian Goldman

// There are a couple of different kind of experiment performed with
// this software, each producing different types of data.  These
// are the drivers of optimization.

#include "Experiments.h"

// A single run of optimization
Record single_run(Random& rand, Configuration& config,
                  evaluation::pointer problem, optimize::pointer solver,
                  int run) {
  size_t limit = config.get<int>("eval_limit");
  // Middle Layer's sit between the problem and the solver, tracking optimization
  shared_ptr<Middle_Layer> recorder(new Middle_Layer(config, problem(config, run)));

  // Start the clock before any evaluations can be performed
  recorder->results.start_clock();
  auto optimizer = solver(rand, std::static_pointer_cast<Evaluator>(recorder), config);
  float good_enough = config.get<float>("fitness_limit");

  // Iterate the optimizer until the solution is reached, the maximum number
  // of evaluations is performed, or the optimizer reaches stagnation
  bool improvement_possible = true;
  while (recorder->best_fitness < good_enough and recorder->evaluations < limit
      and improvement_possible) {
    improvement_possible = optimizer->iterate();
  }
  recorder->results.metadata = optimizer->finalize();
  return recorder->results;
}

// Performs multiple runs of optimization
vector<Record> multirun(Random& rand, Configuration& config,
                        evaluation::pointer problem, optimize::pointer solver) {
  int runs = config.get<int>("runs");
  int verbosity = config.get<int>("verbosity");
  vector<Record> records;
  // perform the desired number of runs
  for (int run = 0; run < runs; run++) {
    records.push_back(single_run(rand, config, problem, solver, run));
    // output partial summaries
    if (verbosity > 0) {
      auto summary = Record::summarize(records, config);
      std::cout << "Run: " << run << " Evals: "
                << records[records.size() - 1].best().second << " MES: "
                << summary[MES] << " MAD: " << summary[MAD] << " FAILURES: "
                << summary[FAILURES] << std::endl;
    }
  }
  return records;
}

// Determine the minimum useful population size
int bisection(Random& rand, Configuration& config, evaluation::pointer problem,
              optimize::pointer solver) {
  int runs = config.get<int>("runs");
  float good_enough = config.get<int>("fitness_limit");
  vector<Record> records;
  // initial bounds
  int min = 0;
  int max = 1;
  size_t failed_on = 0;
  bool success;
  // Double the maximum size until a successful run is found
  do {
    // double the bounds each time
    min = max;
    max *= 2;
    config.set("pop_size", max);
    std::cout << "Pop size: " << max << std::endl;
    success = true;
    // Perform multiple runs, stopping as soon as one of the runs failed to find
    // the global optimum
    for (int i = 0; success and i < runs; i++) {
      int problem_number = ((i + failed_on) % runs);
      std::cout << "\tTrying problem: " << problem_number << std::endl;
      auto results = single_run(rand, config, problem, solver, problem_number);
      if (results.best().first < good_enough) {
        failed_on = problem_number;
        success = false;
      }
    }
  } while (!success);

  // Bisect between the min (unsuccessful) and max (successful) until they meet
  int guess;
  while (min + 1 < max) {
    guess = (max + min) / 2;
    std::cout << "Pop size: " << guess << std::endl;
    config.set("pop_size", guess);
    success = true;
    // Perform multiple runs, stopping as soon as one of the runs failed to find
    // the global optimum
    for (int i = 0; success and i < runs; i++) {
      int problem_number = ((i + failed_on) % runs);
      std::cout << "\tTrying problem: " << problem_number << std::endl;
      auto results = single_run(rand, config, problem, solver, problem_number);
      if (results.best().first < good_enough) {
        failed_on = problem_number;
        success = false;
      }
    }
    // update the bound
    if (success) {
      max = guess;
    } else {
      min = guess;
    }
  }
  return max;
}

// Performs a bisection method where a single success decreases predicted population size,
// and a single failure returns to the previously predicted size.
int fast_bisection(Random& rand, Configuration& config,
                   evaluation::pointer problem, optimize::pointer solver) {
  // initial bounds
  int min = 0;
  int max = 1;
  int result = -1;
  // Double the maximum size until a successful run is found
  do {
    // double the bounds each time
    min = max;
    max *= 2;
    // Perform multiple runs, stopping as soon as one of the runs failed to find
    // the global optimum
    result = recurse(rand, config, problem, solver, min, max);
  } while (result == -1);
  return result;
}

// Recursively drills down from a given size until either it finds the correct
// population size, or you discover that max isn't big enough
int recurse(Random& rand, Configuration& config, evaluation::pointer problem,
            optimize::pointer solver, int min, int max) {
  int runs = config.get<int>("runs");
  float good_enough = config.get<float>("fitness_limit");
  vector<Record> records;
  int result;
  for (int i = 0; i < runs; i++) {
    config.set("pop_size", max);
    std::cout << "Pop size: " << max << " Trying problem: " << i << std::endl;
    auto results = single_run(rand, config, problem, solver, i);
    if (results.best().first < good_enough) {
      // fail this size, you have to get bigger
      std::cout << "\tfailed" << std::endl;
      return -1;
    } else {
      int guess = (max + min) / 2;
      if (min != guess) {
        result = recurse(rand, config, problem, solver, min, guess);
        // If you found the correct population size
        if (result != -1) {
          return result;
        } else {
          // You failed, increase the minimum
          min = guess;
        }
      }
    }
  }
  return max;
}
