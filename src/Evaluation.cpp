// Brian Goldman

// Each evaluation object defines the "evaluate" function for how
// to convert the vector of bool solution into a floating point fitness.

#include "Evaluation.h"
using namespace std;

// Template specialization to convert a string to a function pointer
// Used to allow configurable problem choice
template<>
evaluation::pointer Configuration::get(const string key) {
  return evaluation::lookup[get<string>(key)];
}

// Simple example for how to do evaluation
float OneMax::evaluate(const vector<bool> & solution) {
  float sum = 0;
  for (const bool & bit : solution) {
    sum += bit;
  }
  // Convert to percentage of total
  return sum / solution.size();
}

// Iterate over traps, sum partial scores
float DeceptiveTrap::evaluate(const vector<bool> & solution) {
  int partial;
  int total = 0;

  for (size_t i = 0; i < solution.size(); i += trap_size) {
    // Find the number of bits set in each trap
    partial = 0;
    for (size_t index = i; index < i + trap_size; index++) {
      partial += solution[index];
    }

    // Make it deceptive
    if (partial < trap_size) {
      partial = trap_size - partial - 1;
    }
    total += partial;
  }

  // Convert to percentage of total
  return float_round(float(total) / solution.size(), precision);
}

// Iterate over traps, sum partial scores
float DeceptiveStepTrap::evaluate(const vector<bool> & solution) {
  int partial;
  int total = 0;

  int trap_maximum = (offset + trap_size) / step_size;
  for (size_t i = 0; i < solution.size(); i += trap_size) {
    // Find the number of bits set in each trap
    partial = 0;
    for (size_t index = i; index < i + trap_size; index++) {
      partial += solution[index];
    }

    // Make it deceptive
    if (partial < trap_size) {
      partial = trap_size - partial - 1;
    }

    // Convert it to stepwise
    total += (offset + partial) / step_size;
  }

  // Convert to percentage of total
  float fitness = (float(total) * trap_size) / (solution.size() * trap_maximum);
  return float_round(fitness, precision);
}

// Iterate over the 3 bit traps, sum up their scores
float Deceptive3::evaluate(const vector<bool> & solution) {
  int partial;
  float total = 0;
  for (size_t i = 0; i < solution.size(); i += 3) {
    // Find the number of bits set in each trap
    partial = 0;
    for (size_t index = i; index < i + 3; index++) {
      partial += solution[index];
    }
    // score the trap
    total += scores[partial];
  }

  // Convert to percentage of total
  return float_round((3 * total) / solution.size(), precision);
}

// Iterate over traps, sum partial scores
float Bipolar::evaluate(const vector<bool> & solution) {
  int partial;
  double total = 0;
  int k = trap_size / 2;
  for (size_t i = 0; i < solution.size(); i += trap_size) {
    // Find the number of bits set in each trap
    partial = 0;
    for (size_t index = i; index < i + trap_size; index++) {
      partial += solution[index];
    }
    int score = abs(partial - k);
    if (score == k) {
      total += 1;
    } else if (score == 0) {
      total += 0.9;
    } else if (score == 1) {
      total += 0.8;
    }
  }

  // Convert to percentage of total
  return float_round(total * trap_size / solution.size(), precision);
}

float JumpK::evaluate(const vector<bool> & solution) {
  // Count the number of bits set to one
  size_t count = 0;
  for (const bool bit : solution) {
    count += bit;
  }
  float score;
  if (count == solution.size() or count <= solution.size() - k) {
    score = count + k;
  } else {
    // Push solutions away from all ones if they are closer than k.
    score = solution.size() - count;
  }
  return score / (solution.size() + k);
}

// Attempts to load the problem file, otherwise constructs a new problem
// solves it, and saves it to a problem file
NearestNeighborNK::NearestNeighborNK(Configuration& config, int run_number) {
  k = config.get<int>("k");
  length = config.get<int>("length");
  precision = config.get<int>("precision");
  table.resize(length, vector<float>(2 << k, 0));
  int rng_seed = config.get<int>("problem_seed") + run_number;

  // Build up the filename where this problem is stored
  string filename = config.get<string>("problem_folder");
  filename += +"NearestNeighborNK_";
  filename += config.get<string>("length") + "_";
  filename += config.get<string>("k") + "_";
  filename += to_string(rng_seed) + ".txt";

  ifstream in(filename);
  // If this problem has been tried before
  if (in) {
    // Read in information about the global minimum
    in >> minimum;
    worst.resize(length);
    string temp;
    in >> temp;
    for (int i = 0; i < length; i++) {
      worst[i] = temp[i] == '1';
    }

    // Read in information about the global maximum
    in >> maximum;
    best.resize(length);
    in >> temp;
    for (int i = 0; i < length; i++) {
      best[i] = temp[i] == '1';
    }

    // Read in the fitness table
    for (auto& row : table) {
      for (auto& entry : row) {
        in >> entry;
      }
    }
    in.close();
  } else {
    // Generate the new problem
    Random rand(rng_seed);

    // Generate the table
    auto generator = std::uniform_real_distribution<>(0, 1);
    for (auto& row : table) {
      for (auto& entry : row) {
        entry = make_filable(generator(rand));
      }
    }

    if (config.get<string>("method") == "old") {
      // Find its minimum and maximum
      minimum = make_filable(solve(worst, false));
      maximum = make_filable(solve(best, true));
    } else {
      minimum = make_filable(hammer_solve(worst, false));
      maximum = make_filable(hammer_solve(best, true));
    }
    auto& out = cout;

    out << minimum << " ";
    print(worst, out);

    out << maximum << " ";
    print(best, out);

    throw "STOP";

    for (auto& row : table) {
      for (auto& entry : row) {
        out << entry << " ";
      }
      out << endl;
    }
    //out.close();
  }
}

// Used in finding the minimum / maximum of the generated problem.
// Finds the fitness contribution of setting a particular section in a specific way
float NearestNeighborNK::chunk_fitness(trimap& known, size_t chunk_index,
                                       size_t a, size_t b) {
  // If we know the fitness, return it
  const auto& first = known.find(chunk_index);
  if (first != known.end()) {
    const auto& second = first->second.find(a);
    if (second != first->second.end()) {
      const auto& third = second->second.find(b);
      if (third != second->second.end()) {
        return third->second;
      }
    }
  }

  // Calculate the fitness
  float fitness = 0;
  size_t mask = (2 << k) - 1;

  size_t combined = (a << k) | b;
  combined = (combined << 1) | (a & 1);

  for (size_t g = 0; g < k; g++) {
    size_t value = (combined >> (k - g)) & mask;
    fitness += table[chunk_index * k + g][value];
  }

  // Record the chuck fitness for later
  known[chunk_index][a][b] = fitness;
  return fitness;
}

// Converts a number into a series of bits
void NearestNeighborNK::int_into_bit(size_t src, vector<bool>& dest) {
  for (size_t i = 1; i <= k; i++) {
    dest.push_back((src >> (k - i)) & 1);
  }
}

float NearestNeighborNK::hammer_solve(vector<bool>& solution, bool maximize) {
  // TODO Handle case where N < dependencies
  // TODO Ensure this works for k=0
  const size_t dependencies = 2 * k;
  const size_t patterns = 1 << dependencies;
  const size_t function_mask = (1 << (k+1)) - 1;
  // Pads 1 past N to make the logic simpler
  // TODO pretty sure you only need 2 hammer functions: previous and under construction
  vector<vector<float>> hammer_functions(length+1, vector<float>(patterns, 0));
  vector<vector<char>> best_bit_choice(length, vector<char>(patterns));
  const char TIE = 2;
  // Initialize the table with k last functions already compressed
  for (size_t f=length - k; f < length; f++) {
    // Used to strip off bits not used by this function
    const size_t shift_size = length - f - 1;
    for (size_t pattern=0; pattern < patterns; pattern++) {
      // Leaves only the k+1 bits used by this function
      size_t relevant_bits = (pattern >> shift_size) & function_mask;
      hammer_functions[length][pattern] += table[f][relevant_bits];
    }
  }
  const size_t powk = 1 << k;
  const size_t hammer_function_mask = (1 << dependencies) - 1;
  // "n" is the bit position we are currently trying to remove.
  // You could also say at the end of this loop the problem size will be n.
  for (size_t n=length-1; n >= dependencies; n--) {
    for (size_t left=0; left < powk; left++) {
      for (size_t right=0; right < powk; right++) {
        // quality of setting bit "n" to zero and one.
        float zero_quality = 0, one_quality = 0;
        // The pattern the dependent bits take
        const size_t zero_pattern = (left << (k+1)) | right;
        const size_t one_pattern = zero_pattern | powk;
        // bit "n" only depends on two functions at this point:
        // function_{n-k} and the hammer_function[n+1]
        zero_quality += table[n - k][zero_pattern >> k];
        one_quality += table[n - k][one_pattern >> k];
        zero_quality += hammer_functions[n + 1][zero_pattern & hammer_function_mask];
        one_quality += hammer_functions[n + 1][one_pattern & hammer_function_mask];
        // put the results into the table
        const size_t hammer_pattern = (left << k) | right;
        if ((maximize and zero_quality > one_quality) or
            (not maximize and zero_quality < one_quality)) {
          // Zero was better
          hammer_functions[n][hammer_pattern] = zero_quality;
          best_bit_choice[n][hammer_pattern] = 0;
        } else if (zero_quality == one_quality) {
          // tie
          hammer_functions[n][hammer_pattern] = zero_quality;
          best_bit_choice[n][hammer_pattern] = TIE;
          cout << "A tie happened" << endl;
        } else {
          // One was better
          hammer_functions[n][hammer_pattern] = one_quality;
          best_bit_choice[n][hammer_pattern] = 1;
        }
      }
    }
  }

  // The hammer_functions table is now wide enough to score all possible patterns
  // so we can just copy in the remaining functions
  const size_t k_bits = (1 << k) - 1;
  for (size_t f=0; f < k; f++) {
    // Used to strip off bits not used by this function
    size_t shift_size = k - f - 1;
    for (size_t pattern=0; pattern < patterns; pattern++) {
      // Leaves only the k+1 bits used by this function
      size_t relevant_bits = (pattern >> shift_size) & function_mask;
      // TODO Simplify this, perhaps with 2 loops
      // The handedness of this flips
      size_t hammer_pattern = ((pattern >> k) & k_bits) | ((pattern & k_bits) << k);
      hammer_functions[dependencies][hammer_pattern] += table[f][relevant_bits];
    }
  }

  // Whichever hammer_functions entry for "dependencies" has the highest quality
  // is now the best pattern
  const auto& final_column = hammer_functions[dependencies];
  size_t best_pattern = 0;
  for (size_t pattern=0; pattern < patterns; pattern++) {
    if ((maximize and final_column[pattern] > final_column[best_pattern]) or
        (not maximize and final_column[pattern] < final_column[best_pattern])) {
      // TODO Probably need a list of patterns to deal with ties
      best_pattern = pattern;
    }
  }

  // TODO pretty sure this is broken
  // Extract a global optimum
  solution.resize(length);
  // First "dependencies" bits comes from "best_pattern";
  for (size_t i=0; i < dependencies; i++) {
    solution[i] = (best_pattern >> (dependencies - i - 1)) & 1;
  }
  size_t current_pattern = best_pattern;
  for (size_t i= dependencies; i < length; i++) {
    auto best_bit = best_bit_choice[i][current_pattern];
    if (best_bit == TIE) {
      // TODO Handle this better
      best_bit = 0;
    }
    solution[i] = best_bit;
    // shift out old bit information, add in new bit
    current_pattern = ((current_pattern << 1) & hammer_function_mask) | best_bit;
  }
  return final_column[best_pattern];
}

// Find either the minimum or maximum (depending on the last argument) of the problem
// See the following paper for full explanation:
// "The computational complexity of N-K fitness functions"
// by A. H. Wright, R. K. Thompson, and J. Zhang
float NearestNeighborNK::solve(vector<bool>& solution, bool maximize) {
  size_t numbers = 1 << k;
  trimap known;
  std::unordered_map<size_t,
      std::unordered_map<size_t, std::unordered_map<size_t, size_t> > > partial;
  float current;
  // Iteratively build up partial values and their fitnesses
  for (size_t n = length / k - 1; n > 1; n--) {
    // stores the fitness of settings
    std::unordered_map<size_t, std::unordered_map<size_t, float> > utility;
    // stores the bit values that achieve best partial fitnesses
    std::unordered_map<size_t, std::unordered_map<size_t, size_t> > value;
    // Loop through the pieces
    for (size_t left = 0; left < numbers; left++) {
      for (size_t right = 0; right < numbers; right++) {
        utility[left][right] = -1;
        if (not maximize) {
          utility[left][right] = 2 * length;
        }
        for (size_t middle = 0; middle < numbers; middle++) {
          current = chunk_fitness(known, n - 1, left, middle);
          current += chunk_fitness(known, n, middle, right);
          if ((maximize and utility[left][right] < current)
              or (not maximize and utility[left][right] > current)) {
            utility[left][right] = current;
            value[left][right] = middle;
          }
        }
      }
    }

    // Extract information
    for (size_t left = 0; left < numbers; left++) {
      for (size_t right = 0; right < numbers; right++) {
        known[n - 1][left][right] = utility[left][right];
        partial[n][left][right] = value[left][right];
      }
    }
  }

  // Initialize fitness to be worse than worst possible
  float fitness = -1;
  if (not maximize) {
    fitness = 2 * length;
  }

  // Find the best left and right values
  size_t best_left = 0;
  size_t best_right = 0;
  for (size_t left = 0; left < numbers; left++) {
    for (size_t right = 0; right < numbers; right++) {
      current = chunk_fitness(known, 0, left, right);
      current += chunk_fitness(known, 1, right, left);
      if ((maximize and fitness < current)
          or (not maximize and fitness > current)) {
        fitness = current;
        best_left = left;
        best_right = right;
      }
    }
  }

  // Recreate the optimal string
  solution.clear();
  solution.reserve(length);
  int_into_bit(best_left, solution);
  int_into_bit(best_right, solution);
  size_t last = best_right;
  for (size_t i = 2; i < length / k; i++) {
    last = partial[i][last][best_left];
    int_into_bit(last, solution);
  }

  return fitness;
}

// Use the table to evaluate the quality of the solution
float NearestNeighborNK::evaluate(const vector<bool> & solution) {
  float total = 0;
  for (size_t i = 0; i < solution.size(); i++) {
    // Construct the integer represented by this subset of the solution
    size_t index = 0;
    for (size_t neighbor = i; neighbor <= i + k; neighbor++) {
      index = (index << 1) | solution[neighbor % length];
    }
    // Look up that bit combination in the table
    total += table[i][index];
  }

  // Convert to percentage of total
  float fitness = (total - minimum) / (maximum - minimum);
  // Ensures the best fitness actually gets 1.0
  return float_round(fitness, precision);
}

// Simple evaluation, stop counting when a 0 is found
float LeadingOnes::evaluate(const vector<bool> & solution) {
  for (size_t i = 0; i < solution.size(); i++) {
    if (not solution[i]) {
      return float_round(float(i) / solution.size(), precision);
    }
  }
  // only reach here if all bits are set to 1
  return 1;
}

// Score the solution at iteratively larger block sizes
float HIFF::evaluate(const vector<bool> & solution) {
  // Data structure used to represent values as 0, 1, or -1 (Null)
  int * level = new int[solution.size()];
  int level_length = solution.size();

  // Initialize the level to the current solution
  for (size_t i = 0; i < solution.size(); i++) {
    level[i] = solution[i];
  }
  int power = 1;
  int next_length = level_length >> 1;
  int total = 0;
  int maximum = 0;

  // Keep going while the next level actual has bits in it
  while (next_length > 0) {
    int * next_level = new int[next_length];
    // Construct the next level using the current level
    for (int i = 0; i + 1 < level_length; i += 2) {
      if (level[i] == level[i + 1] and level[i] != -1) {
        // Score points for a correct setting at this level
        total += power;
        next_level[i >> 1] = level[i];
      } else {
        next_level[i >> 1] = -1;
      }
      // Keep track of the maximum possible score
      maximum += power;
    }
    delete[] level;
    level = next_level;
    level_length = next_length;
    next_length = level_length >> 1;
    power <<= 1;
  }
  delete[] level;

  // Convert to percentage of total
  return float(total) / maximum;
}

// Generates the new problem each time its needed, based on
// the problem see and run number
MAXSAT::MAXSAT(Configuration& config, int run_number) {
  size_t length = config.get<int>("length");
  precision = config.get<int>("precision");
  clauses.resize(
      float_round(config.get<float>("clause_ratio") * length, precision));
  signs.resize(clauses.size());

  int rng_seed = config.get<int>("problem_seed") + run_number;
  Random rand(rng_seed);

  // Create the random target solution
  vector<bool> solution = rand_vector(rand, length);

  // Data structure used to select random variables to include in a clause
  vector<int> options(length);
  std::iota(options.begin(), options.end(), 0);

  // Distributions for selecting 3 variables at random
  std::uniform_int_distribution<> dist[] = {
      std::uniform_int_distribution<>(0, length - 1),
      std::uniform_int_distribution<>(1, length - 1),
      std::uniform_int_distribution<>(2, length - 1)
  };

  auto sign_select = std::uniform_int_distribution<>(0, sign_options.size() - 1);

  for (size_t i = 0; i < clauses.size(); i++) {
    // Choose a random sign setting for this clause
    int select = sign_select(rand);

    // For each term in this clause
    for (int c = 0; c < 3; c++) {
      // Choose a random variable to be a term in the clause
      std::swap(options[c], options[dist[c](rand)]);
      clauses[i][c] = options[c];
      // Set the sign to agree / disagree with the solution based on chosen
      // sign option
      signs[i][c] = sign_options[select][c] == solution[options[c]];
    }
  }
}

// Count how many clauses evaluate to true
float MAXSAT::evaluate(const vector<bool> & solution) {
  int total = 0;
  for (size_t i = 0; i < clauses.size(); i++) {
    for (size_t c = 0; c < 3; c++) {
      // if the literal is true, move to the next clause
      if (solution[clauses[i][c]] == signs[i][c]) {
        total++;
        break;
      }
    }
  }

  // Convert to percentage of total
  return float_round(float(total) / clauses.size(), precision);
}

// Read in the problem from a file and set up the evaluation table
IsingSpinGlass::IsingSpinGlass(Configuration& config, int run_number)
    : length(config.get<int>("length")),
      precision(config.get<int>("precision")) {
  int rng_seed = config.get<int>("problem_seed") + run_number;

  // Build up the filename where this problem is stored
  string filename = config.get<string>("problem_folder");
  filename += +"IsingSpinGlass_";
  filename += config.get<string>("ising_type") + "_";
  filename += config.get<string>("length") + "_";
  filename += to_string(rng_seed) + ".txt";
  ifstream in(filename);
  if (!in) {
    throw invalid_argument(
        "IsingSpinGlass data file does not exist: " + filename);
  }
  in >> min_energy;
  string solution_string;
  in >> solution_string;
  int number_of_spins;
  in >> number_of_spins;
  spins.resize(number_of_spins);
  for (auto& spin : spins) {
    for (auto& part : spin) {
      in >> part;
    }
  }
  in.close();
  span = number_of_spins - min_energy;

  // Sanity check
  vector<bool> solution(length);
  for (int i = 0; i < length; i++) {
    solution[i] = solution_string[i] == '1';
  }
  if (evaluate(solution) != 1) {
    throw invalid_argument(
        "IsingSpinGlass data file has inconsistent data: " + filename);
  }
}

// Evaluate using the read in spins
float IsingSpinGlass::evaluate(const vector<bool>& solution) {
  int energy = 0;
  for (const auto& spin : spins) {
    energy -= (bit_to_sign[solution[spin[0]]] * spin[2]
        * bit_to_sign[solution[spin[1]]]);
  }

  // Convert to percentage of total
  return float_round(1 - (energy - min_energy) / span, precision);
}

// Sets up a table of fitness function values based on possible gray
// value settings
Rastrigin::Rastrigin(Configuration& config, int run_number)
    : precision(config.get<int>("precision")),
      converter(
          BinaryToFloat(config.get<int>("bits_per_float"), -5.12, 5.12,
                        precision)) {
  worst = 0;

  // Build a table for all possible function values
  for (const auto& x : converter.possible()) {
    function[x] = 10 + x * x - 10 * cos(2 * PI * x);
    if (worst < function[x]) {
      worst = function[x];
    }
  }
}

// Converts each group of bits into a grey value, then look up the fitness
float Rastrigin::evaluate(const vector<bool>& solution) {
  auto it = solution.begin();
  float x, total = 0;
  int n = 0;

  // Use the converter to turn bits in the solution into float values
  while (it != solution.end()) {
    x = converter.convert(it);
    total += function[x];
    n++;
  }

  // Convert to percentage of total
  total /= (n * worst);
  return float_round(1 - total, precision);
}

// Call the external script file to perform evaluation
float External::evaluate(const vector<bool>& solution) {
  // write the solution to the output file
  ofstream output(out_file);
  print(solution, output);
  output.close();

  // calls the script file
  int error_code = system(script_file.c_str());
  if (error_code) {
    throw invalid_argument("Script file returned non zero success");
  }

  // read the fitness from the input file
  ifstream input(in_file);
  float fitness;
  input >> fitness;
  return fitness;
}
