#include "Evaluation.h"
#include "Configuration.h"
#include "Optimizer.h"
#include "Util.h"
#include "HillClimb.h"
#include <algorithm>
#include <fstream>
#include <unordered_set>

class Factorizer {
 public:
  Factorizer(size_t n);
  void operator()(vector<size_t>& buckets) const;
 private:
  vector<size_t> factors;
};

class Factorial_Fraction {
 public:
  Factorial_Fraction(size_t n=0): upper(n+1, 0), lower(n+1, 0) {}
  void mul_by_factorial(size_t n);
  void div_by_factorial(size_t n);
  void simplify(const Factorizer& factorizer);
  operator double() const;
  friend std::ostream& operator<<(std::ostream& out, const Factorial_Fraction&);
  Factorial_Fraction operator*(const Factorial_Fraction& rhs) const;
  Factorial_Fraction operator/(const Factorial_Fraction& rhs) const;
 private:
  vector<size_t> upper;
  vector<size_t> lower;
  void safe_size(size_t n);
};

struct Bayesian_Tree {
  Bayesian_Tree(size_t i=-1) : index(i), is_leaf(true), left(nullptr), right(nullptr), counts({0, 0}) {};
  Bayesian_Tree(const Bayesian_Tree& rhs);
  ~Bayesian_Tree();
  friend void swap(Bayesian_Tree& lhs, Bayesian_Tree& rhs);
  Bayesian_Tree& operator=(Bayesian_Tree rhs);
  operator string() const;
  friend std::ostream& operator<<(std::ostream& out, const Bayesian_Tree&);

  void add_solution(const vector<bool>& solution);
  void split(size_t new_index);
  void join();

  void set_bit(Random& rand, vector<bool>& solution);

  Factorial_Fraction bde_fraction();
  Factorial_Fraction splitless_bde(size_t new_index);
  size_t index;
  bool is_leaf;
  vector<vector<bool>> solutions;
  Bayesian_Tree* left;
  Bayesian_Tree* right;
  std::array<size_t, 2> counts;
};

class Bayesian_Forest {
 public:
  Bayesian_Forest(size_t n);
  void add_solution(const vector<bool>& solution);
  void build_forest();
  void generate(Random& rand, vector<bool> & solution);
  friend std::ostream& operator<<(std::ostream& out, const Bayesian_Forest&);
 private:
  void build_options(Bayesian_Tree& tree, const Factorizer&, vector<std::tuple<double, Bayesian_Tree*, size_t>>&) const;
  void filter(vector<std::tuple<double, Bayesian_Tree*, size_t>> &, double);
  vector<std::unordered_set<size_t>> prev, post;
  vector<Bayesian_Tree> trees;
  vector<size_t> ordering;
  size_t added;
};

class HBOA : public Optimizer {
 public:
  HBOA(Random& _rand, Evaluator& _evaluator, Configuration& _config);
  bool iterate() override;
  create_optimizer(HBOA);
 private:
  int rtr_nearest(const vector<bool>& solution);
  vector<vector<bool>> solutions;
  vector<float> fitnesses;
  vector<int> selection;
  size_t rtr_size;
};
