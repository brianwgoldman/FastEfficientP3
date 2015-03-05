// Brian Goldman

// Implementation of the Hierarchical Bayesian Optimization Algorithm
// and helper classes based on the publication:
// Pelikan, M. and Goldberg, D. (2006). Hierarchical bayesian optimization algorithm. In Scalable
// Optimization via Probabilistic Modeling, volume 33 of Studies in Computational Intelligence, pages
// 63–90. Springer Berlin Heidelberg.

#include "Evaluation.h"
#include "Configuration.h"
#include "Optimizer.h"
#include "Util.h"
#include "HillClimb.h"
#include <algorithm>
#include <fstream>
#include <unordered_set>

// Helper class used to store and calculate fractions built from factorials
class Factorial_Fraction {
 public:
  Factorial_Fraction(size_t n = 0)
      : upper(n + 1, 0),
        lower(n + 1, 0) {
  }
  // Multiply by n!
  void mul_by_factorial(size_t n);
  // Divide by n!
  void div_by_factorial(size_t n);
  // Removes terms which appear in both numerator and denominator
  void simplify();
  // Convert to real value representation, doing smart operations to avoid rounding errors
  operator double() const;
  // Utility for outputting the fraction
  friend std::ostream& operator<<(std::ostream& out, const Factorial_Fraction&);
  // Multiply two fractions
  Factorial_Fraction operator*(const Factorial_Fraction& rhs) const;
  // Divide two fractions
  Factorial_Fraction operator/(const Factorial_Fraction& rhs) const;
 private:
  // Keep track of how often each number appears in the numerator
  vector<size_t> upper;
  // Keep track of how often each number appears in the denominator
  vector<size_t> lower;
  // Resize upper and lower to fit n!
  void safe_size(size_t n);
};

// Encodes a Bayesian Decision Tree for a single gene
struct Bayesian_Tree {
  Bayesian_Tree(size_t i = -1)
      : index(i),
        is_leaf(true),
        left(nullptr),
        right(nullptr),
        counts( { 0, 0 } ) {
  }
  // Needed to handle dynamic memory copies and assignments using the
  // copy and swap idiom
  Bayesian_Tree(const Bayesian_Tree& rhs);
  ~Bayesian_Tree();
  friend void swap(Bayesian_Tree& lhs, Bayesian_Tree& rhs);
  Bayesian_Tree& operator=(Bayesian_Tree rhs);

  // Utility functions for visualizing the tree
  operator string() const;
  friend std::ostream& operator<<(std::ostream& out, const Bayesian_Tree&);

  // Pass this solution through the decision tree until it reaches the
  // appropriate leaf.
  void add_solution(const vector<bool>* const solution);
  // Split this leaf based on new_index
  void split(size_t new_index);
  // Join the subtrees of this internal node to create a single leaf
  void join();

  // Use the decision tree to set "index" in solution
  void set_bit(Random& rand, vector<bool>& solution);

  // Calculate the BDe score of this tree
  Factorial_Fraction bde_fraction();
  // Calculate the BDe score of this leaf after splitting on new_index
  Factorial_Fraction splitless_bde(size_t new_index);

  // The active index this node uses. For internal nodes this is the splitting index,
  // for leaves this is generating index
  size_t index;
  bool is_leaf;

  // Track all solutions which reach this leaf in the decision tree
  vector<const vector<bool>*> solutions;

  // Children in the decision tree
  Bayesian_Tree* left;
  Bayesian_Tree* right;

  // Keep track of how many solutions have 0 and 1 for "index"
  std::array<size_t, 2> counts;
};

// Collection of trees used to model each index in the solution
class Bayesian_Forest {
 public:
  Bayesian_Forest(size_t n);
  // Provide this solution to each tree for modeling purposes
  void add_solution(const vector<bool>* solution);
  // Construct the entire forest based on the added solutions
  void build_forest(Random & rand);
  // Generate a new solution based on the build model
  void generate(Random& rand, vector<bool> & solution);
  // Utility function for displaying the forest
  friend std::ostream& operator<<(std::ostream& out, const Bayesian_Forest&);
 private:
  // Add all potential splits for a tree to a list of options
  void build_options(Bayesian_Tree& tree,
                     vector<std::tuple<double, Bayesian_Tree*, size_t>>&) const;
  // Remove invalid splits from the list of options based on a threshold
  void filter(vector<std::tuple<double, Bayesian_Tree*, size_t>> &, double);
  // Keep track of tree dependences to prevent the creation of cycles
  // such that prev[i] is all trees which i depends on, directly or indirectly
  // and post[i] is all tree which depend on i, directly or indirectly
  vector<std::unordered_set<size_t>> prev, post;
  // The list of trees, one for each index
  vector<Bayesian_Tree> trees;
  // The order in which trees can be evaluated to ensure all dependencies are met
  vector<size_t> ordering;
  // Keep track of how many solutions have been added to the forest
  size_t added;
};

// The Hierarchical Bayesian Optimization Algorithm implementation
class HBOA : public Optimizer {
 public:
  HBOA(Random& _rand, shared_ptr<Evaluator> _evaluator, Configuration& _config);
  // Perform a generation of search, returns false if convergence is reached
  bool iterate() override;
  create_optimizer(HBOA);
 private:
  // Find the restricted tournament replacement competitor for a given solution.
  int rtr_nearest(const vector<bool>& solution);
  // The population of solutions and their qualities
  vector<vector<bool>> solutions;
  vector<float> fitnesses;

  // Tool used to make RTR more efficients
  vector<int> selection;
  size_t rtr_size;
  // Optional hill climbing
  hill_climb::pointer hc;
  size_t completed_generations;
};
