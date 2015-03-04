#include <cassert>
#include <cfloat>
#include <cmath>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <functional>
#include <random>
#include <string>
#include <vector>
#include <unordered_map>

#define MALE 0
#define FEMALE 1

class Agent;

struct Partnership {
  unsigned partner_id;
  double date_started;
};

struct InitialVals {
  std::mt19937 rng;
  double prob_male = 0.5;
  double prob_hiv_pos = 0.1;
  double prob_tst_pos = 0.4;
  double prob_tb_inf = 0.45;
  double prob_tb_sick = 0.01;
  double prob_hetero = 1.0;
  double current_date = 0.0;
  double max_x_coord = 10.0;
  double max_y_coord = 10.0;
  unsigned last_agent = 0;
  unsigned seed = 0;
  bool verbose = false;

  // proportions for distance calc
  double age_factor = 1.0;
  double orientation_factor = 100.0;
  double tightness_factor = 1.0;
  double distance_factor = 0.1;
  double previous_partner_factor = 500;
};

class Agent {
private:
  InitialVals& initial_vals_;
public:
  Agent(InitialVals & initial_vals) : initial_vals_(initial_vals) {
    std::uniform_real_distribution<double> uni;
    std::uniform_real_distribution<double> uni_x(0, initial_vals_.max_x_coord);
    std::uniform_real_distribution<double> uni_y(0, initial_vals_.max_y_coord);
    id = initial_vals_.last_agent++;
    sex = uni(initial_vals_.rng) < initial_vals_.prob_male ? MALE : FEMALE;
    dob = initial_vals_.current_date;
    tightness = uni(initial_vals.rng);
    hiv_pos = uni(initial_vals_.rng) < initial_vals_.prob_hiv_pos ? true : false;
    tst_pos = uni(initial_vals_.rng) < initial_vals_.prob_tst_pos ? true : false;
    tb_inf = uni(initial_vals_.rng) < initial_vals_.prob_tb_inf ? true : false;
    tb_sick = uni(initial_vals_.rng) < initial_vals_.prob_tb_sick ? true : false;
    hetero = uni(initial_vals_.rng) < initial_vals_.prob_hetero ? 1.0 : 0.0;
    x_coord = uni_x(initial_vals_.rng);
    y_coord = uni_y(initial_vals_.rng);
  }
  double euclidean_distance(const Agent &a)
  {
    double x_d = (x_coord - a.x_coord);
    double y_d = (x_coord - a.y_coord);
    return std::sqrt(x_d * x_d + y_d * y_d);
  }

  double distance(const Agent &a, unsigned partner_count = 0)
  {
    double prev_partner;
    double age_diff;
    double orientation_diff;
    double tightness_diff;
    double distance_diff;

    if (count (partners.begin(), partners.end(), &a) > partner_count)
      prev_partner = initial_vals_.previous_partner_factor;
    else
      prev_partner = 0;

    age_diff = initial_vals_.age_factor * (dob - a.dob);
    if (sex == a.sex)
      orientation_diff = initial_vals_.orientation_factor * (hetero + a.hetero);
    else
      orientation_diff = initial_vals_.orientation_factor *
	((1.0 - hetero) + (1.0 - a.hetero));
    tightness_diff = initial_vals_.tightness_factor * (tightness - a.tightness);
    distance_diff = initial_vals_.distance_factor * euclidean_distance(a);

    return fabs(age_diff) + fabs(orientation_diff) +
      fabs(tightness_diff) + fabs(distance_diff) + prev_partner;
  }

  double cluster_value()
  {
    return initial_vals_.age_factor * dob +
      initial_vals_.orientation_factor * hetero +
      initial_vals_.tightness_factor * tightness;
  }

  unsigned id;
  unsigned sex;
  unsigned dob;
  double tightness;
  double hetero;
  bool hiv_pos;
  bool tst_pos;
  bool tb_inf;
  bool tb_sick;
  double x_coord;
  double y_coord;
  double weight;
  std::vector<Agent* > partners;
};


class Simulation {
private:
  bool free_initial_vals_ = false;
  InitialVals *initial_vals_;
public:
  std::vector<Agent *> agents;
  InitialVals initial_vals;
  unsigned iteration = 0;
  std::vector<unsigned> positions;

  Simulation(InitialVals *initial_vals_parm = NULL)
  {
    if (initial_vals_parm) {
      initial_vals = *initial_vals_parm;
    } else {
      initial_vals_ = (new InitialVals());
      initial_vals = *initial_vals_;
      free_initial_vals_ = true;
    }
    initial_vals.rng.seed(initial_vals.seed);
  }

  ~Simulation() {
    for (auto &a : agents) delete a;
    if (free_initial_vals_) delete initial_vals_;
  }

  void init_population(std::size_t size)
  {
    for (std::size_t i = 0; i < size; ++i) {
      Agent *a = new Agent(initial_vals);
      agents.push_back(a);
    }
  }

  void
  print_agent_ids(std::vector<Agent *> &agents,
		  std::string delim = " ",
		  std::string after = "\n")
  {
    for (auto &a : agents)
      std::cout << a->id << delim;
    std::cout << after;
  }

  void reset()
  {
    for (auto &a : agents) a->partners.clear();
    iteration = 1;
  }

  unsigned
  find_partner_rank(Agent *agent)
  {
    unsigned position = 0;
    double d = agent->distance(*agent->partners.back(), 1);
    for (auto & a : agents) {
      if (a != agent && a != agent->partners.back()) {
	double x = agent->distance(*a);
	if (x < d) ++position;
      }
    }
    return position;
  }

  std::vector<Agent *>::iterator
  closest_pair_match(std::vector<Agent *>::iterator from,
		     std::vector<Agent *>::iterator to)
  {
    double smallest_val = DBL_MAX;
    std::vector<Agent *>::iterator closest_agent = to;

    for (auto it = from + 1; it != to; ++it) {
      if ( (*it)->partners.size() < iteration) {
	double distance = (*from)->distance(**it);
	if (distance < smallest_val) {
	  smallest_val = distance;
	  closest_agent = it;
	}
      }
    }
    return closest_agent;
  }

  std::vector<Agent *>::iterator
  closest_pair_match_n(std::vector<Agent *>::iterator from,
		       std::vector<Agent *>::iterator to,
		       unsigned n)
  {
    double smallest_val = DBL_MAX;
    std::vector<Agent *>::iterator closest_agent = to;

    unsigned j = 0;
    for (auto it = from + 1; j < n && it != to; ++it) {
      if ( (*it)->partners.size() < iteration) {
	double distance = (*from)->distance(**it);
	if (distance < smallest_val) {
	  smallest_val = distance;
	  closest_agent = it;
	}
	++j;
      }
    }
    return closest_agent;
  }

  void
  brute_force_match()
  {
    std::shuffle(agents.begin(), agents.end(), initial_vals.rng);
    for (auto it = agents.begin(); it != agents.end(); ++it) {
      if ( (*it)->partners.size() < iteration) {
	auto partner =  closest_pair_match(it, agents.end());
	if (partner != agents.end()) {
	  (*it)->partners.push_back(*partner);
	  (*partner)->partners.push_back(*it);
	}
      }
    }
  }

  void random_match()
  {
    std::shuffle(agents.begin(), agents.end(), initial_vals.rng);
    std::uniform_real_distribution<double> uni;
    for (auto it = agents.begin(); it < agents.end() - 1; it+=2) {
      (*it)->partners.push_back( *(it + 1) );
      (* (it + 1) )->partners.push_back(*it);
    }
  }

  void random_match_n(unsigned neighbors)
  {
    std::uniform_real_distribution<double> uni;
    std::shuffle(agents.begin(), agents.end(), initial_vals.rng);
    for (auto it = agents.begin(); it < agents.end() - 1; ++it) {
      if ( (*it)->partners.size() < iteration) {
	auto last = agents.end() - it < (neighbors + 1) ?
					agents.end() : it + neighbors + 1;
	auto partner = closest_pair_match_n(it, last, neighbors);
	if (partner != last) {
	  (*it)->partners.push_back(*partner);
	  (*partner)->partners.push_back(*it);
	}
      }
    }
  }

  void weighted_shuffle_match(unsigned neighbors)
  {
    std::uniform_real_distribution<double> uni;
    for (auto & agent: agents) agent->weight =
				 agent->cluster_value() * uni(initial_vals.rng);
    std::sort(agents.rbegin(), agents.rend(),
	      [](Agent *a, Agent *b) {return a->weight < b->weight; });
    for (auto it = agents.begin(); it < agents.end() - 1; ++it) {
      if ( (*it)->partners.size() < iteration) {
	auto last = agents.end() - it < (neighbors + 1) ?
					agents.end() : it + neighbors + 1;
	auto partner = closest_pair_match_n(it, last, neighbors);
	if (partner != last) {
	  (*it)->partners.push_back(*partner);
	  (*partner)->partners.push_back(*it);
	}
      }
    }
  }

  void
  cluster_shuffle_match(unsigned clusters,
			unsigned neighbors)
  {
    unsigned cluster_size = agents.size() / clusters;
    for (auto &a : agents) a->weight = a->cluster_value();
    sort(agents.rbegin(), agents.rend(), [](Agent *a, Agent *b)
	 { return a->weight < b->weight; });
    for (unsigned i = 0; i < clusters; ++i) {
      auto first = agents.begin() + i * cluster_size;
      auto last = first + cluster_size;
      if (last > agents.end()) last = agents.end();
      std::shuffle(first, last, initial_vals.rng);
    }
    for (auto it = agents.begin(); it < agents.end() - 1; ++it) {
      if ( (*it)->partners.size() < iteration) {
	auto last = agents.end() - it < (neighbors + 1) ?
					agents.end() : it + neighbors + 1;
	auto partner = closest_pair_match_n(it, last, neighbors);
	if (partner != last) {
	  (*it)->partners.push_back(*partner);
	  (*partner)->partners.push_back(*it);
	}
      }
    }
  }

  double
  calc_avg_match()
  {
    unsigned total = 0;
    unsigned partnerships = 0;
    unsigned samesex = 0;
    //std::vector<unsigned> positions;
    unsigned denom = agents.size();
    positions.clear();

    for (auto & a: agents) {
      if (a->partners.size() < iteration) {
	--denom;
      } else {
	assert(a->partners.back());
	assert(a->partners.back()->partners.back());
	assert(a->partners.back()->partners.back() == a);
	++partnerships;
	if (a->sex == a->partners.back()->sex) ++samesex;
	unsigned position = find_partner_rank(a);
	positions.push_back(position);
	total += position;
      }
    }
    if (initial_vals.verbose) {
      std::cout << "Number of partnerships: " << partnerships / 2 << std::endl;
      std::cout << "Number same sex partnerships: " << samesex / 2 << std::endl;
      std::cout << "Number of agents without partners: "
		<< agents.size() - denom << std::endl;
    }
    return (double) total / denom;
  }
};

/* Time measuring function taken from:
   http://codereview.stackexchange.com/questions/48872/measuring-execution-time-in-c
*/

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
  template<typename F, typename ...Args>
  static typename TimeT::rep execution(F func, Args&&... args)
  {
    auto start = std::chrono::system_clock::now();

    // Now call the function with all the parameters you need.
    func(std::forward<Args>(args)...);

    auto duration = std::chrono::duration_cast< TimeT>
      (std::chrono::system_clock::now() - start);

    return duration.count();
  }
};

template<class InputIterator> double
median (InputIterator  from, InputIterator to, bool sorted = false)
{
  std::size_t l = std::distance(from, to);
  double median;
  if (sorted == false)
    std::sort(from, to);
  if (l % 2 == 0) {
    median = double (from[l/2] + from[l/2 - 1]) / 2.0;
  } else {
    median = from[l/2];
  }
  return median;
}

void
stats(Simulation &s, const char *description,
      std::function<void(void)> func, unsigned iterations = 1, unsigned run = 0)
{
  double mean;

  s.reset();
  for (unsigned i = 0; i < iterations; ++i, ++s.iteration) {
    auto start = std::chrono::system_clock::now();
    func();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::system_clock::now() - start);

    std::cout << run << ", " << description << ", " << i << ", algorithm time, "
	      << duration.count() << std::endl;

    start = std::chrono::system_clock::now();
    mean = s.calc_avg_match();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::system_clock::now() - start);

    std::cout << run << ", " << description << ", " << i << ", ranker time, "
	      << duration.count() << std::endl;

    std::cout << run << ", " << description  << ", " << i << ", mean, "
	      << mean << std::endl;
    auto mdn = median(s.positions.begin(), s.positions.end());
    std::cout << run << ", " << description  << ", " << i << ", median, "
	      << mdn << std::endl;
  }
}

//////////////////////////////
/* COMMAND LINE PROCESSING */

/*
  Simple command line processing functions taken from:
  http://stackoverflow.com/questions/865668/parse-command-line-arguments
*/

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    {
      return *itr;
    }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

/* END COMMAND LINE PROCESSING */
/////////////////////////////////

void run_tests(std::size_t population = 16,
	       unsigned clusters = 4,
	       unsigned neighbors = 2,
	       unsigned iterations = 1,
	       unsigned seed = 0,
	       unsigned runs = 1,
	       bool verbose = true)
{
  InitialVals initial_vals;

  initial_vals.seed = seed;
  initial_vals.verbose = verbose;

  Simulation s(&initial_vals);

  s.init_population(population);

  if (verbose) {
    std::cout << "Population: " << population << std::endl;
    std::cout << "Clusters: " << clusters << std::endl;
    std::cout << "Neighbors: " << neighbors << std::endl;
  }

  if (verbose) s.print_agent_ids(s.agents);

  for (unsigned i = 0; i < runs; ++i) {
    stats(s, "Random match", [&](){s.random_match();}, iterations, i);
    stats(s, "Random match n", [&](){s.random_match_n(neighbors);},
	  iterations, i);
    stats(s, "Weighted shuffle", [&](){s.weighted_shuffle_match(neighbors); },
	  iterations, i);
    stats(s, "Cluster shuffle", [&](){
	s.cluster_shuffle_match(clusters, neighbors);
      }, iterations, i);
    stats(s, "Brute force", [&](){s.brute_force_match();}, iterations, i);
  }
}

int main(int argc, char *argv[])
{
  unsigned population = 16, clusters = 4, neighbors = 2,
    seed = 0, iterations = 1, runs = 1;
  bool verbose = false;
  char *population_str = getCmdOption(argv, argv + argc, "-p");
  char *neighbors_str = getCmdOption(argv, argv + argc, "-n");
  char *clusters_str = getCmdOption(argv, argv + argc, "-c");
  char *iterations_str = getCmdOption(argv, argv + argc, "-i");
  char *seed_str = getCmdOption(argv, argv + argc, "-s");
  char *runs_str = getCmdOption(argv, argv + argc, "-r");

  if (population_str) {
    population = atoi(population_str);
    if (population > 63) {
      clusters = neighbors = std::round(std::log2(population));
    }
  }
  if (clusters_str)
    clusters = atoi(clusters_str);
  if (neighbors_str)
    neighbors = atoi(neighbors_str);
  if (iterations_str)
    iterations = atoi(iterations_str);
  if (seed_str)
    seed = atoi(seed_str);
  if (runs_str)
    runs = atoi(runs_str);
  if(cmdOptionExists(argv, argv+argc, "-v"))
    verbose = true;

  run_tests(population, clusters, neighbors, iterations, seed, runs, verbose);
}
