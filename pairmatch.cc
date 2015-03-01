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
  double prob_hetero = 0.92;
  double current_date = 0.0;
  double max_x_coord = 10.0;
  double max_y_coord = 10.0;
  unsigned last_agent = 0;
  bool verbose = false;

  // proportions for distance calc
  double age_factor = 1.0;
  double orientation_factor = 1.0;
  double tightness_factor = 1.0;
  double distance_factor = 0.1;
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
    num_partners = 0;
    partner = NULL;
  }
  double euclidean_distance(const Agent &a)
  {
    double x_d = (x_coord - a.x_coord);
    double y_d = (x_coord - a.y_coord);
    return std::sqrt(x_d * x_d + y_d * y_d);
  }

  double distance(const Agent &a)
  {
    double age_diff;
    double orientation_diff;
    double tightness_diff;
    double distance_diff;

    age_diff = initial_vals_.age_factor * (dob - a.dob);
    if (sex == a.sex)
      orientation_diff = initial_vals_.orientation_factor * (hetero + a.hetero);
    else
      orientation_diff = initial_vals_.orientation_factor *
	((1.0 - hetero) + (1.0 - a.hetero));
    tightness_diff = initial_vals_.tightness_factor * (tightness - a.tightness);
    distance_diff = initial_vals_.distance_factor * euclidean_distance(a);

    return std::fabs(age_diff + orientation_diff +
		     tightness_diff + distance_diff);
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
  unsigned num_partners;

  Agent *partner;
};


class Simulation {
public:
  std::vector<Agent *> agents;
  InitialVals initial_vals;

  Simulation(InitialVals *initial_vals_parm = NULL)
  {
    if (initial_vals_parm) initial_vals = *initial_vals_parm;
    initial_vals.rng.seed(0);
  }

  ~Simulation() {
    for (auto &a : agents) delete a;
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

  void reset_partners(std::vector<Agent *> &agents)
  {
    for (auto &a : agents) a->partner = NULL;
  }

  unsigned
  find_partner_rank(std::vector<Agent *> &agents,
		    Agent *agent)
  {
    unsigned position = 0;
    double d = agent->distance(*agent->partner);
    for (auto & a : agents) {
      if (a != agent && a != agent->partner) {
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
      if ( (*it)->partner == NULL) {
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
      if ( (*it)->partner == NULL) {
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
  brute_force_match(std::vector<Agent *> &agents)
  {
    std::shuffle(agents.begin(), agents.end(), initial_vals.rng);
    for (auto it = agents.begin(); it != agents.end(); ++it) {
      if ( (*it)->partner == NULL) {
	auto partner =  closest_pair_match(it, agents.end());
	if (partner != agents.end()) {
	  (*it)->partner = *partner;
	  (*partner)->partner = *it;
	}
      }
    }
  }

  void random_match(std::vector<Agent *> &agents)
  {
    std::shuffle(agents.begin(), agents.end(), initial_vals.rng);
    std::uniform_real_distribution<double> uni;
    for (auto it = agents.begin(); it < agents.end() - 1; it+=2) {
      (*it)->partner = *(it + 1);
      (* (it + 1) )->partner = *it;
    }
  }

  void weighted_shuffle_match(std::vector<Agent *> &agents, unsigned n)
  {
    std::uniform_real_distribution<double> uni;
    for (auto & agent: agents) agent->weight =
				       agent->cluster_value() * uni(initial_vals.rng);
    std::sort(agents.rbegin(), agents.rend(),
	      [](Agent *a, Agent *b) {return a->weight < b->weight; });
    for (auto it = agents.begin(); it < agents.end() - 1; ++it) {
      if ( (*it)->partner == NULL) {
	auto last = agents.end() - it < (n + 1) ? agents.end() : it + n + 1;
	auto partner = closest_pair_match_n(it, last, n);
	if (partner != last) {
	  (*it)->partner = *partner;
	  (*partner)->partner = *it;
	  if (initial_vals.verbose)
	    std::cout << "Match (agent, partner):" << (*it)->id << " "
		      << (*partner)->id << std::endl;
	}
      }
    }
  }

  void
  cluster_shuffle_match(std::vector<Agent *> &agents,
			unsigned clusters,
			unsigned n)
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
      if ( (*it)->partner == NULL) {
	auto last = agents.end() - it < (n + 1) ? agents.end() : it + n + 1;
	auto partner = closest_pair_match_n(it, last, n);
	if (partner != last) {
	  (*it)->partner = *partner;
	  (*partner)->partner = *it;
	  if (initial_vals.verbose)
	    std::cout << "Match (agent, partner): " << (*it)->id << " "
		      << (*partner)->id << std::endl;
	}
      }
    }
  }

  double
  calc_avg_match(std::vector<Agent *> &agents)
  {
    unsigned total = 0;
    unsigned denom = agents.size();
    for (auto & a: agents) {
      if (a->partner == NULL) {
	--denom;
      } else {
	assert(a->partner->partner == a);
	unsigned position = find_partner_rank(agents, a);
	if (initial_vals.verbose)
	  std::cout << "Rank (agent, partner, rank): "
		    << a->id << " " << a->partner->id
		    << " " << position << std::endl;
	total += position;
      }
    }
    if (initial_vals.verbose)
      std::cout << "Number of agents without partners: "
		<< agents.size() - denom << std::endl;
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
	       bool verbose = true)
{
  double ws_avg, rn_avg, cs_avg, bf_avg;
  InitialVals initial_vals;

  initial_vals.verbose = verbose;

  Simulation s(&initial_vals);

  s.init_population(population);

  std::cout << "Population: " << population << std::endl;
  std::cout << "Clusters: " << clusters << std::endl;
  std::cout << "Neighbors: " << neighbors << std::endl;

  if (verbose) s.print_agent_ids(s.agents);

  std::cout << "Algorithm: Weighted shuffle" << std::endl;

  std::cout << "Time taken: " << measure<>::execution( [&]() {
      s.weighted_shuffle_match(s.agents, neighbors);
    }) << std::endl;

  std::cout << "Time taken by ranker: " << measure<>::execution( [&]() {
      ws_avg = s.calc_avg_match(s.agents);
    }) << std::endl;

  std::cout << "Algorithm: Random match" << std::endl;
  s.reset_partners(s.agents);

  std::cout << "Time taken: " << measure<>::execution( [&]() {
      s.random_match(s.agents);
    }) << std::endl;

  std::cout << "Time taken by ranker: " << measure<>::execution( [&]() {
      rn_avg = s.calc_avg_match(s.agents);
    }) << std::endl;

  std::cout << "Algorithm: Cluster shuffle" << std::endl;
  s.reset_partners(s.agents);
  std::cout << "Time taken: " << measure<>::execution( [&]() {
      s.cluster_shuffle_match(s.agents, clusters, neighbors);
    }) << std::endl;

  std::cout << "Time taken by ranker: " << measure<>::execution( [&]() {
      cs_avg = s.calc_avg_match(s.agents);
    }) << std::endl;


  std::cout << "Algorithm: Brute force" << std::endl;
  s.reset_partners(s.agents);
  std::cout << "Time taken: " << measure<>::execution( [&]() {
      s.brute_force_match(s.agents);
    }) << std::endl;

  std::cout << "Time taken by ranker: " << measure<>::execution( [&]() {
      bf_avg = s.calc_avg_match(s.agents);
    }) << std::endl;

  std::cout << "Random match average: " << rn_avg << std::endl;
  std::cout << "Weighted shuffle average: " << ws_avg << std::endl;
  std::cout << "Cluster shuffle average: " << cs_avg << std::endl;
  std::cout << "Brute force average: " << bf_avg << std::endl;
}

int main(int argc, char *argv[])
{
  unsigned population = 16, clusters = 4, neighbors = 2;
  bool verbose = false;
  char *population_str = getCmdOption(argv, argv + argc, "-p");
  char *neighbors_str = getCmdOption(argv, argv + argc, "-n");
  char *clusters_str = getCmdOption(argv, argv + argc, "-c");

  if (population_str)
    population = atoi(population_str);
  if (clusters_str)
    clusters = atoi(clusters_str);
  if (neighbors_str)
    neighbors = atoi(neighbors_str);
  if(cmdOptionExists(argv, argv+argc, "-v")) verbose = true;

  run_tests(population, clusters, neighbors, verbose);
}
