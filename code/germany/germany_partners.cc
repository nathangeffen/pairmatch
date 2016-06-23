#include <algorithm>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <ctime>
#include <unordered_map>
#include <vector>

#include "stats.hh"

std::mt19937 rng;


#define MALE 0
#define FEMALE 1
#define HETEROSEXUAL 0
#define HOMOSEXUAL 1
#define WANTS_TO_BE_SINGLE 0
#define WANTS_TO_BE_PARTNERED 1
#define NUM_RELS 2
#define NUM_AGES 100
#define NUM_SEXES 2
#define NUM_ORIENTATIONS 2

class Agent;

typedef std::unordered_map<const char *, double> ParameterMap;
typedef std::vector<Agent *> AgentVector;


struct Agent {
  unsigned id;
  unsigned age;
  unsigned sex;
  unsigned sexor;
  unsigned rel;
  unsigned page;
  unsigned psex;
  unsigned psexor;
  unsigned desired_age;
  Agent* partner = NULL;
  double weight;

  // A very simple clustering function to get people of same age and
  // sexual orientation nearer each other when clustering
  double cluster_value()
  {
    return (page / 100.0 + age / 100.0) / 2.0 + sexor;
  };
};


void destroy_agents(AgentVector& agents)
{
  for (auto &agent: agents)
    delete agent;
}



static void read_in_csv_file(const char* filename,
			     AgentVector& agents)
{
  std::ifstream infile;
  infile.open (filename, std::ifstream::in);
  if (infile.fail()) {
    std::cerr << "Error opening " << filename << std::endl;
    exit(1);
  }
  unsigned line = 0;
  while (infile)
  {
    try {
      ++line;
      std::string s;
      if (!std::getline( infile, s )) break;
      std::istringstream ss(s);

      if (line == 1) continue;
      Agent* agent = new Agent();
      if (!std::getline( ss, s, ',' )) break;
      agent->id = stol(s);
      if (!std::getline( ss, s, ',' )) break;
      agent->age = stol(s);
      if (!std::getline( ss, s, ',' )) break;
      agent->sex = stol(s);
      if (!std::getline( ss, s, ',' )) break;
      agent->sexor = stol(s);
      if (!std::getline( ss, s, ',' )) break;
      agent->rel = stol(s);
      // pid - skip
      if (!std::getline( ss, s, ',' )) break;
      // page
      if (!std::getline( ss, s, ',' )) break;
      agent->page = stol(s);
      // psex
      if (!std::getline( ss, s, ',' )) break;
      agent->psex = stol(s);
      // psexor
      if (!std::getline( ss, s, ',' )) break;
      agent->psexor = stol(s);
      if (!std::getline( ss, s)) break;
      if (s.size() == 0 || s[0] != 'N') {
	agent->desired_age = stol(s);
	if (agent->desired_age != agent->page) {
	  std::cerr << "Non-matching ages at line " << line << std::endl;
	  exit(1);
	}
      } else {
	agent->desired_age = 0;
      }
      agents.push_back(agent);
    } catch (std::exception &e) {
      std::cerr << "Exception at line " << line << " with message: " << e.what()
		<< std::endl;
      exit(1);
    }
  }

  if (!infile.eof())
  {
    std::cerr << "Error occurred on line: " << line << std::endl;
    exit(1);
  }
}

void clear_partners(AgentVector& agents)
{
  for (auto& agent : agents) {
    agent->partner = NULL;
  }
}

/* Auxiliary matching functions */

static bool check_for_partial_match(const Agent *a, const Agent *b)
{
  if (a->rel == WANTS_TO_BE_PARTNERED &&
      b->rel == WANTS_TO_BE_PARTNERED &&
      a->page == b->age &&
      a->age == b->page &&
      a->sexor == b->sexor) {
    if (a->sexor == HETEROSEXUAL && a->sex != b->sex) {
      return true;
    } else if (a->sexor == HOMOSEXUAL && a->sex == b->sex) {
      return true;
    }
  }
  return false;
}

static bool check_for_match(const Agent *a, const Agent *b)
{
  if (a->partner == false &&
      b->partner == false) {
    return check_for_partial_match(a, b);
  }
  return false;
}


static void make_partner(Agent *a, Agent *b)
{
  a->partner = b;
  b->partner = a;
}

static void check_for_errors(const AgentVector& agents)
{
  unsigned errors = 0;
  for (auto & agent: agents) {
    if (agent->rel == WANTS_TO_BE_SINGLE && agent->partner != NULL) {
      ++errors;
      std::cout << "ERROR - Single in relationship: " << agent->id << std::endl;
      continue;
    }
    if (agent->partner && !check_for_partial_match(agent, agent->partner)) {
      ++errors;
      std::cout << "ERROR - Mismatched partnership: "
		<< agent->id << " - " << agent->partner->id << std::endl;
    }
  }
  if (errors)
    std::cout << "Errors: " << errors << std::endl;
}

static AgentVector::iterator
find_closest_match(AgentVector::iterator from,
			AgentVector::iterator to)
{
  for (auto it = from + 1; it != to; ++it) {
    if (check_for_match(*from, *it)) {
      return it;
    }
  }
  return to;
}

/* Finds first available partner in k nearest neighbours */

void find_partner(AgentVector& agents, const unsigned k)
{
  for (auto it = agents.begin(); it < agents.end() - 1; ++it) {
    if ( (*it)->rel == WANTS_TO_BE_PARTNERED) {
      auto last = (agents.end() - it) < (k + 1) ?
				      agents.end() : it + k + 1;
      auto partner = find_closest_match(it, last);
      if (partner != last) {
	make_partner(*it, *partner);
      }
    }
  }
}

/* Matching Algorithms */


/* Reference partner matching algorithm: Random match. */

/* This algorithm is hopeless but fast. */
void random_match(AgentVector& agents, const ParameterMap& parameters)
{
  std::shuffle(agents.begin(), agents.end(), rng);
  for (size_t i = 0; i < agents.size() - 1; ++i) {
    if (check_for_match(agents[i], agents[i + 1])) {
      make_partner(agents[i], agents[i + 1]);
    }
  }
}


/* Select first matching partner from k nearest neighbours */
void random_k_match(AgentVector& agents, const ParameterMap& parameters)
{
  unsigned k = (unsigned) parameters.at("neighbors");
  std::shuffle(agents.begin(), agents.end(), rng);
  find_partner(agents, k);
}

/* Select first matching partner from weighted shuffle of k nearest
   neightbours.
*/
void weighted_shuffle_match(AgentVector& agents, const ParameterMap& parameters)
{
  std::uniform_real_distribution<double> uni;
  unsigned k = parameters.at("neighbors");
  for (auto & agent: agents)
    agent->weight = agent->cluster_value() * uni(rng);
  std::sort(agents.begin(), agents.end(),
	    [](Agent *a, Agent *b) {return a->weight < b->weight; });
  find_partner(agents, k);
}


/* Select first matching partner from weighted cluster of k nearest neighbours.
 */

void
cluster_shuffle_match(AgentVector& agents, const ParameterMap& parameters)
{
  unsigned k = parameters.at("neighbors");
  unsigned clusters = parameters.at("clusters");
  unsigned cluster_size = agents.size() / clusters;
  for (auto &a : agents) a->weight = a->cluster_value();
  sort(agents.begin(), agents.end(), [](Agent *a, Agent *b)
       { return a->weight < b->weight; });
  for (unsigned i = 0; i < clusters; ++i) {
    auto first = agents.begin() + i * cluster_size;
    auto last = first + cluster_size;
    if (last > agents.end()) last = agents.end();
    std::shuffle(first, last, rng);
  }
  find_partner(agents, k);
}

/*

 */

struct Table {
  size_t start;
  size_t entries;
};

void
distribution_match(AgentVector& agents, const ParameterMap& parameters)
{
  // We are going to match at most k neighbours
  unsigned k = parameters.at("neighbors");

  // Shuffle the agents - O(n)
  std::shuffle(agents.begin(), agents.end(), rng);
  // Make a copy of the agent **POINTERS** - O(n)
  AgentVector& copy_agents = agents;
  // Sort the agent pointers on age, sex, sexor, desired_age O(n log n)
  // Perhaps this can be replaced with distribution_sort but
  // the effort isn't worth it because sorting takes 1/30th of the execution time
  // on a 1m record input file with k=3000.
  std::sort(copy_agents.begin(), copy_agents.end(),
	    [](Agent *a, Agent *b) {
	      if (a->rel < b->rel) return true;
	      if (b->rel < a->rel) return false;
	      if (a->age < b->age) return true;
	      if (b->age < a->age) return false;
	      if (a->sex < b->sex) return true;
	      if (b->sex < a->sex) return false;
	      if (a->sexor < b->sexor) return true;
	      if (b->sexor < a->sexor) return false;
	      if (a->page < b->page) return true;
	      if (b->page < a->page) return false;
	      return false;
	    });
  // We need a distribution table. Initialization O(1)
  Table table[NUM_RELS][NUM_AGES][NUM_SEXES][NUM_ORIENTATIONS] = {0, 0};

  // Populate the table indices - O(n)
  for(auto & agent: copy_agents)
    ++table[agent->rel][agent->age][agent->sex][agent->sexor].entries;
  size_t last_index = 0;
  // This is constant time despite the four loops - O(1)
  for (size_t i = 0; i < NUM_RELS; ++i) {
    for (size_t j = 0; j < NUM_AGES; ++j) {
      for (size_t k = 0; k < NUM_SEXES; ++k) {
	for (size_t l = 0; l < NUM_ORIENTATIONS; ++l) {
	  table[i][j][k][l].start = last_index;
	  last_index += table[i][j][k][l].entries;
	}
      }
    }
  }
  // Now match - O(n)
  for (auto & agent: agents) {
    // Calculate the start and end indices
    size_t rel = agent->rel;
    size_t age = agent->page;
    size_t sex = (agent->sexor == HOMOSEXUAL) ? agent->sex : (!agent->sex);
    size_t sexor = agent->sexor;
    size_t start_index = table[rel][age][sex][sexor].start;
    size_t last_index =
      std::min(start_index + table[rel][age][sex][sexor].entries, agents.size());
    for (size_t i = start_index; (i < last_index) && (i < (start_index + k));
	 ++i) {
      if (agent == copy_agents[i]) continue; // Can't partner yourself
      if (check_for_match(agent, copy_agents[i])) {
	make_partner(agent, copy_agents[i]);
	// Swap with the last entry in this part of the table
	std::swap(copy_agents[i], copy_agents[last_index - 1]);
	--table[rel][age][sex][sexor].entries;
	break;
      }
    }
  }
}

/* Reporting */

static void report_csv(const char *filename, const AgentVector& agents)
{
  std::ofstream fout;
  fout.open (filename, std::ofstream::out);

  fout << "id, age, sex, sexor, rel, pid, page, psex, psexor, page1"
	    << std::endl;
  for (auto & agent: agents) {
    fout << agent->id << ", "
	 << agent->age << ", "
	 << agent->sex << ", "
	 << agent->sexor << ", "
	 << agent->rel << ", "
	 << ((agent->partner) ? agent->partner->id : 0) << ", "
	 << ((agent->partner) ? round(agent->partner->age) : 0) << ", "
	 << ((agent->partner) ? round(agent->partner->sex) : 0) << ", "
	 << ((agent->partner) ? round(agent->partner->sexor) : 0) << ", "
	 << agent->desired_age << std::endl;
  }
  fout.close();
}


static void report_stats(const char *report_name,
			 const unsigned report_num,
			 const AgentVector& agents,
			 const float time_taken,
			 const ParameterMap& parameters)
{
  unsigned agents_to_be_matched = 0;
  unsigned agents_successfully_matched = 0;
  for (auto & agent: agents) {
    if (agent->rel == WANTS_TO_BE_PARTNERED) {
      ++agents_to_be_matched;
      if (agent->partner) {
	++agents_successfully_matched;
      }
    }
  }
  double success_rate = (double)
    agents_successfully_matched / agents_to_be_matched ;
  std::cout << report_name << ", " << report_num << ", "
	    << parameters.at("neighbors")  << ", "
	    << parameters.at("clusters")  << ", "
	    << agents_to_be_matched << ", "
	    << agents_successfully_matched << ", "
	    << success_rate << ", "
	    << time_taken << std::endl;
}

void run_tests(ParameterMap& parameters,
	       const std::string algorithms_to_run,
	       const unsigned number_of_runs,
	       const char* input_file,
	       const bool csv_output)
{
  AgentVector agents;
  std::function<void(AgentVector &, const ParameterMap)> algorithm;
  const char* algorithm_name;
  std::string csv_file_name;
  unsigned neighbors = parameters["neighbors"];
  unsigned clusters = parameters["clusters"];
  read_in_csv_file(input_file, agents);
  std::cout << "alg,run,k,c,tomatch,success,rate,time"
	    << std::endl;
  for (auto & c: algorithms_to_run) {
    switch(c) {
    case 'R':
      algorithm = random_match;
      algorithm_name = "Random";
      break;
    case 'K':
      algorithm = random_k_match;
      algorithm_name = "RandomK";
      break;
    case 'W':
      algorithm = weighted_shuffle_match;
      algorithm_name = "Weighted";
      break;
    case 'C':
      algorithm = cluster_shuffle_match;
      algorithm_name = "Cluster";
      break;
    case 'D':
      algorithm = distribution_match;
      algorithm_name = "Distribution";
      break;
    default:
      std::cerr << "Unrecognised algorithm code: " << c << std::endl;
      exit(1);
    }
    // Resent k size and clusters before new algorithm
    parameters["neighbors"] = neighbors;
    parameters["clusters"] = clusters;
    for (unsigned i = 0; i < number_of_runs; ++i) {
      clear_partners(agents);
      clock_t t = clock();
      algorithm(agents, parameters);
      t = clock() - t;
      float time_taken = (float) t / CLOCKS_PER_SEC;
      check_for_errors(agents);
      if (csv_output) {
	std::string csv_filename =
	  std::string(algorithm_name) + std::string("_") +
	  std::to_string(i) + std::string(".csv");
	report_csv(csv_filename.c_str(), agents);
      }
      report_stats(algorithm_name, i, agents, time_taken, parameters);
      if (parameters.at("varyk") > 0.0)
	parameters["neighbors"] += parameters.at("varyk");
      if (parameters.at("varyc") > 0.0)
	parameters["clusters"] += parameters.at("varyc");
    }
  }
  destroy_agents(agents);
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


int main(int argc, char *argv[])
{
  ParameterMap parameters;
  unsigned seed = 23;
  unsigned runs = 1;
  std::string algorithms = "RKWCD";
  bool csv_output = false;
  unsigned varyk = 0;
  unsigned varyc = 0;
  unsigned k = 350;
  unsigned clusters = 100;

  if (cmdOptionExists(argv, argv + argc, "-h")) {
    std::cout << argv[0] << " options, where options are:\n"
	      << " [-i input file name]\n"
	      << "    (if this option isn't specified, agent_list.csv is used)\n"
	      << " [-s random seed integer]\n"
	      << " [-a ([R][K][W][C][D])+] \n"
	      << "    where:\n"
	      << "    R = Random matching\n"
	      << "    K = Random k matching\n"
	      << "    W = Weighted shuffle matching\n"
	      << "    C = Cluster shuffle matching\n"
	      << "    D = Distribution matching\n"
	      << " [-r number of times each algorithm must be run]\n"
	      << " [-k number of neighbours to search for some algorithms]\n"
	      << " [-c number of clusters for Cluster shuffle algorithm]\n"
	      << " [-vk number to add to k on each run]\n"
	      << " [-vc number to add to clusters on each run]\n"
	      << " [-o] writes CSV output files containing all agents"
	      << std::endl;
    exit(1);
  }

  const char *input_file_str = getCmdOption(argv, argv + argc, "-i");
  const char *seed_str = getCmdOption(argv, argv + argc, "-s");
  const char *algorithms_str = getCmdOption(argv, argv + argc, "-a");
  const char *runs_str = getCmdOption(argv, argv + argc, "-r");
  const char *neighbors_str = getCmdOption(argv, argv + argc, "-k");
  const char *varyk_str = getCmdOption(argv, argv + argc, "-vk");
  const char *varyc_str = getCmdOption(argv, argv + argc, "-vc");
  const char *clusters_str = getCmdOption(argv, argv + argc, "-c");
  if (cmdOptionExists(argv, argv + argc, "-o")) csv_output = true;

  if (!input_file_str) input_file_str = "input_agents.csv";
  if (seed_str) seed = std::stol(std::string(seed_str));
  if (algorithms_str) algorithms = std::string(algorithms_str);
  if (runs_str) runs = std::stol(std::string(runs_str));
  if (neighbors_str)
    parameters["neighbors"] =  (double) std::stol(std::string(neighbors_str));
  else
    parameters["neighbors"] =  (double) k;
  if (varyk_str)
    parameters["varyk"] =  (double) std::stol(std::string(varyk_str));
  else
    parameters["varyk"] =  (double) varyk;
  if (varyc_str)
    parameters["varyc"] =  (double) std::stol(std::string(varyc_str));
  else
    parameters["varyc"] =  (double) varyc;

  if (clusters_str)
    parameters["clusters"] =  (double) std::stol(std::string(clusters_str));
  else
    parameters["clusters"] =  (double) clusters;

  rng.seed(seed);
  run_tests(parameters, algorithms, runs, input_file_str, csv_output);
}
