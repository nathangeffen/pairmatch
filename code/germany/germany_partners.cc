#include <algorithm>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <ctime>
#include <cstdio>
#include <unordered_map>
#include <vector>

#include "CSV_Parser.hh"
#include "sample.hh"
#include "stats.hh"

std::mt19937 rng;

#define FEMALE 0
#define MALE 1
#define HETEROSEXUAL 0
#define HOMOSEXUAL 1
#define WANTS_TO_BE_SINGLE 0
#define WANTS_TO_BE_PARTNERED 1
#define NUM_RELS 2
#define NUM_AGES 101
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

double str_to_num(std::string& s)
{
  boost::replace_all(s, ",", ".");
  boost::replace_all(s, "NA", "0");
  return stod(s);
}

static void read_agent_file(const char* filename,
			    AgentVector& agents)
{
    CSVParser agents_csv(filename, ",", true);
    for (auto& row: agents_csv.string_rows) {
      Agent* agent = new Agent();
      // "id","age","sex","sexor","rel","pid","page","psex","psexor"
      agent->id = str_to_num(row[0]);
      agent->age = str_to_num(row[1]);
      agent->sex = str_to_num(row[2]);
      agent->sexor = str_to_num(row[3]);
      agent->rel = str_to_num(row[4]);
      // skip row[5] pid
      agent->page = str_to_num(row[6]);
      agent->psex = str_to_num(row[7]);
      agent->psexor = str_to_num(row[8]);
      agents.push_back(agent);
    }
}

void clear_partners(AgentVector& agents)
{
  for (auto& agent : agents) {
    agent->partner = NULL;
  }
}

/* Auxiliary matching functions */

void print_vector(const std::vector<double> vec)
{
  for (auto& d: vec)
    std::cout << d << " ";
  std::cout << std::endl;
}

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
  if (a->partner == NULL &&
      b->partner == NULL) {
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
    if (agent->partner && !agent->partner->partner) {
      ++errors;
      std::cout << "ERROR - One way match: "
		<< agent->id << " - " << agent->partner->id << std::endl;
    }
    if (agent->partner && agent->partner && agent->partner->partner != agent) {
      ++errors;
      std::cout << "ERROR - Cuckolded: "
		<< agent->id << " - " << agent->partner->id << std::endl;
    }
  }
  if (errors)
    std::cout << "Errors: " << errors << std::endl;
}

static AgentVector::iterator
find_closest_match(AgentVector::iterator from,
		   AgentVector::iterator to,
		   unsigned& num_comparisons)
{
  for (auto it = from + 1; it != to; ++it) {
    ++num_comparisons;
    if (check_for_match(*from, *it)) {
      return it;
    }
  }
  return to;
}

/* Finds first available partner in k nearest neighbours */

double find_partners(AgentVector& agents, const unsigned k)
{
  unsigned total_k = 0;
  for (auto it = agents.begin(); it < agents.end() - 1; ++it) {
    if ( (*it)->rel == WANTS_TO_BE_PARTNERED) {
      auto last = (agents.end() - it) < (k + 1) ?
				      agents.end() : it + k + 1;
      auto partner = find_closest_match(it, last, total_k);
      if (partner != last) {
	make_partner(*it, *partner);
      }
    }
  }
  return (double) total_k / agents.size();
}

/* Matching Algorithms */


/* Supporting functions for Stefan's initial population creation algorithm. */



std::vector<double> get_col(const DblMatrix& matrix, unsigned col)
{
  std::vector<double> output;
  for (unsigned i = 0; i < matrix.size(); ++i)
    output.push_back(matrix[i][col]);
  return output;
}

std::vector<double> sub_vector(double d,
			       const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d - v[i]);
  return output;
}

std::vector<double> mult_vector(double d,
				const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d * v[i]);
  return output;
}

std::vector<double> mult_vectors(const std::vector<double>& v1,
				 const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] * v2[i]);
  return output;
}


std::vector<double> add_vectors(const std::vector<double>& v1,
				const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] + v2[i]);
  return output;
}

double sum_vector(const std::vector<double>& dbls)
{
  double total = 0.0;
  for (auto d: dbls)
    total += d;
  return total;
}

bool check_gt_0(const std::vector<double>& dbls)
{
  for (auto& d: dbls) {
    if (d > 0.0)
      return true;
  }
  return false;
}

double calc_number_singles(DblMatrix data, unsigned X)
{
  auto ageshare = get_col(data, 1);
  auto femratio = get_col(data, 2);
  auto relm_share = get_col(data, 3);
  auto relw_share = get_col(data, 4);

  auto t1 = sub_vector(1, relw_share);
  auto t2 = mult_vectors(femratio, t1);
  auto t3 = mult_vectors(t2, ageshare);

  auto t4 = sub_vector(1, femratio);;
  auto t5 = sub_vector(1, relm_share);
  auto t6 = mult_vectors(t4, t5);
  auto t7 = mult_vectors(t6, ageshare);

  auto t8 = add_vectors(t3, t7);
  auto t9 = sum_vector(t8);

  double num_agents = X * t9;
  return round(num_agents);
}

void create_singles(AgentVector& agents,
		    const DblMatrix& data,
		    const unsigned S,
		    const std::vector<double>& ageRange,
		    const std::vector<double>& ageShare,
		    const std::vector<double>& femRatio,
		    const std::vector<double>& wswRate,
		    const std::vector<double>& msmRate)
{
  std::uniform_real_distribution<double> uni;
  Sample sample_ageshare(ageShare, &rng);
  for(unsigned i = 0; i < S; ++i) {
    Agent *agent = new Agent();
    // ID
    agent->id = i + 1;
    // Age
    unsigned age = sample_ageshare();
    agent->age = age;
    // Sex
    unsigned sex = uni(rng) < femRatio[age] ? FEMALE : MALE;
    agent->sex = sex;
    // Sexor
    unsigned sexor;
    if (sex == FEMALE) {
      sexor = uni(rng) < wswRate[age] ? HOMOSEXUAL : HETEROSEXUAL;
    } else {
      sexor = uni(rng) < msmRate[age] ? HOMOSEXUAL : HETEROSEXUAL;
    }
    agent->sexor = sexor;
    agent->rel = 0;
    agent->partner = NULL;
    agent->page = 0;
    agent->psex = 0;
    agent->psexor = 0;
    agents.push_back(agent);
  }
}

void create_partners_ACP(AgentVector& agents,
			 const DblMatrix& data,
			 const unsigned fromAgent,
			 const unsigned toAgent,
			 const std::vector<double>& ageRange,
			 const std::vector<double>& ageShare,
			 const std::vector<double>& femRatio,
			 const std::vector<double>& wswRate,
			 const std::vector<double>& msmRate,
			 const DblMatrix& matWW,
			 const DblMatrix& matMW,
			 const DblMatrix& matWM,
			 const DblMatrix& matMM)
{
  std::uniform_real_distribution<double> uni;
  Sample sample_ageshare(ageShare, &rng);
  vector<Sample> sample_matWW(matWW[0].size());
  vector<Sample> sample_matMW(matMW[0].size());
  vector<Sample> sample_matWM(matWM[0].size());
  vector<Sample> sample_matMM(matMM[0].size());
  std::vector<double> placeholder(matMM.size(), 0.000001);

  for (unsigned i = 0; i < matWW[0].size(); ++i) {
    sample_matWW[i].init(get_col(matWW,i), &rng);
    sample_matWM[i].init(get_col(matWM,i), &rng);
    sample_matMW[i].init(get_col(matMW,i), &rng);
    sample_matMM[i].init(get_col(matMM,i), &rng);
  }
  for (unsigned i = fromAgent; i + 1 < toAgent; i+=2) {
    Agent *agent = new Agent();
    // ID
    agent->id = i + 1;
    // Age
    unsigned age = sample_ageshare() + 12;
    agent->age = age;
    // Sex
    unsigned sex = uni(rng) < femRatio[age - 12] ? FEMALE : MALE;
    agent->sex = sex;
    // Sexor
    unsigned sexor;
    if (sex == FEMALE) {
      sexor = uni(rng) < wswRate[age - 12] ? HOMOSEXUAL : HETEROSEXUAL;
    } else {
      sexor = uni(rng) < msmRate[age - 12] ? HOMOSEXUAL : HETEROSEXUAL;
    }
    agent->sexor = sexor;
    // Relationship
    agent->rel = 1;
    Agent* partner = new Agent();
    // ID of partner
    agent->partner = partner;
    // Partner ID
    partner->id = i + 2;
    // Orientation = partner's orientation
    partner->sexor = sexor;
    // Partner sex
    if (sexor == HETEROSEXUAL) {
      if (sex == MALE)
	partner->sex = FEMALE;
      else
	partner->sex = MALE;
    } else {
      partner->sex = sex;
    }
    // Partner age
    if (sex == FEMALE && sexor == HOMOSEXUAL) {
      partner->age  = sample_matWW[age - 12]() + 12;
    } else if (sex == FEMALE && sexor == HETEROSEXUAL) {
      partner->age = sample_matWM[age - 12]() + 12;
    } else if (sex == MALE && sexor == HETEROSEXUAL) {
      partner->age = sample_matMW[age - 12]() + 12;
    } else {
      partner->age = sample_matMM[age - 12]() + 12;
    }
    // Partner in relationship
    partner->rel = 1;
    // Partner's partner
    partner->partner = agent;
    // Preferred age of partner
    agent->page = partner->age;
    partner->page = agent->age;
    // Preferred partner sex
    agent->psex = partner->sex;
    partner->psex = agent->sex;
    // partner sexor
    agent->psexor = partner->sexor;
    partner->psexor = agent->sexor;
    agents.push_back(agent);
    agents.push_back(partner);
  }
}

/* Stefan's initial population creation algorithm. */

void initial_pop(AgentVector& agents, ParameterMap& parameters)
{
  CSVParser data_csv("data.csv", ";", true);
  DblMatrix data = data_csv.convert_all_entries_to_doubles();
  CSVParser singles_csv("singles.csv", ",", true);
  DblMatrix singles = singles_csv.convert_all_entries_to_doubles();
  CSVParser partners_csv("partnersACP.csv", ",", true);
  DblMatrix partners = partners_csv.convert_all_entries_to_doubles();
  CSVParser mm_csv("mat.msm.csv", ";", false);
  DblMatrix mm = mm_csv.convert_all_entries_to_doubles();
  CSVParser ww_csv("mat.wsw.csv", ";", false);
  DblMatrix ww = ww_csv.convert_all_entries_to_doubles();
  CSVParser mw_csv("mat.msw.csv", ";", false);
  DblMatrix mw = mw_csv.convert_all_entries_to_doubles();
  CSVParser wm_csv("mat.wsm.csv", ";", false);
  DblMatrix wm = wm_csv.convert_all_entries_to_doubles();

  unsigned X = parameters.at("agents");
  unsigned S = calc_number_singles(data, X);


  std::vector<double> ageRange;
  for (unsigned i = 12; i <= 100; ++i)
    ageRange.push_back(i);
  auto ageShare = get_col(singles, 1);
  auto femRatio = get_col(singles, 2);
  auto msmRate = get_col(singles, 3);
  auto wswRate = get_col(singles, 4);
  create_singles(agents, data, S, ageRange, ageShare,
		 femRatio, wswRate, msmRate);
  ageShare = femRatio = msmRate = wswRate = {};
  ageShare = get_col(partners, 1);
  femRatio = get_col(partners, 2);
  msmRate = get_col(partners, 3);
  wswRate = get_col(partners, 4);
  create_partners_ACP(agents, data, S, X, ageRange, ageShare, femRatio,
  		      wswRate, msmRate, ww, mw, wm, mm);
  parameters["avg_k"] = 0.0;
}

/* Reference partner matching algorithm: Random match. */

/* This algorithm is hopeless but fast. */
void random_match(AgentVector& agents, ParameterMap& parameters)
{
  std::shuffle(agents.begin(), agents.end(), rng);
  for (size_t i = 0; i < agents.size() - 1; ++i) {
    if (agents[i]->rel == WANTS_TO_BE_PARTNERED &&
	check_for_match(agents[i], agents[i + 1])) {
      make_partner(agents[i], agents[i + 1]);
    }
  }
  parameters["avg_k"] = 1.0;
}


/* Select first matching partner from k nearest neighbours */
void random_k_match(AgentVector& agents, ParameterMap& parameters)
{
  unsigned k = (unsigned) parameters.at("neighbors");
  std::shuffle(agents.begin(), agents.end(), rng);
  parameters["avg_k"] = find_partners(agents, k);
}

/* Select first matching partner from weighted shuffle of k nearest
   neightbours.
*/
void weighted_shuffle_match(AgentVector& agents, ParameterMap& parameters)
{
  std::uniform_real_distribution<double> uni;
  unsigned k = parameters.at("neighbors");
  for (auto & agent: agents)
    agent->weight = agent->cluster_value() * uni(rng);
  std::sort(agents.begin(), agents.end(),
	    [](Agent *a, Agent *b) {return a->weight < b->weight; });
  parameters["avg_k"] = find_partners(agents, k);
}


/* Select first matching partner from weighted cluster of k nearest neighbours.
 */

void
cluster_shuffle_match(AgentVector& agents, ParameterMap& parameters)
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
  parameters["avg_k"] = find_partners(agents, k);
}

/*
  Table needed by Distribution Match
 */

struct Table {
  size_t start;
  size_t entries;
};

/*
  Sort according to distribution, then use distributional knowledge to locate
  partners.
*/

void
distribution_match(AgentVector& agents, ParameterMap& parameters)
{
  // We are going to match at most k neighbours
  unsigned comparisons = 0;
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
    if (agent->partner || agent->rel == WANTS_TO_BE_SINGLE)
      continue;
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
      if (agent == copy_agents[i])
	continue; // Ignore if partnered and can't partner yourself
      ++comparisons;
      if (check_for_match(agent, copy_agents[i])) {
	make_partner(agent, copy_agents[i]);
	// Swap with the last entry in this part of the table
	std::swap(copy_agents[i], copy_agents[last_index - 1]);
	--table[rel][age][sex][sexor].entries;
	break;
      }
    }
  }
  parameters["avg_k"] = (double) comparisons / agents.size();
}

/* Reporting */

static void report_csv(const char *filename, const AgentVector& agents)
{
  FILE* fout = fopen(filename, "w");
  fprintf(fout, "id, age, sex, sexor, rel, pid, page, psex, psexor\n");
  for (auto& agent: agents) {
    fprintf(fout, "%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
	    agent->id, agent->age, agent->sex, agent->sexor, agent->rel,
	    ((agent->partner) ? agent->partner->id : 0),
	    ((agent->partner) ? agent->partner->age : 0),
	    ((agent->partner) ? agent->partner->sex : 0),
	    ((agent->partner) ? agent->partner->sexor : 0));
  }
  fclose(fout);
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
  std::cout << report_name << "," << report_num << ","
	    << parameters.at("neighbors")  << ","
	    << parameters.at("clusters")  << ","
	    << agents_to_be_matched << ","
	    << agents_successfully_matched << ","
	    << success_rate << ","
	    << parameters.at("avg_k") << ","
	    << time_taken << std::endl;
}

void run_tests(ParameterMap& parameters,
	       const std::string algorithms_to_run,
	       const unsigned number_of_runs,
	       const char* input_file,
	       const bool csv_output)
{
  AgentVector agents;
  unsigned j = 0;
  std::function<void(AgentVector &, ParameterMap &)> algorithm;
  const char* algorithm_name;
  std::string csv_file_name;
  unsigned neighbors = parameters.at("neighbors");
  unsigned clusters = parameters.at("clusters");
  unsigned varyt = parameters.at("varyt");
  if (parameters.at("read_agents") == 1.0) {
    read_agent_file(input_file, agents);
    shuffle(agents.begin(), agents.end(), rng);
  }
  std::cout << "alg,run,k,c,tomatch,success,rate,avgk,time"
	    << std::endl;
  for (auto & c: algorithms_to_run) {
    switch(c) {
    case 'I':
      destroy_agents(agents);
      agents = {};
      algorithm = initial_pop;
      algorithm_name = "Initialpop";
      break;
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
    // Reset k size and clusters before new algorithm
    parameters["neighbors"] = neighbors;
    parameters["clusters"] = clusters;
    for (unsigned i = 0; i < number_of_runs; ++i) {
      clear_partners(agents);
      clock_t t = clock();
      algorithm(agents, parameters);
      t = clock() - t;
      float time_taken = (float) t / CLOCKS_PER_SEC;
      check_for_errors(agents);
      report_stats(algorithm_name, i, agents, time_taken, parameters);
      if (csv_output) {
	std::string csv_filename =
	  std::string(algorithm_name) + std::string("_") +
	  std::to_string(j) + std::string("_") +
	  std::to_string(i) + std::string(".csv");
	report_csv(csv_filename.c_str(), agents);
      }
      if (parameters.at("varyk") > 0.0)
	parameters["neighbors"] += parameters.at("varyk");
      if (varyt == 0 || (i + 1) % varyt == 0) {
	if (parameters.at("varyc") > 0.0) {
	  parameters["clusters"] += parameters.at("varyc");
	  if (varyt > 0) {
	    parameters["neighbors"] = neighbors;
	  }
	}
      }
    }
    ++j;
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
  unsigned varyt = 0;
  unsigned k = 350;
  unsigned clusters = 100;
  unsigned num_agents = 100;

  if (cmdOptionExists(argv, argv + argc, "-h")) {
    std::cout << argv[0] << " options, where options are:\n"
	      << " [-i input file name]\n"
	      << "    (if this isn't specified, input_agents.csv is used)\n"
	      << " [-s random seed integer]\n"
	      << " [-a [n,]([I][R][K][W][C][D])+] \n"
	      << "    where:\n"
	      << "    n = how many times to repeat the subsequent string\n"
	      << "    I = Initial population creation matching\n"
	      << "    R = Random matching\n"
	      << "    K = Random k matching\n"
	      << "    W = Weighted shuffle matching\n"
	      << "    C = Cluster shuffle matching\n"
	      << "    D = Distribution matching\n"
	      << " [-r number of times each algorithm must be run]\n"
	      << " [-k number of neighbours to search for some algorithms]\n"
	      << "     (default if left out is 350)\n"
	      << " [-c number of clusters for Cluster shuffle algorithm]\n"
	      << "     (default if left out is 100)\n"
	      << " [-vk number to add to k on each run]\n"
	      << " [-vc number to add to clusters on each run]\n"
	      << " [-vt vary k every run and clusters every specified run]\n"
	      << " [-o] writes CSV output files containing all agents\n"
	      << " [-h] displays this message\n";
    exit(1);
  }

  const char *input_file_str = getCmdOption(argv, argv + argc, "-i");
  const char *seed_str = getCmdOption(argv, argv + argc, "-s");
  const char *algorithms_str = getCmdOption(argv, argv + argc, "-a");
  const char *runs_str = getCmdOption(argv, argv + argc, "-r");
  const char *neighbors_str = getCmdOption(argv, argv + argc, "-k");
  const char *varyk_str = getCmdOption(argv, argv + argc, "-vk");
  const char *varyc_str = getCmdOption(argv, argv + argc, "-vc");
  const char *varyt_str = getCmdOption(argv, argv + argc, "-vt");
  const char *clusters_str = getCmdOption(argv, argv + argc, "-c");
  const char *num_agents_str = getCmdOption(argv, argv + argc, "-n");
  if (cmdOptionExists(argv, argv + argc, "-o")) csv_output = true;

  if (input_file_str)
    parameters["read_agents"] = 1.0;
  else
    parameters["read_agents"] = 0.0;

  if (seed_str) seed = std::stol(std::string(seed_str));
  if (algorithms_str) {
    algorithms = std::string(algorithms_str);
    std::size_t found = algorithms.find_first_of(",");
    if (found != std::string::npos) {
      std::string int_part(algorithms,0, found);
      size_t n = std::stol(int_part);
      std::string alg_part(algorithms,found + 1);
      algorithms = string("");
      for (size_t i = 0; i < n; ++i)
	algorithms += alg_part;
    }
  }
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
  if (varyt_str)
    parameters["varyt"] =  (double) std::stol(std::string(varyt_str));
  else
    parameters["varyt"] =  (double) varyt;
  if (num_agents_str)
    parameters["agents"] = (double) std::stol(std::string(num_agents_str));
  else
    parameters["agents"] = num_agents;

  if (clusters_str)
    parameters["clusters"] =  (double) std::stol(std::string(clusters_str));
  else
    parameters["clusters"] =  (double) clusters;

  rng.seed(seed);
  run_tests(parameters, algorithms, runs, input_file_str, csv_output);
}
