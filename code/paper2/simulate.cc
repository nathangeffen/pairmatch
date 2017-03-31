/**
    stisimulator: simulate.cc
    Purpose: Microsimulation of sexually transmitted infection epidemics

    @author Nathan Geffen
    @version 0.1 23/2/2017
    @license GPL v3
*/

#include <algorithm>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <random>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <cfloat>
#include <cstdint>
#include <sys/time.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "CSV_Parser.hh"
#include "sample.hh"

thread_local std::mt19937 rng;

#define DAY 1.0 / 365;
#define MIN_AGE 12
#define MAX_AGE 100

#define RPM 1
#define RKPM 2
#define CSPM 3

#define HEURISTIC_DISTANCE 0
#define TABLE_DISTANCE 1

#define GRAPH_ACCURACY 10000

#define FEMALE 0
#define MALE 1
#define HOMOSEXUAL 0
#define HETEROSEXUAL 1

/**
    Macro used for unit testing.  Tests if two values are equal, incrementing
    successes if true, else incrementing failures.

    @param x the first value to test
    @param y the second value to test
    @param successes integer to increment if test is successful
    @param failures integer to increment if test fails
*/

#define TESTEQ(x, y, successes, failures)                               \
  do {                                                                  \
    auto _t1 = (x);                                                     \
    auto _t2 = (y);                                                     \
    std::string _t3(#x);                                                \
    std::string _t4(#y);                                                \
    if (_t1 == _t2) {                                                   \
      cout << "PASS:\t" << _t3 << " == " << _t4                         \
           << "\tLine:" << __LINE__ << "\n";                            \
      ++successes;                                                      \
    }                                                                   \
    else {                                                              \
      cout << "FAIL:\t" << _t3 << " == " << _t4                         \
           << "\t" << _t1 << " != " << _t2                              \
           << "\tLine:" << __LINE__ << "\n";                            \
      ++failures;                                                       \
    }                                                                   \
  } while(0)

/**
   Linear algebra functions.
*/

/**
    Gets a column in a matrix.

    @param matrix the matrix to get the column from
    @param col the zero-indexed column in the matrix to return
    @return matrix column as a vector of doubles
*/

std::vector<double> getCol(const DblMatrix& matrix, unsigned col)
{
  std::vector<double> output;
  for (unsigned i = 0; i < matrix.size(); ++i)
    output.push_back(matrix[i][col]);
  return output;
}

/**
    Calculates a vector subtracted from a scalar.

    @param d scalar from which to subtract
    @param v vector to subtract
    @return vector = d - v
*/

std::vector<double> subVector(double d,
                              const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d - v[i]);
  return output;
}

/**
    Calculates a vector multiplied by a scalar.

    @param d scalar to multiply
    @param v vector to multiply
    @return vector = d * v
*/

std::vector<double> multVector(double d,
                               const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d * v[i]);
  return output;
}

/**
    Calculates product of two vectors.

    @param v1 vector to multiply
    @param v2 vector to multiply
    @return vector = v1 X v2
*/

std::vector<double> multVectors(const std::vector<double>& v1,
                                const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] * v2[i]);
  return output;
}

/**
    Calculates sum of two vectors.

    @param v1 vector to add
    @param v2 vector to add
    @return vector = v1 + v2
*/

std::vector<double> addVector(const std::vector<double>& v1,
                              const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] + v2[i]);
  return output;
}

/**
    Calculates sum of elements of vector.

    @param v vector to sum
    @return sum of elements of
*/

double sumVector(const std::vector<double>& v)
{
  double total = 0.0;
  for (auto d: v)
    total += d;
  return total;
}

// End of linear algebra functions
////////////////////////////////////////


/**
    Calculates the number of people who are single in the initial population.

    @param data matrix of probabilities
    @return number of people who are single
*/


double calcNumberSingles(DblMatrix data, unsigned X)
{
  auto ageshare = getCol(data, 1);
  auto femratio = getCol(data, 2);
  auto relm_share = getCol(data, 3);
  auto relw_share = getCol(data, 4);

  auto t1 = subVector(1, relw_share);
  auto t2 = multVectors(femratio, t1);
  auto t3 = multVectors(t2, ageshare);

  auto t4 = subVector(1, femratio);;
  auto t5 = subVector(1, relm_share);
  auto t6 = multVectors(t4, t5);
  auto t7 = multVectors(t6, ageshare);

  auto t8 = addVector(t3, t7);
  auto t9 = sumVector(t8);

  double num_agents = X * t9;
  return round(num_agents);
}


///////////////////////////////////////////////
// Parameter classes and functions

struct ParameterValue {
  double getDbl(size_t index = 0) const
  {
    return values[0];
  };

  inline const bool isSet() const
  {
    return (values[0] != 0.0) ? true : false;
  }
  inline std::string getStr() const
  {
    return stringValue;
  };
  const char* getCStr() const
  {
    return stringValue.c_str();
  }
  const char* description;
  std::vector<double> values;
  bool isString = false;
  std::string stringValue;
};

typedef std::map<std::string,  ParameterValue> ParameterMap;

void printParameters(ParameterMap& parameterMap)
{
  for (auto& p: parameterMap) {
    printf("PARAMETER,%s,%s", p.first.c_str(), p.second.description);
    if (p.second.isString)
      printf(",%s", p.second.stringValue.c_str());
    else
      for (auto& v: p.second.values)
        printf(",%f", v);
    printf("\n");
  }
}

void addParameter(ParameterMap& parameterMap,
                  const char* name,
                  const char* description,
                  std::initializer_list<double> values)
{
  ParameterValue parameterValue;

  parameterValue.description = description;
  parameterValue.values = values;
  parameterMap[name] = parameterValue;
}

void addParameter(ParameterMap& parameterMap,
                  const char* name,
                  const char* description,
                  const char* stringValue)
{
  ParameterValue parameterValue;

  parameterValue.description = description;
  parameterValue.isString = true;
  parameterValue.stringValue = stringValue;
  parameterMap[name] = parameterValue;
}


bool replaceParameter(ParameterMap& parameterMap,
                      const char* name,
                      std::istringstream& line,
                      std::string& errorMessage)
{
  ParameterValue parameterValue;
  auto it = parameterMap.find(name);
  if (it == parameterMap.end()) {
    std::ostringstream errorMessageStream;
    errorMessageStream << "Unknown parameter: " << name;
    errorMessage = errorMessageStream.str();
    return false;
  }
  if (parameterMap.at(name).isString) {
    line >> parameterMap[name].stringValue;
  } else {
    double value;
    std::vector<double> values;

    while (line >> value)
      values.push_back(value);
    if ( (line.fail() && !line.eof()) || values.size() == 0) {
      std::ostringstream errorMessageStream;
      errorMessageStream << "Invalid values for parameter: " << name;
      errorMessage = errorMessageStream.str();
      return false;
    }
    parameterMap[name].values = values;
  }
  return true;
}


/* Add new parameters here. */

void setDefaultParameters(ParameterMap& parameterMap)
{
  addParameter(parameterMap, "SIMULATION_NAME",
               "Name of simulation used in output", "Default");
  addParameter(parameterMap, "NUM_SIMULATIONS",
               "Number of simulations to execute (default is 1)", {1});
  addParameter(parameterMap, "NUM_THREADS",
               "Number of threads (default is one per simulation - 0)", {0});
  addParameter(parameterMap, "NUM_AGENTS", "Number of agents", {100.0});
  addParameter(parameterMap, "MALE_INFECTIVITY", "How easily males are infected",
               {0.1});
  addParameter(parameterMap, "FEMALE_INFECTIVITY",
               "How easily females are infected", {0.2});
  addParameter(parameterMap, "MALE_INFECTIOUSNESS",
               "How likely males are to infect their partners", {0.2});
  addParameter(parameterMap, "FEMALE_INFECTIOUSNESS",
               "How easily females are infected", {0.1});
  addParameter(parameterMap, "BREAKUPINESS",
               "How easily partners breakup", {0.1});
  addParameter(parameterMap, "START_DATE",
               "Start date of simulation", {2017.0});
  addParameter(parameterMap, "END_DATE",
               "End date of simulation", {2018.0});
  addParameter(parameterMap, "TIME_STEP",
               "Time step for each iteration of simulation", {1.0 / 365.0});
  addParameter(parameterMap, "PREV_PARTNER_PENALTY",
               "Previous partner penalty in pair matching", {10.0});
  addParameter(parameterMap, "MEAN_DAYS_RELATIONSHIP",
               "Mean number of days agents are partners", {365.0});
  addParameter(parameterMap, "MEAN_DAYS_SINGLE",
               "Mean number of days agents are single", {365.0});

  addParameter(parameterMap, "HET_MALE_INFECTIVITY",
               "Daily risk of getting infected for male in heterosexual sex",
               {0.001});
  addParameter(parameterMap, "HET_FEMALE_INFECTIVITY",
               "Daily risk of getting infected for female in heterosexual sex",
               {0.005});
  addParameter(parameterMap, "HET_MALE_INFECTIOUSNESS",
               "Daily risk of infecting for male in heterosexual sex",
               {0.005});
  addParameter(parameterMap, "HET_FEMALE_INFECTIOUSNESS",
               "Daily risk of infecting for female in heterosexual sex",
               {0.001});

  addParameter(parameterMap, "HOM_MALE_INFECTIVITY",
               "Daily risk of getting infected for male in homosexual sex",
               {0.005});
  addParameter(parameterMap, "HOM_FEMALE_INFECTIVITY",
               "Daily risk of getting infected for female in homosexual sex",
               {0.0001});
  addParameter(parameterMap, "HOM_MALE_INFECTIOUSNESS",
               "Daily risk of infecting for male in homosexual sex",
               {0.005});
  addParameter(parameterMap, "HOM_FEMALE_INFECTIOUSNESS",
               "Daily risk of infecting for female in homosexual sex",
               {0.0001});

  addParameter(parameterMap, "MATCH_NEIGHBORS",
               "Value of k (neighbors for matching", {300});
  addParameter(parameterMap, "MATCH_CLUSTERS",
               "Number of clusters for cluster shuffle matching", {100});
  addParameter(parameterMap, "MATCH_SCORE_FAIL",
               "Highest score beyond which a match between two agents must fail",
               {2000.00});
  addParameter(parameterMap, "MATCH_SCORE_POOR",
               "Score for which a poor match must be registered",
               {49.50});

  addParameter(parameterMap, "OUTPUT_AGENTS_AFTER_INIT",
               "Print agent info after initialization"
               "(0=no,1=yes)", {0.0});
  addParameter(parameterMap, "OUTPUT_AGENTS_DURING_SIM",
               "Print agent info during simulation"
               "(0=no,1=yes), (frequency of output)", {0.0, 100.0});
  addParameter(parameterMap, "OUTPUT_AGENTS_AT_END",
               "Print agent info at end of simulation "
               "(0=no,1=yes)", {1.0});

  addParameter(parameterMap, "ANALYSIS_AGE_INTERVAL",
               "Length of age intervals for infection analysis", {5.0});
  addParameter(parameterMap, "ANALYZE_AFTER_INIT",
               "Calculate stats after initialization"
               "(0=no,1=yes)", {1.0});
  addParameter(parameterMap, "ANALYZE_DURING_SIM",
               "Calculate stats during simulation"
               "(0=no,1=yes), (frequency of output)", {0.0, 1.0});
  addParameter(parameterMap, "ANALYZE_AT_END",
               "Calculate stats at end of simulation "
               "(0=no,1=yes)", {1.0});
  addParameter(parameterMap, "PRINT_PARAMETERS",
               "Print parameters using in simulation"
               "(0=no,1=yes)", {0.0});

  addParameter(parameterMap, "TIMING", "Output time every x iterations", {0});

  addParameter(parameterMap, "AGENT_DATA_CSV",
               "Agent data for initialization", "data/data.csv");
  addParameter(parameterMap, "SINGLES_DATA_CSV",
               "Agent data for initialization", "data/singles.csv");
  addParameter(parameterMap, "PARTNERS_DATA_CSV",
               "Agent data for initialization", "data/partnersACP.csv");
  addParameter(parameterMap, "MSM_DATA_CSV",
               "Agent data for initialization", "data/mat.msm.csv");
  addParameter(parameterMap, "MSW_DATA_CSV",
               "Agent data for initialization", "data/mat.msw.csv");
  addParameter(parameterMap, "WSW_DATA_CSV",
               "Agent data for initialization", "data/mat.wsw.csv");
  addParameter(parameterMap, "WSM_DATA_CSV",
               "Agent data for initialization", "data/mat.wsm.csv");

  addParameter(parameterMap, "MSM_AGE_DIST_CSV",
               "MSM age preferences CSV file", "data/dist_msm.csv");
  addParameter(parameterMap, "MSW_AGE_DIST_CSV",
               "MSW age preferences CSV file", "data/dist_msw.csv");
  addParameter(parameterMap, "WSW_AGE_DIST_CSV",
               "WSW age preferences CSV file", "data/dist_wsw.csv");
  addParameter(parameterMap, "WSM_AGE_DIST_CSV",
               "WSM age preferences CSV file", "data/dist_wsm.csv");

  addParameter(parameterMap, "INITIAL_INFECTION_RATES_CSV",
               "Initial infection rates in population",
               "data/initial_rates.csv");

  addParameter(parameterMap, "AGE_EVENT", "Execute the aging event", {1});
  addParameter(parameterMap, "INFECT_EVENT", "Execute the infection event", {1});
  addParameter(parameterMap, "BREAKUP_EVENT", "Execute the infection event", {1});
  addParameter(parameterMap, "MATCH_EVENT",
               "Execute the infection event (RPM,RKPM,CSPM)", "RPM");
  addParameter(parameterMap, "DISTANCE_METHOD",
               "Distance calculation to use (HEURISTIC,TABLE)", "HEURISTIC");

  addParameter(parameterMap, "PARTNERS_LOAD_FACTOR",
               "Load factor for partnerships unordered set (advanced)", {0});
  addParameter(parameterMap, "PARTNERS_RESERVE",
               "Capacity for partnerships unordered set (advanced)", {0});

  addParameter(parameterMap, "RANDOM_SEED", "Value to set random seed to",
               {1});
}

void readParameters(std::istream& input,
                    std::vector<ParameterMap>& parameterMaps)
{
  std::string lineString;
  unsigned lineNumber = 0;
  std::string errorMessage;
  parameterMaps.push_back(ParameterMap());
  setDefaultParameters(parameterMaps.back());
  while (std::getline(input, lineString)) {
    ++lineNumber;
    boost::trim(lineString);
    if (lineString.size() == 0) continue;
    if (lineString[0] == '#') continue;
    if (lineString[0] == '-') {
      parameterMaps.push_back(ParameterMap());
      setDefaultParameters(parameterMaps.back());
      continue;
    }
    std::istringstream line(lineString);
    std::string parameterName;
    line >> parameterName;
    if (replaceParameter(parameterMaps.back(), parameterName.c_str(),
                         line, errorMessage) == false) {
      fprintf(stderr,"Error at line: %u\n", lineNumber);
      fprintf(stderr, "%s\n", errorMessage.c_str());
    }
  }
}

// End of Parameter management

////////////////////////////////////
// Agent class and related functions

// Used to keep track of past partnerships

class Agent {
public:
  uint32_t id = 0;
  Agent* partner;
  Agent* infector = NULL;

  float het_infectivity;
  float het_infectiousness;
  float hom_infectivity;
  float hom_infectiousness;
  float age;
  float desired_age;
  float relationship_change_date;
  float binomial_p_relationship_length;
  float binomial_p_relationship_wait;
  float weight;

  unsigned short sex;
  unsigned short sexual_orientation;
  bool infected;
  bool initial_relationship;

  unsigned num_partners = 0;
  unsigned num_infected = 0;


  void setRelationshipLength(const double current_date,
                             const double mean_days_relationship)
  {
    /* This needs to be replaced with Stefan's relationship status equations. */
    double p = (binomial_p_relationship_length +
                partner->binomial_p_relationship_length) / 2.0;
    std::binomial_distribution<int> dist (mean_days_relationship, p);
    double num_days = (double) dist(rng) * DAY;
    /****************************************/
    relationship_change_date = current_date + num_days;
    partner->relationship_change_date = relationship_change_date;
  }
  void setSingleLength(double current_date,
                       double mean_days_single)
  {
    /* This needs to be replaced with Stefan's relationship status equations. */
    std::binomial_distribution<int> dist (mean_days_single,
                                          binomial_p_relationship_wait);
    /****************************************/
    relationship_change_date = current_date + (double) dist(rng) * DAY;
  }

  bool isMatchable(const double current_date) const
  {
    if (partner == NULL && current_date >= relationship_change_date)
      return true;
    else
      return false;
  }

  void print(FILE *f = stdout) const
  {
    fprintf(f, "ID,%7u,Age,%3.2f,Sex,%c,Orientation,%c,Desired,%.0f",
            id, age, (sex == MALE ? 'M' : 'F'),
            (sexual_orientation == HETEROSEXUAL ? 'S' : 'G'), desired_age);
    fprintf(f, ",Risk,%.3f,%.3f,Partner",
            binomial_p_relationship_length,binomial_p_relationship_wait);
    if (partner)
      fprintf(f,",%7u", partner->id);
    else
      fprintf(f, ",      0");
    fprintf(f, ",Date,%.3f,Infection,%d,%.3f,%.3f,%.3f\n",
            relationship_change_date, infected, het_infectiousness,
            hom_infectivity, hom_infectiousness);
  }
};

// We don't use the ostream operator here, preferring the C stdio library
// because it is implemented much faster by the GNU compiler.

ostream& operator<<(ostream& os, const Agent& agent)
{
  os << "ID," << agent.id << ",Age," << agent.age
     << ",Sex," << (agent.sex == MALE ? "M" : "F")
     << ",Orientation," << (agent.sexual_orientation == HETEROSEXUAL
                            ? "S" : "G")
     << ",Desired," << agent.desired_age;
  std::ios_base::fmtflags oldflags = os.flags();
  std::streamsize oldprecision = os.precision();
  os << std::fixed << std::setprecision(3)
     << ",Risk," << agent.binomial_p_relationship_length
     << "," << agent.binomial_p_relationship_wait;
  os.flags (oldflags);
  os.precision (oldprecision);
  if (agent.partner) {
    auto &p = agent.partner;
    os << ",ID," << p->id << ",Age," << p->age
       << ",Sex," << (p->sex == MALE ? "M" : "F")
       << ",Orientation," << (p->sexual_orientation == HETEROSEXUAL
                              ? "S" : "G")
       << ",Desired," << p->desired_age;
  }
  os << ",Date," << agent.relationship_change_date;
  os << ",Infection," << agent.het_infectivity << ","
     << agent.het_infectiousness << ","
     << agent.hom_infectivity << ","
     << agent.hom_infectiousness;
  return os;
}

typedef std::vector<Agent *> AgentVector;

void printAgents(const AgentVector& agents,
                 unsigned simulation_num,
                 double current_date,
                 FILE *out)
{
  for (auto& agent: agents) {
    fprintf(out, "AGENT,%u,%f,", simulation_num, current_date);
    agent->print(out);
  }
}

// End of Agents
////////////////////////////////////////


// Simulation class and supporting routines

// Structure to keep track of partnerships

class Partnerships {
public:
  void insert(const uint32_t id1, const uint32_t id2)
  {
    partnerships.insert(combine(id1, id2));
  }

  bool exists(const uint32_t id1, const uint32_t id2) const
  {
    if (partnerships.find(combine(id1, id2)) == partnerships.end())
      return false;
    else
      return true;
  }
  size_t size() const
  {
    return partnerships.size();
  }

  //private:
  uint64_t combine(uint32_t id1, uint32_t id2) const
  {
    uint64_t A = id1 < id2 ? id1 : id2;
    uint64_t B = id1 > id2 ? id1 : id2;
    uint64_t C = A << 32 | B;
    return C;
  }
  std::unordered_set<uint64_t> partnerships;
};


struct PartnershipScore {
  std::vector<Agent *>::iterator partner;
  double score;
};

// Initialization routines

void setAgentInfectionParameters(Agent& agent,
                                 double het_male_infectivity,
                                 double het_female_infectivity,
                                 double het_male_infectiousness,
                                 double het_female_infectiousness,
                                 double hom_male_infectivity,
                                 double hom_female_infectivity,
                                 double hom_male_infectiousness,
                                 double hom_female_infectiousness)
{
  agent.het_infectivity =
    agent.sex == MALE ? het_male_infectivity : het_female_infectivity;
  agent.hom_infectivity =
    agent.sex == MALE ? hom_male_infectivity : hom_female_infectivity;
  agent.het_infectiousness =
    agent.sex == MALE ? het_male_infectiousness : het_female_infectiousness;
  agent.hom_infectiousness =
    agent.sex == MALE ? hom_male_infectiousness : hom_female_infectiousness;
}

void setInitialInfection(Agent &agent,
                         const std::vector<double>& initialInfectionRatesMSW,
                         const std::vector<double>& initialInfectionRatesMSM,
                         const std::vector<double>& initialInfectionRatesWSM,
                         const std::vector<double>& initialInfectionRatesWSW)
{
  std::uniform_real_distribution<double> uni;
  if (agent.sex == MALE && agent.sexual_orientation == HETEROSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesMSW[(int) agent.age])
      ? true : false;
  } else if (agent.sex == MALE && agent.sexual_orientation == HOMOSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesMSM[(int) agent.age])
      ? true : false;
  } else if (agent.sex == FEMALE && agent.sexual_orientation == HETEROSEXUAL) {
    agent.infected = (uni(rng) < initialInfectionRatesWSM[(int) agent.age])
      ? true : false;
  } else {
    agent.infected = (uni(rng) < initialInfectionRatesWSW[(int) agent.age])
      ? true : false;
  }
}


/****************************/

class Simulation {
public:
  Simulation(const ParameterMap& parameter_map,
             const unsigned simulation_num) :
    parameterMap(parameter_map), simulationNum(simulation_num)
  {
    simulationName = parameterMap.at("SIMULATION_NAME").getStr();
    startDate = parameterMap.at("START_DATE").getDbl();
    endDate = parameterMap.at("END_DATE").getDbl();
    timeStep = parameterMap.at("TIME_STEP").getDbl();
    meanDaysRelationship = parameterMap.at("MEAN_DAYS_RELATIONSHIP").getDbl();
    meanDaysSingle = parameterMap.at("MEAN_DAYS_SINGLE").getDbl();
    currentDate = startDate;
    ageInterval = parameterMap.at("ANALYSIS_AGE_INTERVAL").getDbl();
    neighbors = parameterMap.at("MATCH_NEIGHBORS").getDbl();
    clusters = parameterMap.at("MATCH_CLUSTERS").getDbl();
    failureThresholdScore = parameterMap.at("MATCH_SCORE_FAIL").getDbl();
    poorThresholdScore = parameterMap.at("MATCH_SCORE_POOR").getDbl();
    string s = parameterMap.at("DISTANCE_METHOD").getStr();
    if (s == "HEURISTIC")
      distanceMethod = HEURISTIC_DISTANCE;
    else if (s == "TABLE")
      distanceMethod = TABLE_DISTANCE;
    else {
      fprintf(stderr, "Unknown distance method.\n");
      exit(1);
    }

    size_t max_p = (size_t) parameterMap.at("PARTNERS_RESERVE").getDbl();
    if (max_p > 0) partnerships.partnerships.reserve(max_p);
    auto lf = parameterMap.at("PARTNERS_LOAD_FACTOR").getDbl();
    if (lf > 0) partnerships.partnerships.max_load_factor(lf);

    unsigned seed = parameterMap.at("RANDOM_SEED").getDbl();
    if (seed) rng.seed(seed * simulation_num);
  }
  ~Simulation()
  {
    for (auto& agent: agents) delete agent;
  }
  void simulate()
  {
    unsigned timing = parameterMap.at("TIMING").getDbl();

    // Initialization
    initializeSimulation();
    struct timeval timeBegin, timeEnd;
    double elapsedTime;

    gettimeofday(&timeBegin, NULL);
    initializeAgents();
    gettimeofday(&timeEnd, NULL);
    elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
    printf("%s,TIMING,INIT,%u,%f\n", simulationName.c_str(),
           simulationNum, elapsedTime);
    unsigned outputAgents = parameterMap.at("OUTPUT_AGENTS_AFTER_INIT").getDbl();
    if (outputAgents) printAgents(agents, simulationNum, startDate, stdout);
    /********************/
    /* Count WSW */
    unsigned wswcount = 0;
    for (auto& a: agents) {
      if (a->sex == FEMALE && a->sexual_orientation == HOMOSEXUAL) ++wswcount;
    }
    /*******************/
    unsigned analyzeAgents = parameterMap.at("ANALYZE_AFTER_INIT").getDbl();
    if (analyzeAgents) analysis();
    /* Main loop */

    // Make sure main loop uses integer arithmetic, rather than floats
    // though probably makes little diff.
    unsigned num_iterations = (endDate - startDate) / timeStep;
    outputAgents = parameterMap.at("OUTPUT_AGENTS_DURING_SIM").getDbl(0);
    unsigned outputFrequency =
      parameterMap.at("OUTPUT_AGENTS_DURING_SIM").getDbl(1);
    analyzeAgents = parameterMap.at("ANALYZE_DURING_SIM").getDbl(0);
    unsigned analyze_frequency =
      parameterMap.at("ANALYZE_DURING_SIM").getDbl(1);
    for (unsigned i = 0; i < num_iterations; ++i, currentDate += timeStep) {
      for (auto& e: events) e(this);

      if (timing > 0 && (i + 1) % timing == 0) {
        gettimeofday(&timeEnd, NULL);
        elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
        printf("%s,TIMING,%u,%u,%f\n",
               simulationName.c_str(), i, simulationNum, elapsedTime);
      }
      if (outputAgents > 0 && (i + 1) % outputFrequency == 0)
        printAgents(agents, simulationNum, currentDate, stdout);
      if (analyzeAgents > 0 && (i + 1) % analyze_frequency == 0) analysis();
    }

    gettimeofday(&timeEnd, NULL);
    elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
    printf("%s,TIMING,AFTER,%u,%f\n",
           simulationName.c_str(), simulationNum, elapsedTime);

    /* Wrap up */
    outputAgents = parameterMap.at("OUTPUT_AGENTS_AT_END").getDbl();
    if (outputAgents) printAgents(agents, simulationNum, endDate, stdout);
    analyzeAgents = parameterMap.at("ANALYZE_AT_END").getDbl();
    if (analyzeAgents) analysis();
  }
  void initializeAgents()
  {
    unsigned numIntervals = (MAX_AGE + 1) / ageInterval + 1;
    malesByAge = std::vector<unsigned>(numIntervals, 0);
    femalesByAge = std::vector<unsigned>(numIntervals, 0);
    infectedMalesByAge = std::vector<unsigned>(MAX_AGE + 1, 0);
    infectedFemalesByAge = std::vector<unsigned>(MAX_AGE + 1, 0);

    CSVParser data_csv(parameterMap.at("AGENT_DATA_CSV").stringValue.c_str(),
                       ";", true);
    DblMatrix data = data_csv.convert_all_entries_to_doubles();
    CSVParser singles_csv(parameterMap.at("SINGLES_DATA_CSV").stringValue.c_str(),
                          ",", true);
    DblMatrix singles = singles_csv.convert_all_entries_to_doubles();
    CSVParser partners_csv(parameterMap.
                           at("PARTNERS_DATA_CSV").stringValue.c_str(),
                           ",", true);

    DblMatrix partners = partners_csv.convert_all_entries_to_doubles();
    CSVParser mm_csv(parameterMap.at("MSM_DATA_CSV").stringValue.c_str(), ";",
                     false);
    DblMatrix mm = mm_csv.convert_all_entries_to_doubles();
    CSVParser ww_csv(parameterMap.at("WSW_DATA_CSV").stringValue.c_str(), ";",
                     false);
    DblMatrix ww = ww_csv.convert_all_entries_to_doubles();
    CSVParser mw_csv(parameterMap.at("MSW_DATA_CSV").stringValue.c_str(), ";",
                     false);
    DblMatrix mw = mw_csv.convert_all_entries_to_doubles();
    CSVParser wm_csv(parameterMap.at("WSM_DATA_CSV").stringValue.c_str(), ";",
                     false);
    DblMatrix wm = wm_csv.convert_all_entries_to_doubles();

    CSVParser infectionRatesCsv(parameterMap.at("INITIAL_INFECTION_RATES_CSV").
                                stringValue.c_str(), ",", true);
    DblMatrix initialRatesMatrix =
      infectionRatesCsv.convert_all_entries_to_doubles();
    std::vector<double> initialInfectionRatesMSW(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesMSM(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesWSM(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesWSW(MAX_AGE + 1, 0.0);
    {
      unsigned i = 0;
      for (auto& row : initialRatesMatrix) {
        for (;i <= row[0] && i < MAX_AGE + 1; ++i) {
          assert(i <= MAX_AGE);
          initialInfectionRatesMSW[i] = row[1];
          initialInfectionRatesMSM[i] = row[2];
          initialInfectionRatesWSM[i] = row[3];
          initialInfectionRatesWSW[i] = row[4];
        }
      }
      for (; i <= MAX_AGE; ++i) {
          initialInfectionRatesMSW[i] = 0.0;
          initialInfectionRatesMSM[i] = 0.0;
          initialInfectionRatesWSM[i] = 0.0;
          initialInfectionRatesWSW[i] = 0.0;
      }
    }

    unsigned X = parameterMap.at("NUM_AGENTS").values[0];
    unsigned S = calcNumberSingles(data, X);
    agents.reserve(X);

    std::vector<double> ageRange;
    for (unsigned i = 12; i <= 100; ++i)
      ageRange.push_back(i);
    auto ageShare = getCol(singles, 1);
    auto femRatio = getCol(singles, 2);
    auto msmRate = getCol(singles, 3);
    auto wswRate = getCol(singles, 4);
    createPartners(agents, parameterMap, data, 0, S, ageRange, ageShare, femRatio,
                   wswRate, msmRate, ww, mw, wm, mm,
                   initialInfectionRatesMSW, initialInfectionRatesMSM,
                   initialInfectionRatesWSM, initialInfectionRatesWSW, false);
    ageShare = femRatio = msmRate = wswRate = {};
    ageShare = getCol(partners, 1);
    femRatio = getCol(partners, 2);
    msmRate = getCol(partners, 3);
    wswRate = getCol(partners, 4);
    createPartners(agents, parameterMap, data, S, X, ageRange, ageShare, femRatio,
                   wswRate, msmRate, ww, mw, wm, mm,
                   initialInfectionRatesMSW, initialInfectionRatesMSM,
                   initialInfectionRatesWSM, initialInfectionRatesWSW, true);
    for (unsigned i = 0; i < agents.size(); ++i) {
      unsigned age = agents[i]->age;
      if (agents[i]->sex == MALE) {
        ++numMales;
        ++malesByAge[age / ageInterval];
        if (agents[i]->sexual_orientation == HETEROSEXUAL)
          ++numMsw;
        else
          ++numMsm;
        if (agents[i]->infected) {
          ++numInfectedMales;
          ++infectedMalesByAge[age / ageInterval];
          if (agents[i]->sexual_orientation == HETEROSEXUAL)
            ++numInfectedMsw;
          else
            ++numInfectedMsm;
        }
      } else {
        ++numFemales;
        ++femalesByAge[age /  ageInterval];
        if (agents[i]->sexual_orientation == HETEROSEXUAL)
          ++numWsm;
        else
          ++numWsw;
        if (agents[i]->infected) {
          ++numInfectedFemales;
          ++infectedFemalesByAge[age / ageInterval];
          if (agents[i]->sexual_orientation == HETEROSEXUAL)
            ++numInfectedWsm;
          else
            ++numInfectedWsw;
        }
      }
    }
    agents.shrink_to_fit();
  }

  inline void initAgent(Agent *agent,
                        const unsigned id,
                        Sample& sampleAgeshare,
                        const std::vector<double>& femRatio,
                        const std::vector<double>& wswRate,
                        const std::vector<double>& msmRate,
                        const std::vector<double>& initialInfectionRatesMSW,
                        const std::vector<double>& initialInfectionRatesMSM,
                        const std::vector<double>& initialInfectionRatesWSM,
                        const std::vector<double>& initialInfectionRatesWSW,
                        std::vector<Sample>& sample_matWW,
                        std::vector<Sample>& sample_matMW,
                        std::vector<Sample>& sample_matWM,
                        std::vector<Sample>& sample_matMM)
  {
    std::uniform_real_distribution<double> uni;

    agent->id = id;
    // Age
    unsigned age = sampleAgeshare() + 12;
    agent->age = age;
    // Sex
    unsigned sex = uni(rng) < femRatio[age - 12] ? FEMALE : MALE;
    agent->sex = sex;
    // Sexual_Orientation
    unsigned sexual_orientation;
    if (sex == FEMALE) {
      sexual_orientation = uni(rng) < wswRate[age - 12]
                                      ? HOMOSEXUAL : HETEROSEXUAL;
    } else {
      sexual_orientation = uni(rng) < msmRate[age - 12]
                                      ? HOMOSEXUAL : HETEROSEXUAL;
    }
    agent->sexual_orientation = sexual_orientation;
    setInitialInfection(*agent, initialInfectionRatesMSW,
                        initialInfectionRatesMSM, initialInfectionRatesWSM,
                        initialInfectionRatesWSW);
    // Desire age of partner
    if (sex == FEMALE && sexual_orientation == HOMOSEXUAL)
      agent->desired_age  = sample_matWW[age - 12]() + 12;
    else if (sex == FEMALE && sexual_orientation == HETEROSEXUAL)
      agent->desired_age = sample_matWM[age - 12]() + 12;
    else if (sex == MALE && sexual_orientation == HETEROSEXUAL)
      agent->desired_age = sample_matMW[age - 12]() + 12;
    else
      agent->desired_age = sample_matMM[age - 12]() + 12;

    agent->binomial_p_relationship_length = uni(rng);
    agent->binomial_p_relationship_wait = uni(rng);
    if (agent->binomial_p_relationship_length >= 0.5 &&
        agent->binomial_p_relationship_wait >= 0.5) {
      ++numLongBreakLongPartnership;
    } else if (agent->binomial_p_relationship_length >= 0.5 &&
               agent->binomial_p_relationship_wait < 0.5) {
      ++numShortBreakLongPartnership;
    } else if (agent->binomial_p_relationship_length < 0.5 &&
               agent->binomial_p_relationship_wait >= 0.5) {
      ++numLongBreakShortPartnership;
    } else {
      ++numShortBreakShortPartnership;
    }
  }

  void createPartners(AgentVector& agents,
                      const ParameterMap& parameterMap,
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
                      const DblMatrix& matMM,
                      const std::vector<double>& initialInfectionRatesMSW,
                      const std::vector<double>& initialInfectionRatesMSM,
                      const std::vector<double>& initialInfectionRatesWSM,
                      const std::vector<double>& initialInfectionRatesWSW,
                      bool initial_relation)
  {
    std::uniform_real_distribution<double> uni;
    double het_male_infectivity =
      parameterMap.at("HET_MALE_INFECTIVITY").getDbl();
    double het_male_infectiousness =
      parameterMap.at("HET_MALE_INFECTIOUSNESS").getDbl();
    double het_female_infectivity =
      parameterMap.at("HET_FEMALE_INFECTIVITY").getDbl();
    double het_female_infectiousness =
      parameterMap.at("HET_FEMALE_INFECTIOUSNESS").getDbl();
    double hom_male_infectivity =
      parameterMap.at("HOM_MALE_INFECTIVITY").getDbl();
    double hom_male_infectiousness =
      parameterMap.at("HOM_MALE_INFECTIOUSNESS").getDbl();
    double hom_female_infectivity =
      parameterMap.at("HOM_FEMALE_INFECTIVITY").getDbl();
    double hom_female_infectiousness =
      parameterMap.at("HOM_FEMALE_INFECTIOUSNESS").getDbl();
    Sample sample_ageshare(ageShare, &rng);
    vector<Sample> sample_matWW(matWW[0].size());
    vector<Sample> sample_matMW(matMW[0].size());
    vector<Sample> sample_matWM(matWM[0].size());
    vector<Sample> sample_matMM(matMM[0].size());
    std::vector<double> placeholder(matMM.size(), 0.000001);

    for (unsigned i = 0; i < matWW[0].size(); ++i) {
      sample_matWW[i].init(getCol(matWW,i), &rng);
      sample_matWM[i].init(getCol(matWM,i), &rng);
      sample_matMW[i].init(getCol(matMW,i), &rng);
      sample_matMM[i].init(getCol(matMM,i), &rng);
    }


    unsigned increment = initial_relation ? 2 : 1;
    for (unsigned i = fromAgent; (i + increment - 1) < toAgent; i += increment) {
      Agent *agent = new Agent();
      initAgent(agent, i + 1,
                sample_ageshare,
                femRatio,
                wswRate,
                msmRate,
                initialInfectionRatesMSW,
                initialInfectionRatesMSM,
                initialInfectionRatesWSM,
                initialInfectionRatesWSW,
                sample_matWW,
                sample_matMW,
                sample_matWM,
                sample_matMM);
    setAgentInfectionParameters(*agent,
                                het_male_infectivity,
                                het_female_infectivity,
                                het_male_infectiousness,
                                het_female_infectiousness,
                                hom_male_infectivity,
                                hom_female_infectivity,
                                hom_male_infectiousness,
                                hom_female_infectiousness);
      if (initial_relation) {
        // Relationship
        agent->initial_relationship = true;
        Agent* partner = new Agent();
        initAgent(partner, i + 2,
                  sample_ageshare,
                  femRatio,
                  wswRate,
                  msmRate,
                  initialInfectionRatesMSW,
                  initialInfectionRatesMSM,
                  initialInfectionRatesWSM,
                  initialInfectionRatesWSW,
                  sample_matWW,
                  sample_matMW,
                  sample_matWM,
                  sample_matMM);
        // Orientation = partner's orientation
        partner->sexual_orientation = agent->sexual_orientation;
        // Partner sex
        if (agent->sexual_orientation == HETEROSEXUAL) {
          if (agent->sex == MALE)
            partner->sex = FEMALE;
          else
            partner->sex = MALE;
        } else {
          partner->sex = agent->sex;
        }
        partner->age = agent->desired_age;
        // Partner in relationship
        partner->initial_relationship = true;
        // Preferred age of partner
        partner->desired_age = agent->age;
        // partner infection risk parameters
        setInitialInfection(*partner, initialInfectionRatesMSW,
                            initialInfectionRatesMSM, initialInfectionRatesWSM,
                            initialInfectionRatesWSW);
        setAgentInfectionParameters(*partner,
                                    het_male_infectivity,
                                    het_female_infectivity,
                                    het_male_infectiousness,
                                    het_female_infectiousness,
                                    hom_male_infectivity,
                                    hom_female_infectivity,
                                    hom_male_infectiousness,
                                    hom_female_infectiousness);
        makePartner(agent, partner, distance(agent, partner),
                    meanDaysRelationship);
        // Correct relationship time because this is in the middle of relationship
        std::uniform_real_distribution<double>
          uni2(currentDate, agent->relationship_change_date);
        agent->relationship_change_date = uni2(rng);
        partner->relationship_change_date = agent->relationship_change_date;
      } else {
        agent->setSingleLength(currentDate, meanDaysSingle);
        // Correct single time because this is in the middle of relationship
        std::uniform_real_distribution<double>
          uni2(currentDate, agent->relationship_change_date);
        agent->relationship_change_date = uni2(rng);
        agent->partner = NULL;
      }

      agents.push_back(agent);
      if (agent->partner)
        agents.push_back(agent->partner);
    }
  }

  void initializeSimulation()
  {
    // Age distribution matrices
    msmAgeDist = CSVParser(parameterMap.at("MSM_AGE_DIST_CSV").stringValue.
                           c_str(), ",", true).convert_all_entries_to_doubles();
    wswAgeDist = CSVParser(parameterMap.at("WSW_AGE_DIST_CSV").stringValue.
                           c_str(), ",", true).convert_all_entries_to_doubles();
    mswAgeDist = CSVParser(parameterMap.at("MSW_AGE_DIST_CSV").stringValue.
                           c_str(), ",", true).convert_all_entries_to_doubles();
    wsmAgeDist = CSVParser(parameterMap.at("WSM_AGE_DIST_CSV").stringValue.
                           c_str(), ",", true).convert_all_entries_to_doubles();
    setEvents();
  }


  void trackRiskFactors(Agent* agent)
  {
    if (agent->sex == MALE) {
      ++numInfectedMales;
      ++infectedMalesByAge[ (int) agent->age / ageInterval];
      if (agent->sexual_orientation == HETEROSEXUAL)
        ++numInfectedMsw;
      else
        ++numInfectedMsm;
    } else {
      ++numInfectedFemales;
      ++infectedFemalesByAge[ (int) agent->age / ageInterval];
      if (agent->sexual_orientation == HETEROSEXUAL)
        ++numInfectedWsm;
      else
        ++numInfectedWsw;
    }

    if (agent->binomial_p_relationship_length >= 0.5 &&
        agent->binomial_p_relationship_wait >= 0.5) {
      ++numInfectedLongBreakLongPartnership;
    } else if (agent->binomial_p_relationship_length >= 0.5 &&
               agent->binomial_p_relationship_wait < 0.5) {
      ++numInfectedShortBreakLongPartnership;
    } else if (agent->binomial_p_relationship_length < 0.5 &&
               agent->binomial_p_relationship_wait >= 0.5) {
      ++numInfectedLongBreakShortPartnership;
    } else {
      ++numInfectedShortBreakShortPartnership;
    }
  }

  AgentVector getUnmatchedAgents()
  {
    AgentVector matingPool;
    for (auto& agent : agents) {
      if (agent->isMatchable(currentDate))
        matingPool.push_back(agent);
    }
    return matingPool;
  }

  AgentVector getShuffledUnmatchedAgents()
  {
    AgentVector matingPool = getUnmatchedAgents();
    shuffle(matingPool.begin(), matingPool.end(), rng);
    // Remove back agent if odd
    if (matingPool.size() % 2 == 1) matingPool.pop_back();
    return matingPool;
  }

  PartnershipScore
  closestPairMatchN(std::vector<Agent *>::iterator from,
		    std::vector<Agent *>::iterator to)
  {
    double smallest_val = DBL_MAX;
    std::vector<Agent *>::iterator closest_agent = to;

    for (auto it = from + 1; it != to; ++it) {
      if ( (*it)->partner == NULL) {
        double d = distance(*from, *it);
        if (d < smallest_val) {
          smallest_val = d;
          closest_agent = it;
        }
      }
    }
    PartnershipScore partnershipScore = {closest_agent, smallest_val};
    return partnershipScore;
  }


  void analysis()
  {
    double prevalence = (double) (numInfectedMales + numInfectedFemales) /
      (numMales + numFemales);
    double malePrevalence = (double) numInfectedMales / numMales;
    double femalePrevalence = (double) numInfectedFemales / numFemales;
    double msmPrevalence = (double) numInfectedMsm / numMsm;
    double wswPrevalence = (double) numInfectedWsw / numWsw;
    printf("%s,ANALYSIS,INFECTED,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate,
           numInfectedMales + numInfectedFemales);
    printf("%s,ANALYSIS,PREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(), simulationNum, currentDate, prevalence);
    printf("%s,ANALYSIS,MALEPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, malePrevalence);
    printf("%s,ANALYSIS,FEMALEPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, femalePrevalence);
    printf("%s,ANALYSIS,MSMPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, msmPrevalence);
    printf("%s,ANALYSIS,WSWPREVALENCE,%d,%.3f,%f\n",
           simulationName.c_str(),simulationNum, currentDate, wswPrevalence);
    printf("%s,ANALYSIS,PARTNERSHIPS,%d,%.3f,%u\n",
           simulationName.c_str(),simulationNum, currentDate, totalPartnerships);
    printf("%s,ANALYSIS,MSMPARTNERSHIPS,%d,%.3f,%u\n",
           simulationName.c_str(),simulationNum, currentDate,
           totalMsmPartnerships);
    printf("%s,ANALYSIS,WSWPARTNERSHIPS,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate,
           totalWswPartnerships);

    double rateShortBreakShortPartnership = (double)
      numInfectedShortBreakShortPartnership / numShortBreakShortPartnership;
    double rateShortBreakLongPartnership = (double)
      numInfectedShortBreakLongPartnership /  numShortBreakLongPartnership;
    double rateLongBreakShortPartnership = (double)
      numInfectedLongBreakShortPartnership / numLongBreakShortPartnership;
    double rateLongBreakLongPartnership = (double)
      numInfectedLongBreakLongPartnership / numLongBreakLongPartnership;
    printf("%s,ANALYSIS,SHORTBREAK_SHORTPARTNERSHIP,%d,%.3f,%.3f\n",
           simulationName.c_str(), simulationNum, currentDate,
           rateShortBreakShortPartnership);
    printf("%s,ANALYSIS,SHORTBREAK_LONGPARTNERSHIP,%d,%.3f,%.3f\n",
           simulationName.c_str(), simulationNum, currentDate,
           rateShortBreakLongPartnership);
    printf("%s,ANALYSIS,LONGBREAK_SHORTPARTNERSHIP,%d,%.3f,%.3f\n",
           simulationName.c_str(), simulationNum, currentDate,
           rateLongBreakShortPartnership);
    printf("%s,ANALYSIS,LONGBREAK_LONGPARTNERSHIP,%d,%.3f,%.3f\n",
           simulationName.c_str(), simulationNum, currentDate,
           rateLongBreakLongPartnership);

    for (unsigned i = 0; i < malesByAge.size(); ++i) {
      ostringstream ssmale, ssfemale;
      if (malesByAge[i] == 0)
        ssmale << "NA";
      else
        ssmale << std::setprecision(3)
               << ( (double) infectedMalesByAge[i] / malesByAge[i] );
      if (femalesByAge[i] == 0)
        ssfemale << "NA";
      else
        ssfemale << std::setprecision(3)
                 << ((double) infectedFemalesByAge[i] / femalesByAge[i]);
      printf("%s,ANALYSIS,MALE_AGE_%03u-%03u,%d,%.3f,%s\n",
             simulationName.c_str(),
             i * ageInterval, i * ageInterval + ageInterval - 1,
             simulationNum, currentDate, ssmale.str().c_str());
      printf("%s,ANALYSIS,FEMALE_AGE_%03u-%03u,%d,%.3f,%s\n",
             simulationName.c_str(),
             i * ageInterval, i * ageInterval + ageInterval - 1,
             simulationNum, currentDate, ssfemale.str().c_str());
    }
    printf("%s,ANALYSIS,SCORE,%d,%.3f,%f\n",
           simulationName.c_str(), simulationNum, currentDate,
           totalPartnershipScore / totalPartnerships);
    printf("%s,ANALYSIS,FAILED,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate, failedMatches);
    printf("%s,ANALYSIS,POOR,%d,%.3f,%u\n",
           simulationName.c_str(), simulationNum, currentDate, poorMatches);
  }


  void makePartner(Agent* a, Agent *b,
                   const double score,
                   const double mean_days_relationship)
  {
    assert(a->partner == NULL);
    assert(b->partner == NULL);

    if (score < failureThresholdScore) {
      if (score > poorThresholdScore) ++poorMatches;
      totalPartnershipScore += score;
      partnerships.insert(a->id, b->id);
      a->partner = b;
      b->partner = a;
      a->setRelationshipLength(currentDate, mean_days_relationship);
      ++a->num_partners;
      ++b->num_partners;
      ++totalPartnerships;
      if (a->sex == b->sex) {
        if (a->sex == MALE)
          ++totalMsmPartnerships;
        else
          ++totalWswPartnerships;
      } else {
        ++totalMswPartnerships;
      }
    } else {
      ++failedMatches;
    }
  }


  double distance(const Agent *a, const Agent *b) const
  {
    return distanceMethod == HEURISTIC_DISTANCE
      ? heuristicDistance(a, b) : tableDistance(a, b);
  }

  double heuristicDistance(const Agent *a, const Agent *b) const
  {
    double score = 0.0;

    score +=  (fabs(a->desired_age - b->age) +
               fabs(b->desired_age - a->age)) / 2.0;
    if (a->sexual_orientation != b->sexual_orientation) {
      score += 50.0;
    } else if (a->sexual_orientation == HETEROSEXUAL) {
      if (a->sex == b->sex)
        score += 50.0;
    } else if (a->sex != b->sex) {
      score += 50.0;
    }
    if (partnerships.exists(a->id, b->id))
        score += 50.0;
    return score;
  }

  double tableDistance(const Agent *a, const Agent *b) const
  {
    double score = 0.0;
    unsigned a_age = a->age - MIN_AGE;
    unsigned b_age = b->age - MIN_AGE;

    if (a->sex == MALE and b->sex == FEMALE)
      score += (mswAgeDist[a_age][b_age] + wsmAgeDist[b_age][a_age]) * 25;
    else if (a->sex == FEMALE and b->sex == MALE)
      score += (wsmAgeDist[a_age][b_age] + mswAgeDist[b_age][a_age]) * 25;
    else if (a->sex == MALE and b->sex == MALE)
      score += (msmAgeDist[a_age][b_age] + msmAgeDist[b_age][a_age]) * 25;
    else if (a->sex == FEMALE and b->sex == FEMALE)
      score += (wswAgeDist[a_age][b_age] + wswAgeDist[b_age][a_age]) * 25;
    else {
      fprintf(stderr, "Error: Unknown sex combination.\n");
      exit(1);
    }

    if (a->sexual_orientation != b->sexual_orientation) {
      score += 50.0;
    } else if (a->sexual_orientation == HETEROSEXUAL) {
      if (a->sex == b->sex)
        score += 50.0;
    } else if (a->sex != b->sex) {
      score += 50.0;
    }

    if (partnerships.exists(a->id, b->id))
        score += 50.0;

    return score;
  }

  double clusterValue(const Agent *a) const
  {
    return ( (double) a->desired_age / 100.0 + a->age / 100.0) / 2.0
      + a->sexual_orientation;
  }

  void setEvents();

  std::string simulationName;
  AgentVector agents;
  Partnerships partnerships;
  ParameterMap parameterMap;
  unsigned simulationNum;
  unsigned ageInterval;
  unsigned distanceMethod;
  double startDate;
  double endDate;
  double timeStep;
  double meanDaysRelationship;
  double meanDaysSingle;
  double currentDate;
  double failureThresholdScore;
  double poorThresholdScore;
  double totalPartnershipScore = 0.0;
  unsigned clusters;
  unsigned neighbors;
  unsigned totalBreakups = 0;
  unsigned totalPartnerships = 0;
  unsigned totalMswPartnerships = 0;
  unsigned totalMsmPartnerships = 0;
  unsigned totalWswPartnerships = 0;
  unsigned numMales = 0;
  unsigned numFemales = 0;
  unsigned numMsm = 0;
  unsigned numMsw = 0;
  unsigned numWsm = 0;
  unsigned numWsw = 0;
  unsigned numInfectedMales = 0;
  unsigned numInfectedFemales = 0;
  unsigned numInfectedMsm = 0;
  unsigned numInfectedMsw = 0;
  unsigned numInfectedWsm = 0;
  unsigned numInfectedWsw = 0;
  unsigned numShortBreakShortPartnership = 0;
  unsigned numShortBreakLongPartnership = 0;
  unsigned numLongBreakShortPartnership = 0;
  unsigned numLongBreakLongPartnership = 0;
  unsigned numInfectedShortBreakShortPartnership = 0;
  unsigned numInfectedShortBreakLongPartnership = 0;
  unsigned numInfectedLongBreakShortPartnership = 0;
  unsigned numInfectedLongBreakLongPartnership = 0;
  unsigned poorMatches = 0;
  unsigned failedMatches = 0;
  std::vector<unsigned> malesByAge;
  std::vector<unsigned> femalesByAge;
  std::vector<unsigned> infectedMalesByAge;
  std::vector<unsigned> infectedFemalesByAge;
  DblMatrix mswAgeDist;
  DblMatrix wsmAgeDist;
  DblMatrix msmAgeDist;
  DblMatrix wswAgeDist;
  std::vector< std::function<void(Simulation*)> > events;
};


void ageEvent(Simulation* simulation)
{
  for (auto& agent : simulation->agents) agent->age += simulation->timeStep;
}

void infectEvent(Simulation* simulation)
{
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (auto& agent: simulation->agents) {
    if (agent->partner && agent->infected == false &&
        agent->partner->infected == true) {
      double risk_infection;
      if (agent->partner->sex == agent->sex)
        risk_infection = (agent->hom_infectivity +
                          agent->partner->hom_infectiousness) / 2.0;
      else
        risk_infection = (agent->het_infectivity +
                          agent->partner->het_infectiousness) / 2.0;
      if (dist(rng) < risk_infection) {
        agent->infected = true;
        agent->infector = agent->partner;
        ++agent->partner->num_infected;
        simulation->trackRiskFactors(agent);
      }
    }
  }
}

void breakupEvent(Simulation* simulation)
{
  unsigned breakups = 0;
  for (auto& agent: simulation->agents) {
    if (agent->partner && simulation->currentDate >
        agent->relationship_change_date) {
      Agent* partner = agent->partner;
      agent->partner = NULL;
      partner->partner = NULL;
      agent->setSingleLength(simulation->currentDate,
                             simulation->meanDaysSingle);
      partner->setSingleLength(simulation->currentDate,
                               simulation->meanDaysSingle);
      ++breakups;
    }
  }
  simulation->totalBreakups += breakups;
}

void randomMatchEvent(Simulation* simulation)
{
  // struct timeval timeBegin, timeEnd;
  // double elapsedTime;
  // gettimeofday(&timeBegin, NULL);

  AgentVector unmatchedAgents = simulation->getShuffledUnmatchedAgents();

  if (unmatchedAgents.size() > 0)
    for (size_t i = 0; i < unmatchedAgents.size()  - 1; i += 2)
      simulation->makePartner(unmatchedAgents[i], unmatchedAgents[i + 1],
                              simulation->distance(unmatchedAgents[i],
                                       unmatchedAgents[i + 1]),
                              simulation->meanDaysRelationship);

  // gettimeofday(&timeEnd, NULL);
  // elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
  // printf("%s,PM_TIMING,RPM,%u,%f\n", simulation->simulationName.c_str(),
  //          simulation->simulationNum, elapsedTime);
  // printf("%s,PM_SIZE,RPM,%u,%lu\n", simulation->simulationName.c_str(),
  //        simulation->simulationNum, unmatchedAgents.size());
}

void randomKMatchEvent(Simulation *simulation)
{
  // struct timeval timeBegin, timeEnd;
  // double elapsedTime;
  // gettimeofday(&timeBegin, NULL);


  std::uniform_real_distribution<double> uni;
  AgentVector unmatchedAgents = simulation->getShuffledUnmatchedAgents();

  if(unmatchedAgents.size()) {
    for (auto it = unmatchedAgents.begin(); it < unmatchedAgents.end() - 1;
         ++it) {
      if ( (*it)->partner == NULL) {
        auto last =  (unmatchedAgents.end() - it) < (simulation->neighbors + 1)
          ? unmatchedAgents.end() : it + simulation->neighbors + 1;
        auto partnershipScore = simulation->closestPairMatchN(it, last);
        if (partnershipScore.partner != last)
          simulation->makePartner(*it, *partnershipScore.partner,
                                  partnershipScore.score,
                                  simulation->meanDaysRelationship);
      }
    }
  }

  // gettimeofday(&timeEnd, NULL);
  // elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
  // printf("%s,PM_TIMING,RKPM,%u,%f\n", simulation->simulationName.c_str(),
  //          simulation->simulationNum, elapsedTime);
  // printf("%s,PM_SIZE,RKPM,%u,%lu\n", simulation->simulationName.c_str(),
  //        simulation->simulationNum, unmatchedAgents.size());
}

void clusterShuffleMatchEvent(Simulation* simulation)
{
  // struct timeval timeBegin, timeEnd;
  // double elapsedTime;
  // gettimeofday(&timeBegin, NULL);

  AgentVector unmatchedAgents = simulation->getShuffledUnmatchedAgents();
  uint64_t cluster_size = unmatchedAgents.size() / simulation->clusters;
  for (auto& a : unmatchedAgents) a->weight = simulation->clusterValue(a);
  sort(unmatchedAgents.rbegin(), unmatchedAgents.rend(), [](Agent *a, Agent *b)
       { return a->weight < b->weight; });
  for (uint64_t i = 0; i < simulation->clusters; ++i) {
    auto first = unmatchedAgents.begin() + i * cluster_size;
    auto last = first + cluster_size;
    if (last > unmatchedAgents.end()) last = unmatchedAgents.end();
    std::shuffle(first, last, rng);
  }
  if(unmatchedAgents.size()) {
    for (auto it = unmatchedAgents.begin(); it < unmatchedAgents.end() - 1;
         ++it) {
      if ( (*it)->partner == NULL) {
        auto last = (unmatchedAgents.end() - it) < (simulation->neighbors + 1)
          ? unmatchedAgents.end() : it + simulation->neighbors + 1;
        auto partnershipScore = simulation->closestPairMatchN(it, last);
        if (partnershipScore.partner != last)
          simulation->makePartner(*it, *partnershipScore.partner,
                      partnershipScore.score, simulation->meanDaysRelationship);
      }
    }
  }
  // gettimeofday(&timeEnd, NULL);
  // elapsedTime = timeEnd.tv_sec - timeBegin.tv_sec;
  // printf("%s,PM_TIMING,CSPM,%u,%f\n", simulation->simulationName.c_str(),
  //          simulation->simulationNum, elapsedTime);
  // printf("%s,PM_SIZE,CSPM,%u,%lu-%lu\n", simulation->simulationName.c_str(),
  //        simulation->simulationNum, unmatchedAgents.size(),
  //        unmatchedAgents.capacity());
  // auto &p = simulation->partnerships.partnerships;
  // printf("%s,PM_SET,CSPM,%u,%.2f-%.2f-%lu-%lu\n",
  //        simulation->simulationName.c_str(),
  //        simulation->simulationNum,
  //        p.load_factor(), p.max_load_factor(),
  //        p.bucket_count(), p.max_bucket_count());

}


/******************** Call Blossom V reference algorithm ******************/

void graphPairs(const char *graph,
                const Simulation* simulation,
                AgentVector& agents)
{
  FILE *f = fopen(graph, "w");

  uint64_t vertices = (uint64_t) agents.size();
  uint64_t edges = vertices * (vertices - 1) / 2;
  fprintf(f, "%lu %lu\n", vertices, edges);
  for (uint64_t i = 0; i < agents.size(); ++i) {
    for (uint64_t j = i + 1; j < agents.size(); ++j) {
      double d = simulation->distance(agents[i],agents[j]);
      fprintf(f, "%lu %lu %.0f\n", i, j, d * GRAPH_ACCURACY);
    }
  }
  fclose(f);
}

void blossomVMatchEvent(Simulation* simulation)
{
  AgentVector agents = simulation->getShuffledUnmatchedAgents();

  if (agents.size() == 0) return;

  std::ostringstream ss;
  ss << std::this_thread::get_id();
  std::string suffix = ss.str() + std::string(".txt");
  std::string graph_file = std::string("bv_graph_") + suffix;
  std::string blossom_out_file = std::string("bv_out_") + suffix;

  graphPairs(graph_file.c_str(), simulation, agents);
  std::string command("./blossom5 -e ");
  command += graph_file;
  command +=  std::string(" -w ");
  command += blossom_out_file;
  command += std::string(" > /dev/null");
  fflush(stdout);
  if(system(command.c_str()) != 0) {
    std::cerr << "Error executing Blossom V." << std::endl;
    exit(1);
  }

  FILE *f = fopen(blossom_out_file.c_str(), "r");
  uint64_t from, to;
  if (fscanf(f, "%lu %lu\n", &from, &to) != 2) {
    fprintf(stderr, "Error reading graph file.\n");
    exit(1);
  }
  for (size_t i = 0; i < agents.size() / 2; ++i) {
    if (fscanf(f, "%lu %lu\n", &from, &to) != 2) {
      fprintf(stderr, "Error reading graph file.\n");
      exit(1);
    }
    simulation->makePartner(agents[from], agents[to],
                            simulation->distance(agents[from], agents[to]),
                            simulation->meanDaysRelationship);
  }
  fclose(f);
}

/*********************************************************************/

void Simulation::setEvents()
{
  if (parameterMap.at("AGE_EVENT").isSet()) events.push_back(ageEvent);
  if (parameterMap.at("INFECT_EVENT").isSet()) events.push_back(infectEvent);
  if (parameterMap.at("BREAKUP_EVENT").isSet()) events.push_back(breakupEvent);

  string s = parameterMap.at("MATCH_EVENT").getStr();
  if (s == "RPM") {
    events.push_back(randomMatchEvent);
  } else if (s == "RKPM") {
    events.push_back(randomKMatchEvent);
  } else if (s == "CSPM") {
    events.push_back(clusterShuffleMatchEvent);
  }  else if (s == "BLOSSOMV") {
    events.push_back(blossomVMatchEvent);
  } else {
    fprintf(stderr, "Unknown matching algorithm.\n");
    exit(1);
  }
}

/***************************/

void runTests(ParameterMap& parameterMap)
{
  unsigned successes = 0, failures = 0;

  // Partnerships
  Partnerships partnerships;
  partnerships.insert(12995, 271);
  partnerships.insert(3, 23994);
  TESTEQ(partnerships.exists(12995, 271), true, successes, failures);
  TESTEQ(partnerships.exists(271, 12995), true, successes, failures);
  TESTEQ(partnerships.exists(3, 23994), true, successes, failures);
  TESTEQ(partnerships.exists(23994, 3), true, successes, failures);
  TESTEQ(partnerships.exists(30, 23994), false, successes, failures);
  TESTEQ(partnerships.exists(23994, 30), false, successes, failures);
  TESTEQ(partnerships.exists(12995, 272), false, successes, failures);
  TESTEQ(partnerships.size() == 2, true, successes, failures);

  Simulation simulation(parameterMap, 0);
  simulation.initializeSimulation();
  simulation.initializeAgents();

  TESTEQ(simulation.agents.size(), parameterMap.at("NUM_AGENTS").getDbl(),
         successes, failures);

  simulation.agents[0]->age = 20;
  simulation.agents[0]->sex = MALE;
  simulation.agents[0]->desired_age = 50;
  simulation.agents[0]->sexual_orientation = HETEROSEXUAL;

  simulation.agents[1]->age = 30;
  simulation.agents[1]->sex = MALE;
  simulation.agents[1]->desired_age = 25;
  simulation.agents[1]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[2]->age = 40;
  simulation.agents[2]->sex = FEMALE;
  simulation.agents[2]->desired_age = 30;
  simulation.agents[2]->sexual_orientation = HETEROSEXUAL;

  simulation.agents[3]->age = 25;
  simulation.agents[3]->sex = FEMALE;
  simulation.agents[3]->desired_age = 25;
  simulation.agents[3]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[4]->age = 30;
  simulation.agents[4]->sex = FEMALE;
  simulation.agents[4]->desired_age = 20;
  simulation.agents[4]->sexual_orientation = HOMOSEXUAL;

  simulation.agents[5]->age = 30;
  simulation.agents[5]->sex = MALE;
  simulation.agents[5]->desired_age = 20;
  simulation.agents[5]->sexual_orientation = HOMOSEXUAL;

  double d1 = simulation.heuristicDistance(simulation.agents[0],
                                          simulation.agents[2]);
  TESTEQ(d1, 10.0, successes, failures);

  double d2 = simulation.heuristicDistance(simulation.agents[2],
                                           simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[0],
                                   simulation.agents[3]);
  TESTEQ(d1, 65.0, successes, failures);
  d2 = simulation.heuristicDistance(simulation.agents[3],
                                    simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[0],
                                   simulation.agents[1]);
  TESTEQ(d1, 62.5, successes, failures);
  d2 = simulation.heuristicDistance(simulation.agents[1],
                                   simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);
  d1 = simulation.heuristicDistance(simulation.agents[1],
                                   simulation.agents[2]);
  TESTEQ(d1, 57.5, successes, failures);
  d2 = simulation.heuristicDistance(simulation.agents[2],
                                    simulation.agents[1]);
  TESTEQ(d1, d2, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[1],
                                   simulation.agents[3]);
  TESTEQ(d1, 52.5, successes, failures);
  d1 = simulation.heuristicDistance(simulation.agents[2],
                                   simulation.agents[3]);
  TESTEQ(d1, 60, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[3],
                                   simulation.agents[4]);
  TESTEQ(d1, 5, successes, failures);

  d1 = simulation.heuristicDistance(simulation.agents[1],
                                   simulation.agents[5]);
  TESTEQ(d1, 7.5, successes, failures);

  d1 = simulation.tableDistance(simulation.agents[0],
                               simulation.agents[3]);
  TESTEQ(fabs(d1 - 92.85787107631695) < 0.00000000001, true, successes, failures);
  d2 = simulation.tableDistance(simulation.agents[3],
                                simulation.agents[0]);
  TESTEQ(d1, d2, successes, failures);
  d1 = simulation.tableDistance(simulation.agents[3],
                                simulation.agents[4]);
  TESTEQ(fabs(d1 - 21.9126355725) < 0.00000001, true, successes, failures);
  d2 = simulation.tableDistance(simulation.agents[4],
                                simulation.agents[3]);
  TESTEQ(d1, d2, successes, failures);
  d1 = simulation.tableDistance(simulation.agents[1],
                                simulation.agents[5]);
  TESTEQ(fabs(d1 - 0.644530208167975) < 0.00000001, true, successes, failures);
  d2 = simulation.tableDistance(simulation.agents[5],
                                simulation.agents[1]);
  TESTEQ(d1, d2, successes, failures);
  printf("Successes: %u. Failures: %u.\n", successes, failures);
}

/***************************/

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

void callSimulation(ParameterMap& parameterMap, const unsigned simulationNum)
{

  Simulation(parameterMap, simulationNum).simulate();
}


int main(int argc, char *argv[])
{
  if (cmdOptionExists(argv, argv + argc, "-h")) {
    printf("%s options, where options are:\n"
           "-h: help - Print this message.\n"
           "-f filename: Use filename as input.\n"
           "-s integer: Random seed (0 - use time).\n"
           "-t: run tests.\n",
           argv[0]);
    ParameterMap parameterMap;
    setDefaultParameters(parameterMap);
    printParameters(parameterMap);
    exit(1);
  }

  const char *seed_str = getCmdOption(argv, argv + argc, "-s");
  std::vector<ParameterMap> parameterMaps;
  const char *input_file_str = getCmdOption(argv, argv + argc, "-f");
  if (input_file_str) {
    std::ifstream infile;
    infile.open (input_file_str, std::ifstream::in);
    if (infile.fail()) {
      fprintf(stderr, "Error opening %s\n", input_file_str);
      exit(1);
    }
    readParameters(infile, parameterMaps);
    infile.close();
  }

  if (cmdOptionExists(argv, argv + argc, "-t")) {
    ParameterMap parameterMap;
    setDefaultParameters(parameterMap);
    runTests(parameterMap);
  }
  for (auto& parameterMap: parameterMaps) {
    if (seed_str) {
      unsigned seed = atoi(seed_str);
      if (seed == 0) // Use time
        parameterMap["RANDOM_SEED"].values = { (double) time(NULL) };
      else
        parameterMap["RANDOM_SEED"].values = { (double) seed};
    }
    if (parameterMap.at("PRINT_PARAMETERS").getDbl() == 1)
      printParameters(parameterMap);
    else {
      unsigned numThreads = parameterMap.at("NUM_THREADS").getDbl();
      unsigned numSimulations = parameterMap.at("NUM_SIMULATIONS").getDbl();
      if (numThreads == 0) numThreads = numSimulations;
      unsigned simulationsRun = 0;
      while(simulationsRun < numSimulations) {
        if (simulationsRun + numThreads > numSimulations)
          numThreads = numSimulations - simulationsRun;
        std::thread t[numThreads];
        for (unsigned i = 0; i < numThreads; ++i)
          t[i] = std::thread(callSimulation, std::ref(parameterMap),
                             simulationsRun + i);
        for (unsigned i = 0; i < numThreads; ++i)
          t[i].join();
        simulationsRun += numThreads;
      }
    }
  }
}
