#include <cstdio>
#include <ctime>

#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <random>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "CSV_Parser.hh"
#include "sample.hh"

thread_local std::mt19937 rng;

#define DAY 1.0 / 365;
#define MAX_AGE 100

#define BFPM 1
#define RPM 2
#define CSPM 3

#define FEMALE 0
#define MALE 1
#define HETEROSEXUAL 0
#define HOMOSEXUAL 1



// Linear algebra functions for Stefan's initial population creation algorithm.

std::vector<double> getCol(const DblMatrix& matrix, unsigned col)
{
  std::vector<double> output;
  for (unsigned i = 0; i < matrix.size(); ++i)
    output.push_back(matrix[i][col]);
  return output;
}

std::vector<double> subVector(double d,
			       const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d - v[i]);
  return output;
}

std::vector<double> multVector(double d,
				const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d * v[i]);
  return output;
}

std::vector<double> multVectors(const std::vector<double>& v1,
				 const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] * v2[i]);
  return output;
}


std::vector<double> addVector(const std::vector<double>& v1,
				const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] + v2[i]);
  return output;
}

double sumVector(const std::vector<double>& dbls)
{
  double total = 0.0;
  for (auto d: dbls)
    total += d;
  return total;
}

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


// End of linear algebra functions
////////////////////////////////////////


///////////////////////////////////////////////
// Parameter classes and functions

struct ParameterValue {
  double getDbl(size_t index = 0) const
  {
    return values[0];
  };
  std::string getStr() const
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
  addParameter(parameterMap, "MATCHING_ALGORITHM",
               "1 for BFPM, 2 for RPM, 3 for CSPM", {BFPM});
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
               "(0=no,1=yes)", {1.0});

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

  addParameter(parameterMap, "AGE_EVENT", "Execute the aging event", {0});
  addParameter(parameterMap, "INFECT_EVENT", "Execute the infection event", {1});
  addParameter(parameterMap, "BREAKUP_EVENT", "Execute the infection event", {1});
  addParameter(parameterMap, "RANDOM_MATCH_EVENT",
               "Execute the infection event", {1});
}

void readParameters(std::istream& input, ParameterMap& parameterMap)
{
  std::string lineString;
  unsigned lineNumber = 0;
  std::string errorMessage;
  while (std::getline(input, lineString)) {
    ++lineNumber;
    boost::trim(lineString);
    if (lineString.size() == 0) continue;
    if (lineString[0] == '#') continue;
    std::istringstream line(lineString);
    std::string parameterName;
    line >> parameterName;
    if (replaceParameter(parameterMap, parameterName.c_str(),
                         line, errorMessage) == false) {
      fprintf(stderr,"Error at line: %u\n", lineNumber);
      fprintf(stderr, "%s\n", errorMessage.c_str());
    }
  }
}

// End of Parameter management

////////////////////////////////////
// Agent class and related functions

class Agent;
typedef std::vector<Agent *> AgentVector;

class Agent {
public:
  unsigned id = 0;
  Agent* partner;
  AgentVector past_partners;
  double het_infectivity;
  double het_infectiousness;
  double hom_infectivity;
  double hom_infectiousness;
  bool infected;
  bool initial_relationship;
  double age;
  unsigned sex;
  unsigned sexual_orientation;
  unsigned desired_age;
  double relationship_change_date;
  double binomial_p_relationship_length;
  double binomial_p_relationship_wait;

  void setRelationshipLength(const double current_date,
                             const double mean_days_relationship)
  {
    double p = (binomial_p_relationship_length +
                partner->binomial_p_relationship_length) / 2.0;
    std::binomial_distribution<int> dist (mean_days_relationship, p);
    relationship_change_date = current_date + (double) dist(rng) * DAY;
    partner->relationship_change_date = relationship_change_date;
  }
  void setSingleLength(double current_date,
                       double mean_days_single)
  {
    std::binomial_distribution<int> dist (mean_days_single,
                                          binomial_p_relationship_wait);
    relationship_change_date = current_date + (double) dist(rng) * DAY;
  }

  bool isMatchable(const double current_date) const
  {
    if (partner == NULL && current_date >= relationship_change_date)
      return true;
    else
      return false;
  }

  void makePartner(Agent* agent,
                   const double current_date,
                   const double mean_days_relationship)
  {
    assert(partner == NULL);
    assert(agent->partner == NULL);
    partner = agent;
    agent->partner = this;
    setRelationshipLength(current_date, mean_days_relationship);
  }

  void print(FILE *f = stdout)
  {
    fprintf(f, "ID,%7u,Age,%3.2f,Sex,%c,Orientation,%c,Desired,%u",
            id, age, (sex == MALE ? 'M' : 'F'),
            (sexual_orientation == HETEROSEXUAL ? 'S' : 'G'), desired_age);
    fprintf(f, ",Risk,%.3f,%.3f,Partner",
            binomial_p_relationship_length,binomial_p_relationship_wait);
    if (partner)
      fprintf(f,",%7u", partner->id);
      // auto &p = partner;
      // fprintf(f, ",ID,%u,Age,%.2f,Sex,%c,Orientation,%c,Desired,%u",
      //         p->id, p->age, (p->sex == MALE ? 'M' : 'F'),
      //         (p->sexual_orientation == HETEROSEXUAL ? 'S' : 'G'),
      //         p->desired_age);
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

// End of Agents
////////////////////////////////////////

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

// Simulation class and supporting routines

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
                         const std::vector<double>& initialInfectionRatesMale,
                         const std::vector<double>& initialInfectionRatesFemale)
{
  std::uniform_real_distribution<double> uni;
  agent.infected = (agent.sex == MALE) ?
    (uni(rng) < initialInfectionRatesMale[(int) agent.age] ? true : false) :
    (uni(rng) < initialInfectionRatesFemale[(int) agent.age] ? true : false);
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
                    const std::vector<double>& initialInfectionRatesMale,
                    const std::vector<double>& initialInfectionRatesFeMale,
                    bool initial_relation)
{
  std::uniform_real_distribution<double> uni;
  double start_date = parameterMap.at("START_DATE").getDbl();
  double mean_days_relationship =
    parameterMap.at("MEAN_DAYS_RELATIONSHIP").getDbl();
  double mean_days_single =
    parameterMap.at("MEAN_DAYS_SINGLE").getDbl();
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
    // ID
    agent->id = i + 1;
    // Age
    unsigned age = sample_ageshare() + 12;
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
    setInitialInfection(*agent, initialInfectionRatesMale,
                        initialInfectionRatesFeMale);
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
      // ID of partner
      agent->partner = partner;
      // Partner ID
      partner->id = i + 2;
      // Orientation = partner's orientation
      partner->sexual_orientation = sexual_orientation;
      // Partner sex
      if (sexual_orientation == HETEROSEXUAL) {
        if (sex == MALE)
          partner->sex = FEMALE;
        else
          partner->sex = MALE;
      } else {
        partner->sex = sex;
      }
      partner->age = agent->desired_age;
      // Partner in relationship
      partner->initial_relationship = true;
      // Partner's partner
      partner->partner = agent;
      // Preferred age of partner
      partner->desired_age = agent->age;
      // partner infection risk parameters
      setInitialInfection(*partner, initialInfectionRatesMale,
                          initialInfectionRatesFeMale);
      setAgentInfectionParameters(*partner,
                                  het_male_infectivity,
                                  het_female_infectivity,
                                  het_male_infectiousness,
                                  het_female_infectiousness,
                                  hom_male_infectivity,
                                  hom_female_infectivity,
                                  hom_male_infectiousness,
                                  hom_female_infectiousness);

    } else {
      agent->partner = NULL;
    }
    if (agent->partner)
      agent->setRelationshipLength(start_date, mean_days_relationship);
    else
      agent->setSingleLength(start_date, mean_days_single);
    agents.push_back(agent);
    if (agent->partner)
      agents.push_back(agent->partner);
  }
}

class Simulation {
public:
  Simulation(const ParameterMap& parameter_map,
             const unsigned simulation_num) :
    parameterMap(parameter_map), simulationNum(simulation_num)
  {
    startDate = parameterMap.at("START_DATE").getDbl();
    endDate = parameterMap.at("END_DATE").getDbl();
    timeStep = parameterMap.at("TIME_STEP").getDbl();
    meanDaysRelationship = parameterMap.at("MEAN_DAYS_RELATIONSHIP").getDbl();
    meanDaysSingle = parameterMap.at("MEAN_DAYS_SINGLE").getDbl();
    currentDate = startDate;
    ageInterval = parameterMap.at("ANALYSIS_AGE_INTERVAL").getDbl();
  }
  ~Simulation()
  {
    for (auto& agent: agents)
      delete agent;
  }
  void simulate()
  {
    // Events to execute
    bool execAgeEvent = (parameterMap.at("AGE_EVENT").getDbl() == 1)
      ? true : false;
    bool execInfectEvent = (parameterMap.at("INFECT_EVENT").getDbl() == 1)
      ? true : false;
    bool execBreakupEvent = (parameterMap.at("BREAKUP_EVENT").getDbl() == 1)
      ? true : false;
    bool execRandomMatchEvent = (parameterMap.at("RANDOM_MATCH_EVENT").
                                 getDbl() == 1) ? true : false;
    // Initialization
    clock_t t = clock();
    initializeAgents();
    t = clock() - t;
    printf("TIMING,Initialization,%ld,%f\n",t,((float)t)/CLOCKS_PER_SEC);
    unsigned outputAgents = parameterMap.at("OUTPUT_AGENTS_AFTER_INIT").getDbl();
    if (outputAgents)
      printAgents(agents, simulationNum, startDate, stdout);
    unsigned analyzeAgents = parameterMap.
      at("ANALYZE_AFTER_INIT").getDbl();
    if (analyzeAgents)
      analysis();
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

      /* Execute events */
      if (execAgeEvent) ageEvent();
      if (execInfectEvent) infectionEvent();
      if (execBreakupEvent) breakupEvent();
      if (execRandomMatchEvent)
        randomMatchEvent();

      if (outputAgents > 0 && (i + 1) % outputFrequency == 0)
        printAgents(agents, simulationNum, currentDate, stdout);
      if (analyzeAgents > 0 && (i + 1) % analyze_frequency == 0)
        analysis();
    }

    /* Wrap up */
    outputAgents = parameterMap.at("OUTPUT_AGENTS_AT_END").getDbl();
    if (outputAgents)
      printAgents(agents, simulationNum, endDate, stdout);
    analyzeAgents = parameterMap.at("ANALYZE_AT_END").getDbl();
    if (analyzeAgents)
      analysis();
  }
  void initializeAgents()
  {
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
    std::vector<double> initialInfectionRatesMale(MAX_AGE + 1, 0.0);
    std::vector<double> initialInfectionRatesFemale(MAX_AGE + 1, 0.0);

    unsigned i = 0;
    for (auto& row : initialRatesMatrix) {
      for (;i <= row[0] && i < MAX_AGE + 1; ++i) {
        assert(i <= MAX_AGE);
        initialInfectionRatesMale[i] = row[1];
        initialInfectionRatesFemale[i] = row[2];
      }
    }
    for (; i <= MAX_AGE; ++i) {
      initialInfectionRatesMale[i] = 0.0;
      initialInfectionRatesFemale[i] = 0.0;
    }

    unsigned X = parameterMap.at("NUM_AGENTS").values[0];
    unsigned S = calcNumberSingles(data, X);

    std::vector<double> ageRange;
    for (unsigned i = 12; i <= 100; ++i)
      ageRange.push_back(i);
    auto ageShare = getCol(singles, 1);
    auto femRatio = getCol(singles, 2);
    auto msmRate = getCol(singles, 3);
    auto wswRate = getCol(singles, 4);
    createPartners(agents, parameterMap, data, 0, S, ageRange, ageShare, femRatio,
                   wswRate, msmRate, ww, mw, wm, mm, initialInfectionRatesMale,
                   initialInfectionRatesFemale, false);
    ageShare = femRatio = msmRate = wswRate = {};
    ageShare = getCol(partners, 1);
    femRatio = getCol(partners, 2);
    msmRate = getCol(partners, 3);
    wswRate = getCol(partners, 4);
    createPartners(agents, parameterMap, data, S, X, ageRange, ageShare, femRatio,
                   wswRate, msmRate, ww, mw, wm, mm, initialInfectionRatesMale,
                   initialInfectionRatesFemale, true);
    malesByAge.resize( (MAX_AGE + 1) / ageInterval, 0);
    femalesByAge.resize( (MAX_AGE + 1) / ageInterval, 0);
    infectedMalesByAge.resize( (MAX_AGE + 1) / ageInterval, 0);
    infectedFemalesByAge.resize( (MAX_AGE + 1) / ageInterval, 0);
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
        ++femalesByAge[age / ageInterval];
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
  }
  void ageEvent()
  {
    for (auto& agent : agents)
      agent->age += timeStep;
  }
  void infectionEvent()
  {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (auto& agent: agents) {
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
        }
      }
    }
  }

  void breakupEvent()
  {
    for (auto& agent: agents) {
      if (agent->partner && currentDate > agent->relationship_change_date) {
        Agent* partner = agent->partner;
        agent->partner = NULL;
        partner->partner = NULL;
        agent->past_partners.push_back(partner);
        partner->past_partners.push_back(agent);
        agent->setSingleLength(currentDate, meanDaysSingle);
        partner->setSingleLength(currentDate, meanDaysSingle);
      }
    }
  }

  void randomMatchEvent()
  {
    AgentVector unmatchedAgents;
    // Find unmatched agents due for relationship
    for (auto& agent : agents) {
      if (agent->isMatchable(currentDate))
        unmatchedAgents.push_back(agent);
    }
    shuffle(unmatchedAgents.begin(), unmatchedAgents.end(), rng);
    if (unmatchedAgents.size() > 0)
      for (size_t i = 0; i < unmatchedAgents.size()  - 1; i += 2)
        unmatchedAgents[i]->makePartner(unmatchedAgents[i + 1],
                                        currentDate,
                                        meanDaysRelationship);
  }

  void analysis()
  {
    double prevalence = (double) (numInfectedMales + numInfectedFemales) /
      (numMales + numFemales);
    double malePrevalence = (double) numInfectedMales / numMales;
    double femalePrevalence = (double) numInfectedFemales / numFemales;
    double msmPrevalence = (double) numInfectedMsm / numMsm;
    printf("ANALYSIS-TOT,%d,%f,pop,%u,infected,%u,prev,%.3f,"
           "maleprev,%.3f,femaleprev,%.3f,"
           "msmprev,%.3f\n",
           simulationNum, currentDate, (numMales + numFemales),
           (numInfectedMales + numInfectedFemales),
           prevalence, malePrevalence,
           femalePrevalence, msmPrevalence);
    ostringstream ss;
    ss.precision(3);
    for (unsigned i = 0; i < malesByAge.size(); ++i) {
      double malePrev, femalePrev;
      if (malesByAge[i] == 0)
        malePrev = -1.0;
      else
        malePrev = (double) infectedMalesByAge[i] / malesByAge[i];
      if (femalesByAge[i] == 0)
        femalePrev = -1.0;
      else
        femalePrev = (double) infectedFemalesByAge[i] / femalesByAge[i];

      ss << "," << i * ageInterval + ageInterval << ","
         << malePrev << ","
         << femalePrev;
    }
    printf("ANALYSIS-AGE,%d,%f%s\n",
           simulationNum, currentDate, ss.str().c_str());
  }
private:
  AgentVector agents;
  ParameterMap parameterMap;
  unsigned simulationNum;
  unsigned ageInterval;
  double startDate;
  double endDate;
  double timeStep;
  double meanDaysRelationship;
  double meanDaysSingle;
  double currentDate;
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
  vector<unsigned> malesByAge;
  vector<unsigned> femalesByAge;
  vector<unsigned> infectedMalesByAge;
  vector<unsigned> infectedFemalesByAge;
};

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



int main(int argc, char *argv[])
{
  ParameterMap parameterMap;

  setDefaultParameters(parameterMap);
  if (cmdOptionExists(argv, argv + argc, "-h")) {
    printf("%s options, where options are:\n"
           "-h: help - Print this message.\n"
           "-f filename: Use filename as input.\n\n"
           "Parameters\n"
           "**********\n",
           argv[0]);
    printParameters(parameterMap);
    exit(1);
  }

  const char *input_file_str = getCmdOption(argv, argv + argc, "-f");
  if (input_file_str) {
    std::ifstream infile;
    infile.open (input_file_str, std::ifstream::in);
    if (infile.fail()) {
      fprintf(stderr, "Error opening %s\n", input_file_str);
      exit(1);
    }
    readParameters(infile, parameterMap);
    infile.close();
  } else {
    readParameters(std::cin, parameterMap);
  }
  if (parameterMap.at("PRINT_PARAMETERS").getDbl() == 1)
    printParameters(parameterMap);
  Simulation(parameterMap, 0).simulate();
}
