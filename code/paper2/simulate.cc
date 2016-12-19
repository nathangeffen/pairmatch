#include <algorithm>
#include <fstream>
#include <initializer_list>
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

#define BFPM 1
#define RPM 2
#define CSPM 3

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

const unsigned SEXOR_DIM = NUM_AGES;
const unsigned SEX_DIM = NUM_ORIENTATIONS * SEXOR_DIM;
const unsigned AGE_DIM = NUM_SEXES * SEX_DIM;
const unsigned REL_DIM = NUM_AGES * AGE_DIM;
const unsigned NUM_BUCKETS = NUM_RELS * REL_DIM;

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
    std::cout << p.first << "\t" << p.second.description << "\t";
    if (p.second.isString)
      std::cout << p.second.stringValue;
    else
      for (auto& v: p.second.values)
        std::cout << v << "\t";
    std::cout << std::endl;
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


void setDefaultParameters(ParameterMap& parameterMap)
{
  addParameter(parameterMap, "NUM_AGENTS", "Number of agents", {10000.0});
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
      std::cerr << "Error at line: " << lineNumber << std::endl;
      std::cerr << errorMessage << std::endl;
    }
  }
}

/****************************/
/* Simulation data */

class Agent {
public:
  unsigned id;
  Agent* partner;
  std::vector<Agent*> past_partners;
  double infectivity;
  double infectiousness;
  double breakupiness;
  bool infected;
  unsigned rel;
  double age;
  unsigned sex;
  unsigned sexor;
  unsigned page;
  unsigned psex;
  unsigned psexor;
};

typedef std::vector<Agent *> AgentVector;

/* Initialization routines */


/* Supporting functions for Stefan's initial population creation algorithm. */



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
    sample_matWW[i].init(getCol(matWW,i), &rng);
    sample_matWM[i].init(getCol(matWM,i), &rng);
    sample_matMW[i].init(getCol(matMW,i), &rng);
    sample_matMM[i].init(getCol(matMM,i), &rng);
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

void initializeAgents(AgentVector& agents,
                      const ParameterMap& parameterMap)
{
  CSVParser data_csv(parameterMap.at("AGENT_DATA_CSV").stringValue.c_str(),
                     ";", true);
  DblMatrix data = data_csv.convert_all_entries_to_doubles();
  CSVParser singles_csv(parameterMap.at("SINGLES_DATA_CSV").stringValue.c_str(),
                        ",", true);
  DblMatrix singles = singles_csv.convert_all_entries_to_doubles();
  CSVParser partners_csv(parameterMap.at("PARTNERS_DATA_CSV").stringValue.c_str(),
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

  unsigned X = parameterMap.at("NUM_AGENTS").values[0];
  unsigned S = calcNumberSingles(data, X);


  std::vector<double> ageRange;
  for (unsigned i = 12; i <= 100; ++i)
    ageRange.push_back(i);
  auto ageShare = getCol(singles, 1);
  auto femRatio = getCol(singles, 2);
  auto msmRate = getCol(singles, 3);
  auto wswRate = getCol(singles, 4);
  create_singles(agents, data, S, ageRange, ageShare,
		 femRatio, wswRate, msmRate);
  ageShare = femRatio = msmRate = wswRate = {};
  ageShare = getCol(partners, 1);
  femRatio = getCol(partners, 2);
  msmRate = getCol(partners, 3);
  wswRate = getCol(partners, 4);
  create_partners_ACP(agents, data, S, X, ageRange, ageShare, femRatio,
  		      wswRate, msmRate, ww, mw, wm, mm);
}

/* Simulation routines */

void simulate(const ParameterMap& parameterMap)
{
  AgentVector agents;
  initializeAgents(agents, parameterMap);

  double start_date = parameterMap.at("START_DATE").getDbl();
  double end_date = parameterMap.at("END_DATE").getDbl();
  double time_step = parameterMap.at("TIME_STEP").getDbl();
  for (double current_date = start_date; current_date < end_date;
       current_date += time_step) {

  }
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



int main(int argc, char *argv[])
{
  ParameterMap parameterMap;

  setDefaultParameters(parameterMap);
  if (cmdOptionExists(argv, argv + argc, "-h")) {
    std::cout << argv[0] << " options, where options are:\n"
              << "-h: help - Print this message.\n"
              << "-f filename: Use filename as input.\n\n"
              << "Parameters\n"
              << "**********\n";
    printParameters(parameterMap);
    exit(1);
  }

  const char *input_file_str = getCmdOption(argv, argv + argc, "-f");
  if (input_file_str) {
    std::ifstream infile;
    infile.open (input_file_str, std::ifstream::in);
    if (infile.fail()) {
      std::cerr << "Error opening " << input_file_str << std::endl;
      exit(1);
    }
    readParameters(infile, parameterMap);
    infile.close();
  } else {
    readParameters(std::cin, parameterMap);
  }
  simulate(parameterMap);
  printParameters(parameterMap);
}
