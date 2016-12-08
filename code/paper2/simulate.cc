#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "CSV_Parser.hh"

#define BFPM 1
#define RPM 2
#define CSPM 3

struct ParameterValue {
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
  addParameter(parameterMap, "MSM_CSV",
               "MSM age preferences CSV file", "data/dist_msm.csv");
  addParameter(parameterMap, "MSW_CSV",
               "MSW age preferences CSV file", "data/dist_msw.csv");
  addParameter(parameterMap, "WSW_CSV",
               "WSW age preferences CSV file", "data/dist_wsw.csv");
  addParameter(parameterMap, "WSM_CSV",
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
  double age;
  int sex;
  unsigned preferred_age;
  int preferred_sex;
  bool wants_relationship;
};

/* Simulation routines */

void simulate(ParameterMap& parameterMap)
{

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
