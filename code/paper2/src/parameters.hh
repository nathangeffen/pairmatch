#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <cassert>
#include <exception>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <random>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "common.hh"

#define NO_PARMS 0
#define ALL_PARMS 1
#define VARYING_PARMS 2

enum RangeType {
  NONE = 0,
  THREE_PARM,
  COMMA_SEPARATED
};

class BadParameter: public std::exception
{
public:
  explicit BadParameter(const char* parameter, const char* additionalMsg = "")
  {
    msg.append(additionalMsg);
    msg.append(" for ");
    msg.append(parameter);
  }

  virtual const char* what() const throw (){
    return msg.c_str();
  }
protected:
  std::string msg;
};

class UnknownParameter: public BadParameter
{
public:
  explicit UnknownParameter(const char* parameter, const char* additionalMsg = "")
    : BadParameter(parameter, additionalMsg)
  {}
protected:
  std::string msg = "Unknown parameter ";
};


class ConfigurationError: public std::exception
{
public:
  explicit ConfigurationError(const char* message, unsigned lineNo)
  {
    msg.append(message);
    msg.append(" at line ");
    msg.append(std::to_string(lineNo));
  }

  virtual ~ConfigurationError() throw (){}

  virtual const char* what() const throw (){
    return msg.c_str();
  }

protected:
  std::string msg;
};


std::vector<double> setRange(const double, const double to, const double);
std::vector<double> setRangeUniRand(const size_t, const double, const double);


/**
   Class to handle a parameter, which is a pair consisting of a string key
   associated with either a vector of doubles or a string.
*/

struct ParameterValue {
  inline double dbl(size_t index = 0) const
  {
    return values[index];
  };
  inline std::string str() const
  {
    return strValue;
  };
  inline const char* c_str() const
  {
    return str().c_str();
  }
  inline const bool isSet() const
  {
    return (values[0] != 0.0) ? true : false;
  }

  const char* description;
  std::vector<double> values;
  bool isString = false;
  std::string strValue;
  bool isVaryingRange = false;
  std::string parentParameter;
};

/**
   Keeps track of all parameters and indicates which ones are varying.
*/

class ParameterMap  : public std::unordered_map<std::string,  ParameterValue>  {
public:

  /**
     Creates all parameters and sets them to their default parameters. Users of
     this code should always be able to consult the list of parameters by
     examining this method.

     @param parameterMap to create
  */

  ParameterMap(bool setDefaults = true) :
    std::unordered_map<std::string,  ParameterValue>()
  {
    if (setDefaults) setDefaultParameters();
  }

  void setDefaultParameters()
  {
    addParameter("SIMULATION_NAME",
                 "Name of simulation used in output", "Default");
    addParameter("NUM_SIMULATIONS",
                 "Number of simulations to execute (default is 1)", {1});
    addParameter("NUM_THREADS",
                 "Number of threads (default is one per simulation - 0)", {0});
    addParameter("NUM_AGENTS", "Number of agents", {100.0});

    addParameter("START_DATE",
                 "Start date of simulation", {2017.0});
    addParameter("END_DATE",
                 "End date of simulation", {2018.0});
    addParameter("TIME_STEP",
                 "Time step for each iteration of simulation", {DAY});
    addParameter("PREV_PARTNER_PENALTY",
                 "Previous partner penalty in pair matching", {10.0});

    addParameter("SHAPE_REL_CSV",
                 "File of shapes for partnership", "data/Rel_shape.csv");

    addParameter("SCALE_REL_CSV",
                 "File of scales for partnership", "data/Rel_scale.csv");

    addParameter("HET_MALE_INFECTIOUSNESS",
                 "Daily risk of infecting for male in heterosexual sex",
                 {0.005});
    addParameter("HET_FEMALE_INFECTIOUSNESS",
                 "Daily risk of infecting for female in heterosexual sex",
                 {0.001});

    addParameter("HOM_MALE_INFECTIOUSNESS",
                 "Daily risk of infecting for male in homosexual sex",
                 {0.005});
    addParameter("HOM_FEMALE_INFECTIOUSNESS",
                 "Daily risk of infecting for female in homosexual sex",
                 {0.0001});

    addParameter("MATCH_NEIGHBORS",
                 "Value of k (neighbors for matching", {300});
    addParameter("MATCH_CLUSTERS",
                 "Number of clusters for cluster shuffle matching", {100});
    addParameter("MATCH_SCORE_FAIL",
                 "Highest score beyond which a match between two agents fails",
                 {2000.00});
    addParameter("MATCH_SCORE_POOR",
                 "Score for which a poor match must be registered",
                 {49.50});

    addParameter("OUTPUT_AGENTS_AFTER_INIT",
                 "Print agent info after initialization"
                 "(0=no,1=yes)", {0.0});
    addParameter("OUTPUT_AGENTS_DURING_SIM",
                 "Print agent info during simulation"
                 "(0=no,1=yes), (frequency of output)", {0.0, 100.0});
    addParameter("OUTPUT_AGENTS_AT_END",
                 "Print agent info at end of simulation "
                 "(0=no,1=yes)", {0.0});

    addParameter("ANALYSIS_AGE_INTERVAL",
                 "Length of age intervals for infection analysis", {5.0});
    addParameter("ANALYZE_AFTER_INIT",
                 "Calculate stats after initialization"
                 "(0=no,1=yes)", {1.0});
    addParameter("ANALYZE_DURING_SIM",
                 "Calculate stats during simulation"
                 "(frequency of output, 0 for never)", {0});
    addParameter("ANALYZE_AT_END",
                 "Calculate stats at end of simulation "
                 "(0=no,1=yes)", {1.0});
    addParameter("PRINT_PARAMETERS",
                 "Print parameters using in simulation"
                 "(0=no,1=all,2=Varying only)", {VARYING_PARMS});


    addParameter("OUTPUT_INIT_TIMING",
                 "Output time every x iterations", {0});
    addParameter("OUTPUT_TIMING_DURING",
                 "Output time every x iterations", {0});
    addParameter("OUTPUT_TIMING_AFTER",
                 "Output time every x iterations", {1});
    addParameter("AGENT_DATA_CSV",
                 "Agent data for initialization", "data/data.csv");
    addParameter("SINGLES_DATA_CSV",
                 "Agent data for initialization", "data/singles.csv");
    addParameter("PARTNERS_DATA_CSV",
                 "Agent data for initialization", "data/partnersACP.csv");
    addParameter("MSM_DATA_CSV",
                 "Agent data for initialization", "data/mat.msm.csv");
    addParameter("MSW_DATA_CSV",
                 "Agent data for initialization", "data/mat.msw.csv");
    addParameter("WSW_DATA_CSV",
                 "Agent data for initialization", "data/mat.wsw.csv");
    addParameter("WSM_DATA_CSV",
                 "Agent data for initialization", "data/mat.wsm.csv");

    addParameter("MSM_AGE_DIST_CSV",
                 "MSM age preferences CSV file", "data/dist_msm.csv");
    addParameter("MSW_AGE_DIST_CSV",
                 "MSW age preferences CSV file", "data/dist_msw.csv");
    addParameter("WSW_AGE_DIST_CSV",
                 "WSW age preferences CSV file", "data/dist_wsw.csv");
    addParameter("WSM_AGE_DIST_CSV",
                 "WSM age preferences CSV file", "data/dist_wsm.csv");

    addParameter("INITIAL_INFECTION_RATES_CSV",
                 "Initial infection rates in population",
                 "data/initial_rates.csv");
    addParameter("SCALE_INITIAL_RATES",
                 "Scale the initial infection rates by this value", {1.0});

    addParameter("AGE_EVENT", "Execute the aging event", {1});
    addParameter("INFECT_EVENT", "Execute the infection event", {1});
    addParameter("BREAKUP_EVENT",
                 "Execute the infection event", {1});
    addParameter("OUTPUT_NUM_BREAKUPS",
                 "Print number of breakups on every time step", {1});
    addParameter("MATCH_EVENT",
                 "Execute the infection event (RPM,RKPM,CSPM)", "RPM");
    addParameter("OUTPUT_NUM_MATINGPOOL",
                 "Print number in mating pool on every time step", {1});
    addParameter("DISTANCE_METHOD",
                 "Distance calculation to use (HEURISTIC,TABLE)", "HEURISTIC");

    addParameter("PARTNERS_LOAD_FACTOR",
                 "Load factor for partnerships unordered set (advanced)", {0});
    addParameter("PARTNERS_RESERVE",
                 "Capacity for partnerships unordered set (advanced)", {0});

    addParameter("RANDOM_SEED", "Value to set random seed to",
                 {1});

    // Parameters that may need to be fitted
    addParameter("SHAPE_SINGLE_PERIOD_INITIAL",
                 "Weibull shape of single period at begin of simulation",
                 {0.85});
    addParameter("SCALE_SINGLE_PERIOD_INITIAL",
                 "Weibull scale of single period at begin of simulation",
                 {1200.0});
    addParameter("SHAPE_SINGLE_PERIOD_DURING",
                 "Weibull shape of single period during simulation", {0.4});
    addParameter("SCALE_SINGLE_PERIOD_DURING",
                 "Weibull scale of single period during simulation", {100.0});
    addParameter("MEAN_SINGLE_PERIOD",
                 "Mean difference from expected single period", {0});
    addParameter("SD_SINGLE_PERIOD",
                 "Standard deviation of single period", {3.0});

    addParameter("SCALE_RELATIONSHIP_PERIOD_INITIAL",
                 "Multiple relationship scale parameters by this", {1.0 / 3.0});
    addParameter("SCALE_RELATIONSHIP_PERIOD_DURING",
                 "Multiple relationship scale parameters by this", {1.0});
    addParameter("MEAN_RELATIONSHIP_PERIOD",
                 "Mean difference from expected relationship period.", {0});
    addParameter("SD_RELATIONSHIP_PERIOD",
                 "Standard deviation of relationship period.", {3.0});

    addParameter("PROB_INFECTED_IF_PARTNER",
                 "Probability on initialization of an agent being infected "
                 "if its partner is infected", {0.5});
  }

  /**
     Adds a parameter containing one or more doubles to a ParameterMap.

     @param parameterMap ParameterMap to be added to
     @param name key of the parameter
     @param description of the parameter
     @param values of the parameter
  */

  void addParameter(const char* name,
                    const char* description,
                    const std::initializer_list<double> values)
  {
    ParameterValue parameterValue;

    parameterValue.description = description;
    parameterValue.values = values;
    (*this)[name] = parameterValue;
  }

  /**
     Adds a parameter containing a string to a ParameterMap.

     @param parameterMap ParameterMap to be added to
     @param name key of the parameter
     @param description of the parameter
     @param strValue of the parameter
  */

  void addParameter(const char* name,
                    const char* description,
                    const char* strValue)
  {
    ParameterValue parameterValue;

    parameterValue.description = description;
    parameterValue.isString = true;
    parameterValue.strValue = strValue;
    (*this)[name] = parameterValue;
  }

  const ParameterValue& at ( const std::string& s ) const {
    try {
      return std::unordered_map<std::string,  ParameterValue>::at(s);
    } catch (const std::out_of_range& oor) {
      throw UnknownParameter(s.c_str());
    }
  }

  /**
     Prints the parameter keys and values in a ParameterMap.

     @param printDescription if true
     @param varyingParametersOnly if true prints only parameters that vary
     on each simulation
  */
  void print(const std::string& prefix,
             const bool printDescription = true,
             const bool varyingParametersOnly = false)
  {
    for (auto& p: (*this)) {
      std::ostringstream oss;
      if (varyingParametersOnly && p.second.isVaryingRange == false) continue;
      oss << prefix;
      oss << "PARAMETER," << p.first;
      if (printDescription) {
        oss << "\"" << p.second.description << "\"";
      }
      if (p.second.isString) {
        oss << "," << p.second.strValue;
      } else {
        for (auto& v: p.second.values) oss << "," << v;
      }
      printf("%s\n", oss.str().c_str());
    }
  }

  /**
     Replace the value of a parameter with values in a line from a stringstream.

     @param parameterMap ParameterMap to be added to
     @param name key of the parameter
     @param line containing the new values
  */

  void replaceParameter(const char* name,
                        std::istringstream& line,
                        const RangeType rangeType,
                        const bool varyWithPrevious)
  {
    auto it = this->find(name);
    if (it == this->end()) {
      throw UnknownParameter(name);
    }
    if (this->at(name).isString) {
      line >> (*this)[name].strValue;
    } else {
      double value;
      std::vector<double> values;

      while (line >> value)
        values.push_back(value);

      if ( (line.fail() && !line.eof()) || values.size() == 0) {
        throw BadParameter(name, "Invalid values");
      }
      if (rangeType == NONE) {
        (*this)[name].values = values;
      } else {
        (*this)[name].isVaryingRange = true;
        if (rangeType == THREE_PARM) {
          assert(values.size() == 3);
          (*this)[name].values = setRange(values[0], values[1], values[2]);
        } else if (rangeType == COMMA_SEPARATED) {
          (*this)[name].values = values;
        } else {
          throw std::runtime_error("Unknown range type");
        }
        if (varyWithPrevious) {
          (*this)[name].parentParameter =
            this->varyingParameters.back();
        } else {
          (*this)[name].parentParameter = "";
        }
        this->varyingParameters.push_back(name);
      }
    }
  }


  std::vector<std::string>  varyingParameters;
};


void readParameters(std::istream& input,
                    std::vector<ParameterMap>& parameterMaps);

#endif
