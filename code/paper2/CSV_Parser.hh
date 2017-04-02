#ifndef CSVPARSER_HH
#define CSVPARSER_HH

#include <iostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

extern "C" {
#include "csvparser.h"
}

typedef std::vector< std::vector<double > > DblMatrix;

class CSVParser {
public:
  CSVParser(const char* filename,
	    const char* delim = ",",
	    const bool has_header = true) {
    csvparser_ = CsvParser_new(filename, delim, has_header);
    if (has_header) {
      header_ = CsvParser_getHeader(csvparser_);
      if (!header_) {
	std::cerr << "Error getting csv header in " << filename << std::endl;
	exit(1);
      }
    }
    CsvRow *row;
    while ((row = CsvParser_getRow(csvparser_)) ) {
      rows_.push_back(row);
      std::vector< std::string > str_row;
      const char **rowFields = CsvParser_getFields(row);
      for (int i = 0 ; i < CsvParser_getNumFields(row) ; i++)
	str_row.push_back(rowFields[i]);
      string_rows.push_back(str_row);
    }
  };
  DblMatrix convert_all_entries_to_doubles() {
    DblMatrix double_rows;
    for (auto& r: string_rows) {
      std::vector<double> double_row;
      for (auto& s: r) {
	boost::replace_all(s, ",", ".");
	boost::replace_all(s, "NA", "0");
	double_row.push_back(stod(s));
      }
      double_rows.push_back(double_row);
    }
    return double_rows;
  };
  ~CSVParser() {
    for (auto& row: rows_)
      CsvParser_destroy_row(row);
    CsvParser_destroy(csvparser_);
  };
  std::vector< std::vector<std::string> > string_rows;
private:
  CsvRow *header_ = NULL;
  std::vector<CsvRow *> rows_;

  CsvParser *csvparser_;
};

void print_dbl_matrix(const DblMatrix&);

#endif
