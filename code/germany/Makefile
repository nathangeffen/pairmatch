CXX = g++

CXXFLAGS = -Wall -std=c++11
DEVFLAGS  = -g
RELFLAGS = -O3
LDFLAGS =

# the build target executable:
SOURCES = germany_partners.cc
ADD_SOURCES = csvparser.c ransampl.c CSV_Parser.cc
OBJECTS = germany_partners.o csvparser.o ransampl.o CSV_Parser.o
EXECUTABLE = germany_partners

all: $(SOURCES) $(EXECUTABLE)-dev

$(EXECUTABLE)-dev: $(OBJECTS)
	$(CXX) $(DEVFLAGS) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)-dev

germany_partners.o: germany_partners.cc
	$(CXX) -c $(DEVFLAGS) $(CXXFLAGS) germany_partners.cc -o germany_partners.o

CSV_Parser.o: CSV_Parser.cc
	$(CXX) -c $(DEVFLAGS) $(CXXFLAGS) CSV_Parser.cc -o CSV_Parser.o

csvparser.o: csvparser.c
	$(CXX) -c $(DEVFLAGS) $(CXXFLAGS) csvparser.c -o csvparser.o

ransampl.o: ransampl.c
	$(CXX) -c $(DEVFLAGS) $(CXXFLAGS) ransampl.c -o ransampl.o


release: clean
	$(CXX) $(RELFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $(EXECUTABLE)-rel $(SOURCES) $(ADD_SOURCES)

clean:
	rm -f $(EXECUTABLE)-rel $(EXECUTABLE)-dev *.o
