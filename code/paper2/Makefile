# the compiler: gcc for C program, define as g++ for C++
CXX = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CXXFLAGS  = -g -Wall -std=c++11
RELFLAGS = -O3 -std=c++11
EXTRAFLAGS =

# the build target executable:
TARGET = simulate

all: $(TARGET)

$(TARGET): $(TARGET).cc CSV_Parser.o CSV_Parser.hh
	$(CXX) $(CXXFLAGS) $(EXTRAFLAGS) -o $(TARGET) $(TARGET).cc

CSV_Parser.o: CSV_Parser.cc CSV_Parser.hh
	$(CXX) $(CXXFLAGS) -c CSV_Parser.cc csvparser.c

release: clean
	$(CXX) $(RELFLAGS) $(EXTRAFLAGS) -o $(TARGET) $(TARGET).cc CSV_Parser.cc csvparser.c

clean:
	$(RM) $(TARGET)
