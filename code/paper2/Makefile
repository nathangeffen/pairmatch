# the compiler: gcc for C program, define as g++ for C++
CXX = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
DEBUGFLAGS = -DDEBUG -g
CXXFLAGS  = -Wall -std=c++11 -pthread
RELFLAGS = -O3
EXTRAFLAGS =

OBJDIR = bin/
SRCDIR = src/

INCLUDES_CC = $(addprefix $(SRCDIR), simulate.hh agent.hh CSV_Parser.hh parameters.hh sample.hh linear.hh common.hh symboltable.hh)
INCLUDES_C = csvparser.h ransampl.h
INCLUDES = INCLUDES_CC INCLUDES_C
SOURCES_CC = main.cc simulate.cc agent.cc parameters.cc CSV_Parser.cc linear.cc
SOURCES_C = csvparser.c ransampl.c
SOURCES = $(addprefix $(SRCDIR), $(SOURCES_CC)) $(addprefix $(SRCDIR), $(SOURCES_C))
OBJS = $(addprefix $(OBJDIR), $(SOURCES_CC:.cc=.o)) $(addprefix $(OBJDIR), $(SOURCES_C:.c=.o))
TARGET_NAME = simulate
TARGET = $(OBJDIR)$(TARGET_NAME)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) $(EXTRAFLAGS) -o $(TARGET) $(OBJS)

$(OBJDIR)main.o: $(SRCDIR)main.cc $(SRCDIR)simulate.hh
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)main.cc -o $@

$(OBJDIR)simulate.o: $(SRCDIR)simulate.cc $(INCLUDES_CC)
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)simulate.cc -o $@

$(OBJDIR)agent.o: $(SRCDIR)agent.cc $(SRCDIR)common.hh
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)agent.cc -o $@

$(OBJDIR)parameters.o: $(SRCDIR)parameters.cc $(SRCDIR)parameters.hh $(SRCDIR)common.hh
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)parameters.cc -o $@

$(OBJDIR)linear.o: $(SRCDIR)linear.cc $(SRCDIR)linear.hh
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)linear.cc -o $@

$(OBJDIR)CSV_Parser.o: $(SRCDIR)CSV_Parser.cc $(SRCDIR)CSV_Parser.hh
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)CSV_Parser.cc -o $@

$(OBJDIR)csvparser.o: $(SRCDIR)csvparser.c $(SRCDIR)csvparser.h
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)csvparser.c -o $@

$(OBJDIR)ransampl.o: $(SRCDIR)ransampl.c $(SRCDIR)ransampl.h
	$(CXX) $(DEBUGFLAGS) $(CXXFLAGS) -c $(SRCDIR)ransampl.c -o $@

release: clean
	$(CXX) $(RELFLAGS) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

clean:
	$(RM) $(OBJS) $(TARGET)
