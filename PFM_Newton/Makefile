# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# the build target executable:
TARGET = main
HEADERS +=

SOURCES += main.cpp 

OBJECTS = $(SOURCES:.cpp=.o)

# define the executable file 
MAIN = main


.PHONY: depend clean

all:    $(MAIN)
	@echo  "\nuloha ok\n"

$(MAIN): $(OBJECTS) 
	$(CXX) $(CPPFLAGS) $(HEADERS) -o $(MAIN) $(OBJECTS) $(LD_FLAGS)
	$(RM) *.o

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXX_FLAGS) -c -o $@ $<

clean:
	$(RM) *.o *~ $(MAIN) *.txt


depend: $(SOURCES)
	makedepend $(HEADERS) $^

# DO NOT DELETE THIS LINE -- make depend needs it

