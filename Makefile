#########################
# customise these
CXXFILES := main.cpp readArgs.cpp
PROG := KmerCounter
CXXFLAGS := -Wall -Wextra -O3 -std=c++11
LDFLAGS :=
########################

# -MMD generates dependencies while compiling
CXXFLAGS += -MMD
CXX := g++
CC := g++

OBJFILES := $(CXXFILES:.cpp=.o)
DEPFILES := $(CXXFILES:.cpp=.d)

$(PROG) : $(OBJFILES)
	$(LINK.o) $(LDFLAGS) -o $@ $^

clean :
	rm -f $(PROG) $(OBJFILES) $(DEPFILES)

-include $(DEPFILES)