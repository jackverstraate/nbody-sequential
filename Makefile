CXX := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -pedantic

all: nbody

nbody: nbody.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f nbody