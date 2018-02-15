CXX = g++

debug_or_optimize = -O1
CXXFLAGS = -Wall -Werror -pedantic --std=c++11 $(debug_or_optimize)

all: test simulation

debug: debug_or_optimize = -g
debug: test simulation

test: ./test_particle

output_pos.txt: 
	touch output_pos.txt

test_particle: particle.cpp test_particle.cpp unit_test_framework.cpp array_helpers.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

simulation: LWS_driver.cpp simulation.cpp particle.cpp array_helpers.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@


.SUFFIXES:

.PHONY: clean
clean:
	rm -vf *.out *.exe simulation
