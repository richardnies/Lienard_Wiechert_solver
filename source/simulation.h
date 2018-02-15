#ifndef GUARD_simulation
#define GUARD_simulation

#include "particle.h"

#include <iostream>
#include <vector>

class Simulation {
public:
	Simulation(const char* filename);

	~Simulation();

	void run();

	std::ostream& print_pos(std::ostream& os);
	
private:
	std::vector<Particle*> particles;

	double dt_sim;
	int total_steps_sim;
	int current_steps_sim;
};

#endif