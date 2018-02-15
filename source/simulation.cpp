#include "simulation.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace std;

Simulation::Simulation(const char* filename): current_steps_sim(1)
{
	ifstream inFile;
	try {
		inFile.open(filename);
	} catch(const ifstream::failure& e) {
		throw runtime_error("Could not open file.");
	}

	inFile >> dt_sim >> total_steps_sim;

	double m, q;
	array<double, 3> pos, vel;

	while(inFile >> m >> q >> pos[0] >> pos[1] >> pos[2] 
		>> vel[0] >> vel[1] >> vel[2])
	{
		particles.push_back(new Particle(m, q, pos, vel, total_steps_sim, dt_sim));
	}

	cout << "Initialised " << particles.size() << " particles." << endl;

	inFile.close();
}


Simulation::~Simulation()
{
	while (!particles.empty())
	{
		delete particles.back();
		particles.pop_back();
	}
}


void Simulation::run()
{
	for (int step_nr = 0; step_nr < total_steps_sim; step_nr++)
	{
		for (int part_nr = 0; part_nr < (int)particles.size(); ++part_nr)
		{
			// find E and B fields from every other particle
			array<double, 3> E_field = {0, 0, 0};
			array<double, 3> B_field = {0, 0, 0};
			double t_0;

			for (int part_nr_inter = 0; part_nr_inter < (int)particles.size();
					 ++part_nr_inter)
			{
				// no self-interaction
				if (part_nr_inter == part_nr)
					continue;

				t_0 = particles[part_nr]->find_retarded_time(*particles[part_nr_inter]);

				particles[part_nr]->compute_E_B(*particles[part_nr_inter], t_0, 
							E_field, B_field);
			}

			particles[part_nr]->advance(E_field, B_field);
		}

		cout << "Finished step " << step_nr+1 << " of " << total_steps_sim << endl;
	}
}



std::ostream& Simulation::print_pos(std::ostream& os)
{
	os << setw(OUTPUT_WIDTH) << "TIME";
	for (int part_nr = 0; part_nr < (int)particles.size(); ++part_nr)
	{
		os << setw(OUTPUT_WIDTH) << "\tx_" << part_nr << 
					setw(OUTPUT_WIDTH) << "\ty_" << part_nr << 
					setw(OUTPUT_WIDTH) << "\tz_" << part_nr;
	}
	os << "\n";


	for (int step_nr = 0; step_nr < total_steps_sim; step_nr++)
	{
		os << scientific;
		os.precision(PRECISION);

		os << setw(OUTPUT_WIDTH) << dt_sim * step_nr;

		for (int part_nr = 0; part_nr < (int)particles.size(); ++part_nr)
		{
			particles[part_nr]->print_pos(os, step_nr);
		}

		os << "\n";
	}

	return os;
}