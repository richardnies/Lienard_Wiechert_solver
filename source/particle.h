#ifndef GUARD_Particle
#define GUARD_Particle

#include <array>
#include<iostream>

const double C_LIGHT = 299792458.;
const double INV_C_LIGHT = 1. / C_LIGHT;
const double EPS_0 = 8.854187817e-12;
const double PI = 3.14159265359;

const int OUTPUT_WIDTH = 12;
const int PRECISION = 5;

const int NR_NEWTON_STEPS = 3;

class Particle
{
private:
	void copy_particle(const Particle& other);
	std::pair<int, bool> find_time
		(int lower_index, int upper_index, const double time) const;

public:

	Particle(const double m_in, const double q_in, 
			const std::array<double, 3>& startPos, const std::array<double, 3>& startVel, 
			 const int total_steps_in, const double dt_in);

	Particle(const Particle& otherParticle);

	Particle& operator=(const Particle& rhs);

	~Particle();

	double find_retarded_time(const Particle& other) const;

	void compute_E_B(const Particle& other, double t_0_prime, 
				std::array<double, 3>& E_field, std::array<double, 3>& B_field) const;

	void advance(const std::array<double, 3>& E_field, 
							const std::array<double, 3>& B_field);


	std::array<double, 3> getPos_step(const int timestep) const;
	std::array<double, 3> getVel_step(const int timestep) const;


	std::array<double, 3> getPos(const double time) const;
	std::array<double, 3> getVel(const double time) const;
	std::array<double, 3> getAccel(const double time) const;

	std::ostream& print_pos(std::ostream& os, int timestep) const;

private:
	// particle position at timesteps n*dt (starting at 0)
	std::array<double*, 3> pos;
	// particle velocities at timesteps (n-1/2)*dt (starting at 0)
	std::array<double*, 3> vel;
	// Particle mass in units of electron mass
	double m;
	// Particle charge in units of proton charge
	double q;
	// all times for which particle and velocity of the particle are available
	double* all_times;
	// total number of steps of the simulation
	int total_steps;
	// current number of steps in the simulation (1 step for t = 0!)
	int current_steps;
	// timestep size
	double dt;
};

#endif