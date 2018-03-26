#include "particle.h"
#include "array_helpers.h"

#include <cmath>
#include <iomanip>
#include <utility>
#include <stdexcept>

using namespace std;

Particle::Particle(const double m_in, const double q_in, 
		const array <double, 3>& startPos,
		const array <double, 3>& startVel,  const int total_steps_in, const double dt_in)
	 				: m(9.11e-31 * m_in), q(1.6e-19 * q_in), total_steps(total_steps_in),
	 				current_steps(1), dt(dt_in)
{
	for (int i = 0; i < 3; ++i)
	{
		pos[i] = new double[total_steps+1];
		vel[i] = new double[total_steps+1];

		pos[i][0] = startPos[i];
		vel[i][0] = startVel[i];
	}

	all_times = new double[total_steps+1];
	for (int i = 0; i < total_steps+1; ++i)
		all_times[i] = i * dt;

}

Particle::Particle(const Particle& other)
{
	for (int i = 0; i < 3; ++i)
	{
		pos[i] = new double[other.total_steps+1];
		vel[i] = new double[other.total_steps+1];
	}

	all_times = new double[other.total_steps+1];

	copy_particle(other);
}

Particle& Particle::operator=(const Particle& rhs)
{
	if (this == &rhs)
		return *this;
	copy_particle(rhs);
	return *this;
}

Particle::~Particle()
{
	for (int i = 0; i < 3; ++i)
	{
		delete[] pos[i];
		delete[] vel[i];
	}

	delete[] all_times;
}


array<double, 3> Particle::getPos_step(const int timestep) const
{

	if (timestep >= current_steps)
		throw invalid_argument("asked position at timestep > current_step");

	if (timestep < 0)
		return { pos[0][0], pos[1][0], pos[2][0] };

	return { pos[0][timestep], pos[1][timestep], pos[2][timestep] };
}


array<double, 3> Particle::getVel_step(const int timestep) const
{
	if (timestep >= current_steps)
		throw invalid_argument("asked velocity at timestep > current_step");

	if (timestep < 0)
		return { pos[0][0], pos[1][0], pos[2][0] };

	return { vel[0][timestep], vel[1][timestep], vel[2][timestep] };
}




// REQUIRES: time is smaller than time at current timestep
// EFFECTS: outputs position at given time by interpolating between two
// nearest timsteps (if exact timestep is not present)
array<double, 3> Particle::getPos(const double time) const
{
	if (time > all_times[current_steps-1])
		throw invalid_argument("asked position at time t > simulation time");

	if (time <= 0 || current_steps <= 1) {
		return getPos_step(0);
	} else {

		pair<int, bool> time_interp = 
				find_time(0, current_steps-1, time);

		int index_time = time_interp.first;
		bool unique_time = time_interp.second;

		

		if (unique_time == true) {
			return getPos_step(index_time);
			
		} else {
			double time_prev = all_times[index_time];
			double time_next = all_times[index_time+1];
	
			double frac_prev = (time_next - time) / dt;
			double frac_next = (time - time_prev) / dt;
	
			return array_add( array_mult_scalar( getPos_step(index_time), frac_prev), 
										  array_mult_scalar( getPos_step(index_time+1), frac_next) );
		}
	}
}

// REQUIRES: time is smaller than time at current timestep
// EFFECTS: outputs velocity at given time by interpolating between two
// nearest timsteps (if exact timestep is not present)
array<double, 3> Particle::getVel(const double time) const
{
	if (time > all_times[current_steps-1])
		throw invalid_argument("asked velocity at time t > simulation time");

	if (time <= 0 || current_steps <= 1) {
		return getVel_step(0);
	} else {

		// take into account that velocities are stored at (n-1/2)*dt
		pair<int, bool> time_interp = 
				find_time(0, current_steps-1, time+dt/2);

		int index_time = time_interp.first;
		bool unique_time = time_interp.second;

		

		if (unique_time == true) {
			return getVel_step(index_time);

		} else {
			double time_prev = all_times[index_time];
			double time_next = all_times[index_time+1];
	
			double frac_prev = (time_next - time) / dt;
			double frac_next = (time - time_prev) / dt;
	
			return array_add( array_mult_scalar( getVel_step(index_time), frac_prev), 
										  array_mult_scalar( getVel_step(index_time+1), frac_next) );
		}
	}
}
// REQUIRES: time is smaller than time at current timestep
// EFFECTS: outputs acceleration at given time
array<double, 3> Particle::getAccel(const double time) const
{
	if (time > all_times[current_steps-1])
		throw invalid_argument("asked acceleration at time t > simulation time");

	if (time <= 0 || current_steps <= 1) {
		return {0., 0. , 0.};
	} else {

		// take into account that velocities are stored at (n-1/2)*dt
		pair<int, bool> time_interp = 
				find_time(0, current_steps-1, time+dt/2);

		int index_time = time_interp.first;

		if (index_time == current_steps-1) 
		{
			return array_mult_scalar(array_substract(getVel_step(index_time), 
																							getVel_step(index_time-1)),
															1. / dt);
			
		} else {
			return array_mult_scalar(array_substract(getVel_step(index_time+1), 
																							getVel_step(index_time)),
															1. / dt);
		}
	}
}

double Particle::find_retarded_time(const Particle& other) const
{
	if (&other == this)
		throw invalid_argument("tried to compute self-interaction.");
	
	// use newton's algorithm to find root t_0
	array<double, 3> x_curr = getPos_step(current_steps-1);
	double 					 t_curr = all_times[current_steps-1];

	// first guess for t_0
	double t_0 = t_curr - array_norm(array_substract(x_curr, 
											other.getPos_step(current_steps-1))) * INV_C_LIGHT;

	for (int i = 0; i < NR_NEWTON_STEPS; i++)
	{
		array<double, 3> x_diff = array_substract(x_curr, other.getPos(t_0));
		double norm_x_diff = array_norm(x_diff);

		t_0 -= (t_0 - t_curr + norm_x_diff * INV_C_LIGHT) / 
			(1. + array_dot_product(x_diff, getVel(t_0)) / (C_LIGHT * norm_x_diff)  );
	}

	return t_0;
}


void Particle::compute_E_B(const Particle& other, double t_0, 
			std::array<double, 3>& E_field, std::array<double, 3>& B_field) const
{
	if (&other == this)
		throw invalid_argument("tried to compute self-interaction.");

	// calculate E and B via Lienard Wiechert potentials
	array<double, 3> r_0 = other.getPos(t_0);
	array<double, 3> beta_0 = array_mult_scalar(other.getVel(t_0), INV_C_LIGHT);
	array<double, 3> a_0 = other.getAccel(t_0);


	array<double, 3> diff_x_r = array_substract(getPos_step(current_steps-1), r_0);
	double inv_norm_diff_x_r = 1. / array_norm(diff_x_r);

	array<double, 3> n_0 = array_mult_scalar(diff_x_r, inv_norm_diff_x_r);
	double inv_K = 1. / (1. - array_dot_product(n_0, beta_0));
	array<double, 3> n_0_min_beta_0 = array_substract(n_0, beta_0);


	array<double, 3> E_field_first_term = array_mult_scalar( n_0_min_beta_0,
						(1. - array_norm_sq(beta_0)) * pow(inv_norm_diff_x_r, 2) );

	array<double, 3> E_field_second_term = array_mult_scalar(
		 	array_cross_product(n_0, array_cross_product(n_0_min_beta_0, a_0)),
		 	pow(INV_C_LIGHT, 2) * inv_norm_diff_x_r );

	array<double, 3> E_change = array_mult_scalar( 
			array_add( E_field_first_term, E_field_second_term) , 
			q * pow(inv_K, 3) / ( 4 * PI * EPS_0) );

	array<double, 3> B_change = array_mult_scalar( 
			array_cross_product(n_0, E_change), INV_C_LIGHT);


	E_field = array_add(E_field, E_change);
	B_field = array_add(B_field, B_change);
}


void Particle::advance(const array<double, 3>& E_field, 
											 const array<double, 3>& B_field)
{
	if (current_steps-1 == total_steps)
		throw invalid_argument("total number of steps exceeded");

	array<double, 3> vel_curr = getVel_step(current_steps-1);
  double gamma_curr = 1. / sqrt(1 - pow(INV_C_LIGHT, 2) * array_norm_sq(vel_curr));
  array<double, 3> mom_curr_non_norm = array_mult_scalar(vel_curr, gamma_curr*m);

	// relativistic Boris particle pusher

  array<double, 3> delta_mom_half_step_non_norm = array_mult_scalar(E_field, dt * q / 2.);

  array<double, 3> mom_half_step = array_mult_scalar( array_add(mom_curr_non_norm, delta_mom_half_step_non_norm),
       1./(m * C_LIGHT));

  double root = q*dt / (2. * m * sqrt(1. + array_norm_sq(mom_half_step)) );

  array<double, 3> tau = array_mult_scalar(B_field, root);

  double tau_scalar = 1. / (1 + array_norm_sq(tau));

  array<double, 3> mom_prime_1 = array_mult_scalar(mom_half_step, 1. - array_norm_sq(tau));
  array<double, 3> mom_prime_2 = array_mult_scalar(tau, array_dot_product(tau, mom_half_step) * 2);
  array<double, 3> mom_prime_3 = array_mult_scalar(array_cross_product(tau, mom_half_step), -2.);
  array<double, 3> mom_prime = array_mult_scalar(
    array_add(mom_prime_1, array_add(mom_prime_2, mom_prime_3)), tau_scalar);

  array<double, 3> new_mom = array_add(mom_prime, 
    array_mult_scalar(delta_mom_half_step_non_norm, 1./ (m*C_LIGHT)));

  double inv_gamma = 1. / sqrt(1. + array_norm_sq(new_mom));

  array<double, 3> new_vel = array_mult_scalar(new_mom, C_LIGHT * inv_gamma);

  vel[0][current_steps] = new_vel[0];
  vel[1][current_steps] = new_vel[1];
  vel[2][current_steps] = new_vel[2];

  pos[0][current_steps] = pos[0][current_steps-1] + new_vel[0] * dt;
  pos[1][current_steps] = pos[1][current_steps-1] + new_vel[1] * dt;
  pos[2][current_steps] = pos[2][current_steps-1] + new_vel[2] * dt;



	// double gamma = 1./ sqrt(1. - array_norm_sq(vel_curr) * pow(INV_C_LIGHT, 2) );

	// array<double, 3> p_min = array_add( array_mult_scalar(vel_curr, gamma * m),
	// 			array_mult_scalar(E_field, q * dt / 2.) );

	// array<double, 3> t_vec = array_mult_scalar(B_field, q * dt / (2. * gamma * m));

	// array<double, 3> s_vec = array_mult_scalar(t_vec, 2. / (1 + array_norm_sq(t_vec)) );

	// array<double, 3> p_bar = array_add(p_min, array_cross_product(p_min, t_vec));

	// array<double, 3> p_plus = array_add(p_min, array_cross_product(p_bar, s_vec));

	// array<double, 3> p_new = array_add(p_plus, array_mult_scalar(E_field, q * dt / 2.));


	// double inv_gamma_m = 1. / ( m * sqrt(1 + array_norm_sq(p_new) / pow(m * C_LIGHT, 2)) );

	// vel[0][current_steps] = p_new[0] * inv_gamma_m;
	// vel[1][current_steps] = p_new[1] * inv_gamma_m;
	// vel[2][current_steps] = p_new[2] * inv_gamma_m;

	// pos[0][current_steps] = pos[0][current_steps-1] + vel[0][current_steps] * dt;
	// pos[1][current_steps] = pos[1][current_steps-1] + vel[1][current_steps] * dt;
	// pos[2][current_steps] = pos[2][current_steps-1] + vel[2][current_steps] * dt;


	current_steps++;
}

std::ostream& Particle::print_pos(std::ostream& os, int timestep) const
{
	os << scientific;
	os.precision(PRECISION);

	os << "\t" << setw(OUTPUT_WIDTH) << pos[0][timestep] << "\t" 
						 << setw(OUTPUT_WIDTH) << pos[1][timestep] << "\t" 
						 << setw(OUTPUT_WIDTH) << pos[2][timestep];

	return os;
}



// 
// HELPER FUNCTIONS
// 

void Particle::copy_particle(const Particle& other)
{

	m = other.m;
	q = other.q;
	total_steps = other.total_steps;
	current_steps = other.current_steps;

	for (int j = 0; j < current_steps; ++j)
	{
		for (int i = 0; i < 3; ++i)
		{
			pos[i][j] = other.pos[i][j];
			vel[i][j] = other.vel[i][j];
		}
	}

	for (int i = 0; i < total_steps; ++i)
		all_times[i] = other.all_times[i];
}

// EFFECTS: returns the index of the time in all_times which is equal to time (true)
// or closest to time and smaller (false)
pair<int, bool> Particle::find_time(int lower_index, int upper_index, double time) const
{
	double lower_time = all_times[lower_index];
	double upper_time = all_times[upper_index];

	if (lower_time == time || lower_index == upper_index)
		return pair<int, bool>(lower_index, true);
	else if (upper_time == time)
		return pair<int, bool>(upper_index, true);
	else if (upper_index - lower_index == 1)
		return pair<int, bool>(lower_index, false);
	else {
		int middle_index = (lower_index + upper_index) / 2;

		if (time < all_times[middle_index])
			return find_time(lower_index, middle_index, time);
		else
			return find_time(middle_index, upper_index, time);
	}
}