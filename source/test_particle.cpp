#include "particle.h"
#include "unit_test_framework.h"

#include <iostream>
#include <array>

using namespace std;




TEST(test_ctor_getPos_getVel) {
	double const eps = 1e-12;
	const int nr_steps = 1;

	const array<double, 3> startPos = {0, 0, 3.};
	const array<double, 3> startVel = {0, 1., 0};
	
	double m = 1;
	double q = 1;
	double dt_dummy = 1;

	Particle positron(m, q, startPos, startVel, nr_steps, dt_dummy);

	double neg_time = -1;

	array<double, 3> pos = positron.getPos(neg_time);
	array<double, 3> vel = positron.getVel(neg_time);

	ASSERT_ALMOST_EQUAL(pos[0], 0, eps);
	ASSERT_ALMOST_EQUAL(pos[1], 0, eps);
	ASSERT_ALMOST_EQUAL(pos[2], 3., eps);
	ASSERT_ALMOST_EQUAL(vel[0], 0, eps);
	ASSERT_ALMOST_EQUAL(vel[1], 1., eps);
	ASSERT_ALMOST_EQUAL(vel[2], 0, eps);
}

TEST(test_advance_no_fields) {
	const double eps = 1e-12;
	const int nr_steps = 2;

	const array<double, 3> startPos = {0, 0, 3.};
	const array<double, 3> startVel = {0, 1., 0};

	const array<double, 3> E = {0, 0, 0.};
	const array<double, 3> B = {0, 0., 0};

	const double dt = 1;

	const double m = 1;
	const double q = 1;

	Particle positron(m, q, startPos, startVel, nr_steps, dt);
	positron.advance(E, B);

	double time = dt;

	array<double, 3> pos = positron.getPos(time);
	array<double, 3> vel = positron.getVel(time);

	ASSERT_ALMOST_EQUAL(pos[0], 0, eps);
	ASSERT_ALMOST_EQUAL(pos[1], 1, eps);
	ASSERT_ALMOST_EQUAL(pos[2], 3., eps);
	ASSERT_ALMOST_EQUAL(vel[0], 0, eps);
	ASSERT_ALMOST_EQUAL(vel[1], 1., eps);
	ASSERT_ALMOST_EQUAL(vel[2], 0, eps);
}

TEST(test_advance_E_field) {
	const double eps = 1e-12;
	const int nr_steps = 3;

	const array<double, 3> startPos = {0, 0, 0.};
	const array<double, 3> startVel = {0, 1., 0};

	const array<double, 3> E = {1, 2, -1.};
	const array<double, 3> B = {0, 0., 0};

	const double dt = 1;

	const double m = 1. / (9.11e-31);
	const double q = 1. / (1.6e-19);

	Particle test_part(m, q, startPos, startVel, nr_steps, dt);
	test_part.advance(E, B);
	test_part.advance(E, B);

	array<double, 3> pos1 = test_part.getPos(dt);
	array<double, 3> vel1 = test_part.getVel(dt/2);
	array<double, 3> pos2 = test_part.getPos(2*dt);
	array<double, 3> vel2 = test_part.getVel(3*dt/2);

	ASSERT_ALMOST_EQUAL(pos1[0], 1, eps);
	ASSERT_ALMOST_EQUAL(pos1[1], 3, eps);
	ASSERT_ALMOST_EQUAL(pos1[2], -1., eps);
	ASSERT_ALMOST_EQUAL(vel1[0], 1, eps);
	ASSERT_ALMOST_EQUAL(vel1[1], 3., eps);
	ASSERT_ALMOST_EQUAL(vel1[2], -1, eps);

	ASSERT_ALMOST_EQUAL(pos2[0], 3, eps);
	ASSERT_ALMOST_EQUAL(pos2[1], 8, eps);
	ASSERT_ALMOST_EQUAL(pos2[2], -3., eps);
	ASSERT_ALMOST_EQUAL(vel2[0], 2, eps);
	ASSERT_ALMOST_EQUAL(vel2[1], 5., eps);
	ASSERT_ALMOST_EQUAL(vel2[2], -2, eps);
}

TEST_MAIN()
