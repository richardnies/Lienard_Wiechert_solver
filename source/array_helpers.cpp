#include "array_helpers.h"

#include <cmath>

using namespace std;

array<double, 3> array_add(array<double, 3> array1, array<double, 3> array2)
{
	array<double, 3> array_out;

	for (int i = 0; i < 3; i++)
		array_out[i] = array1[i] + array2[i];

	return array_out;
}

array<double, 3> array_substract(array<double, 3> array1, array<double, 3> array2)
{
	array<double, 3> array_out;

	for (int i = 0; i < 3; i++)
		array_out[i] = array1[i] - array2[i];

	return array_out;
}

double array_norm(array<double, 3> arr)
{
	return sqrt( array_norm_sq(arr) );
}

double array_norm_sq(array<double, 3> arr)
{
	return pow(arr[0], 2) + pow(arr[1], 2) + pow(arr[2], 2);
}

double array_dot_product(array<double, 3> arr1, array<double, 3> arr2)
{
	return arr1[0]*arr2[0] + arr1[1]*arr2[1] + arr1[2]*arr2[2];
}

array<double, 3> array_cross_product(array<double, 3> arr1, array<double, 3> arr2)
{
	return { arr1[1] * arr2[2] - arr1[2] * arr2[1], 
					 arr1[2] * arr2[0] - arr1[0] * arr2[2], 
					 arr1[0] * arr2[1] - arr1[1] * arr2[0]  };
}

array<double, 3> array_mult_scalar(array<double, 3> arr, double mult)
{
	return { arr[0] * mult, arr[1] * mult, arr[2] * mult };
}
