#ifndef GUARD_array_helpers
#define GUARD_array_helpers 

#include <array>

std::array<double, 3> array_add(std::array<double, 3> array1, std::array<double, 3> array2);

std::array<double, 3> array_substract(std::array<double, 3> array1, std::array<double, 3> array2);

double array_norm(std::array<double, 3> arr);

double array_norm_sq(std::array<double, 3> arr);

double array_dot_product(std::array<double, 3> arr1, std::array<double, 3> arr2);

std::array<double, 3> array_cross_product(std::array<double, 3> arr1, std::array<double, 3> arr2);

std::array<double, 3> array_mult_scalar(std::array<double, 3> arr, double mult);


#endif