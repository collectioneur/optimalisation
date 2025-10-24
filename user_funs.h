#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
double ff1(double);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix l2_dvdt(double, matrix, matrix, matrix = NAN);
double target_f_l2(double);
void f_l2_print(double);
matrix target_func_l3(matrix, matrix = NAN, matrix = NAN);

