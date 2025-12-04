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
matrix target_func_real_l3(double, matrix,  matrix = NAN, matrix = NAN);
matrix Q_real_l3(matrix, matrix = NAN, matrix = NAN);
matrix ff4T(matrix, matrix = NAN, matrix = NAN);  // testowa funkcja celu dla lab4
matrix ball_motion_l4(double, matrix, matrix = NAN, matrix = NAN);  // równania ruchu piłki
matrix ff4R(matrix, matrix = NAN, matrix = NAN);  // funkcja celu dla problemu rzeczywistego lab4
void simulate_ball_flight(double, double);  // symulacja lotu piłki dla sprawdzenia

