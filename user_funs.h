#pragma once

#include "ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
double ff1(double);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix l2_dvdt(double, matrix, matrix, matrix = NAN);
double target_f_l2(double);
void f_l2_print(double);
matrix target_func_l3(matrix, matrix = NAN, matrix = NAN);
matrix target_func_real_l3(double, matrix, matrix = NAN, matrix = NAN);
matrix Q_real_l3(matrix, matrix = NAN, matrix = NAN);
matrix ff4T(matrix, matrix = NAN, matrix = NAN);                   // testowa funkcja celu dla lab4
matrix ball_motion_l4(double, matrix, matrix = NAN, matrix = NAN); // równania ruchu piłki
matrix ff4R(matrix, matrix = NAN, matrix = NAN);                   // funkcja celu dla problemu rzeczywistego lab4
void simulate_ball_flight(double, double);                         // symulacja lotu piłki dla sprawdzenia
// --- Funkcje pomocnicze do Lab 4 (Funkcje Kary) ---
// Zaimplementowane w opt_alg.cpp, wywoływane w main.cpp
matrix penalty_objective_function(matrix x, matrix ud1, matrix ud2);
matrix penalty_objective_function_lab4R(matrix x, matrix ud1, matrix ud2);

// --- Funkcje do Lab 5 (Metody Gradientowe) ---
// Funkcja testowa
matrix ff5T(matrix x, matrix ud1, matrix ud2); // Funkcja celu
matrix gf5T(matrix x, matrix ud1, matrix ud2); // Gradient
matrix Hf5T(matrix x, matrix ud1, matrix ud2); // Hesjan
matrix ff5T_multi(matrix x, matrix ud1, matrix ud2); // Funkcja wielokryterialna testowa

// Problem rzeczywisty (Regresja logistyczna)
matrix ff5R(matrix theta, matrix ud1, matrix ud2); // Funkcja kosztu
matrix gf5R(matrix theta, matrix ud1, matrix ud2); // Gradient kosztu
matrix hypothesis(matrix theta, matrix X);

// Funkcje pomocnicze do wczytywania danych (Lab 5)
// Wymaga #include <string> na początku pliku opt_alg.h
matrix read_matrix_from_file(std::string filename, int rows, int cols);

matrix ff6T(matrix x, matrix ud1, matrix ud2);
