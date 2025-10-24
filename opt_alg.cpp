#include"opt_alg.h"
#include <vector>
#include <cmath>
#include <limits>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(double(*ff)(double), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		//Tu wpisz kod funkcji

		int i = 0;
		vector<double> X = { x0, x0 + d };

		double f0 = ff(X[0]);
		double f1 = ff(X[1]);

		if (f1 == f0) {
			p[0] = X[0];
			p[1] = X[1];
			return p;
		}

		if (f1 > f0) {
			d = -d;
			X[1] = X[0] + d;
			f1 = ff(X[1]);
			if (f1 >= f0) {
				p[0] = X[1];
				p[1] = X[0] - d;
				return p;
			}
		}

		while (f1 < f0) {
			if (i >= Nmax) {
				throw ("Przekroczono maksymalna liczbe iteracji");
			}
			
			i += 1;
			X.push_back(X[0] + pow(alpha, i) * d);

			f0 = f1;
			f1 = ff(X[i + 1]);
		}

		if (d > 0)
			p[0] = X[i - 1], p[1] = X[i + 1];
		else
			p[0] = X[i + 1], p[1] = X[i - 1];

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

double* fib(double(*ff)(double), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		std::vector<unsigned long long> F = {0, 1};
		//calculate fibonacci numbers based on a, b and epsilon 
		while (F.back() < static_cast<unsigned long long>((b - a) / epsilon))
			F.push_back(F[F.size() - 1] + F[F.size() - 2]);
		int N = static_cast<int>(F.size()) - 1;
		
		int iterations = 0;

		double x1 = a + (double)F[N - 2] / F[N] * (b - a);
		double x2 = a + (double)F[N - 1] / F[N] * (b - a);

		double f1 = ff(x1);
		double f2 = ff(x2);

		//minimize range [a, b] using fibonacci numbers
		for (int k = 1; k <= N - 2; ++k)
    	{
			if (f1 > f2)
			{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + (double)F[N - k - 1] / F[N - k] * (b - a);
				f2 = ff(x2);
			}
			else
			{
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + (double)F[N - k - 2] / F[N - k] * (b - a);
				f1 = ff(x1);
			}
			iterations++;
		}
		cout<<"Number of iterations: "<< iterations << endl;
		double* xmin_val = new double((x1 + x2) / 2.0);
		return xmin_val;
	}
	catch (string ex_info)
	{
			throw ("double* fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		
		int iterations = 0;
		int i = 0;
		double a_i = a, b_i = b;
		double c_i = (a + b) / 2.0;
		double d_i = 0.0, d_prev = 0.0;
		
		if (!(a_i < c_i && c_i < b_i)) {
			Xopt.flag = 0;
			return Xopt;
		}
		
		do {

			matrix x_a(1, 1), x_b(1, 1), x_c(1, 1);
			x_a(0) = a_i;
			x_b(0) = b_i;
			x_c(0) = c_i;
			
			matrix f_a = ff(x_a, ud1, ud2);
			matrix f_b = ff(x_b, ud1, ud2);
			matrix f_c = ff(x_c, ud1, ud2);
			
			double l = f_a(0) * (b_i * b_i - c_i * c_i) + 
					   f_b(0) * (c_i * c_i - a_i * a_i) + 
					   f_c(0) * (a_i * a_i - b_i * b_i);
			
			double m = f_a(0) * (b_i - c_i) + 
					   f_b(0) * (c_i - a_i) + 
					   f_c(0) * (a_i - b_i);
			
			if (m <= 0) {
				Xopt.flag = 0;
				return Xopt;
			}
			
			d_prev = d_i;
			

			d_i = 0.5 * l / m;
			
			if (a_i < d_i && d_i < c_i) {
				matrix x_d(1, 1);
				x_d(0) = d_i;
				matrix f_d = ff(x_d, ud1, ud2);
				
				if (f_d(0) < f_c(0)) {
					double new_a = a_i;
					double new_c = d_i;
					double new_b = c_i;
					
					a_i = new_a;
					c_i = new_c;
					b_i = new_b;
				} else {
					double new_a = d_i;
					double new_c = c_i;
					double new_b = b_i;
					
					a_i = new_a;
					c_i = new_c;
					b_i = new_b;
				}
			} else if (c_i < d_i && d_i < b_i) {
				matrix x_d(1, 1);
				x_d(0) = d_i;
				matrix f_d = ff(x_d, ud1, ud2);
				
				if (f_d(0) < f_c(0)) {
					double new_a = c_i;
					double new_c = d_i;
					double new_b = b_i;
					
					a_i = new_a;
					c_i = new_c;
					b_i = new_b;
				} else {
					double new_a = a_i;
					double new_c = c_i;
					double new_b = d_i;
					
					a_i = new_a;
					c_i = new_c;
					b_i = new_b;
				}
			} else {
				Xopt.flag = 0;
				return Xopt;
			}
			
			i++;
			
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
			
			iterations++;
		} while ((b_i - a_i) >= epsilon && (i == 1 || std::fabs(d_i - d_prev) >= gamma)); 
		
		
		Xopt.x = matrix(1, 1);
		Xopt.x(0) = d_i;
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.flag = 1; 

		cout<<"Number of iterations: "<< iterations << endl;
		
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
