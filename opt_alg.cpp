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

int expansion_calls = 0; 

double* expansion(double(*ff)(double), double x0, double d, double alpha, int Nmax)
{
    try {
        double* p = new double[2]{0, 0};
        int i = 0;
        expansion_calls = 0;
        std::vector<double> X = {x0, x0 + d};

        double f0 = ff(X[0]); expansion_calls++;
        double f1 = ff(X[1]); expansion_calls++;

        if (f1 == f0) {
            p[0] = X[0];
            p[1] = X[1];
            return p;
        }

        if (f1 > f0) {
            d = -d;
            X[1] = X[0] + d;
            f1 = ff(X[1]); expansion_calls++;
            if (f1 >= f0) {
                p[0] = X[1];
                p[1] = X[0] - d;
                return p;
            }
        }

        while (f1 < f0) {
            if (i >= Nmax) break;
            i++;
            X.push_back(X[0] + pow(alpha, i) * d);
            f0 = f1;
            f1 = ff(X[i + 1]); expansion_calls++;
        }

        if (d > 0) {
            p[0] = X[i - 1];
            p[1] = X[i + 1];
        } else {
            p[0] = X[i + 1];
            p[1] = X[i - 1];
        }

        return p;
    }
    catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

int fib_calls = 0;

double* fib(double(*ff)(double), double a, double b, double epsilon)
{
    try {
        std::vector<unsigned long long> F = {0, 1};
        while (F.back() < static_cast<unsigned long long>((b - a) / epsilon))
            F.push_back(F[F.size() - 1] + F[F.size() - 2]);
        int N = static_cast<int>(F.size()) - 1;

        double x1 = a + (double)F[N - 2] / F[N] * (b - a);
        double x2 = a + (double)F[N - 1] / F[N] * (b - a);

        double f1 = ff(x1); fib_calls++;
        double f2 = ff(x2); fib_calls++;

        for (int k = 1; k <= N - 2; ++k) {
            if (f1 > f2) {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + (double)F[N - k - 1] / F[N - k] * (b - a);
                f2 = ff(x2); fib_calls++;
            } else {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + (double)F[N - k - 2] / F[N - k] * (b - a);
                f1 = ff(x1); fib_calls++;
            }
        }

        double* xmin_val = new double((x1 + x2) / 2.0);
        return xmin_val;
    }
    catch (string ex_info) {
        throw ("double* fib(...):\n" + ex_info);
    }
}

int lag_calls = 0;
double* lag(double(*ff)(double), double a, double b, double epsilon, double gamma, int Nmax)
{
    try {
        lag_calls = 0;
        int i = 0;
        double a_i = a, b_i = b;
        double c_i = (a + b) / 2.0;
        double d_i = 0.0, d_prev = 0.0;

        if (!(a_i < c_i && c_i < b_i))
            throw string("Nieprawidlowy przedzial poczatkowy");

        do {
            double f_a = ff(a_i); lag_calls++;
            double f_b = ff(b_i); lag_calls++;
            double f_c = ff(c_i); lag_calls++;

            double numerator = ((c_i - a_i) * (c_i - a_i)) * (f_c - f_b)
                             - ((c_i - b_i) * (c_i - b_i)) * (f_c - f_a);
            double denominator = 2.0 * ((c_i - a_i) * (f_c - f_b)
                             - (c_i - b_i) * (f_c - f_a));

            if (fabs(denominator) < 1e-12)
                break;

            d_prev = d_i;
            d_i = c_i - numerator / denominator;
			
            if (d_i <= a_i) d_i = a_i + 1e-8;
            if (d_i >= b_i) d_i = b_i - 1e-8;

            double f_d = ff(d_i); lag_calls++;

            if (d_i < c_i) {
                if (f_d < f_c) b_i = c_i;
                else a_i = d_i;
            } else {
                if (f_d < f_c) a_i = c_i;
                else b_i = d_i;
            }

            c_i = d_i;
            i++;

            if (lag_calls > Nmax)
                throw string("Przekroczono maksymalna liczbe wywolan funkcji celu (Nmax)");

        } while ((b_i - a_i) >= epsilon && (i == 1 || fabs(d_i - d_prev) >= gamma));

        double* result = new double(d_i);
        return result;
    }
    catch (string ex_info) {
        throw ("double* lag(...):\n" + ex_info);
    }
}




solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xb;
		matrix x = x0;
		//Tu wpisz kod funkcji
		while (s > epsilon)
		{
			Xb = x;
			x = HJ_trial(ff, Xb, s, ud1, ud2).x;

			//cout << "HJ expansion step: " << x << " f: " << ff(x, ud1, ud2) << endl;
			//cout << "Current step size s: " << s << endl;

			if (ff(x, ud1, ud2) < ff(Xb.x, ud1, ud2))
			{
				while (ff(x, ud1, ud2) < ff(Xb.x, ud1, ud2) && s > epsilon)
				{
					matrix xb_ = Xb.x;
					Xb = x;
					x = 2 * Xb.x - xb_;
					x = HJ_trial(ff, x, s, ud1, ud2).x;
					if (solution::f_calls > Nmax)
					{
						throw string("Przekroczono maksymalna liczbe wywolan funkcji celu (Nmax)");
					}
				}
				x = Xb.x;
			}
			else{
				s *= alpha;
			}
			if (solution::f_calls > Nmax)
			{
				throw string("Przekroczono maksymalna liczbe wywolan funkcji celu (Nmax)");
			}

		}
		

		return Xb;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

matrix dir[4] = {
	matrix(2, new double[2]{0, 1}),
	matrix(2, new double[2]{0, -1}),
	matrix(2, new double[2]{1, 0}),
	matrix(2, new double[2]{-1, 0})
};

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji
		for (int j = 0; j < 4; ++j)
		{
			if (ff(XB.x + dir[j] * s, ud1, ud2) < ff(XB.x, ud1, ud2))
			{
				XB.x = XB.x + dir[j] * s;
			}
			else if (ff(XB.x - dir[j] * s, ud1, ud2) < ff(XB.x, ud1, ud2))
			{
				XB.x = XB.x - dir[j] * s;
			}
		}
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
		
		int n = 2;
		
		int i = 0;
		matrix d = ident_mat(n);
		matrix lambda(n, 1, 0.0);
		matrix p(n, 1, 0.0);
		matrix s = s0;
		matrix s_initial = s0;
		matrix xB = x0;
		
		solution temp_sol;
		temp_sol.x = xB;
		temp_sol.fit_fun(ff, ud1, ud2);
		double f_xB = temp_sol.y(0);

		do {
			for (int j = 0; j < n; j++) {
				matrix test_point = xB;
				for (int k = 0; k < n; k++) {
					test_point(k) += s(j) * d(k, j);
				}
				
				temp_sol.x = test_point;
				temp_sol.fit_fun(ff, ud1, ud2);
				double f_test = temp_sol.y(0);
				
				if (f_test < f_xB) {
					xB = test_point;
					f_xB = f_test;
					lambda(j) += s(j);
					s(j) = alpha * s(j);
				} else {
					s(j) = -beta * s(j);
					p(j) += 1;
				}
			}

			i++;
			
			double max_step = 0.0;
			for (int j = 0; j < n; j++) {
				if (std::fabs(s(j)) > max_step) {
					max_step = std::fabs(s(j));
				}
			}
			
			if (solution::f_calls > Nmax) {
				Xopt.x = xB;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				return Xopt;
			}
			
			if (max_step < epsilon) {
				break;
			}
			
		} while (i < 100);
		

		Xopt.x = xB;
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.flag = 1;
		
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