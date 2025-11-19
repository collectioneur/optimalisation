#include"user_funs.h"
#include<cmath>

#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

double ff1(double x)				// funkcja celu dla przypadku testowego 1D
{
	return -cos(0.1 * x) * exp(-pow(0.1 * x - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x, 2);
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)		
{
	matrix y;
	double val = x(0);	
	y = -cos(0.1 * val) * exp(-pow(0.1 * val - 2 * M_PI, 2)) + 0.002 * pow(0.1 * val, 2);
	return y;
}
//

// Lab 2
matrix l2_dvdt(double t,matrix Y, matrix ud1, matrix ud2)
{
	double P_A = 2;
	double V_A = Y(0);
	double P_B = 1;
	double V_B = Y(1);
	double T_A = 95;
	double T_B = Y(2);
	double T_B_in = 20;
	double F_B_in = 0.01; // m3/s
	double D_a = ud1(0);
	double D_b = 0.00365665;

	double a = 0.98;
	double b = 0.63;
	double g = 9.81;

	matrix dY = matrix(3, 1);
	dY(0) = -a*b*D_a*sqrt(2*g*(V_A)/P_A);
	dY(1) = -a*b*D_b*sqrt(2*g*(V_B)/P_B) + a*b*D_a*sqrt(2*g*(V_A)/P_A) + F_B_in;

	double F_A_out = a*b*D_a*sqrt(2*g*(V_A)/P_A);
	dY(2) = F_B_in*(T_B_in - T_B)/V_B + F_A_out*(T_A - T_B)/V_B;
	return dY;
}

matrix f_l2(double x, double t) {
	x *= 0.0001;
	return solve_ode(l2_dvdt, 0, 1, t, matrix(3, new double[3]{5, 1, 20.0}), matrix(1, new double[1]{x}))[1];
}

double f_l2_max(double x)
{
	int t = 500;
	matrix res = f_l2(x, t);
	double ans = 0.0;
	for(int i = 0; i < t; ++i)
	{
		if(res(i, 2) > ans) {
			ans = res(i, 2);
		}
	}
	return ans;
}

void f_l2_print(double x)
{
	int t = 500;
	cout << f_l2(x, t) << endl;
	return;
}

double target_f_l2(double x)
{
	double ans = f_l2_max(x);

	double target_temp = 50.0;
	ans = pow(ans - target_temp, 2);
	return ans;
}

matrix target_func_l3(matrix x, matrix ud1, matrix ud2)
{
	matrix y(1,1);
	double x1 = x(0);
	double x2 = x(1);
	y(0) = x1*x1 + x2*x2 - cos(2.5*M_PI*x1) - cos(2.5*M_PI*x2) + 2.0;
	return y;
}

matrix target_func_real_l3(double t, matrix Y, matrix ud1, matrix ud2)
{
    matrix dY(2,1);

    double alpha = Y(0);
    double omega = Y(1);
    double l = 2.0;
    double m_r = 1.0;  
    double m_w = 5.0;  
    double b = 0.25;

    double I = (1.0/3.0) * m_r * l * l + m_w * l * l;

    double k1 = ud2(0);
    double k2 = ud2(1);

    double alpha_ref = M_PI;
    double omega_ref = 0.0;

    double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

    dY(0) = omega;
    dY(1) = (M - b * omega) / I;

    return dY;
}

matrix Q_real_l3(matrix x, matrix ud1, matrix ud2)
{	
    matrix y(1, 1);
    
    double k1 = x(0);
    double k2 = x(1);
	
    double t0 = 0.0;
    double dt = 0.1;
    double t_end = 100.0;
    
    matrix Y0(2, 1);
    Y0(0) = 0.0;  
    Y0(1) = 0.0;  
    
    matrix k_params(2, 1);
    k_params(0) = k1;
    k_params(1) = k2;

    matrix* Y = solve_ode(target_func_real_l3, t0, dt, t_end, Y0, NAN, k_params);
    
    double alpha_ref = M_PI;
    double omega_ref = 0.0;
	
    int n = get_len(Y[0]);
    double Q = 0.0;
    
    for (int i = 0; i < n; ++i)
    {
        double alpha_t = Y[1](i, 0);  
        double omega_t = Y[1](i, 1);  
        
        double M_t = k1 * (alpha_ref - alpha_t) + k2 * (omega_ref - omega_t);
        
        double term1 = 10.0 * pow(alpha_ref - alpha_t, 2);
        double term2 = pow(omega_ref - omega_t, 2);
        double term3 = pow(M_t, 2);

        Q += (term1 + term2 + term3) * dt;
    }
    
    y(0) = Q;

    delete[] Y;
    
    return y;
}
	

