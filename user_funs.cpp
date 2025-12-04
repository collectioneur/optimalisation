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

// Lab 4 - testowa funkcja celu
matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
    matrix y(1, 1);
    double x1 = x(0);
    double x2 = x(1);
    
    // f(x1, x2) = sin(π*sqrt((x1/π)² + (x2/π)²)) / (π*sqrt((x1/π)² + (x2/π)²))
    double norm_squared = (x1/M_PI)*(x1/M_PI) + (x2/M_PI)*(x2/M_PI);
    double norm = sqrt(norm_squared);
    
    if (norm < 1e-10) {
        // Gdy norm → 0, sin(π*norm)/(π*norm) → 1
        y(0) = 1.0;
    } else {
        double argument = M_PI * norm;
        y(0) = sin(argument) / argument;
    }
    
    return y;
}

// Lab 4 - Problem rzeczywisty - równania ruchu piłki z efektem Magnusa
matrix ball_motion_l4(double t, matrix Y, matrix ud1, matrix ud2)
{
    // Y(0) = x - pozycja pozioma
    // Y(1) = y - pozycja pionowa  
    // Y(2) = vx - prędkość pozioma
    // Y(3) = vy - prędkość pionowa
    
    // ud2(0) = v0x - początkowa prędkość pozioma
    // ud2(1) = omega - rotacja piłki
    
    matrix dY(4, 1);
    
    // Parametry fizyczne
    double m = 0.6;        // masa piłki [kg]
    double r = 0.12;       // promień piłki [m]
    double g = 9.81;       // przyspieszenie ziemskie [m/s²]
    double C = 0.47;       // współczynnik oporu
    double rho = 1.2;      // gęstość powietrza [kg/m³]
    double S = M_PI * r * r;  // powierzchnia przekroju piłki
    double omega = ud2(1);    // rotacja piłki [rad/s]
    
    double x = Y(0);
    double y = Y(1);
    double vx = Y(2);
    double vy = Y(3);
    
    // Siły oporu powietrza
    double Dx = 0.5 * C * rho * S * vx * fabs(vx);
    double Dy = 0.5 * C * rho * S * vy * fabs(vy);
    
    // Siły Magnusa
    double F_Mx = rho * vy * omega * M_PI * r * r * r;
    double F_My = rho * vx * omega * M_PI * r * r * r;
    
    // Równania ruchu
    dY(0) = vx;  // dx/dt = vx
    dY(1) = vy;  // dy/dt = vy
    dY(2) = (-Dx - F_Mx) / m;  // m * dvx/dt = -Dx - F_Mx
    dY(3) = (-m * g - Dy - F_My) / m;  // m * dvy/dt = -mg - Dy - F_My
    
    return dY;
}

// Funkcja celu dla problemu rzeczywistego - maksymalizacja x_end
matrix ff4R(matrix params, matrix ud1, matrix ud2)
{
    matrix y(1, 1);
    
    double v0x = params(0);  // początkowa prędkość pozioma
    double omega = params(1);  // rotacja piłki
    
    // Warunki początkowe: x0=0, y0=100, vx0=v0x, vy0=0
    matrix Y0(4, 1);
    Y0(0) = 0.0;    // x0
    Y0(1) = 100.0;  // y0
    Y0(2) = v0x;    // vx0
    Y0(3) = 0.0;    // vy0
    
    matrix motion_params(2, 1);
    motion_params(0) = v0x;
    motion_params(1) = omega;
    
    // Symulacja lotu piłki
    double t0 = 0.0;
    double dt = 0.01;
    double t_end = 7.0;
    
    matrix* Y = solve_ode(ball_motion_l4, t0, dt, t_end, Y0, NAN, motion_params);
    
    int n = get_len(Y[0]);
    double x_end = 0.0;
    double x_at_y50 = 0.0;
    bool found_y50 = false;
    
    // Znajdź x_end (miejsce gdzie piłka uderza w ziemię) i x dla y=50
    for (int i = 1; i < n; i++) {
        double y_curr = Y[1](i, 1);  // aktualna pozycja y
        double y_prev = Y[1](i-1, 1);  // poprzednia pozycja y
        double x_curr = Y[1](i, 0);  // aktualna pozycja x
        double x_prev = Y[1](i-1, 0);  // poprzednia pozycja x
        
        // Sprawdź przecięcie z y = 50
        if (!found_y50 && y_prev >= 50.0 && y_curr <= 50.0) {
            // Interpolacja liniowa
            double t_interp = (50.0 - y_prev) / (y_curr - y_prev);
            x_at_y50 = x_prev + t_interp * (x_curr - x_prev);
            found_y50 = true;
        }
        
        // Sprawdź uderzenie w ziemię (y <= 0)
        if (y_curr <= 0.0 && y_prev > 0.0) {
            // Interpolacja liniowa do znalezienia dokładnego miejsca uderzenia
            double t_interp = (0.0 - y_prev) / (y_curr - y_prev);
            x_end = x_prev + t_interp * (x_curr - x_prev);
            break;
        }
    }
    
    // Funkcja celu: maksymalizujemy x_end, więc minimalizujemy -x_end
    // Dodajemy karę za naruszenie ograniczenia x ∈ [3, 7] dla y = 50
    double penalty = 0.0;
    if (found_y50) {
        if (x_at_y50 < 3.0) {
            penalty = 1000.0 * (3.0 - x_at_y50) * (3.0 - x_at_y50);
        } else if (x_at_y50 > 7.0) {
            penalty = 1000.0 * (x_at_y50 - 7.0) * (x_at_y50 - 7.0);
        }
    } else {
        penalty = 10000.0;  // Duża kara jeśli nie znaleziono y=50
    }
    
    y(0) = -x_end + penalty;  // Minimalizujemy -x_end + penalty
    
    delete[] Y;
    
    return y;
}

// Funkcja do sprawdzenia poprawności implementacji
void simulate_ball_flight(double v0x, double omega)
{
    cout << "=== Symulacja lotu piłki ===" << endl;
    cout << "v0x = " << v0x << " m/s, omega = " << omega << " rad/s" << endl;
    
    matrix Y0(4, 1);
    Y0(0) = 0.0;    // x0
    Y0(1) = 100.0;  // y0
    Y0(2) = v0x;    // vx0
    Y0(3) = 0.0;    // vy0
    
    matrix motion_params(2, 1);
    motion_params(0) = v0x;
    motion_params(1) = omega;
    
    double t0 = 0.0;
    double dt = 0.01;
    double t_end = 7.0;
    
    matrix* Y = solve_ode(ball_motion_l4, t0, dt, t_end, Y0, NAN, motion_params);
    
    int n = get_len(Y[0]);
    double x_end = 0.0;
    double x_at_y50 = 0.0;
    double t_at_y50 = 0.0;
    bool found_y50 = false;
    
    // Znajdź wyniki
    for (int i = 1; i < n; i++) {
        double t_curr = Y[0](i);
        double y_curr = Y[1](i, 1);
        double y_prev = Y[1](i-1, 1);
        double x_curr = Y[1](i, 0);
        double x_prev = Y[1](i-1, 0);
        double t_prev = Y[0](i-1);
        
        // Sprawdź przecięcie z y = 50
        if (!found_y50 && y_prev >= 50.0 && y_curr <= 50.0) {
            double t_interp = (50.0 - y_prev) / (y_curr - y_prev);
            x_at_y50 = x_prev + t_interp * (x_curr - x_prev);
            t_at_y50 = t_prev + t_interp * (t_curr - t_prev);
            found_y50 = true;
        }
        
        // Sprawdź uderzenie w ziemię
        if (y_curr <= 0.0 && y_prev > 0.0) {
            double t_interp = (0.0 - y_prev) / (y_curr - y_prev);
            x_end = x_prev + t_interp * (x_curr - x_prev);
            double t_end_actual = t_prev + t_interp * (t_curr - t_prev);
            cout << "x_end ≈ x(" << t_end_actual << "s) ≈ " << x_end << " m" << endl;
            break;
        }
    }
    
    if (found_y50) {
        cout << "x ≈ " << x_at_y50 << " m dla y ≈ 50 m (t ≈ " << t_at_y50 << " s)" << endl;
    }
    
    delete[] Y;
}