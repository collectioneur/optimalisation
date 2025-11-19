/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include "ode_solver.h"
#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();
void lab1_all();

int main()
{
	try
	{
		lab3();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}
void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;						 					// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1() {
    ofstream out("wyniki_lab1.csv");
    out << "i;alpha;x0;a;b;it_exp;x_fib;f_fib;it_fib;min_fib;x_lag;f_lag;it_lag;min_lag\n";

    double eps = 1e-4;
    int Nmax = 1000;
    double alphas[3] = {1.2, 1.5, 2.0};

    srand(time(NULL));

    for (int k = 0; k < 3; ++k) {
        double alpha = alphas[k];

        for (int i = 0; i < 100; ++i) {
            double x0 = -100 + (rand() % 200);

            //metoda ekspansji
            double* p = expansion(ff1, x0, 1.0, alpha, Nmax);
            extern int expansion_calls;
            int iter_exp = expansion_calls;

            //metoda Fibonacciego
            double* fib_res = fib(ff1, p[0], p[1], eps);
            extern int fib_calls;
            int iter_fib = fib_calls;
            double f_fib = ff1(*fib_res);

            //metoda Lagrange’a (nowa wersja double*)
            double* lag_res = lag(ff1, p[0], p[1], eps, 1e-5, Nmax);
            extern int lag_calls;
            int iter_lag = lag_calls;
            double f_lag = ff1(*lag_res);

            //klasyfikacja minimum 
            string min_fib_type = (fabs(f_fib + 0.9211) < 1e-3) ? "globalne" : "lokalne";
            string min_lag_type = (fabs(f_lag + 0.9211) < 1e-3) ? "globalne" : "lokalne";

            //zapis do csv
            out << i + 1 << ";"
                << alpha << ";"
                << x0 << ";"
                << p[0] << ";"
                << p[1] << ";"
                << iter_exp << ";"
                << *fib_res << ";"
                << f_fib << ";"
                << iter_fib << ";"
                << min_fib_type << ";"
                << *lag_res << ";"
                << f_lag << ";"
                << iter_lag << ";"
                << min_lag_type << "\n";

            delete[] p;
            delete fib_res;
            delete lag_res;
        }
    }
    out.close();
}

void lab1_all() {
    ofstream out("wyniki_lab1_all.csv");
    out << "Opis;alpha;a;b;x_fib;f_fib;it_fib;min_fib;x_lag;f_lag;it_lag;min_lag\n";

    double eps = 1e-4;
    int Nmax = 1000;
    double a = -100;
    double b = 100;
    double alphas[3] = {1.2, 1.5, 2.0};

    for (int k = 0; k < 3; ++k) {
        double alpha = alphas[k];
        cout << "=== ALPHA = " << alpha << " ===\n";

        //bez ekspansji
        double* fib_res_noexp = fib(ff1, a, b, eps);
        extern int fib_calls;
        int iter_fib_noexp = fib_calls;
        double f_fib_noexp = ff1(*fib_res_noexp);

        double* lag_res_noexp = lag(ff1, a, b, eps, 1e-5, Nmax);
        extern int lag_calls;
        int iter_lag_noexp = lag_calls;
        double f_lag_noexp = ff1(*lag_res_noexp);

        string min_fib_noexp = (fabs(f_fib_noexp + 0.9211) < 1e-3) ? "globalne" : "lokalne";
        string min_lag_noexp = (fabs(f_lag_noexp + 0.9211) < 1e-3) ? "globalne" : "lokalne";

        out << "Fibonacci bez ekspansji;" << alpha << ";"
            << a << ";" << b << ";"
            << *fib_res_noexp << ";" << f_fib_noexp << ";" << iter_fib_noexp << ";" << min_fib_noexp << ";"
            << *lag_res_noexp << ";" << f_lag_noexp << ";" << iter_lag_noexp << ";" << min_lag_noexp << "\n";

        cout << "Bez ekspansji: x_fib=" << *fib_res_noexp << ", f_fib=" << f_fib_noexp
             << ", x_lag=" << *lag_res_noexp << ", f_lag=" << f_lag_noexp << "\n";

        //z ekspansją
        double x0 = -100 + (rand() % 200);
        double* p = expansion(ff1, x0, 1.0, alpha, Nmax);
        double a_exp = p[0];
        double b_exp = p[1];

        double* fib_res_exp = fib(ff1, a_exp, b_exp, eps);
        extern int fib_calls;
        int iter_fib_exp = fib_calls;
        double f_fib_exp = ff1(*fib_res_exp);

        double* lag_res_exp = lag(ff1, a_exp, b_exp, eps, 1e-5, Nmax);
        extern int lag_calls;
        int iter_lag_exp = lag_calls;
        double f_lag_exp = ff1(*lag_res_exp);

        string min_fib_exp = (fabs(f_fib_exp + 0.9211) < 1e-3) ? "globalne" : "lokalne";
        string min_lag_exp = (fabs(f_lag_exp + 0.9211) < 1e-3) ? "globalne" : "lokalne";

        out << "Fibonacci z ekspansją;" << alpha << ";"
            << a_exp << ";" << b_exp << ";"
            << *fib_res_exp << ";" << f_fib_exp << ";" << iter_fib_exp << ";" << min_fib_exp << ";"
            << *lag_res_exp << ";" << f_lag_exp << ";" << iter_lag_exp << ";" << min_lag_exp << "\n";

        cout << "Z ekspansją (α=" << alpha << "): a=" << a_exp << ", b=" << b_exp
             << ", x_fib=" << *fib_res_exp << ", f_fib=" << f_fib_exp
             << ", x_lag=" << *lag_res_exp << ", f_lag=" << f_lag_exp << "\n\n";

        delete fib_res_noexp;
        delete lag_res_noexp;
        delete fib_res_exp;
        delete lag_res_exp;
        delete p;
    }

    out.close();
    cout << "=== Zapisano wyniki do wyniki_lab1_all.csv ===\n";
}

void lab2()
{
    double x0 = 20.0;
    double d = 5.0;
    double alpha = 1.5;
    double epsilon = 1e-2;

    double* interval = expansion(target_f_l2, x0, d, alpha, 100);
    cout << "Found interval: [" << interval[0] << ", " << interval[1] << "]" << endl;

    double* fibo = fib(target_f_l2, interval[0], interval[1], epsilon);
    cout << "Fibonacci: D_A* = " << *fibo << ", y* = " << target_f_l2(*fibo) << endl;

    double* lagr = lag(target_f_l2, interval[0], interval[1], epsilon, 1e-5, 1000);
    cout << "Lagrange: D_A* = " << *lagr << ", y* = " << target_f_l2(*lagr) << endl;

    cout << "Simulation with optimal D_A = " << *lagr << endl;
    f_l2_print(*lagr);

    delete[] interval;
    delete fibo;
    delete lagr;
}

void lab3()
{
    matrix x0(2, new double[2]{1.0, 0.0});
    matrix s0(2, new double[2]{1.0, 1.0});
    double alpha = 1.2, beta = 0.8, epsilon = 1e-4;
    int Nmax = 200;

    //rosenbrock testowa funkcja - podobnie zrobić z metodą HJ
    solution::clear_calls();
    solution result = Rosen(target_func_l3, x0, s0, alpha, beta, epsilon, Nmax);

    
    cout << "Wywołanie testowe dla metody Rosenbrocka" << endl;
    cout << "x = [" << result.x(0) << ", " << result.x(1) << "]" << endl;
    cout << "f(x) = " << result.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl << endl;


    //test porpawności dla funkcji Q(Q_real_l3)
    cout << "Test poprawności dla Q" << endl;

    matrix test_k(2, new double[2]{5.0, 5.0});
    matrix Q = Q_real_l3(test_k, NAN, test_k);

    cout << "Q(k1=5, k2=5) = " << Q(0) << endl;

    cout << "Oczekiwana wartość(ze sprawka): 775.229" << endl << endl;

    //rozwozanie problemu rzeczywistego -suzkanie najlepszego k1 i k2 - metoda rosenbrocka - podobnie z HJ
    cout << "Optymalizacja problemu rzeczywistego metoda Rosenbrocka" << endl;
    matrix k0(2, new double[2]{10.0, 10.0});
    matrix s0_real(2, new double[2]{2.0, 2.0});
    epsilon = 1e-2;
    Nmax = 100;

    solution::clear_calls();
    solution result_real = Rosen(Q_real_l3, k0, s0_real, alpha, beta, epsilon, Nmax);

    cout << "k1 = " << result_real.x(0) << endl;
    cout << "k2 = " << result_real.x(1) << endl;
    cout << "Q(k1,k2) = " << result_real.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl << endl;

    //symulcja z czasem z otpymalnymi parametrami
    cout << "Symulacja z optymalnymi parametrami" << endl;

    matrix Y0(2,1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix k_opt(2,1);
    k_opt(0) = result_real.x(0);
    k_opt(1) = result_real.x(1);

    matrix* Y = solve_ode(target_func_real_l3, 0, 0.1, 100, Y0, NAN, k_opt);

    ofstream Sout("symulacja_lab3.csv");
    Sout << "t,alpha,omega\n";  
    int n = get_len(Y[0]);
    
    for (int i = 0; i < n; ++i)
    {
        Sout << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 1) << "\n";
    }
    Sout.close();
    cout << "zapisano do symulacja_lab3.csv" << endl;

    delete[] Y;
}


void lab4()
{
	cout << HJ(target_func_l3, matrix(2, new double[2]{-0.1, 0.2}), 1.0, 0.3, 1e-2, 1000, NAN, NAN).x << endl;
}

void lab5()
{

}

void lab6()
{

}