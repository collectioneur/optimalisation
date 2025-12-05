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
		lab5();
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

    //rosenbrock testowa funkcja
    solution::clear_calls();
    solution result = Rosen(target_func_l3, x0, s0, alpha, beta, epsilon, Nmax);

    
    cout << "Wywołanie testowe dla metody Rosenbrocka" << endl;
    cout << "x = [" << result.x(0) << ", " << result.x(1) << "]" << endl;
    cout << "f(x) = " << result.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl << endl;

    //test HJ dla funkcji testowej
    solution::clear_calls();
    solution result_HJ_test = HJ(target_func_l3, x0, 1.0, 0.3, epsilon, Nmax);
    
    cout << "Wywołanie testowe dla metody Hooke-Jeeves" << endl;
    cout << "x = [" << result_HJ_test.x(0) << ", " << result_HJ_test.x(1) << "]" << endl;
    cout << "f(x) = " << result_HJ_test.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl << endl;


    //test porpawności dla funkcji Q(Q_real_l3)
    cout << "Test poprawności dla Q" << endl;

    matrix test_k(2, new double[2]{5.0, 5.0});
    matrix Q = Q_real_l3(test_k, NAN, test_k);

    cout << "Q(k1=5, k2=5) = " << Q(0) << endl;

    cout << "Oczekiwana wartość(ze sprawka): 775.229" << endl << endl;

    //rozwozanie problemu rzeczywistego -suzkanie najlepszego k1 i k2 - metoda rosenbrocka
    cout << "Optymalizacja problemu rzeczywistego metoda Rosenbrocka" << endl;
    matrix k0(2, new double[2]{10.0, 10.0});
    matrix s0_real(2, new double[2]{2.0, 2.0});
    epsilon = 1e-2;
    Nmax = 100;

    solution::clear_calls();
    solution result_real_rosen = Rosen(Q_real_l3, k0, s0_real, alpha, beta, epsilon, Nmax);

    cout << "k1 = " << result_real_rosen.x(0) << endl;
    cout << "k2 = " << result_real_rosen.x(1) << endl;
    cout << "Q(k1,k2) = " << result_real_rosen.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl << endl;

    //rozwozanie problemu rzeczywistego - metoda HJ
    cout << "Optymalizacja problemu rzeczywistego metoda Hooke-Jeeves" << endl;
    solution::clear_calls();
    solution result_real_HJ = HJ(Q_real_l3, k0, 2.0, 0.3, epsilon, Nmax);

    cout << "k1 = " << result_real_HJ.x(0) << endl;
    cout << "k2 = " << result_real_HJ.x(1) << endl;
    cout << "Q(k1,k2) = " << result_real_HJ.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl << endl;

    //symulcja z czasem z optymalnymi parametrami - metoda Rosenbrocka
    cout << "Symulacja z optymalnymi parametrami (Rosenbrock)" << endl;

    matrix Y0(2,1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix k_opt_rosen(2,1);
    k_opt_rosen(0) = result_real_rosen.x(0);
    k_opt_rosen(1) = result_real_rosen.x(1);

    matrix* Y_rosen = solve_ode(target_func_real_l3, 0, 0.1, 100, Y0, NAN, k_opt_rosen);

    //symulcja z czasem z optymalnymi parametrami - metoda HJ
    cout << "Symulacja z optymalnymi parametrami (Hooke-Jeeves)" << endl;

    matrix k_opt_HJ(2,1);
    k_opt_HJ(0) = result_real_HJ.x(0);
    k_opt_HJ(1) = result_real_HJ.x(1);

    matrix* Y_HJ = solve_ode(target_func_real_l3, 0, 0.1, 100, Y0, NAN, k_opt_HJ);

    //zapis wyników do pliku csv z obiema metodami
    ofstream Sout("symulacja_lab3.csv");
    Sout << "t,alpha_rosen,omega_rosen,alpha_HJ,omega_HJ\n";  
    int n_rosen = get_len(Y_rosen[0]);
    int n_HJ = get_len(Y_HJ[0]);
    int n = (n_rosen < n_HJ) ? n_rosen : n_HJ;
    
    for (int i = 0; i < n; ++i)
    {
        Sout << Y_rosen[0](i) << "," 
             << Y_rosen[1](i, 0) << "," << Y_rosen[1](i, 1) << ","
             << Y_HJ[1](i, 0) << "," << Y_HJ[1](i, 1) << "\n";
    }
    Sout.close();
    cout << "zapisano do symulacja_lab3.csv" << endl;

    delete[] Y_rosen;
    delete[] Y_HJ;
}


void lab4()
{
    // cout << "=== LAB 4 - Testowa funkcja celu z ograniczeniami ===" << endl << endl;
    
    // // Utworzenie pliku CSV dla wyników
    // ofstream results_file("lab4_results.csv");
    // results_file << "a,i,x1_start,x2_start,x1_ext,x2_ext,r_ext,y_ext,calls_ext,x1_int,x2_int,r_int,y_int,calls_int\n";
    
    // ofstream stats_file("lab4_statistics.csv");
    // stats_file << "a,penalty_type,avg_x1,avg_x2,avg_r,avg_y,avg_calls,success_rate\n";
    
    // // Parametry optymalizacji
    // double epsilon = 1e-3;
    // int Nmax = 10000;
    // double values_a[] = {4.0, 4.4934, 5.0};
    // int num_tests = 100; // Полные 100 тестов для каждого параметра a
    
    // srand(time(NULL));
    
    // // Test funkcji celu
    // cout << "Test funkcji celu f(x1, x2):" << endl;
    // matrix test_point(2, new double[2]{1.0, 1.0});
    // matrix f_test = ff4T(test_point, NAN, NAN);
    // cout << "f(1, 1) = " << f_test(0) << endl << endl;
    
    // // Przeprowadź 100 testów dla każdej wartości a
    // for (int k = 0; k < 3; k++) {
    //     double a = values_a[k];
    //     cout << "=== Testowanie dla a = " << a << " ===" << endl;
        
    //     // Statystyki dla zewnętrznej funkcji kary
    //     double sum_x1_ext = 0.0, sum_x2_ext = 0.0, sum_f_ext = 0.0, sum_calls_ext = 0.0, sum_r_ext = 0.0;
    //     int success_count_ext = 0;
        
    //     // Statystyki dla wewnętrznej funkcji kary
    //     double sum_x1_int = 0.0, sum_x2_int = 0.0, sum_f_int = 0.0, sum_calls_int = 0.0, sum_r_int = 0.0;
    //     int success_count_int = 0;
        
    //     for (int i = 0; i < num_tests; i++) {
    //         // Generowanie losowego punktu startowego w obszarze dopuszczalnym dla wewnętrznej funkcji kary
    //         matrix x0(2, 1);
    //         bool valid_point = false;
    //         int attempts = 0;
            
    //         while (!valid_point && attempts < 1000) {
    //             // Generuj punkt w większym obszarze dla bardziej różnorodnych punktów startowych
    //             double r = (a - 0.5) * ((double)rand() / RAND_MAX + 0.2); // Więcej różnorodności w promieniu
    //             double theta = 2.0 * M_PI * ((double)rand() / RAND_MAX);
    //             x0(0) = r * cos(theta) + 1.3 + 0.5 * ((double)rand() / RAND_MAX); // Większa losowość
    //             x0(1) = r * sin(theta) + 1.3 + 0.5 * ((double)rand() / RAND_MAX);
                
    //             // Sprawdź ograniczenia (punkt musi być bezpiecznie w środku dla wewnętrznej funkcji kary)
    //             if (x0(0) >= 1.1 && x0(1) >= 1.1 && sqrt(x0(0)*x0(0) + x0(1)*x0(1)) <= a - 0.2) {
    //                 valid_point = true;
    //             }
    //             attempts++;
    //         }
            
    //         if (!valid_point) {
    //             cout << "Nie można znaleźć poprawnego punktu startowego dla testu " << i+1 << endl;
    //             continue;
    //         }
            
    //         // Test zewnętrznej funkcji kary - bezpośrednio sym_NM
    //         matrix penalty_params_ext(3, 1);
    //         penalty_params_ext(0) = 0.0;  // zewnętrzna funkcja kary
    //         penalty_params_ext(1) = a;    // parametr a  
    //         penalty_params_ext(2) = 10.0; // współczynnik c (zmniejszony dla więcej iteracji)
            
    //         solution::clear_calls();
    //         matrix pen_params_ext(2, 1);
    //         pen_params_ext(0) = 0.0;  // внешняя функция кары
    //         pen_params_ext(1) = a;    // параметр a
    //         solution result_ext = pen(ff4T, x0, 1.0, 10.0, epsilon, Nmax, pen_params_ext, NAN);
            
    //         double r_ext = sqrt(result_ext.x(0)*result_ext.x(0) + result_ext.x(1)*result_ext.x(1));
    //         bool converged_ext = (result_ext.flag == 1);
    //         int calls_ext = solution::f_calls;
            
    //         // Aktualizuj statystyki zewnętrznej funkcji kary
    //         if (converged_ext) {
    //             sum_x1_ext += result_ext.x(0);
    //             sum_x2_ext += result_ext.x(1);
    //             sum_f_ext += result_ext.y(0);
    //             sum_calls_ext += calls_ext;
    //             sum_r_ext += r_ext;
    //             success_count_ext++;
    //         }
            
    //         // Test wewnętrznej funkcji kary - tym samym punktem startowym
    //         matrix penalty_params_int(3, 1);
    //         penalty_params_int(0) = 1.0;  // wewnętrzna funkcja kary
    //         penalty_params_int(1) = a;    // parametr a
    //         penalty_params_int(2) = 0.1;  // współczynnik c (zmniejszony dla więcej iteracji)
            
    //         solution::clear_calls();
    //         matrix pen_params_int(2, 1);
    //         pen_params_int(0) = 1.0;  
    //         pen_params_int(1) = a;    
    //         solution result_int = pen(ff4T, x0, 1.0, 0.1, epsilon, Nmax, pen_params_int, NAN);
            
    //         double r_int = sqrt(result_int.x(0)*result_int.x(0) + result_int.x(1)*result_int.x(1));
    //         bool converged_int = (result_int.flag == 1);
    //         int calls_int = solution::f_calls;
            
    //         // Aktualizuj statystyki wewnętrznej funkcji kary
    //         if (converged_int) {
    //             sum_x1_int += result_int.x(0);
    //             sum_x2_int += result_int.x(1);
    //             sum_f_int += result_int.y(0);
    //             sum_calls_int += calls_int;
    //             sum_r_int += r_int;
    //             success_count_int++;
    //         }
            
    //         // Zapisz jedną linię z obydwoma wynikami
    //         results_file << a << "," << i+1 << "," << x0(0) << "," << x0(1) << ",";
            
    //         if (converged_ext) {
    //             results_file << result_ext.x(0) << "," << result_ext.x(1) << "," << r_ext << "," 
    //                        << result_ext.y(0) << "," << calls_ext;
    //         } else {
    //             results_file << "NA,NA,NA,NA,NA";
    //         }
            
    //         results_file << ",";
            
    //         if (converged_int) {
    //             results_file << result_int.x(0) << "," << result_int.x(1) << "," << r_int << "," 
    //                        << result_int.y(0) << "," << calls_int;
    //         } else {
    //             results_file << "NA,NA,NA,NA,NA";
    //         }
            
    //         results_file << "\n";
            
    //         if ((i+1) % 20 == 0) {
    //             cout << "Ukończono " << i+1 << "/" << num_tests << " testów dla a=" << a << endl;
    //         }
    //     }
        
    //     // Zapisz statystyki do pliku
    //     if (success_count_ext > 0) {
    //         stats_file << a << ",external," 
    //                   << sum_x1_ext/success_count_ext << "," << sum_x2_ext/success_count_ext << ","
    //                   << sum_r_ext/success_count_ext << "," << sum_f_ext/success_count_ext << ","
    //                   << sum_calls_ext/success_count_ext << "," << (double)success_count_ext/num_tests << "\n";
    //     }
        
    //     if (success_count_int > 0) {
    //         stats_file << a << ",internal," 
    //                   << sum_x1_int/success_count_int << "," << sum_x2_int/success_count_int << ","
    //                   << sum_r_int/success_count_int << "," << sum_f_int/success_count_int << ","
    //                   << sum_calls_int/success_count_int << "," << (double)success_count_int/num_tests << "\n";
    //     }
        
    //     cout << "Ukończono wszystkie testy dla a=" << a << endl;
    //     cout << "Zewnętrzna funkcja kary: " << success_count_ext << "/" << num_tests << " udanych" << endl;
    //     cout << "Wewnętrzna funkcja kary: " << success_count_int << "/" << num_tests << " udanych\n" << endl;
    // }
    
    // // Test metody Nelder-Mead dla funkcji testowej bez ograniczeń
    // cout << "=== Test metody sympleks Nelder-Mead (bez ograniczeń) ===" << endl;
    // matrix x0_nm(2, new double[2]{2.0, 2.0});
    
    // solution::clear_calls();
    // solution result_nm = sym_NM(ff4T, x0_nm, 0.5, 1.0, 0.5, 2.0, 0.5, epsilon, Nmax, NAN, NAN);
    
    // cout << "Punkt startowy: [" << x0_nm(0) << ", " << x0_nm(1) << "]" << endl;
    // cout << "Wynik: x* = [" << result_nm.x(0) << ", " << result_nm.x(1) << "]" << endl;
    // cout << "f(x*) = " << result_nm.y(0) << endl;
    // cout << "Liczba wywołań: " << solution::f_calls << endl;
    // cout << "r = " << sqrt(result_nm.x(0)*result_nm.x(0) + result_nm.x(1)*result_nm.x(1)) << endl;
    
    // // Zamknij pliki CSV
    // results_file.close();
    // stats_file.close();
    // cout << "Zapisano wyniki do lab4_results.csv i lab4_statistics.csv\n" << endl;
    
    // === PROBLEM RZECZYWISTY ===
    cout << "=== PROBLEM RZECZYWISTY - Lot piłki z efektem Magnusa ===" << endl;

    // ---------------- Parametry metody funkcji kary ----------------
    const double penaltyC0           = 2.0;      // początkowy współczynnik kary
    const double penaltyScale        = 20.0;     // współczynnik skalowania c
    const double penaltyEps          = 1e-10;     // dokładność metody kary
    const int    penaltyMaxCalls     = 10000;     // maks. liczba wywołań funkcji celu
    const int    penaltyMaxOuterIter = 50;       // maks. liczba pętli zewnętrznych (aktualizacji c)

    // Punkt startowy: [v0x, omega]
    matrix x0_real(2, new double[2]{0.0, 0.0});

    // ---------------- Optymalizacja (funkcja kary + sym_NM) ----------------
    matrix x_curr_real = x0_real;
    double c_curr      = penaltyC0;

    // solution::clear_calls();

    for (int outer = 0; outer < penaltyMaxOuterIter &&
                        solution::f_calls < penaltyMaxCalls; ++outer)
    {
        // bieżące parametry kary
        matrix penalty_params_real(1, 1);
        penalty_params_real(0) = c_curr;

        // minimalizacja F(x) = f(x) + c * S(x)
        solution res = sym_NM(
            penalty_objective_function_lab4R,
            x_curr_real,          // punkt startowy dla tej iteracji
            1.0, 1.0, 0.5, 2.0,   // parametry N-M (alfa, beta, gamma, delta)
            0.5,                  // sigma
            penaltyEps / 10.0,    // dokładność wewnętrznego N-M
            penaltyMaxCalls / penaltyMaxOuterIter,  // limit wywołań na 1 iterację
            penalty_params_real,  // parametry funkcji kary
            NAN
        );

        // warunek stopu metody funkcji kary
        if (norm(res.x - x_curr_real) < penaltyEps) {
            x_curr_real = res.x;
            break;
        }

        x_curr_real = res.x;
        c_curr *= penaltyScale;
    }

        // ---------------- Wyniki optymalizacji ----------------
    cout << "\n--- Wyniki optymalizacji ---" << endl;
    cout << "Optymalne parametry:" << endl;
    cout << "v0x*   = " << x_curr_real(0) << " m/s" << endl;
    cout << "omega* = " << x_curr_real(1) << " rad/s" << endl;

    // rzeczywista wartość funkcji celu (bez kary, minimalizowaliśmy -x_end)
    // solution::clear_calls();
    matrix final_result = ff4R(x_curr_real, NAN, NAN);
    double x_end_opt    = -final_result(0);
    int f_calls_opt     = solution::f_calls;  // liczba wywołań funkcji celu dla optymalizacji

    cout << "x_end* = " << x_end_opt << " m" << endl;
    cout << "Liczba wywołań funkcji celu: " << f_calls_opt << endl;

    // =================== Symulacja trajektorii ===================

    // warunki początkowe: [x0, y0, vx0, vy0]
    matrix Y0_opt(4, 1);
    Y0_opt(0) = 0.0;               // x0
    Y0_opt(1) = 100.0;             // y0
    Y0_opt(2) = x_curr_real(0);    // vx0 = v0x*
    Y0_opt(3) = 0.0;               // vy0

    // parametry ruchu: [v0x, omega]
    matrix motion_params_opt(2, 1);
    motion_params_opt(0) = x_curr_real(0);
    motion_params_opt(1) = x_curr_real(1);

    // rozwiązywanie ODE
    matrix* Y_opt = solve_ode(
        ball_motion_l4,
        0.0,    // t0
        0.01,   // dt
        7.0,    // t_end
        Y0_opt,
        NAN,
        motion_params_opt
    );

    // =================== Obliczenie x* dla y = 50 m ===================
    double x_at_y50 = 0.0;
    int n_opt = get_len(Y_opt[0]);

    for (int i = 0; i < n_opt - 1; ++i) {
        double y1 = Y_opt[1](i, 1);
        double y2 = Y_opt[1](i + 1, 1);

        if ((y1 >= 50.0 && y2 <= 50.0) || (y1 <= 50.0 && y2 >= 50.0)) {
            double t1 = Y_opt[0](i);
            double t2 = Y_opt[0](i + 1);
            double x1 = Y_opt[1](i, 0);
            double x2 = Y_opt[1](i + 1, 0);

            // interpolacja czasu dla y = 50
            double t_interp = t1 + (50.0 - y1) / (y2 - y1) * (t2 - t1);
            // interpolacja x po czasie
            x_at_y50 = x1 + (t_interp - t1) / (t2 - t1) * (x2 - x1);
            break;
        }
    }

    // =================== CSV z wynikami optymalizacji ===================
    ofstream results_file("lab4_problem_rzeczywisty.csv");
    results_file << "v0x(0),omega(0),v0x*,omega*,xend*,x* dla y = 50m,Liczba wywołań funkcji celu\n";
    results_file << x0_real(0) << "," << x0_real(1) << ","
                << x_curr_real(0) << "," << x_curr_real(1) << ","
                << x_end_opt << "," << x_at_y50 << "," << f_calls_opt << "\n";
    results_file.close();

    // =================== CSV z trajektorią ===================
    ofstream trajectory_file("lab4_trajectory.csv");
    trajectory_file << "t,x,y\n";

    for (int i = 0; i < n_opt; ++i) {
        if (Y_opt[1](i, 1) >= 0.0) { // tylko nad ziemią
            trajectory_file << Y_opt[0](i)   << ","
                            << Y_opt[1](i,0) << ","  // x
                            << Y_opt[1](i,1) << "\n";  // y
        }
    }
    trajectory_file.close();

    delete[] Y_opt;
}

void lab5()
{
    cout << "=== LAB 5 - Metody Gradientowe ===" << endl;
    
    // --- CZĘŚĆ A: Funkcja Testowa ---
    cout << "\n--- Czesc A: Funkcja testowa ---" << endl;
    
    // Konfiguracja
    double steps[] = {0.05, 0.25, 0.0}; // 0.0 oznacza krok zmienny
    string step_names[] = {"0.05", "0.25", "Variable"};
    int Nmax = 1000;
    double epsilon = 1e-4;
    
    ofstream res_test("lab5_test_results.csv");
    res_test << "Method,Step,SuccessRate,AvgIter,AvgX1,AvgX2,AvgF\n";

    // Wykonaj 100 losowań dla każdego wariantu
    for (int s = 0; s < 3; ++s) {
        double h = steps[s];
        
        // Statystyki dla SD, CG, Newton
        double stats[3][4] = {0}; // [Method][Success, IterSum, x1Sum, x2Sum]
        int attempts = 100;

        srand(1234); // Stałe ziarno dla powtarzalności porównań między krokami

        for (int i = 0; i < attempts; ++i) {
            // Losowy punkt startowy [-2, 2]
            matrix x0(2, 1);
            x0(0) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
            x0(1) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
            
            solution sol;

            // 1. Steepest Descent (SD)
            solution::clear_calls();
            try {
                sol = SD(ff5T, gf5T, x0, h, epsilon, Nmax, NAN, NAN);
                // Sprawdź czy znalazł minimum globalne (około -1.0316)
                // Punkt globalny to ok. (1.74, -0.87) lub (-1.74, 0.87) f ~ 0.29
                // Instrukcja mówi o minimum globalnym. Funkcja camel ma ich kilka.
                // Przyjmijmy sukces jeśli flaga=1
                if (sol.flag == 1) {
                    stats[0][0]++;
                    stats[0][1] += solution::f_calls; // Lub liczba iteracji wewnątrz algorytmu
                    stats[0][2] += sol.x(0);
                    stats[0][3] += sol.x(1);
                }
            } catch(...) {}

            // 2. Conjugate Gradient (CG)
            solution::clear_calls();
            try {
                sol = CG(ff5T, gf5T, x0, h, epsilon, Nmax, NAN, NAN);
                if (sol.flag == 1) {
                    stats[1][0]++;
                    stats[1][1] += solution::f_calls;
                    stats[1][2] += sol.x(0);
                    stats[1][3] += sol.x(1);
                }
            } catch(...) {}

            // 3. Newton
            solution::clear_calls();
            try {
                sol = Newton(ff5T, gf5T, Hf5T, x0, h, epsilon, Nmax, NAN, NAN);
                if (sol.flag == 1) {
                    stats[2][0]++;
                    stats[2][1] += solution::f_calls;
                    stats[2][2] += sol.x(0);
                    stats[2][3] += sol.x(1);
                }
            } catch(...) {}
        }
        
        // Zapisz wyniki uśrednione
        string methods[] = {"SD", "CG", "Newton"};
        for(int m=0; m<3; ++m) {
            double success = stats[m][0];
            if(success > 0) {
                res_test << methods[m] << "," << step_names[s] << ","
                             << (success/attempts)*100 << "%,"
                             << stats[m][1]/success << ","
                             << stats[m][2]/success << ","
                             << stats[m][3]/success << ","
                             << ff5T(matrix(2, new double[2]{stats[m][2]/success, stats[m][3]/success}), NAN, NAN)(0)
                             << "\n";
            } else {
                res_test << methods[m] << "," << step_names[s] << ",0%,0,0,0,0\n";
            }
        }
    }
    res_test.close();
    cout << "Wyniki testow zapisano do lab5_test_results.csv" << endl;


    // --- CZĘŚĆ B: Problem Rzeczywisty (Regresja) ---
    cout << "\n--- Czesc B: Problem rzeczywisty ---" << endl;
    
    // Wczytanie danych
    // Uwaga: Zakładam strukturę plików zgodną z instrukcją: XData.txt (3x100), YData.txt (1x100)
    // Jeśli pliki nie istnieją, należy je utworzyć lub podstawić dane.
    // XData w instrukcji ma [1, x1, x2]^T
    
    try {
        matrix X = read_matrix_from_file("XData.txt", 3, 100);
        matrix Y = read_matrix_from_file("YData.txt", 1, 100);
        
        cout << "Dane wczytane pomyslnie." << endl;
        
        matrix theta0(3, 1, 0.0); // Punkt startowy [0,0,0]
        double steps_real[] = {0.01, 0.001, 0.0001};
        
        ofstream res_real("lab5_real_results.csv");
        res_real << "Step,Theta0,Theta1,Theta2,Cost,Calls,Accuracy\n";
        
        for (int i = 0; i < 3; ++i) {
            double h = steps_real[i];
            solution::clear_calls();
            
            // Używamy metody CG zgodnie z instrukcją
            solution sol = CG(ff5R, gf5R, theta0, h, 1e-4, 10000, X, Y);
            
            // Obliczamy dokładność klasyfikacji
            matrix h_theta = hypothesis(sol.x, X);
            int correct = 0;
            int m = 100;
            for(int j=0; j<m; ++j) {
                int pred = (h_theta(0, j) >= 0.5) ? 1 : 0;
                if(pred == (int)Y(0, j)) correct++;
            }
            double acc = (double)correct / m * 100.0;
            
            cout << "Krok: " << h << ", Koszt: " << sol.y(0) << ", Acc: " << acc << "%" << endl;
            
            res_real << h << "," 
                     << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
                     << sol.y(0) << "," << solution::f_calls << "," << acc << "\n";
        }
        res_real.close();
        
    } catch (string ex) {
        cout << "Problem z czescia B (brak plikow?): " << ex << endl;
        cout << "Upewnij sie, ze pliki XData.txt i YData.txt sa w katalogu roboczym." << endl;
    }
}

void lab6()
{

}