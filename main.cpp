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
#include "opt_alg.h"
#include <sstream>
#include <vector>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>

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
        // double l_test = 0.5;   // 500 mm = 0.5 m
        // double d_test = 0.025; // 25 mm = 0.025 m
        // double P = 2000.0;
        // double E = 120e9;
        // double rho = 8920.0;
        // double m = (M_PI * pow(d_test, 2) * l_test * rho) / 4.0;
        // double u = (64.0 * P * pow(l_test, 3)) / (3.0 * E * M_PI * pow(d_test, 4));
        // double sigma = (32.0 * P * l_test) / (M_PI * pow(d_test, 3));

        // cout << "=== WERYFIKACJA WZOROW (DANE Z KONSPEKTU) ===" << endl;
        // cout << "Zadane l = " << l_test * 1000 << " mm" << endl;
        // cout << "Zadane d = " << d_test * 1000 << " mm" << endl;
        // cout << "---------------------------------------------" << endl;
        // cout << "Masa (oczekiwana: ~2.19 kg):       " << m << " kg" << endl;
        // cout << "Ugiecie (oczekiwana: ~36.22 mm):   " << u * 1000.0 << " mm" << endl;
        // cout << "Naprezenie (oczekiwana: ~651.9 MPa): " << sigma / 1e6 << " MPa" << endl;
        // cout << "---------------------------------------------" << endl;
        // lab5();
        lab6();
    }
    catch (string EX_INFO)
    {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl
             << endl;
    }
    return 0;
}
void lab0()
{
    // Funkcja testowa
    double epsilon = 1e-2;            // dok�adno��
    int Nmax = 10000;                 // maksymalna liczba wywo�a� funkcji celu
    matrix lb(2, 1, -5), ub(2, 1, 5), // dolne oraz g�rne ograniczenie
        a(2, 1);                      // dok�adne rozwi�zanie optymalne
    solution opt;                     // rozwi�zanie optymalne znalezione przez algorytm
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a); // wywo�anie procedury optymalizacji
    cout << opt << endl
         << endl;            // wypisanie wyniku
    solution::clear_calls(); // wyzerowanie licznik�w

    // Wahadlo
    Nmax = 1000;                                        // dok�adno��
    epsilon = 1e-2;                                     // maksymalna liczba wywo�a� funkcji celu
    lb = 0, ub = 5;                                     // dolne oraz g�rne ograniczenie
    double teta_opt = 1;                                // maksymalne wychylenie wahad�a
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt); // wywo�anie procedury optymalizacji
    cout << opt << endl
         << endl;            // wypisanie wyniku
    solution::clear_calls(); // wyzerowanie licznik�w

    // Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1),                            // Y0 zawiera warunki pocz�tkowe
        MT = matrix(2, new double[2]{m2d(opt.x), 0.5});  // MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT); // rozwi�zujemy r�wnanie r�niczkowe
    ofstream Sout("symulacja_lab0.csv");                 // definiujemy strumie� do pliku .csv
    Sout << hcat(Y[0], Y[1]);                            // zapisyjemy wyniki w pliku
    Sout.close();                                        // zamykamy strumie�
    Y[0].~matrix();                                      // usuwamy z pami�ci rozwi�zanie RR
    Y[1].~matrix();
}

void lab1()
{
    ofstream out("wyniki_lab1.csv");
    out << "i;alpha;x0;a;b;it_exp;x_fib;f_fib;it_fib;min_fib;x_lag;f_lag;it_lag;min_lag\n";

    double eps = 1e-4;
    int Nmax = 1000;
    double alphas[3] = {1.2, 1.5, 2.0};

    srand(time(NULL));

    for (int k = 0; k < 3; ++k)
    {
        double alpha = alphas[k];

        for (int i = 0; i < 100; ++i)
        {
            double x0 = -100 + (rand() % 200);

            // metoda ekspansji
            double *p = expansion(ff1, x0, 1.0, alpha, Nmax);
            extern int expansion_calls;
            int iter_exp = expansion_calls;

            // metoda Fibonacciego
            double *fib_res = fib(ff1, p[0], p[1], eps);
            extern int fib_calls;
            int iter_fib = fib_calls;
            double f_fib = ff1(*fib_res);

            // metoda Lagrange’a (nowa wersja double*)
            double *lag_res = lag(ff1, p[0], p[1], eps, 1e-5, Nmax);
            extern int lag_calls;
            int iter_lag = lag_calls;
            double f_lag = ff1(*lag_res);

            // klasyfikacja minimum
            string min_fib_type = (fabs(f_fib + 0.9211) < 1e-3) ? "globalne" : "lokalne";
            string min_lag_type = (fabs(f_lag + 0.9211) < 1e-3) ? "globalne" : "lokalne";

            // zapis do csv
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

void lab1_all()
{
    ofstream out("wyniki_lab1_all.csv");
    out << "Opis;alpha;a;b;x_fib;f_fib;it_fib;min_fib;x_lag;f_lag;it_lag;min_lag\n";

    double eps = 1e-4;
    int Nmax = 1000;
    double a = -100;
    double b = 100;
    double alphas[3] = {1.2, 1.5, 2.0};

    for (int k = 0; k < 3; ++k)
    {
        double alpha = alphas[k];
        cout << "=== ALPHA = " << alpha << " ===\n";

        // bez ekspansji
        double *fib_res_noexp = fib(ff1, a, b, eps);
        extern int fib_calls;
        int iter_fib_noexp = fib_calls;
        double f_fib_noexp = ff1(*fib_res_noexp);

        double *lag_res_noexp = lag(ff1, a, b, eps, 1e-5, Nmax);
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

        // z ekspansją
        double x0 = -100 + (rand() % 200);
        double *p = expansion(ff1, x0, 1.0, alpha, Nmax);
        double a_exp = p[0];
        double b_exp = p[1];

        double *fib_res_exp = fib(ff1, a_exp, b_exp, eps);
        extern int fib_calls;
        int iter_fib_exp = fib_calls;
        double f_fib_exp = ff1(*fib_res_exp);

        double *lag_res_exp = lag(ff1, a_exp, b_exp, eps, 1e-5, Nmax);
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

    double *interval = expansion(target_f_l2, x0, d, alpha, 100);
    cout << "Found interval: [" << interval[0] << ", " << interval[1] << "]" << endl;

    double *fibo = fib(target_f_l2, interval[0], interval[1], epsilon);
    cout << "Fibonacci: D_A* = " << *fibo << ", y* = " << target_f_l2(*fibo) << endl;

    double *lagr = lag(target_f_l2, interval[0], interval[1], epsilon, 1e-5, 1000);
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

    // rosenbrock testowa funkcja
    solution::clear_calls();
    solution result = Rosen(target_func_l3, x0, s0, alpha, beta, epsilon, Nmax);

    cout << "Wywołanie testowe dla metody Rosenbrocka" << endl;
    cout << "x = [" << result.x(0) << ", " << result.x(1) << "]" << endl;
    cout << "f(x) = " << result.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl
         << endl;

    // test HJ dla funkcji testowej
    solution::clear_calls();
    solution result_HJ_test = HJ(target_func_l3, x0, 1.0, 0.3, epsilon, Nmax);

    cout << "Wywołanie testowe dla metody Hooke-Jeeves" << endl;
    cout << "x = [" << result_HJ_test.x(0) << ", " << result_HJ_test.x(1) << "]" << endl;
    cout << "f(x) = " << result_HJ_test.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl
         << endl;

    // test porpawności dla funkcji Q(Q_real_l3)
    cout << "Test poprawności dla Q" << endl;

    matrix test_k(2, new double[2]{5.0, 5.0});
    matrix Q = Q_real_l3(test_k, NAN, test_k);

    cout << "Q(k1=5, k2=5) = " << Q(0) << endl;

    cout << "Oczekiwana wartość(ze sprawka): 775.229" << endl
         << endl;

    // rozwozanie problemu rzeczywistego -suzkanie najlepszego k1 i k2 - metoda rosenbrocka
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
    cout << "Calls: " << solution::f_calls << endl
         << endl;

    // rozwozanie problemu rzeczywistego - metoda HJ
    cout << "Optymalizacja problemu rzeczywistego metoda Hooke-Jeeves" << endl;
    solution::clear_calls();
    solution result_real_HJ = HJ(Q_real_l3, k0, 2.0, 0.3, epsilon, Nmax);

    cout << "k1 = " << result_real_HJ.x(0) << endl;
    cout << "k2 = " << result_real_HJ.x(1) << endl;
    cout << "Q(k1,k2) = " << result_real_HJ.y(0) << endl;
    cout << "Calls: " << solution::f_calls << endl
         << endl;

    // symulcja z czasem z optymalnymi parametrami - metoda Rosenbrocka
    cout << "Symulacja z optymalnymi parametrami (Rosenbrock)" << endl;

    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix k_opt_rosen(2, 1);
    k_opt_rosen(0) = result_real_rosen.x(0);
    k_opt_rosen(1) = result_real_rosen.x(1);

    matrix *Y_rosen = solve_ode(target_func_real_l3, 0, 0.1, 100, Y0, NAN, k_opt_rosen);

    // symulcja z czasem z optymalnymi parametrami - metoda HJ
    cout << "Symulacja z optymalnymi parametrami (Hooke-Jeeves)" << endl;

    matrix k_opt_HJ(2, 1);
    k_opt_HJ(0) = result_real_HJ.x(0);
    k_opt_HJ(1) = result_real_HJ.x(1);

    matrix *Y_HJ = solve_ode(target_func_real_l3, 0, 0.1, 100, Y0, NAN, k_opt_HJ);

    // zapis wyników do pliku csv z obiema metodami
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
    cout << "=== LAB 4 - Optymalizacja funkcji wielu zmiennych metodami gradientowymi ===" << endl
         << endl;

    // ========== CZĘŚĆ A: TESTOWA FUNKCJA CELU ==========
    cout << "--- Czesc A: Testowa funkcja celu ---" << endl;

    // Parametry
    double epsilon = 1e-4;
    int Nmax = 10000;
    double steps[] = {0.05, 0.25, 0.0}; // 0.0 oznacza krok zmienny (złoty podział)
    string step_names[] = {"0.05", "0.25", "Variable"};
    string methods[] = {"SD", "CG", "Newton"};
    int num_tests = 100;

    // Pliki wynikowe
    ofstream table1("lab4_tabela1.csv"); // Tabela 1: wszystkie wyniki
    table1 << "Metoda,Krok,Test,x1_start,x2_start,x1_opt,x2_opt,f_opt,Iteracje,Flaga\n";

    ofstream table2("lab4_tabela2.csv"); // Tabela 2: wartości średnie (tylko dla minimum globalnego)
    table2 << "Metoda,Krok,Srednia_x1,Srednia_x2,Srednia_f,Srednia_iteracji,Liczba_sukcesow\n";

    srand(time(NULL));

    // Dla każdej metody
    for (int m = 0; m < 3; ++m)
    {
        // Dla każdej długości kroku
        for (int s = 0; s < 3; ++s)
        {
            double h0 = steps[s];

            // Statystyki dla wartości średnich (tylko dla minimum globalnego)
            double sum_x1 = 0.0, sum_x2 = 0.0, sum_f = 0.0, sum_iter = 0.0;
            int count_global = 0;

            // Wykonaj 100 optymalizacji
            for (int i = 0; i < num_tests; ++i)
            {
                // Losowy punkt startowy z przedziału [-2, 2] x [-2, 2]
                matrix x0(2, 1);
                x0(0) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
                x0(1) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);

                solution sol;
                solution::clear_calls();

                try
                {
                    if (m == 0)
                    { // SD
                        sol = SD(ff5T, gf5T, x0, h0, epsilon, Nmax, NAN, NAN);
                    }
                    else if (m == 1)
                    { // CG
                        sol = CG(ff5T, gf5T, x0, h0, epsilon, Nmax, NAN, NAN);
                    }
                    else
                    { // Newton
                        sol = Newton(ff5T, gf5T, Hf5T, x0, h0, epsilon, Nmax, NAN, NAN);
                    }

                    // Zapisz do tabeli 1
                    table1 << methods[m] << "," << step_names[s] << "," << (i + 1) << ","
                           << x0(0) << "," << x0(1) << ","
                           << sol.x(0) << "," << sol.x(1) << ","
                           << sol.y(0) << "," << solution::f_calls << ","
                           << sol.flag << "\n";

                    // Sprawdź czy znaleziono minimum globalne (flaga = 1 oznacza sukces)
                    if (sol.flag == 1)
                    {
                        sum_x1 += sol.x(0);
                        sum_x2 += sol.x(1);
                        sum_f += sol.y(0);
                        sum_iter += solution::f_calls;
                        count_global++;
                    }
                }
                catch (string ex)
                {
                    table1 << methods[m] << "," << step_names[s] << "," << (i + 1) << ","
                           << x0(0) << "," << x0(1) << ","
                           << "ERROR,ERROR,ERROR,ERROR,ERROR\n";
                }

                if ((i + 1) % 20 == 0)
                {
                    cout << "Metoda: " << methods[m] << ", Krok: " << step_names[s]
                         << ", Ukończono: " << (i + 1) << "/" << num_tests << endl;
                }
            }

            // Zapisz wartości średnie do tabeli 2
            if (count_global > 0)
            {
                table2 << methods[m] << "," << step_names[s] << ","
                       << sum_x1 / count_global << ","
                       << sum_x2 / count_global << ","
                       << sum_f / count_global << ","
                       << sum_iter / count_global << ","
                       << count_global << "\n";
            }
            else
            {
                table2 << methods[m] << "," << step_names[s] << ",0,0,0,0,0\n";
            }
        }
    }

    table1.close();
    table2.close();
    cout << "Zapisano wyniki do lab4_tabela1.csv i lab4_tabela2.csv" << endl;

    // ========== CZĘŚĆ B: PROBLEM RZECZYWISTY (KLASYFIKACJA LOGISTYCZNA) ==========
    cout << "\n--- Czesc B: Problem rzeczywisty (klasyfikacja logistyczna) ---" << endl;

    // Wczytanie danych
    try
    {
        matrix X = read_matrix_from_file("XData.txt", 3, 100);
        matrix Y = read_matrix_from_file("YData.txt", 1, 100);

        cout << "Dane wczytane pomyslnie." << endl;

        // Punkt startowy: θ^(0) = [0, 0, 0]
        matrix theta0(3, 1, 0.0);
        double steps_real[] = {0.01, 0.001, 0.0001};

        // Tabela 3: wyniki optymalizacji
        ofstream table3("lab4_tabela3.csv");
        table3 << "Krok,Theta0,Theta1,Theta2,Koszt,Iteracje,P_theta\n";

        for (int i = 0; i < 3; ++i)
        {
            double h = steps_real[i];
            solution::clear_calls();

            // Optymalizacja metodą CG (zgodnie z wymaganiami)
            solution sol = CG(ff5R, gf5R, theta0, h, epsilon, Nmax, X, Y);

            // Oblicz procent poprawnie zaklasyfikowanych przypadków P(θ*)
            matrix h_theta = hypothesis(sol.x, X);
            int m = 100;
            int correct = 0;
            for (int j = 0; j < m; ++j)
            {
                int pred = (h_theta(0, j) >= 0.5) ? 1 : 0;
                if (pred == (int)Y(0, j))
                {
                    correct++;
                }
            }
            double P_theta = (double)correct / m * 100.0;

            // Zapisz do tabeli 3
            table3 << h << ","
                   << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
                   << sol.y(0) << "," << solution::f_calls << ","
                   << P_theta << "\n";

            cout << "Krok: " << h << ", Koszt: " << sol.y(0)
                 << ", P(θ*): " << P_theta << "%" << endl;
        }

        table3.close();
        cout << "Zapisano wyniki do lab4_tabela3.csv" << endl;

        // Znajdź najlepszy przypadek (najwyższe P(θ*))
        ifstream best_file("lab4_tabela3.csv");
        string line;
        getline(best_file, line); // nagłówek

        double best_P = 0.0;
        matrix best_theta(3, 1);

        while (getline(best_file, line))
        {
            stringstream ss(line);
            string token;
            vector<string> tokens;

            while (getline(ss, token, ','))
            {
                tokens.push_back(token);
            }

            if (tokens.size() >= 7)
            {
                double P = stod(tokens[6]);
                if (P > best_P)
                {
                    best_P = P;
                    best_theta(0) = stod(tokens[1]);
                    best_theta(1) = stod(tokens[2]);
                    best_theta(2) = stod(tokens[3]);
                }
            }
        }
        best_file.close();

        // Zapisz dane do wykresu (granica klasyfikacji)
        ofstream plot_data("lab4_wykres_dane.csv");
        plot_data << "x1,x2,y,decision_boundary\n";

        int m = 5;

        // Granica klasyfikacji: h_θ(x) = 0.5 => θ₀ + θ₁x₁ + θ₂x₂ = 0
        // Jeśli θ₂ ≠ 0: x₂ = -(θ₀ + θ₁x₁) / θ₂
        for (int j = 0; j < m; ++j)
        {
            double x1 = X(1, j);
            double x2 = X(2, j);
            double y = Y(0, j);

            // Oblicz x2 dla granicy klasyfikacji dla tego x1
            double x2_boundary = 0.0;
            if (fabs(best_theta(2)) > 1e-10)
            {
                x2_boundary = -(best_theta(0) + best_theta(1) * x1) / best_theta(2);
            }

            plot_data << x1 << "," << x2 << "," << y << "," << x2_boundary << "\n";
        }

        plot_data.close();
        cout << "Zapisano dane do wykresu do lab4_wykres_dane.csv" << endl;
    }
    catch (string ex)
    {
        cout << "Błąd w części B: " << ex << endl;
        cout << "Upewnij się, że pliki XData.txt i YData.txt są w katalogu roboczym." << endl;
    }

    cout << "\n=== LAB 4 zakończone ===" << endl;
}

// TO NIE LAB 5!!!!
// void lab5()
// {
//     cout << "=== LAB 5 - Metody Gradientowe ===" << endl;

//     // --- CZĘŚĆ A: Funkcja Testowa ---
//     cout << "\n--- Czesc A: Funkcja testowa ---" << endl;

//     // Konfiguracja
//     double steps[] = {0.05, 0.25, 0.0}; // 0.0 oznacza krok zmienny
//     string step_names[] = {"0.05", "0.25", "Variable"};
//     int Nmax = 1000;
//     double epsilon = 1e-4;

//     ofstream res_test("lab5_test_results.csv");
//     res_test << "Method,Step,SuccessRate,AvgIter,AvgX1,AvgX2,AvgF\n";

//     // Wykonaj 100 losowań dla każdego wariantu
//     for (int s = 0; s < 3; ++s)
//     {
//         double h = steps[s];

//         // Statystyki dla SD, CG, Newton
//         double stats[3][4] = {0}; // [Method][Success, IterSum, x1Sum, x2Sum]
//         int attempts = 100;       // Zmniejsz do 10 jesli test 'Variable' trwa za dlugo

//         srand(1234);

//         cout << "Rozpoczynanie testow dla kroku: " << step_names[s] << endl;

//         for (int i = 0; i < attempts; ++i)
//         {
//             // Logowanie postępu co 10 prób, żeby nie zalewać konsoli
//             if (i % 10 == 0)
//                 cout << "  Probka " << i << " / " << attempts << "..." << endl;

//             matrix x0(2, 1);
//             x0(0) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
//             x0(1) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);

//             // 1. Steepest Descent (SD)
//             solution::clear_calls();
//             try
//             {
//                 // TWORZYMY NOWĄ ZMIENNĄ sol DLA KAŻDEJ METODY (Zapobiega błędom)
//                 solution sol = SD(ff5T, gf5T, x0, h, epsilon, Nmax, NAN, NAN);
//                 if (sol.flag == 1)
//                 {
//                     stats[0][0]++;
//                     stats[0][1] += solution::f_calls + solution::g_calls + solution::H_calls;
//                     stats[0][2] += sol.x(0);
//                     stats[0][3] += sol.x(1);
//                 }
//             }
//             catch (...)
//             {
//             }

//             // 2. Conjugate Gradient (CG)
//             solution::clear_calls();
//             try
//             {
//                 solution sol = CG(ff5T, gf5T, x0, h, epsilon, Nmax, NAN, NAN);
//                 if (sol.flag == 1)
//                 {
//                     stats[1][0]++;
//                     stats[1][1] += solution::f_calls + solution::g_calls + solution::H_calls;
//                     stats[1][2] += sol.x(0);
//                     stats[1][3] += sol.x(1);
//                 }
//             }
//             catch (...)
//             {
//             }

//             // 3. Newton
//             solution::clear_calls();
//             try
//             {
//                 solution sol = Newton(ff5T, gf5T, Hf5T, x0, h, epsilon, Nmax, NAN, NAN);
//                 if (sol.flag == 1)
//                 {
//                     stats[2][0]++;
//                     stats[2][1] += solution::f_calls + solution::g_calls + solution::H_calls;
//                     stats[2][2] += sol.x(0);
//                     stats[2][3] += sol.x(1);
//                 }
//             }
//             catch (...)
//             {
//             }
//         }

//         // Zapisz wyniki uśrednione
//         string methods[] = {"SD", "CG", "Newton"};
//         for (int m = 0; m < 3; ++m)
//         {
//             double success = stats[m][0];
//             if (success > 0)
//             {
//                 res_test << methods[m] << "," << step_names[s] << ","
//                          << (success / attempts) * 100 << "%,"
//                          << stats[m][1] / success << ","
//                          << stats[m][2] / success << ","
//                          << stats[m][3] / success << ","
//                          << ff5T(matrix(2, new double[2]{stats[m][2] / success, stats[m][3] / success}), NAN, NAN)(0)
//                          << "\n";
//             }
//             else
//             {
//                 res_test << methods[m] << "," << step_names[s] << ",0%,0,0,0,0\n";
//             }
//         }
//         cout << "Zakonczono krok: " << step_names[s] << "\n----------------" << endl;
//     }
//     res_test.close();
//     cout << "Wyniki testow czesci A zapisano do lab5_test_results.csv" << endl;

//     // --- CZĘŚĆ B: Problem Rzeczywisty (Regresja) ---
//     cout << "\n--- Czesc B: Problem rzeczywisty ---" << endl;

//     try
//     {
//         // Upewnij się, że masz te pliki w folderze z projektem!
//         matrix X = read_matrix_from_file("XData.txt", 3, 100);
//         matrix Y = read_matrix_from_file("YData.txt", 1, 100);

//         cout << "Dane wczytane pomyslnie." << endl;

//         matrix theta0(3, 1, 0.0);
//         double steps_real[] = {0.01, 0.001, 0.0001};

//         ofstream res_real("lab5_real_results.csv");
//         res_real << "Step,Theta0,Theta1,Theta2,Cost,Calls,Accuracy\n";

//         for (int i = 0; i < 3; ++i)
//         {
//             double h = steps_real[i];
//             solution::clear_calls();

//             // Zwiekszamy limit iteracji dla problemu rzeczywistego
//             solution sol = CG(ff5R, gf5R, theta0, h, 1e-4, 20000, X, Y);

//             matrix h_theta = hypothesis(sol.x, X);
//             int correct = 0;
//             int m = 100;
//             for (int j = 0; j < m; ++j)
//             {
//                 int pred = (h_theta(0, j) >= 0.5) ? 1 : 0;
//                 if (pred == (int)Y(0, j))
//                     correct++;
//             }
//             double acc = (double)correct / m * 100.0;

//             cout << "Krok: " << h << ", Koszt: " << sol.y(0) << ", Acc: " << acc << "%" << endl;

//             res_real << h << ","
//                      << sol.x(0) << "," << sol.x(1) << "," << sol.x(2) << ","
//                      << sol.y(0) << "," << solution::f_calls << "," << acc << "\n";
//         }
//         res_real.close();
//         cout << "Wyniki testow czesci B zapisano do lab5_real_results.csv" << endl;
//     }
//     catch (string ex)
//     {
//         cout << "Problem z czescia B (brak plikow?): " << ex << endl;
//     }
// }

void lab5()
{
    cout << "=== LAB 5 - Funkcja testowa wielokryterialna ===" << endl;
    
    try
    {
        ofstream plik_test("lab5_test_wielokryterialna.csv");
        if (!plik_test.is_open())
            throw string("Nie udalo sie otworzyc pliku lab5_test_wielokryterialna.csv!");
        
        plik_test << "x1(0),x2(0),";
        plik_test << "x1* (a=1),x2* (a=1),f1* (a=1),f2* (a=1),Liczba wywołań (a=1),";
        plik_test << "x1* (a=10),x2* (a=10),f1* (a=10),f2* (a=10),Liczba wywołań (a=10),";
        plik_test << "x1* (a=100),x2* (a=100),f1* (a=100),f2* (a=100),Liczba wywołań (a=100)\n";

        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dist_x(-5.0, 5.0);
        
        double a_values[] = {1.0, 10.0, 100.0};
        
        cout << "Rozpoczynam optymalizacje funkcji testowej..." << endl;

        for (int i = 0; i <= 100; ++i)
        {
            double w = i / 100.0;

            double x1_0 = dist_x(gen);
            double x2_0 = dist_x(gen);
            
            plik_test << fixed << setprecision(6) << x1_0 << "," << x2_0;

            for (int j = 0; j < 3; ++j)
            {
                double a = a_values[j];
                
                matrix x_start(2, 1);
                x_start(0) = x1_0;
                x_start(1) = x2_0;
                
                matrix ud1(2, 1);
                ud1(0, 0) = w;
                ud1(1, 0) = a;
                matrix ud2; 

                solution::clear_calls();
                solution res = Powell(ff5T_multi, x_start, 1e-6, 10000, ud1, ud2);
                
                int total_calls = solution::f_calls;

                double x1_opt = res.x(0);
                double x2_opt = res.x(1);
                double f1_opt = pow(x1_opt, 2) + pow(x2_opt, 2);
                double f2_opt = pow(x1_opt - a, 2) + pow(x2_opt - a, 2);

                plik_test << "," << x1_opt << "," << x2_opt << "," 
                         << f1_opt << "," << f2_opt << "," << total_calls;
            }
            
            plik_test << "\n";
            
            if (i % 20 == 0)
            {
                cout << "Postep: w = " << w << " zakonczone" << endl;
            }
        }
        
        plik_test.close();
        cout << "Zakonczono funkcje testowa. Wyniki w 'lab5_test_wielokryterialna.csv'." << endl;
    }
    catch (string ex)
    {
        cerr << "Blad w czesci testowej: " << ex << endl;
    }
    cout << "\n=== LAB 5 - Problem rzeczywisty (belka) ===" << endl;
    
    const double P_force = 2000.0; 
    const double E_modulus = 120e9;  
    const double rho_density = 8920.0; 

    const double u_max = 2.5e-3;   
    const double sigma_max = 300e6;

    const double l_min = 0.2; 
    const double l_max = 1.0; 
    const double d_min = 0.01;
    const double d_max = 0.05;
    try
    {
        ofstream plik("wyniki_lab5.csv");
        if (!plik.is_open())
            throw string("Nie udalo sie otworzyc pliku do zapisu!");

        plik << "w;l_opt[mm];d_opt[mm];Masa[kg];Ugiecie[mm];Naprezenie[MPa];FunkcjaCelu\n";

        ofstream plik_rzecz("lab5_problem_rzeczywisty.csv");
        if (!plik_rzecz.is_open())
            throw string("Nie udalo sie otworzyc pliku lab5_problem_rzeczywisty.csv!");

        plik_rzecz << "l(0) [mm],d(0) [mm],l* [mm],d*[mm],masa* [kg],ugięcie* [mm],naprężenie* [Mpa],Liczba wywołań funkcji celu\n";

        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dist_l(l_min, l_max);
        uniform_real_distribution<double> dist_d(d_min, d_max);

        cout << "Rozpoczynam obliczenia dla 101 punktow..." << endl;

        matrix ud1(2, 1);
        matrix ud2;

        for (int i = 0; i <= 100; ++i)
        {
            double w = i / 100.0;
            ud1(0, 0) = w;

            matrix x_curr(2, 1);
            x_curr(0, 0) = dist_l(gen);
            x_curr(1, 0) = dist_d(gen);

            double l_start = x_curr(0, 0);
            double d_start = x_curr(1, 0);

            double c = 100.0;      
            double dc = 1.5; 
            int penalty_steps = 15;

            solution res;

            solution::clear_calls();

            for (int k = 0; k < penalty_steps; ++k)
            {
                ud1(1, 0) = c; 

                res = Powell(ff5R, x_curr, 1e-6, 2000, ud1, ud2);

                x_curr = res.x; 
                c *= dc;      
            }

            int total_calls = solution::f_calls;

            double l_opt = x_curr(0, 0);
            double d_opt = x_curr(1, 0);

            double mass = (M_PI * pow(d_opt, 2) * l_opt * rho_density) / 4.0;
            double defl = (64.0 * P_force * pow(l_opt, 3)) / (3.0 * E_modulus * M_PI * pow(d_opt, 4));
            double stress = (32.0 * P_force * l_opt) / (M_PI * pow(d_opt, 3));
            double f_val = w * mass + (1.0 - w) * defl;

            if (i % 10 == 0)
            {
                cout << "Postep: w = " << w << " -> Masa: " << mass << " kg, Ugiecie: " << defl * 1000 << " mm" << endl;
            }

            plik << fixed << setprecision(6)
                 << w << ";"
                 << l_opt * 1000.0 << ";"
                 << d_opt * 1000.0 << ";"
                 << mass << ";"
                 << defl * 1000.0 << ";" 
                 << stress / 1e6 << ";" 
                 << f_val << "\n";

            plik_rzecz << fixed << setprecision(6)
                       << l_start * 1000.0 << ","
                       << d_start * 1000.0 << ","
                       << l_opt * 1000.0 << ","
                       << d_opt * 1000.0 << ","
                       << mass << ","
                       << defl * 1000.0 << ","
                       << stress / 1e6 << ","
                       << total_calls << "\n";
        }

        plik.close();
        plik_rzecz.close();
        cout << "Zakonczono. Wyniki zapisano w 'wyniki_lab5.csv' i 'lab5_problem_rzeczywisty.csv'." << endl;
    }
    catch (string ex)
    {
        cerr << "Blad krytyczny: " << ex << endl;
    }
}

void lab6()
{
    cout << "=== LAB 6 - Algorytmy Ewolucyjne ===" << endl;
    srand(time(NULL));

    cout << "\n--- Czesc A: Testowa funkcja celu ---" << endl;
    
    int N = 2;
    matrix lb(N, 1, -5.0);
    matrix ub(N, 1, 5.0);
    double epsilon = 1e-6; 
    int Nmax = 10000;
    int pop_size_mi = 20;   
    int pop_size_lambda = 40; 

    double sigma_values[] = {0.01, 0.1, 1.0, 10.0, 100.0};
    int repetitions = 100;

    ofstream tab1("lab6_tabela1.csv");
    tab1 << "Początkowa wartość zakresu mutacji;Lp.;x1*;x2*;y*;Liczba wywołań funkcji celu;Minimum globalne [tak/nie]\n";

    ofstream tab2("lab6_tabela2.csv");
    tab2 << "Początkowa wartość zakresu mutacji;x1*;x2*;y*;Liczba wywołań funkcji celu;Liczba minimów globalnych\n";

    for (int s_idx = 0; s_idx < 5; ++s_idx)
    {
        double sigma_val = sigma_values[s_idx];
        matrix sigma0(N, 1, sigma_val);
        
        double sum_f = 0.0;
        double sum_calls = 0.0;
        double sum_x1 = 0.0;
        double sum_x2 = 0.0;
        int success_count = 0;

        cout << "Przetwarzanie sigma = " << sigma_val << "..." << endl;

        for (int i = 0; i < repetitions; ++i)
        {
            solution::clear_calls();
            
            solution res = EA(ff6T, N, lb, ub, pop_size_mi, pop_size_lambda, sigma0, epsilon, Nmax, NAN, NAN);
            
            string status_str = (res.flag == 1) ? "tak" : "nie";

            tab1 << sigma_val << ";" << (i+1) << ";"
                 << res.x(0) << ";" << res.x(1) << ";" 
                 << res.y(0) << ";" << solution::f_calls << ";"
                 << status_str << "\n";

            if (res.y(0) < 0.05) 
            {
                sum_f += res.y(0);
                sum_calls += solution::f_calls;
                sum_x1 += res.x(0);
                sum_x2 += res.x(1);
                success_count++;
            }
        }

        double avg_f = (success_count > 0) ? sum_f / success_count : 0.0;
        double avg_calls = (success_count > 0) ? sum_calls / success_count : 0.0;
        double avg_x1 = (success_count > 0) ? sum_x1 / success_count : 0.0;
        double avg_x2 = (success_count > 0) ? sum_x2 / success_count : 0.0;

        tab2 << sigma_val << ";" 
             << avg_x1 << ";" << avg_x2 << ";" 
             << avg_f << ";" << avg_calls << ";" 
             << success_count << "\n";
    }

    tab1.close();
    tab2.close();
    cout << "Wyniki czesci A zapisano." << endl;

    cout << "\n--- Czesc B: Problem rzeczywisty ---" << endl;

    ifstream file("lab6_experiment.txt");
    if (!file.is_open()) {
        cout << "Error: Nie mozna otworzyc pliku lab6_experiment.txt!" << endl;
        return;
    }

    vector<double> x1_vec, x2_vec;
    double val1, val2;
    char sep;

    while (file >> val1 >> sep >> val2 >> sep) {
        x1_vec.push_back(val1);
        x2_vec.push_back(val2);
    }
    file.close();

    int n_points = x1_vec.size();
    matrix ref_data(n_points, 2);
    for (int i = 0; i < n_points; ++i) {
        ref_data(i, 0) = x1_vec[i];
        ref_data(i, 1) = x2_vec[i];
    }
    cout << "Wczytano " << n_points << " punktow referencyjnych." << endl;

    matrix lb_real(2, 1, 0.1);
    matrix ub_real(2, 1, 3.0);
    matrix sigma0_real(2, 1, 0.5);

    solution::clear_calls();
    solution sol_real = EA(ff6R, 2, lb_real, ub_real, 20, 40, sigma0_real, 1e-4, 5000, ref_data, NAN);

    ofstream tab3("lab6_tabela3.csv");
    tab3 << "b1*;b2*;y*;Liczba wywołań funkcji celu\n";
    tab3 << sol_real.x(0) << ";" << sol_real.x(1) << ";" << sol_real.y(0) << ";" << solution::f_calls << "\n";
    tab3.close();

    matrix Y0(4, 1);
    matrix *Y_opt = solve_ode(df6, 0, 0.1, (n_points - 1) * 0.1, Y0, sol_real.x, NAN);
    
    ofstream sim_out("lab6_symulacja.csv");
    sim_out << "t;x1_ref;x2_ref;x1_opt;x2_opt\n";
    for (int i = 0; i < n_points; ++i) {
        sim_out << i * 0.1 << ";" << ref_data(i, 0) << ";" << ref_data(i, 1) << ";"
                << Y_opt[1](i, 0) << ";" << Y_opt[1](i, 2) << "\n";
    }
    sim_out.close();
    delete[] Y_opt;
}