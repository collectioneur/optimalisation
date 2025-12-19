void lab4()
{
    cout << "=== LAB 4 - Optymalizacja funkcji wielu zmiennych metodami gradientowymi ===" << endl << endl;
    
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
    
    // Wartość minimum globalnego (przybliżona) - funkcja ma kilka minimów lokalnych
    // Sprawdzamy czy wartość funkcji jest bliska minimum globalnemu
    double f_global_min = 0.0; // Przybliżona wartość minimum globalnego
    
    srand(time(NULL));
    
    // Dla każdej metody
    for (int m = 0; m < 3; ++m) {
        // Dla każdej długości kroku
        for (int s = 0; s < 3; ++s) {
            double h0 = steps[s];
            
            // Statystyki dla wartości średnich (tylko dla minimum globalnego)
            double sum_x1 = 0.0, sum_x2 = 0.0, sum_f = 0.0, sum_iter = 0.0;
            int count_global = 0;
            
            // Wykonaj 100 optymalizacji
            for (int i = 0; i < num_tests; ++i) {
                // Losowy punkt startowy z przedziału [-2, 2] x [-2, 2]
                matrix x0(2, 1);
                x0(0) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
                x0(1) = -2.0 + 4.0 * ((double)rand() / RAND_MAX);
                
                solution sol;
                solution::clear_calls();
                
                try {
                    if (m == 0) { // SD
                        sol = SD(ff5T, gf5T, x0, h0, epsilon, Nmax, NAN, NAN);
                    } else if (m == 1) { // CG
                        sol = CG(ff5T, gf5T, x0, h0, epsilon, Nmax, NAN, NAN);
                    } else { // Newton
                        sol = Newton(ff5T, gf5T, Hf5T, x0, h0, epsilon, Nmax, NAN, NAN);
                    }
                    
                    // Zapisz do tabeli 1
                    table1 << methods[m] << "," << step_names[s] << "," << (i+1) << ","
                           << x0(0) << "," << x0(1) << ","
                           << sol.x(0) << "," << sol.x(1) << ","
                           << sol.y(0) << "," << solution::f_calls << ","
                           << sol.flag << "\n";
                    
                    // Sprawdź czy znaleziono minimum globalne (flaga = 1 oznacza sukces)
                    if (sol.flag == 1) {
                        // Uznajemy za minimum globalne jeśli flaga = 1
                        // (można dodać dodatkowe sprawdzenie wartości funkcji)
                        sum_x1 += sol.x(0);
                        sum_x2 += sol.x(1);
                        sum_f += sol.y(0);
                        sum_iter += solution::f_calls;
                        count_global++;
                    }
                    
                } catch (string ex) {
                    table1 << methods[m] << "," << step_names[s] << "," << (i+1) << ","
                           << x0(0) << "," << x0(1) << ","
                           << "ERROR,ERROR,ERROR,ERROR,ERROR\n";
                }
                
                if ((i+1) % 20 == 0) {
                    cout << "Metoda: " << methods[m] << ", Krok: " << step_names[s] 
                         << ", Ukończono: " << (i+1) << "/" << num_tests << endl;
                }
            }
            
            // Zapisz wartości średnie do tabeli 2
            if (count_global > 0) {
                table2 << methods[m] << "," << step_names[s] << ","
                       << sum_x1 / count_global << ","
                       << sum_x2 / count_global << ","
                       << sum_f / count_global << ","
                       << sum_iter / count_global << ","
                       << count_global << "\n";
            } else {
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
    try {
        matrix X = read_matrix_from_file("XData.txt", 3, 100);
        matrix Y = read_matrix_from_file("YData.txt", 1, 100);
        
        cout << "Dane wczytane pomyslnie." << endl;
        
        // Punkt startowy: θ^(0) = [0, 0, 0]
        matrix theta0(3, 1, 0.0);
        double steps_real[] = {0.01, 0.001, 0.0001};
        
        // Tabela 3: wyniki optymalizacji
        ofstream table3("lab4_tabela3.csv");
        table3 << "Krok,Theta0,Theta1,Theta2,Koszt,Iteracje,P_theta\n";
        
        for (int i = 0; i < 3; ++i) {
            double h = steps_real[i];
            solution::clear_calls();
            
            // Optymalizacja metodą CG (zgodnie z wymaganiami)
            solution sol = CG(ff5R, gf5R, theta0, h, epsilon, Nmax, X, Y);
            
            // Oblicz procent poprawnie zaklasyfikowanych przypadków P(θ*)
            matrix h_theta = hypothesis(sol.x, X);
            int m = 100;
            int correct = 0;
            for (int j = 0; j < m; ++j) {
                int pred = (h_theta(0, j) >= 0.5) ? 1 : 0;
                if (pred == (int)Y(0, j)) {
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
        // Wczytaj wyniki i znajdź najlepszy
        ifstream best_file("lab4_tabela3.csv");
        string line;
        getline(best_file, line); // nagłówek
        
        double best_P = 0.0;
        matrix best_theta(3, 1);
        
        while (getline(best_file, line)) {
            stringstream ss(line);
            string token;
            vector<string> tokens;
            
            while (getline(ss, token, ',')) {
                tokens.push_back(token);
            }
            
            if (tokens.size() >= 7) {
                double P = stod(tokens[6]);
                if (P > best_P) {
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
        
        // Granica klasyfikacji: h_θ(x) = 0.5 => θ₀ + θ₁x₁ + θ₂x₂ = 0
        // Jeśli θ₂ ≠ 0: x₂ = -(θ₀ + θ₁x₁) / θ₂
        for (int j = 0; j < m; ++j) {
            double x1 = X(1, j);
            double x2 = X(2, j);
            double y = Y(0, j);
            
            // Oblicz x2 dla granicy klasyfikacji dla tego x1
            double x2_boundary = 0.0;
            if (fabs(best_theta(2)) > 1e-10) {
                x2_boundary = -(best_theta(0) + best_theta(1) * x1) / best_theta(2);
            }
            
            plot_data << x1 << "," << x2 << "," << y << "," << x2_boundary << "\n";
        }
        
        plot_data.close();
        cout << "Zapisano dane do wykresu do lab4_wykres_dane.csv" << endl;
        
    } catch (string ex) {
        cout << "Błąd w części B: " << ex << endl;
        cout << "Upewnij się, że pliki XData.txt i YData.txt są w katalogu roboczym." << endl;
    }
    
    cout << "\n=== LAB 4 zakończone ===" << endl;
}

