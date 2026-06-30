#include <bits/stdc++.h>
using namespace std;
#include <fstream>
#include <string>

//Budowa macierzy w formie trzech wektorów 
void build_system(int N, vector<double>& a, vector<double>& d, vector<double>& c, vector<double>& bvec) {
    a.assign(N, 0.0); //poddiagonala
    d.assign(N, 0.0); //diagonala
    c.assign(N, 0.0); //naddiagonala
    bvec.assign(N, 0.0); //wektor b
    double h = 2.0 / (N - 1); //wartości z zadania
    double h2 = h * h;
    d[0] = 1.0; bvec[0] = 1.0; //boki macierzy A i wetkora b
    d[N-1] = 1.0; bvec[N-1] = 1.0;

    for (int i = 1; i < N-1; ++i) {
        a[i] = 1.0 / h2;
        d[i] = -2.0 / h2; //Wypełnienie macierzy wartościami
        c[i] = 1.0 / h2;
        bvec[i] = 0.0;
    }
}

// Eliminacja gaussa razem z budową macierzy A.
// Samo zbudowanie macierzy ma złożoność O(N^2) a 
// Eliminacja ma złożoność O(N^3)
vector<double> solve_dense_gauss(const vector<double>& a_tr, const vector<double>& d_tr, const vector<double>& c_tr, const vector<double>& bvec) {
    int N = (int)bvec.size();
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    // Budowa macierzy A na podstawie trzech wektorów (a_tr,d_tr,c_tr)
    for (int i = 0; i < N; ++i) {
        A[i][i] = d_tr[i];
        if (i > 0) A[i][i-1] = a_tr[i];
        if (i < N-1) A[i][i+1] = c_tr[i];
    }
    // Kopia wektora b, tak aby operacje eliminacji nie modyfikowały oryginalnego bvec
    vector<double> b = bvec; //Kopia wektora b O(N)

    // Główna pętla eliminacji Gaussa z częściowym pivotowaniem
    for (int i = 0; i < N; ++i) {
        // Znajdź wiersz z największym elementem w kolumnie i (od i do N-1)
        int piv = i;
        double maxv = fabs(A[i][i]);
        for (int r = i+1; r < N; ++r) {
            double v = fabs(A[r][i]);
            if (v > maxv) { maxv = v; piv = r; }
        }
        // Jeśli największy element nie jest na przekątnej, zamień wiersze
        if (piv != i) {
            swap(A[piv], A[i]);    // zamiana całych wierszy w macierzy A
            swap(b[piv], b[i]);    // i odpowiadających elementów wektora b
        }
        // Częściowe pivotowanie poprawia stabilność numeryczną (unika dzielenia przez małe liczby)
        double diag = A[i][i];
        if (fabs(diag) < 1e-18) throw runtime_error("Zero pivot in Gaussian elimination");

        // Eliminacja: dla każdego wiersza poniżej i anulujemy element w kolumnie i
        for (int r = i+1; r < N; ++r) {
            double fac = A[r][i] / diag;            // mnożnik 
            for (int col = i; col < N; ++col)       // aktualizujemy pozostałą część wiersza
                A[r][col] -= fac * A[i][col];
            b[r] -= fac * b[i];                     // i aktualizujemy odpowiadający element b
        }
    }
    // backsubstitution
    vector<double> x(N, 0.0);
    // Po zakończeniu eliminacji macierz A jest górnotrójkątna; rozwiązujemy U x = b
    for (int i = N-1; i >= 0; --i) {
        double s = b[i];
        // Odejmujemy wkład znanych już zmiennych x[j] dla j>i
        for (int j = i+1; j < N; ++j) s -= A[i][j] * x[j];
        // Dzielimy przez element diagonalny U[i][i] aby uzyskać x[i]
        x[i] = s / A[i][i];
    }
    return x;
}

//Algorytm thomasa dla dobrze uwarunkowanych macierzy trójdiagonlanych
//daje złożoność czasową O(N)
vector<double> solve_thomas(const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& bvec) {
    int N = (int)bvec.size();
    vector<double> cp(N, 0.0), dp(N, 0.0);
    vector<double> x(N, 0.0);
    // Inicjalizacja zmodyfikowanych współczynników c' i d' (cp i dp)
    // cp[0] = c[0] / d[0], dp[0] = b[0] / d[0] (jeśli d[0] != 0)
    cp[0] = (d[0] == 0.0) ? 0.0 : (c[0] / d[0]);  // c' dla pierwszego wiersza
    dp[0] = (d[0] == 0.0) ? 0.0 : (bvec[0] / d[0]); // d' dla pierwszego wiersza

    // forward substitution - przekształcamy układ do postaci górnotrójkątnej zapisanej przez cp i dp
    for (int i = 1; i < N; ++i) {
        // element diagonalny zmodyfikowany przez poprzednie kroki
        double den = d[i] - a[i] * cp[i-1];
        // Jeśli den jest bliskie zeru, algorytm zawodzi
        if (fabs(den) < 1e-18) throw runtime_error("Zero pivot in Thomas algorithm");

        // cp[i] to nowy współczynnik naddiagonalny dla wiersza i
        // Dla ostatniego wiersza nie ma naddiagonalnego elementu, więc cp[N-1]=0
        cp[i] = (i == N-1) ? 0.0 : (c[i] / den);

        // dp[i] to zmodyfikowany element prawej strony po uwzględnieniu a[i]*dp[i-1]
        dp[i] = (bvec[i] - a[i] * dp[i-1]) / den;
    }

    // backsubstitution - odzyskujemy x od końca (x[N-1] = dp[N-1], potem x[i] = dp[i] - cp[i]*x[i+1])
    x[N-1] = dp[N-1];
    for (int i = N-2; i >= 0; --i) x[i] = dp[i] - cp[i] * x[i+1];
    return x;
}

double l2norm(const vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x*x;
    return sqrt(s);
}
// Metoda iteracyjna Gauss-Seidel omega = 1
// Funkcja iteracyjnie przybliża rozwiązanie układu trójdiagonalnego A x = b
vector<double> solve_sor_tridiag(const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& bvec, double tol=1e-12, int maxit=100000, double omega=1.0) {
    int N = (int)bvec.size();               // rozmiar układu
    vector<double> x(N, 0.0);               // wektor przybliżony, inicjalizowany zerami

    // Główna pętla iteracyjna: wykonujemy kolejne iteracje aż do osiągnięcia
    // tolerancji (maxdiff < tol) lub do przekroczenia maxit.
    for (int iter = 0; iter < maxit; ++iter) {
        double maxdiff = 0.0;                // maksymalna zmiana elementu w tej iteracji

        // Przechodzimy po wszystkich komponentach x 
        for (int i = 0; i < N; ++i) {
            double sigma = 0.0;              // suma wkładów A[i][j]*x[j] dla j != i

            // jeśli istnieje element podprzekątnej, dodajemy jego wkład
            if (i > 0) sigma += a[i] * x[i-1];

            // jeśli istnieje element nadprzekątnej, dodajemy jego wkład
            // tutaj x[i+1] to jeszcze stara wartość z tej samej iteracji
            if (i < N-1) sigma += c[i] * x[i+1];

            // rozwiązujemy równanie i-tego wiersza względem x_i przy założeniu
            // że pozostałe składowe są znane: x_new = (b_i - sigma) / A[i][i]
            double x_new = (bvec[i] - sigma) / d[i];

            // przechowujemy starą wartość, by móc policzyć zmianę
            double x_old = x[i];

            x[i] = (1.0 - omega) * x_old + omega * x_new;

            // aktualizujemy maksymalną zmianę w wektorze (do kryterium stopu)
            maxdiff = max(maxdiff, fabs(x[i] - x_old));
        }

        // Jeśli największa zmiana w tej iteracji jest mniejsza niż tolerancja,
        // uznajemy, że iteracje zbiegnęły i przerywamy pętlę.
        if (maxdiff < tol) break;
    }

    return x; // zwracamy przybliżone rozwiązanie
}

int main(int argc, char** argv) {
        // Tryby uruchomienia programu:
        // 1) Normalne: ./N01.exe N 
        // 2) Tryb wsadowy: ./N01.exe --batch [points] [minN] [maxN] [out.csv] np .\N01.exe --batch 25 10 1000 times_N01.csv
        if (argc >= 2 && string(argv[1]) == "--batch") {
            int points = 25;
            long long minN = 10;
            long long maxN = 100000;
            string out = "times_N01.csv";
            if (argc >= 3) points = atoi(argv[2]);
            if (argc >= 4) minN = atoll(argv[3]);
            if (argc >= 5) maxN = atoll(argv[4]);
            if (argc >= 6) out = argv[5];

            ofstream ofs(out);
            if (!ofs) {
                cerr << "Failed to open output file: " << out << "\n";
                return 1;
            }
            ofs << "N,dense_time,thomas_time,sor_time\n";

            // Generujemy listę rozmiarów N rozłożonych logarytmicznie
            // między minN a maxN. Następnie
            // dla każdej unikalnej wartości N (>1) budujemy układ równań
            // i mierzymy czas trzech metod: eliminacji Gaussa (gęsta macierz),
            // algorytmu Thomasa (trójdiagonalny) oraz iteracyjnego SOR.
            // Wyniki (N oraz czasy) zapisywane są do pliku CSV, a postęp

            vector<long long> Ns;
            for (int i = 0; i < points; ++i) {
                double t = (double)i / (points - 1);
                double val = exp(log((double)minN) * (1 - t) + log((double)maxN) * t);
                long long ni = max( (long long)2, (long long) (val + 0.5) );
                if (Ns.empty() || ni != Ns.back()) Ns.push_back(ni);
            }

            for (long long Ni : Ns) {
                int Ni_int = (int)Ni;
                if (Ni_int < 2) continue;
                vector<double> a,d,c,bvec;
                build_system(Ni_int, a,d,c,bvec);

                auto t0 = chrono::high_resolution_clock::now();
                vector<double> x_dense = solve_dense_gauss(a,d,c,bvec);
                auto t1 = chrono::high_resolution_clock::now();
                double t_dense = chrono::duration<double>(t1-t0).count();

                t0 = chrono::high_resolution_clock::now();
                vector<double> x_thomas = solve_thomas(a,d,c,bvec);
                t1 = chrono::high_resolution_clock::now();
                double t_thomas = chrono::duration<double>(t1-t0).count();

                t0 = chrono::high_resolution_clock::now();
                vector<double> x_sor = solve_sor_tridiag(a,d,c,bvec,1e-12, 100000, 1.0);
                t1 = chrono::high_resolution_clock::now();
                double t_sor = chrono::duration<double>(t1-t0).count();

                ofs << Ni_int << "," << setprecision(12) << t_dense << "," << t_thomas << "," << t_sor << "\n";
                cerr << "Done N=" << Ni_int << "\n";
            }
            ofs.close();
            cerr << "Batch complete, written to: " << out << "\n";
            return 0;
        }

        int N = 1000;
        if (argc >= 2) N = atoi(argv[1]);
        if (N < 2) { cerr << "N must be >= 2\n"; return 1; }
        // Można zapisać pełne rozwiązanie do CSV przy wywołaniu: <program> N --save-solution filename
        bool save_solution = false;
        string solution_fname;
        if (argc >= 3) {
            string arg2 = argv[2];
            if (arg2 == string("--save-solution") && argc >= 4) {
                save_solution = true;
                solution_fname = argv[3];
            }
        }

    vector<double> a, d, c, bvec;
    build_system(N, a, d, c, bvec);
    
    // Uruchamianie eliminacji gaussa i pomiar czasu
    auto t0 = chrono::high_resolution_clock::now();
    vector<double> x_dense = solve_dense_gauss(a,d,c,bvec);
    auto t1 = chrono::high_resolution_clock::now();
    double t_dense = chrono::duration<double>(t1-t0).count();

    // Uruchamianie Thomasa i pomiar czasu
    t0 = chrono::high_resolution_clock::now();
    vector<double> x_thomas = solve_thomas(a,d,c,bvec);
    t1 = chrono::high_resolution_clock::now();
    double t_thomas = chrono::duration<double>(t1-t0).count();

    // Uruchamianie Gauss-Seidel i pomiar czasu
    t0 = chrono::high_resolution_clock::now();
    vector<double> x_sor = solve_sor_tridiag(a,d,c,bvec,1e-12, 100000, 1.0);
    t1 = chrono::high_resolution_clock::now();
    double t_sor = chrono::duration<double>(t1-t0).count();

    // Compare solutions
    vector<double> diff_dt(N), diff_ds(N);
    for (int i = 0; i < N; ++i) {
        diff_dt[i] = x_dense[i] - x_thomas[i];
        diff_ds[i] = x_dense[i] - x_sor[i];
    }

    cout << fixed << setprecision(8);
    cout << "N = " << N << "\n";
    cout << "Dense time: " << t_dense << " s\n";
    cout << "Thomas time: " << t_thomas << " s\n";
    cout << "Gauss-Seidel time: " << t_sor << " s\n";
    cout << "L2(diff dense-thomas): " << l2norm(diff_dt) << "\n"; //Różnica między wartościami dla eliminacji gaussa i algorytmu thomasa
    cout << "L2(diff gauss-sor):     " << l2norm(diff_ds) << "\n"; //Różnica między wartościami dla eliminacji gaussa i metody gauss-seidel
    cout << "Sample solution (first 10 entries):\n";
    int M = min(N, 10);
    for (int i = 0; i < M; ++i) cout << "i=" << i << "  " << x_dense[i] << "  " << x_thomas[i] << "  " << x_sor[i] << "\n";

    cout << "CSV_TIMES," << N << "," << t_dense << "," << t_thomas << "," << t_sor << "\n";

    // Automatyczne zapisywanie dla N >= 1000 
    bool will_save = save_solution || (N >= 1000 && !save_solution);
    string auto_fname;
    if (save_solution) auto_fname = solution_fname;
    else auto_fname = string("solution_N") + to_string(N) + string(".csv");

    if (will_save) {
        ofstream ofs(auto_fname);
        if (!ofs) {
            cerr << "Failed to open solution file for writing: " << auto_fname << "\n";
        } else {
            ofs << "x,u\n";
            double h = 2.0 / (N - 1);
            for (int i = 0; i < N; ++i) {
                double x = i * h;
                ofs << setprecision(8) << x << "," << x_dense[i] << "\n";
            }
            ofs.close();
            cout << "Saved solution to: " << auto_fname << "\n";
        }
    }

    return 0;
}
