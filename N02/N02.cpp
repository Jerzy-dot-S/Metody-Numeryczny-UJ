#include <bits/stdc++.h>
using namespace std;
// Przenośna definicja PI (M_PI nie jest standardowo gwarantowane)
const double PI = acos(-1.0);

void build_system(int N, vector<double>& a, vector<double>& d, vector<double>& c, vector<double>& b, double& corner_val) {
    // Budujemy trójwektorową reprezentację macierzy T (pod-diagonal, diagonal, nad-diagonal)
    // oraz wektor prawej strony b. a[0] oraz c[N-1] pozostają nieużywane.
    a.assign(N, 0.0);
    d.assign(N, 0.0);
    c.assign(N, 0.0);
    b.assign(N, 0.0);
    double h = 2.0 / (N - 1);
    double h2 = h * h;
    double off = 1.0 / h2;
    double diag = -2.0 / h2;
    corner_val = off; // wartości elementów rogowych macierzy A

    for (int i = 0; i < N; ++i) {
        d[i] = diag;
        if (i > 0) a[i] = off;
        if (i < N-1) c[i] = off;
        b[i] = cos(4.0 * PI * (i) / (double)N); // indeksowanie: i odpowiada (n-1)
    }
}

// Eliminacja gaussa razem z budową macierzy A.
// Samo zbudowanie macierzy ma złożoność O(N^2) a 
// Eliminacja ma złożoność O(N^3)
vector<double> solve_dense_gauss(const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& b, double corner) {
    int N = (int)b.size();
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        A[i][i] = d[i];
        if (i > 0) A[i][i-1] = a[i];
        if (i < N-1) A[i][i+1] = c[i];
    }
    // elementy rogowe macierzy (dodane osobno)
    A[0][N-1] = corner;
    A[N-1][0] = corner;

    vector<double> rhs = b;

    // Główna pętla eliminacji Gaussa z częściowym pivotowaniem
    for (int i = 0; i < N; ++i) {
        int piv = i;
        double maxv = fabs(A[i][i]);
        for (int r = i+1; r < N; ++r) {
            double v = fabs(A[r][i]);
            if (v > maxv) { maxv = v; piv = r; }
        }
        if (piv != i) {
            swap(A[piv], A[i]);
            swap(rhs[piv], rhs[i]);
        }
        double diag = A[i][i];
        if (fabs(diag) < 1e-18) throw runtime_error("Zero pivot in dense Gaussian");
        for (int r = i+1; r < N; ++r) {
            double fac = A[r][i] / diag;
            for (int col = i; col < N; ++col) A[r][col] -= fac * A[i][col];
            rhs[r] -= fac * rhs[i];
        }
    }
    vector<double> x(N,0.0);
    for (int i = N-1; i >= 0; --i) {
        double s = rhs[i];
        for (int j = i+1; j < N; ++j) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return x;
}

// Algorytm Thomasa dla macierzy trójdiagonalnej T: rozwiązuje T x = b w czasie O(N).
// Dla działania shermana morissona
vector<double> solve_thomas(const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& b) {
    int N = (int)b.size();
    vector<double> cp(N, 0.0), dp(N, 0.0), x(N, 0.0);
    if (fabs(d[0]) < 1e-18) throw runtime_error("Zero pivot in Thomas");
    cp[0] = (N==1) ? 0.0 : c[0] / d[0];
    dp[0] = b[0] / d[0];
    for (int i = 1; i < N; ++i) {
        double den = d[i] - a[i] * cp[i-1];
        if (fabs(den) < 1e-18) throw runtime_error("Zero pivot in Thomas forward");
        cp[i] = (i==N-1) ? 0.0 : c[i] / den;
        dp[i] = (b[i] - a[i] * dp[i-1]) / den;
    }
    x[N-1] = dp[N-1];
    for (int i = N-2; i >= 0; --i) x[i] = dp[i] - cp[i] * x[i+1];
    return x;
}

// Rozwiązanie macierzy z rogami stosując dwa razy wzór shermana morissona
vector<double> sherman_morisson(const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& b, double corner) {
    int N = (int)b.size();
    // Krok 0: oblicz T^{-1} * b, T^{-1} * e1, T^{-1} * eN
    vector<double> y = solve_thomas(a,d,c,b); // T^{-1} b
    vector<double> e1(N,0.0), eN(N,0.0);
    e1[0] = 1.0; eN[N-1] = 1.0;
    vector<double> Tinve1 = solve_thomas(a,d,c,e1); // T^{-1} e1
    vector<double> TinveN = solve_thomas(a,d,c,eN); // T^{-1} eN

    // Pierwszy krok Sherman–Morrison: B = T + corner * (e1 eN^T)
    // B^{-1} = T^{-1} - (corner * T^{-1} e1 * (eN^T T^{-1})) / (1 + corner * eN^T T^{-1} e1)
    double alpha1 = Tinve1[N-1]; // eN^T T^{-1} e1
    double beta1 = y[N-1];      // eN^T T^{-1} b
    double denom1 = 1.0 + corner * alpha1;
    if (fabs(denom1) < 1e-18) throw runtime_error("Singular in Sherman-Morrison step 1");
    vector<double> x1(N,0.0);
    double fac1 = (corner * beta1) / denom1;
    for (int i = 0; i < N; ++i) x1[i] = y[i] - fac1 * Tinve1[i]; // B^{-1} b

    // Potrzebujemy teraz B^{-1} * eN (do drugiego kroku)
    // B^{-1} eN = T^{-1} eN - (corner * T^{-1} e1) * (eN^T T^{-1} eN) / denom1
    double gamma = TinveN[N-1]; // eN^T T^{-1} eN
    vector<double> Binv_eN(N,0.0);
    double coeff = (corner * gamma) / denom1;
    for (int i = 0; i < N; ++i) Binv_eN[i] = TinveN[i] - coeff * Tinve1[i];

    // Drugi krok Sherman–Morrison: A = B + corner * (eN e1^T)
    // A^{-1} = B^{-1} - (corner * B^{-1} eN * (e1^T B^{-1})) / (1 + corner * e1^T B^{-1} eN)
    double beta2 = x1[0]; // e1^T B^{-1} b
    double alpha2 = Binv_eN[0]; // e1^T B^{-1} eN
    double denom2 = 1.0 + corner * alpha2;
    if (fabs(denom2) < 1e-18) throw runtime_error("Singular in Sherman-Morrison step 2");
    double fac2 = (corner * beta2) / denom2;
    vector<double> x(N,0.0);
    for (int i = 0; i < N; ++i) x[i] = x1[i] - fac2 * Binv_eN[i];
    return x;
}

// Metoda iteracyjna Gauss-Seidel omega = 1
vector<double> solve_sor_full(const vector<double>& a, const vector<double>& d, const vector<double>& c, const vector<double>& b, double corner, double tol=1e-12, int maxit=100000, double omega=1.0) {
    int N = (int)b.size();
    vector<double> x(N,0.0);
    for (int iter = 0; iter < maxit; ++iter) {
        double maxdiff = 0.0;
        for (int i = 0; i < N; ++i) {
            double sigma = 0.0;
            if (i > 0) sigma += a[i] * x[i-1];
            if (i < N-1) sigma += c[i] * x[i+1];
            // wkład elementów rogowych do sumy sigma
            if (i == 0) sigma += corner * x[N-1];
            if (i == N-1) sigma += corner * x[0];
            double x_new = (b[i] - sigma) / d[i];
            double x_old = x[i];
            x[i] = (1.0 - omega) * x_old + omega * x_new;
            maxdiff = max(maxdiff, fabs(x[i] - x_old));
        }
        if (maxdiff < tol) break;
    }
    return x;
}

int main(int argc, char** argv) {
    // Główna funkcja programu.
    // Program obsługuje dwa tryby pracy:
    // 1) Tryb wsadowy: --batch [points] [minN] [maxN] [out.csv] — generuje listę
    //    rozmiarów N i mierzy czasy trzech metod zapisując
    //    wyniki do pliku CSV.
    // 2) Tryb pojedynczy: .\N02.exe [matrix_size] — mierzy czasy dla
    //    pojedynczego N, wypisuje krótkie podsumowanie i zapisuje
    //    pełne rozwiązanie do pliku CSV od N > 1000
    if (argc >= 2 && string(argv[1]) == "--batch") {
        int points = 25;
        long long minN = 10;
        long long maxN = 100000;
        string out = "times_N02.csv";
        if (argc >= 3) points = atoi(argv[2]);
        if (argc >= 4) minN = atoll(argv[3]);
        if (argc >= 5) maxN = atoll(argv[4]);
        if (argc >= 6) out = argv[5];

        ofstream ofs(out);
        if (!ofs) {
            cerr << "Failed to open output file: " << out << "\n";
            return 1;
        }
        ofs << "N,dense_time,sherman_morrison_time,sor_time\n";

            // wygeneruj liczby całkowite rozłożone logarytmicznie między minN a maxN
        vector<long long> Ns;
        for (int i = 0; i < points; ++i) {
            double t = (double)i / (points - 1);
            double val = exp(log((double)minN) * (1 - t) + log((double)maxN) * t);
            long long ni = max( (long long)2, (long long) (val + 0.5) );
            if (Ns.empty() || ni != Ns.back()) Ns.push_back(ni);
        }

        // Dla każdej wartości N w liście: budujemy układ, uruchamiamy trzy
        // metody i mierzymy czas ich wykonania. Wyniki zapisujemy do CSV.
        for (long long Ni : Ns) {
            int Ni_int = (int)Ni;
            if (Ni_int < 2) continue;
            vector<double> a,d,c,b;
            double corner;
            build_system(Ni_int, a,d,c,b, corner);

            using clk = chrono::high_resolution_clock;
            double t_dense = -1.0, t_sherman = -1.0, t_sor = -1.0;
            // 1) Eliminacja Gaussa na pełnej macierzy — O(N^3)
            try {
                auto t0 = clk::now();
                vector<double> x_dense = solve_dense_gauss(a,d,c,b,corner);
                auto t1 = clk::now();
                t_dense = chrono::duration<double>(t1-t0).count();
            } catch (const exception &e) {
                cerr << "Dense solve failed for N=" << Ni_int << ": " << e.what() << "\n";
            }
            // 2) Sherman–Morrison koszt O(N).
            try {
                auto t0 = clk::now();
                vector<double> x_sherman = sherman_morisson(a,d,c,b,corner);
                auto t1 = clk::now();
                t_sherman = chrono::duration<double>(t1-t0).count();
            } catch (const exception &e) {
                cerr << "Sherman-Morrison solve failed for N=" << Ni_int << ": " << e.what() << "\n";
            }
            // 3) Metoda iteracyjna SOR Może być wolniejsza koło O(N^2)
            try {
                auto t0 = clk::now();
                vector<double> x_sor = solve_sor_full(a,d,c,b,corner,1e-12,100000,1.0);
                auto t1 = clk::now();
                t_sor = chrono::duration<double>(t1-t0).count();
            } catch (const exception &e) {
                cerr << "SOR solve failed for N=" << Ni_int << ": " << e.what() << "\n";
            }

            ofs << Ni_int << "," << setprecision(12) << t_dense << "," << t_sherman << "," << t_sor << "\n";
            cerr << "Done N=" << Ni_int << "\n";
        }
        ofs.close();
        cerr << "Batch complete, written to: " << out << "\n";
        return 0;
    }

    // Tryb pojedynczy: jeśli podano N jako argument, użyj go, domyślnie N=1000.
    int N = 1000;
    if (argc >= 2) N = atoi(argv[1]);
    if (N < 2) { cerr << "N must be >= 2\n"; return 1; }

    vector<double> a,d,c,b;
    double corner;
    build_system(N, a,d,c,b, corner);

    using clk = chrono::high_resolution_clock;

    double t_dense = -1.0, t_sherman = -1.0, t_sor = -1.0;
    vector<double> x_dense, x_sherman, x_sor;
    // Pojedyncze pomiary: Eliminacja gaussa
    try {
        auto t0 = clk::now();
        x_dense = solve_dense_gauss(a,d,c,b,corner);
        auto t1 = clk::now();
        t_dense = chrono::duration<double>(t1-t0).count();
    } catch (const exception &e) {
        cerr << "Dense solve failed for N=" << N << ": " << e.what() << "\n";
    }
    // Pojedyncze pomiary: Sherman-Morrison
    try {
        auto t0 = clk::now();
        x_sherman = sherman_morisson(a,d,c,b,corner);
        auto t1 = clk::now();
        t_sherman = chrono::duration<double>(t1-t0).count();
    } catch (const exception &e) {
        cerr << "Sherman-Morrison solve failed for N=" << N << ": " << e.what() << "\n";
    }
    // Pojedyncze pomiary: SOR
    try {
        auto t0 = clk::now();
        x_sor = solve_sor_full(a,d,c,b,corner,1e-12,100000,1.0);
        auto t1 = clk::now();
        t_sor = chrono::duration<double>(t1-t0).count();
    } catch (const exception &e) {
        cerr << "SOR solve failed for N=" << N << ": " << e.what() << "\n";
    }

    // Wypis podsumowania: czasy i krótka próbka rozwiązań
    cout << fixed << setprecision(12);
    cout << "N = " << N << "\n";
    cout << "Dense time: " << t_dense << " s\n";
    cout << "Sherman-Morrison time: " << t_sherman << " s\n";
    cout << "SOR time: " << t_sor << " s\n";
    cout << "Sample solution (first 10 entries):\n";
    int M = min(N, 10);
    for (int i = 0; i < M; ++i) cout << "i=" << i << "  " << x_dense[i] << "  " << x_sherman[i] << "  " << x_sor[i] << "\n";
    cout << "CSV_TIMES," << N << "," << t_dense << "," << t_sherman << "," << t_sor << "\n";

    bool save_solution = false;
    string solution_fname;
    if (argc >= 3 && string(argv[2]) == "dump") {
        save_solution = true;
        solution_fname = string("solution_N") + to_string(N) + string(".csv");
    }

    // Automatyczne zapisywanie dla N >= 1000 
    // Automatyczne zapisywanie rozwiązania dla N >= 1000, jeśli nie wyłączono jawnie
    // Automatyczne zapisanie pełnego rozwiązania dla N >= 1000 lub gdy wywołano
    // program z opcją dump
    bool will_save = save_solution || (N >= 1000 && !save_solution);
    if (will_save) {
        double h = 2.0 / (N - 1);
        string fname = save_solution ? solution_fname : (string("solution_N") + to_string(N) + string(".csv"));
        ofstream fout(fname);
        if (!fout) {
            cerr << "Failed to open solution file for writing: " << fname << "\n";
        } else {
            fout << "x,dense,sherman_morrison,sor\n";
            for (int i = 0; i < N; ++i) {
                double xcoord = (double)i * h;
                fout << xcoord << "," << x_dense[i] << "," << x_sherman[i] << "," << x_sor[i] << "\n";
            }
            fout.close();
            cout << "Saved solution to: " << fname << "\n";
        }
    }
    return 0;
}
