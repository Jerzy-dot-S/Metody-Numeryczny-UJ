// N05: Numeryczne obliczanie calki
//   I = \int_{-1}^{1} exp(x^2) / sqrt(1 - x^2) dx
// za pomoca zlozonej metody trapezow, Simpsona i reguly 3/8,
// z iteracyjnym zageszczaniem podprzedzialow (bez ponownego
// liczenia wartosci funkcji na istniejacych wezach).
//
// Podstawienie: x = sin(t), t in [-pi/2, pi/2]
// dx = cos(t) dt, sqrt(1-x^2) = sqrt(1-sin^2 t) = cos(t) (na tym przedziale >=0)
// => calka I = \int_{-pi/2}^{pi/2} exp(sin^2 t) dt.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

// Funkcja podcalkowa po podstawieniu x = sin(t)
static inline double g(double t)
{
    double s = std::sin(t);
    return std::exp(s * s);
}

// Oblicza przy danym n (liczba podprzedzialow) przyblizenia calki
// metoda trapezow, Simpsona i reguly 3/8 na podstawie tablicy f[0..n]
// (f[i] = g(t_i), t_i rownomierne na [a,b]).
static void compute_integrals(double a, double b,
                              const std::vector<double>& f,
                              int n,
                              double& T, double& S, double& R38)
{
    double h = (b - a) / static_cast<double>(n);

    // Zlozona metoda trapezow
    double sumT = 0.5 * (f[0] + f[n]);
    for (int i = 1; i < n; ++i)
        sumT += f[i];
    T = h * sumT;

    // Zlozona metoda Simpsona (n parzyste, tutaj zawsze n wielokrotnosc 2)
    double sum_odd = 0.0;
    double sum_even = 0.0;
    for (int i = 1; i < n; ++i)
    {
        if (i % 2 == 1)
            sum_odd += f[i];
        else
            sum_even += f[i];
    }
    S = (h / 3.0) * (f[0] + f[n] + 4.0 * sum_odd + 2.0 * sum_even);

    // Zlozona regula 3/8 (wymaga, aby n bylo wielokrotnoscia 3)
    double sum_mult3 = 0.0;      // indeksy podzielne przez 3 (bez 0 i n)
    double sum_not_mult3 = 0.0;  // pozostale indeksy wewnetrzne
    for (int i = 1; i < n; ++i)
    {
        if (i % 3 == 0)
            sum_mult3 += f[i];
        else
            sum_not_mult3 += f[i];
    }
    R38 = (3.0 * h / 8.0) * (f[0] + f[n] + 3.0 * sum_not_mult3 + 2.0 * sum_mult3);
}

int main()
{
#ifdef _WIN32
    // Wymuszenie UTF-8 w konsoli Windows
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    const double a = -0.5 * M_PI;
    const double b =  0.5 * M_PI;
    const double tol = 1e-6;

    // Zaczynamy od n = 6 (wielokrotnosc 2 i 3); przy kazdym kroku n <- 2n
    int n = 6;
    const int maxRefinements = 20; // zabezpieczenie

    // Inicjalna siatka i wartosci funkcji
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; ++i)
    {
        double t = a + (b - a) * static_cast<double>(i) / static_cast<double>(n);
        f[i] = g(t);
    }

    double T_prev, S_prev, R_prev;
    double T_curr, S_curr, R_curr;

    compute_integrals(a, b, f, n, T_prev, S_prev, R_prev);

    bool doneT = false, doneS = false, doneR = false;

    for (int level = 0; level < maxRefinements; ++level)
    {
        // Zigeszczamy: n_new = 2*n, nowa siatka zawiera stare wezly + nowe w srodkach
        int n_new = 2 * n;
        std::vector<double> f_new(n_new + 1);

        // Kopiujemy stare wartosci na parzyste indeksy
        for (int i = 0; i <= n; ++i)
        {
            f_new[2 * i] = f[i];
        }
        // Liczymy nowe punkty (niepowtarzajace sie wezly)
        for (int i = 1; i < n_new; i += 2)
        {
            double t = a + (b - a) * static_cast<double>(i) / static_cast<double>(n_new);
            f_new[i] = g(t);
        }

        f.swap(f_new);
        n = n_new;

        compute_integrals(a, b, f, n, T_curr, S_curr, R_curr);

        if (!doneT && std::fabs(T_curr - T_prev) < tol)
            doneT = true;
        if (!doneS && std::fabs(S_curr - S_prev) < tol)
            doneS = true;
        if (!doneR && std::fabs(R_curr - R_prev) < tol)
            doneR = true;

        T_prev = T_curr;
        S_prev = S_curr;
        R_prev = R_curr;

        if (doneT && doneS && doneR)
            break;
    }

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  metoda trapezow : " << T_prev << "\n";
    std::cout << "  metoda Simpsona : " << S_prev << "\n";
    std::cout << "  regula 3/8      : " << R_prev << "\n";
    std::cout << "Uzyta liczba podprzedzialow (n) = " << n << "\n";

    return 0;
}
