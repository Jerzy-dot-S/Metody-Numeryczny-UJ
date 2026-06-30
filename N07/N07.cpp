// N07: Rozwiązanie nieliniowego układu równań różnicowych
//
// Zadanie:
//  - znaleźć przybliżone rozwiązanie {u_n} spełniające dla n = 1..N-1:
//      -(u_{n-1} - 2 u_n + u_{n+1}) / h^2 + 2 u_n (u_n^2 - 1) = 0,
//    z warunkami brzegowymi u_0 = 0, u_N = 1,
//    gdzie h = 20 / (N - 1).
//  - N ~ 100–1000, dokładność ~ 1e-6.
//  - narysować wykres rozwiązania u(x_n), x_n = n h.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#endif

// Rozwiązywanie układu trójdiagonalnego metodą Thomasa:
//  lower[i]  - element pod diagonalą (i = 1..m-1)
//  diag[i]   - elementy diagonalne (i = 0..m-1)
//  upper[i]  - element nad diagonalą (i = 0..m-2)
//  rhs[i]    - prawa strona
//  solution  - wynik (rozmiar m)
static bool solve_tridiagonal(const std::vector<double> &lower,
                              const std::vector<double> &diag,
                              const std::vector<double> &upper,
                              const std::vector<double> &rhs,
                              std::vector<double> &solution)
{
	int m = static_cast<int>(diag.size());
	if (m == 0) return false;

	std::vector<double> c(m, 0.0);
	std::vector<double> d(m, 0.0);

	// Modyfikacja pierwszego wiersza
	if (std::fabs(diag[0]) < 1e-14) return false;
	c[0] = (m > 1) ? upper[0] / diag[0] : 0.0;
	d[0] = rhs[0] / diag[0];

	// Eliminacja w przód
	for (int i = 1; i < m; ++i)
	{
		double denom = diag[i] - lower[i] * c[i - 1];
		if (std::fabs(denom) < 1e-14) return false;
		c[i] = (i < m - 1) ? upper[i] / denom : 0.0;
		d[i] = (rhs[i] - lower[i] * d[i - 1]) / denom;
	}

	// Podstawianie wstecz
	solution.assign(m, 0.0);
	solution[m - 1] = d[m - 1];
	for (int i = m - 2; i >= 0; --i)
	{
		solution[i] = d[i] - c[i] * solution[i + 1];
	}

	return true;
}

int main()
{
	#ifdef _WIN32
	// Wymuszenie UTF-8 w konsoli Windows, żeby polskie znaki wyświetlały się poprawnie
	SetConsoleOutputCP(CP_UTF8);
	SetConsoleCP(CP_UTF8);
	#endif

	int N;
	std::printf("Podaj N (100..1000): ");
	if (std::scanf("%d", &N) != 1 || N < 3)
	{
		std::fprintf(stderr, "Błędne N.\n");
		return 1;
	}

	// Dla zadania: h = 20 / (N - 1), przedział długości 20.
	double h = 20.0 / (N - 1);
	int M = N - 1; // liczba niewiadomych: u_1..u_{N-1}

	// Wektor u[0..N]; u_0 = 0, u_N = 1, reszta będzie iterowana.
	std::vector<double> u(N + 1, 0.0);
	u[0] = 0.0;
	u[N] = 1.0;
	for (int n = 1; n < N; ++n)
	{
		// Początkowe przybliżenie – prosta od 0 do 1
		u[n] = static_cast<double>(n) / N;
	}

	const double tol = 1e-6;
	const int maxNewtonIters = 50;

	std::vector<double> F(M, 0.0);
	std::vector<double> lower(M, 0.0), diag(M, 0.0), upper(M, 0.0);
	std::vector<double> delta(M, 0.0);

	int iter;
	for (iter = 0; iter < maxNewtonIters; ++iter)
	{
		// Oblicz residual F i macierz Jacobiego (trójdiagonalna)
		double inv_h2 = 1.0 / (h * h);

		for (int n = 1; n <= N - 1; ++n)
		{
			int i = n - 1; // indeks w wektorach F, lower, diag, upper
			double un = u[n];
			double unm1 = u[n - 1];
			double unp1 = u[n + 1];

			// F_n(u) = -(u_{n-1} - 2u_n + u_{n+1}) / h^2 + 2 u_n (u_n^2 - 1)
			double termFD = -(unm1 - 2.0 * un + unp1) * inv_h2;
			double termNL = 2.0 * un * (un * un - 1.0);
			F[i] = termFD + termNL;

			// Pochodne: dF_n/du_{n-1}, dF_n/du_n, dF_n/du_{n+1}
			double dFd_unm1 = -inv_h2;
			double dFd_un = 2.0 * inv_h2;
			double dFd_unp1 = -inv_h2;
			double dNLd_un = 6.0 * un * un - 2.0; // pochodna 2u(u^2-1)

			diag[i] = dFd_un + dNLd_un;
			lower[i] = (n > 1) ? dFd_unm1 : 0.0;
			upper[i] = (n < N - 1) ? dFd_unp1 : 0.0;
		}

		// Rozwiązujemy J * delta = -F
		std::vector<double> rhs(M, 0.0);
		for (int i = 0; i < M; ++i)
			rhs[i] = -F[i];

		if (!solve_tridiagonal(lower, diag, upper, rhs, delta))
		{
			std::fprintf(stderr, "Błąd: nie udało się rozwiązać układu liniowego w iteracji Newtona.\n");
			break;
		}

		// Aktualizacja u i sprawdzenie zbieżności
		double maxDelta = 0.0;
		for (int n = 1; n <= N - 1; ++n)
		{
			int i = n - 1;
			u[n] += delta[i];
			maxDelta = std::max(maxDelta, std::fabs(delta[i]));
		}

		if (maxDelta < tol)
			break;
	}

	if (iter == maxNewtonIters)
	{
		std::fprintf(stderr, "Ostrzeżenie: Newton nie osiągnął zadanej dokładności w %d iteracjach.\n", maxNewtonIters);
	}

	// Zapis rozwiązania do pliku: kolumny x, u(x)
	FILE *f = std::fopen("N07_solution.dat", "w");
	if (!f)
	{
		std::perror("Nie można otworzyć pliku N07_solution.dat");
		return 1;
	}

	std::fprintf(f, "# x\tu(x)\n");
	for (int n = 0; n <= N; ++n)
	{
		double x = n * h;
		std::fprintf(f, "% .10f\t% .10f\n", x, u[n]);
	}
	std::fclose(f);

	std::printf("Zakończono po %d iteracjach Newtona. Wynik zapisany w pliku N07_solution.dat.\n", iter + 1);

	// Automatyczne uruchomienie gnuplota z plikiem plot_gnuplot.plt (jak w N04)
	int rc = std::system("gnuplot plot_gnuplot.plt");
	if (rc == 0)
	{
		std::printf("Wygenerowano wykres: N07_solution.png\n");
	}
	else
	{
		std::printf("Nie znaleziono gnuplota lub rysowanie nie powiodło się.\n");
		std::printf("Możesz samodzielnie uruchomić:\n");
		std::printf("  gnuplot plot_gnuplot.plt\n");
	}

	return 0;
}