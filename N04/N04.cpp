// N04: Interpolacja Lagrange'a funkcji f(x) = 1/(1+x^2) na przedziale [-5,5]
// Dwa zbiory węzłów:
// (a) Równoodległe: X_n = -5 + 10*n/N, n = 0..N
// (b) Węzły cosinusowe: X_n = 5*cos(n*pi/N), n = 0..N
// Program przeszukuje N, aby zminimalizować maksymalny błąd na gęstej siatce,
// a następnie zapisuje dane i wypisuje najlepsze N dla każdego zbioru węzłów.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <fstream>
#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#endif

using std::cout;
using std::endl;
using std::vector;

static inline double f(double x) {
	return 1.0 / (1.0 + x * x);
}

// Funkcja do obliczania p(x0) używając wzory interpolacyjnego Lagrange'a
static double lagrange(const vector<double>& x, const vector<double>& y, double x0) {
	const size_t n = x.size();
	// jeśli x0 jest dokładnie jednym z węzłów, zwracamy odpowiadającą wartość
	for (size_t j = 0; j < n; ++j) {
		if (x0 == x[j]) return y[j];
	}
	double px = 0.0;
	for (size_t j = 0; j < n; ++j) {
		double Lj = 1.0;
		for (size_t k = 0; k < n; ++k) {
			if (k == j) continue;
			Lj *= (x0 - x[k]) / (x[j] - x[k]);
		}
		px += Lj * y[j];
	}
	return px;
}

//Funkcja do generowania siatki do błędu. Sprawdzamy błąd
//dla każdego punktu w przedziale i stąd wyliczamy przybliżony 
//maksymalny błąd
static vector<double> make_grid(double a, double b, int m) {
	// tworzymy wektor m punktów równomiernie rozmieszczonych w [a,b]
	vector<double> g;
	g.reserve(static_cast<size_t>(m));
	double h = (b - a) / (m - 1); // krok siatki
	for (int i = 0; i < m; ++i) g.push_back(a + i * h);
	return g;
}

//Funkcja do ustawienia punktów między którymi znajdujemy interpolacje
//dla X_n = -5 + 10*n/N (węzły równoodległe)
static vector<double> make_wezly_rownolegle(int N, double a, double b) {
	int points = N + 1;
	vector<double> x(points);
	double h = (b - a) / (points - 1);
	for (int n = 0; n < points; ++n) x[n] = a + h * n;
	return x;
}
	
	//Funkcja do ustawienia punktów między którymi znajdujemy interpolacje
	//dla X_n = 5*cos(n*pi/N) (węzły cosinusowe)
static vector<double> make_wezly_cos(int N) {
	int points = N + 1;
	vector<double> x(points);
	const double pi = std::acos(-1.0);
	double denom = static_cast<double>(N);
	for (int n = 0; n < points; ++n) {
		x[n] = 5.0 * std::cos(pi * n / denom);
	}
	return x;
}

//Funckja do liczenia przybliżonego maksymalnego błędu.
static double max_error_on_grid(const vector<double>& grid,const vector<double>& x,const vector<double>& y) {
	double maxerr = 0.0;
	for (double gx : grid) {
		double fa = lagrange(x, y, gx);
		double err = std::fabs(f(gx) - fa);
		if (err > maxerr) maxerr = err;
	}
	return maxerr;
}

int main(int argc, char** argv) {
	// Parametry zadania
	const double a = -5.0, b = 5.0;

#ifdef _WIN32
	// Wymuszamy UTF-8 dla polskich znaków
	SetConsoleOutputCP(CP_UTF8);
	SetConsoleCP(CP_UTF8);
#endif

	int gridSize = 4001; // Siatka do oceny błędu
	int Nmin = 1;        // od którego N liczymy
	int Nmax = 60;       // maksymalny badany stopień wielomianu

	auto grid = make_grid(a, b, gridSize);

	// Przejście po N, obliczenie maksymalnego błędu dla każdego zbioru węzłów
	std::ofstream errs("errors.dat");
	errs << std::fixed << std::setprecision(12);
	errs << "# N\tErrRownolegle\tErrCos" << '\n';

	double bestErrEqui = std::numeric_limits<double>::infinity();
	int bestNEqui = Nmin;

	double bestErrCos = std::numeric_limits<double>::infinity();
	int bestNCos = Nmin;

	//Główna pętla dla znajdywania najlepszego N
	for (int N = Nmin; N <= Nmax; ++N) {

		// Węzły równoodległe
		auto xe = make_wezly_rownolegle(N, a, b);
		vector<double> ye(xe.size());
		for (size_t i = 0; i < xe.size(); ++i) ye[i] = f(xe[i]);
		double errEqui = max_error_on_grid(grid, xe, ye);

		// Węzły cosinusowe
		double errCos = std::numeric_limits<double>::infinity();
		if (N >= 3) {
			auto xc = make_wezly_cos(N);
			vector<double> yc(xc.size());
			for (size_t i = 0; i < xc.size(); ++i) yc[i] = f(xc[i]);
			errCos = max_error_on_grid(grid, xc, yc);
		}

		errs << N << '\t' << errEqui << '\t' << errCos << '\n';

		if (errEqui < bestErrEqui) {
			bestErrEqui = errEqui;
			bestNEqui = N;
		}
		if (errCos < bestErrCos) {
			bestErrCos = errCos;
			bestNCos = N;
		}
	}
	errs.close();

	cout << "Najlepsze N (węzły równoodległe): " << bestNEqui << ", max error = " << std::setprecision(12) << bestErrEqui << endl;
	cout << "Najlepsze N (węzły cosinusowe): " << bestNCos << ", max error = " << std::setprecision(12) << bestErrCos << endl;

	//Generowanie danych dla najlepszych N
	std::ofstream out("out_best.dat");
	out << std::fixed << std::setprecision(12);
	out << "# x\tf(x)\tP_equidistant\tP_cosineSpec" << '\n';

	//Interpolanty dla najlepszych N
	auto xe_best = make_wezly_rownolegle(bestNEqui, a, b);
	vector<double> ye_best(xe_best.size());
	for (size_t i = 0; i < xe_best.size(); ++i) ye_best[i] = f(xe_best[i]);

	auto xc_best = make_wezly_cos(bestNCos);
	vector<double> yc_best(xc_best.size());
	for (size_t i = 0; i < xc_best.size(); ++i) yc_best[i] = f(xc_best[i]);

	//Wartości do zapisania dla wykresu
	for (double gx : grid) {
		double fx = f(gx);
		double pe = lagrange(xe_best, ye_best, gx);
		double pc = lagrange(xc_best, yc_best, gx);
		out << gx << '\t' << fx << '\t' << pe << '\t' << pc << '\n';
	}
	out.close();

	// Zapisz węzły dla najlepszych N 
	{
		std::ofstream ne("nodes_rowno.dat");
		ne << std::fixed << std::setprecision(12);
		ne << "# x\tf(x)" << '\n';
		for (size_t i = 0; i < xe_best.size(); ++i) ne << xe_best[i] << '\t' << ye_best[i] << '\n';
	}
	{
		std::ofstream nc("nodes_cos.dat");
		nc << std::fixed << std::setprecision(12);
		nc << "# x\tf(x)" << '\n';
		for (size_t i = 0; i < xc_best.size(); ++i) nc << xc_best[i] << '\t' << yc_best[i] << '\n';
	}

	cout << "Dane zapisane do plików out_best.dat i errors.dat" << endl;

	std::string scriptPath = "plot_gnuplot.plt";

	int rc = -1;
	if (!scriptPath.empty()) {
		std::string cmd = std::string("gnuplot ") + scriptPath;
		rc = std::system(cmd.c_str());
	}
	if (rc == 0) {
		cout << "Wygenerowano wykresy: plot_best.png, plot_errors.png" << endl;
	} else {
		cout << "Nie znaleziono gnuplota lub rysowanie nie powiodło się. Spróbuj: gnuplot plot_gnuplot.plt" << endl;
	}
	return 0;
}

