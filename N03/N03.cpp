#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#endif


// Prosta iteracja potęgowa dla największego wektora własnego
std::vector<double> power_iteration(const std::vector<std::vector<double>>& A, int num_iterations) {
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");

	std::mt19937 rng(12345);
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	std::vector<double> b(n);
	for (int i = 0; i < n; ++i) b[i] = dist(rng);

	std::vector<double> b1(n);
	for (int it = 0; it < num_iterations; ++it) {
		for (int i = 0; i < n; ++i) {
			double s = 0.0;
			for (int j = 0; j < n; ++j) s += A[i][j] * b[j];
			b1[i] = s;
		}
		double norm = 0.0;
		for (int i = 0; i < n; ++i) norm += b1[i] * b1[i];
		norm = std::sqrt(norm);
		for (int i = 0; i < n; ++i) b[i] = b1[i] / norm;
	}
	return b;
}

// Drugi co do wielkości wektor własny – iteracja potęgowa z ortogonalizacją do e1
std::vector<double> power_iteration_second(const std::vector<std::vector<double>>& A,
	const std::vector<double>& e1,
	int num_iterations) {
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");
	if (static_cast<int>(e1.size()) != n) throw std::invalid_argument("e1 has wrong size");

	std::mt19937 rng(67890);
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	std::vector<double> y(n);
	for (int i = 0; i < n; ++i) y[i] = dist(rng);

	// y_1 prostopadły do e1: y <- y - e1 (e1^T y)
	double proj = 0.0;
	for (int i = 0; i < n; ++i) proj += e1[i] * y[i];
	for (int i = 0; i < n; ++i) y[i] -= proj * e1[i];
	double norm = 0.0;
	for (int i = 0; i < n; ++i) norm += y[i] * y[i];
	norm = std::sqrt(norm);
	for (int i = 0; i < n; ++i) y[i] /= norm;

	std::vector<double> z(n);
	for (int it = 0; it < num_iterations; ++it) {
		// z_k = A y_k
		for (int i = 0; i < n; ++i) {
			double s = 0.0;
			for (int j = 0; j < n; ++j) s += A[i][j] * y[j];
			z[i] = s;
		}
		// Reortogonalizacja: z_k <- z_k - e1 (e1^T z_k)
		double proj_z = 0.0;
		for (int i = 0; i < n; ++i) proj_z += e1[i] * z[i];
		for (int i = 0; i < n; ++i) z[i] -= proj_z * e1[i];
		// Normalizacja y_{k+1} = z_k / ||z_k||
		double norm_z = 0.0;
		for (int i = 0; i < n; ++i) norm_z += z[i] * z[i];
		norm_z = std::sqrt(norm_z);
		for (int i = 0; i < n; ++i) y[i] = z[i] / norm_z;
	}
	return y;
}

// Rozkład LU z częściowym wyborem elementu głównego
void lu_factor(std::vector<std::vector<double>>& A, std::vector<int>& piv) {
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");

	piv.resize(n);
	for (int i = 0; i < n; ++i) piv[i] = i;

	for (int k = 0; k < n; ++k) {
		int pivot = k;
		double maxAbs = std::abs(A[k][k]);
		for (int i = k + 1; i < n; ++i) {
			if (std::abs(A[i][k]) > maxAbs) {
				maxAbs = std::abs(A[i][k]);
				pivot = i;
			}
		}
		if (maxAbs < 1e-14) throw std::runtime_error("Singular matrix in lu_factor");
		if (pivot != k) {
			std::swap(A[k], A[pivot]);
			std::swap(piv[k], piv[pivot]);
		}

		double diag = A[k][k];
		for (int i = k + 1; i < n; ++i) {
			double factor = A[i][k] / diag;
			A[i][k] = factor; // element L
			for (int j = k + 1; j < n; ++j) {
				A[i][j] -= factor * A[k][j]; // aktualizacja U
			}
		}
	}
}

// Rozwiązywanie układu LU x = b 
std::vector<double> lu_solve(const std::vector<std::vector<double>>& LU, const std::vector<int>& piv, const std::vector<double>& b) {
	int n = static_cast<int>(LU.size());
	if (n == 0 || static_cast<int>(LU[0].size()) != n) throw std::invalid_argument("LU must be square");
	if (static_cast<int>(b.size()) != n) throw std::invalid_argument("b has wrong size");
	if (static_cast<int>(piv.size()) != n) throw std::invalid_argument("piv has wrong size");

	// Zastosowanie permutacji P do b
	std::vector<double> x(n);
	for (int i = 0; i < n; ++i) x[i] = b[piv[i]];

	// Rozwiązanie L y = P b (L ma jedynki na diagonali)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j) x[i] -= LU[i][j] * x[j];
	}

	// Rozwiązanie U x = y
	for (int i = n - 1; i >= 0; --i) {
		for (int j = i + 1; j < n; ++j) x[i] -= LU[i][j] * x[j];
		x[i] /= LU[i][i];
	}

	return x;
}

// Iloraz Rayleigha do przybliżenia wartości własnej dla wektora v
double rayleigh(const std::vector<std::vector<double>>& A, const std::vector<double>& v) {
	int n = static_cast<int>(v.size());
	std::vector<double> Av(n);
	for (int i = 0; i < n; ++i) {
		double s = 0.0;
		for (int j = 0; j < n; ++j) s += A[i][j] * v[j];
		Av[i] = s;
	}
	double num = 0.0, den = 0.0;
	for (int i = 0; i < n; ++i) { num += v[i] * Av[i]; den += v[i] * v[i]; }
	return num / den;
}

// Iteracyjna metoda Rayleigha dla pojedynczej pary własnej
void rayleigh_iteration(const std::vector<std::vector<double>>& A,
	std::vector<double>& v,
	double& lambda,
	int max_iters,
	double tol)
{
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");
	if (static_cast<int>(v.size()) != n) throw std::invalid_argument("v has wrong size");

	// Normalizacja wektora startowego
	double norm = 0.0;
	for (int i = 0; i < n; ++i) norm += v[i] * v[i];
	norm = std::sqrt(norm);
	for (int i = 0; i < n; ++i) v[i] /= norm;

	for (int it = 0; it < max_iters; ++it) {
		lambda = rayleigh(A, v);

		// Budujemy przesuniętą macierz A - lambda I
		std::vector<std::vector<double>> Ashift = A;
		for (int i = 0; i < n; ++i) Ashift[i][i] -= lambda;

		// Rozkład LU(A - lambda I)
		std::vector<std::vector<double>> LU = Ashift;
		std::vector<int> piv;
		lu_factor(LU, piv);

		// Rozwiązujemy (A - lambda I) w = v
		std::vector<double> w = lu_solve(LU, piv, v);

		// Normalizacja nowego wektora
		double norm_w = 0.0;
		for (int i = 0; i < n; ++i) norm_w += w[i] * w[i];
		norm_w = std::sqrt(norm_w);
		for (int i = 0; i < n; ++i) v[i] = w[i] / norm_w;

		// Sprawdzamy residuum ||A v - lambda v||
		double res2 = 0.0;
		for (int i = 0; i < n; ++i) {
			double s = 0.0;
			for (int j = 0; j < n; ++j) s += A[i][j] * v[j];
			s -= lambda * v[i];
			res2 += s * s;
		}
		if (std::sqrt(res2) < tol) break;
	}

	// Ostateczne odświeżenie wartości własnej
	lambda = rayleigh(A, v);
}

// Odwrotna iteracja potęgowa z ortogonalizacją względem wcześniej znalezionych wektorów
std::vector<double> inverse_power_iteration_orth(const std::vector<std::vector<double>>& A,
	const std::vector<std::vector<double>>& prev_eigvecs,
	int num_iterations)
{
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");

	// Jednorazowy rozkład LU macierzy A
	std::vector<std::vector<double>> LU = A;
	std::vector<int> piv;
	lu_factor(LU, piv);

	std::mt19937 rng(98765);
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	std::vector<double> y(n);
	for (int i = 0; i < n; ++i) y[i] = dist(rng);

	// Ortogonalizacja względem wcześniejszych wektorów własnych
	for (const auto& e : prev_eigvecs) {
		double proj = 0.0;
		for (int i = 0; i < n; ++i) proj += e[i] * y[i];
		for (int i = 0; i < n; ++i) y[i] -= proj * e[i];
	}

	double norm = 0.0;
	for (int i = 0; i < n; ++i) norm += y[i] * y[i];
	norm = std::sqrt(norm);
	for (int i = 0; i < n; ++i) y[i] /= norm;

	for (int it = 0; it < num_iterations; ++it) {
		std::vector<double> z = lu_solve(LU, piv, y);
		// Reortogonalizacja względem wcześniejszych wektorów
		for (const auto& e : prev_eigvecs) {
			double proj = 0.0;
			for (int i = 0; i < n; ++i) proj += e[i] * z[i];
			for (int i = 0; i < n; ++i) z[i] -= proj * e[i];
		}
		double norm_z = 0.0;
		for (int i = 0; i < n; ++i) norm_z += z[i] * z[i];
		norm_z = std::sqrt(norm_z);
		for (int i = 0; i < n; ++i) y[i] = z[i] / norm_z;
	}
	return y;
}

// Dwie największe wartości własne metodą Rayleigha
std::vector<double> find_two_largest_rayleigh(const std::vector<std::vector<double>>& A,
	int iters_power,
	int iters_rayleigh,
	double tol)
{
	// Najpierw przybliżamy dwa największe wektory własne metodą potęgową,
	// a następnie doszlifowujemy je metodą Rayleigha.
	std::vector<double> lambdas(2);
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");

	// Przybliżenia wektorów własnych z metody potęgowej
	std::vector<double> v1 = power_iteration(A, iters_power);
	std::vector<double> v2 = power_iteration_second(A, v1, iters_power);

	double lambda1 = 0.0;
	rayleigh_iteration(A, v1, lambda1, iters_rayleigh, tol);
	lambdas[0] = lambda1;

	double lambda2 = 0.0;
	rayleigh_iteration(A, v2, lambda2, iters_rayleigh, tol);
	lambdas[1] = lambda2;

	return lambdas;
}

// Cztery najmniejsze wartości własne metodą Rayleigha
std::vector<double> find_four_smallest_rayleigh(const std::vector<std::vector<double>>& A,int iters_inverse,int iters_rayleigh,double tol)
{
	// Najpierw znajdujemy przybliżenia czterech najmniejszych wartości własnych
	// odwrotną metodą potęgową z ortogonalizacją, a następnie poprawiamy je
	// iteracją Rayleigha.
	std::vector<double> lambdas;
	lambdas.reserve(4);
	std::vector<std::vector<double>> eigvecs; // wektory własne do ortogonalizacji
	int n = static_cast<int>(A.size());
	if (n == 0 || static_cast<int>(A[0].size()) != n) throw std::invalid_argument("A must be square");

	for (int k = 0; k < 4; ++k) {
		// Przybliżony wektor własny k-tej najmniejszej wartości własnej
		std::vector<double> v = inverse_power_iteration_orth(A, eigvecs, iters_inverse);
		double lambda_k = 0.0;
		// Doszlifowanie metodą Rayleigha
		rayleigh_iteration(A, v, lambda_k, iters_rayleigh, tol);
		lambdas.push_back(lambda_k);
		eigvecs.push_back(v);
	}

	return lambdas;
}

// Cztery najmniejsze wartości własne odwrotną metodą potęgową
std::vector<double> find_four_smallest_power(const std::vector<std::vector<double>>& A,
	int iters_inverse)
{
	std::vector<double> lambdas;
	lambdas.reserve(4);
	std::vector<std::vector<double>> eigvecs;

	for (int k = 0; k < 4; ++k) {
		std::vector<double> v = inverse_power_iteration_orth(A, eigvecs, iters_inverse);
		double lambda_k = rayleigh(A, v);
		lambdas.push_back(lambda_k);
		eigvecs.push_back(v);
	}

	return lambdas;
}

//Konstrukcja macierzy A
std::vector<std::vector<double>> generateA(int N) {
	if (N < 2) {
		throw std::invalid_argument("N must be at least 2");
	}

	const double h = 20.0 / static_cast<double>(N - 1);
	const double inv_h2 = 1.0 / (h * h);

	std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));

	for (int i = 0; i < N; ++i) {
		const int i_math = i + 1; 
		const double x_i = i_math * h; 
		const double V_i = (x_i - 10.0) * (x_i - 10.0);

		// Diagonal
		A[i][i] = -2.0 * inv_h2 + V_i;

		// Off-diagonals
		if (i + 1 < N) {
			A[i][i + 1] = inv_h2;
		}
		if (i - 1 >= 0) {
			A[i][i - 1] = inv_h2;
		}
	}

	return A;
}

int main() {
	try {
		std::cout << std::fixed << std::setprecision(6);

		#ifdef _WIN32
		SetConsoleOutputCP(CP_UTF8);
		SetConsoleCP(CP_UTF8);
		#endif

		// Parametry liczby iteracji
		int N = 500; // Wielkość macierzy 
		int itersPowerDominant   = 100;   // iteracje metody potęgowej (największa wartość własna)
		int itersPowerSecond     = 200;   // iteracje dla drugiej wartości własnej (potęgowa)
		int itersInverseSmall    = 60;    // iteracje odwrotnej metody potęgowej (najmniejsze wartości własne)
		int itersRayleighInit    = 60;    // ile iteracji potęgowej przed startem iteracji Rayleigha
		int itersRayleighMax     = 200;  // maksymalna liczba iteracji metody Rayleigha
		int itersPowerCompare    = 100;   // iteracje potęgowej w porównaniu (ten sam wektor startowy)
		int itersRayleighCompare = 50;    // maks. iteracji Rayleigha w porównaniu
		double tolRayleigh       = 1e-10; // tolerancja dla residuum w metodzie Rayleigha
		
		auto A = generateA(N);
		// Metoda potęgowa (zwykła)
		auto vN = power_iteration(A, itersPowerDominant);
		double lambdaN = rayleigh(A, vN);
		auto vN2 = power_iteration_second(A, vN, itersPowerSecond);
		double lambdaN2 = rayleigh(A, vN2);
		// Cztery najmniejsze odwrotną metodą potęgową
		auto four_small_power = find_four_smallest_power(A, itersInverseSmall);

		// Metoda Rayleigha
		auto two_big   = find_two_largest_rayleigh(A, itersRayleighInit, itersRayleighMax, tolRayleigh);
		auto four_small = find_four_smallest_rayleigh(A, itersInverseSmall, itersRayleighMax, tolRayleigh);

		std::cout << "\nLarge A (N=" << N << ") eigenvalues comparison:\n";
		std::cout << "  Largest:       power = " << lambdaN
			          << ",  Rayleigh = " << two_big[0] << "\n";
		std::cout << "  Second largest: power = " << lambdaN2
			          << ",  Rayleigh = " << two_big[1] << "\n";
		std::cout << "  Four smallest by |lambda|:\n";
		for (std::size_t i = 0; i < four_small.size() && i < four_small_power.size(); ++i) {
			std::cout << "    #" << (i + 1)
				      << ": power = " << four_small_power[i]
				      << ",  Rayleigh = " << four_small[i] << "\n";
		}

		// Dodatkowe: porównanie zbieżności dla największej wartości własnej,
		// przy tym samym losowym wektorze startowym dla obu metod.
		{
			int nDim = static_cast<int>(A.size());
			std::mt19937 rng_same(11111);
			std::uniform_real_distribution<double> dist_same(0.0, 1.0);
			std::vector<double> v0(nDim);
			for (int i = 0; i < nDim; ++i) v0[i] = dist_same(rng_same);
			// normalizacja v0
			double norm0 = 0.0;
			for (int i = 0; i < nDim; ++i) norm0 += v0[i] * v0[i];
			norm0 = std::sqrt(norm0);
			for (int i = 0; i < nDim; ++i) v0[i] /= norm0;

			std::vector<double> v_power = v0;
			std::vector<double> v_rayleigh = v0;

			// Kilka iteracji metody potęgowej z tym samym startem 
			for (int it = 0; it < itersPowerCompare; ++it) {
				std::vector<double> tmp(nDim);
				for (int i = 0; i < nDim; ++i) {
					double s = 0.0;
					for (int j = 0; j < nDim; ++j) s += A[i][j] * v_power[j];
					tmp[i] = s;
				}
				double norm = 0.0;
				for (int i = 0; i < nDim; ++i) norm += tmp[i] * tmp[i];
				norm = std::sqrt(norm);
				for (int i = 0; i < nDim; ++i) v_power[i] = tmp[i] / norm;
			}
			double lambda_power_same = rayleigh(A, v_power);

			// Iteracja Rayleigha z tym samym startem
			double lambda_rayleigh_same = 0.0;
			rayleigh_iteration(A, v_rayleigh, lambda_rayleigh_same, itersRayleighCompare, tolRayleigh);

			std::cout << "\nPorównanie od tego samego losowego wektora startowego (największa wart. własna):\n";
			std::cout << "  Power   lambda ≈ " << lambda_power_same << "\n";
			std::cout << "  Rayleigh lambda ≈ " << lambda_rayleigh_same << "\n";
		}

		return 0;
	} catch (const std::exception &ex) {
		std::cerr << "Error: " << ex.what() << "\n";
		return 1;
	}
}

