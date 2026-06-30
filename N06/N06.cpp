#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include <utility>
#include <iostream>

double horner(double x, int n, const std::vector<double> &P)
{
	double result = P[n];
	for (int i = n - 1; i >= 0; i--)
	{
		result = result * x + P[i];
	}
	return result;
}

double dhorner(double x, int n, const std::vector<double> &P)
{
	double result = P[n] * n;
	for (int i = n - 1; i >= 1; i--)
	{
		result = result * x + P[i] * i;
	}
	return result;
}

double d2horner(double x, int n, const std::vector<double> &P)
{
	double result = P[n] * n * (n - 1);
	for (int i = n - 1; i >= 2; i--)
	{
		result = result * x + P[i] * i * (i - 1);
	}
	return result;
}

// Complex Horner using pairs (Re, Im)
void hornerZ(double &resRe, double &resIm, double xRe, double xIm, int n, const std::vector<double> &P)
{
	resRe = P[n];
	resIm = 0.0;
	for (int i = n - 1; i >= 0; i--)
	{
		double origResRe = resRe;
		// (res * x)
		resRe = (resRe * xRe - resIm * xIm);
		resIm = (resIm * xRe + origResRe * xIm);
		// + P[i]
		resRe += P[i];
	}
}

void dhornerZ(double &resRe, double &resIm, double a, double b, int n, const std::vector<double> &P)
{
	resRe = P[n] * n;
	resIm = 0.0;
	for (int i = n - 1; i >= 1; i--)
	{
		double origResRe = resRe;
		resRe = (resRe * a - resIm * b);
		resIm = (resIm * a + origResRe * b);
		resRe += P[i] * i;
	}
}

void d2hornerZ(double &resRe, double &resIm, double a, double b, int n, const std::vector<double> &P)
{
	resRe = P[n] * n * (n - 1);
	resIm = 0.0;
	for (int i = n - 1; i >= 2; i--)
	{
		double origResRe = resRe;
		resRe = (resRe * a - resIm * b);
		resIm = (resIm * a + origResRe * b);
		resRe += P[i] * i * (i - 1);
	}
}

// Laguerre's method for complex numbers; returns iteration count and final root in (resRe, resIm)
int LaguerreZ(double &resRe, double &resIm, double xRe, double xIm, int n, const std::vector<double> &P)
{
	int loopCnt = 0;
	const int maxIters = 200;
	double maxA = 0.0;
	for (int i = 0; i <= n; i++)
		maxA = std::max(maxA, std::fabs(P[i]));

	if (maxA == 0.0) { resRe = xRe; resIm = xIm; return loopCnt; }

	while (true)
	{
		double hRe, hIm;
		hornerZ(hRe, hIm, xRe, xIm, n, P);
		if (std::sqrt(hRe * hRe + hIm * hIm) < DBL_EPSILON * maxA)
			break;

		double dhRe, dhIm;
		dhornerZ(dhRe, dhIm, xRe, xIm, n, P);

		// G = dh / h
		double denomh = hRe * hRe + hIm * hIm;
		double GRe = (dhRe * hRe + dhIm * hIm) / denomh;
		double GIm = (dhIm * hRe - dhRe * hIm) / denomh;

		// G^2
		double G2Re = GRe * GRe - GIm * GIm;
		double G2Im = 2.0 * GRe * GIm;

		// d2h
		double d2hRe, d2hIm;
		d2hornerZ(d2hRe, d2hIm, xRe, xIm, n, P);

		// d2h/h
		double GRe2 = (d2hRe * hRe + d2hIm * hIm) / denomh;
		double GIm2 = (d2hIm * hRe - d2hRe * hIm) / denomh;

		// H = G^2 - d2h/h
		double HRe = G2Re - GRe2;
		double HIm = G2Im - GIm2;

		// v = (n-1)*(n*H - G^2)
		double vRe = (n - 1) * (n * HRe - G2Re);
		double vIm = (n - 1) * (n * HIm - G2Im);

		// sqrt(v) via magnitude decomposition
		double modv = std::sqrt(vRe * vRe + vIm * vIm);
		double sRe = std::sqrt((modv + vRe) / 2.0);
		double sIm = std::sqrt((modv - vRe) / 2.0);
		if (vIm < 0) sIm = -sIm;

		// denominator choice with larger modulus
		double denom1Re = GRe + sRe;
		double denom1Im = GIm + sIm;
		double denom2Re = GRe - sRe;
		double denom2Im = GIm - sIm;

		double denomRe, denomIm;
		double mod1 = denom1Re * denom1Re + denom1Im * denom1Im;
		double mod2 = denom2Re * denom2Re + denom2Im * denom2Im;
		if (mod1 > mod2)
		{
			denomRe = denom1Re; denomIm = denom1Im;
		}
		else
		{
			denomRe = denom2Re; denomIm = denom2Im;
		}

		// a = n / denom
		double denoma = denomRe * denomRe + denomIm * denomIm;
		double aRe = n * denomRe / denoma;
		double aIm = -n * denomIm / denoma;

		double xmod = std::sqrt(xRe * xRe + xIm * xIm);
		double amod = std::sqrt(aRe * aRe + aIm * aIm);
		if (amod < 2 * DBL_EPSILON * std::max(1.0, xmod))
		{
			xRe -= 0.5 * aRe;
			xIm -= 0.5 * aIm;
			break;
		}

		xRe -= aRe;
		xIm -= aIm;
		loopCnt++;
		if (loopCnt > maxIters) break; // safeguard
	}

	resRe = xRe;
	resIm = xIm;
	return loopCnt;
}

// Real synthetic division by (x - r)
std::vector<double> deflateLinear(const std::vector<double> &P, int n, double r)
{
	std::vector<double> Q(n);
	Q[n - 1] = P[n];
	for (int i = n - 2; i >= 0; --i)
	{
		Q[i] = P[i + 1] + r * Q[i + 1];
	}
	return Q;
}

// Real synthetic division by quadratic x^2 + a x + b
std::vector<double> deflateQuadratic(const std::vector<double> &P, int n, double a, double b)
{
	std::vector<double> S(n - 1, 0.0); // degree becomes n-2 -> size n-1 entries [0..n-2]
	// S[k] corresponds to coefficient for x^k in the quotient
	S[n - 2] = P[n];
	if (n - 3 >= 0) S[n - 3] = P[n - 1] - a * S[n - 2];
	for (int k = n - 4; k >= 0; --k)
	{
		S[k] = P[k + 2] - a * S[k + 1] - b * S[k + 2];
	}
	return S;
}

int main()
{
	// Wielomian: 243x^7 − 486x^6 + 783x^5 − 990x^4 + 558x^3 − 28x^2 − 72x + 16
	const int deg = 7;
	std::vector<double> P = { 16.0, -72.0, -28.0, 558.0, -990.0, 783.0, -486.0, 243.0 };
	int n = deg;

	std::vector<std::pair<double, double>> roots;
	const double tolIm = 1e-12;

	while (n > 0)
	{
		// Find one root with initial guess (0, 0)
		double rRe = 0.0, rIm = 0.0;
		LaguerreZ(rRe, rIm, 0.0, 0.0, n, P);

		if (std::fabs(rIm) <= tolIm)
		{
			// Real root; deflate linearly
			roots.emplace_back(rRe, 0.0);
			auto Q = deflateLinear(P, n, rRe);
			P.assign(Q.begin(), Q.end());
			n -= 1;
		}
		else
		{
			// Complex root; find its conjugate for stable real deflation
			double r2Re = 0.0, r2Im = 0.0;
			LaguerreZ(r2Re, r2Im, rRe, -rIm, n, P);

			// Average conjugate pair to stabilize
			double ReAvg = 0.5 * (rRe + r2Re);
			double ImAvg = 0.5 * (rIm - r2Im); // since r2Im ~ -rIm

			roots.emplace_back(ReAvg, ImAvg);
			roots.emplace_back(ReAvg, -ImAvg);

			// Quadratic factor: x^2 - 2*ReAvg*x + (ReAvg^2 + ImAvg^2)
			double a = -2.0 * ReAvg;
			double b = ReAvg * ReAvg + ImAvg * ImAvg;
			auto S = deflateQuadratic(P, n, a, b);
			P.assign(S.begin(), S.end());
			n -= 2;
		}
	}

	std::cout.setf(std::ios::fixed); 
	std::cout.precision(12);
	std::cout << "Miejsca zerowe (Laguerre):\n";
	for (size_t i = 0; i < roots.size(); ++i)
	{
		double re = roots[i].first;
		double im = roots[i].second;
		if (std::fabs(im) <= tolIm)
			std::cout << "x_" << (i + 1) << " = " << re << "\n";
		else
			std::cout << "x_" << (i + 1) << " = " << re << (im >= 0 ? " + " : " - ") << std::fabs(im) << " i" << "\n";
	}

	return 0;
}

