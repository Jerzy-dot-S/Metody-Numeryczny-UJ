N06 — Roots of a degree-7 polynomial via Laguerre's method

Polynomial:
  243 x^7 − 486 x^6 + 783 x^5 − 990 x^4 + 558 x^3 − 28 x^2 − 72 x + 16 = 0

Implementation:
- Uses Laguerre's method in complex arithmetic to find one root at a time
- Deflates the polynomial after each root using synthetic division (ascending coefficients)
- Prints roots and residuals |P(r)| to confirm accuracy

Build (Windows, PowerShell):
  g++ -O2 -std=c++17 -o N06/N06.exe N06/N06.cpp

Run:
  .\N06\N06.exe

Notes:
- Coefficients are stored in ascending order [a0, a1, ..., a7]
- Small imaginary parts (< 1e-10) are printed as 0 for readability
- If any root seems inaccurate, rerun with different initial seeds or increase max iterations
