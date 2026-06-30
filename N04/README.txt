N04 — Interpolacja Lagrange'a funkcji f(x) = 1/(1 + x^2)

Program oblicza interpolację Lagrange'a na przedziale [-5, 5] z użyciem dwóch
rodzajów węzłów:
(a) węzły równoodległe:   X_n = -5 + 10*n/(N-1),  n = 0..N
(b) węzły cosinusowe (typ Czebyszewa): X_n = 5*cos(n*pi/N), n = 0..N

Jak widać na wykresie pierwszy sposób wyboru wezłów jest dla tego przypadku gorszy.
Dla małego N czyli w tym zadaniu 7 interpolacja X_n = -5 + 10*n/(N-1) przyjmuje
najmniej błędne wartości i nadal jest daleko od prawdziwej funkcji. Dla dużych N 
z kolei zaczynają pojawiać się oscylacje Rungego w granicach przedziału które dla 
N = 100 osiągają wartość 1.36457868797e+15. Natomiast dla węzłów Czebyszewa bład zmniejsza
się, z tego co sprawdzałem to aż do N = 190, i interpolacja jest tak dokładna że przy 
rozdzielczości wykresu funkcje pokrywają się(Z tego powodu dobranie kolorów, można zmienić w plot_gnuplot.plt).
Na wykresie plot_errors.png jest porównianie błędów obu wyboru węzłów, gdzie 
węzły Czebyszewa są na tyle lepsze że widnieją jak prosta linia w tej rozdzielczości.
Dopiero w pliku errors.dat można zobaczyć że wartości tych błędów są bardzo małe.

Kompilacja (Windows, PowerShell):
  g++ -O2 -std=c++23 -o N04.exe N04.cpp

Pliki wyjściowe (w bieżącym katalogu):
- out_best.dat          — kolumny: x, f(x), P_equidistant(x), P_cosineSpec(x)
- errors.dat            — kolumny: N, max error (equidistant), max error (cosine spec)
- nodes_equidistant.dat — węzły i wartości f(x) dla najlepszego N (równoodległe)
- nodes_cosineSpec.dat  — węzły i wartości f(x) dla najlepszego N (cosinusowe)

Rysowanie wykresów:
- Program automatycznie próbuje uruchomić gnuplota i wygenerować pliki:
  - plot_best.png   — funkcja oryginalna oraz oba interpolanty + węzły
  - plot_errors.png — maksymalny błąd w zależności od N dla obu metod
- Jeśli gnuplot nie jest zainstalowany lub nie jest w PATH, można go uruchomić
  ręcznie z katalogu N04 poleceniem:
  gnuplot plot_gnuplot.plt
