Projekt Numerki - krótki opis i instrukcje
=========================================

Opis
----
Program zawiera trzy implementacje metod obliczania układów równań 
Są to eliminacja gaussa, metoda thomasa i metoda sor. Metoda thomasa jest szybsza ale działa tylko dla macierzy trójdiagonalnych
Przy uruchamianiu programu przy samej metodzie thomasa udało mi się obliczyć macierz N = 100000000 w czasie 1.4 sekundy gdzie przy 
eliminacji gaussa prawdopodobnie czekałbym 10000000000000000 dłużej więc różnica czasu jest ogromna.

Kompilacja
----------

g++ -O2 -std=c++23 -Wall -Wextra -o N01.exe N01.cpp

```

Uruchomienie - pojedyncze N
---------------------------
- Uruchom pojedynczo i zapisz pełne rozwiązanie (N=1000):

.\N01.exe 1000 --save-solution solution_N1000.csv

Uruchomienie - tryb wsadowy (generacja danych)
---------------------------------------------
- Przykład generowania tabeli czasów (log‑rozłożone N):

.\N01.exe --batch 25 10 1000 times_N01.csv

Wyniki zostaną zapisane w `times_*.csv` z nagłówkiem np. N,dense_time,sherman_morrison_time,sor_time
Wartości -1 oznaczają, że dana metoda nie powiodła się i
skrypty gnuplot ignorują je jako brakujące dane.

Gnuplot — generowanie wykresów
--------------------------------
Po zapisaniu danych do pliku z wynikami można stworzyć wykres poprzez.

gnuplot plot_gnuplot.plt
