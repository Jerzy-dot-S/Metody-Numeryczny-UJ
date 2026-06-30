Projekt Numerki - krótki opis i instrukcje
=========================================

Opis
----
Program zawiera trzy implementacje metod obliczania układów równań. Wzór shermana-morrisona jest najszybszą metodą lecz
działa tylko na dany typ macierzy a eliminacja gaussa i metoda sor są wolniejsze ale działają dla wszystkich macierzy.
Niestety czas oczekiwania na stworzenie wykresu większego niż 1000 jest za długi więc nie udało mi sie sprawdzić dalej ale
wzór shermana morrisona powinien być 10000000000000000 szybszy niż eliminacja gaussa

Kompilacja
----------

g++ -O2 -std=c++23 -Wall -Wextra -o N02.exe N02.cpp

```

Uruchomienie - pojedyncze N
---------------------------
- Uruchom pojedynczo i zapisz pełne rozwiązanie (N=1000):

.\N01.exe 1000 --save-solution solution_N1000.csv
.\N02.exe 1000 dump

Uruchomienie - tryb wsadowy (generacja danych)
---------------------------------------------
- Przykład generowania tabeli czasów (log‑rozłożone N):

.\N01.exe --batch 25 10 1000 times_N01.csv
.\N02.exe --batch 25 10 1000 times_N02.csv


Wyniki zostaną zapisane w `times_*.csv` z nagłówkiem np. N,dense_time,sherman_morrison_time,sor_time
Wartości -1 oznaczają, że dana metoda nie powiodła się i
skrypty gnuplot ignorują je jako brakujące dane.

Gnuplot — generowanie wykresów
--------------------------------
Po zapisaniu danych do pliku z wynikami można stworzyć wykres poprzez.

gnuplot plot_gnuplot.plt
