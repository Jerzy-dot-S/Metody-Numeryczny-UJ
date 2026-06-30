Program znajduję wartości własne macierzy opisanej w zadaniu dwoma sposobami.
Niestety ciężko mi powiedzieć która metoda jest lepsza, przy wysokich iteracjach
z dokładnością sześciu liczb po przecinku wektory policzone innymi metodami mogą się
pokrywać, ale z tego co wiem przy odpowiednim wektorze początkowym metoda rayleigha jest 
szybsza i lepsza lecz nie znam dokładnych wartości macierzy dla porównania. 
Natomiast przy słabo ustawionym wektorze początkowym metoda rayleigha dąży do najbliżej 
wartości własnej a standardowa metoda potęgowa zawsze leci do najwyżej co do modułu.


Działanie programu:
- Dla wybranego N (ustawionego w main, domyślnie N = 500) program:
  - wypisuje największą i drugą co do wielkości wartość własną oszacowaną przez
    - metodę potęgową,
    - iterację Rayleigha (po starcie z metody potęgowej),
    aby można było je porównać;
  - wypisuje cztery najmniejsze (modułowo) wartości własne, wyznaczone przez
    - odwrotną metodę potęgową + iloraz Rayleigha (kolumna „power”),
    - odwrotną metodę potęgową + pełną iterację Rayleigha (kolumna „Rayleigh”).
- Dodatkowo drukowany jest krótki eksperyment porównujący zbieżność
  - metody potęgowej
  - i iteracji Rayleigha
  przy tym samym losowym wektorze startowym.
  Można zmienić ilość iteracji oraz wielkość macierzy.

Kompilacja (Windows, PowerShell):
  g++ -O2 -std=c++23 -o N03.exe N03.cpp

Uruchomienie:
  .\N03.exe


