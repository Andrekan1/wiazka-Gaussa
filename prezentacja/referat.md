# Wiązka Gaussa i mody wyższych rzędów: od wzoru do modelu w COMSOL

**Referat do prezentacji (ok. 15 min, dwóch prelegentów)**

- **Prelegent 1 – teoria:** slajdy 1–5 i 14, ok. 6,5 min.
- **Prelegent 2 – COMSOL:** slajdy 6–13, ok. 8 min.
- Kolejność slajdów, czasy i obrazy: [plan_prezentacji.md](plan_prezentacji.md).

**Źródła.** Korzystamy wyłącznie z zawartości tego repozytorium (README, modele `.mph`, obrazy w `zdj/` i `laguerre/`) oraz z samego COMSOL-a 6.4: dokumentacji modułu Wave Optics, przykładu z biblioteki aplikacji i własnego eksperymentu (slajd 11). Pełna lista jest w aneksie F.

Tekst pod nagłówkami slajdów to wypowiedź prelegenta. Aneksy na końcu nie są do czytania. To materiał pomocniczy: tabele, liczby i instrukcje odtworzenia.

---

## Część I – Teoria (Prelegent 1)

### Slajd 1. Wstęp *(0:00–0:30)*

Dzień dobry. Opowiemy o projekcie z laboratorium 4 „Propagacja wiązek światła strukturyzowanego”. Zbudowaliśmy w COMSOL-u trójwymiarowe modele wiązki gaussowskiej wyższego rzędu w dwóch wersjach: Hermite'a-Gaussa i Laguerre'a-Gaussa. Najpierw wyjaśnię, co oznaczają wzory, na których opierają się modele. Potem kolega pokaże, jak te wzory trafiają do programu, które parametry można zmieniać i jak sam program może zafałszować wynik przez tzw. dyspersję numeryczną.

### Slajd 2. Wiązka Gaussa: wzór i znaczenie symboli *(0:30–2:00)*

Wiązka laserowa nie jest falą płaską: jej energia jest skupiona wokół osi. Najprostszy opis takiej wiązki to wiązka Gaussa, rozwiązanie równania falowego w przybliżeniu przyosiowym:

```math
E(r,z) = E_0\,\frac{w_0}{w(z)}\,\exp\!\left(-\frac{r^2}{w^2(z)}\right)
\exp\!\left(-ikz - ik\frac{r^2}{2R(z)} + i\zeta(z)\right)
```

Każdy symbol ma swoją rolę:

- $`E_0`$ to amplituda. Równania są liniowe, więc zmiana $`E_0`$ tylko skaluje wynik.
- $`w_0`$ to promień przewężenia. W odległości $`w_0`$ od osi pole spada do 1/e, a natężenie do 1/e².
- $`w(z)=w_0\sqrt{1+(z/z_0)^2}`$ to promień wiązki w odległości z od przewężenia.
- $`z_0=\pi w_0^2/\lambda`$ to zasięg Rayleigha, czyli odległość, na której wiązka poszerza się √2 razy. Zależy od kwadratu $`w_0`$: wiązka dwa razy szersza „trzyma się” cztery razy dłużej.
- $`R(z)`$ to promień krzywizny frontu falowego. W przewężeniu jest nieskończony, czyli front jest płaski.
- $`\zeta(z)=\arctan(z/z_0)`$ to faza Gouya: dodatkowe przesunięcie fazy względem fali płaskiej. Przy przejściu przez całe ognisko wynosi łącznie π.
- $`k = 2\pi n/\lambda`$ to liczba falowa w ośrodku o współczynniku załamania n.

Ten sam wzór, z fazą Gouya, podaje dokumentacja COMSOL-a. Ostrzega też, że jest to wzór przyosiowy, który traci dokładność, gdy promień wiązki zbliża się do długości fali. U nas $`w_0 = 10\lambda`$, więc przybliżenie jest dobre.

### Slajd 3. Skale w naszym modelu *(2:00–3:00)*

Wstawmy liczby z modeli: λ = 1 µm i $`w_0`$ = 10 µm. Zasięg Rayleigha wynosi wtedy 314 µm, a kąt rozbieżności 32 mrad. Obszar obliczeniowy ma jednak tylko h = 3 µm długości, około sto razy mniej niż $`z_0`$. Na końcu obszaru wiązka jest szersza zaledwie o 0,005%, a faza Gouya wynosi 9,5 mrad, czyli około pół stopnia.

Wynikają z tego dwa wnioski. Po pierwsze, w naszych modelach widać głównie poprzeczną strukturę modów i to, że zachowuje się ona w trakcie propagacji. Rozmywania dyfrakcyjnego nie widać. Po drugie, wejście modelu leży dokładnie w przewężeniu. Tam $`w = w_0`$, $`R`$ jest nieskończone, a $`\zeta = 0`$, więc z całego wzoru zostaje sam profil poprzeczny. To jedyna rzecz, którą podajemy programowi. Rozszerzanie się wiązki, krzywiznę frontu i fazę Gouya COMSOL musi „wyliczyć sam”.

### Slajd 4. Mody Hermite'a-Gaussa *(3:00–4:15)*

Mody wyższych rzędów powstają przez pomnożenie profilu Gaussa przez wielomiany. W modach Hermite'a-Gaussa są to wielomiany Hermite'a, osobno w x i w y:

```math
E_{l,m} = E_0\,\frac{w_0}{w}\,
H_l\!\left(\tfrac{\sqrt{2}\,x}{w}\right)H_m\!\left(\tfrac{\sqrt{2}\,y}{w}\right)
e^{-(x^2+y^2)/w^2}\,
e^{-ikz - ik\frac{x^2+y^2}{2R} + i(l+m+1)\zeta}
```

Wielomian $`H_l`$ ma l miejsc zerowych, więc mod HG(l, m) ma l linii węzłowych w poprzek osi x i m w poprzek osi y. Daje to **(l+1)·(m+1) plam**, a sąsiednie plamy mają przeciwne znaki pola. Faza Gouya rośnie z rzędem jako (l+m+1)ζ: dla HG(4,4) jest dziewięć razy większa niż dla wiązki podstawowej. W COMSOL-u wielomiany $`H_0`$…$`H_4`$ siedzą w funkcjach analitycznych `Hx` i `Hy`, które wybierają wielomian zależnie od parametrów `lx` i `my`.

Jedna uwaga praktyczna: wzór nie jest unormowany. Mody wyższe mają więc w modelu dużo większą amplitudę: od ok. 1 V/m dla HG(0,0) do ok. 220 V/m dla HG(4,4). Żeby wszystkie mody niosły tę samą moc, trzeba podzielić amplitudę przez $`\sqrt{2^{l+m}\,l!\,m!}`$.

### Slajd 5. Mody Laguerre'a-Gaussa, moment pędu i interferogramy *(4:15–6:00)*

Druga rodzina to mody Laguerre'a-Gaussa, zapisane we współrzędnych walcowych:

```math
E_{p,m} = E_0\,\frac{w_0}{w}\left(\tfrac{\sqrt{2}\,r}{w}\right)^{|m|}
L_p^{|m|}\!\left(\tfrac{2r^2}{w^2}\right)e^{-r^2/w^2}\,
e^{-ikz - ik\frac{r^2}{2R} + im\varphi + i(2p+|m|+1)\zeta}
```

Indeks radialny p mówi, ile jest ciemnych pierścieni: mod ma **p+1 jasnych pierścieni**. Indeks azymutalny m to ładunek topologiczny. Dla m ≠ 0 pole przy osi rośnie jak $`r^{|m|}`$, więc na osi jest ciemny punkt, w którym faza jest nieokreślona. To jest **wir optyczny**.

Najciekawszy jest czynnik $`e^{im\varphi}`$. Na jednym obiegu wokół osi faza zmienia się o m·2π, więc front falowy nie jest płaski, tylko skręcony w m splecionych helis, jak korkociąg. Taka wiązka niesie **orbitalny moment pędu m·ħ na foton**, a w modelu możemy go policzyć całką:

```math
\frac{\langle L_z\rangle}{\hbar}=\frac{\int E^*\left(-i\,\partial_\varphi\right)E\,dA}{\int|E|^2\,dA}.
```

Fazy nie da się zobaczyć kamerą, dlatego w eksperymencie nakłada się wiązkę wirową na wiązkę odniesienia. Z pochyloną falą płaską powstają proste prążki, które przy osi tworzą **widelec**: po jednej stronie osi jest o |m| prążków więcej. Ze współosiową falą sferyczną powstaje **spirala o |m| ramionach**. Oba obrazy pozwalają odczytać m bez żadnych obliczeń.

Dwie uwagi praktyczne. W COMSOL-u kąt azymutalny liczymy jako `atan2(y,x)` i kolejność argumentów ma znaczenie. Poza tym COMSOL używa konwencji czasowej $`e^{+i\omega t}`$, więc fizycznie nasz mod z $`e^{+im\varphi}`$ ma moment $`L_z=-m\hbar`$ względem osi +z. Liczba +m, którą policzymy, to ładunek topologiczny zespolonej amplitudy.

Przekazuję głos koledze, który pokaże, jak to wszystko liczy COMSOL.

---

## Część II – COMSOL (Prelegent 2)

### Slajd 6. Jak COMSOL liczy wiązkę: metoda obwiedni *(6:00–7:15)*

Najprostszym podejściem byłoby rozwiązać równania Maxwella bezpośrednio dla pola E. Dokumentacja COMSOL-a mówi jednak, że do dobrego opisu fali potrzeba około pięciu elementów drugiego rzędu na każdą długość fali. Nasz walec ma średnicę 100 długości fali, więc taka siatka szybko robi się ogromna.

Dlatego używamy interfejsu **Electromagnetic Waves, Beam Envelopes** (w skrócie ewbe). Pole zapisuje się w nim jako iloczyn obwiedni i znanej z góry fazy:

```math
\mathbf{E}(\mathbf{r}) = \mathbf{E}_1(\mathbf{r})\,e^{-i\,\mathbf{k}_1\cdot\mathbf{r}},
\qquad
(\nabla - i\mathbf{k}_1)\times\big[(\nabla - i\mathbf{k}_1)\times\mathbf{E}_1\big] - k^2\mathbf{E}_1 = 0,\quad k = n k_0 .
```

To nadal dokładne równanie Maxwella, bez przybliżenia przyosiowego. Niewiadomą jest jednak tylko obwiednia $`\mathbf{E}_1`$. Jeśli wektor $`\mathbf{k}_1`$ dobrze opisuje prawdziwą falę, obwiednia zmienia się wolno i elementy wzdłuż z mogą być dużo dłuższe niż λ. W obu modelach wybieramy wariant jednokierunkowy i $`\mathbf{k}_1 = (0, 0, k_0)`$, gdzie $`k_0`$ to w COMSOL-u zmienna `ewbe.k0`.

### Slajd 7. Budowa modelu krok po kroku *(7:15–8:30)*

Na schemacie widać cały model, a w drzewie COMSOL-a to dziewięć kroków:

1. **Parameters**: `lam`, `w0` równe dziesięciu długościom fali, `E0`, indeksy modu, liczby elementów siatki i współczynnik załamania.
2. **Functions**: funkcje analityczne z wielomianami `Hx`, `Hy` albo `L`.
3. **Variables**: współrzędne bezwymiarowe `xp`, `rp` i kąt `phi`.
4. **Geometry**: walec o promieniu 50 µm i długości 3 µm.
5. **Materials**: współczynnik załamania.
6. **Physics** (ewbe): na wejściu *Matched Boundary Condition* z polem padającym, czyli tu wpisujemy profil modu (wyrażenie widać na slajdzie). Na wyjściu ten sam warunek bez pola padającego, z boku PEC.
7. **Mesh**: siatka zamiatana wzdłuż z.
8. **Study**: *Wavelength Domain* plus *Parametric Sweep* po indeksach modu.
9. **Results**: przekroje, izopowierzchnie, animacje i całki.

Najważniejsze: narzucamy tylko profil w z = 0, a resztę liczy program. Sprawdziliśmy to na modzie HG(2,3). Faza Gouya odczytana z rozwiązania wyniosła 0,05737 rad, a teoria daje 0,05729 rad. Różnica to 0,13%, chociaż tej fazy nigdzie nie wpisaliśmy.

### Slajd 8. Parametry i ich wpływ *(8:30–9:15)*

Tabela na slajdzie zbiera parametry, które warto zmieniać:

- `lx`, `my` (HG) oraz `m`, `l` = p (LG) zmieniają kształt modu: liczbę plam, pierścieni i ramion spirali.
- `w0` zmienia szerokość wiązki i, z kwadratem, zasięg Rayleigha.
- `lam` zmienia długość fali. W modelu LG wymiary są podane w długościach fali, więc zmiana `lam` tylko przeskalowuje cały obraz. W modelu HG długość obszaru jest stała (3 µm), więc zmienia się stosunek h/z₀.
- `n_real` i `n_im` zmieniają materiał. Tu trzeba uważać na wektor $`\mathbf{k}_1`$, o czym za chwilę.
- `Nl` i rozmiar siatki decydują o dokładności i o czasie obliczeń.

### Slajd 9. Wyniki HG *(9:15–10:05)*

W animacji wszystkich 25 modów zawsze widać (l+1)·(m+1) plam o naprzemiennych znakach, tak jak przewiduje teoria. Profil wzdłuż trzech mikrometrów się nie zmienia, bo jesteśmy daleko przed zasięgiem Rayleigha.

Kolega zbadał też siatkę na modzie HG(3,3). Na najgrubszych ustawieniach plamy są kanciaste i zniekształcone. Od ustawienia *Fine* obraz się stabilizuje. Wzdłuż osi przy jednym elemencie widać tylko gładkie „rury” bez żadnej struktury. Od dwóch-trzech elementów pojawia się struktura w kierunku propagacji, a od około siedmiu kolejne obrazy prawie się nie różnią.

### Slajd 10. Wyniki LG *(10:05–11:15)*

Dla m = 3 i p = 0…3 widać kolejno 1, 2, 3 i 4 jasne pierścienie oraz ciemny wir na osi. Interferogram z pochyloną falą daje widelec, w którym trzy prążki kończą się na osi. Interferogram z falą sferyczną daje spiralę o trzech ramionach. Powierzchnia stałej fazy to trzy splecione helisy, które w animacji obracają się wraz z upływem okresu optycznego.

Na koniec liczymy moment pędu całką po przekroju, za pomocą operatorów całkowania `intop_in` i `intop_out`. Wychodzi od 2,97 do 2,93 zamiast 3. Ta sama całka dla wzoru analitycznego daje dokładnie 3,0000, więc różnica nie pochodzi z całkowania, tylko z dyskretyzacji pola na siatce. Po 3 µm propagacji wartość spada jeszcze o 0,3–0,6%, a moc zachowuje się z dokładnością do 0,02%. Fizycznie moment pędu w jednorodnym ośrodku jest stały, więc ten spadek to błąd numeryczny. Prowadzi nas to do dyspersji numerycznej.

### Slajd 11. Dyspersja numeryczna: eksperyment w COMSOL-u *(11:15–12:45)*

Dyspersja numeryczna polega na tym, że fala na siatce biegnie z nieco inną liczbą falową niż prawdziwa. Błąd fazy narasta z każdą długością fali, a zbyt rzadka siatka zaniża też amplitudę. W metodzie obwiedni siatka musi odwzorować tylko obwiednię. Gdy $`k_1 = n k_0`$, obwiednia jest prawie stała. Gdy nie, oscyluje jak $`e^{-i(nk_0-k_1)z}`$ i trzeba ją rozdzielić siatką.

Sprawdziliśmy to na kopii modelu HG. Wzięliśmy wiązkę podstawową, n = 1 albo 1,5 i $`k_1`$ równe 0, $`k_0`$ albo $`1{,}5\,k_0`$, oraz od 1 do 8 elementów wzdłuż z. Na wyjściu porównaliśmy fazę na osi z teorią.

- **Gdy $`k_1 = n k_0`$**, błąd fazy wynosi 0,0007 rad niezależnie od liczby elementów, nawet przy jednym.
- **Gdy $`k_1`$ jest niedopasowane**, błąd maleje z siatką od 3 rad do 0,012 rad przy ośmiu elementach. To około pięciu elementów na „długość fali” obwiedni, czyli dokładnie tyle, ile zaleca dokumentacja.
- **Gdy $`k_1 = 0`$**, czyli liczymy pełne pole, osiem elementów to mniej niż trzy na λ. Fala ma złą długość i amplitudę 0,6 zamiast 1. To jest dyspersja numeryczna w czystej postaci.

Zaskoczyło nas jeszcze jedno: przy złym $`k_1`$ amplituda dąży do $`2k_1/(k_1+nk_0)`$, czyli 0,8 albo 1,2 zamiast 1, i zagęszczanie siatki tego nie naprawia. Dlatego serię z różnymi n traktujemy jakościowo, a poprawka to `k1 = n_real*ewbe.k0`.

### Slajd 12. Zasoby i powtarzalność *(12:45–13:25)*

Wszystko zrobiliśmy w jednym module, Wave Optics, z interfejsem ewbe. Pełna lista węzłów COMSOL-a użytych w każdym zadaniu jest w aneksie B referatu. Najcięższy jest model LG: prawie milion stopni swobody, 2 min 45 s obliczeń na cztery mody i około 24 GB pamięci. Model HG liczy 25 modów w półtorej minuty.

Żeby projekt dało się powtórzyć, w repozytorium są oba modele i opis każdego kroku. Model HG jest zapisany bez rozwiązań: wystarczy kliknąć *Compute*. Model LG ma rozwiązania, dlatego waży 98 MB, tuż pod limitem GitHuba. Opisaliśmy też pułapki, na które trafiliśmy: siatkę sterowaną fizyką, która nie przebudowywała się w sweepie, i eksport danych, który wymagał układu wierszowego.

### Slajd 13. Co jeszcze można zmienić *(13:25–14:00)*

Kilka pomysłów, które łatwo sprawdzić w tych modelach:

- wydłużyć obszar do $`z_0`$, żeby zobaczyć poszerzenie √2 i fazę Gouya π/4·(N+1);
- zmniejszyć $`w_0`$ do kilku λ i zobaczyć, gdzie przestaje działać przybliżenie przyosiowe;
- złożyć HG(1,0) + i·HG(0,1) i dostać wir LG z m = 1;
- odwrócić znak m i zobaczyć, że widelec i spirala się odwracają;
- dodać polaryzację kołową.

Szczegóły są w aneksie G.

### Slajd 14. Podsumowanie *(Prelegent 1, 14:00–14:30)*

Podsumowując:

- wzory na mody Hermite'a-Gaussa i Laguerre'a-Gaussa wystarczy podać COMSOL-owi tylko na wejściu, w przewężeniu, a resztę program liczy sam, łącznie z fazą Gouya;
- metoda obwiedni pozwala to zrobić małym kosztem, pod warunkiem że wektor falowy obwiedni pasuje do ośrodka, bo inaczej wraca dyspersja numeryczna i błąd warunku brzegowego;
- modele LG pokazują wir optyczny, jego interferogramy i moment pędu równy trzem jednostkom ħ na foton.

Dziękujemy, chętnie odpowiemy na pytania.

---

## Aneksy (materiał pomocniczy, nie do czytania)

### A. Symbole ze wzorów i ich odpowiedniki w COMSOL

| Symbol | Znaczenie | W COMSOL | Wartość w modelach | Na co wpływa |
|---|---|---|---|---|
| λ | długość fali w próżni | parametr `lam`; badanie *Wavelength Domain* | 1 µm | skala całego problemu; `ewbe.k0` = 2π/λ |
| $`k_0`$, k | liczba falowa w próżni / w ośrodku | `ewbe.k0`; k = n·k₀ liczy fizyka | 6,28·10⁶ rad/m | szybkość oscylacji fazy |
| $`\mathbf{k}_1`$ | wektor falowy obwiedni | ewbe → *Wave Vectors*: `(0,0,ewbe.k0)` | (0, 0, k₀) | co musi rozdzielić siatka (slajd 11) |
| $`w_0`$ | promień przewężenia | parametr `w0` = 10·`lam` | 10 µm | szerokość wiązki, $`z_0 \propto w_0^2`$ |
| $`z_0`$ | zasięg Rayleigha | zmienna `z0 = k*w0^2/2` (HG) | 314 µm | jak szybko wiązka się rozszerza |
| w(z), R(z), ζ(z) | promień, krzywizna, faza Gouya | `w` (HG); R i ζ nie są wpisane | w(0) = w₀ | **wynikają z rozwiązania** |
| $`E_0`$ | amplituda | parametr `E0` | 1 V/m | tylko skala (problem liniowy) |
| √2x/w₀, √2r/w₀ | współrzędne bezwymiarowe | zmienne `xp`, `yp`, `rp` | – | argumenty wielomianów |
| φ | kąt azymutalny | `phi = atan2(y,x)` | – | kierunek skręcenia wiru |
| $`H_l`$, $`H_m`$ | wielomiany Hermite'a | funkcje `Hx(xp)`, `Hy(yp)`, wybór przez `lx`, `my` | rzędy 0…4 | liczba plam (l+1)(m+1) |
| $`L_p^{\lvert m\rvert}`$ | uogólnione wielomiany Laguerre'a | funkcja `L(rp^2, abs(m))`, wybór przez `l` = p | p = 0…3 | liczba pierścieni p+1 |
| m | ładunek topologiczny | parametr `m`, czynnik `exp(i*m*phi)` | 3 | OAM, widelec, spirala |
| n, κ | współczynnik załamania i ekstynkcji | *Materials* → *Refractive index*: `n_real`, `n_im` | 1 i 0 | długość fali w ośrodku i tłumienie |
| E, $`\mathbf{E}_1`$ | pole całkowite / obwiednia | `ewbe.Ex`, `ewbe.normE` / `E1x` | – | wielkości do wykresów i całek |
| h, R_dom | długość i promień obszaru | *Geometry* → *Cylinder* | 3 µm, 50 µm | ile propagacji widać, koszt |
| N_T, N_L | elementy poprzeczne / podłużne | *Mesh* (ewbe): `Nt`, `Nl` | Nl = 5 (HG), 10 (LG) | dokładność, pamięć, czas |

### B. Zasoby COMSOL użyte w poszczególnych zadaniach

| Zadanie | Węzły i funkcje COMSOL |
|---|---|
| Model podstawowy HG (`wzf_l04_01_hermite_gauss.mph`) | Model Wizard 3D; *Wave Optics → Electromagnetic Waves, Beam Envelopes* (jednokierunkowo, wektor falowy k₁); *Parameters*; *Analytic* (`Hx`, `Hy` z `if()`); *Variables*; *Cylinder*; *Material* z węzłem *Refractive index*; *Matched Boundary Condition* ×2 (z polem padającym / bez); *Perfect Electric Conductor*; siatka sterowana fizyką (*Swept*, *Extremely fine*); *Wavelength Domain*; *Parametric Sweep* |
| Wszystkie rzędy HG (`zdj/orders/`) | *Parametric Sweep* po `lx`, `my` ∈ {0…4} (*All combinations*, 25 rozwiązań); wykresy 3D modułu pola i Ex; *Export → Animation* (GIF po parametrze) |
| Rozmiar siatki (`zdj/mesh_topology/`) | *Mesh → Element size*: 9 predefiniowanych rozmiarów od *Extremely coarse* do *Extremely fine*; *Export → Image* |
| Elementy wzdłuż z (`zdj/mesh_sideways/`) | siatka sterowana fizyką ewbe: liczba elementów poprzecznych N_T = 50, podłużnych N_L = `Nl` = 1…10 |
| Ośrodek stratny (`zdj/n_real/`) | *Material → Refractive index* (n = `n_real`, κ = `n_im`); sweep n przy stałym κ; animacje GIF |
| Ośrodek anizotropowy (`zdj/orders/anisotropy/`) | *Refractive index* jako tensor (wartości nie są zapisane w pliku w repo); sweep `lx`, `my`; animacje GIF |
| Model LG (`wzf_l04_02_laguerre_gauss.mph`) | jak HG, plus funkcja *Analytic* `L(x,m)`; zmienne `r`, `phi = atan2(y,x)`, `rp`; *Parametric Sweep* po `l` = p ∈ {0…3}; wykresy *Slice* (moduł pola, faza `arg(ewbe.Ex)` z tablicą kolorów *Cyclic*); *Animation* (Player + GIF) |
| Prezentacja LG (`laguerre/lg_oam/`) | *Slice* z wyrażeniem interferogramu (fala odniesienia jako wyrażenie analityczne), *Filter*, tablica *GrayScale*; *Isosurface* okna fazy z filtrem i *Color Expression*; animacja z *Dynamic data extension* (obrót fazy); *Definitions → Integration* (`intop_in`, `intop_out`); *Derived Values → Global Evaluation*; *Table → Export* (CSV); *Study → Update Solution* (nowe operatory bez ponownego liczenia) |
| Eksperyment z dyspersją (`prezentacja/`) | wektor falowy ewbe jako wyrażenie `kfac*ewbe.k0`; siatka przełączona na *User-controlled* (*Swept* + *Distribution* = `Nl`); *Parametric Sweep* (`n_real`, `kfac`, `Nl`); *Datasets → Cut Point 3D* i *Cut Line 3D*; *Derived Values → Point Evaluation*; *Export → Table* i *Export → Data* (siatka punktów, układ wierszowy) |

### C. Zasoby obliczeniowe

Wszystko liczyliśmy na komputerze z procesorem 6-rdzeniowym (Ryzen 5 3600) i 32 GB RAM, w COMSOL 6.4.

| Model / badanie | Siatka | Stopnie swobody | Czas | Pamięć |
|---|---|---|---|---|
| HG, 25 modów (zapisane w pliku) | *Extremely fine*, Nl = 5: 6786 trójkątów w przekroju, 33 930 pryzm | nie zapisano w logu | 1 min 35 s | – |
| LG, 4 mody (p = 0…3) | *Extremely fine*, N_T = N_L = 10 | 992 682 | 2 min 45 s | 23,8 GB (fizyczna), 31,7 GB (wirtualna) |
| Eksperyment z dyspersją, 24 przypadki | *Normal* (1152 trójkąty w przekroju), Nl = 1…8 | 22 250 / 38 668 / 71 504 / 137 176 dla Nl = 1 / 2 / 4 / 8 | 1 min 42 s | – |

Pamięć rośnie szybko z liczbą warstw siatki. Model LG z 10 warstwami zajął 24 z 32 GB. Dlatego w eksperymencie zmniejszyliśmy siatkę poprzeczną do *Normal*. Wpływa to tylko na stały „próg” błędu (0,0007 rad), a nie na zależność od Nl.

### D. Jak odtworzyć modele

1. Potrzebny jest COMSOL Multiphysics 6.4 z modułem Wave Optics.
2. **HG:** otwórz `wzf_l04_01_hermite_gauss.mph` i uruchom *Study 1 → Compute*. Liczy się sweep `lx`, `my` = 0…4. Zmiana listy w *Parametric Sweep* daje inne mody.
3. **LG:** otwórz `wzf_l04_02_laguerre_gauss.mph`. Rozwiązania są w pliku, więc wykresy, animacje i *Derived Values* działają od razu. Po dodaniu nowych operatorów w *Definitions* trzeba wykonać *Study → Update Solution*. Ponowne liczenie wymaga ok. 24 GB RAM. Na słabszym komputerze zmniejsz rozmiar siatki albo liczbę warstw.
4. **Budowa od zera (skrót):**
   1. *Model Wizard* → 3D → *Wave Optics → Electromagnetic Waves, Beam Envelopes* → *Wavelength Domain*.
   2. *Parameters* jak w aneksie A.
   3. Funkcje *Analytic* z wielomianami.
   4. Zmienne `r`, `xp`, `yp`, `rp`, `phi`.
   5. *Cylinder* (oś z).
   6. *Material* z *Refractive index*.
   7. W ewbe: *Unidirectional*, k₁ = `(0,0,ewbe.k0)`.
   8. *Matched Boundary Condition* na z = 0 z polem padającym (wyrażenie profilu jako składowa x).
   9. *Matched Boundary Condition* na z = h bez pola padającego.
   10. Siatka *Swept*.
   11. *Parametric Sweep* po indeksach modu.
5. **Czego nie da się odtworzyć z plików:** ustawień serii `zdj/mesh_sideways`, `zdj/n_real` i `zdj/orders/anisotropy`. Liczono je na zmodyfikowanych wersjach modelu, które nie trafiły do repozytorium. Zostały tylko obrazy. Na przyszłość warto zapisywać każdą serię jako osobne *Study* w tym samym pliku.
6. **Pułapki, na które trafiliśmy:**
   - siatka sterowana fizyką nie przebudowała się w sweepie po `Nl`, więc wszystkie przypadki liczyły się na tej samej siatce, aż przełączyliśmy ją na *User-controlled* z *Distribution* = `Nl`;
   - *Export → Data* dla rozwiązań na różnych siatkach wymaga ustawienia *Parameter layout: Rows*;
   - nowe operatory całkowania po obliczeniu wymagają *Update Solution*;
   - GitHub odrzuca pliki powyżej 100 MB, a plik LG z rozwiązaniami ma 98 MB.

### E. Eksperyment z dyspersją numeryczną: ustawienia i dane

**Plik:** kopia modelu HG, `COMSOL/dyspersja/hg_dyspersja.mph`, poza repozytorium (79,5 MB z rozwiązaniami). Oryginał w repo nie był zmieniany.

**Zmiany względem `wzf_l04_01_hermite_gauss.mph`:**

1. *Parameters*: `lx` = `my` = 0 (wiązka podstawowa) i nowy parametr `kfac` (k₁ = `kfac`·k₀).
2. ewbe → *Wave Vectors*: składowa z wektora k₁ = `kfac*ewbe.k0`.
3. *Mesh*:
   - rozmiar *Normal*;
   - typ sekwencji *User-controlled*;
   - *Swept 1 → Distribution 1*: liczba elementów = `Nl`.
4. *Parametric Sweep* (*All combinations*):
   - `n_real` ∈ {1; 1,5};
   - `kfac` ∈ {0; 1; 1,5};
   - `Nl` ∈ {1, 2, 4, 8}.
5. *Cut Point 3D* w (0, 0, 3 µm) i *Point Evaluation* z dwoma wyrażeniami:
   - błąd fazy `arg(ewbe.Ex*exp(i*n_real*ewbe.k0*z))`, od którego odejmujemy teoretyczną fazę Gouya $`\arctan(h/z_0)`$, gdzie $`z_0=\pi n w_0^2/\lambda`$ (9,55 mrad dla n = 1 i 6,37 mrad dla n = 1,5);
   - amplituda `abs(ewbe.Ex)/E0`.
6. *Cut Line 3D* od (0,0,0) do (0,0,3 µm) i *Export → Data* dla `real(ewbe.Ex)` oraz `abs(ewbe.Ex)`, w 301 punktach.

**Dane i wykresy w repo:**
- dane: [dane/dyspersja_numeryczna.csv](dane/dyspersja_numeryczna.csv) i [dane/pole_na_osi.txt](dane/pole_na_osi.txt);
- skrypt rysujący: [dane/wykresy_dyspersja.py](dane/wykresy_dyspersja.py);
- wykresy: `img/dyspersja_*.png`.

| Przypadek | „Długość fali” obwiedni λ/\|n − k₁/k₀\| | \|Δφ\| przy Nl = 1 / 2 / 4 / 8 [rad] | Amplituda przy Nl = 1 / 2 / 4 / 8 |
|---|---|---|---|
| n = 1, k₁ = k₀ (dopasowany) | ∞ | 0,00074 / 0,00075 / 0,00074 / 0,00074 | 0,997 / 0,997 / 0,997 / 0,997 |
| n = 1,5, k₁ = 1,5·k₀ (dopasowany) | ∞ | 0,00074 / 0,00074 / 0,00074 / 0,00074 | 0,997 / 0,997 / 0,997 / 0,997 |
| n = 1,5, k₁ = k₀ | 2 µm | 3,04 / 0,259 / 0,073 / 0,0116 | 0,315 / 0,432 / 0,788 / 0,820 |
| n = 1, k₁ = 1,5·k₀ | 2 µm | 2,77 / 0,779 / 0,022 / 0,016 | 0,371 / 0,397 / 1,218 / 1,242 |
| n = 1, k₁ = 0 (pełne pole) | 1 µm | 1,10 / 0,57 / 0,56 / 1,20 | 0,148 / 0,084 / 0,045 / 0,600 |
| n = 1,5, k₁ = 0 (pełne pole) | 0,67 µm | 1,88 / 2,19 / 1,20 / 2,97 | 0,102 / 0,042 / 0,487 / 0,104 |

**Interpretacja.** Błąd fazy przy dopasowanym k₁ nie zależy od Nl. Stała resztka 0,0007 rad i amplituda 0,997 zamiast ≈ 1 pochodzą z siatki poprzecznej *Normal*.

Przy niedopasowanym k₁ amplituda zbiega do $`2k_1/(k_1+nk_0)`$, czyli 0,8 i 1,2, z dokładnością ok. 3%. To nasza interpretacja: warunek *Matched Boundary Condition* zachowuje się tak, jakby był dopasowany do fali o liczbie falowej k₁. Nie znaleźliśmy tego wprost w dokumentacji.

### F. Źródła (tylko repozytorium i COMSOL)

1. README tego repozytorium: wstęp teoretyczny, opis modeli i wyników, wzory HG i LG.
2. Modele `wzf_l04_01_hermite_gauss.mph` i `wzf_l04_02_laguerre_gauss.mph`: ustawienia, logi solvera (stopnie swobody, czas, pamięć).
3. Obrazy w `zdj/` i `laguerre/` (wyniki zespołu).
4. COMSOL Multiphysics 6.4, *Wave Optics Module User's Guide*:
   - s. 29–30: wymagana rozdzielczość siatki, ok. 5 elementów 2. rzędu na λ, i wyjątek dla wydzielonej szybkiej fazy;
   - s. 161: *Matched Boundary Condition*;
   - s. 211–214: interfejs ewbe, sformułowanie jednokierunkowe i wektor k₁, siatka sterowana fizyką;
   - s. 219–220: równanie dla obwiedni;
   - s. 256–257: wiązka Gaussa we wzorach COMSOL-a, z fazą Gouya, i granice przybliżenia przyosiowego.
5. COMSOL Application Library, *Orbital Angular Momentum Beam* (`Wave_Optics_Module/Beam_Propagation/orbital_angular_momentum`): wzór LG i sposób zadania wiązki w ewbe.
6. Własny eksperyment w COMSOL 6.4 (aneks E).

### G. Co jeszcze można zmienić w modelach

| Pomysł | Co zmienić w COMSOL | Czego się spodziewać (z teorii) |
|---|---|---|
| Zobaczyć dyfrakcję i pełną fazę Gouya | *Cylinder*: h = z₀ ≈ 314 µm. Przy k₁ = k₀ obwiednia zmienia się w skali z₀, więc wystarczy kilkanaście–kilkadziesiąt warstw siatki | w(z₀) = √2·w₀; faza Gouya (N+1)·π/4, gdzie N = l+m albo 2p+\|m\| |
| Ognisko wewnątrz obszaru | pełny wzór z w(z), R(z), ζ(z) na wejściu przesuniętym o −z_f; dla modu podstawowego można też użyć wbudowanej opcji *Gaussian beam* w *Matched Boundary Condition* (położenie ogniska p₀) | wiązka zwęża się do ogniska, potem rozszerza; skok fazy Gouya przy ognisku |
| Granica przybliżenia przyosiowego | `w0` = 2·`lam` albo mniej; gęstsza siatka poprzeczna | odstępstwa od wzorów przyosiowych (ostrzega przed tym dokumentacja COMSOL-a) |
| HG → LG | pole padające `HG10 + i·HG01`: zmienne z H₁(xp)·H₀(yp) i H₀(xp)·H₁(yp) | pierścień z wirem m = 1; całka OAM ≈ 1 |
| Odwrócony wir | `m` = −3 | widelec i spirala odwrócone; OAM ≈ −3 |
| Spinowy moment pędu | pole padające `(Ex, i*Ex, 0)/sqrt(2)` (polaryzacja kołowa) | wiązka z polaryzacją kołową (spin) i dodatkowo z wirem (orbitalny moment pędu) |
| Poprawna seria n, κ | ewbe: `k1 = (0,0,n_real*ewbe.k0)`; Nl tak, żeby na drogę zaniku δ = λ/(2πκ) wypadało kilka elementów | zanik amplitudy jak exp(−κk₀z) bez artefaktów brzegowych |
| Parametry interferogramów | kąt α, promień R_ref i amplituda odniesienia w wyrażeniu *Slice* | okres prążków λ/sin α; inna gęstość ramion spirali |
