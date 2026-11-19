# Plan prezentacji: „Wiązka Gaussa i mody wyższych rzędów w COMSOL”

**Czas:** ok. 15 min, potem pytania.

**Prelegenci:**
- **P1 – teoria:** slajdy 1–5 i 14, ok. 6,5 min.
- **P2 – COMSOL:** slajdy 6–13, ok. 8 min.

Pełny tekst wypowiedzi jest w [referat.md](referat.md). Każdy slajd poniżej odsyła do odpowiedniego fragmentu referatu, a aneksy referatu są ściągą na pytania.

Ścieżki do obrazów są podane względem folderu `prezentacja/`. GIF-y działają w PowerPoincie i Google Slides. W PDF-ie trzeba wstawić pojedynczą klatkę.

## Podział czasu

| Blok | Slajdy | Prelegent | Czas | Narastająco |
|---|---|---|---|---|
| Wstęp i wiązka Gaussa | 1–3 | P1 | 3:00 | 3:00 |
| Mody HG i LG, OAM | 4–5 | P1 | 3:00 | 6:00 |
| Metoda obwiedni i budowa modelu | 6–8 | P2 | 3:15 | 9:15 |
| Wyniki HG i LG | 9–10 | P2 | 2:00 | 11:15 |
| Dyspersja numeryczna | 11 | P2 | 1:30 | 12:45 |
| Zasoby, powtarzalność, pomysły | 12–13 | P2 | 1:15 | 14:00 |
| Podsumowanie | 14 | P1 | 0:30 | 14:30 |

Zapas to ok. 30 s. Jeśli brakuje czasu, skracamy slajd 13 do jednego zdania i odsyłamy do aneksu G.

---

## Slajdy

### 1. Tytuł i cel *(P1, 0:00–0:30)*
- **Na slajdzie:**
  - tytuł, autorzy;
  - jedno zdanie celu: „Mody HG i LG w COMSOL: od wzoru do modelu i do jego ograniczeń numerycznych”.
- **Obraz:** `../laguerre/lg_sweep_l.gif` (pierścienie LG jako tło) albo `../laguerre/lg_oam/korkociag_faza_l0.png`.
- **Powiedzieć:** referat, slajd 1.

### 2. Wiązka Gaussa: wzór i symbole *(P1, 0:30–2:00)*
- **Na slajdzie:**
  - wzór na wiązkę Gaussa;
  - pod nim lista symboli: E₀, w₀, w(z), z₀, R(z), ζ(z), k, każdy z jednym zdaniem znaczenia;
  - ramka: „wzór przyosiowy – dobry, gdy w₀ ≫ λ (u nas w₀ = 10λ)”.
- **Obraz:** opcjonalnie pierwsza klatka `../zdj/orders/3x3_all_orders_hermite_enorm.gif`, czyli HG(0,0).
- **Powiedzieć:** referat, slajd 2. Mówić znaczenie symboli, nie czytać wzoru.

### 3. Skale w naszym modelu *(P1, 2:00–3:00)*
- **Na slajdzie:** tabelka liczb:
  - λ = 1 µm, w₀ = 10 µm;
  - z₀ = 314 µm, θ₀ = 32 mrad;
  - h = 3 µm (≈ z₀/100);
  - w(h)/w₀ − 1 = 0,005%, ζ(h) = 9,5 mrad.
- **Wniosek na slajdzie:** „Na wejściu (przewężenie): w = w₀, R = ∞, ζ = 0 → podajemy tylko profil; resztę liczy COMSOL”.
- **Obraz:** `img/schemat_modelu.png` (można go też zostawić na slajd 7).
- **Powiedzieć:** referat, slajd 3.

### 4. Mody Hermite'a-Gaussa *(P1, 3:00–4:15)*
- **Na slajdzie:**
  - wzór HG;
  - H₀…H₄;
  - „(l+1)(m+1) plam, sąsiednie z przeciwnym znakiem”;
  - „faza Gouya (l+m+1)ζ”;
  - „brak normalizacji → 1 V/m dla HG(0,0) do ~220 V/m dla HG(4,4)”.
- **Obraz:** `../zdj/orders/3x3_all_orders_hermite_ex.gif` (znaki plam) albo `…_enorm.gif`.
- **Powiedzieć:** referat, slajd 4.

### 5. Mody Laguerre'a-Gaussa, OAM, interferogramy *(P1, 4:15–6:00)*
- **Na slajdzie:**
  - wzór LG;
  - p+1 pierścieni, wir na osi;
  - e^{imφ} → helikalny front, OAM = mħ na foton;
  - wzór całki OAM;
  - mała ramka o znaku: konwencja e^{+iωt} → fizycznie L_z = −mħ.
- **Obrazy:**
  - `../laguerre/lg_sweep_l.gif`;
  - `../laguerre/lg_oam/interferogram_widelec_l0.png`;
  - `../laguerre/lg_oam/interferogram_spirala_l0.png`.
- **Powiedzieć:** referat, slajd 5. Na końcu przekazać głos: „Teraz kolega pokaże, jak to liczy COMSOL”.

### 6. Metoda obwiedni (ewbe) *(P2, 6:00–7:15)*
- **Na slajdzie:**
  - E = E₁·e^(−ik₁·r);
  - równanie dla E₁;
  - porównanie: pełne pole wymaga ~5 elementów 2. rzędu na λ (dokumentacja COMSOL, s. 30), a obwiednia tylko tylu, ile trzeba do zmian E₁;
  - nasze ustawienia: *Unidirectional*, k₁ = (0, 0, ewbe.k0).
- **Obraz:** `img/schemat_modelu.png`.
- **Powiedzieć:** referat, slajd 6.

### 7. Budowa modelu krok po kroku *(P2, 7:15–8:30)*
- **Na slajdzie:**
  - drzewo modelu w 9 krokach (Parameters → Functions → Variables → Geometry → Materials → ewbe + warunki brzegowe → Mesh → Study → Results);
  - wyróżnione wyrażenie pola padającego LG;
  - wynik testu: faza Gouya HG(2,3) = 0,05737 rad vs teoria 0,05729 rad.
- **Obraz:** `img/schemat_modelu.png`, jeśli nie był na slajdzie 6. Ewentualnie zrzut ekranu drzewa modelu z COMSOL-a, do zrobienia przed prezentacją.
- **Powiedzieć:** referat, slajd 7.
- **Opcjonalne demo na żywo (ryzykowne, ~30 s):** otworzyć `wzf_l04_02_laguerre_gauss.mph` (ma rozwiązania) i kliknąć *Matched Boundary Condition 1*, żeby pokazać wyrażenie. Plik ładuje się długo (98 MB), więc otworzyć go **przed** prezentacją.

### 8. Parametry i ich wpływ *(P2, 8:30–9:15)*
- **Na slajdzie:** tabela 5 wierszy: parametr → gdzie w COMSOL → efekt.
  - `lx`, `my`, `m`, `l`: kształt modu;
  - `w0`: szerokość i z₀;
  - `lam`: skala; w LG tylko przeskalowanie;
  - `n_real`, `n_im`: ośrodek, uwaga na k₁;
  - `Nl` i rozmiar siatki: dokładność i koszt.
- **Źródło treści:** referat, aneks A (skrócona wersja).
- **Powiedzieć:** referat, slajd 8.

### 9. Wyniki HG *(P2, 9:15–10:05)*
- **Na slajdzie:** trzy obrazy z podpisami:
  - wszystkie mody;
  - siatka zbyt gruba vs wystarczająca;
  - liczba elementów wzdłuż z.
- **Obrazy:**
  - `../zdj/orders/3x3_all_orders_hermite_enorm.gif`;
  - `../zdj/mesh_topology/3x3ExtremelyCoarseEnorm.png` obok `../zdj/mesh_topology/3x3FinerEnorm.png`;
  - `../zdj/mesh_sideways/3x3_NI1_50.png` obok `../zdj/mesh_sideways/3x3_NI10_50.png`.
- **Powiedzieć:** referat, slajd 9.

### 10. Wyniki LG *(P2, 10:05–11:15)*
- **Na slajdzie:**
  - pierścienie, widelec i spirala w jednym rzędzie;
  - pod nimi korkociąg i wykres OAM;
  - liczby: OAM 2,97 → 2,93 (wejście), kontrola 3,0000, moc zachowana z dokładnością do 0,02%.
- **Obrazy:**
  - `../laguerre/lg_oam/interferogram_widelec_sweep_l.gif`;
  - `../laguerre/lg_oam/interferogram_spirala_sweep_l.gif`;
  - `../laguerre/lg_oam/korkociag_rotacja_fazy.gif`;
  - `../laguerre/lg_oam/oam_na_foton.png`.
- **Powiedzieć:** referat, slajd 10. Zakończyć przejściem: „spadek OAM o 0,3–0,6% to błąd numeryczny → dyspersja numeryczna”.

### 11. Dyspersja numeryczna: eksperyment *(P2, 11:15–12:45)*
- **Na slajdzie:**
  - definicja w jednym zdaniu;
  - „siatka musi rozdzielić obwiednię ~ e^(−i(nk₀−k₁)z)”;
  - trzy wnioski:
    - k₁ = n·k₀ → błąd niezależny od Nl;
    - k₁ niedopasowane → potrzeba ~5 el. na „długość fali” obwiedni;
    - pełne pole przy < 3 el./λ → zła długość fali i amplituda;
  - dodatkowa obserwacja: amplituda → 2k₁/(k₁+nk₀) przy złym k₁.
- **Obrazy:**
  - główny: `img/dyspersja_pole_n1.png` (najbardziej obrazowy);
  - drugi: `img/dyspersja_blad_fazy.png`;
  - zapasowe: `img/dyspersja_pole_n15.png`, `img/dyspersja_amplituda.png` (do pytań albo jako slajd 11b, jeśli zostanie czas).
- **Powiedzieć:** referat, slajd 11. Zakończyć zdaniem o serii n i κ w projekcie i o poprawce `k1 = n_real*ewbe.k0`.

### 12. Zasoby i powtarzalność *(P2, 12:45–13:25)*
- **Na slajdzie:**
  - lewa kolumna, zasoby COMSOL:
    - Wave Optics → ewbe;
    - Matched BC;
    - Swept mesh;
    - Wavelength Domain + Parametric Sweep;
    - Integration, Global/Point Evaluation;
    - Animation/Data export.
  - prawa kolumna, koszt:
    - LG: 992 682 stopni swobody, 2:45 min, 23,8 GB RAM;
    - HG: 25 modów w 1:35 min.
  - stopka: „wszystkie kroki odtworzenia i pułapki – aneks D referatu”.
- **Źródło treści:** referat, aneksy B, C, D.
- **Powiedzieć:** referat, slajd 12.

### 13. Co jeszcze można zmienić *(P2, 13:25–14:00)*
- **Na slajdzie:** 4–5 punktów z aneksu G:
  - h = z₀ (dyfrakcja, Gouya π/4·(N+1));
  - mniejsze w₀ (granica przybliżenia przyosiowego);
  - HG10 + i·HG01 = LG z m = 1;
  - m → −m;
  - polaryzacja kołowa.
- **Powiedzieć:** referat, slajd 13.

### 14. Podsumowanie *(P1, 14:00–14:30)*
- **Na slajdzie:** trzy punkty z referatu (slajd 14), potem „Dziękujemy – pytania?”.
- **Obraz:** `../laguerre/lg_oam/korkociag_rotacja_fazy.gif` jako tło.

---

## Przygotowanie przed wystąpieniem

- [ ] Wstawić obrazy w slajdy. GIF-y sprawdzić w trybie pokazu.
- [ ] Wzory przepisać do edytora równań (LaTeX z referatu można wkleić do PowerPointa jako równanie).
- [ ] Jeśli będzie demo: otworzyć `wzf_l04_02_laguerre_gauss.mph` w COMSOL-u przed startem (ładowanie 98 MB) i ustawić widok na *Matched Boundary Condition 1*.
- [ ] Opcjonalnie zrobić zrzut drzewa modelu w COMSOL-u na slajd 7.
- [ ] Próba z zegarem. P1 kończy przed 6:00, P2 przed 14:00.
- [ ] Mieć pod ręką referat z aneksami A–G jako ściągę na pytania.

## Ściąga liczb

| Wielkość | Wartość |
|---|---|
| λ, w₀, h, promień obszaru | 1 µm, 10 µm, 3 µm, 50 µm |
| z₀, θ₀, ζ(h) | 314 µm, 32 mrad, 9,55 mrad |
| Faza Gouya HG(2,3): symulacja / teoria | 0,05737 / 0,05729 rad (0,13%) |
| OAM (m = 3): wejście / wyjście / kontrola analityczna | 2,970–2,928 / 2,960–2,909 / 3,0000 |
| P_wyj / P_wej (LG) | ≥ 0,99983 |
| LG: stopnie swobody, czas, RAM | 992 682; 2 min 45 s; 23,8 GB |
| Dyspersja, k₁ = n·k₀ | \|Δφ\| = 0,0007 rad dla każdego Nl |
| Dyspersja, n = 1,5, k₁ = k₀ | \|Δφ\| = 3,04 → 0,26 → 0,073 → 0,012 rad (Nl = 1 → 8) |
| Dyspersja, pełne pole, n = 1, Nl = 8 | amplituda 0,60 zamiast 1; zła długość fali |

## Prawdopodobne pytania z sali

1. **Dlaczego obszar ma tylko 3 µm?**
   Koszt: model LG z 10 warstwami zajął już 24 GB RAM. Przy h ≪ z₀ widzimy strukturę modu, a nie dyfrakcję. Dzięki metodzie obwiedni dłuższy obszar (h ≈ z₀) nie wymagałby proporcjonalnie gęstszej siatki (aneks G).
2. **Skoro nie wpisaliście fazy Gouya, skąd ona w modelu?**
   Wynika z rozwiązania równań Maxwella. Test HG(2,3) daje 0,05737 rad wobec 0,05729 rad z teorii.
3. **Dlaczego OAM to 2,93–2,97, a nie 3?**
   To dyskretyzacja pola. Ta sama całka dla wzoru analitycznego daje 3,0000, a błąd rośnie z p, bo struktura radialna jest drobniejsza.
4. **Czym różni się *Matched Boundary Condition* od *Scattering BC* czy PML?**
   Według dokumentacji COMSOL-a (s. 161) *Matched BC* uwzględnia znany kierunek fali i mało odbija także przy dużych kątach. *Scattering BC* działa dobrze przy fali prawie prostopadłej do brzegu. PML to dodatkowa warstwa pochłaniająca, która wymaga własnej siatki. U nas P_wyj/P_wej ≥ 0,9998, więc odbicia są pomijalne.
5. **Czy wyniki serii z różnymi n i κ są błędne?**
   Są jakościowe: k₁ = k₀ i za mało elementów wzdłuż z (slajd 11). Poprawka to `k1 = n_real*ewbe.k0` i gęstsza siatka.
6. **Dlaczego mod z e^{+imφ} ma L_z = −mħ?**
   COMSOL używa konwencji e^{+iωt}, a literatura fizyczna e^{−iωt}. Amplitudy w tych konwencjach są sprzężone, więc znak się odwraca.
7. **Dlaczego maksimum |E| HG(4,4) to ~220 V/m?**
   Wielomiany Hermite'a nie są unormowane. Do porównań mocy trzeba podzielić amplitudę przez √(2^{l+m} l! m!).
8. **Dlaczego siatka zamiatana (*Swept*)?**
   To domyślny typ dla ewbe w 3D. Warstwy idą wzdłuż kierunku propagacji, więc łatwo osobno kontrolować rozdzielczość poprzeczną (N_T) i podłużną (N_L).
9. **Co to jest dyspersja numeryczna, jednym zdaniem?**
   Fala na siatce biegnie z trochę inną liczbą falową niż prawdziwa, więc błąd fazy narasta z każdą długością fali. Dlatego potrzeba ~5 elementów 2. rzędu na λ, chyba że, jak w metodzie obwiedni, szybką fazę wydzielimy z góry.
