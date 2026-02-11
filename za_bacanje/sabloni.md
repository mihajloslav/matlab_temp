# MATLAB Šabloni - Vodič za Učenje

## 1. Funkcije za crtanje grafika

### 2D Grafici
```matlab
plot(x, y, 'format_linije', 'ime_svojstva', vrednost_svojstva)  % linijski 2D
% format_linije: 'r-' (crvena linija), 'b--' (plava isprekidana), 'go' (zeleni krugovi)
% svojstva: 'LineWidth', 2 ili 'MarkerSize', 10

bar         % histogram
stem        % diskretni - za diskretne signale (tačke sa linijama)
stairs      % stepnasti - za step funkcije
pie         % kružni - procenat, distribucija
```

### 3D Grafici
```matlab
plot3       % linijski 3D
```

### Crtanje funkcija
```matlab
fplot(@(x)funkcija_od_X, granice, specifikacija_linije)
ylim([-1.5, 1.5])
plottools
```

### Labele i naslovi
```matlab
xlabel('x')
ylabel('y')
title('Prigušeni kosinus')
text(0.45, 0.7, 'y=exp(-x).*cos(6*pi*x)')
```

### Multipli grafici
```matlab
figure
subplot(3, 1, 1);
```

### Sistemski odzivi
```matlab
step(gz, gs)       % step odziv - odziv sistema na jedinični skok
impulse(gz, gs)    % impulsni odziv - odziv sistema na impuls
% gz - brojilac prenosne funkcije
% gs - imenilac prenosne funkcije
```

---

## 2. Specijalne funkcije

### Heaviside funkcija
```matlab
defaultParam = sympref('HeavisideAtOrigin', 1);  % postavlja vrednosti na 0 i 1
h = heaviside(5)
sympref('default');
```

### Dirac delta funkcija
```matlab
d = dirac(x)       % DELTA funkcija
```

### Kronecker delta
```matlab
k = -10:1:10;
k_delta = kroneckerDelta(sym(k));  % KVAZI DELTA
```

### Pravougaoni i trouglasti impulsi
```matlab
rectangularPulse(a, b, x)      % pravougaoni impuls
triangularPulse(a, b, c, x)    % trouglasti impuls
```

---

## 3. Prenosne funkcije

### Kontinualno vreme
**Osnovna forma:** Y(s) = G(s)U(s)
- **U(s)** - ulaz sistema (input)
- **Y(s)** - izlaz sistema (output)
- **G(s)** - prenosna funkcija sistema

```matlab
sys = tf(num, den)  % prenosna funkcija u kontinualnom vremenu
% num - brojilac (numerator) - koeficijenti polinoma u brojiocu
% den - imenilac (denominator) - koeficijenti polinoma u imeniocu
% Primer: G(s) = (s+1)/(s^2+3s+2) => num=[1 1], den=[1 3 2]
```

### Diskretno vreme
```matlab
sys = tf(num, den, Ts)  % sa periodom odabiranja Ts u sekundama
                        % Ts = -1 ili Ts = [] ako ne želimo specifičnu vrednost
```

### Konverzija iz nula/polovi/pojačanja
```matlab
% Zadavanje sistema preko nula, polova i pojačanja
z = [];              % nule (zeros) - koreni brojioca
p = [0 -1 -2];      % polovi (poles) - koreni imenioca
k = 7;              % pojačanje (gain) - konstanta ispred
[num1 den1] = zp2tf(z, p, k)  % zero-pole-gain to transfer function

% Faktorizovana forma: G(s) = k * (s-z1)(s-z2)... / (s-p1)(s-p2)...
```

---

## 4. Veze sistema

### Redna veza
```matlab
[ns ds] = series(num1, den1, num2, den2)
sys1 = tf(ns, ds)
```

### Paralelna veza
```matlab
[np, dp] = parallel(num1, den1, num2, den2)
sys2 = tf(np, dp)
```

### Povratna sprega
```matlab
[nf df] = feedback(num1, den1, num2, den2, -1)
% -1 za negativnu povratnu spregu (najčešće)
% +1 za pozitivnu povratnu spregu
% G1 u direktnoj grani, G2 u povratnoj grani
sys3 = tf(nf, df)
```

---

## 5. Mason formula

**Primena:** Računanje prenosne funkcije za složene blok dijagrame sa više petlji.

```matlab
[num, den] = mason('zadatak1.net', 1, 7)
% 'zadatak1.net' - fajl sa opisom blok dijagrama
% 1 - broj ulaznog čvora
% 7 - broj izlaznog čvora

num = str2sym(num);  % konvertuj string u simbolički izraz
den = str2sym(den);

Gs = simplify(num/den)  % pojednostavi prenosnu funkciju
pretty(Gs)              % lepo formatiran prikaz
```

---

## 6. Model u prostoru stanja

### Kontinualno vreme
```matlab
sys = ss(A, B, C, D)
```
**Reprezentacija prostora stanja:**
- **ẋ = Ax + Bu** (jednačina stanja)
- **y = Cx + Du** (jednačina izlaza)

Gde:
- **A** (n×n) - matrica sistema (dinamika)
- **B** (n×m) - matrica ulaza (uticaj ulaza na stanja)
- **C** (p×n) - matrica izlaza (koji izlazi se mere)
- **D** (p×m) - matrica direktne veze (obično nula)
- **n** - broj stanja, **m** - broj ulaza, **p** - broj izlaza

### Diskretno vreme
```matlab
sys = ss(A, B, C, D, Ts)
```

---

## 7. Konverzije modela

### Funkcije za konverziju
```matlab
tf()       % kreira prenosnu funkciju ili prevodi model iz prostora stanja
ss()       % kreira model u prostoru stanja ili prevodi prenosnu funkciju

ss2tf()    % kreira prenosnu funkciju preko reprezentacije u prostoru stanja
tf2ss()    % kreira model u prostoru stanja za datu prenosnu funkciju
tf2zp()    % kreira prenosnu funkciju preko nula i polova
ss2zp()    % kreira prenosnu funkciju preko nula i polova za model u prostoru stanja

sysd = c2d(sys, Ts, 'zoh')    % prevodi kontinualan sistem u diskretan ('zoh', 'foh')
sysc = d2c(sysd, 'zoh')       % prevodi diskretan sistem u kontinualan

size()     % prikazuje broj izlaza i ulaza sistema
```

---

## 8. Osobine modela sistema

### Upravljivost
**Značenje:** Da li možemo dovesti sistem iz bilo kog početnog stanja u bilo koje željeno stanje pomoću ulaza.

```matlab
matricaUpravljivosti = ctrb(A, B)  % [B AB A²B ... Aⁿ⁻¹B]
brNeupravljivihStanja = length(A) - rank(matricaUpravljivosti)
% rank = broj linearno nezavisnih kolona
% Ako je rank(matricaUpravljivosti) = n, sistem je potpuno upravljiv
```

### Osmotrivost
**Značenje:** Da li možemo odrediti sva unutrašnja stanja sistema samo na osnovu merenja izlaza.

```matlab
matricaOsmotrivosti = obsv(A, C)  % [C; CA; CA²; ...; CAⁿ⁻¹]
brNeosmotrivihStanja = length(A) - rank(matricaOsmotrivosti)
% Ako je rank(matricaOsmotrivosti) = n, sistem je potpuno osmotriv
```

---

## 9. Stabilnost sistema

### Provera stabilnosti - Step odziv
```matlab
num = [1 1]; 
den = [1 -1 -2];
t = 0:0.1:5;
step(num, den, t)
% Posmatraj da li sistem ide u stabilno stanje
```

---

## 10. Stabilnost kontinualnih sistema

**Pravilo:** Vremenski kontinualan LTI sistem je asimptotski stabilan ako i samo ako sve sopstvene vrednosti matrice A imaju **negativne realne delove** (leže u levoj poluraini).

**Zašto?** Sopstvene vrednosti određuju ponašanje sistema:
- **Re(λ) < 0** → eksponencijalno opadanje → **STABILAN**
- **Re(λ) = 0** → oscilacije konstantne amplitude → **granično stabilan**
- **Re(λ) > 0** → eksponencijalni rast → **NESTABILAN**

### Metoda 1: Sopstvene vrednosti
```matlab
A = [-1 1; -4 -5];
sv = eig(A)
% Proveri: real(sv) < 0
```

### Metoda 2: Karakteristična jednačina
```matlab
syms s
jd = det(s*eye(2) - A)
solve(jd)
```

### Metoda 3: Routh kriterijum
```matlab
syms lambda, eps
p = det(lambda*eye(2) - A)
ra = routh(sym2poly(p), eps);
r = double(ra)
```
**Pravilo:** Ako su sve vrednosti **prve kolone istog znaka**, sistem je **STABILAN**.

**Routh kriterijum:** Bez računanja korena, proverava stabilnost preko prve kolone Routh tabele.
- Sve pozitivne → stabilan
- Sve negativne → stabilan
- Mešoviti znaci → nestabilan (ima polove u desnoj poluraini)

### OUOI stabilnost (prenosna funkcija)
**OUOI = Ograničen Ulaz Ograničen Izlaz (BIBO - Bounded Input Bounded Output)**

Vremenski kontinualan LTI sistem je OUOI stabilan ako i samo ako svi **polovi prenosne funkcije G(s)** imaju negativne realne delove.

**Značenje:** Ako je ulaz ograničen, izlaz će takođe biti ograničen (neće ići u beskonačnost).

#### Metoda 1: Simbolička
```matlab
syms s
Gs = (s^2 - 2*s + 1)/(s^3 + 4*s^2 + 5*s + 2)
[num den] = numden(Gs)
polovi = solve(den)
% Proveri: real(polovi) < 0
```

#### Metoda 2: Numerička
```matlab
num = [1 -2 1];
den = [1 4 5 2];
[z p k] = tf2zp(num, den)
% Proveri: real(p) < 0
```

---

## 11. Stabilnost diskretnih sistema

**Pravilo:** Vremenski diskretan LTI sistem je asimptotski stabilan ako i samo ako su sve sopstvene vrednosti matrice A **po modulu manje od 1** (nalaze se u jediničnom krugu).

**Zašto jedinični krug?** U diskretnom vremenu:
- **|λ| < 1** → sekvenca se smanjuje → **STABILAN**
- **|λ| = 1** → sekvenca konstantna → **granično stabilan**
- **|λ| > 1** → sekvenca raste → **NESTABILAN**

### Metoda 1: Sopstvene vrednosti
```matlab
A = [-5 -2 -1; 4 0 0; 0 1 0];
sv = eig(A)
% Proveri: abs(sv) < 1
```

### Metoda 2: Karakteristična jednačina
```matlab
syms z
jd = det(z*eye(3) - A)
% jd = z^3 + 5*z^2 + 8*z + 4
polovi = solve(jd)
% Proveri: abs(polovi) < 1
```

### Metoda 3: Jury kriterijum
```matlab
syms z
p = det(z*eye(3) - A)
% p = z^3 + 5*z^2 + 8*z + 4
coef = sym2poly(p)
[J C] = jury(coef);
all(sign(C) == 1)  % Ako je rezultat 1 (true), sistem je stabilan

% Jury kriterijum - diskretni analog Routh kriterijuma
% Proverava da li su svi polovi unutar jediničnog kruga
% bez eksplicitnog računanja korena
```

### OUOI stabilnost (prenosna funkcija)
Vremenski diskretan LTI sistem je OUOI stabilan ako i samo ako su svi **polovi prenosne funkcije G(z)** po modulu manji od 1.

```matlab
num = [1 -2 1];
den = [1 4 5 2];
[z p k] = tf2zp(num, den)

abs(p) < 1  % Ako je svuda 1 (true), sistem je stabilan
```

---

## Rezime - Kriterijumi stabilnosti

| Tip sistema | Domen | Uslov stabilnosti |
|-------------|-------|-------------------|
| Kontinualan (prostor stanja) | s | Real(λ) < 0 (leva poluravan) |
| Kontinualan (prenosna funkcija) | s | Real(polovi) < 0 |
| Diskretan (prostor stanja) | z | \|λ\| < 1 (jedinični krug) |
| Diskretan (prenosna funkcija) | z | \|polovi\| < 1 |

