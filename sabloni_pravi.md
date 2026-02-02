# MATLAB Šabloni - Kompletni vodič za zadatke iz Teorije Sistema

---

## 📊 1. CRTANJE GRAFIKA

### 1.1 Osnovni 2D Grafici

#### Linijski grafikon (plot)
```matlab
% Osnovni linijski grafikon
x = [2 2.5 3 3.5];
y = [3 5 7 9];
plot(x, y)
xlabel('x'); ylabel('y'); title('Naslov');

% Sa formatiranjem linije
plot(x, y, 'r--', 'LineWidth', 2)  % crvena isprekidana linija debljine 2

% Format: 'boja+tip_linije+oznaka'
% Boje: 'r' crvena, 'b' plava, 'g' zelena, 'k' crna, 'y' žuta, 'w' bela
% Linije: '-' puna, '--' isprekidana, ':' tačkasta, '-.' tačka-crtica
% Oznake: 'o' krug, 'x' iks, '*' zvezda, '+' plus, 'v' trougao dole, '^' trougao gore
```

#### Više funkcija na jednom grafiku
```matlab
t = 0:pi/100:2*pi;
y1 = sin(t);
y2 = sin(t-0.5);
y3 = sin(t-1);
plot(t, y1, t, y2, t, y3)
legend('sin(t)', 'sin(t-0.5)', 'sin(t-1)');
```

#### Polinomi - polyval
```matlab
% Crtanje polinoma p(x) = x^5 - 12.1x^4 + 40.5x^3 - 17x^2 - 71.5x + 35.8
p = [1 -12.1 40.5 -17 -71.5 35.8];
x = -1.5:0.1:6.7;
y = polyval(p, x);
plot(x, y)
xlabel('x'); ylabel('y');
```

### 1.2 Podgrafici (subplot)
```matlab
figure
subplot(3, 1, 1);  % 3 reda, 1 kolona, pozicija 1
plot(t, y1);
title('Prvi grafik');

subplot(3, 1, 2);
plot(t, y2, 'r--');
title('Drugi grafik');

subplot(3, 1, 3);
plot(t, y3, 'g-.');
title('Treći grafik');
```

### 1.3 Crtanje funkcija (fplot)
```matlab
% Osnovni fplot
fplot(@(x) x^2 + 4*sin(2*x) - 1, [-2.5, 2.5])

% Sa hold on za više funkcija
fplot(@(x) x^2, [1 4.5])
hold on
fplot(@(x) 2^x, [1 4.5])
hold off

% Dodavanje teksta i elemenata - koristi plottools
plottools
```

### 1.4 Specijalni grafici

#### Histogram (bar)
```matlab
% Vertikalni histogram
god = [2016:2020];
pro = [27194, 28962, 30730, 34202, 26165];
bar(god, pro, 'r');
xlabel('Godina'); ylabel('Prodaja');

% Horizontalni histogram
barh(god, pro, 'r');
```

#### Diskretni grafikon (stem)
```matlab
x = -1.5:0.05:1.5;
p = [1 0 -3 0 2 0 0];
stem(x, polyval(p, x))
xlabel('x'); ylabel('y');
```

#### Stepenasti grafikon (stairs)
```matlab
stairs(x, polyval(p, x))
```

#### Kružni grafikon (pie)
```matlab
ocene = [64 23 18 13 11 7];
pie(ocene);
title('Ocene studenata');
```

### 1.5 3D Grafici
```matlab
x = [1 2 3 4];
y = [1 4 9 16];
z = [1 8 27 64];
plot3(x, y, z, 'o', 'LineWidth', 2)
grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
```

---

## 🌊 2. SIGNALI I KARAKTERISTIČNE FUNKCIJE

### 2.1 Heaviside funkcija (jedinična odskočna)

#### Kontinualno vreme
```matlab
clear all
syms t

% Postavljanje vrednosti na 0 i 1 (preporučeno za TS)
defaultParam = sympref('HeavisideAtOrigin', 1);

h = heaviside(t);
fplot(h, [-10, 10])
ylim([-1, 2]);
xlabel('t'); ylabel('h(t)');
title('Jedinicna odskocna funkcija');

% Vratiti na default
sympref('default');
```

#### Diskretno vreme
```matlab
k = -10:1:10;
defaultParam = sympref('HeavisideAtOrigin', 1);
h = heaviside(k);
stem(k, h);
xlabel('k'); ylabel('h(k)');
sympref('default');
```

#### Pomerena Heaviside funkcija
```matlab
syms t
% h(t-2) - pomerena za 2 u desno
fplot(heaviside(t-2), [-10, 10])

% Kombinacija: h(t) - h(t-2)
fplot(heaviside(t) - heaviside(t-2), [-10, 10])
```

#### Laplasova i Z transformacija
```matlab
syms t
h = heaviside(t);
laplace(h)  % Rezultat: 1/s

syms k
defaultParam = sympref('HeavisideAtOrigin', 1);
h = heaviside(k);
ztrans(h)   % Rezultat: z/(z-1)
sympref('default');
```

### 2.2 Dirac delta funkcija

#### Kontinualno vreme
```matlab
syms t
delta = dirac(t);

% Veza sa Heaviside: δ(t) = dh(t)/dt
diff(heaviside(t))  % Rezultat: dirac(t)

% Laplasova transformacija
laplace(delta)  % Rezultat: 1
```

### 2.3 Kronecker delta (kvazi-delta) - diskretno

```matlab
% Diskretni ekvivalent Dirac funkcije
k = -10:1:10;
k_delta = kroneckerDelta(sym(k));
stem(k, double(k_delta));
xlabel('k'); ylabel('δ(k)');
title('Kvazi Dirakova funkcija (Kronecker delta)');

% Z transformacija
syms k
k_delta = kroneckerDelta(k);
ztrans(k_delta)  % Rezultat: 1
```

### 2.4 Pravougaoni impuls

#### Kontinualno vreme
```matlab
syms t
f = rectangularPulse(-1, 1, t);  % od -1 do 1
fplot(f)
ylim([-0.5, 1.5]);
xlabel('t'); ylabel('Π(t)');

% Alternativno - razlika dve Heaviside funkcije
% fplot(heaviside(t+1) - heaviside(t-1))
```

#### Diskretno vreme
```matlab
k = -5:1:5;
defaultParam = sympref('HeavisideAtOrigin', 1);
f = rectangularPulse(-2, 2, k);
stem(k, f)
xlabel('k'); ylabel('f(k)');
sympref('default');
```

### 2.5 Trougaoni impuls

#### Kontinualno vreme
```matlab
syms t
fplot(triangularPulse(t), [-3, 3])  % default: -1, 0, 1
ylim([-0.5, 1.5]);
xlabel('t'); ylabel('f(t)');

% Sa eksplicitnim granicama a, b, c
% triangularPulse(a, b, c, t)
```

#### Diskretno vreme
```matlab
k = -3:1:3;
f = triangularPulse(-2, 0, 2, k);
stem(k, f)
xlabel('k'); ylabel('f(k)');
```

### 2.6 Diskretizacija (Shannon-Nyquist)

```matlab
clear all
% Teorema: fs >= 2*fg (učestanost odabiranja >= 2*granična učestanost)

t = 0:0.001:1;
freq = 15;  % granična učestanost 15 Hz
f = sin(2*pi*freq*t);

% Analogi signal
subplot(3, 1, 1);
plot(t, f);
xlabel('Vreme'); ylabel('Amplituda');
title('Analogni signal');

% Semplovan signal na fs = 2*freq (minimalno po teoremi)
freqS = 30;
n = 0:1/freqS:1;
f1 = sin(2*pi*freq*n);
subplot(3, 1, 2);
stem(n, f1);
xlabel('Vreme'); ylabel('Amplituda');
title('Semplovan signal na fs=2*freq');

% Rekonstruisan signal
subplot(3, 1, 3);
plot(n, f1);
xlabel('Vreme'); ylabel('Amplituda');
title('Rekonstruisan signal');
```

---

## 🔄 3. PRENOSNE FUNKCIJE

### 3.1 Kontinualno vreme - tf()

#### Osnovna forma
```matlab
% G(s) = (s+1)/(s^2+3s+2)
num = [1 1];      % brojilac: s + 1
den = [1 3 2];    % imenilac: s^2 + 3s + 2
Gs = tf(num, den)

% Ili korišćenjem simboličke promenljive
s = tf('s');
Gs = (s+1)/(s^2 + 3*s + 2);
```

#### MIMO sistemi (više ulaza/izlaza)
```matlab
% Sistem sa 2 ulaza i 1 izlaz
num1 = [2, 0]; den1 = [4, 1];
num2 = [3]; den2 = [5, 1];
Gs1 = tf(num1, den1);
Gs2 = tf(num2, den2);
Gs = [Gs1, Gs2]  % horizontalno spajanje = više ulaza

% Sistem sa 1 ulaz i 2 izlaza
num = {[1 1]; 1};
den = {[1 2 2]; [1 0]};
Gs = tf(num, den)  % vertikalno = više izlaza

% Provera dimenzija
size(Gs)
```

### 3.2 Diskretno vreme - tf(num, den, Ts)

```matlab
% G(z) = (z+0.5)/(z^2+1.5z+2), Ts = 0.4s
numd = [1, 0.5];
dend = [1, 1.5, 2];
Tsd = 0.4;
Gz = tf(numd, dend, Tsd)

% Ili sa simboličkom promenljivom
z = tf('z', 0.4);
Gz = (z+0.5)/(z^2 + 1.5*z + 2);
```

### 3.3 Zadavanje preko nula, polova i pojačanja (zp2tf)

```matlab
% G(s) = k*(s-z1)(s-z2)../((s-p1)(s-p2)..)
z = [];              % nule (zeros) - koreni brojioca
p = [0 -1 -2];       % polovi (poles) - koreni imenioca
k = 7;               % pojačanje (gain)

[num, den] = zp2tf(z, p, k)
Gs = tf(num, den)
```

---

## 🔗 4. VEZE SISTEMA

### 4.1 Redna (serijska) veza - series()

```matlab
% G1 i G2 vezani redom => G = G1 * G2
num1 = [7]; den1 = [1 3 2 0];
num2 = [1 3]; den2 = [1 5];

[ns, ds] = series(num1, den1, num2, den2)
sys = tf(ns, ds)
```

### 4.2 Paralelna veza - parallel()

```matlab
% G1 i G2 vezani paralelno => G = G1 + G2
[np, dp] = parallel(num1, den1, num2, den2)
sys = tf(np, dp)
```

### 4.3 Povratna sprega - feedback()

```matlab
% G1 u direktnoj grani, G2 u povratnoj
% -1 za negativnu povratnu spregu (najčešće)
% +1 za pozitivnu povratnu spregu
% G = G1/(1 ± G1*G2)

[nf, df] = feedback(num1, den1, num2, den2, -1)
sys = tf(nf, df)
```

### 4.4 Kompleksni primer - kombinovane veze

```matlab
% Primer: negativna povratna sprega G1,H1 -> redna sa G2 -> 
% paralelna sa G3 -> negativna povratna sa H2

% Negativna povratna sprega G1, H1
[n, d] = feedback([0 0 5], [1 1 0], [1 0], [1 4], -1);

% Redna veza sa G2
[n, d] = series(n, d, [0 2], [1 0]);

% Paralelna veza sa G3
[n, d] = parallel(n, d, [2], [1]);

% Negativna povratna sprega sa H2
[n, d] = feedback(n, d, [5 0], [1 2], -1);

Gs = tf(n, d)
```

---

## 📐 5. MASON FORMULA (za složene blok dijagrame)

### 5.1 Priprema .net fajla

```
% Fajl: zadatak.net
% Format: broj_grane početni_čvor krajnji_čvor pojačanje
1 1 2 1
2 2 3 1/(s+4)
3 3 4 -1
4 4 2 1/(s+1)
5 3 5 1
6 5 4 -1
7 5 6 s+1
8 6 7 1
```

### 5.2 Poziv Mason funkcije

```matlab
% Preuzmi mason.m sa MathWorks File Exchange
% Stavi u Current Folder

[num, den] = mason('zadatak.net', 1, 7)
% Argumenti: ime_fajla, početni_čvor, krajnji_čvor

% Konverzija u simboličke promenljive
num = str2sym(num);
den = str2sym(den);

% Pojednostavljivanje
Gs = simplify(num/den)
pretty(Gs)  % Lepo formatiran prikaz
```

### 5.3 Primer za diskretne sisteme

```matlab
% Za diskretne sisteme koristi z u .net fajlu
% Primer: 1/z, 2/z, itd.

[num, den] = mason('zadatak_z.net', 1, 7);
num = str2sym(num);
den = str2sym(den);
Gz = simplify(num/den)
```

---

## 🎛️ 6. MODEL U PROSTORU STANJA

### 6.1 Kontinualni sistem - ss()

#### Osnovna forma
```matlab
% dx/dt = Ax + Bu
% y = Cx + Du

A = [0, 1; -4, -2];
B = [0; -2];
C = [1, 0];
D = [0];

sys = ss(A, B, C, D)

% Pristup matricama
sys.A
sys.B
sys.C
sys.D
```

#### Kreiranje od prenosne funkcije
```matlab
% Iz prenosne funkcije u prostor stanja
Gs = tf([1 1], [1 2 2]);
sys = ss(Gs)

% Ili eksplicitno
[A, B, C, D] = tf2ss([1 1], [1 2 2])
```

### 6.2 Diskretni sistem - ss(A, B, C, D, Ts)

```matlab
A = [-1, 0; 0, -1];
B = [1; 0];
C = eye(2);  % matrica identiteta 2x2
D = zeros(2, 1);
Ts = 0.4;

sys = ss(A, B, C, D, Ts)
```

### 6.3 Konverzija iz prostora stanja u prenosnu funkciju

#### Metoda 1: tf(sys)
```matlab
sys = ss(A, B, C, D);
Gs = tf(sys)  % Prikaže prenosnu funkciju
```

#### Metoda 2: ss2tf() - za koeficijente
```matlab
% Za SISO sisteme
[num, den] = ss2tf(A, B, C, D)

% Za MIMO - zadaj koji ulaz (poslednji argument)
[num, den] = ss2tf(A, B, C, D, 1)  % prvi ulaz
[num, den] = ss2tf(A, B, C, D, 2)  % drugi ulaz
```

#### Metoda 3: Simbolička - preko inverzije
```matlab
syms s
Gs = C*inv(s*eye(3) - A)*B + D
Gs = simplify(Gs)
```

### 6.4 Konverzija preko nula i polova

```matlab
% Iz prenosne funkcije u nule/polovi
[z, p, k] = tf2zp(num, den)

% Iz prostora stanja u nule/polovi
[z, p, k] = ss2zp(A, B, C, D)
```

---

## 🔀 7. KONVERZIJE IZMEĐU KONTINUALNIH I DISKRETNIH SISTEMA

### 7.1 Kontinualan → Diskretan (c2d)

```matlab
% Osnovni primer
Gs = tf([1 -1], [1 4 5]);
Ts = 0.1;

% Zero-Order Hold (ZOH) - najčešća metoda
Gz = c2d(Gs, Ts, 'zoh')

% First-Order Hold (FOH)
Gz = c2d(Gs, Ts, 'foh')

% Poređenje odziva
step(Gs, Gz)
impulse(Gs, Gz)
```

### 7.2 Diskretan → Kontinualan (d2c)

```matlab
% Samo ZOH metoda podržana
sysc = d2c(Gz, 'zoh')
```

---

## ✅ 8. OSOBINE MODELA SISTEMA

### 8.1 Upravljivost - ctrb()

```matlab
% Sistem je upravljiv ako je rank(matriceUpravljivosti) = n
% Matrica upravljivosti: [B AB A²B ... Aⁿ⁻¹B]

A = [1, 1; 4, -2];
B = [1, -1; 1, -1];

matricaUpravljivosti = ctrb(A, B)
brNeupravljivihStanja = length(A) - rank(matricaUpravljivosti)

if brNeupravljivihStanja == 0
    disp('Sistem je upravljiv')
else
    fprintf('Sistem ima %d neupravljivih stanja\n', brNeupravljivihStanja)
end
```

### 8.2 Osmotrivost - obsv()

```matlab
% Sistem je osmotriv ako je rank(matriceOsmotrivosti) = n
% Matrica osmotrivosti: [C; CA; CA²; ...; CAⁿ⁻¹]

A = [1, 1; 4, -2];
C = [1, 0; 0, 1];

matricaOsmotrivosti = obsv(A, C)
brNeosmotrivihStanja = length(A) - rank(matricaOsmotrivosti)

if brNeosmotrivihStanja == 0
    disp('Sistem je osmotriv')
else
    fprintf('Sistem ima %d neosmotrivih stanja\n', brNeosmotrivihStanja)
end
```

---

## 📈 9. ODZIVI SISTEMA

### 9.1 Step odziv - step()

```matlab
% Odziv na jediničnu odskočnu funkciju
num = [1 1];
den = [1 4 5];
Gs = tf(num, den);

% Osnovni step odziv
step(Gs)

% Sa vremenskim vektorom
t = 0:0.1:5;
step(Gs, t)

% Step odziv za više sistema
step(Gs1, Gs2, Gs3)
```

### 9.2 Impulsni odziv - impulse()

```matlab
% Odziv na Dirac delta funkciju
impulse(Gs)
impulse(Gs, t)

% Poređenje kontinualnog i diskretnog
impulse(Gs, Gz)
```

---

## 🎯 10. STABILNOST KONTINUALNIH SISTEMA

### 10.1 Asimptotska stabilnost - Sopstvene vrednosti

**Uslov:** Sve sopstvene vrednosti matrice A moraju imati **Re(λ) < 0** (leva poluravan)

#### Metoda 1: eig()
```matlab
A = [-1, 1; -4, -5];
sv = eig(A)

% Provera stabilnosti
if all(real(sv) < 0)
    disp('Sistem je asimptotski stabilan')
else
    disp('Sistem NIJE asimptotski stabilan')
end
```

#### Metoda 2: Rešavanje karakteristične jednačine
```matlab
syms s
jd = det(s*eye(2) - A)
polovi = solve(jd)

% Provera
if all(real(polovi) < 0)
    disp('Sistem je stabilan')
end
```

### 10.2 Asimptotska stabilnost - Routh kriterijum

**Uslov:** Svi elementi **prve kolone** Routh tabele moraju biti **istog znaka**

```matlab
% Preuzmi routh.m sa MathWorks File Exchange

syms lambda, eps
A = [-1, 1; -4, -5];
p = det(lambda*eye(2) - A)

% Formiranje Routh tabele
ra = routh(sym2poly(p), eps);
r = double(ra)

% Provera prve kolone
prva_kolona = r(:, 1);
if all(prva_kolona > 0) || all(prva_kolona < 0)
    disp('Sistem je stabilan (svi elementi prve kolone istog znaka)')
else
    disp('Sistem NIJE stabilan (mešoviti znaci u prvoj koloni)')
end
```

### 10.3 OUOI (BIBO) stabilnost - Polovi prenosne funkcije

**Uslov:** Svi polovi prenosne funkcije moraju imati **Re(p) < 0**

#### Metoda 1: Simbolička
```matlab
syms s
Gs = (s^2 - 2*s + 1)/(s^3 + 4*s^2 + 5*s + 2);
[num, den] = numden(Gs);
polovi = solve(den)

if all(real(polovi) < 0)
    disp('Sistem je OUOI stabilan')
end
```

#### Metoda 2: Numerička - tf2zp()
```matlab
num = [1, -2, 1];
den = [1, 4, 5, 2];
[z, p, k] = tf2zp(num, den)

if all(real(p) < 0)
    disp('Sistem je OUOI stabilan')
else
    disp('Sistem NIJE OUOI stabilan')
end
```

### 10.4 Grafička analiza - step odziv
```matlab
% Stabilan sistem se vraća u stacionarno stanje
num = [1, 1];
den = [1, 5, 16];
step(num, den)
% Ako grafik konvergira → STABILAN

% Nestabilan sistem
den_nestabilan = [1, -1, -2];
step(num, den_nestabilan)
% Ako grafik ide u beskonačnost → NESTABILAN
```

---

## 🎲 11. STABILNOST DISKRETNIH SISTEMA

### 11.1 Asimptotska stabilnost - Sopstvene vrednosti

**Uslov:** Sve sopstvene vrednosti matrice A moraju biti **|λ| < 1** (jedinični krug)

```matlab
A = [-5, -2, -1; 4, 0, 0; 0, 1, 0];
sv = eig(A)

% Provera
if all(abs(sv) < 1)
    disp('Sistem je asimptotski stabilan')
else
    disp('Sistem NIJE asimptotski stabilan')
    fprintf('Broj sopstvenih vrednosti van jediničnog kruga: %d\n', sum(abs(sv) >= 1))
end
```

### 11.2 Asimptotska stabilnost - Jury kriterijum

**Uslov:** Prvi elementi u **neparnim redovima** Jury tabele moraju biti **pozitivni**

```matlab
% Preuzmi jury.m sa MathWorks File Exchange

syms z
A = [-5, -2, -1; 4, 0, 0; 0, 1, 0];
p = det(z*eye(3) - A)

coef = sym2poly(p);
[J, C] = jury(coef);

% J - kompletna Jury tabela
% C - elementi prve kolone na neparnim pozicijama

if all(sign(C) == 1)
    disp('Sistem je stabilan (svi elementi C pozitivni)')
else
    disp('Sistem NIJE stabilan')
    fprintf('Broj polova van jediničnog kruga: %d\n', sum(sign(C) ~= 1))
end
```

### 11.3 OUOI stabilnost - Polovi prenosne funkcije

**Uslov:** Svi polovi prenosne funkcije moraju biti **|p| < 1**

```matlab
num = [1, 3];
den = [1, 4, 5, 2];
[z, p, k] = tf2zp(num, den)

% Provera
moduli = abs(p);
if all(moduli < 1)
    disp('Sistem je OUOI stabilan')
else
    disp('Sistem NIJE OUOI stabilan')
    fprintf('Polovi van jediničnog kruga: %d\n', sum(moduli >= 1))
end
```

---

## 📋 12. REZIME - KRITERIJUMI STABILNOSTI

| **Tip sistema** | **Domen** | **Uslov stabilnosti** | **Funkcija** |
|-----------------|-----------|----------------------|--------------|
| **Kontinualan (prostor stanja)** | s | Re(λ) < 0 | `eig(A)` |
| **Kontinualan (Routh)** | s | Prva kolona istog znaka | `routh()` |
| **Kontinualan (OUOI)** | s | Re(polovi) < 0 | `tf2zp()` |
| **Diskretan (prostor stanja)** | z | \|λ\| < 1 | `eig(A)` |
| **Diskretan (Jury)** | z | C svi pozitivni | `jury()` |
| **Diskretan (OUOI)** | z | \|polovi\| < 1 | `tf2zp()` |

---

## 🔧 13. KORISNI TRIKOVI I SAVETI

### 13.1 Rad sa simboličkim promenljivim
```matlab
syms s t z k x

% Konverzija string → simbolička
num_str = '2*s + 1';
num_sym = str2sym(num_str);

% Konverzija simbolička → koeficijenti polinoma
p = s^3 + 5*s^2 + 8*s + 4;
coef = sym2poly(p);  % [1 5 8 4]

% Pojednostavljivanje
expr = (s^2 - 1)/(s - 1);
simplify(expr)  % s + 1

% Lepo formatiran ispis
pretty(expr)
```

### 13.2 Provera dimenzija sistema
```matlab
sys = tf(num, den);
size(sys)  % Prikazuje: "Transfer function with X outputs and Y inputs"
```

### 13.3 Učitavanje podataka iz Excel-a
```matlab
data = xlsread('filmovi.xlsx', 'Sheet1', 'B2:D32');
% Argumenti: ime_fajla, list, opseg_ćelija
```

### 13.4 Sortiranje podataka
```matlab
% Sortiranje matrice po prvoj koloni (opadajuće)
sortedData = sortrows(data, 1, 'descend');

% Po više kolona
sortedData = sortrows(data, [1 2], {'descend', 'ascend'});
```

### 13.5 Aproksimacija krivom (polyfit)
```matlab
x = [1, 2, 3, 4, 5];
y = [2.1, 3.9, 6.2, 8.1, 9.8];

% Aproksimacija polinomom drugog reda
koef = polyfit(x, y, 2);

% Generisanje tačaka za crtanje
xp = min(x):0.1:max(x);
yp = polyval(koef, xp);

plot(x, y, 'ro', xp, yp, 'b-');
legend('Podaci', 'Aproksimacija');
```

### 13.6 Rad sa objektima tf i ss
```matlab
% Pristup elementima
sys = tf([1 2], [1 3 2]);
sys.Numerator    % {[1 2]}
sys.Denominator  % {[1 3 2]}

% Za prostor stanja
sys_ss = ss(A, B, C, D);
sys_ss.A
sys_ss.B
```

---

## 📚 14. KOMPLETNI PRIMERI ZADATAKA

### Primer 1: Analiza stabilnosti kontinualnog sistema
```matlab
% Dat je sistem: G(s) = (s+1)/(s^3 + 4s^2 + 5s + 2)

% 1. Kreiranje prenosne funkcije
num = [1, 1];
den = [1, 4, 5, 2];
Gs = tf(num, den);

% 2. Konverzija u prostor stanja
sys = ss(Gs);

% 3. Provera upravljivosti i osmotrivosti
Ct = ctrb(sys.A, sys.B);
unco = length(sys.A) - rank(Ct);
fprintf('Broj neupravljivih stanja: %d\n', unco);

O = obsv(sys.A, sys.C);
unob = length(sys.A) - rank(O);
fprintf('Broj neosmotrivih stanja: %d\n', unob);

% 4. Asimptotska stabilnost (sopstvene vrednosti)
sv = eig(sys.A);
if all(real(sv) < 0)
    disp('Sistem je asimptotski stabilan');
end

% 5. OUOI stabilnost (polovi)
[z, p, k] = tf2zp(num, den);
if all(real(p) < 0)
    disp('Sistem je OUOI stabilan');
end

% 6. Grafička analiza
figure;
subplot(2,1,1);
step(Gs);
title('Step odziv');
subplot(2,1,2);
impulse(Gs);
title('Impulsni odziv');
```

### Primer 2: Diskretizacija i analiza diskretnog sistema
```matlab
% Kontinualni sistem
Gs = tf([1, -1], [1, 4, 5]);

% Diskretizacija
Ts = 0.1;
Gz_zoh = c2d(Gs, Ts, 'zoh');
Gz_foh = c2d(Gs, Ts, 'foh');

% Poređenje
figure;
step(Gs, Gz_zoh, Gz_foh);
legend('Kontinualan', 'ZOH', 'FOH');

% Analiza stabilnosti diskretnog
sys_d = ss(Gz_zoh);
sv_d = eig(sys_d.A);

fprintf('Sopstvene vrednosti diskretnog sistema:\n');
disp(sv_d);

if all(abs(sv_d) < 1)
    disp('Diskretan sistem je stabilan');
end

% Jury kriterijum
syms z
p = det(z*eye(length(sys_d.A)) - sys_d.A);
coef = sym2poly(p);
[J, C] = jury(coef);

if all(sign(C) == 1)
    disp('Potvrđeno Jury kriterijumom: sistem je stabilan');
end
```

### Primer 3: Kompleksni blok dijagram sa Mason formulom
```matlab
% 1. Kreiranje .net fajla (uradi van MATLAB-a)
% Fajl: sistem.net
% 1 1 2 1
% 2 2 3 1/(s+4)
% 3 3 4 -1
% 4 4 2 1/(s+1)
% 5 3 5 1
% 6 5 4 -1
% 7 5 6 s+1
% 8 6 7 1

% 2. Primena Mason formule
[num, den] = mason('sistem.net', 1, 7);

% 3. Konverzija i pojednostavljivanje
num = str2sym(num);
den = str2sym(den);
Gs = simplify(num/den);

fprintf('Prenosna funkcija sistema:\n');
pretty(Gs);

% 4. Kreiranje tf objekta
[num_vec, den_vec] = numden(Gs);
Gs_tf = tf(sym2poly(num_vec), sym2poly(den_vec));

% 5. Analiza
step(Gs_tf);
[z, p, k] = tf2zp(Gs_tf.Numerator{1}, Gs_tf.Denominator{1});
fprintf('\nPolovi sistema:\n');
disp(p);
```

---

## 🎓 15. DODATNI RESURSI

### Korisni linkovi:
- **Mason.m**: https://www.mathworks.com/matlabcentral/fileexchange/22-mason-m
- **Routh.m**: https://www.mathworks.com/matlabcentral/fileexchange/58-routh-m
- **Jury.m**: https://www.mathworks.com/matlabcentral/fileexchange/13904-jury

### Često korišćene funkcije:
```matlab
help funkcija      % Pomoć za funkciju
doc funkcija       % Detaljna dokumentacija
clear all          % Briše sve promenljive
clc                % Čisti Command Window
close all          % Zatvara sve prozore sa grafikonima
```

---

**Kraj šablona - Srećno sa zadacima! 🚀**
