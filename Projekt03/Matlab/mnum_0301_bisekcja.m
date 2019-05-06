%todo dodac ograniczenie na ilosc iteracji

% usuwanie wartosci skrajnych
%   [wart indx] = max(cond_sym(:,rozmiar));
%   iteracje_sym_bezprzes(indx,rozmiar) = 0;

% warunek stopu
%   while tol(fun(x)) > 0.00001 & i < 200+1

% % wyznaczanie pierw f. kwadratowej
% function [x1 x2] = pierw_f_kwadratowej(mala_macierz)
%     a = 1;
%     b = -(mala_macierz(1,1)+mala_macierz(2,2));
%     c = (mala_macierz(1,1)*mala_macierz(2,2))-(mala_macierz(2,1)*mala_macierz(1,2));
%     
%     x1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
%     x2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
%     if abs(x2) > abs(x1)
%         x1 = x2;
%     end
%    %drugi pierwiastek ze wzorów Viete'a
%     x2 = ((-b)/a) - x1;
% end
%% 
% 
% 
% start

clc;
clear;
x = linspace(0,15,100);
%% 
% wykres

figure
y = wartosc_funkcji(x);
plot(x,y,'b-',[0 15], [0 0], 'k--');
legend({'f(x)', 'y=0'},'Location','southwest');
title('wykres funkcji f(x)=2,3*sin(x)+4*ln(x+2)-11');
%% 
% oszacowanie miejsc zerowych na podstawie rysunku (skrypt prof Tatjeskiego 
% mowi, zeby estymowac na poodstawie rysunku)

% wartosc_funkcji(7.25)
% wartosc_funkcji(8)
% wartosc_funkcji(9)
% wartosc_funkcji(12.5)
%% 
% przedzialy poczatkowe

    przedzialy = [1 8; 8 10; 10 15];

for i=1:3
    %wybor i-tego przedzialu
    sprawdzenie_przedzialu(przedzialy(i,:));
end
%% 
% metoda bisekcji

x_bisekcja = bisekcja(przedzialy);
x_bisekcja
y_bisekcja = wartosc_funkcji(x_bisekcja)
fprintf('najwiekszy blad zera przy metodzie bisekcji %f\n', najwieksze_zero(x_bisekcja));

figure
plot(x,y,'b-',[0 15], [0 0], 'k--', x_bisekcja, wartosc_funkcji(x_bisekcja), 'go');
legend({'f(x)', 'y=0', 'm. zerowe met. bisekcji'},'Location','southwest');
title('wykres funkcji f(x)=2,3*sin(x)+4*ln(x+2)-11');
%% 
% metoda siecznych

x_sieczne = met_siecznych(przedzialy);
x_sieczne
y_sieczne = wartosc_funkcji(x_sieczne)
fprintf('najwiekszy blad zera przy metodzie siecznych %f\n', najwieksze_zero(x_sieczne));

figure
plot(x,y,'b-',[0 15], [0 0], 'k--', x_sieczne, wartosc_funkcji(x_sieczne), 'go');
legend({'f(x)', 'y=0', 'm. zerowe met. siecznych'},'Location','southwest');
title('wykres funkcji f(x)=2,3*sin(x)+4*ln(x+2)-11');
%% 
% *funkcje pomocnicze glowne*
% 
% metoda bisekcji

function x = bisekcja(przedzialy)
    dokladnosc_zer = 0.1;
    wielkosc_przedzialu = 0.1;
    
    ilosc_pierwiastkow = size(przedzialy,1);
    x = zeros(ilosc_pierwiastkow);
    x = wektor(x);
    for i = 1:ilosc_pierwiastkow
        c = 0;
        while (abs(wartosc_funkcji(c))>dokladnosc_zer | (przedzialy(i,2)-przedzialy(i,1))>wielkosc_przedzialu)
            [c przedzialy(i,:)] = polowienie_przedzialu(przedzialy(i,:));
        end
        x(i) = c;
    end
end
%% 
% polowienie przedzialow dla metody bisekcji

function [c nowy_przedzial] = polowienie_przedzialu(przedzial)
    c = (przedzial(1)+przedzial(2))/2;
    if(sprawdzenie_przedzialu([przedzial(1) c]) == 1)
        nowy_przedzial = [przedzial(1) c]; 
    else
        nowy_przedzial = [c przedzial(2)];
    end
end
%% 
% metoda siecznych

function x = met_siecznych(przedzialy)
    dokladnosc_zer = 0.1;
    wielkosc_przedzialu = 0.1;

    ilosc_pierwiastkow = size(przedzialy,1);
    x = zeros(ilosc_pierwiastkow);
    x = wektor(x);
    for i = 1:ilosc_pierwiastkow
        c = 0;
        while (abs(wartosc_funkcji(c))>dokladnosc_zer | (przedzialy(i,2)-przedzialy(i,1))>wielkosc_przedzialu)
            [c przedzialy(i,:)] = nowy_sieczny_przedzial(przedzialy(i,:), c);
        end
        x(i) = c;
    end
end
%% 
% zmniejszenie przedzialow dla metody siecznych

function [d nowy_przedzial] = nowy_sieczny_przedzial(przedzial, c)
    d = wyznacz_zero_f_liniowej(przedzial);
    if(przedzial(2) == c)
        nowy_przedzial = [d przedzial(2)]; 
    else
        nowy_przedzial = [przedzial(1) d]; 
    end
end
%% 
% wyznaczenie zera f. liniowej wyznaczonej na podstawie punktow z koncow 
% przedzialu na podstawie wyliczonego analitycznie wzoru

function c = wyznacz_zero_f_liniowej(przedzial)
    x1 = przedzial(1);
    x2 = przedzial(2);
    y1 = wartosc_funkcji(x1);
    y2 = wartosc_funkcji(x2);
    c = (x2*y1-x1*y2)/(y1-y2);
end
%% 
% *funkcje pomocnicze dodatkowe*
% 
% wyznaczenie wyjsc dla podanych x-ow i zadanej funkcji

function y = wartosc_funkcji(x)
    [temp rozmiar_x] = size(x);
    if rozmiar_x == 1
        x = x';
        rozmiar_x = temp;
    end
    y = zeros(rozmiar_x);
    y = y(:,1);
    for i = 1:rozmiar_x
        y(i,1) = 2.3*sin(x(1,i))+4*log(x(1,i)+2)-11;
    end
end
%% 
% sprawdzenie czy w podanym przedziale jest miejsce zerowe

function result = sprawdzenie_przedzialu(przedzial)
    if wartosc_funkcji(przedzial(1))*wartosc_funkcji(przedzial(2)) < 0
        result = 1;
    else
        result = 0;
    end
end
%% 
% wartosc najwiekszego bledu

function y = najwieksze_zero(x)
    y = max(abs(wartosc_funkcji(x)));
end
%% 
% funkcja wektoryzujaca macierz diagonalna

function w = wektor(A)
    rozmiar = size(A);
    for i = 1:rozmiar
        w(i,1) = A(i,i);
    end
end
%% 
% 
% 
% <https://www.matemaks.pl/program-do-rysowania-wykresow-funkcji.html https://www.matemaks.pl/program-do-rysowania-wykresow-funkcji.html>
% 
% 2,3*sin(x)+4*log(x+2)-11
% 
%