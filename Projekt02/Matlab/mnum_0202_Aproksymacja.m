%% 
% program

clc;
clear;

dane = [-5 -5.4606;-4 -3.8804;-3 -1.9699;-2 -1.6666;-1 -0.0764;0 -0.3971;1 -1.0303;2 -4.5483;3 -11.528;4 -21.6417;5 -34.4458];

funkcja = uklad_rownan_normalnych(dane, 1);

x = linspace(-5,5, 100);
y = fun(funkcja, x);
plot(dane(:,1),dane(:,2),'bo', x, y, 'r-')
%A
%% 
% wyswietlenie bledow

% blad_res_a
%% 
% *funkcje pomocnicze*
% 
% uklad rownan normalnych

function wspolczynniki = uklad_rownan_normalnych(dane, st_wielomianu)
   % wyznaczanie macierzy Grama - <przeksztalcenie_i,przeksztalcenie_j>
    st_wielomianu = st_wielomianu + 1;
    macierz_Grama = zeros(st_wielomianu);
    [r_wiersze, r_kolumny] = size(dane);
    for i = 1:st_wielomianu
        for j = 1:st_wielomianu
            for k = 1:r_wiersze
                macierz_Grama(i,j) = macierz_Grama(i,j) + (dane(k,1))^(i+j-2);
            end
        end
    end
    
   % wektor prawej strony
    prawa_strona = zeros(st_wielomianu);
    prawa_strona = prawa_strona(:,1);
    for i = 1:st_wielomianu
        for k = 1:r_wiersze
            prawa_strona(i) = prawa_strona(i) + (dane(k,1))^(i-1)*dane(k,2);
        end
    end
    
    wspolczynniki = gauss(macierz_Grama, prawa_strona);
end
%% 
% wyznaczenie wyjsc dla podanych x-ow i zadanej funkcji

function y = fun(funkcja, x)
    rozmiar_x = size(x);
    st_wielomianu = size(funkcja);
    y = zeros(rozmiar_x);
    y = y(:,1);
    
    for i = 1:rozmiar_x
        for j = 1:st_wielomianu
            y(i) = y(i) + funkcja(i)*(x(i))^(j-1);
        end
    end
end
%% 
% algorytm obliczania wartosci wlasnych metoda QR bez przesuniec

function A = qr_bezprzesuniec(A)
    while tolerancja(A) > 0.00001
        [Q R] = qr_rozklad(A);
        A = R * Q;
    end
    A = wektor(A);
end
%% 
% algorytm obliczania wartosci wlasnych metoda QR z przesunieciami

function qrp = qr_przesuniecia(A)
    qrp = 0;
end
%% 
% rozklad QR dla macierzy niekwadratowych

function [Q R] = qr_rozklad_niekwadrat(A)
     [r_wiersze r_kolumny] = size(A);
     Q = zeros(r_wiersze);
     if r_wiersze > r_kolumny
        R = eye(r_wiersze);
        Q = eye(r_wiersze);
     else
        R = eye(r_kolumny);
        Q = eye(r_wiersze);
     end
    %Gram-Schmidt
     for i = 1:r_kolumny
         Q(:,i) = A(:,i);
         for j = 1:(i-1)
             R(j,i) = mydot(Q(:,j),A(:,i))/mydot(Q(:,j),Q(:,j));
             Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
         end
     end
     Q = Q(1:r_wiersze,1:r_kolumny);
     %normalizacja
     N = zeros(r_wiersze);
     for i = 1:r_kolumny
         N(i,i) = norm(Q(:,i));
         Q(:,i) = Q(:,i)/N(i,i);
     end
     R = N*R;
     
     if r_wiersze > r_kolumny
        R = R(1:r_kolumny,1:r_kolumny);
     else
        R = R(1:r_wiersze,1:r_wiersze);
     end
end
%% 
% autorska implementacja matlabowej funkcji dot()

function md = mydot(A,B)
    rozmiar = size(A);
    md = 0;
    for i = 1:rozmiar
        md = md + A(i)*B(i);
    end
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
% sprawdzenie tolerancji

function tol = tolerancja(A)
    rozmiar = size(A);
    A = abs(A);
    tol = 0;
    for i = 1:rozmiar
        if max(A(i,i+1:end)) > tol
            tol = max(A(i,i+1:end));
        end
        if max(A(i,1:i-1)) > tol
            tol = max(A(i,1:i-1));
        end
    end
end
%% 
% tworzenie macierzy symetrycznej o zadanym rozmiarze

function mac_sym = macierz_symetryczna(rozmiar)
    mac_sym = randi([0 50],rozmiar,rozmiar);
    mac_sym = mac_sym + mac_sym';
end
%% 
% tworzenie macierzy niesymetrycznej o zadanym rozmiarze

function mac_nsym = macierz_niesymetryczna(rozmiar)
    mac_nsym = randi([0 100],rozmiar,rozmiar);
end
%% 
% norma residuum

function nr = norma_residuum(wspolczynniki, x, rozw)
    residuum = wspolczynniki*x - rozw;
    nr = norm(residuum);
end
%% 
% algorytm Gaussa z projektu1

function rozw = gauss(wspolczynniki, rozw)
    rozmiar = size(rozw);
    pierwotne_b = rozw;
    pierwotne_wspolczynniki = wspolczynniki;
    for i = 1:(rozmiar-1)
        podmacierz = wspolczynniki(i:rozmiar, i:rozmiar);
        [el_glowny, ind_glowny] = max(podmacierz(:));
        [ind_wiersz, ind_kolumna] = ind2sub(size(podmacierz), ind_glowny);
        wspolczynniki([i,ind_wiersz+i-1],:) = wspolczynniki([ind_wiersz+i-1,i],:);
        pierwotne_wspolczynniki([i,ind_wiersz+i-1],:) = pierwotne_wspolczynniki([i,ind_wiersz+i-1],:);
        rozw([i,ind_wiersz+i-1],:) = rozw([ind_wiersz+i-1,i],:);
        pierwotne_b([i,ind_kolumna+i-1], :) = pierwotne_b([i,ind_kolumna+i-1], :);
        wspolczynniki(:,[i,ind_kolumna+i-1]) = wspolczynniki(:,[ind_kolumna+i-1,i]);
        pierwotne_wspolczynniki(:,[i,ind_kolumna+i-1]) = pierwotne_wspolczynniki(:,[ind_kolumna+i-1,i]);
        a = wspolczynniki(i,i);
        for l = (i+1):rozmiar
            wspolczynniki(l,i) = wspolczynniki(l,i) / a;
            wspolczynniki(l, (i+1):end) = wspolczynniki(l, (i+1):end) - wspolczynniki(l,i) * wspolczynniki(i, (i+1):end);
            rozw(l) = rozw(l) - wspolczynniki(l,i) * rozw(i);
        end
    end
    
    % U * x = rozw -> wyznaczenie x
    for i = (rozmiar:-1:1)
        % Ax=b; -> odejmowanie wartosci znanych ze strony rozwiazan
        for j = (i+1):(rozmiar)
            rozw(i) = rozw(i) - wspolczynniki((i),(j))*wspolczynniki((j),(j));
        end
        wspolczynniki(i,i) = rozw(i)/wspolczynniki(i,i);
    end
    
    %przepisanie rozwiazan do macierzy rozw
    for i = (1:rozmiar)
        rozw(i) = wspolczynniki(i,i);
    end
end