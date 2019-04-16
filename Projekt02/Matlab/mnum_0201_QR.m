%% 
% program

clc;
clear;
% for rozmiar = [5 10 20]
%     for i = 1:30
%         macierz_symetryczna(rozmiar);
%         macierz_niesymetryczna(rozmiar);
%     end
% end
%A = [1 1 9;2 -1 0; -2 4 0];
%A = [1 1;2 -1; -2 4];
%A = [1 1 2; -1 -2 4];
A = macierz_symetryczna(5);
    
[q r] = qr_rozklad_niekwadrat(A)
[Q R] = qr(A)
%eig(A)
%qr_bezprzesuniec(A)
%% 
% wyswietlenie bledow

% blad_res_a
%% 
% *funkcje pomocnicze*
% 
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
%% rozklad QR dla macierzy niekwadratowych

function [Q R] = qr_rozklad_niekwadrat(A)
     [r_wiersze r_kolumny] = size(A);
     Q = zeros(r_wiersze);
     if r_wiersze > r_kolumny
        R = eye(r_wiersze);
     else
        R = eye(r_kolumny);
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
% rozklad QR dla macierzy kwadratowych

function [Q R] = qr_rozklad(A)
     rozmiar = size(A);
     Q = zeros(rozmiar);
     R = eye(rozmiar);
    %Gram-Schmidt
     for i = 1:rozmiar
         Q(:,i) = A(:,i);
         for j = 1:(i-1)
             R(j,i) = mydot(Q(:,j),A(:,i))/mydot(Q(:,j),Q(:,j));
             Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
         end
     end
    %normalizacja
     N = zeros(rozmiar);
     for i = 1:rozmiar
         N(i,i) = norm(Q(:,i));
         Q(:,i) = Q(:,i)/N(i,i);
     end
     R = N*R;
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