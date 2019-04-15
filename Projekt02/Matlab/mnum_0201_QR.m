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
A = [1 1 9;2 -1 0; -2 4 0];
%[q r] = qr(A)
[Q R] = gram_schmidt2(A)
[Q R] = gram_schmidt(A)
orth(A)
%% 
% wyswietlenie bledow

% blad_res_a
%% 
% *funkcje pomocnicze*
% 
% algorytm obliczania wartosci wlasnych metoda QR bez przesuniec

function qrb = qr_bezprzesuniec(A)
     rozmiar = size(A);
     Q = zeros(rozmiar);
     % 1
     Q(:,1) = A(:,1);
     % 2
     x = Q(:,1).*A(:,2)/(Q(:,1).*Q(:,1))
     Q(:,2) = A(:,2) - x(:,2).*Q(:,1);
     % 3
     x = x + Q(:,2)/(Q(:,2).*Q(:,2));
     Q(:,3) = A(:,3) - x(:,2).*A(:,3).*Q(:,1) - x(:,3).*A(:,3).*Q(:,2);
%      for i = 1:rozmiar
%          B[:,1]
%      end
    qrb = Q;
end
%% 
% algorytm obliczania wartosci wlasnych metoda QR z przesunieciami

function qrp = qr_przesuniecia(A)
    qrp = 0;
end
%% 
% algorytm Grama-Schmidta 2

function [Q R] = gram_schmidt2(A)
     rozmiar = size(A);
     Q = zeros(rozmiar);
     R = eye(rozmiar);
     % 1
     Q(:,1) = A(:,1);
     % 2
     R(1,2) = mydot(Q(:,1),A(:,2))/mydot(Q(:,1),Q(:,1));
     Q(:,2) = A(:,2) - R(1,2)*Q(:,1);
     % 3
     R(1,3) = mydot(Q(:,1),A(:,3))/mydot(Q(:,1),Q(:,1));
     R(2,3) = mydot(Q(:,2),A(:,3))/mydot(Q(:,2),Q(:,2));
     Q(:,3) = A(:,3) - R(1,3)*Q(:,1) - R(2,3)*Q(:,2);
%      for i = 1:rozmiar
%          B[:,1]
%      end
end
%% *algorytm Grama-Schmidta*

function [Q R] = gram_schmidt(A)
     rozmiar = size(A);
     Q = zeros(rozmiar);
     R = eye(rozmiar);
     
%      % 1
%      Q(:,1) = A(:,1);
%      % 2
%      R(1,2) = mydot(Q(:,1),A(:,2))/mydot(Q(:,1),Q(:,1));
%      Q(:,2) = A(:,2) - R(1,2)*Q(:,1);
%      % 3
%      R(1,3) = mydot(Q(:,1),A(:,3))/mydot(Q(:,1),Q(:,1));
%      R(2,3) = mydot(Q(:,2),A(:,3))/mydot(Q(:,2),Q(:,2));
%      Q(:,3) = A(:,3) - R(1,3)*Q(:,1) - R(2,3)*Q(:,2);

     for i = 1:rozmiar
         Q(:,i) = A(:,i);
         for j = 1:(i-1)             
             R(j,i) = mydot(Q(:,j),A(:,i))/mydot(Q(:,j),Q(:,j));
             Q(:,i) = Q(:,i) - R(j,i)*Q(:,j);
         end
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