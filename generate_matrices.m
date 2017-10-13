clear all;
close all;

%********************
%   rate 1 --> 1/2  *
%   rate 2 --> 2/3B *
%   rate 3 --> 3/4A *
%   rate 4 --> 3/4B *
%   rate 5 --> 5/6  *
%********************

%Parameters*************************
rate = 1; %Code rate               *
n = 576; %Codeword length 2304     *
%***********************************

switch rate
    case 1, R = 1/2; load Hb_r1_2.mat;
    case 2, R = 2/3; load Hb_r2_3B.mat;
    case 3, R = 3/4; load Hb_r3_4A.mat;
    case 4, R = 3/4; load Hb_r3_4B.mat;
    case 5, R = 5/6; load Hb_r5_6.mat;
end     
k = R*n; %Infoword length
m = n - k; %Number of parity checks
z = n/24; %Expansion factor
 
%Compute the base matrix for the given rate and codeword length
Hbm = p(Hbm0,z);

%Derive the parity check matrix
H = zeros(m,n);
for i = 0 : size(Hbm,1)-1
    for j = 0 : size(Hbm,2)-1
        H(i*z+1 : (i+1)*z , j*z+1 : (j+1)*z) = expand(Hbm(i+1,j+1),z);
    end
end

%Compute the generating matrix
B = H(: , 1 : k);
C = H(: , k+1 : end);
while(gfrank(C,2) < m)
    H_perm = H(:,randperm(n));
    C = H_perm(: , k+1 : end);
    B = H_perm(: , 1 : k);
end

C_1 = invGF2(C);
A = mod(C_1*B,2);
G_prime = [eye(k) ; A]; %generating matrix (systematic form [Ik ; -A])
H_prime = [A eye(m)];

%test
if(H*G_prime ~= zeros(m,k))
    display('ERROR: matrix G is not correct 1');
end

u = randi([0 1],k,1);
c = mod(G_prime*u,2);

if(mod(H*c,2) ~= zeros(m,1))
    display('ERROR: matrix G is not correct');
end

G = G_prime;

save(strcat('matrices/r',num2str(rate),'n',num2str(n),'.mat'),'H','A','n','k');

function Y = expand(x,z)
    if(x<0)
        Y = zeros(z);
    else
        Y = circshift(eye(z),x,1);
    end
end

function Y = p(Hbm0,z)
    
    Y = zeros(size(Hbm0));
    for i = 1 : size(Hbm0,1)
        for j = 1 : size(Hbm0,2)
            if(Hbm0(i,j)<=0)
                Y(i,j) = Hbm0(i,j);
            else
                Y(i,j) = floor(Hbm0(i,j)*z/96);
            end
        end
    end
    
end

%Find the inverse of a binary matrix using the gauss elimination method
function C_1 = invGF2(C)

    [m n] = size(C);
    if(m ~= n)
        display('The matrix is not square');
    end
    
    A = [C eye(m)];
    
    for i = 1 : m
        it = i;
        while(A(i,i) == 0 && it <= m)
            if(A(it,i) ~= 0)
                tmp = 1 : m;
                tmp(i) = it;
                tmp(it) = i;
                A = A(tmp,:);
            end
            it = it + 1;
        end

        for l = [1 : i-1 i+1 : m]
            if(A(l,i) ~= 0)
                A(l,:) = mod(A(i,:)+A(l,:),2);
            end
        end
    end
    C_1 = A(:,m+1:end);
end

