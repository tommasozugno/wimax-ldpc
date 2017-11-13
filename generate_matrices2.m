clear all;
close all;

%********************
%   rate 1 --> 1/2  *
%   rate 2 --> 2/3B *
%   rate 3 --> 3/4A *
%   rate 4 --> 3/4B * <- Error phy is not an identity
%   rate 5 --> 5/6  *
%********************

%Parameters*************************
rate = 5; %Code rate               *
n = 576; %Codeword length 2304     *
%***********************************

switch rate
    case 1, R = 1/2; load matrices/Hb_r1_2.mat;
    case 2, R = 2/3; load matrices/Hb_r2_3B.mat;
    case 3, R = 3/4; load matrices/Hb_r3_4A.mat;
    %case 4, R = 3/4; load matrices/Hb_r3_4B.mat;
    case 5, R = 5/6; load matrices/Hb_r5_6.mat;
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

%Find the submatrices
A = gf(H(1:m-z , 1:k),1);
B = gf(H(1:m-z , k+1:k+z),1);
T = gf(H(1:m-z , k+z+1:end),1);
C = gf(H(m-z+1:end , 1:k),1);
D = gf(H(m-z+1:end , k+1:k+z),1);
E = gf(H(m-z+1:end , k+z+1:end),1);

%test
phy = E*T^-1*B+D;
if(~isequal(phy,gf(eye(z),1)))
    disp('Error')
end

M1 = E*T^-1*A+C;
M2 = T^-1*A;
M3 = T^-1*B;

%test
u = randi([0 1],k,1);

p1t = M1*u;
p2t = M2*u+M3*p1t;

c = [u ; p1t ; p2t];

if(~isequal(gf(H,1)*c,gf(zeros(m,1),1)))
    display('ERROR: matrix G is not correct');
end

save(strcat('matrices/r',num2str(rate),'n',num2str(n),'.mat'),'M1','M2','M3','H','n','k');

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

