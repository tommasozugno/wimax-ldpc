%Script that generates the matrices G and H given the corresponding alist
%file. The matrix G in systematic form is obtained by using the gauss-
%eliminaition method.

clear all;
close all;

load matrices/alist96.33.964.mat;

H = alisttoH(alist);
H = H.';
n = alist(1,1); m = alist(1,2); k = alist(1,3);

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
H_prime = [A eye(n-k)];

%test
if(H*G_prime ~= zeros(m,m))
    display('ERROR: matrix G is not correct 1');
end

u = randi([0 1],k,1);
c = mod(G_prime*u,2);

if(mod(H*c,2) ~= zeros(m,1))
    display('ERROR: matrix G is not correct');
end

save('matrices/96.33.964.mat','H','A','n','k');


%Reference: http://www.inference.org.uk/mackay/codes/alist.html
%Function to obtain H from the corresponding alist file.
%You have to remove the 3th and the 4th row from the original alist file.
function H = alisttoH(alist)

    H = zeros(alist(1,1),alist(1,2));
    for i = 1 : size(H,1)
        for j = 1 : alist(2,1)
            H(i,alist(i+2,j)) = 1;
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