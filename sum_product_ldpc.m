clear all;
close all; 

Nit = 1;

H = [1 0 1 0 1 0 1 0 ; 1 0 0 1 0 1 0 1 ; 0 1 1 0 0 1 1 0 ; 0 1 0 1 1 0 0 1];
r = [-2.5467 0.2358 -1.3929 -3.0287 -1.8290 -1.1768 -1.9434 -0.1152];

g = -2*r;

%initialization
c_hat = zeros(1,length(r));
mu_hf = H.';
mu_fh = H;
mu_hg = zeros(1,length(g));
for i = 1 : size(mu_hf,1)
    mu_hf(i,:) = mu_hf(i,:)*g(i);
end

for it = 1 : Nit
    %check nodes update
    for i = 1 : size(mu_fh,1)
        for j = 1 : size(mu_fh,2)

            if(mu_fh(i,j) ~= 0)
                tmp1 = 0; tmp2 = 1;
                for l = 1 : size(mu_hf,1)
                    if(l ~= j && mu_hf(l,i) ~= 0)
                        tmp1 = tmp1 + phy_tilde(abs(mu_hf(l,i)));
                        tmp2 = tmp2*sign(mu_hf(l,i));
                    end
                end
                mu_fh(i,j) = phy_tilde(tmp1)*tmp2;
            end

        end
    end

    %variable nodes update
    for i = 1 : size(mu_hf,1)
        for j = 1 : size(mu_hf,2)

            if(mu_hf(i,j) ~= 0)
                tmp1 = 0;
                for l = 1 : size(mu_fh,1)
                    if(l ~= i)
                        tmp1 = tmp1 + mu_fh(l,i);
                    end
                end
                mu_hf(i,j) = tmp1 + g(i);
            end
        end
    end

    for i = 1 : length(g)
        mu_hg(i) = sum(mu_hf(i,:));
    end

    %marginalization
    for i = 1 : length(g)
        if(mu_hg(i)+g(i) >= 0)
            c_hat(i) = 0;
        else
            c_hat(i) = 1;
        end
    end
    c_hat
end




function y = phy_inv(x)
    y = tanh(0.5*x);
end

function y = phy_tilde(x)
    y = -log(phy_inv(x));
end