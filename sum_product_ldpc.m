%Reference: http://webstaff.itn.liu.se/~eribe/tne066/lab/ldpc_harq/Tutorial%20%96%20the%20sum-product%20algorithm.pdf
clear all;
close all; 

%% Encoder
n = 8; k = 4; %Codeword length and info bits
H = [1 0 1 0 1 0 1 0 ; 1 0 0 1 0 1 0 1 ; 0 1 1 0 0 1 1 0 ; 0 1 0 1 1 0 0 1]; %Parity check matrix
G = [eye(k) ; 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 1 1 0 1]; %Generetor matrix
u = randi([0,1],[k,1]); %Word
c = mod(G*u , 2); %Codeword

%% BPSK
c_mod = 2*c-1; %Modulate

%% Channel
SNR_dB = 0;
SNR = 10^(SNR_dB/10);
sigmaw = sqrt(1/SNR); %Noise variance

w = randn(size(c_mod)); %AWGN noise
r = c_mod + w*sigmaw; %Received vector
%r = [-2.5467 0.2358 -1.3929 -3.0287 -1.8290 -1.1768 -1.9434 -0.1152];

%% Decoder
Nit = 20; %Number of iterations on the graph
g = -2*r/(sigmaw^2); %LLR leaf nodes

%initialization
c_hat = zeros(length(r),1); %estimated codeword
mu_hf = H.'; %messages from variable to check
mu_fh = H; %messages from check to variable
mu_hg = zeros(1,length(g)); %messages from variable to leaf
for i = 1 : size(mu_hf,1)
    mu_hf(i,:) = mu_hf(i,:)*g(i);
end

it = 0; stopp = 0;
while(it < Nit && stopp == 0)
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
        mu_hg(i) = sum(mu_fh(:,i));
    end

    %marginalization
    for i = 1 : length(g)
        if(mu_hg(i)+g(i) >= 0)
            c_hat(i) = 0;
        else
            c_hat(i) = 1;
        end
    end
    
    if(sum(mod(H*c_hat,2)) == 0)
        %stopp = 1;
    end
    
    it = it + 1;
end

u_hat = c_hat(1:k); %Estimated message
err = sum(u_hat ~= u);


function y = phy_inv(x)
    y = tanh(0.5*x);
end

function y = phy_tilde(x)
    y = -log(phy_inv(x));
end