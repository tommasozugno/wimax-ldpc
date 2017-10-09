clear all;
close all;

save_results = 1;

load matrices/96.33.964.mat

Nit = 20; %Number of iterations on the graph

%Generate data
Npck = 100;
u = randi([0,1],[Npck*k,1]);

%Encode
c = [];
for i = 0 : Npck-1
    c = [c ; mod(G*u(k*i+1:(i+1)*k) , 2)];
end

%Modulate
c_mod = 2*c-1;

%Channel
SNR_dB = 1: 1 : 5;
SNR = 10.^(SNR_dB/10);
sigmaw = sqrt(1./SNR); %Noise variance
w = randn(length(c_mod),1);

tic
for snr = 1 : length(SNR)
    r = c_mod + w*sigmaw(snr);
    
    u_hat = [];
    for npck = 0 : Npck-1
        u_hat = [u_hat ; decode2(r(npck*n+1 : (npck+1)*n),sigmaw(snr),H,k,Nit)];
    end
    
    err(snr) = sum(u_hat ~= u);
end
toc

Pbit = err/(Npck*k);

%% - Results

%save results
if(save_results)
    save(strcat('results/',datestr(clock)),'Pbit','SNR_dB','Npck','Nit');
end

%Uncoded BER
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB,Pbit,'k-',SNR_dB,Pbit_uncoded,'b--','LineWidth',2)
axis([min(SNR_dB) max(SNR_dB) 1e-7 1e0])
hleg = legend('Simulation','Uncoded BER');
set(hleg,'position',[0.15 0.13 0.32 0.15]);
xlabel('$E_b/N_0$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');