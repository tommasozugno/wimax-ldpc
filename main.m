clear all;
close all;

save_results = 0;
load matrices/96.33.964.mat

%Parameters  **************************************************************
Nit = 20; %Number of iterations on the graph
Max_npck = 1000; %Maximum number of packets
Th_err = 100; %Error threshold
SNR_dB = 1: 1 : 5; %SNR range in dB
SNR = 10.^(SNR_dB/10); %Linear SNR range
sigmaw = sqrt(1./SNR); %Noise variance range
%**************************************************************************

%Initialization
err = zeros(1,length(SNR_dB));
Npck = zeros(1,length(SNR_dB));

tic
for npck = 1 : Max_npck
    
    %Generate info data
    u = randi([0,1],[k,1]);
    
    %Encoding
    c = [u ; mod(A*u,2)];
    
    %Modulation
    c_mod = c*2-1;
    
    %Generate noise samples
    w = randn(n,1);
    
    for snr = 1 : length(SNR_dB)
        
        if(err(snr) < Th_err)
            %Received vector
            r = c_mod + sigmaw(snr)*w;

            %Decoding
            u_hat = decode2(r, sigmaw(snr), H, k, Nit);

            %Count errors
            err(snr) = err(snr) + sum(u_hat ~= u);

            %Count packets
            Npck(snr) = Npck(snr) + 1;
        end
    end
end

%Time per packet
time = toc/sum(Npck);

%BER
Pbit = err./(Npck*k);

%Warning flag = 1 if err < th_err
warn = err < Th_err;

%% - Results

%save results
if(save_results)
    save(strcat('results/',datestr(clock)),'Pbit','SNR_dB','warn','Npck','err','Nit','time','Max_npck','Th_err');
end

%Uncoded BER
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB,Pbit,'k-',SNR_dB,Pbit_uncoded,'b--',SNR_dB,Pbit.*warn,'rx','LineWidth',2)
axis([min(SNR_dB) max(SNR_dB) 1e-7 1e0])
hleg = legend('Simulation','Uncoded BER');
set(hleg,'position',[0.15 0.13 0.32 0.15]);
xlabel('$E_b/N_0$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');