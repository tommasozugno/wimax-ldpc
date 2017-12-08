clear all;
close all;

save_results = 0;

%********************
%   rate 1 --> 1/2  *
%   rate 2 --> 2/3B *
%   rate 3 --> 3/4A *
%   rate 4 --> 3/4B *
%   rate 5 --> 5/6  *
%********************

rate = 5; %Code rate
n = 576;  %Codeword length

load(strcat('matrices/r',num2str(rate),'n',num2str(n)));
%load matrices/96.33.964.mat
H = sparse(H);

M1 = double(M1.x);
M2 = double(M2.x);
M3 = double(M3.x);

%Parameters  **************************************************************
Nit = 50; %Number of iterations on the graph
Max_npck = 1000; %Maximum number of packets
Th_err = 100; %Error threshold
SNR_dB = [13 15]; %SNR range in dB
SNR = 10.^(SNR_dB/10); %Linear SNR range
sigmaw = sqrt(1./SNR); %Noise variance range
Q = 6; %Modulation order
%**************************************************************************

%Initialization
err = zeros(1,length(SNR_dB));
Npck = zeros(1,length(SNR_dB));
time = zeros(1,length(SNR_dB));

%Arguments for decodeBICM
if(Q > 1)
    C = logical(de2bi(0:2^Q-1,'left-msb'));
    for i = 1 : size(C,1)
        d(i) = modulate(C(i,:),Q);
    end
end

for npck = 1 : Max_npck
    
    if(~all(err >= Th_err) ) %to stop the execution if the min number
                             %of errors is reached for all SNRs
    
        %Generate info data
        u = randi([0,1],[k,1]);

        %Encoding
        p1t = mod(M1*u,2);
        p2t = mod(M2*u+M3*p1t,2);
        c = [u ; p1t ; p2t];

        %Modulation
        c_mod = modulate(c.', Q);


        %Generate noise samples
        if(Q > 1)
            %QPSK, 16QAM, 64QAM
            w = randn(size(c_mod)) + 1i*randn(size(c_mod));
        else
            %BPSK
            w = randn(n,1);
        end

        for snr = 1 : length(SNR_dB)

            if(err(snr) < Th_err)

                if(Q > 1)
                    %QPSK, 16QAM, 64QAM
                    %Received vector
                    r = c_mod + w*sigmaw(snr);
                    %Decoding
                    u_hat = decodeBICM(r, sigmaw(snr), H, k, Nit, Q, C, d);
                else
                    %BPSK
                    %Received vector
                    r = c_mod + sigmaw(snr)*w;
                    %Decoding
                    u_hat = decode2(r, sigmaw(snr), H, k, Nit);
                end

                %tic

                %time(snr) = time(snr) + toc;

                %Count errors
                err(snr) = err(snr) + sum(u_hat ~= u);

                %Count packets
                Npck(snr) = Npck(snr) + 1;
            end
        end
    end
end

%Average time to decode a packet
%time = time./Npck;

%BER
Pbit = err./(Npck*k);

%Warning flag = 1 if err < th_err
warn = err < Th_err;

%% - Results

%save results
if(save_results)
    save(strcat('results/',datestr(clock)),'Pbit','SNR_dB','Npck','err','Nit','time','Max_npck','Th_err','rate','n');
end

if(false)
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
end