%% Confronto Nit, n = 576 r = 5
clear all;
close all;

load results/Nit/Copy_of_Nit1.mat;
Pbit1 = Pbit;
SNR_dB1 = SNR_dB;

load results/Nit/Copy_of_Nit10.mat;
Pbit2 = Pbit;
SNR_dB2 = SNR_dB;

load results/Nit/Copy_of_Nit50.mat;
Pbit3 = Pbit;
SNR_dB3 = SNR_dB;

load results/Nit/Copy_of_Nit100.mat;
Pbit4 = Pbit;
SNR_dB4 = SNR_dB;


%Uncoded BER
SNR_dBU = [5 : 7 , 7.2];
SNR = 10.^(SNR_dBU/10);
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB1,Pbit1,'k-',SNR_dB2,Pbit2,'rs-',SNR_dB3,Pbit3,'go-',SNR_dB4,Pbit4,'b+-',SNR_dBU,Pbit_uncoded,'b--','LineWidth',2,'MarkerSize',10)
axis([5 7.2 1e-7 1e-1])
hleg = legend('Nit = 1, n = 576, r = 5/6',...
              'Nit = 10, n = 576, r = 5/6',...
              'Nit = 50, n = 576, r = 5/6',...
              'Nit = 100, n = 576, r = 5/6','Uncoded BER');
set(hleg,'position',[0.15 0.13 0.32 0.15]);
xlabel('$E_b/N_0$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');
    
%% - Confronto rate, n = 576 Nit = 50

clear all;

load results/r/Copy_of_r1.mat;
Pbit1 = Pbit;
SNR_dB1 = SNR_dB;

load results/r/Copy_of_r2.mat;
Pbit2 = Pbit;
SNR_dB2 = SNR_dB;

load results/r/Copy_of_r3.mat;
Pbit3 = Pbit;
SNR_dB3 = SNR_dB;

load results/r/Copy_of_r5.mat;
Pbit5 = Pbit;
SNR_dB5 = SNR_dB;

%Uncoded BER
SNR_dBU = [1 : 7 , 7.2];
SNR = 10.^(SNR_dBU/10);
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB1,Pbit1,'k-',SNR_dB2,Pbit2,'go-',SNR_dB3,Pbit3,'b+-',SNR_dB5,Pbit5,'rs-',SNR_dBU,Pbit_uncoded,'b--','LineWidth',2,'MarkerSize',10)
axis([1 7.2 1e-5 1e-1])
hleg = legend('R = 1/2, n = 576, N_{it} = 50'...
              ,'R = 2/3 B, n = 576, N_{it} = 50'...
              ,'R = 3/4 A, n = 576, N_{it} = 50'...
              ,'R = 5/6, n = 576, N_{it} = 50'...
              ,'Uncoded BER');
set(hleg,'position',[0.15 0.13 0.32 0.15]);
xlabel('$E_b/N_0$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');
    
%% - Confronto n, R=3 Nit = 50

clear all;

load results/n/n576;
Pbit1 = Pbit;
SNR_dB1 = SNR_dB;

load results/n/Copy_of_n1344_prova.mat;
Pbit2 = Pbit;
SNR_dB2 = SNR_dB;

load results/n/Copy_of_n2304.mat;
Pbit3 = Pbit;
SNR_dB3 = SNR_dB;

%Uncoded BER
SNR_dBU = [1 : 7 , 7.2];
SNR = 10.^(SNR_dBU/10);
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB1,Pbit1,'k-',SNR_dB2,Pbit2,'go-',SNR_dB3,Pbit3,'rs-',SNR_dBU,Pbit_uncoded,'b--','LineWidth',2,'MarkerSize',10)
axis([3 6 1e-5 1e-1])
hleg = legend('n = 576, r = 3/4 A, Nit = 50'...
              ,'n = 1344, r = 3/4 A, Nit = 50'...
              ,'n = 2304, r = 3/4 A, Nit = 50','Uncoded BER');
set(hleg,'position',[0.15 0.13 0.32 0.15]);
xlabel('$E_b/N_0$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');
    
    
%% - Confronto BICM, n = 576 R = 5

clear all;

load results/r/Copy_of_r5.mat;
Pbit1 = Pbit;
SNR_dB1 = SNR_dB;

% load results/BICM/qpsk.mat;
% Pbit2 = Pbit;
% SNR_dB2 = SNR_dB;

load results/BICM/qpsk.mat;
Pbit3 = Pbit;
SNR_dB3 = SNR_dB;


%Uncoded BER
SNR_dBU = [1 : 7 , 7.2];
SNR = 10.^(SNR_dBU/10);
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB1,Pbit1,'k-',SNR_dB3,Pbit3,'ro-',SNR_dBU,Pbit_uncoded,'b--','LineWidth',2,'MarkerSize',10)
axis([1 10 1e-5 1])
hleg = legend('BPSK','QPSK');
set(hleg,'position',[0.15 0.13 0.32 0.15]);
xlabel('$E_b/N_0$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');
    


    

