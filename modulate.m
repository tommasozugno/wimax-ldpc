function y = modulate(c, Q)

    if(Q == 1) %BPSK
        y = c*2-1;
    else %QPSK, 16QAM, 64QAM
        for i = 1 : length(c)/Q
            tmp(i) = bi2de(c(Q*(i-1)+1:Q*i),'left-msb');
        end
        y = qammod(tmp,2^Q,'UnitAveragePower',true); %constellation energy 1

    end
    
end   