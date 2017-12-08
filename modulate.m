function y = modulate(c, Q)

    switch Q
        
        case 1  %BPSK
            y = c*2-1;
            
        case 2 %QPSK
            for i = 1 : length(c)/2
                tmp(i) = bi2de(c(2*i-1:2*i),'left-msb');
            end
            y = qammod(tmp,4);
            
        case 4 %16QAM
            for i = 1 : length(c)/4
                tmp(i) = bi2de(c(4*i-3:4*i),'left-msb');
            end
            y = qammod(tmp,16,'UnitAveragePower',true)*sqrt(4);
            
        case 6 %64QAM
            for i = 1 : length(c)/6
                tmp(i) = bi2de(c(6*i-5:6*i),'left-msb');
            end
            y = qammod(tmp,64,'UnitAveragePower',true)*sqrt(6);
            
    end
end   