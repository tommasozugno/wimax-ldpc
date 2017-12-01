function y = modulate(c, useQPSK)

    if(useQPSK)
        %QPSK
        if(mod(length(c),2) ~= 0)
            display('Error: c length is not multiple of 2');
            c = [c , zeros(1,2-mod(length(c),2))]; %zero padding
        end

        c_mod = zeros(length(c)/2,1);
        for i = 0 : length(c)/2-1
            if(c(2*i+1) == 0)
                c_mod(i+1) = 1;
            else
                c_mod(i+1) = -1;
            end

            if(c(2*(i+1)) == 0)
                c_mod(i+1) = c_mod(i+1) + 1i;
            else
                c_mod(i+1) = c_mod(i+1) - 1i;
            end
        end
        y = c_mod;
    else
        %BPSK
        y = c*2-1;
    end
end   