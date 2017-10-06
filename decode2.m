function u_hat = decode2(r,sigmaw,H,k,Nit)

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
        tmp1 = arrayfun(@phy_tilde,abs(mu_hf));
    
        tmp2 = sum(tmp1,1).';

        tmp3 = (tmp2*ones(1,size(H,2))).*H-tmp1.';

        tmp4 = (mu_hf>=0)*2-1;

        tmp5 = prod(tmp4,1).';

        mu_fh = arrayfun(@phy_tilde,tmp3).*(tmp5*ones(1,size(H,2)).*tmp4.');

        %variable nodes update
        tmp = sum(mu_fh).' + g;
        mu_hf = (tmp*ones(1,size(H,1))).*(H.') - mu_fh.';
        
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
    u_hat = c_hat(1:k);
end

function y = phy_tilde(x)
    %y = -log(phy_inv(x));
    if(x<=0)
        y = 0;
    else
        if(x<10^-2)
            y = 6;
        else
            y = -log(tanh(0.5*x));
        end
    end
end