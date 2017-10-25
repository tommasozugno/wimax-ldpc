function u_hat = decode2(r,sigmaw,H,k,Nit)

    m = size(H,1); %number of check nodes
    n = size(H,2); %number of variable nodes

    g = -2*r/(sigmaw^2); %LLR leaf nodes

    %initialization
    c_hat = zeros(n,1); %estimated codeword
    mu_hf = H.'; %messages from variable to check
    for i = 1 : size(mu_hf,1)
        mu_hf(i,:) = mu_hf(i,:)*g(i);
    end
    
    it = 0; stopp = 0;
    while(it < Nit && stopp == 0)
        
        %check nodes update
        tmp1 = phy_tilde(abs(mu_hf));
        tmp3 = ((sum(tmp1,1).')*ones(1,n)).*H-tmp1.';
        tmp4 = (mu_hf>=0)*2-1;
        tmp5 = prod(tmp4,1).';
        mu_fh =  phy_tilde(tmp3).*(tmp5*ones(1,n).*tmp4.'); %messages from check to variable

        %variable nodes update
        tmp = sum(mu_fh).' + g;
        mu_hf = (tmp*ones(1,m)).*(H.') - mu_fh.';
        
        mu_hg = sum(mu_fh,1); %messages from variable to leaf

        %marginalization
        c_hat = (mu_hg.'+g)<0;
       
        if(sum(mod(H*c_hat,2)) == 0)
            stopp = 1;
        end
        
        it = it + 1;
    end
    u_hat = c_hat(1:k);
end

function y = phy_tilde(x)

    ind = x>0 & x<10^-5;
    ind2 = x>=10^-5 & x<50;
    ind3 = x>=50;
    y = sparse(size(x,1),size(x,2));
    y(ind) = 12.5;
    y(ind2) = -log(tanh(0.5*x(ind2)));
    y(ind3) = 0;    
end