function u_hat = decodeBICM(r,sigmaw,H,k,Nit,Q,C,d)

    m = size(H,1); %number of check nodes
    n = size(H,2); %number of variable nodes

    %initialization
    
    %leaf nodes
    g = zeros(2^Q,n/Q);
    for i = 1 : n/Q
        for j = 1 : 2^Q
            g(j,i) = exp(-(abs(r(i)-d(j))^2) / (2*sigmaw^2));
        end
    end
    
    c_hat = zeros(n,1); %estimated codeword
    
    %messages from variable to conform
    mu_hw = zeros(n, n/Q);
    for i = 1 : n/Q
        mu_hw(2*i-1:2*i, i) = [1 1];
    end
    
    %messages from conform to variable
    mu_wh = zeros(n/Q, n);
     
    mu_hf = H.'; %messages from variable to check
    mu_fh = H; %messages from check to variable
    
    it = 0; stopp = 0;
    while(it < Nit && stopp == 0)
        
        %Messages from conform to variable
        for i = 1 : n/Q
            for q = 1 : Q
                tmp = zeros(2^Q,1);
                for p = 1 : Q
                    if(p~=q)
                        tmp = tmp + g(:,i).*(~C(:,p)*exp(mu_hw(2*(i-1)+p,i))+C(:,p));
                    end
                end
                mu_wh(i,2*(i-1)+q) = log( sum(tmp(~C(:,q))) / sum(tmp(C(:,q))) );
            end
        end
        
        %messages from variable to check
        for i = 1 : size(mu_hf,1)
            for j = 1 : size(mu_hf,2)

                if(mu_hf(i,j) ~= 0)
                    tmp1 = 0;
                    for l = 1 : size(mu_fh,1)
                        if(l ~= j)
                            tmp1 = tmp1 + mu_fh(l,i);
                        end
                    end
                    mu_hf(i,j) = tmp1 + mu_wh(ceil(i/Q),i);
                end
            end
        end
        
        %check nodes update
        tmp1 = phy_tilde(abs(mu_hf));
        tmp3 = ((sum(tmp1,1).')*ones(1,n)).*H-tmp1.';
        tmp4 = (mu_hf>=0)*2-1;
        tmp5 = prod(tmp4,1).';
        mu_fh =  phy_tilde(tmp3).*(tmp5*ones(1,n).*tmp4.'); %messages from check to variable

        % Messages from variable to conform, bit
        for i = 1 : n
            mu_hw(i,ceil(i/Q)) = sum(mu_fh(:,i));
        end
        
        %marginalization
        for i = 1 : n
           if(mu_hw(i,ceil(i/Q)) + mu_wh(ceil(i/Q),i) >= 0)
               c_hat(i) = 0;
           else
               c_hat(i) = 1;
           end
        end
       
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