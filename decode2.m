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
        for i = 1 : size(mu_fh,1)
            for j = 1 : size(mu_fh,2)

                if(mu_fh(i,j) ~= 0)
                    tmp1 = 0; tmp2 = 1;
                    for l = 1 : size(mu_hf,1)
                        if(l ~= j && mu_hf(l,i) ~= 0)
                            tmp1 = tmp1 + phy_tilde(abs(mu_hf(l,i)));
                            tmp2 = tmp2*sign(mu_hf(l,i));
                        end
                    end
                    mu_fh(i,j) = phy_tilde(tmp1)*tmp2;
                end

            end
        end

        %variable nodes update
        for i = 1 : size(mu_hf,1)
            for j = 1 : size(mu_hf,2)

                if(mu_hf(i,j) ~= 0)
                    tmp1 = 0;
                    for l = 1 : size(mu_fh,1)
                        if(l ~= i)
                            tmp1 = tmp1 + mu_fh(l,i);
                        end
                    end
                    mu_hf(i,j) = tmp1 + g(i);
                end
            end
        end

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

function y = phy_inv(x)
    y = tanh(0.5*x);
end

function y = phy_tilde(x)
    y = -log(phy_inv(x));
end