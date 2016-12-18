function pchange = pexeest_p(times, cens, tchange, mono)

% a function produces the P-value at the changepoints "tchange"

tchange = sort(tchange);
nchange = length(tchange);
[time_die,ttot,deaths] = totaltest(times,cens);
ntime   = length(time_die);


if nchange >= 1,
    % find the index
    indchange = zeros(nchange,1);
    for j = 1:nchange,
        for i = 1:ntime,
            if abs(time_die(i)-tchange(j)) < 0.00001, % for round off error
                indchange(j) = i;
            end
        end
    end

    if length(unique(indchange)) < length(indchange),
        % there is a need to look for unique variables
        [indchange, I, J] = unique(indchange,'first'); 
        tchange = tchange(I);
        nchange = length(tchange);
    end
    
    % estimate the piecewise exponential parameter 
    %  lambda1-lambda_nchange
    % compute all the changepoint pvalues
    lamest  = zeros(nchange+1,1);
    E = zeros(nchange,1);
    for j = 1:nchange,
        if j == 1,
            ttot1 = sum(ttot(1:indchange(1)));
            d1    = sum(deaths(1:indchange(1)));
            ttot2 = sum(ttot(indchange(1)+1:indchange(2)));
            d2    = sum(deaths(indchange(1)+1:indchange(2)));
            [a2,pval] = exact_pvalue(ttot1,ttot2,d1,d2,mono);
            pchange(1) = pval;
        elseif j == nchange,
            ttot1 = sum(ttot(indchange(j-1)+1:indchange(j)));
            d1 = sum(deaths(indchange(j-1)+1:indchange(j)));
            ttot2 = sum(ttot(indchange(j)+1:end));
            d2 = sum(deaths(indchange(j)+1:end));    
            [a2,pval] = exact_pvalue(ttot1,ttot2,d1,d2,mono);
            pchange(j) = pval;   
        else
            ttot1 = sum(ttot(indchange(j-1)+1:indchange(j)));
            d1 = sum(deaths(indchange(j-1)+1:indchange(j)));
            ttot2 = sum(ttot(indchange(j)+1:indchange(j+1)));
            d2 = sum(deaths(indchange(j)+1:indchange(j+1)));           
            [a2,pval] = exact_pvalue(ttot1,ttot2,d1,d2,mono);
            pchange(j) = pval;       
        end
    end

%     for k = 1:length(tx),
%         if tx(k) < tchange(1), % the first piece
%             quan(k) = exp(-tx(k)/lamest(1));
%         elseif tx(k) < tchange(nchange), % piece 2 -- nchange 
%             for j = 2:nchange,
%                 if tx(k) >= tchange(j-1),
%                     if tx(k) < tchange(j),
%                         quan(k) = E(j-1)*exp(-(tx(k)-tchange(j-1))/lamest(j));
%                     end
%                 end
%             end
%         else % piece nchange+1
%             quan(k) = E(nchange)* exp(-(tx(k)-tchange(nchange))/...
%                 lamest(nchange+1));
%         end
%     end
else % nchange < 1
    pchange = [];
end
