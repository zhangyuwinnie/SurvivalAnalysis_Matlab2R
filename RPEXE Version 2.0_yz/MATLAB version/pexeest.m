function [quan lamest] = pexeest(times, cens, tchange, tx)

% a function produces the the quantile estimate at tx, 
% when a piecewise exponential distribution is fitted to 
% (times,cens) cens = 0 for censored, cens = 1 for uncensored.
% the chance point is tchange
% lamest is the estimated parameters

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

    if length(unique(indchange))<length(indchange),
        % there is a need to look for unique variables
        [indchange, I, J] = unique(indchange,'first'); 
        tchange = tchange(I);
        nchange = length(tchange);
    end
    
    % estimate the piecewise exponential parameter lambda1-lambda_nchange
    % compute all the changepoint quantile
    lamest  = zeros(nchange+1,1);
    E = zeros(nchange,1);
    for j = 1:nchange,
        if j == 1,
            lamest(j) = sum(ttot(1:indchange(j)))/sum(deaths(1:indchange(j)));
            E(j) = exp(-tchange(j)/lamest(j));
        else
            lamest(j) = sum(ttot(indchange(j-1)+1:indchange(j)))/...
                sum(deaths(indchange(j-1)+1:indchange(j)));
            E(j) = E(j-1)*exp(-(tchange(j)-tchange(j-1))/lamest(j));
        end
    end
    lamest(nchange+1) = sum(ttot(indchange(nchange)+1:ntime))/...
               sum(deaths(indchange(nchange)+1:ntime));
    for k = 1:length(tx),
        if tx(k) < tchange(1), % the first piece
            quan(k) = exp(-tx(k)/lamest(1));
        elseif tx(k) < tchange(nchange), % piece 2 -- nchange 
            for j = 2:nchange,
                if tx(k) >= tchange(j-1),
                    if tx(k) < tchange(j),
                        quan(k) = E(j-1)*exp(-(tx(k)-tchange(j-1))/lamest(j));
                    end
                end
            end
        else % piece nchange+1
            quan(k) = E(nchange)* exp(-(tx(k)-tchange(nchange))/...
                lamest(nchange+1));
        end
    end
else % nchange < 1
    lamest = sum(ttot)/sum(deaths);
    for k = 1:length(tx),
        quan(k) = exp(-tx(k)/lamest);
    end
end
% 
% tchange
% E
% lamest