function loglik = gamllik(structtime, structttot, structdeaths,...
                    time_die, ttot, deaths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% A function computing the log likelihood from the gamma distribution under
%       an order restriction reduction 
% Inputs
%       structtime   == times under restriction
%       structttot   == ttots between each time point and the previous 
%                           time point (or 0) under restriction
%       structdeaths == number of deaths corresponding to struct.ttot
%       time_die      == all possible times to make the cuts
%       ttot          == ttots corresponding to time_die
%       deaths        == dealthes corresponding to ttot
% Outputs      
%       loglik        == log of the likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the Gamma parameter
for j = 1: length(structtime),
    structgamindi(j) = structttot(j)/structdeaths(j);
end
% the likelihood
% get the indices of the cut time:
structindi = zeros(length(structtime),1);
for j = 1:length(structtime),
    for jj = 1:length(time_die),
        if time_die(jj) == structtime(j);
            structindi(j) = jj;
        end
    end
end
% set the scale parameter for Gamma distribution of the ttot
structgampar = zeros(length(time_die),1);
for ii = 1:length(time_die),
    for j = 1:length(structtime),
        if ii <= structindi(j),
            if j ==1, 
                structgampar(ii) = structgamindi(j);
            else 
                if ii > structindi(j-1),
                    structgampar(ii) = structgamindi(j);
                end
            end
        end
    end
end
loglik = 0;
for ii = 1 : length(structgampar),
    loglik = loglik + log(gampdf(ttot(ii)/deaths(ii),deaths(ii),...
         structgampar(ii)));
end


