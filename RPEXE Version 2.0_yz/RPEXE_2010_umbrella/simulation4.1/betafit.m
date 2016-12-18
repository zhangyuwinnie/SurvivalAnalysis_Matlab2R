function dis = betafit(par,per5,per10)
% This function calcultes the discrepancy between the 
% beta fits and the (5, 10) percentiles.
% 
% Inputs 
%       per5: 5 percent quantile
%       per10: 10 percent quantile
%       par: 2 by 1 parameter vector of the beta distribution.
% Output
%       dis: log of the sum squared discrepancy between the given 
%           props (0.05, 0.1) and the beta proportions


dis = log((betacdf(per5,par(1),par(2))-0.05)^2) + log((betacdf(per10,par(1),par(2))-0.1)^2);


