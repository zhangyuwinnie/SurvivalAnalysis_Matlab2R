function t = ggaminv(p,beta,lambda,sigma)
% The inverse CDF function of the generalized gamma distribution. 
% Use the paramters (beta,lambda,sigma) to compute the inverse cdf at p.
% Input 
%   p: the probability at which inverse cdf shall be computed 
%   beta, lambda, sigma: parameters of the genralized gamma distribution
% Output
%   t: inverse cdf of the generalized gamma at p

t = (lambda^(-2)*gaminv(p,lambda^(-2),1))^(sigma/lambda)*exp(beta);



