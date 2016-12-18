function lowerbd = Exp_Med_lower_bound(parest,n,alp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program returns the lower bound of the parameter of exponential
% distribution
% Inputs:
%      parest == estimated parameter from exponential distribution
%      n      == sample size
%      alp    == alpha level of the test
% Output: 
%      lowerbd== lower bound of the one-sided confidence interval with 
%                level alp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nume = 2*n*parest;
dem  = gaminv(1-alp,n,2);
lowerbd   = nume/dem;




