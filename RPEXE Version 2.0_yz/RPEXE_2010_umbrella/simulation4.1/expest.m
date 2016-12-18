function [est ci] = expest(ttot,d,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a function returns the estimate and the confidence interval
% inputs:
%   ttot = total time on test
%   d    = the number of events
%   alpha= the size of the test
% outputs
%   est  = estimated model parameter 
%   ci   = confidence interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

est = ttot/d;

ci  = [2*est*d/gaminv(1-alpha/2,d,2) 2*est*d/gaminv(alpha/2,d,2)];
