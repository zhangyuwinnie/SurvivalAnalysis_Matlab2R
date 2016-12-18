function [lam,dea] = plotpexe_driver3piece(ostime,oscensor,t1,t2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% given t1, t2, 
% compute lam1, lam2, lam3
%       using ostime and oscensor
% Then plot the PEXE overlaid with Kaplan-Meier estimate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[time_die,ttot,deaths] = totaltest(ostime,oscensor);
% compute lam1 - 3
% t1 = 1.6438;
% t2 = 29.0301;
% index for lam1-3 
% this "+0.00005" is to deal with the round off error
ind1  = find(time_die<t1+0.00005);
ind12 = find(time_die<t2+0.00005);
ind2  = setdiff(ind12, ind1);
ind3  = find(time_die>=t2+0.00005);
lam1  = sum(ttot(ind1))/sum(deaths(ind1));
lam2  = sum(ttot(ind2))/sum(deaths(ind2));
lam3  = sum(ttot(ind3))/sum(deaths(ind3));
lam   = [lam1;lam2;lam3];
dea   = [sum(deaths(ind1)); sum(deaths(ind2)); sum(deaths(ind3))];
xaxis = [0:0.001:1]*max(ostime);
y2    = pexe3piece(lam1,lam2,lam3,t1,t2,xaxis);
[xpart,ypart] = km_overlay(ostime, oscensor, xaxis, y2);