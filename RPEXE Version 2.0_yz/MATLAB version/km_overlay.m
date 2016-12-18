function [xpart,ypart] = km_overlay(time, censor, x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the Kaplan-Meier curve overlaied 
% with a parametric estimate.                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xpart,ypart] = km(time, censor);

hold on;
plot(x,y,'--r','LineWidth',2);
hold off;






