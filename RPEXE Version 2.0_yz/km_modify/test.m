addpath('..\RPEXE Version 1.0 (03032014)\MATLAB version');
addpath('..\data');

[bstnum] = xlsread('km_test_data.xlsx');
% Three columns in the number: indicator of the validation, censor
%   indicator, time to drfs; 
% Two columns in the text: er status ("P" or "N"), 
%   chemo prediction: ("Rx Sensitive" or "Rx Insensitive");
d_time = bstnum(:,1);
censor = bstnum(:,2);

% All training;
figure(1);
[xpart1 ypart1] = km_modify(d_time,censor, 1);
xlabel('Years');
ylabel('Survival probability');