
% Include the MATLAB files' directory
cd '/Users/ganghan/Desktop/work/GH work/topics/Testing Exponential/code/RPEXE Version 2.0 (June2014)';
addpath('./MATLAB version/');



% Load the data
data_ex1 = textread('data.ex1.txt');
oscensor      = data_ex1(:,1);
ostime        = data_ex1(:,2);

% Plot the data
figure(1)
[xpart,ypart] = km(ostime, oscensor, 0);
figure(2)
[xpart,ypart] = km2(ostime, oscensor, 0);
figure(3)
[xpart,ypart] = km_log(ostime, oscensor, 0);
figure(4)
[xpart,ypart] = km_log2(ostime, oscensor, 0);


% If you run original code you get the same results as version 1;

pexeout1     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',1);
pexeout1
[pexeout1.times pexeout1.pvalues]

pexeout2     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',3);
pexeout2
[pexeout2.times pexeout2.pvalues]
pexeout2.struct
pexeout2.struct(1)
[pexeout2.struct(1).time pexeout2.struct(1).ttot pexeout2.struct(1).deaths]
pexeout2.struct(2)
[pexeout2.struct(2).time pexeout2.struct(2).ttot pexeout2.struct(2).deaths]




% Below are new codes that make plots and computes p-values at a fixed time;


% Fit RPEXE with decreasing failure rate assumption, critical value 0.05, 
% decreasing failure rate
pexeout0     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',1,...
    'Criticalp', 0.05);
pexeout0
[pexeout0.times pexeout0.pvalues]
[pexeout0.times_c pexeout0.pvalues_c]



% Fit RPEXE with decreasing failure rate assumption, critical value 0.5,
% decreasing failure rate
pexeout1     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',1,...
    'Criticalp', 0.5);
pexeout1
[pexeout1.times pexeout1.pvalues]
[pexeout1.times_c pexeout1.pvalues_c]



% Fit RPEXE with monotonic failure rate assumption, critical value 0.7,
% monotonic failure rate
pexeout2     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',3,...
    'Criticalp', 0.7);
pexeout2
[pexeout2.times pexeout2.pvalues]
[pexeout2.times_c pexeout2.pvalues_c]



% Fit RPEXE with decreasing then increasing failure rate assumption, 
% critical value 0.7 
% decreasing then increasing failure rate
pexeout3     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',5,...
    'Criticalp', 0.7);
pexeout3
[pexeout3.times pexeout3.pvalues]
[pexeout3.times_c pexeout3.pvalues_c]




% Fit RPEXE with decreasing failure rate assumption, critical value 1e-8, 
% decreasing failure rate; The result is only one exponential fit;
pexeout4     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',1,...
    'Criticalp', 1e-8);
pexeout4
[pexeout4.times pexeout4.pvalues]
%[pexeout4.times_c pexeout4.pvalues_c]

