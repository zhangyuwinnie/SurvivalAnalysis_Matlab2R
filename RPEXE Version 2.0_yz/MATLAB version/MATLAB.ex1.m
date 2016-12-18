
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


% Fit RPEXE with decreasing failure rate assumption 
pexeout1     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',1);
pexeout1
[pexeout1.times pexeout1.pvalues]
% pexeout1 = 
%       trend: 'Decreasing failure rate'
%       times: [9x1 double]
%     pvalues: [9x1 double]
% ans =
%     1.8932    0.9780
%    23.4603    0.9608
%    47.8548    0.7539
%    12.7096    0.7045
%    24.7096    0.6876
%     0.0849    0.6013
%     2.3863    0.3060
%    18.5945    0.0852
%    28.0301    0.0000    


% Fit RPEXE with monotonic failure rate assumption 
pexeout2     = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',3);
pexeout2
[pexeout2.times pexeout2.pvalues]
pexeout2.struct
pexeout2.struct(1)
[pexeout2.struct(1).time pexeout2.struct(1).ttot pexeout2.struct(1).deaths]
pexeout2.struct(2)
[pexeout2.struct(2).time pexeout2.struct(2).ttot pexeout2.struct(2).deaths]
% pexeout2 = 
%      struct: [1x2 struct]
%       trend: 'Monotone failure rate'
%       times: [9x1 double]
%     pvalues: [9x1 double]
% ans =
%     1.8932    0.9783
%    23.4603    0.9608
%    47.8548    0.7539
%    12.7096    0.7045
%    24.7096    0.6877
%     0.0849    0.6010
%     2.3863    0.3060
%    18.5945    0.0852
%    28.0301    0.0000
% ans = 
% 1x2 struct array with fields:
%     time
%     ttot
%     deaths
%     loglik
% ans = 
%       time: [10x1 double]
%       ttot: [10x1 double]
%     deaths: [10x1 double]
%     loglik: -543.9703
% ans =
%    1.0e+03 *
%     0.0001    0.0151    0.0020
%     0.0019    0.3021    0.0250
%     0.0024    0.0733    0.0060
%     0.0127    1.0907    0.0760
%     0.0186    0.3450    0.0220
%     0.0235    0.1992    0.0100
%     0.0247    0.0412    0.0020
%     0.0280    0.0997    0.0040
%     0.0479    0.3699    0.0060
%     0.0511    0.0825    0.0010
% ans = 
%       time: [2x1 double]
%       ttot: [2x1 double]
%     deaths: [2x1 double]
%     loglik: -559.4657
% ans =
%    1.0e+03 *
%     0.0006    0.1130    0.0050
%     0.0511    2.5057    0.1490

