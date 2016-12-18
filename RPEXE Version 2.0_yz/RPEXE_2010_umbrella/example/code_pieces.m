
% [input1 txt1] = xlsread('gang.xls');
% oscensor      = input1(:,2);
% ostime        = input1(:,3);


% Overall Survival, without the first month
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella\example');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\commonly used functions');

% [input1 txt1] = xlsread('os_short.xls');
% oscensor      = input1(:,1);
% ostime        = input1(:,2);

[input1 txt1] = xlsread('gang.xls');
oscensor      = input1(:,2);
ostime        = input1(:,3);
labless1      = find(ostime<1);
ostime(labless1)   = [];
oscensor(labless1) = [];
%ostime   = ostime -1;


pexeout0      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',0);
[pexeout0.times pexeout0.pvalues]


pexeout1      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',1);
[pexeout1.times pexeout1.pvalues]
pexeout1
% decreasing failure rate
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
% pexeout1 = 
%       trend: 'Decreasing failure rate'
%       times: [9x1 double]
%     pvalues: [9x1 double]
    

[time_die,ttot,deaths] = totaltest(ostime-1,oscensor);
plotspacinghazard(time_die,ttot,deaths)

pexeout2      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',2);
[pexeout2.times pexeout2.pvalues]
pexeout2
% increasing failure rate
% ans =
%     0.6438    0.4659

t1 = 28.0301;
[lam,dea] = plotpexe_driver2piece(ostime-1,oscensor,t1);
% lam =
%    14.7370
%    64.6215
% dea =
%    147
%      7


% pexeout3      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',3);
% [pexeout3.times pexeout3.pvalues]
% pexeout3
% % monotone failure rate
% % pexeout3      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',3);
% % [pexeout3.times pexeout3.pvalues]
% % ans =
% %     1.8932    0.9780
% %    23.4603    0.9608
% %    47.8548    0.7539
% %    12.7096    0.7045
% %    24.7096    0.6876
% %     0.0849    0.6013
% %     2.3863    0.3060
% %    18.5945    0.0852
% %    28.0301    0.0000
% 
% 
% pexeout4      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',4);
% [pexeout4.times pexeout4.pvalues]
% pexeout4
% % increasing then decreasing failure rate
% % ans =
% %    23.4603    0.9608
% %    47.8548    0.7539
% %     5.2795    0.7208
% %    24.7096    0.6876
% %    12.7096    0.6500
% %     5.3452    0.4152
% %     0.6438    0.1855
% %    18.5945    0.1441
% %     5.3781    0.0044
% %     5.4110    0.3086
% %    28.0301    0.0000
% % pexeout4 = 
% %      struct: [1x139 struct]
% %     changet: 5.3781
% %       trend: 'Increasing-decreasing failure rate'
% %       times: [11x1 double]
% %     pvalues: [11x1 double]
% 
% pexeout5      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',5);
% [pexeout5.times pexeout5.pvalues]
% pexeout5
% % decreasing then increasing failure rate
% % ans =
% %     1.8932    0.9780
% %    23.4603    0.9608
% %    12.7096    0.7045
% %    24.7096    0.6876
% %    35.8548    0.6152
% %     0.0849    0.6013
% %    42.2329    0.5430
% %     2.3863    0.3060
% %    18.5945    0.0852
% %    28.0301    0.0000
% % pexeout5 = 
% %      struct: [1x139 struct]
% %     changet: 35.8548
% %       trend: 'Decreasing-increasing failure rate'
% %       times: [10x1 double]
% %     pvalues: [10x1 double]
% 
% 
% pexeout6      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',6);
% [pexeout6.times pexeout6.pvalues]
% pexeout6
% % increasing then decreasing failure rate, without including the peak point
% % ans =
% %    23.4603    0.9608
% %    47.8548    0.7539
% %     5.2795    0.7208
% %    24.7096    0.6876
% %    12.7096    0.6500
% %     0.6438    0.1977
% %    18.5945    0.1441
% %     5.3452    0.0053
% %     5.4110    0.3086
% %    28.0301    0.0000
% % pexeout6 = 
% %      struct: [1x139 struct]
% %     changet: 5.3781
% %       trend: 'Increasing-decreasing failure rate'
% %       times: [10x1 double]
% %     pvalues: [10x1 double]
% 
% 
% pexeout7      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Trend',7);
% [pexeout7.times pexeout7.pvalues]
% pexeout7
% % decreasing then increasing failure rate, without including the peak point
% % ans =
% %     1.8932    0.9780
% %    23.4603    0.9608
% %    12.7096    0.7045
% %    24.7096    0.6876
% %     0.0849    0.6013
% %    42.2329    0.5430
% %     2.3863    0.3060
% %    18.5945    0.0852
% %    28.0301    0.0000
% % pexeout7 = 
% %      struct: [1x139 struct]
% %     changet: 31.3178
% %       trend: 'Decreasing-increasing failure rate'
% %       times: [9x1 double]
% %     pvalues: [9x1 double]







% Overall Survival, With the first Month
[input1 txt1] = xlsread('gang.xls');
oscensor      = input1(:,2);
ostime        = input1(:,3);
% labless1      = find(ostime<1);
% ostime(labless1)   = [];
% oscensor(labless1) = [];
% ostime   = ostime -1;


% Kaplan-Meier plot;
[xpart ypart] = km(ostime,oscensor);

% Plot the instantaneous hazard
[time_die,ttot,deaths] = totaltest(ostime,oscensor);
plotspacinghazard(time_die,ttot,deaths)



pexeout0      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',0);
[pexeout0.times pexeout0.pvalues]


pexeout1      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',1);
[pexeout1.times pexeout1.pvalues]
pexeout1
pexeout1
% ans =
%    24.4603    0.9608
%    13.7096    0.8018
%    48.8548    0.7539
%    25.7096    0.6876
%    19.5945    0.1543
%    29.0301    0.0000
% pexeout1 = 
%       trend: 'Decreasing failure rate'
%       times: [6x1 double]
%     pvalues: [6x1 double]
    

[time_die,ttot,deaths] = totaltest(ostime,oscensor);
pexeout2      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',2);
[pexeout2.times pexeout2.pvalues]
pexeout2
% ans =
%     1.3808    0.8258
%     1.0521    0.4156
%     1.6438    0.0133
% pexeout2 = 
%       trend: 'Increasing failure rate'
%       times: [3x1 double]
%     pvalues: [3x1 double]


pexeout3      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',3);
[pexeout3.times pexeout3.pvalues]
pexeout3
% ans =
%    24.4603    0.9608
%    13.7096    0.8018
%    48.8548    0.7539
%    25.7096    0.6876
%    19.5945    0.1543
%    29.0301    0.0000
% pexeout3 = 
%      struct: [1x2 struct]
%       trend: 'Monotone failure rate'
%       times: [6x1 double]
%     pvalues: [6x1 double]


pexeout4      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',4);
[pexeout4.times pexeout4.pvalues]
pexeout4
% ans =
%    24.4603    0.9608
%     1.3808    0.8258
%    48.8548    0.7539
%     6.2795    0.7208
%    25.7096    0.6876
%    13.7096    0.6500
%     1.0521    0.4156
%     6.3452    0.4152
%    19.5945    0.1441
%     6.3781    0.0053
%     6.4110    0.1723
%     1.6438    0.0026
%    29.0301    0.0000
% pexeout4 = 
%      struct: [1x142 struct]
%     changet: 6.3781
%       trend: 'Increasing-decreasing failure rate'
%       times: [13x1 double]
%     pvalues: [13x1 double]

% given t1, t2, 
% compute lam1, lam2, lam3
%       using ostime and oscensor
% [time_die,ttot,deaths] = totaltest(ostime,oscensor);
% % compute lam1 - 3
% t1 = 1.6438;
% t2 = 29.0301;
% % index for lam1-3 
% ind1  = find(time_die<t1+0.00005);
% ind12 = find(time_die<t2+0.00005);
% ind2  = setdiff(ind12, ind1);
% ind3  = find(time_die>=t2+0.00005);
% lam1  = sum(ttot(ind1))/sum(deaths(ind1));
% lam2  = sum(ttot(ind2))/sum(deaths(ind2));
% lam3  = sum(ttot(ind3))/sum(deaths(ind3));
% xaxis = [0:0.01:1]*length(ostime);
% y2    = pexe3piece(lam1,lam2,lam3,t1,t2,xaxis);
% [xpart,ypart] = km_overlay(ostime, oscensor, xaxis, y2);
% 

t1 = 1.6438;
t2 = 29.0301;
[lam,dea] = plotpexe_driver3piece(ostime,oscensor,t1,t2);
% lam =
%    36.5342
%    14.4604
%    64.6215
% dea =
%      8
%    142
%      7

figure(1);
plotspacinghazard(time_die,ttot,deaths)

figure(2);
plotspacinghazard_overlay(time_die,ttot,deaths,[t1 t2],lam)

figure(3);
plotspacinghazard_bynumber(time_die,ttot,deaths)

figure(4);
plotspacinghazard_bynumber_overlay(time_die,ttot,deaths,[dea(1) dea(1)+dea(2)],lam)


pexeout5      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',5);
[pexeout5.times pexeout5.pvalues]
pexeout5
% ans =
%    24.4603    0.9608
%    13.7096    0.8018
%    25.7096    0.6876
%    36.8548    0.6152
%    43.2329    0.5430
%    19.5945    0.1543
%    29.0301    0.0000
% pexeout5 = 
%      struct: [1x142 struct]
%     changet: 36.8548
%       trend: 'Decreasing-increasing failure rate'
%       times: [7x1 double]
%     pvalues: [7x1 double]

pexeout6      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',6);
[pexeout6.times pexeout6.pvalues]
pexeout6
% ans =
%    24.4603    0.9608
%     1.3808    0.8258
%    48.8548    0.7539
%     6.2795    0.7208
%    25.7096    0.6876
%    13.7096    0.6500
%     1.0521    0.4156
%    19.5945    0.1441
%     6.3452    0.0067
%     6.4110    0.1723
%     1.6438    0.0026
%    29.0301    0.0000
% pexeout6 = 
%      struct: [1x142 struct]
%     changet: 6.3781
%       trend: 'Increasing-decreasing failure rate'
%       times: [12x1 double]
%     pvalues: [12x1 double]


pexeout7      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Trend',7);
[pexeout7.times pexeout7.pvalues]
pexeout7
% ans =
%    24.4603    0.9608
%    13.7096    0.8018
%    25.7096    0.6876
%    43.2329    0.5430
%    19.5945    0.1543
%    29.0301    0.0000
% pexeout7 = 
%      struct: [1x142 struct]
%     changet: 36.8548
%       trend: 'Decreasing-increasing failure rate'
%       times: [6x1 double]
%     pvalues: [6x1 double]















% PFS, without the first month
[input1 txt1] = xlsread('gang.xls');
pfscensor      = input1(:,4);
pfstime        = input1(:,5);
labless1      = find(pfstime<1);
pfstime(labless1)   = [];
pfscensor(labless1) = [];
%pfstime   = pfstime -1;


% [input2 txt2] = xlsread('pfs.xls');
% pfscensor      = input2(:,1);
% pfstime        = input2(:,2);

pfsout0      = RPEXEv1('EventTime',pfstime-1,'Censor',pfscensor,'Trend',0);
[pfsout0.times pfsout0.pvalues]


pfsout1      = RPEXEv1('EventTime',pfstime-1,'Censor',pfscensor,'Trend',1);
[pfsout1.times pfsout1.pvalues]
pfsout1
% ans =
%    12.5452    0.8331
%     0.1507    0.7269
%    37.2027    0.6974
%     0.7753    0.4658
%    17.8712    0.3488
%    24.4795    0.0236
%     0.3808    0.0124
%    11.9863    0.0000
% pfsout1 = 
%       trend: 'Decreasing failure rate'
%       times: [8x1 double]
%     pvalues: [8x1 double]   


t1 = 11.9863;
[lam,dea] = plotpexe_driver2piece(pfstime-1,pfscensor,t1);
% lam =
%     5.7166
%    28.5597
% dea =
%    153
%     13


[time_die,ttot,deaths] = totaltest(pfstime-1,pfscensor);
plotspacinghazard(time_die,ttot,deaths)
% 
% pfsout2      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',2);
% [pfsout2.times pfsout2.pvalues]
% pfsout2
% % ans =
% %     0.0521    0.8188
% % pfsout2 = 
% %       trend: 'Increasing failure rate'
% %       times: 0.0521
% %     pvalues: 0.8188
% 
% pfsout3      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',3);
% [pfsout3.times pfsout3.pvalues]
% pfsout3
% pfsout3
% % ans =
% %    12.5452    0.8331
% %     0.1507    0.7269
% %    37.2027    0.6974
% %     0.7753    0.4658
% %    17.8712    0.3488
% %    24.4795    0.0236
% %     0.3808    0.0124
% %    11.9863    0.0000
% % pfsout3 = 
% %      struct: [1x2 struct]
% %       trend: 'Monotone failure rate'
% %       times: [8x1 double]
% %     pvalues: [8x1 double]
% 
% 
% 
% pfsout4      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',4);
% [pfsout4.times pfsout4.pvalues]
% pfsout4
% % ans =
% %    12.5452    0.8331
% %     0.0521    0.7714
% %    37.2027    0.6974
% %     0.7753    0.4658
% %    17.8712    0.3488
% %     0.0849    0.2234
% %     0.3479    0.0276
% %    24.4795    0.0236
% %     0.3808    0.0124
% %    11.9863    0.0000
% % pfsout4 = 
% %      struct: [1x126 struct]
% %     changet: 0.3479
% %       trend: 'Increasing-decreasing failure rate'
% %       times: [10x1 double]
% %     pvalues: [10x1 double]
% 
% pfsout5      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',5);
% [pfsout5.times pfsout5.pvalues]
% pfsout5
% % ans =
% %    12.5452    0.8331
% %     0.1507    0.7269
% %    34.4082    0.6784
% %     0.7753    0.4658
% %    17.8712    0.3488
% %    24.4795    0.0300
% %     0.3808    0.0124
% %    11.9863    0.0000
% % pfsout5 = 
% %      struct: [1x126 struct]
% %     changet: 24.4795
% %       trend: 'Decreasing-increasing failure rate'
% %       times: [8x1 double]
% %     pvalues: [8x1 double]
% 
% 
% pfsout6      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',6);
% [pfsout6.times pfsout6.pvalues]
% pfsout6
% % ans =
% %    12.5452    0.8331
% %     0.0521    0.7714
% %    37.2027    0.6974
% %     0.7753    0.4658
% %    17.8712    0.3488
% %     0.0849    0.1038
% %    24.4795    0.0236
% %     0.3808    0.0124
% %    11.9863    0.0000
% % pfsout6 = 
% %      struct: [1x126 struct]
% %     changet: 0.3479
% %       trend: 'Increasing-decreasing failure rate'
% %       times: [9x1 double]
% %     pvalues: [9x1 double]
% 
% 
% 
% pfsout7      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',7);
% [pfsout7.times pfsout7.pvalues]
% pfsout7
% % ans =
% %    12.5452    0.8331
% %     0.1507    0.7269
% %     0.7753    0.4658
% %    17.8712    0.0361
% %     0.3808    0.0124
% %    11.9863    0.0000
% % pfsout7 = 
% %      struct: [1x126 struct]
% %     changet: 24.4795
% %       trend: 'Decreasing-increasing failure rate'
% %       times: [6x1 double]
% %     pvalues: [6x1 double]
% % 
% 









% PFS, with the first month
[input1 txt1]  = xlsread('gang.xls');
pfscensor      = input1(:,4);
pfstime        = input1(:,5);
% labless1      = find(pfstime<1);
% pfstime(labless1)   = [];
% pfscensor(labless1) = [];
% pfstime   = pfstime -1;


% [input2 txt2] = xlsread('pfs.xls');
% pfscensor      = input2(:,1);
% pfstime        = input2(:,2);


% Kaplan-Meier plot;
[xpart ypart] = km(pfstime,pfscensor);

% Plot the instantaneous hazard
[time_die,ttot,deaths] = totaltest(pfstime,pfscensor);
plotspacinghazard(time_die,ttot,deaths);



pfsout0      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',0);
[pfsout0.times pfsout0.pvalues]


pfsout1      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',1);
[pfsout1.times pfsout1.pvalues]
pfsout1
% ans =
%    13.5452    0.8331
%    38.2027    0.6974
%    18.8712    0.3488
%    25.4795    0.0236
%    12.9863    0.0000
% pfsout1 = 
%       trend: 'Decreasing failure rate'
%       times: [5x1 double]
%     pvalues: [5x1 double]    

[time_die,ttot,deaths] = totaltest(pfstime,pfscensor);

pfsout2      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',2);
[pfsout2.times pfsout2.pvalues]
pfsout2
% ans =
%     0.6575    0.9882
%     0.8877    0.7069
%     0.1973    0.4686
%     1.0521    0.0004
% pfsout2 = 
%       trend: 'Increasing failure rate'
%       times: [4x1 double]
%     pvalues: [4x1 double]


pfsout3      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',3);
[pfsout3.times pfsout3.pvalues]
pfsout3
% ans =
%    13.5452    0.8331
%    38.2027    0.6974
%    18.8712    0.3488
%    25.4795    0.0236
%    12.9863    0.0000
% pfsout3 = 
%      struct: [1x2 struct]
%       trend: 'Monotone failure rate'
%       times: [5x1 double]
%     pvalues: [5x1 double]


pfsout4      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',4);
[pfsout4.times pfsout4.pvalues]
pfsout4
% ans =
%     0.6575    0.9882
%    13.5452    0.8331
%     0.8877    0.7069
%    38.2027    0.6974
%     1.0849    0.5259
%     0.1973    0.4686
%     1.7753    0.4658
%    18.8712    0.3488
%     1.3479    0.0410
%    25.4795    0.0236
%     1.3808    0.0053
%     1.0521    0.0000
%    12.9863    0.0000
% pfsout4 = 
%      struct: [1x134 struct]
%     changet: 1.3479
%       trend: 'Increasing-decreasing failure rate'
%       times: [13x1 double]
%     pvalues: [13x1 double]


t1 =  1.0521;
t2 = 12.9863;
[lam,dea] = plotpexe_driver3piece(pfstime,pfscensor,t1,t2);


% Plot the instantaneous hazard
[time_die,ttot,deaths] = totaltest(pfstime,pfscensor);
% plotspacinghazard_overlay(time_die,ttot,deaths,[t1 t2],lam);
% 
% plotspacinghazard_bynumber(time_die,ttot,deaths)


figure(1);
plotspacinghazard(time_die,ttot,deaths)

figure(2);
plotspacinghazard_overlay(time_die,ttot,deaths,[t1 t2],lam)

figure(3);
plotspacinghazard_bynumber(time_die,ttot,deaths)

figure(4);
plotspacinghazard_bynumber_overlay(time_die,ttot,deaths,[dea(1) dea(1)+dea(2)],lam)



% lam =
%    20.7196
%     5.6950
%    28.5597
% dea =
%      9
%    152
%     13



pfsout5      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',5);
[pfsout5.times pfsout5.pvalues]
pfsout5
% ans =
%    13.5452    0.8331
%    35.4082    0.6784
%    18.8712    0.3488
%    25.4795    0.0300
%    12.9863    0.0000
% pfsout5 = 
%      struct: [1x134 struct]
%     changet: 25.4795
%       trend: 'Decreasing-increasing failure rate'
%       times: [5x1 double]
%     pvalues: [5x1 double]


pfsout6      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',6);
[pfsout6.times pfsout6.pvalues]
pfsout6
% ans =
%     0.6575    0.9882
%    13.5452    0.8331
%     0.8877    0.7069
%    38.2027    0.6974
%     0.1973    0.4686
%     1.7753    0.4658
%     1.0521    0.3849
%    18.8712    0.3488
%    25.4795    0.0236
%     1.3808    0.0037
%     1.0849    0.0000
%    12.9863    0.0000
% pfsout6 = 
%      struct: [1x134 struct]
%     changet: 1.3479
%       trend: 'Increasing-decreasing failure rate'
%       times: [12x1 double]
%     pvalues: [12x1 double]



pfsout7      = RPEXEv1('EventTime',pfstime,'Censor',pfscensor,'Trend',7);
[pfsout7.times pfsout7.pvalues]
pfsout7
% ans =
%    13.5452    0.8331
%    18.8712    0.0361
%    12.9863    0.0000
% pfsout7 = 
%      struct: [1x134 struct]
%     changet: 25.4795
%       trend: 'Decreasing-increasing failure rate'
%       times: [3x1 double]
%     pvalues: [3x1 double]



















% %%% four regresion lines:
% X  = [1 20;1 50;1 200;1 800];
% y1 = [log(0.0016) log(0.00044) log(0.00008) log(0.00001)]';
% y2 = [log(0.0048) log(0.0028) log(0.0017) log(0.0010)]';
% y3 = [log(0.0110) log(0.0068) log(0.0038) log(0.0026)]';
% y4 = [log(0.0049) log(0.0029) log(0.0018) log(0.0011)]';
% 
% beta1 = X\y1
% beta2 = X\y2
% beta3 = X\y3
% beta4 = X\y4
% 
% 
% [time_die,ttot,deaths] = totaltest(ostime-1,oscensor);
% [time_die,ttot,deaths]
% ind = zeros(length(ttot),1);
% for i = 1:length(ttot),
%     if time_die(i) < 28.0302,
%         ind(i) = 1;
%     end
% end
% lam1 = sum(ttot(find(ind)))/sum(deaths(find(ind)));
% lam2 = sum(ttot(find(1-ind)))/sum(deaths(find(1-ind)));
%     
% ind1 = zeros(length(ttot),1);
% ind2 = zeros(length(ttot),1);
% t1   = 18.595;
% t2   = 28.031; 
% for i = 1:length(ttot),
%     if time_die(i) < t1,
%         ind1(i) = 1;
%     elseif time_die(i) < t2,
%         ind2(i) = 1;
%     end
% end
% lam1 = sum(ttot(find(ind1)))/sum(deaths(find(ind1)))
% lam2 = sum(ttot(find(ind2)))/sum(deaths(find(ind2)))
%     
%     
% 
% [input1 txt1] = xlsread('gang.xls');
% oscensor      = input1(:,2);
% ostime        = input1(:,3);
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_Oct09');
% pexeout1      = RPEXEv1('EventTime',ostime-1,'Censor',oscensor,'Monotone',1);
% [pexeout1.times pexeout1.pvalues]
% 
% 
% [input1 txt1] = xlsread('gang.xls');
% oscensor      = input1(:,2);
% ostime        = input1(:,3);
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_Oct09');
% pexeout1      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Monotone',3);
% [pexeout1.times pexeout1.pvalues]
% 
% 
% [input1 txt1] = xlsread('gang.xls');
% oscensor      = input1(:,2);
% ostime        = input1(:,3);
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_Oct09');
% pexeout1      = RPEXEv1('EventTime',ostime,'Censor',oscensor,'Monotone',3);
% [pexeout1.times pexeout1.pvalues]


