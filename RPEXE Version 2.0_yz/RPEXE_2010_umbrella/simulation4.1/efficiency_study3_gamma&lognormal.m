Ns    = 100;
comb1 = [0.5 240];
comb2 = [0.8 150];
comb3 = [0.95 126.3];
comb4 = [1.2 100];
comb5 = [1.5 80];
comb6 = [3 40];
comb7 = [10 12];

res_gam1 = gamtest(comb1(1),comb1(2),Ns);
res_gam2 = gamtest(comb2(1),comb2(2),Ns);
res_gam3 = gamtest(comb3(1),comb3(2),Ns);
res_gam4 = gamtest(comb4(1),comb4(2),Ns);
res_gam5 = gamtest(comb5(1),comb5(2),Ns);
res_gam6 = gamtest(comb6(1),comb6(2),Ns);
res_gam7 = gamtest(comb7(1),comb7(2),Ns);

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(0.5,240)'},'A1:A1')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam1,'A2:C16')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam1(:,3)),'D2:D16')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(0.8,150)'},'A17:A17')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam2,'A18:C32')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam2(:,3)),'D18:D32')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(0.95,126.3)'},'A33:A33')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam3,'A34:C48')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam3(:,3)),'D34:D48')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(1.2,100)'},'A49:A49')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam4,'A50:C64')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam4(:,3)),'D50:D64')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(1.5,80)'},'A65:A65')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam5,'A66:C80')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam5(:,3)),'D66:D80')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(3,40)'},'A81:A81')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam6,'A82:C96')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam6(:,3)),'D82:D96')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',{'gamma par=(10,12)'},'A97:A97')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',res_gam7,'A98:C112')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gamest1.xls',rankavec(res_gam7(:,3)),'D98:D112')





% 
% 
% %%% Compare the Monte Carlo variances of four predictors in 
% %%% estimating the quantile in the day 100.
% 
% 
% % 1, Kaplan-Meier binomial distribution
% % 2, Kaplan-Meier averaging over length 30
% % 3, Exponential all the way
% % 4, Exponential up to a point 150
% 
% % There are more adds on in this version of the code.
% % 5, Kaplan-Meier averaging over length 10
% % 6, gamma
% % 7, lognormal
% % 8, weibull
% % 9, Exponential up to a point 60
% % 10, Exponential up to a point 90
% % 11, Exponential up to a point 120
% % 12, Exponential up to a point 180
% % 13, Exponential up to a point 210
% % 14, RPEXE monotone restriction with alpha^star = 0.005
% 
% 
% % 15, Gompertz
% % 16, generalized gamma
% 
% % We fit the functions by looking maximizing the likelihood functions.
% % There is no censored observations in any of the cases
% 
% % set the number of iteration samples
% N = 5000;
% 
% kme1 = zeros(N,1);
% kme2 = zeros(N,1);
% exp1 = zeros(N,1);
% exp2 = zeros(N,1);
% kme3 = zeros(N,1);
% gam1 = zeros(N,1);
% wei1 = zeros(N,1);
% ln1  = zeros(N,1);
% exp9 = zeros(N,1);
% exp10 = zeros(N,1);
% exp11 = zeros(N,1);
% exp12 = zeros(N,1);
% exp13 = zeros(N,1);
% exp14 = zeros(N,1);
% exp15 = zeros(N,1);
% 
% rpexrec = zeros(N,1);
% 
% Ns   = 30;
% crival = exp(-4.2331-0.3938*log(Ns));
% 
% for iter = 1:N,
%     
% % Step 1, generate data from exponential distribution 120, with sample size 50
% exprd = sort(gamrnd(0.5,240,Ns,1));
% 
% % Step 2, estimate the quantile at 100 days using four approaches
% % 1, KME
% qan1   = length(find(exprd>100))/length(exprd);
% 
% % 2, KME with averaging 85-115 days window
% tmpq   = intersect(find(exprd>85),find(exprd<115)); 
% np     = length(tmpq);
% if np <2, 
%     qan2   = qan1;
% else
%     qan2   = 0;
%         % calculate the area
%     for i = 1:(np-1),
%         qan2 = qan2 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
%             length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan2 = qan2/(exprd(tmpq(np))-exprd(tmpq(1)));
% end
% 
% % 3, exponential using all the data
% expar1 = mean(exprd);
% stdpar1= std(exprd);
% qan3 = exp(-100/expar1);
% 
% 
% % 4, exponential using data regarding whatever larger than 150 as censored
% expar2 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
%     /length(find(exprd<=150));
% qan4 = exp(-100/expar2);
% 
% 
% % 5, KME with 10 days window from days 95 to 105.
% tmpq2   = intersect(find(exprd>95),find(exprd<105)); 
% np2      = length(tmpq2);
% if np2 < 2, 
%     qan5   = qan1;
% else
%     qan5   = 0;
%         % calculate the area
%     for i = 1:(np2-1),
%         qan5 = qan5 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
%             length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan5 = qan5/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
% end
% 
% 
% % 6, gamma distribution
% % find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% % and b0 = ave(x)/a0;
% lgavex = log(mean(exprd));
% avelgx = mean(log(exprd));
% a0     = 0.5/(lgavex-avelgx);
% %options = optimset('Display','off','Simplex','on');
% options = optimset('Display','off','LargeScale','on',...
%                    'TolFun',.0001, 'Algorithm','active-set');
% [a fval]   = fmincon(@(a) gamlike([a; mean(exprd)/a],exprd),a0,...
%     [],[],[],[],0.00001,10000,[],options); 
% qan6 = 1-gamcdf(100,a,mean(exprd)/a);
% 
% 
% % 7, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
% qan7 = 1-wblcdf(100,wei(1),wei(2));
% 
% 
% % 8, lognormal distribution
% % m     = expar1;
% % v     = stdpar1;
% % lmu   = log(m^2/sqrt(v+m^2));
% % lsigma= sqrt(log(v/m^2 + 1)); 
% % [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
% %     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
% logparhat = lognfit(exprd);
% qan8 = 1-logncdf(100,logparhat(1),logparhat(2));
% % lnorm
% % qan8
% 
% % try the exponential distribution with multiple ratio of time and 
% %  exponential parameter
% % 150/120 = 1.25
% % try ratio = 0.5, 0.75,  1, 1.5, 1.75, 2
% %    time   = 60,   90, 120, 180, 210,  240
% 
% % 9, exponential using data regarding whatever larger than 60 as censored
% expar9 = (sum(exprd(find(exprd<=60)))+length(find(exprd>60))*60)...
%     /length(find(exprd<=60));
% qan9 = exp(-100/expar9);
% 
% 
% % 10, exponential using data regarding whatever larger than 90 as censored
% expar10 = (sum(exprd(find(exprd<=90)))+length(find(exprd>90))*90)...
%     /length(find(exprd<=90));
% qan10 = exp(-100/expar10);
% 
% 
% % 11, exponential using data regarding whatever larger than 120 as censored
% expar11 = (sum(exprd(find(exprd<=120)))+length(find(exprd>120))*120)...
%     /length(find(exprd<=120));
% qan11 = exp(-100/expar11);
% 
% 
% % 12, exponential using data regarding whatever larger than 180 as censored
% expar12 = (sum(exprd(find(exprd<=180)))+length(find(exprd>180))*180)...
%     /length(find(exprd<=180));
% qan12 = exp(-100/expar12);
% 
% 
% % 13, exponential using data regarding whatever larger than 210 as censored
% expar13 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
%     /length(find(exprd<=210));
% qan13 = exp(-100/expar13);
% 
% 
% % 14, pexe estimate
% [time_die14,ttot14,deaths14] = totaltest(exprd,ones(Ns,1));
% lendis = length(ttot14);
% [qan14 lamest14] = pexeest(exprd, ones(Ns,1), time_die14(1:lendis-1), 100);
% 
% 
% % 15, rpexe estimate
% % Estimate the first change point
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
% pexeout     = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
% if min(pexeout.pvalues) < crival, % if the change point is found
%     tchange = pexeout.times(find(pexeout.pvalues < crival));
%     [qan15 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
%     rpexrec(iter) = min(pexeout.pvalues);
% else 
%     expar15 = expar1;
%     qan15   = qan3;
% end
% 
% 
% kme1(iter) = qan1;
% kme2(iter) = qan2;
% exp1(iter) = qan3;
% exp2(iter) = qan4;
% kme3(iter) = qan5;
% gam1(iter) = qan6;
% wei1(iter) = qan7;
% ln1(iter)  = qan8;
% exp9(iter) = qan9;
% exp10(iter) = qan10;
% exp11(iter) = qan11;
% exp12(iter) = qan12;
% exp13(iter) = qan13;
% exp14(iter) = qan14;
% exp15(iter) = qan15;
% end
% 
% 
% % 
% quanval = 1-gamcdf(100,0.5,240);
% 
% % compute the roon mean squared prediction error
% rmspe1 = sqrt(sum((kme1-quanval)'*(kme1-quanval)/N))
% rmspe2 = sqrt(sum((kme2-quanval)'*(kme2-quanval)/N))
% rmspe3 = sqrt(sum((exp1-quanval)'*(exp1-quanval)/N))
% rmspe4 = sqrt(sum((exp2-quanval)'*(exp2-quanval)/N))
% rmspe5 = sqrt(sum((kme3-quanval)'*(kme3-quanval)/N))
% rmspe6 = sqrt(sum((gam1-quanval)'*(gam1-quanval)/N))
% rmspe7 = sqrt(sum((wei1-quanval)'*(wei1-quanval)/N))
% rmspe8 = sqrt(sum((ln1-quanval)'*(ln1-quanval)/N))
% rmspe9 = sqrt(sum((exp9-quanval)'*(exp9-quanval)/N))
% rmspe10 = sqrt(sum((exp10-quanval)'*(exp10-quanval)/N))
% rmspe11 = sqrt(sum((exp11-quanval)'*(exp11-quanval)/N))
% rmspe12 = sqrt(sum((exp12-quanval)'*(exp12-quanval)/N))
% rmspe13 = sqrt(sum((exp13-quanval)'*(exp13-quanval)/N))
% rmspe14 = sqrt(sum((exp14-quanval)'*(exp14-quanval)/N))
% rmspe15 = sqrt(sum((exp15-quanval)'*(exp15-quanval)/N))
% 
% 
% res30_gam1 = [[var(kme1) mean(kme1)-quanval rmspe1]
% [var(kme2) mean(kme2)-quanval rmspe2]
% [var(kme3) mean(kme3)-quanval rmspe5]
% [var(gam1) mean(gam1)-quanval rmspe6]
% [var(wei1) mean(wei1)-quanval rmspe7]
% [var(ln1) mean(ln1)-quanval rmspe8]
% [var(exp1) mean(exp1)-quanval rmspe3]
% [var(exp9) mean(exp9)-quanval rmspe9]
% [var(exp10) mean(exp10)-quanval rmspe10]
% [var(exp11) mean(exp11)-quanval rmspe11]
% [var(exp2) mean(exp2)-quanval rmspe4]
% [var(exp12) mean(exp12)-quanval rmspe12]
% [var(exp13) mean(exp13)-quanval rmspe13]
% [var(exp14) mean(exp14)-quanval rmspe14]
% [var(exp15) mean(exp15)-quanval rmspe15]]
% 
% 
% 
% 
% 
% % 
% % 
% % subplot(2,2,1);
% % hist(kme1);
% % title('Kaplan-Meier Estimate');
% % xlim([0.1,0.9]);
% % subplot(2,2,2);
% % hist(kme2);
% % title('Kaplan-Meier Estimate by averaging over window of 30');
% % xlim([0.1,0.9]);
% % subplot(2,2,3);
% % hist(exp1);
% % title('Exponential Estimate');
% % xlim([0.1,0.9]);
% % subplot(2,2,4);
% % hist(exp2);
% % title('Exponential Estimate with censoring at 150');
% % xlim([0.1,0.9]);
% % 
% % % compute the roon mean squared prediction error
% % rmspe1 = sqrt(sum((kme1-0.4346)'*(kme1-0.4346)/5000))
% % rmspe2 = sqrt(sum((kme2-0.4346)'*(kme2-0.4346)/5000))
% % rmspe3 = sqrt(sum((exp1-0.4346)'*(exp1-0.4346)/5000))
% % rmspe4 = sqrt(sum((exp2-0.4346)'*(exp2-0.4346)/5000))
% % 
% % 
% % [var(kme1) mean(kme1)-0.4346]
% % [var(kme2) mean(kme2)-0.4346]
% % [var(exp1) mean(exp1)-0.4346]
% % [var(exp2) mean(exp2)-0.4346]
% 
% 
% 
% 
% 
% 
% % We fit the functions by looking maximizing the likelihood functions.
% % There is no censored observations in any of the cases
% 
% % set the number of iteration samples
% N = 5000;
% 
% kme1 = zeros(N,1);
% kme2 = zeros(N,1);
% exp1 = zeros(N,1);
% exp2 = zeros(N,1);
% kme3 = zeros(N,1);
% gam1 = zeros(N,1);
% wei1 = zeros(N,1);
% ln1  = zeros(N,1);
% exp9 = zeros(N,1);
% exp10 = zeros(N,1);
% exp11 = zeros(N,1);
% exp12 = zeros(N,1);
% exp13 = zeros(N,1);
% exp14 = zeros(N,1);
% exp15 = zeros(N,1);
% 
% rpexrec = zeros(N,1);
% 
% Ns   = 30;
% crival = exp(-4.2331-0.3938*log(Ns));
% 
% for iter = 1:N,
%     
% % Step 1, generate data from exponential distribution 120, with sample size 50
% exprd = sort(gamrnd(0.9,133.3,Ns,1));
% 
% % Step 2, estimate the quantile at 100 days using four approaches
% % 1, KME
% qan1   = length(find(exprd>100))/length(exprd);
% 
% % 2, KME with averaging 85-115 days window
% tmpq   = intersect(find(exprd>85),find(exprd<115)); 
% np     = length(tmpq);
% if np <2, 
%     qan2   = qan1;
% else
%     qan2   = 0;
%         % calculate the area
%     for i = 1:(np-1),
%         qan2 = qan2 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
%             length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan2 = qan2/(exprd(tmpq(np))-exprd(tmpq(1)));
% end
% 
% % 3, exponential using all the data
% expar1 = mean(exprd);
% stdpar1= std(exprd);
% qan3 = exp(-100/expar1);
% 
% 
% % 4, exponential using data regarding whatever larger than 150 as censored
% expar2 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
%     /length(find(exprd<=150));
% qan4 = exp(-100/expar2);
% 
% 
% % 5, KME with 10 days window from days 95 to 105.
% tmpq2   = intersect(find(exprd>95),find(exprd<105)); 
% np2      = length(tmpq2);
% if np2 < 2, 
%     qan5   = qan1;
% else
%     qan5   = 0;
%         % calculate the area
%     for i = 1:(np2-1),
%         qan5 = qan5 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
%             length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan5 = qan5/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
% end
% 
% 
% % 6, gamma distribution
% % find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% % and b0 = ave(x)/a0;
% lgavex = log(mean(exprd));
% avelgx = mean(log(exprd));
% a0     = 0.5/(lgavex-avelgx);
% %options = optimset('Display','off','Simplex','on');
% options = optimset('Display','off','LargeScale','on',...
%                    'TolFun',.0001, 'Algorithm','active-set');
% [a fval]   = fmincon(@(a) gamlike([a; mean(exprd)/a],exprd),a0,...
%     [],[],[],[],0.00001,10000,[],options); 
% qan6 = 1-gamcdf(100,a,mean(exprd)/a);
% 
% 
% % 7, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
% qan7 = 1-wblcdf(100,wei(1),wei(2));
% 
% 
% % 8, lognormal distribution
% % m     = expar1;
% % v     = stdpar1;
% % lmu   = log(m^2/sqrt(v+m^2));
% % lsigma= sqrt(log(v/m^2 + 1)); 
% % [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
% %     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
% logparhat = lognfit(exprd);
% qan8 = 1-logncdf(100,logparhat(1),logparhat(2));
% % lnorm
% % qan8
% 
% % try the exponential distribution with multiple ratio of time and 
% %  exponential parameter
% % 150/120 = 1.25
% % try ratio = 0.5, 0.75,  1, 1.5, 1.75, 2
% %    time   = 60,   90, 120, 180, 210,  240
% 
% % 9, exponential using data regarding whatever larger than 60 as censored
% expar9 = (sum(exprd(find(exprd<=60)))+length(find(exprd>60))*60)...
%     /length(find(exprd<=60));
% qan9 = exp(-100/expar9);
% 
% 
% % 10, exponential using data regarding whatever larger than 90 as censored
% expar10 = (sum(exprd(find(exprd<=90)))+length(find(exprd>90))*90)...
%     /length(find(exprd<=90));
% qan10 = exp(-100/expar10);
% 
% 
% % 11, exponential using data regarding whatever larger than 120 as censored
% expar11 = (sum(exprd(find(exprd<=120)))+length(find(exprd>120))*120)...
%     /length(find(exprd<=120));
% qan11 = exp(-100/expar11);
% 
% 
% % 12, exponential using data regarding whatever larger than 180 as censored
% expar12 = (sum(exprd(find(exprd<=180)))+length(find(exprd>180))*180)...
%     /length(find(exprd<=180));
% qan12 = exp(-100/expar12);
% 
% 
% % 13, exponential using data regarding whatever larger than 210 as censored
% expar13 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
%     /length(find(exprd<=210));
% qan13 = exp(-100/expar13);
% 
% 
% % 14, pexe estimate
% [time_die14,ttot14,deaths14] = totaltest(exprd,ones(Ns,1));
% lendis = length(ttot14);
% [qan14 lamest14] = pexeest(exprd, ones(Ns,1), time_die14(1:lendis-1), 100);
% 
% 
% % 15, rpexe estimate
% % Estimate the first change point
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
% pexeout     = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
% if min(pexeout.pvalues) < crival, % if the change point is found
%     tchange = pexeout.times(find(pexeout.pvalues < crival));
%     [qan15 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
%     rpexrec(iter) = min(pexeout.pvalues);
% else 
%     expar15 = expar1;
%     qan15   = qan3;
% end
% 
% 
% kme1(iter) = qan1;
% kme2(iter) = qan2;
% exp1(iter) = qan3;
% exp2(iter) = qan4;
% kme3(iter) = qan5;
% gam1(iter) = qan6;
% wei1(iter) = qan7;
% ln1(iter)  = qan8;
% exp9(iter) = qan9;
% exp10(iter) = qan10;
% exp11(iter) = qan11;
% exp12(iter) = qan12;
% exp13(iter) = qan13;
% exp14(iter) = qan14;
% exp15(iter) = qan15;
% end
% 
% 
% % 
% quanval = 1-gamcdf(100,0.9,133);
% 
% % compute the roon mean squared prediction error
% rmspe1 = sqrt(sum((kme1-quanval)'*(kme1-quanval)/N))
% rmspe2 = sqrt(sum((kme2-quanval)'*(kme2-quanval)/N))
% rmspe3 = sqrt(sum((exp1-quanval)'*(exp1-quanval)/N))
% rmspe4 = sqrt(sum((exp2-quanval)'*(exp2-quanval)/N))
% rmspe5 = sqrt(sum((kme3-quanval)'*(kme3-quanval)/N))
% rmspe6 = sqrt(sum((gam1-quanval)'*(gam1-quanval)/N))
% rmspe7 = sqrt(sum((wei1-quanval)'*(wei1-quanval)/N))
% rmspe8 = sqrt(sum((ln1-quanval)'*(ln1-quanval)/N))
% rmspe9 = sqrt(sum((exp9-quanval)'*(exp9-quanval)/N))
% rmspe10 = sqrt(sum((exp10-quanval)'*(exp10-quanval)/N))
% rmspe11 = sqrt(sum((exp11-quanval)'*(exp11-quanval)/N))
% rmspe12 = sqrt(sum((exp12-quanval)'*(exp12-quanval)/N))
% rmspe13 = sqrt(sum((exp13-quanval)'*(exp13-quanval)/N))
% rmspe14 = sqrt(sum((exp14-quanval)'*(exp14-quanval)/N))
% rmspe15 = sqrt(sum((exp15-quanval)'*(exp15-quanval)/N))
% 
% 
% res30_gam2 = [[var(kme1) mean(kme1)-quanval rmspe1]
% [var(kme2) mean(kme2)-quanval rmspe2]
% [var(kme3) mean(kme3)-quanval rmspe5]
% [var(gam1) mean(gam1)-quanval rmspe6]
% [var(wei1) mean(wei1)-quanval rmspe7]
% [var(ln1) mean(ln1)-quanval rmspe8]
% [var(exp1) mean(exp1)-quanval rmspe3]
% [var(exp9) mean(exp9)-quanval rmspe9]
% [var(exp10) mean(exp10)-quanval rmspe10]
% [var(exp11) mean(exp11)-quanval rmspe11]
% [var(exp2) mean(exp2)-quanval rmspe4]
% [var(exp12) mean(exp12)-quanval rmspe12]
% [var(exp13) mean(exp13)-quanval rmspe13]
% [var(exp14) mean(exp14)-quanval rmspe14]
% [var(exp15) mean(exp15)-quanval rmspe15]]
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % We fit the functions by looking maximizing the likelihood functions.
% % There is no censored observations in any of the cases
% 
% % set the number of iteration samples
% N = 5000;
% 
% kme1 = zeros(N,1);
% kme2 = zeros(N,1);
% exp1 = zeros(N,1);
% exp2 = zeros(N,1);
% kme3 = zeros(N,1);
% gam1 = zeros(N,1);
% wei1 = zeros(N,1);
% ln1  = zeros(N,1);
% exp9 = zeros(N,1);
% exp10 = zeros(N,1);
% exp11 = zeros(N,1);
% exp12 = zeros(N,1);
% exp13 = zeros(N,1);
% exp14 = zeros(N,1);
% exp15 = zeros(N,1);
% 
% rpexrec = zeros(N,1);
% 
% Ns   = 30;
% crival = exp(-4.2331-0.3938*log(Ns));
% 
% for iter = 1:N,
%     
% % Step 1, generate data from exponential distribution 120, with sample size 50
% exprd = sort(gamrnd(1.25,96,Ns,1));
% 
% % Step 2, estimate the quantile at 100 days using four approaches
% % 1, KME
% qan1   = length(find(exprd>100))/length(exprd);
% 
% % 2, KME with averaging 85-115 days window
% tmpq   = intersect(find(exprd>85),find(exprd<115)); 
% np     = length(tmpq);
% if np <2, 
%     qan2   = qan1;
% else
%     qan2   = 0;
%         % calculate the area
%     for i = 1:(np-1),
%         qan2 = qan2 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
%             length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan2 = qan2/(exprd(tmpq(np))-exprd(tmpq(1)));
% end
% 
% % 3, exponential using all the data
% expar1 = mean(exprd);
% stdpar1= std(exprd);
% qan3 = exp(-100/expar1);
% 
% 
% % 4, exponential using data regarding whatever larger than 150 as censored
% expar2 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
%     /length(find(exprd<=150));
% qan4 = exp(-100/expar2);
% 
% 
% % 5, KME with 10 days window from days 95 to 105.
% tmpq2   = intersect(find(exprd>95),find(exprd<105)); 
% np2      = length(tmpq2);
% if np2 < 2, 
%     qan5   = qan1;
% else
%     qan5   = 0;
%         % calculate the area
%     for i = 1:(np2-1),
%         qan5 = qan5 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
%             length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan5 = qan5/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
% end
% 
% 
% % 6, gamma distribution
% % find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% % and b0 = ave(x)/a0;
% lgavex = log(mean(exprd));
% avelgx = mean(log(exprd));
% a0     = 0.5/(lgavex-avelgx);
% %options = optimset('Display','off','Simplex','on');
% options = optimset('Display','off','LargeScale','on',...
%                    'TolFun',.0001, 'Algorithm','active-set');
% [a fval]   = fmincon(@(a) gamlike([a; mean(exprd)/a],exprd),a0,...
%     [],[],[],[],0.00001,10000,[],options); 
% qan6 = 1-gamcdf(100,a,mean(exprd)/a);
% 
% 
% % 7, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
% qan7 = 1-wblcdf(100,wei(1),wei(2));
% 
% 
% % 8, lognormal distribution
% % m     = expar1;
% % v     = stdpar1;
% % lmu   = log(m^2/sqrt(v+m^2));
% % lsigma= sqrt(log(v/m^2 + 1)); 
% % [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
% %     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
% logparhat = lognfit(exprd);
% qan8 = 1-logncdf(100,logparhat(1),logparhat(2));
% % lnorm
% % qan8
% 
% % try the exponential distribution with multiple ratio of time and 
% %  exponential parameter
% % 150/120 = 1.25
% % try ratio = 0.5, 0.75,  1, 1.5, 1.75, 2
% %    time   = 60,   90, 120, 180, 210,  240
% 
% % 9, exponential using data regarding whatever larger than 60 as censored
% expar9 = (sum(exprd(find(exprd<=60)))+length(find(exprd>60))*60)...
%     /length(find(exprd<=60));
% qan9 = exp(-100/expar9);
% 
% 
% % 10, exponential using data regarding whatever larger than 90 as censored
% expar10 = (sum(exprd(find(exprd<=90)))+length(find(exprd>90))*90)...
%     /length(find(exprd<=90));
% qan10 = exp(-100/expar10);
% 
% 
% % 11, exponential using data regarding whatever larger than 120 as censored
% expar11 = (sum(exprd(find(exprd<=120)))+length(find(exprd>120))*120)...
%     /length(find(exprd<=120));
% qan11 = exp(-100/expar11);
% 
% 
% % 12, exponential using data regarding whatever larger than 180 as censored
% expar12 = (sum(exprd(find(exprd<=180)))+length(find(exprd>180))*180)...
%     /length(find(exprd<=180));
% qan12 = exp(-100/expar12);
% 
% 
% % 13, exponential using data regarding whatever larger than 210 as censored
% expar13 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
%     /length(find(exprd<=210));
% qan13 = exp(-100/expar13);
% 
% 
% % 14, pexe estimate
% [time_die14,ttot14,deaths14] = totaltest(exprd,ones(Ns,1));
% lendis = length(ttot14);
% [qan14 lamest14] = pexeest(exprd, ones(Ns,1), time_die14(1:lendis-1), 100);
% 
% 
% % 15, rpexe estimate
% % Estimate the first change point
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
% pexeout     = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
% if min(pexeout.pvalues) < crival, % if the change point is found
%     tchange = pexeout.times(find(pexeout.pvalues < crival));
%     [qan15 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
%     rpexrec(iter) = min(pexeout.pvalues);
% else 
%     expar15 = expar1;
%     qan15   = qan3;
% end
% 
% 
% kme1(iter) = qan1;
% kme2(iter) = qan2;
% exp1(iter) = qan3;
% exp2(iter) = qan4;
% kme3(iter) = qan5;
% gam1(iter) = qan6;
% wei1(iter) = qan7;
% ln1(iter)  = qan8;
% exp9(iter) = qan9;
% exp10(iter) = qan10;
% exp11(iter) = qan11;
% exp12(iter) = qan12;
% exp13(iter) = qan13;
% exp14(iter) = qan14;
% exp15(iter) = qan15;
% end
% 
% 
% % 
% quanval = 1-gamcdf(100,1.25,96);
% 
% % compute the roon mean squared prediction error
% rmspe1 = sqrt(sum((kme1-quanval)'*(kme1-quanval)/N))
% rmspe2 = sqrt(sum((kme2-quanval)'*(kme2-quanval)/N))
% rmspe3 = sqrt(sum((exp1-quanval)'*(exp1-quanval)/N))
% rmspe4 = sqrt(sum((exp2-quanval)'*(exp2-quanval)/N))
% rmspe5 = sqrt(sum((kme3-quanval)'*(kme3-quanval)/N))
% rmspe6 = sqrt(sum((gam1-quanval)'*(gam1-quanval)/N))
% rmspe7 = sqrt(sum((wei1-quanval)'*(wei1-quanval)/N))
% rmspe8 = sqrt(sum((ln1-quanval)'*(ln1-quanval)/N))
% rmspe9 = sqrt(sum((exp9-quanval)'*(exp9-quanval)/N))
% rmspe10 = sqrt(sum((exp10-quanval)'*(exp10-quanval)/N))
% rmspe11 = sqrt(sum((exp11-quanval)'*(exp11-quanval)/N))
% rmspe12 = sqrt(sum((exp12-quanval)'*(exp12-quanval)/N))
% rmspe13 = sqrt(sum((exp13-quanval)'*(exp13-quanval)/N))
% rmspe14 = sqrt(sum((exp14-quanval)'*(exp14-quanval)/N))
% rmspe15 = sqrt(sum((exp15-quanval)'*(exp15-quanval)/N))
% 
% 
% res30_gam25 = [[var(kme1) mean(kme1)-quanval rmspe1]
% [var(kme2) mean(kme2)-quanval rmspe2]
% [var(kme3) mean(kme3)-quanval rmspe5]
% [var(gam1) mean(gam1)-quanval rmspe6]
% [var(wei1) mean(wei1)-quanval rmspe7]
% [var(ln1) mean(ln1)-quanval rmspe8]
% [var(exp1) mean(exp1)-quanval rmspe3]
% [var(exp9) mean(exp9)-quanval rmspe9]
% [var(exp10) mean(exp10)-quanval rmspe10]
% [var(exp11) mean(exp11)-quanval rmspe11]
% [var(exp2) mean(exp2)-quanval rmspe4]
% [var(exp12) mean(exp12)-quanval rmspe12]
% [var(exp13) mean(exp13)-quanval rmspe13]
% [var(exp14) mean(exp14)-quanval rmspe14]
% [var(exp15) mean(exp15)-quanval rmspe15]]
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % We fit the functions by looking maximizing the likelihood functions.
% % There is no censored observations in any of the cases
% 
% % set the number of iteration samples
% N = 5000;
% 
% kme1 = zeros(N,1);
% kme2 = zeros(N,1);
% exp1 = zeros(N,1);
% exp2 = zeros(N,1);
% kme3 = zeros(N,1);
% gam1 = zeros(N,1);
% wei1 = zeros(N,1);
% ln1  = zeros(N,1);
% exp9 = zeros(N,1);
% exp10 = zeros(N,1);
% exp11 = zeros(N,1);
% exp12 = zeros(N,1);
% exp13 = zeros(N,1);
% exp14 = zeros(N,1);
% exp15 = zeros(N,1);
% 
% rpexrec = zeros(N,1);
% 
% Ns   = 30;
% crival = exp(-4.2331-0.3938*log(Ns));
% 
% for iter = 1:N,
%     
% % Step 1, generate data from exponential distribution 120, with sample size 50
% exprd = sort(gamrnd(2,60,Ns,1));
% 
% % Step 2, estimate the quantile at 100 days using four approaches
% % 1, KME
% qan1   = length(find(exprd>100))/length(exprd);
% 
% % 2, KME with averaging 85-115 days window
% tmpq   = intersect(find(exprd>85),find(exprd<115)); 
% np     = length(tmpq);
% if np <2, 
%     qan2   = qan1;
% else
%     qan2   = 0;
%         % calculate the area
%     for i = 1:(np-1),
%         qan2 = qan2 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
%             length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan2 = qan2/(exprd(tmpq(np))-exprd(tmpq(1)));
% end
% 
% % 3, exponential using all the data
% expar1 = mean(exprd);
% stdpar1= std(exprd);
% qan3 = exp(-100/expar1);
% 
% 
% % 4, exponential using data regarding whatever larger than 150 as censored
% expar2 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
%     /length(find(exprd<=150));
% qan4 = exp(-100/expar2);
% 
% 
% % 5, KME with 10 days window from days 95 to 105.
% tmpq2   = intersect(find(exprd>95),find(exprd<105)); 
% np2      = length(tmpq2);
% if np2 < 2, 
%     qan5   = qan1;
% else
%     qan5   = 0;
%         % calculate the area
%     for i = 1:(np2-1),
%         qan5 = qan5 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
%             length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan5 = qan5/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
% end
% 
% 
% % 6, gamma distribution
% % find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% % and b0 = ave(x)/a0;
% lgavex = log(mean(exprd));
% avelgx = mean(log(exprd));
% a0     = 0.5/(lgavex-avelgx);
% %options = optimset('Display','off','Simplex','on');
% options = optimset('Display','off','LargeScale','on',...
%                    'TolFun',.0001, 'Algorithm','active-set');
% [a fval]   = fmincon(@(a) gamlike([a; mean(exprd)/a],exprd),a0,...
%     [],[],[],[],0.00001,10000,[],options); 
% qan6 = 1-gamcdf(100,a,mean(exprd)/a);
% 
% 
% % 7, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
% qan7 = 1-wblcdf(100,wei(1),wei(2));
% 
% 
% % 8, lognormal distribution
% % m     = expar1;
% % v     = stdpar1;
% % lmu   = log(m^2/sqrt(v+m^2));
% % lsigma= sqrt(log(v/m^2 + 1)); 
% % [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
% %     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
% logparhat = lognfit(exprd);
% qan8 = 1-logncdf(100,logparhat(1),logparhat(2));
% % lnorm
% % qan8
% 
% % try the exponential distribution with multiple ratio of time and 
% %  exponential parameter
% % 150/120 = 1.25
% % try ratio = 0.5, 0.75,  1, 1.5, 1.75, 2
% %    time   = 60,   90, 120, 180, 210,  240
% 
% % 9, exponential using data regarding whatever larger than 60 as censored
% expar9 = (sum(exprd(find(exprd<=60)))+length(find(exprd>60))*60)...
%     /length(find(exprd<=60));
% qan9 = exp(-100/expar9);
% 
% 
% % 10, exponential using data regarding whatever larger than 90 as censored
% expar10 = (sum(exprd(find(exprd<=90)))+length(find(exprd>90))*90)...
%     /length(find(exprd<=90));
% qan10 = exp(-100/expar10);
% 
% 
% % 11, exponential using data regarding whatever larger than 120 as censored
% expar11 = (sum(exprd(find(exprd<=120)))+length(find(exprd>120))*120)...
%     /length(find(exprd<=120));
% qan11 = exp(-100/expar11);
% 
% 
% % 12, exponential using data regarding whatever larger than 180 as censored
% expar12 = (sum(exprd(find(exprd<=180)))+length(find(exprd>180))*180)...
%     /length(find(exprd<=180));
% qan12 = exp(-100/expar12);
% 
% 
% % 13, exponential using data regarding whatever larger than 210 as censored
% expar13 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
%     /length(find(exprd<=210));
% qan13 = exp(-100/expar13);
% 
% 
% % 14, pexe estimate
% [time_die14,ttot14,deaths14] = totaltest(exprd,ones(Ns,1));
% lendis = length(ttot14);
% [qan14 lamest14] = pexeest(exprd, ones(Ns,1), time_die14(1:lendis-1), 100);
% 
% 
% % 15, rpexe estimate
% % Estimate the first change point
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
% pexeout     = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
% if min(pexeout.pvalues) < crival, % if the change point is found
%     tchange = pexeout.times(find(pexeout.pvalues < crival));
%     [qan15 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
%     rpexrec(iter) = min(pexeout.pvalues);
% else 
%     expar15 = expar1;
%     qan15   = qan3;
% end
% 
% 
% kme1(iter) = qan1;
% kme2(iter) = qan2;
% exp1(iter) = qan3;
% exp2(iter) = qan4;
% kme3(iter) = qan5;
% gam1(iter) = qan6;
% wei1(iter) = qan7;
% ln1(iter)  = qan8;
% exp9(iter) = qan9;
% exp10(iter) = qan10;
% exp11(iter) = qan11;
% exp12(iter) = qan12;
% exp13(iter) = qan13;
% exp14(iter) = qan14;
% exp15(iter) = qan15;
% end
% 
% 
% % 
% quanval = 1-gamcdf(100,2,60);
% 
% % compute the roon mean squared prediction error
% rmspe1 = sqrt(sum((kme1-quanval)'*(kme1-quanval)/N))
% rmspe2 = sqrt(sum((kme2-quanval)'*(kme2-quanval)/N))
% rmspe3 = sqrt(sum((exp1-quanval)'*(exp1-quanval)/N))
% rmspe4 = sqrt(sum((exp2-quanval)'*(exp2-quanval)/N))
% rmspe5 = sqrt(sum((kme3-quanval)'*(kme3-quanval)/N))
% rmspe6 = sqrt(sum((gam1-quanval)'*(gam1-quanval)/N))
% rmspe7 = sqrt(sum((wei1-quanval)'*(wei1-quanval)/N))
% rmspe8 = sqrt(sum((ln1-quanval)'*(ln1-quanval)/N))
% rmspe9 = sqrt(sum((exp9-quanval)'*(exp9-quanval)/N))
% rmspe10 = sqrt(sum((exp10-quanval)'*(exp10-quanval)/N))
% rmspe11 = sqrt(sum((exp11-quanval)'*(exp11-quanval)/N))
% rmspe12 = sqrt(sum((exp12-quanval)'*(exp12-quanval)/N))
% rmspe13 = sqrt(sum((exp13-quanval)'*(exp13-quanval)/N))
% rmspe14 = sqrt(sum((exp14-quanval)'*(exp14-quanval)/N))
% rmspe15 = sqrt(sum((exp15-quanval)'*(exp15-quanval)/N))
% 
% 
% res30_gam3 = [[var(kme1) mean(kme1)-quanval rmspe1]
% [var(kme2) mean(kme2)-quanval rmspe2]
% [var(kme3) mean(kme3)-quanval rmspe5]
% [var(gam1) mean(gam1)-quanval rmspe6]
% [var(wei1) mean(wei1)-quanval rmspe7]
% [var(ln1) mean(ln1)-quanval rmspe8]
% [var(exp1) mean(exp1)-quanval rmspe3]
% [var(exp9) mean(exp9)-quanval rmspe9]
% [var(exp10) mean(exp10)-quanval rmspe10]
% [var(exp11) mean(exp11)-quanval rmspe11]
% [var(exp2) mean(exp2)-quanval rmspe4]
% [var(exp12) mean(exp12)-quanval rmspe12]
% [var(exp13) mean(exp13)-quanval rmspe13]
% [var(exp14) mean(exp14)-quanval rmspe14]
% [var(exp15) mean(exp15)-quanval rmspe15]]
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%% Compare the Monte Carlo variances of four predictors in 
% %%% estimating the quantile in the day 100.
% 
% 
% % 1, Kaplan-Meier binomial distribution
% % 2, Kaplan-Meier averaging over length 30
% % 3, Exponential all the way
% % 4, Exponential up to a point 150
% 
% % There are more adds on in this version of the code.
% % 5, Kaplan-Meier averaging over length 10
% % 6, gamma
% % 7, lognormal
% % 8, weibull
% % 9, Exponential up to a point 60
% % 10, Exponential up to a point 90
% % 11, Exponential up to a point 120
% % 12, Exponential up to a point 180
% % 13, Exponential up to a point 210
% % 14, RPEXE monotone restriction with alpha^star = 0.005
% 
% 
% % 15, Gompertz
% % 16, generalized gamma
% 
% % We fit the functions by looking maximizing the likelihood functions.
% % There is no censored observations in any of the cases
% 
% % set the number of iteration samples
% N = 5000;
% 
% kme1 = zeros(N,1);
% kme2 = zeros(N,1);
% exp1 = zeros(N,1);
% exp2 = zeros(N,1);
% kme3 = zeros(N,1);
% gam1 = zeros(N,1);
% wei1 = zeros(N,1);
% ln1  = zeros(N,1);
% exp9 = zeros(N,1);
% exp10 = zeros(N,1);
% exp11 = zeros(N,1);
% exp12 = zeros(N,1);
% exp13 = zeros(N,1);
% exp14 = zeros(N,1);
% exp15 = zeros(N,1);
% 
% rpexrec = zeros(N,1);
% 
% Ns   = 30;
% crival = exp(-4.2331-0.3938*log(Ns));
% 
% for iter = 1:N,
%     
% % Step 1, generate data from exponential distribution 120, with sample size 50
% exprd = sort(gamrnd(10,12,Ns,1));
% 
% % Step 2, estimate the quantile at 100 days using four approaches
% % 1, KME
% qan1   = length(find(exprd>100))/length(exprd);
% 
% % 2, KME with averaging 85-115 days window
% tmpq   = intersect(find(exprd>85),find(exprd<115)); 
% np     = length(tmpq);
% if np <2, 
%     qan2   = qan1;
% else
%     qan2   = 0;
%         % calculate the area
%     for i = 1:(np-1),
%         qan2 = qan2 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
%             length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan2 = qan2/(exprd(tmpq(np))-exprd(tmpq(1)));
% end
% 
% % 3, exponential using all the data
% expar1 = mean(exprd);
% stdpar1= std(exprd);
% qan3 = exp(-100/expar1);
% 
% 
% % 4, exponential using data regarding whatever larger than 150 as censored
% expar2 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
%     /length(find(exprd<=150));
% qan4 = exp(-100/expar2);
% 
% 
% % 5, KME with 10 days window from days 95 to 105.
% tmpq2   = intersect(find(exprd>95),find(exprd<105)); 
% np2      = length(tmpq2);
% if np2 < 2, 
%     qan5   = qan1;
% else
%     qan5   = 0;
%         % calculate the area
%     for i = 1:(np2-1),
%         qan5 = qan5 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
%             length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
%     end
%         % calculate the weighted average of the quantile
%     qan5 = qan5/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
% end
% 
% 
% % 6, gamma distribution
% % find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% % and b0 = ave(x)/a0;
% lgavex = log(mean(exprd));
% avelgx = mean(log(exprd));
% a0     = 0.5/(lgavex-avelgx);
% %options = optimset('Display','off','Simplex','on');
% options = optimset('Display','off','LargeScale','on',...
%                    'TolFun',.0001, 'Algorithm','active-set');
% [a fval]   = fmincon(@(a) gamlike([a; mean(exprd)/a],exprd),a0,...
%     [],[],[],[],0.00001,10000,[],options); 
% qan6 = 1-gamcdf(100,a,mean(exprd)/a);
% 
% 
% % 7, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
% qan7 = 1-wblcdf(100,wei(1),wei(2));
% 
% 
% % 8, lognormal distribution
% % m     = expar1;
% % v     = stdpar1;
% % lmu   = log(m^2/sqrt(v+m^2));
% % lsigma= sqrt(log(v/m^2 + 1)); 
% % [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
% %     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
% logparhat = lognfit(exprd);
% qan8 = 1-logncdf(100,logparhat(1),logparhat(2));
% % lnorm
% % qan8
% 
% % try the exponential distribution with multiple ratio of time and 
% %  exponential parameter
% % 150/120 = 1.25
% % try ratio = 0.5, 0.75,  1, 1.5, 1.75, 2
% %    time   = 60,   90, 120, 180, 210,  240
% 
% % 9, exponential using data regarding whatever larger than 60 as censored
% expar9 = (sum(exprd(find(exprd<=60)))+length(find(exprd>60))*60)...
%     /length(find(exprd<=60));
% qan9 = exp(-100/expar9);
% 
% 
% % 10, exponential using data regarding whatever larger than 90 as censored
% expar10 = (sum(exprd(find(exprd<=90)))+length(find(exprd>90))*90)...
%     /length(find(exprd<=90));
% qan10 = exp(-100/expar10);
% 
% 
% % 11, exponential using data regarding whatever larger than 120 as censored
% expar11 = (sum(exprd(find(exprd<=120)))+length(find(exprd>120))*120)...
%     /length(find(exprd<=120));
% qan11 = exp(-100/expar11);
% 
% 
% % 12, exponential using data regarding whatever larger than 180 as censored
% expar12 = (sum(exprd(find(exprd<=180)))+length(find(exprd>180))*180)...
%     /length(find(exprd<=180));
% qan12 = exp(-100/expar12);
% 
% 
% % 13, exponential using data regarding whatever larger than 210 as censored
% expar13 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
%     /length(find(exprd<=210));
% qan13 = exp(-100/expar13);
% 
% 
% % 14, pexe estimate
% [time_die14,ttot14,deaths14] = totaltest(exprd,ones(Ns,1));
% lendis = length(ttot14);
% [qan14 lamest14] = pexeest(exprd, ones(Ns,1), time_die14(1:lendis-1), 100);
% 
% 
% % 15, rpexe estimate
% % Estimate the first change point
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
% addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
% pexeout     = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
% if min(pexeout.pvalues) < crival, % if the change point is found
%     tchange = pexeout.times(find(pexeout.pvalues < crival));
%     [qan15 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
%     rpexrec(iter) = min(pexeout.pvalues);
% else 
%     expar15 = expar1;
%     qan15   = qan3;
% end
% 
% 
% kme1(iter) = qan1;
% kme2(iter) = qan2;
% exp1(iter) = qan3;
% exp2(iter) = qan4;
% kme3(iter) = qan5;
% gam1(iter) = qan6;
% wei1(iter) = qan7;
% ln1(iter)  = qan8;
% exp9(iter) = qan9;
% exp10(iter) = qan10;
% exp11(iter) = qan11;
% exp12(iter) = qan12;
% exp13(iter) = qan13;
% exp14(iter) = qan14;
% exp15(iter) = qan15;
% end
% 
% 
% % 
% quanval = 1-gamcdf(100,10,12);
% 
% % compute the roon mean squared prediction error
% rmspe1 = sqrt(sum((kme1-quanval)'*(kme1-quanval)/N))
% rmspe2 = sqrt(sum((kme2-quanval)'*(kme2-quanval)/N))
% rmspe3 = sqrt(sum((exp1-quanval)'*(exp1-quanval)/N))
% rmspe4 = sqrt(sum((exp2-quanval)'*(exp2-quanval)/N))
% rmspe5 = sqrt(sum((kme3-quanval)'*(kme3-quanval)/N))
% rmspe6 = sqrt(sum((gam1-quanval)'*(gam1-quanval)/N))
% rmspe7 = sqrt(sum((wei1-quanval)'*(wei1-quanval)/N))
% rmspe8 = sqrt(sum((ln1-quanval)'*(ln1-quanval)/N))
% rmspe9 = sqrt(sum((exp9-quanval)'*(exp9-quanval)/N))
% rmspe10 = sqrt(sum((exp10-quanval)'*(exp10-quanval)/N))
% rmspe11 = sqrt(sum((exp11-quanval)'*(exp11-quanval)/N))
% rmspe12 = sqrt(sum((exp12-quanval)'*(exp12-quanval)/N))
% rmspe13 = sqrt(sum((exp13-quanval)'*(exp13-quanval)/N))
% rmspe14 = sqrt(sum((exp14-quanval)'*(exp14-quanval)/N))
% rmspe15 = sqrt(sum((exp15-quanval)'*(exp15-quanval)/N))
% 
% 
% res30_gam4 = [[var(kme1) mean(kme1)-quanval rmspe1]
% [var(kme2) mean(kme2)-quanval rmspe2]
% [var(kme3) mean(kme3)-quanval rmspe5]
% [var(gam1) mean(gam1)-quanval rmspe6]
% [var(wei1) mean(wei1)-quanval rmspe7]
% [var(ln1) mean(ln1)-quanval rmspe8]
% [var(exp1) mean(exp1)-quanval rmspe3]
% [var(exp9) mean(exp9)-quanval rmspe9]
% [var(exp10) mean(exp10)-quanval rmspe10]
% [var(exp11) mean(exp11)-quanval rmspe11]
% [var(exp2) mean(exp2)-quanval rmspe4]
% [var(exp12) mean(exp12)-quanval rmspe12]
% [var(exp13) mean(exp13)-quanval rmspe13]
% [var(exp14) mean(exp14)-quanval rmspe14]
% [var(exp15) mean(exp15)-quanval rmspe15]]
% 
% 
% 
% 
% 
% 
