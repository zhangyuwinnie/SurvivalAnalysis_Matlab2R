%%% Compare the Monte Carlo variances of four predictors in 
%%% estimating the quantile in the day 100.


% 1, Kaplan-Meier binomial distribution
% 2, Kaplan-Meier averaging over length 30
% 3, Exponential all the way
% 4, Exponential up to a point

% There are more adds on in this version of the code.
% 5, Kaplan-Meier averaging over length 10
% 6, gamma
% 7, lognormal
% 8, weibull
% 9, Gompertz
% 10, generalized gamma

% We fit the functions by looking maximizing the likelihood functions.
% There is no censored observations in any of the cases

% set the number of iteration samples
N = 5000;

kme1 = zeros(N,1);
kme2 = zeros(N,1);
exp1 = zeros(N,1);
exp2 = zeros(N,1);
kme3 = zeros(N,1);
gam1 = zeros(N,1);
wei1 = zeros(N,1);
ln1  = zeros(N,1);

Ns   = 15;
for iter = 1:N,
    
% Step 1, generate data from exponential distribution 120, with sample size 50
exprnd = sort(gamrnd(1,120,Ns,1));

% Step 2, estimate the quantile at 100 days using four approaches
% 1, KME
qan1   = length(find(exprnd>100))/length(exprnd);

% 2, KME with averaging 85-115 days window
tmpq   = intersect(find(exprnd>85),find(exprnd<115)); 
np     = length(tmpq);
if np <2, 
    qan2   = qan1;
else
    qan2   = 0;
        % calculate the area
    for i = 1:(np-1),
        qan2 = qan2 + (exprnd(tmpq(i+1))-exprnd(tmpq(i)))*...
            length(find(exprnd>exprnd(tmpq(i+1))))/length(exprnd);
    end
        % calculate the weighted average of the quantile
    qan2 = qan2/(exprnd(tmpq(np))-exprnd(tmpq(1)));
end

% 3, exponential using all the data
expar1 = mean(exprnd);
stdpar1= std(exprnd);
qan3 = exp(-100/expar1);


% 4, exponential using data regarding whatever larger than 150 as censored
expar2 = (sum(exprnd(find(exprnd<=150)))+length(find(exprnd>150))*150)...
    /length(find(exprnd<=150));
qan4 = exp(-100/expar2);


% 5, KME with 10 days window from days 95 to 105.
tmpq2   = intersect(find(exprnd>95),find(exprnd<105)); 
np2      = length(tmpq2);
if np2 < 2, 
    qan5   = qan1;
else
    qan5   = 0;
        % calculate the area
    for i = 1:(np2-1),
        qan5 = qan5 + (exprnd(tmpq2(i+1))-exprnd(tmpq2(i)))*...
            length(find(exprnd>exprnd(tmpq2(i+1))))/length(exprnd);
    end
        % calculate the weighted average of the quantile
    qan5 = qan5/(exprnd(tmpq2(np2))-exprnd(tmpq2(1)));
end


% 6, gamma distribution
% find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% and b0 = ave(x)/a0;
lgavex = log(mean(exprnd));
avelgx = mean(log(exprnd));
a0     = 0.5/(lgavex-avelgx);
%options = optimset('Display','off','Simplex','on');
options = optimset('Display','off','LargeScale','on',...
                   'TolFun',.0001, 'Algorithm','active-set');
[a fval]   = fmincon(@(a) gamlike([a; mean(exprnd)/a],exprnd),a0,...
    [],[],[],[],0.00001,10000,[],options); 
qan6 = 1-gamcdf(100,a,mean(exprnd)/a);



% 7, weibull distribution
wa0     = expar1;
wb0     = 1;
[wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprnd),...
    [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
qan7 = 1-wblcdf(100,wei(1),wei(2));



% 8, lognormal distribution
m     = expar1;
v     = stdpar1;
lmu   = log(m^2/sqrt(v+m^2));
lsigma= sqrt(log(v/m^2 + 1)); 
[lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprnd),...
    [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
qan8 = 1-logncdf(100,lnorm(1),lnorm(2));
% lnorm
% qan8

kme1(iter) = qan1;
kme2(iter) = qan2;
exp1(iter) = qan3;
exp2(iter) = qan4;
kme3(iter) = qan5;
gam1(iter) = qan6;
wei1(iter) = qan7;
ln1(iter)  = qan8;
end;


% true quantile = exp(-100/120)=0.4346
subplot(4,2,1);
hist(kme1);
title('Kaplan-Meier Estimate');
xlim([0.1,0.9]);
subplot(4,2,2);
hist(kme2);
title('Kaplan-Meier Estimate by averaging over window of 30');
xlim([0.1,0.9]);
subplot(4,2,3);
hist(exp1);
title('Exponential Estimate');
xlim([0.1,0.9]);
subplot(4,2,4);
hist(exp2);
title('Exponential Estimate with censoring at 150');
xlim([0.1,0.9]);
subplot(4,2,5);
hist(kme3);
title('Kaplan-Meier Estimate by averaging over window of 10');
xlim([0.1,0.9]);
subplot(4,2,6);
hist(gam1);
title('Gamma Estimate');
xlim([0.1,0.9]);
subplot(4,2,7);
hist(wei1);
title('Weibull Estimate');
xlim([0.1,0.9]);
subplot(4,2,8);
hist(ln1);
title('Lognormal Estimate');
xlim([0.1,0.9]);


% compute the roon mean squared prediction error
rmspe1 = sqrt(sum((kme1-0.4346)'*(kme1-0.4346)/N))
rmspe2 = sqrt(sum((kme2-0.4346)'*(kme2-0.4346)/N))
rmspe3 = sqrt(sum((exp1-0.4346)'*(exp1-0.4346)/N))
rmspe4 = sqrt(sum((exp2-0.4346)'*(exp2-0.4346)/N))
rmspe5 = sqrt(sum((kme3-0.4346)'*(kme3-0.4346)/N))
rmspe6 = sqrt(sum((gam1-0.4346)'*(gam1-0.4346)/N))
rmspe7 = sqrt(sum((wei1-0.4346)'*(wei1-0.4346)/N))
rmspe8 = sqrt(sum((ln1-0.4346)'*(ln1-0.4346)/N))



[var(kme1) mean(kme1)-0.4346 rmspe1]
[var(kme2) mean(kme2)-0.4346 rmspe2]
[var(exp1) mean(exp1)-0.4346 rmspe3]
[var(exp2) mean(exp2)-0.4346 rmspe4]
[var(kme3) mean(kme3)-0.4346 rmspe5]
[var(gam1) mean(gam1)-0.4346 rmspe6]
[var(wei1) mean(wei1)-0.4346 rmspe7]
[var(ln1) mean(ln1)-0.4346   rmspe8]



% 
% 
% subplot(2,2,1);
% hist(kme1);
% title('Kaplan-Meier Estimate');
% xlim([0.1,0.9]);
% subplot(2,2,2);
% hist(kme2);
% title('Kaplan-Meier Estimate by averaging over window of 30');
% xlim([0.1,0.9]);
% subplot(2,2,3);
% hist(exp1);
% title('Exponential Estimate');
% xlim([0.1,0.9]);
% subplot(2,2,4);
% hist(exp2);
% title('Exponential Estimate with censoring at 150');
% xlim([0.1,0.9]);
% 
% % compute the roon mean squared prediction error
% rmspe1 = sqrt(sum((kme1-0.4346)'*(kme1-0.4346)/5000))
% rmspe2 = sqrt(sum((kme2-0.4346)'*(kme2-0.4346)/5000))
% rmspe3 = sqrt(sum((exp1-0.4346)'*(exp1-0.4346)/5000))
% rmspe4 = sqrt(sum((exp2-0.4346)'*(exp2-0.4346)/5000))
% 
% 
% [var(kme1) mean(kme1)-0.4346]
% [var(kme2) mean(kme2)-0.4346]
% [var(exp1) mean(exp1)-0.4346]
% [var(exp2) mean(exp2)-0.4346]















