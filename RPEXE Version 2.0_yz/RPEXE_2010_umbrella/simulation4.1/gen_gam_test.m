function res_gam = gen_gam_test(beta,lambda,sigma,Ns)

% N = 10000, Gamma distribution, using a set of estimates

% set the number of iteration samples
N = 10000;

kme1 = zeros(N,1);
kme2 = zeros(N,1);
exp1 = zeros(N,1);
exp2 = zeros(N,1);
kme3 = zeros(N,1);
gam1 = zeros(N,1);
wei1 = zeros(N,1);
ln1  = zeros(N,1);
exp9 = zeros(N,1);
exp10 = zeros(N,1);
exp11 = zeros(N,1);
exp12 = zeros(N,1);
exp13 = zeros(N,1);
exp14 = zeros(N,1);
exp15 = zeros(N,1);
exp16 = zeros(N,1);
exp17 = zeros(N,1);
exp18 = zeros(N,1);
exp19 = zeros(N,1);
exp20 = zeros(N,1);

rpexrec = zeros(N,1);

crivalmono5  = exp(-4.2331-0.3938*log(Ns))
crivalmono10 = exp(-3.4481-0.3847*log(Ns))

crivalum5  = exp(-1.2234-2.7044*log(Ns))
crivalum10 = exp(-1.1951-2.0238*log(Ns))

% approximate 20% by the same proportion assumption
crivalmono20 = crivalmono10 * crivalmono10/crivalmono5;
crivalum20   = crivalum10 * crivalum10/crivalum5;


for iter = 1:N,
    
% Step 1, generate data from exponential distribution 120, with sample size 50
exprd = sort(gamrnd(1/(lambda^2),lambda^2,Ns,1));
for k = 1:Ns,
    exprd(k) = exp(beta)*exprd(k)^(sigma/lambda);
end

% Step 2, estimate the quantile at 100 days using four approaches
% 1, KME
qan1   = length(find(exprd>100))/length(exprd);

% 2, KME with averaging 85-115 days window
tmpq   = intersect(find(exprd>85),find(exprd<115)); 
np     = length(tmpq);
if np <2, 
    qan2   = qan1;
else
    qan2   = 0;
        % calculate the area
    for i = 1:(np-1),
        qan2 = qan2 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
            length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
    end
        % calculate the weighted average of the quantile
    qan2 = qan2/(exprd(tmpq(np))-exprd(tmpq(1)));
end

% 3, exponential using all the data
expar1 = mean(exprd);
stdpar1= std(exprd);
qan3 = exp(-100/expar1);


% 4, exponential using data regarding whatever larger than 150 as censored
expar2 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
    /length(find(exprd<=150));
qan4 = exp(-100/expar2);


% 5, KME with 10 days window from days 95 to 105.
tmpq2   = intersect(find(exprd>95),find(exprd<105)); 
np2      = length(tmpq2);
if np2 < 2, 
    qan5   = qan1;
else
    qan5   = 0;
        % calculate the area
    for i = 1:(np2-1),
        qan5 = qan5 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
            length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
    end
        % calculate the weighted average of the quantile
    qan5 = qan5/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
end


% 6, gamma distribution
% find out a0 = 0.5/(log(ave(x))-ave(log(x)))
% and b0 = ave(x)/a0;
% lgavex = log(mean(exprd));
% avelgx = mean(log(exprd));
% a0     = 0.5/(lgavex-avelgx);
% %options = optimset('Display','off','Simplex','on');
% options = optimset('Display','off','LargeScale','on',...
%                    'TolFun',.0001, 'Algorithm','active-set');
% [a fval]   = fmincon(@(a) gamlike([a; mean(exprd)/a],exprd),a0,...
%     [],[],[],[],0.00001,10000,[],options); 
gmfit = gamfit(exprd);
qan6 = 1-gamcdf(100,gmfit(1),gmfit(2));




% 7, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
weifit = wblfit(exprd);
qan7 = 1-wblcdf(100,weifit(1),weifit(2));




% 8, lognormal distribution
% m     = expar1;
% v     = stdpar1;
% lmu   = log(m^2/sqrt(v+m^2));
% lsigma= sqrt(log(v/m^2 + 1)); 
% [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
%     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
logparhat = lognfit(exprd);
qan8 = 1-logncdf(100,logparhat(1),logparhat(2));
% lnorm
% qan8

% try the exponential distribution with multiple ratio of time and 
%  exponential parameter
% 150/120 = 1.25
% try ratio = 0.5, 0.75,  1, 1.5, 1.75, 2
%    time   = 60,   90, 120, 180, 210,  240

% 9, exponential using data regarding whatever larger than 60 as censored
expar9 = (sum(exprd(find(exprd<=60)))+length(find(exprd>60))*60)...
    /length(find(exprd<=60));
qan9 = exp(-100/expar9);


% 10, exponential using data regarding whatever larger than 90 as censored
expar10 = (sum(exprd(find(exprd<=90)))+length(find(exprd>90))*90)...
    /length(find(exprd<=90));
qan10 = exp(-100/expar10);


% 11, exponential using data regarding whatever larger than 120 as censored
expar11 = (sum(exprd(find(exprd<=120)))+length(find(exprd>120))*120)...
    /length(find(exprd<=120));
qan11 = exp(-100/expar11);


% 12, exponential using data regarding whatever larger than 180 as censored
expar12 = (sum(exprd(find(exprd<=180)))+length(find(exprd>180))*180)...
    /length(find(exprd<=180));
qan12 = exp(-100/expar12);


% 13, exponential using data regarding whatever larger than 210 as censored
expar13 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
    /length(find(exprd<=210));
qan13 = exp(-100/expar13);


% 14, pexe estimate
[time_die14,ttot14,deaths14] = totaltest(exprd,ones(Ns,1));
lendis = length(ttot14);
[qan14 lamest14] = pexeest(exprd, ones(Ns,1), time_die14(1:lendis-1), 100);


% Estimate the change point from mono and umbrella assumptions
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
pexeoutmono  = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
pexeoutum    = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',4);
pexeoutum2    = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',5);

% 15, rpexe estimate monotone, p-value = 0.05
if min(pexeoutmono.pvalues) < crivalmono5, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono5));
    [qan15 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
    %rpexrec(iter) = min(pexeoutmono.pvalues);
else 
    expar15 = expar1;
    qan15   = qan3;
end
% 16, rpexe estimate montone, p-value = 0.1
if min(pexeoutmono.pvalues) < crivalmono10, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono10));
    [qan16 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar16 = expar1;
    qan16   = qan3;
end
% 16_5, rpexe estimate montone, p-value = 0.1
if min(pexeoutmono.pvalues) < crivalmono20, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono20));
    [qan16_5 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar16_5 = expar1;
    qan16_5   = qan3;
end

% 17, rpexe estimate umbrella, down and up, p-value = 0.05
if min(pexeoutum.pvalues) < crivalum5, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum5));
    [qan17 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar17 = expar1;
    qan17   = qan3;
end
% 18, rpexe estimate umbrella, down and up, p-value = 0.1
if min(pexeoutum.pvalues) < crivalum10, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum10));
    [qan18 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar18 = expar1;
    qan18   = qan3;
end
% 18_5, rpexe estimate umbrella, down and up, p-value = 0.1
if min(pexeoutum.pvalues) < crivalum20, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum20));
    [qan18_5 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar18_5 = expar1;
    qan18_5   = qan3;
end

% 19, rpexe estimate umbrella, up and down, p-value = 0.05
if min(pexeoutum2.pvalues) < crivalum5, % if the change point is found
    tchange = pexeoutum2.times(find(pexeoutum2.pvalues < crivalum5));
    [qan19 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar19 = expar1;
    qan19   = qan3;
end
% 20, rpexe estimate umbrella, up and down, p-value = 0.1
if min(pexeoutum2.pvalues) < crivalum10, % if the change point is found
    tchange = pexeoutum2.times(find(pexeoutum2.pvalues < crivalum10));
    [qan20 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar20 = expar1;
    qan20   = qan3;
end
% 20_5, rpexe estimate umbrella, up and down, p-value = 0.1
if min(pexeoutum2.pvalues) < crivalum10, % if the change point is found
    tchange = pexeoutum2.times(find(pexeoutum2.pvalues < crivalum20));
    [qan20_5 lamest] = pexeest(exprd, ones(Ns,1), tchange, 100);
else 
    expar20_5 = expar1;
    qan20_5   = qan3;
end


kme1(iter) = qan1;
kme2(iter) = qan2;
exp1(iter) = qan3;
exp2(iter) = qan4;
kme3(iter) = qan5;
gam1(iter) = qan6;
wei1(iter) = qan7;
ln1(iter)  = qan8;
exp9(iter) = qan9;
exp10(iter) = qan10;
exp11(iter) = qan11;
exp12(iter) = qan12;
exp13(iter) = qan13;
exp14(iter) = qan14;
exp15(iter) = qan15;
exp16(iter) = qan16;
exp16_5(iter) = qan16_5;
exp17(iter) = qan17;
exp18(iter) = qan18;
exp18_5(iter) = qan18_5;
exp19(iter) = qan19;
exp20(iter) = qan20;
exp20_5(iter) = qan20_5;
end


% 
quanval = 1-gamcdf(lambda^(-2)*(exp(-beta)*100)^(lambda/sigma),lambda^(-2),1);

        
% compute the roon mean squared prediction error
rmspe1 = sqrt(sum((kme1-quanval)'*(kme1-quanval)/N));
rmspe2 = sqrt(sum((kme2-quanval)'*(kme2-quanval)/N));
rmspe3 = sqrt(sum((exp1-quanval)'*(exp1-quanval)/N));
rmspe4 = sqrt(sum((exp2-quanval)'*(exp2-quanval)/N));
rmspe5 = sqrt(sum((kme3-quanval)'*(kme3-quanval)/N));
rmspe6 = sqrt(sum((gam1-quanval)'*(gam1-quanval)/N));
rmspe7 = sqrt(sum((wei1-quanval)'*(wei1-quanval)/N));
rmspe8 = sqrt(sum((ln1-quanval)'*(ln1-quanval)/N));
rmspe9 = sqrt(sum((exp9-quanval)'*(exp9-quanval)/N));
rmspe10 = sqrt(sum((exp10-quanval)'*(exp10-quanval)/N));
rmspe11 = sqrt(sum((exp11-quanval)'*(exp11-quanval)/N));
rmspe12 = sqrt(sum((exp12-quanval)'*(exp12-quanval)/N));
rmspe13 = sqrt(sum((exp13-quanval)'*(exp13-quanval)/N));
rmspe14 = sqrt(sum((exp14-quanval)'*(exp14-quanval)/N));
rmspe15 = sqrt(sum((exp15-quanval)'*(exp15-quanval)/N));
rmspe16 = sqrt(sum((exp16-quanval)'*(exp16-quanval)/N));
rmspe16_5 = sqrt(sum((exp16_5-quanval)'*(exp16_5-quanval)/N));
rmspe17 = sqrt(sum((exp17-quanval)'*(exp17-quanval)/N));
rmspe18 = sqrt(sum((exp18-quanval)'*(exp18-quanval)/N));
rmspe18_5 = sqrt(sum((exp18_5-quanval)'*(exp18_5-quanval)/N));
rmspe19 = sqrt(sum((exp19-quanval)'*(exp19-quanval)/N));
rmspe20 = sqrt(sum((exp20-quanval)'*(exp20-quanval)/N));
rmspe20_5 = sqrt(sum((exp20_5-quanval)'*(exp20_5-quanval)/N));

res_gam = [[var(kme1) mean(kme1)-quanval rmspe1]
[var(kme3) mean(kme3)-quanval rmspe5]
[var(kme2) mean(kme2)-quanval rmspe2]
[var(gam1) mean(gam1)-quanval rmspe6]
[var(wei1) mean(wei1)-quanval rmspe7]
[var(ln1) mean(ln1)-quanval rmspe8]
[var(exp1) mean(exp1)-quanval rmspe3]
[var(exp13) mean(exp13)-quanval rmspe13]
[var(exp12) mean(exp12)-quanval rmspe12]
[var(exp2) mean(exp2)-quanval rmspe4]
[var(exp11) mean(exp11)-quanval rmspe11]
% [var(exp9) mean(exp9)-quanval rmspe9]
% [var(exp10) mean(exp10)-quanval rmspe10]
[var(exp14) mean(exp14)-quanval rmspe14]
[var(exp15) mean(exp15)-quanval rmspe15]
[var(exp16) mean(exp16)-quanval rmspe16]
[var(exp16_5) mean(exp16_5)-quanval rmspe16_5]
[var(exp17) mean(exp17)-quanval rmspe17]
[var(exp18) mean(exp18)-quanval rmspe18]
[var(exp18_5) mean(exp18_5)-quanval rmspe18_5]
[var(exp19) mean(exp19)-quanval rmspe19]
[var(exp20) mean(exp20)-quanval rmspe20]
[var(exp20_5) mean(exp20_5)-quanval rmspe20_5]];

