function res_gam = gen_gam_test_new2(beta,lambda,sigma,Ns,ts,Nsim)
% A function for estimating the predictive accuracy of various predictions
% for ganeralized gamma distribution.
% Input variables
%   beta, lambda, and sigma: three parameters in the generalized gamma
%                           distribution. 
%   Ns: sample size at each iteration
%   ts: time for estimating the quantile
%   Nsim: number of simulation
% Output variables
%   res_gam: bias, variance, and the rmspe of the estimators.

% Reference paper
%
% Parametric survival analysis and taxonomy of hazard functions
% for the generalized gamma distribution
% Christopher Cox, Haitao Chu, Michael F. Schneider and Alvaro Munoz,
% STATISTICS IN MEDICINE
% Statist. Med. 2007; 26:4352–4374


crivalmono5  = exp(-4.2331-0.3938*log(Ns));
crivalmono10 = exp(-3.4481-0.3847*log(Ns));
% approximate 20% using the quantiles of beta distribution
% par = fmincon(@(par) betafit(par,crivalmono5,crivalmono10),[1 10], ...
%     [],[],[],[],[0.001 0.001],[]);
parmono = [0.8639   9.6603];
crivalmono25 = betainv(0.25,parmono(1), parmono(2)); 
crivalmono50 = betainv(0.50,parmono(1), parmono(2)); 
crivalmono75 = betainv(0.75,parmono(1), parmono(2)); 
crivalmono100 = 1;


crivalum5  = exp(-1.2234-2.7044*log(Ns));
crivalum10 = exp(-1.1951-2.0238*log(Ns));
% approximate 20% using the quantiles of beta distribution
% options = optimset('Display','off','Simplex','on',...
%                    'TolFun',.000001, 'Algorithm','active-set');
% par = fmincon(@(par) betafit(par,crivalum5,crivalum10),[0.2 1.43], ...
%     [],[],[],[],[0.001 0.001],[],[],options);
% betacdf(crivalum5,par(1),par(2))
% betacdf(crivalum10,par(1),par(2))
parum = [0.2856 2.9087];
crivalum25 = betainv(0.25,parum(1), parum(2)); 
crivalum50 = betainv(0.50,parum(1), parum(2)); 
crivalum75 = betainv(0.75,parum(1), parum(2)); 
crivalum100 = 1;



% set the number of iteration samples
N = Nsim;

kme_raw = zeros(N,1);
kme_10  = zeros(N,1);
kme_30  = zeros(N,1);
exp_raw = zeros(N,1);
exp_210 = zeros(N,1);
exp_150 = zeros(N,1);
gam     = zeros(N,1);
wei     = zeros(N,1);
logn    = zeros(N,1);
pexe    = zeros(N,1);
umrpexe40 = zeros(N,1);
umrpexe10 = zeros(N,1);
umrpexe5  = zeros(N,1);
monorpexe40 = zeros(N,1);
monorpexe10 = zeros(N,1);
monorpexe5  = zeros(N,1);


for iter = 1:N,
    
% Step 1, generate data from generalized gamma distribution, with sample
% size Ns
exprd = sort(gamrnd(1/(lambda^2),1,Ns,1));
for k = 1:Ns,
    exprd(k) = exp(beta)*(exprd(k)*lambda^2)^(sigma/lambda);
end

% Step 2, estimate the quantile at ts days using four approaches

% 1, KME
qkme_raw   = length(find(exprd>ts))/length(exprd);

% 2, KME with 10 days window from days ts-5 to ts+5.
tmpq2  = intersect(find(exprd>(ts-5)),find(exprd<(ts+5))); 
np2    = length(tmpq2);
if np2 < 2, 
    qkme_10   = qkme_raw;
else
    qkme_10   = 0;
        % calculate the area
    for i = 1:(np2-1),
        qkme_10 = qkme_10 + (exprd(tmpq2(i+1))-exprd(tmpq2(i)))*...
            length(find(exprd>exprd(tmpq2(i+1))))/length(exprd);
    end
        % calculate the weighted average of the quantile
    qkme_10 = qkme_10/(exprd(tmpq2(np2))-exprd(tmpq2(1)));
end


% 3, KME with averaging ts-15 -- ts+15 days window
tmpq   = intersect(find(exprd>ts-15),find(exprd<ts+15)); 
np     = length(tmpq);
if np <2, 
    qkme_30   = qkme_raw;
else
    qkme_30   = 0;
        % calculate the area
    for i = 1:(np-1),
        qkme_30 = qkme_30 + (exprd(tmpq(i+1))-exprd(tmpq(i)))*...
            length(find(exprd>exprd(tmpq(i+1))))/length(exprd);
    end
        % calculate the weighted average of the quantile
    qkme_30 = qkme_30/(exprd(tmpq(np))-exprd(tmpq(1)));
end


% 4, exponential using all the data
expar1 = mean(exprd);
stdpar1= std(exprd);
qexp_raw = exp(-ts/expar1);


% 5, exponential using data regarding whatever larger than 210 as censored
expar2 = (sum(exprd(find(exprd<=210)))+length(find(exprd>210))*210)...
    /length(find(exprd<=210));
qexp_210 = exp(-ts/expar2);


% 6, exponential using data regarding whatever larger than 150 as censored
expar3 = (sum(exprd(find(exprd<=150)))+length(find(exprd>150))*150)...
    /length(find(exprd<=150));
qexp_150 = exp(-ts/expar3);


% 7, gamma distribution
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
qgam  = 1-gamcdf(ts,gmfit(1),gmfit(2));




% 8, weibull distribution
% wa0     = expar1;
% wb0     = 1;
% [wei fval]   = fmincon(@(wei) wbllike([wei(1);wei(2)],exprd),...
%     [wa0; wb0],[],[],[],[],[0.001;0.001],[1000;1000],[],options); 
weifit = wblfit(exprd);
qwei   = 1-wblcdf(ts,weifit(1),weifit(2));


% 9, lognormal distribution
% m     = expar1;
% v     = stdpar1;
% lmu   = log(m^2/sqrt(v+m^2));
% lsigma= sqrt(log(v/m^2 + 1)); 
% [lnorm fval] = fmincon(@(lnorm) lognlike([lnorm(1);lnorm(2)],exprd),...
%     [m;v],[],[],[],[],[-1000;0.0001],[1000;1000],[],options); 
logparhat = lognfit(exprd);
qlogn  = 1-logncdf(ts,logparhat(1),logparhat(2));


% 10, pexe estimate
[time_die,ttot,deaths] = totaltest(exprd,ones(Ns,1));
lendis = length(ttot);
[qpexe lamest] = pexeest(exprd, ones(Ns,1), time_die(1:lendis-1), ts);




% Estimate the change point from mono and umbrella assumptions
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
pexeoutmono  = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
pexeoutum    = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',4);


% 11, rpexe estimate umbrella, down and up, p-value = 0.1
if min(pexeoutum.pvalues) < crivalum10, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum10));
    [qumrpexe10 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe10 = qexp_raw;
end

% 12, rpexe estimate umbrella, down and up, p-value = 0.25
if min(pexeoutum.pvalues) < crivalum25, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum25));
    [qumrpexe25 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe25 = qexp_raw;
end

% 13, rpexe estimate umbrella, down and up, p-value = 0.5
if min(pexeoutum.pvalues) < crivalum50, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum50));
    [qumrpexe50 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe50 = qexp_raw;
end

% 14, rpexe estimate umbrella, down and up, p-value = 0.75
if min(pexeoutum.pvalues) < crivalum75, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum75));
    [qumrpexe75 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe75 = qexp_raw;
end

% 15, rpexe estimate umbrella, down and up, p-value = 1
if min(pexeoutum.pvalues) < crivalum100, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum100));
    [qumrpexe100 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe100 = qexp_raw;
end

% 16, rpexe estimate montone, p-value = 0.1
if min(pexeoutmono.pvalues) < crivalmono10, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono10));
    [qmonorpexe10 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe10 = qexp_raw;
end

% 17, rpexe estimate monotone, p-value = 0.25
if min(pexeoutmono.pvalues) < crivalmono25, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono25));
    [qmonorpexe25 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe25 = qexp_raw;   
end

% 18, rpexe estimate monotone, p-value = 0.5
if min(pexeoutmono.pvalues) < crivalmono50, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono50));
    [qmonorpexe50 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe50 = qexp_raw;   
end

% 19, rpexe estimate monotone, p-value = 0.75
if min(pexeoutmono.pvalues) < crivalmono75, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono75));
    [qmonorpexe75 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe75 = qexp_raw;   
end


% 20, rpexe estimate monotone, p-value = 1.00
if min(pexeoutmono.pvalues) < crivalmono100, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono100));
    [qmonorpexe100 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe100 = qexp_raw;   
end

kme_raw(iter) = qkme_raw;
kme_10(iter)  = qkme_10;
kme_30(iter)  = qkme_30;
exp_raw(iter) = qexp_raw;
exp_210(iter) = qexp_210;
exp_150(iter) = qexp_150;
gam(iter)     = qgam;
wei(iter)     = qwei;
logn(iter)    = qlogn;
pexe(iter)    = qpexe;
umrpexe10(iter) = qumrpexe10;
umrpexe25(iter) = qumrpexe25;
umrpexe50(iter) = qumrpexe50;
umrpexe75(iter) = qumrpexe75;
umrpexe100(iter) = qumrpexe100;
monorpexe10(iter) = qmonorpexe10;
monorpexe25(iter) = qmonorpexe25;
monorpexe50(iter) = qmonorpexe50;
monorpexe75(iter) = qmonorpexe75;
monorpexe100(iter) = qmonorpexe100;
end

% 
if lambda >= 0, 
    quanval = 1-gamcdf(lambda^(-2)*(exp(-beta)*ts)^(lambda/sigma),lambda^(-2),1);
else
    quanval = gamcdf(lambda^(-2)*(exp(-beta)*ts)^(lambda/sigma),lambda^(-2),1);
end
% compute the roon mean squared prediction error

rmspekme_raw = sqrt(sum((kme_raw-quanval)'*(kme_raw-quanval)/N));
rmspekme_10 = sqrt(sum((kme_10-quanval)'*(kme_10-quanval)/N));
rmspekme_30 = sqrt(sum((kme_30-quanval)'*(kme_30-quanval)/N));
rmspeexp_raw = sqrt(sum((exp_raw-quanval)'*(exp_raw-quanval)/N));
rmspeexp_210 = sqrt(sum((exp_210-quanval)'*(exp_210-quanval)/N));
rmspeexp_150 = sqrt(sum((exp_150-quanval)'*(exp_150-quanval)/N));
rmspegam = sqrt(sum((gam-quanval)'*(gam-quanval)/N));
rmspewei = sqrt(sum((wei-quanval)'*(wei-quanval)/N));
rmspelogn = sqrt(sum((logn-quanval)'*(logn-quanval)/N));
rmspepexe = sqrt(sum((pexe-quanval)'*(pexe-quanval)/N));
rmspeumrpexe10 = sqrt(sum((umrpexe10-quanval)'*(umrpexe10-quanval)/N));
rmspeumrpexe25 = sqrt(sum((umrpexe25-quanval)'*(umrpexe25-quanval)/N));
rmspeumrpexe50 = sqrt(sum((umrpexe50-quanval)'*(umrpexe50-quanval)/N));
rmspeumrpexe75 = sqrt(sum((umrpexe75-quanval)'*(umrpexe75-quanval)/N));
rmspeumrpexe100 = sqrt(sum((umrpexe100-quanval)'*(umrpexe100-quanval)/N));
rmspemonorpexe10 = sqrt(sum((monorpexe10-quanval)'*(monorpexe10-quanval)/N));
rmspemonorpexe25 = sqrt(sum((monorpexe25-quanval)'*(monorpexe25-quanval)/N));
rmspemonorpexe50 = sqrt(sum((monorpexe50-quanval)'*(monorpexe50-quanval)/N));
rmspemonorpexe75 = sqrt(sum((monorpexe75-quanval)'*(monorpexe75-quanval)/N));
rmspemonorpexe50 = sqrt(sum((monorpexe50-quanval)'*(monorpexe50-quanval)/N));

res_gam = [[var(kme_raw) mean(kme_raw)-quanval rmspekme_raw max(abs(kme_raw-quanval)) mean(abs(kme_raw-quanval))]
[var(kme_10) mean(kme_10)-quanval rmspekme_10 max(abs(kme_10-quanval)) mean(abs(kme_10-quanval))]
[var(kme_30) mean(kme_30)-quanval rmspekme_30 max(abs(kme_30-quanval)) mean(abs(kme_30-quanval))]
[var(exp_raw) mean(exp_raw)-quanval rmspeexp_raw max(abs(exp_raw-quanval)) mean(abs(exp_raw-quanval))]
[var(exp_210) mean(exp_210)-quanval rmspeexp_210 max(abs(exp_210-quanval)) mean(abs(exp_210-quanval))]
[var(exp_150) mean(exp_150)-quanval rmspeexp_150 max(abs(exp_150-quanval)) mean(abs(exp_150-quanval))]
[var(gam) mean(gam)-quanval rmspegam max(abs(gam-quanval)) mean(abs(gam-quanval))]
[var(wei) mean(wei)-quanval rmspewei max(abs(wei-quanval)) mean(abs(wei-quanval))]
[var(logn) mean(logn)-quanval rmspelogn max(abs(logn-quanval)) mean(abs(logn-quanval))]
[var(pexe) mean(pexe)-quanval rmspepexe max(abs(pexe-quanval)) mean(abs(pexe-quanval))]
[var(umrpexe10) mean(umrpexe10)-quanval rmspeumrpexe10 max(abs(umrpexe10-quanval)) mean(abs(umrpexe10-quanval))]
[var(umrpexe25) mean(umrpexe25)-quanval rmspeumrpexe25 max(abs(umrpexe25-quanval)) mean(abs(umrpexe25-quanval))]
[var(umrpexe50) mean(umrpexe50)-quanval rmspeumrpexe50 max(abs(umrpexe50-quanval)) mean(abs(umrpexe50-quanval))]
[var(umrpexe75) mean(umrpexe75)-quanval rmspeumrpexe75 max(abs(umrpexe75-quanval)) mean(abs(umrpexe75-quanval))]
[var(umrpexe50) mean(umrpexe100)-quanval rmspeumrpexe100 max(abs(umrpexe100-quanval)) mean(abs(umrpexe100-quanval))]
[var(monorpexe10) mean(monorpexe10)-quanval rmspemonorpexe10 max(abs(monorpexe10-quanval)) mean(abs(monorpexe10-quanval))]
[var(monorpexe25) mean(monorpexe25)-quanval rmspemonorpexe25 max(abs(monorpexe25-quanval)) mean(abs(monorpexe25-quanval))]
[var(monorpexe50) mean(monorpexe50)-quanval rmspemonorpexe50 max(abs(monorpexe50-quanval)) mean(abs(monorpexe50-quanval))]
[var(monorpexe75) mean(monorpexe75)-quanval rmspemonorpexe75 max(abs(monorpexe75-quanval)) mean(abs(monorpexe75-quanval))]
[var(monorpexe100) mean(monorpexe100)-quanval rmspemonorpexe100 max(abs(monorpexe100-quanval)) mean(abs(monorpexe100-quanval))]];



