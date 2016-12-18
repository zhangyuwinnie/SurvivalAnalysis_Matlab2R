function res_gam = gen_gam_test_new(beta,lambda,sigma,Ns,ts,Nsim)
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



% set the number of iteration samples
N = Nsim;

kme_raw     = zeros(N,1);
exp_raw     = zeros(N,1);
pexe        = zeros(N,1);
        
umrpexe100  = zeros(N,1);
umrpexe70   = zeros(N,1);
umrpexe40   = zeros(N,1);
umrpexe20   = zeros(N,1);
umrpexe10   = zeros(N,1);
umrpexe5    = zeros(N,1);

monorpexe100= zeros(N,1);
monorpexe70 = zeros(N,1);
monorpexe40 = zeros(N,1);
monorpexe20 = zeros(N,1);
monorpexe10 = zeros(N,1);
monorpexe5  = zeros(N,1);

Ns = 50;

crivalmono5  = exp(-4.2331-0.3938*log(Ns));
crivalmono10 = exp(-3.4481-0.3847*log(Ns));
% approximate 20% using the quantiles of beta distribution
% par = fmincon(@(par) betafit(par,crivalmono5,crivalmono10),[1 10], ...
%     [],[],[],[],[0.001 0.001],[]);
parmono = [0.8639   9.6603];
crivalmono20 = betainv(0.2,parmono(1), parmono(2)); 
crivalmono40 = betainv(0.4,parmono(1), parmono(2)); 
crivalmono70 = betainv(0.7,parmono(1), parmono(2)); 
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
crivalum20 = betainv(0.2,parum(1), parum(2)); 
crivalum40 = betainv(0.4,parum(1), parum(2)); 
crivalum70 = betainv(0.7,parum(1), parum(2)); 
crivalum100 = 1;


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


% 2, exponential using all the data
expar1 = mean(exprd);
stdpar1= std(exprd);
qexp_raw = exp(-ts/expar1);


% 3, pexe estimate
[time_die,ttot,deaths] = totaltest(exprd,ones(Ns,1));
lendis = length(ttot);
[qpexe lamest] = pexeest(exprd, ones(Ns,1), time_die(1:lendis-1), ts);


% Estimate the change point from mono and umbrella assumptions
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
pexeoutmono  = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',3);
pexeoutum    = RPEXEv1('EventTime',exprd,'Censor',ones(Ns,1),'Trend',4);


% 4, rpexe estimate umbrella, down and up, p-value = 1
if min(pexeoutum.pvalues) < crivalum100, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum100));
    [qumrpexe100 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe100 = qexp_raw;
end


% 5, rpexe estimate umbrella, down and up, p-value = 0.7
if min(pexeoutum.pvalues) < crivalum70, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum70));
    [qumrpexe70 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe70 = qexp_raw;
end


% 6, rpexe estimate umbrella, down and up, p-value = 0.4
if min(pexeoutum.pvalues) < crivalum40, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum40));
    [qumrpexe40 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe40 = qexp_raw;
end


% 7, rpexe estimate umbrella, down and up, p-value = 0.2
if min(pexeoutum.pvalues) < crivalum20, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum20));
    [qumrpexe20 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe20 = qexp_raw;
end


% 8, rpexe estimate umbrella, down and up, p-value = 0.1
if min(pexeoutum.pvalues) < crivalum10, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum10));
    [qumrpexe10 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe10 = qexp_raw;
end


% 9, rpexe estimate umbrella, down and up, p-value = 0.05
if min(pexeoutum.pvalues) < crivalum5, % if the change point is found
    tchange = pexeoutum.times(find(pexeoutum.pvalues < crivalum5));
    [qumrpexe5 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qumrpexe5 = qexp_raw;
end


% 10, rpexe estimate montone, p-value = 1
if min(pexeoutmono.pvalues) < crivalmono100, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono100));
    [qmonorpexe100 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe100 = qexp_raw;
end


% 11, rpexe estimate montone, p-value = 0.7
if min(pexeoutmono.pvalues) < crivalmono70, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono70));
    [qmonorpexe70 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe70 = qexp_raw;
end


% 12, rpexe estimate montone, p-value = 0.4
if min(pexeoutmono.pvalues) < crivalmono40, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono40));
    [qmonorpexe40 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe40 = qexp_raw;
end


% 13, rpexe estimate montone, p-value = 0.2
if min(pexeoutmono.pvalues) < crivalmono20, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono20));
    [qmonorpexe20 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe20 = qexp_raw;
end


% 14, rpexe estimate montone, p-value = 0.1
if min(pexeoutmono.pvalues) < crivalmono10, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono10));
    [qmonorpexe10 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe10 = qexp_raw;
end


% 15, rpexe estimate monotone, p-value = 0.05
if min(pexeoutmono.pvalues) < crivalmono5, % if the change point is found
    tchange = pexeoutmono.times(find(pexeoutmono.pvalues < crivalmono5));
    [qmonorpexe5 lamest] = pexeest(exprd, ones(Ns,1), tchange, ts);
else 
    qmonorpexe5 = qexp_raw;   
end

kme_raw(iter) = qkme_raw;
exp_raw(iter) = qexp_raw;
pexe(iter)    = qpexe;
umrpexe100(iter) = qumrpexe100;
umrpexe70(iter) = qumrpexe70;
umrpexe40(iter) = qumrpexe40;
umrpexe20(iter) = qumrpexe20;
umrpexe10(iter) = qumrpexe10;
umrpexe5(iter)  = qumrpexe5;
monorpexe100(iter) = qmonorpexe100;
monorpexe70(iter) = qmonorpexe70;
monorpexe40(iter)  = qmonorpexe40;
monorpexe20(iter) = qmonorpexe20;
monorpexe10(iter) = qmonorpexe10;
monorpexe5(iter)  = qmonorpexe5;
end

% 
if lambda >= 0, 
    quanval = 1-gamcdf(lambda^(-2)*(exp(-beta)*ts)^(lambda/sigma),lambda^(-2),1);
else
    quanval = gamcdf(lambda^(-2)*(exp(-beta)*ts)^(lambda/sigma),lambda^(-2),1);
end
% compute the roon mean squared prediction error

rmspekme_raw = sqrt(sum((kme_raw-quanval)'*(kme_raw-quanval)/N));
rmspeexp_raw = sqrt(sum((exp_raw-quanval)'*(exp_raw-quanval)/N));
rmspepexe = sqrt(sum((pexe-quanval)'*(pexe-quanval)/N));
rmspeumrpexe100= sqrt(sum((umrpexe100-quanval)'*(umrpexe100-quanval)/N));
rmspeumrpexe70 = sqrt(sum((umrpexe70-quanval)'*(umrpexe70-quanval)/N));
rmspeumrpexe40 = sqrt(sum((umrpexe40-quanval)'*(umrpexe40-quanval)/N));
rmspeumrpexe20 = sqrt(sum((umrpexe20-quanval)'*(umrpexe20-quanval)/N));
rmspeumrpexe10 = sqrt(sum((umrpexe10-quanval)'*(umrpexe10-quanval)/N));
rmspeumrpexe5  = sqrt(sum((umrpexe5-quanval)'*(umrpexe5-quanval)/N));
rmspemonorpexe100= sqrt(sum((monorpexe100-quanval)'*(monorpexe100-quanval)/N));
rmspemonorpexe70 = sqrt(sum((monorpexe70-quanval)'*(monorpexe70-quanval)/N));
rmspemonorpexe40 = sqrt(sum((monorpexe40-quanval)'*(monorpexe40-quanval)/N));
rmspemonorpexe20 = sqrt(sum((monorpexe20-quanval)'*(monorpexe20-quanval)/N));
rmspemonorpexe10 = sqrt(sum((monorpexe10-quanval)'*(monorpexe10-quanval)/N));
rmspemonorpexe5  = sqrt(sum((monorpexe5-quanval)'*(monorpexe5-quanval)/N));


res_gam = [[var(kme_raw) mean(kme_raw)-quanval rmspekme_raw max(abs(kme_raw-quanval)) mean(abs(kme_raw-quanval))]
[var(exp_raw) mean(exp_raw)-quanval rmspeexp_raw max(abs(exp_raw-quanval)) mean(abs(exp_raw-quanval))]
[var(pexe) mean(pexe)-quanval rmspepexe max(abs(pexe-quanval)) mean(abs(pexe-quanval))]
[var(umrpexe100) mean(umrpexe100)-quanval rmspeumrpexe100 max(abs(umrpexe100-quanval)) mean(abs(umrpexe100-quanval))]
[var(umrpexe70) mean(umrpexe70)-quanval rmspeumrpexe70 max(abs(umrpexe70-quanval)) mean(abs(umrpexe70-quanval))]
[var(umrpexe40) mean(umrpexe40)-quanval rmspeumrpexe40 max(abs(umrpexe40-quanval)) mean(abs(umrpexe40-quanval))]
[var(umrpexe20) mean(umrpexe20)-quanval rmspeumrpexe20 max(abs(umrpexe20-quanval)) mean(abs(umrpexe20-quanval))]
[var(umrpexe10) mean(umrpexe10)-quanval rmspeumrpexe10 max(abs(umrpexe10-quanval)) mean(abs(umrpexe10-quanval))]
[var(umrpexe5) mean(umrpexe5)-quanval rmspeumrpexe5 max(abs(umrpexe5-quanval)) mean(abs(umrpexe5-quanval))]
[var(monorpexe100) mean(monorpexe100)-quanval rmspemonorpexe100 max(abs(monorpexe100-quanval)) mean(abs(monorpexe100-quanval))]
[var(monorpexe70) mean(monorpexe70)-quanval rmspemonorpexe70 max(abs(monorpexe70-quanval)) mean(abs(monorpexe70-quanval))]
[var(monorpexe40) mean(monorpexe40)-quanval rmspemonorpexe40 max(abs(monorpexe40-quanval)) mean(abs(monorpexe40-quanval))]
[var(monorpexe20) mean(monorpexe20)-quanval rmspemonorpexe20 max(abs(monorpexe20-quanval)) mean(abs(monorpexe20-quanval))]
[var(monorpexe10) mean(monorpexe10)-quanval rmspemonorpexe10 max(abs(monorpexe10-quanval)) mean(abs(monorpexe10-quanval))]
[var(monorpexe5) mean(monorpexe5)-quanval rmspemonorpexe5 max(abs(monorpexe5-quanval)) mean(abs(monorpexe5-quanval))]];



