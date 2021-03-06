%%% Compare the Monte Carlo variances of four predictors in 
%%% estimating the quantile in the day 100.


% 1, Kaplan-Meier binomial distribution
% 2, Kaplan-Meier averaging some observations
% 3, Exponential all the way
% 4, Exponential up to a point


kme1 = zeros(5000,1);
kme2 = zeros(5000,1);
exp1 = zeros(5000,1);
exp2 = zeros(5000,1);

Ns   = 15;

for iter = 1:5000,
    
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
qan3 = exp(-100/expar1);
% 4, exponential using data regarding whatever larger than 150 as censored
expar2 = (sum(exprnd(find(exprnd<=150)))+length(find(exprnd>150))*150)...
    /length(find(exprnd<=150));
qan4 = exp(-100/expar2);

kme1(iter) = qan1;
kme2(iter) = qan2;
exp1(iter) = qan3;
exp2(iter) = qan4;

end;

% true quantile = exp(-100/120)=0.4346

subplot(2,2,1);
hist(kme1);
title('Kaplan-Meier Estimate');
xlim([0.1,0.9]);
subplot(2,2,2);
hist(kme2);
title('Kaplan-Meier Estimate by averaging over window of 30');
xlim([0.1,0.9]);
subplot(2,2,3);
hist(exp1);
title('Exponential Estimate');
xlim([0.1,0.9]);
subplot(2,2,4);
hist(exp2);
title('Exponential Estimate with censoring at 150');
xlim([0.1,0.9]);

% compute the roon mean squared prediction error
rmspe1 = sqrt(sum((kme1-0.4346)'*(kme1-0.4346)/5000))
rmspe2 = sqrt(sum((kme2-0.4346)'*(kme2-0.4346)/5000))
rmspe3 = sqrt(sum((exp1-0.4346)'*(exp1-0.4346)/5000))
rmspe4 = sqrt(sum((exp2-0.4346)'*(exp2-0.4346)/5000))


[var(kme1) mean(kme1)-0.4346]
[var(kme2) mean(kme2)-0.4346]
[var(exp1) mean(exp1)-0.4346]
[var(exp2) mean(exp2)-0.4346]















