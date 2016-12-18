function [xpart,ypart] = km_csolid(time, censor, plotcens)
% Function km plots the kaplan meier curve 
% Given times and censor (0 = censored; 1 = uncensored)
% returns x's and y's on kaplan meier curve

% compute realt and deaths
tmptime = time;
for i=1:length(censor),
    if censor(i)==0,
        tmptime(i) = 0;
    end
end
timesort = sort(tmptime);
%timesort(1)=[];
realt = unique(timesort);

%%% input %%%
%%% realt:  sorted times when death occur,
%%% deaths: corresponding number of deaths at the sorted times
%%% time:   all the times, including censor and death,
%%% censor: a vector indicating censored or not at the given times in time
%%%         censored: censor() = 0; noncensored: censor() =1;

%%% output %%
%%% the label on the kaplan-mejer curve, which is corresponding to each
%%% point on realt and the plot of the kaplan-meier curve.


ldea    = length(realt);
pos_km  = zeros(ldea,1);
% at each time point, the key is to compute ni and di,
% the kaplan meier estimate has the form
%       est_s(t) = prod_{ti<t}((ni-di)/ni);
%
% get the values of the first item in pos_km;

% compute two sequences ni and di
niseq   = zeros(ldea,1);
diseq   = zeros(ldea,1);
difnidi = zeros(ldea,1);
for i   = 1:ldea,
    numberp = ones(length(censor),1);
    for ii = 1:length(censor),
        if time(ii) <= realt(i),
            numberp(ii) = 0;
            if censor(ii) == 1,
                diseq(i) = diseq(i)+1; % this is the cumulative deaths
            end
        end
        niseq(i) = sum(numberp);
    end
end

% substract the cumulative death
for j = 2:ldea,
    diseq(ldea+2-j) = diseq(ldea+2-j) - diseq(ldea+2-j-1);
end

difnidi = niseq - diseq;
elekm   = zeros(ldea,1);

for i = 1:ldea,
    elekm(i) = difnidi(i)/niseq(i);
    pos_km(i) = prod(elekm(1:i));
end

% given pos_km, plot kaplan-meier curve
upperx              = zeros(ldea+1,1);
lowerx              = zeros(ldea+1,1);
upperx(1)           = 0;
upperx(2:(ldea+1))  = realt;
lowerx(1:ldea)      = realt;
lowerx(1+ldea)      = max(time);
uppery              = zeros(ldea+1,1);
lowery              = zeros(ldea+1,1);
uppery(1:2)         = 1;
uppery(3:ldea+1)      = pos_km(1:ldea-1);
lowery(1:ldea)      = pos_km;
lowery(1+ldea)      = pos_km(ldea);

xpart               = zeros(2*ldea+2,1);
ypart               = zeros(2*ldea+2,1);
xpartl              = lowerx(1:ldea);
xpartu              = upperx(2:ldea+1);
ypartl              = lowery(1:ldea);
ypartu              = uppery(2:ldea+1);

xpart = zeros(2*ldea+2,1);
ypart = zeros(2*ldea+2,1);
% save the data in xpart and ypart
xpart(1) = upperx(1);
xpart(2) = upperx(2);
for i = 2:ldea,
    xpart(2*(i-1)+1)= lowerx(i-1);
    xpart(2*(i-1)+2)= upperx(i+1);
end
xpart(2*ldea+1) = lowerx(ldea);
xpart(2*ldea+2) = lowerx(ldea+1);

ypart(1) = uppery(1);
ypart(2) = uppery(2);
for i = 2:ldea,
    ypart(2*(i-1)+1)= lowery(i-1);
    ypart(2*(i-1)+2)= uppery(i+1);
end
ypart(2*ldea+1) = lowery(ldea);
ypart(2*ldea+2) = lowery(ldea+1);


plot(xpart,ypart,'c','LineWidth',3);
hold on;

%plot censored data 
if plotcens == 1,
    for i = 1:length(censor),
        tcen = 0;
        if censor(i)== 0,
            tcen    = time(i);
        end
        % find the value
        k = 1;
        while xpart(k)<tcen,
            k = k+1;
        end
        if tcen~=0,
            plot(tcen,ypart(k),'oc');
        end
    end
end

hold off;


