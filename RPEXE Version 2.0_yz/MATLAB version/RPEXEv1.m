function pexeout = RPEXEv1(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% 'EventTime' == A sequence of times where the events occur 
% 'Censor'    == A sequence of dichotomous values indicating 
%                censored or not (0=censored and 1=not censored)
% 'CutTimes'  == A vector of unique, sorted, possible times to 
%                make the cuts. Default is sorted (from small to 
%                large) event times
%                Default == 'EventTime'
% 'Trend'     == An input having indicating the monotonicity assumption
%                -- 0: no monotonic assumption
%                -- 1: failure rate is decreasing over time
%                -- 2: failure rate is increasing over time
%                -- 3: monotonic failure rate
%                -- 4: failure rate is increasing and then decreasing
%                -- 5: failure rate is decreasing and then increasing
%                -- 6: failure rate is increasing and then decreasing with 
%                      the peak removed first
%                -- 7: failure rate is decreasing and then increasing with 
%                      the peak removed first
%                Default == 0
% 'Criticalp'  == The critical (naive) p-value cutoff where all p-values 
%                in the backward elimination that are lower than this 
%                will be regarded as being significant. The prediction of 
%                the survival probability will be made on 100 equally
%                spaced time points within the range of the event times 
%                based on the piecewise exponential estimate determined by 
%                all the changepoints. 
%                Default == -1 (equivalent to NA).
% Output: 
% pexeout.times   ==  times to make the cuts 
% pexeout.pvalues ==  pvalues correspond to the times
% pexeout.times_c ==  critical times to make the cuts
% pexeout.pvalues_c ==  critical p-values that are smaller than the 
% pexeout.trend   ==  trend information
% pexeout.struct  ==  structure information for multiple order restrictions
% pexeout.changet ==  change point in time for umbrella alternatives.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cuttimes_default = 0;
monotone_default = 0;
criticalp_default = 0;

%%%% reset parameter inputs
paramPairs  = varargin(1:end);
for k       = 1:2:length(paramPairs)
  param     = lower(paramPairs{k});
  
%%%% disp(param)
  if (~ischar(param))
    error('Optional parameter names must be strings');
  end
  value     = paramPairs{k+1};

%%%  disp(value) 
  switch (param)
  case 'eventtime'
    times       = value;
    
  case 'censor'
    censoring   = value;
    
  case 'cuttimes'
    cuttimes    = value;
    cuttimes_default = 1;
    
  case 'trend'
     monotone           = value;
     monotone_default   = 1;
     
  case 'criticalp'
     criticalp          = value;
     criticalp_default = 1;
     
  otherwise
      error(['Unrecognized option ' param '.']);
  end % end switch
end % end k loop on data setup

% Compute the time of the death(in increasing order), the total time on
% test, and the number of deaths corresponding to ttot. 
% These quantities will be used later. 
[time_die,ttot,deaths]  = totaltest(times,censoring);

% Set the default values
if cuttimes_default == 0,
    cuttimes   = time_die(1:(length(time_die)-1));  
end

if monotone_default == 0,
    monotone = 0;
end

if  criticalp_default == 0,
    criticalp = -1;
end

if monotone     == 0,
    [ts,pvalues]    = loopcuts(times,censoring,cuttimes,monotone); 
    pexeout.trend   = 'No order restriction';    
elseif monotone == 1,
    [time2, ttot2, deaths2] = pava_dfr(time_die,ttot,deaths);
    cuttimes_trend          = time2(1:(length(time2)-1));    
    cuttimes                = intersect(cuttimes_trend, cuttimes);
    [ts,pvalues]    = loopcuts(times,censoring,cuttimes,monotone); 
    pexeout.trend   = 'Decreasing failure rate';
elseif monotone == 2,    
    [time2, ttot2, deaths2] = pava_ifr(time_die,ttot,deaths);
    cuttimes_trend          = time2(1:(length(time2)-1));    
    cuttimes                = intersect(cuttimes_trend, cuttimes);
    [ts,pvalues]    = loopcuts(times,censoring,cuttimes,monotone);
    pexeout.trend   = 'Increasing failure rate';    
elseif monotone == 3,  % monotonic failure rate
    [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
    cuttimes_trend  = time2(1:(length(time2)-1));    
    cuttimes = intersect(cuttimes_trend, cuttimes);
    if struct(1).loglik > struct(2).loglik,
        monotone = 1;
    else 
        monotone = 2;
    end
    [ts,pvalues]    = loopcuts(times,censoring,cuttimes,monotone);     
    pexeout.struct  = struct;
    pexeout.trend   = 'Monotone failure rate';        
elseif monotone == 4,  % increasing then decreasing failure rate
    [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
%     indx
%     label
    cuttimes_trend  = time2(1:(length(time2)-1));    
    cuttimes        = intersect(cuttimes_trend, cuttimes);
    changetime      = time_die(indx);    
    indsmall        = find(cuttimes< changetime);
    indlarge        = find(cuttimes>changetime);
    cuttimes        = [cuttimes(indsmall);changetime;cuttimes(indlarge)];
    mono            = [2*ones(length(indsmall),1);0;ones(length(indlarge),1)];
    [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono);  
    pexeout.struct  = struct; 
    pexeout.changet = changetime; 
    pexeout.trend   = 'Increasing-decreasing failure rate'; 
    
elseif monotone == 5,   % decreasing then increasing failure rate
    [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
    cuttimes_trend  = time2(1:(length(time2)-1));    
    cuttimes        = intersect(cuttimes_trend, cuttimes);
    changetime      = time_die(indx);    
    indsmall        = find(cuttimes< changetime);
    indlarge        = find(cuttimes>changetime);
    cuttimes        = [cuttimes(indsmall);changetime;cuttimes(indlarge)];
    mono            = [ones(length(indsmall),1);0;2*ones(length(indlarge),1)];
    [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono);      
                        
%     [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
%     cuttimes_trend  = time2(1:(length(time2)-1));    
%     cuttimes        = intersect(cuttimes_trend, cuttimes);
%     changetime      = time_die(indx);
%     indsmall        = find(cuttimes< changetime);
%     indlarge        = find(cuttimes>=changetime);
%     [ts1,pvalues1]  = loopcuts(times,censoring,cuttimes(indsmall),1); 
%                         % decreasing Failure Rate
%     [ts2,pvalues2]  = loopcuts(times,censoring,cuttimes(indlarge),2); 
%                         % increasing Failure Rate
%     ts              = [ts1;ts2];
%     pvalues         = [pvalues1;pvalues2];
    pexeout.struct  = struct; 
    pexeout.changet = changetime; 
    pexeout.trend   = 'Decreasing-increasing failure rate';       
    
elseif monotone == 6,   % increasing then decreasing failure rate, no peak  
    [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
    cuttimes_trend  = time2(1:(length(time2)-1));    
    cuttimes        = intersect(cuttimes_trend, cuttimes);
    changetime      = time_die(indx);    
    indsmall        = find(cuttimes<changetime);
    indlarge        = find(cuttimes>changetime);
    cuttimes        = [cuttimes(indsmall);cuttimes(indlarge)];
    mono            = [2*ones(length(indsmall),1);ones(length(indlarge),1)];
    [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono);  
    pexeout.struct  = struct; 
    pexeout.changet = changetime; 
    pexeout.trend   = 'Increasing-decreasing failure rate';

%     
%     [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
%     cuttimes_trend  = time2(1:(length(time2)-1));    
%     cuttimes        = intersect(cuttimes_trend, cuttimes);
%     changetime      = time_die(indx);    
%     indsmall        = find(cuttimes< changetime);
%     indlarge        = find(cuttimes>changetime);
%     cuttimes        = [cuttimes(indsmall);changetime;cuttimes(indlarge)];
%     mono            = [2*ones(length(indsmall),1);0;ones(length(indlarge),1)];
%     [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono);  
%     pexeout.struct  = struct; 
%     pexeout.changet = changetime; 
%     pexeout.trend   = 'Increasing-decreasing failure rate'; 
    
elseif monotone == 7,   % decreasing then increasing failure rate, no peak 
    [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
    cuttimes_trend  = time2(1:(length(time2)-1));    
    cuttimes        = intersect(cuttimes_trend, cuttimes);
%     cuttimes_trend
%     cuttimes
    changetime      = time_die(indx);    
    indsmall        = find(cuttimes<changetime);
    indlarge        = find(cuttimes>changetime);
%     changetime
    cuttimes        = [cuttimes(indsmall);cuttimes(indlarge)];
%     cuttimes
    mono            = [ones(length(indsmall),1);2*ones(length(indlarge),1)];
%     mono
    [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono); 
%     [ts pvalues]
    pexeout.struct  = struct; 
    pexeout.changet = changetime; 
    pexeout.trend   = 'Decreasing-increasing failure rate';

end

validp = find(pvalues<1);
pexeout.times       = ts(validp);
pexeout.pvalues     = pvalues(validp);



%% The following part is new in the RPEXE version 2.
% 1. In the new version the program will determine location of the 
% change point given a naive p-value.
% 2. The program is able to plot the estimated survival function 
% using the change point determined by the naive p-value and overlay 
% the estimated piecewise exponential estimate with the Kaplan-Meier curve.

if criticalp_default == 1,
    
% compute the grid of 100 times that is equally spaced between 0 and the
% maximum of the event time; 
tmax = max(times);
t100 = [0.01:0.01:1]*tmax;

% Compute the predicted value using the grid and use it as 
if min(pexeout.pvalues) < criticalp, % if the change point is found
        % critical sentence: when any p-value is less than the critical 
        %     value, keep all the times at that point.
    tchange = pexeout.times(min(find(pexeout.pvalues < criticalp)):end);
    [pred100 lamest] = pexeest(times, censoring, tchange, t100);
    % save critical time and p-values;
    tchange = sort(tchange);
    % pexeout.pvalues_c = ...
    % pexeout.pvalues(min(find(pexeout.pvalues < criticalp)):end);
    % compute the critical value at the changetime of tchange;
    [predc lamest] = pexeest(times, censoring, tchange, tchange);
    %pchange = pexeest_p(times, censoring, tchange, mono);  
    % program mono
    if monotone < 4,
        mono = monotone * ones(length(tchange),1);
    else
        indsmall        = find(tchange < changetime);
        indlarge        = find(tchange > changetime);
        if monotone == 4,
        if ismember(changetime, tchange) == 1,
            % order restriction changepoint in the critical times set
            tchange         = [tchange(indsmall);changetime;tchange(indlarge)];
            mono            = [2*ones(length(indsmall),1);0;ones(length(indlarge),1)];
        else 
            % order restriction changepoint not in the critical times set
            tchange         = [tchange(indsmall);tchange(indlarge)];
            mono            = [2*ones(length(indsmall),1);ones(length(indlarge),1)]; 
        end
        elseif monotone == 5,
        if ismember(changetime, tchange) == 1,
            % order restriction changepoint in the critical times set
            tchange         = [tchange(indsmall);changetime;tchange(indlarge)];
            mono            = [ones(length(indsmall),1);0;2*ones(length(indlarge),1)];
        else
            % order restriction changepoint not in the critical times set
            tchange         = [tchange(indsmall);tchange(indlarge)];
            mono            = [ones(length(indsmall),1);2*ones(length(indlarge),1)];        
        end
        elseif monotone == 6,    
        tchange         = [tchange(indsmall);tchange(indlarge)];
        mono            = [2*ones(length(indsmall),1);ones(length(indlarge),1)]; 
        elseif monotone == 7,    
        tchange         = [tchange(indsmall);tchange(indlarge)];
        mono            = [ones(length(indsmall),1);2*ones(length(indlarge),1)];         
        end
    end
    pexeout.times_c = tchange;
    pexeout.lamest  = lamest; 
    pvalall = loopcuts_onestep(times,censoring,tchange,mono);
    pexeout.pvalues_c = pvalall;
    % If an output has p-value 
else
    % Note that if we set the changepoint to be the maximum time in the
    % record, then the function will return the exponential estimate;
    %[pred100 lamest] = pexeest(times, censoring, max(times), t100);
    expar1 = sum(times)/sum(censoring);
    pred100 = exp(-t100/expar1);
    pexeout.lamest  = expar1; 
    
end

% Plot the Kaplan-Meier curve overlaid with the estimates;
% original scale;
figure(11);
[xpart,ypart] = km(times, censoring, 1);
% Note that 0 can be changed to 1 show censoring data;
hold on;
plot(t100,pred100,'--r','LineWidth',2);
if min(pexeout.pvalues) < criticalp,
    plot(tchange,predc,'^r','LineWidth',4);
end
hold off;


% log scale;
figure(12);
[xpart,ypart] = km_log(times, censoring, 1);
% Note that 0 can be changed to 1 show censoring data;
hold on;
plot(t100,log(pred100),'--r','LineWidth',2);
if min(pexeout.pvalues) < criticalp, 
    plot(tchange,log(predc),'^r','LineWidth',4);
end
hold off;


figure(13);
[xpart,ypart] = km(times, censoring, 0);
% Note that 0 can be changed to 1 show censoring data;
hold on;
plot(t100,pred100,'--r','LineWidth',2);
if min(pexeout.pvalues) < criticalp,
    plot(tchange,predc,'^r','LineWidth',4);
end
ylim([0 1]);
legend('Kaplan-Meier estimate','Exponential estimate')
hold off;

% Save the data for makng the
pexeout.plotdatakme.times = times;
pexeout.plotdatakme.censoring = censoring;
pexeout.plotdatapexe.t100 = t100;
pexeout.plotdatapexe.pred100 = pred100;
if min(pexeout.pvalues) < criticalp, 
    pexeout.plotdatapexe.tchange = tchange;
    pexeout.plotdatapexe.predc = predc;
end

end



