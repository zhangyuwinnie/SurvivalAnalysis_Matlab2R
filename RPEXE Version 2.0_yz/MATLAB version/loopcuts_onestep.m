function pvalall = loopcuts_onestep(time,censor,cuttimes,mono)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% using a loop format to find out the times to make the cuts
%%% from large pvalues to small pvalues and the lambdas in between
%%%
%%% Input: 
%%%        time: A sequence of time 
%%%        censor: and censor (0=censored and 1=not censored)
%%%        cuttimes: unique, sorted, possible times to make the cuts, 
%%%              including 0 and the ending time
%%%        mono: a vector indicating the type of the test, 
%%%              of the size of "cuttimes"
%%%             mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
%%%                  == 1: 1-sides hypothesis: H0:lam1>=lam2; H1:lam1 < lam2
%%%                        decreasing failure rate constraint
%%%                  == 2: 1-sides hypothesis: H0:lam1<=lam2; H1:lam1 > lam2
%%%                        increasing failure rate constraint
%%% Output: 
%%%        ts: the times where the cuts shall be made 
%%%        pvalues: the p values for deleting each cutting times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Sort cuttimes and find cutlen
cuttimes = unique(sort(cuttimes));
cutlen   = length(cuttimes);

% prepare the death time, ttot, and the number of deaths 
[time_die,ttot,deaths] = totaltest(time,censor);
              
% Based on cuttimes, find out the corresponding totaltime on test and 
% the number of deaths.

%%% Set up ts and pvalues
if cutlen == 1,
    for j = 1:length(time_die),
        if time_die(j) == cuttimes,
            ttot1 = sum(ttot(1:j));
            ttot2 = sum(ttot(j+1:end));
            death1 = sum(deaths(1:j));
            death2 = sum(deaths(j+1:end));
            %p = lrtpvalue(ttot1,ttot2,death1,death2,n);
            [a2,p] = exact_pvalue(ttot1,ttot2,death1,death2,mono);
        end
    end
    allt = cuttimes;
    pvalall = p;
else 
    %%% initialize the matrix storing ttots and deaths and p values
    ttotvec = zeros(cutlen+1,1);
    deavec  = zeros(cutlen+1,1);
    allt    = cuttimes;
    pvalall = zeros(cutlen,1);
    monall  = mono;
        for ii = 1:cutlen,
            for j = 1:length(time_die),
                if time_die(j) == cuttimes(ii),
                   if ii == 1,
                       ttotvec(ii) = sum(ttot(1:j));
                        deavec(ii)  = sum(deaths(1:j));
                   else 
                       ttotvec(ii) = sum(ttot(1:j)) - sum(ttotvec(1:(ii-1)));
                       deavec(ii)  = sum(deaths(1:j)) - sum(deavec(1:(ii-1)));
                   end
                end
            end
        end
    %%% do the last item in ttotvec and deavec
    ttotvec(cutlen+1) = sum(ttot)-sum(ttotvec(1:cutlen));
    deavec(cutlen+1)  = sum(deaths)-sum(deavec(1:cutlen)); 
    %%% Compute ts(1) and pvalues(1)
    for k = 1:length(allt),
        ttot1 = ttotvec(k);
        ttot2 = ttotvec(k+1);
        dea1  = deavec(k);
        dea2  = deavec(k+1);
        [a2,pvalall(k)] = exact_pvalue(ttot1,ttot2,dea1,dea2,monall(k));
    end
end
