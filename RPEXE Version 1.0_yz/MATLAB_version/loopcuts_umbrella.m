function [ts, pvalues] = loopcuts_umbrella(time,censor,cuttimes,mono)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% using a loop format to find out the times to make the cuts
%%% from large pvalues to small pvalues and the lambdas in between
%%%
%%% Input: 
%%%        time: A sequence of time 
%%%        censor: and censor (0=censored and 1=not censored)
%%%        cuttimes: unique, sorted, possible times to make the cuts, 
%%%              including 0 and the ending time
%%%        mono: indicate the type of the test, a vector of the size of "cuttimes"
%%%            mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
%%%                 == 1: 1-sides hypothesis: H0:lam1>=lam2; H1:lam1 < lam2
%%%                         decreasing failure rate constraint
%%%                 == 2: 1-sides hypothesis: H0:lam1<=lam2; H1:lam1 > lam2
%%%                         increasing failure rate constraint
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
ts      = zeros(cutlen,1);
pvalues = zeros(cutlen,1);
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
    ts = cuttimes;
    pvalues = p;
else 
    for i = 1:cutlen,
        if i == 1,
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
            maxp = max(pvalall);
            maxplab = 0;
            for k1 = 1:length(pvalall),
                if maxp == pvalall(k1),
                    maxplab = k1;
                end
            end
            %pvalall;
            %%% save ts(i) and pvalues(i) 
            ts(i) = allt(maxplab);
            pvalues(i) = pvalall(maxplab);
            allt(maxplab) = [];
            pvalall(maxplab) = [];   
            monall(maxplab) = [];
            %%% merge centain items on ttotvec and deavec
            ttotvec(maxplab) = ttotvec(maxplab) + ttotvec(maxplab+1);
            ttotvec(maxplab+1) = [];
            deavec(maxplab)  = deavec(maxplab) + deavec(maxplab+1);
            deavec(maxplab+1) = [];
        else
            %%% Use position maxplab, two vectors:(allt,pvalall), and
            %%% ttotvec, deavec to compute ts(i) and pvalues(i)
            if maxplab == 1,
                %p = lrtpvalue(ttotvec(1),ttotvec(2),deavec(1),...
                %        deavec(2),n);
                [a2,p] = exact_pvalue(ttotvec(1),ttotvec(2),deavec(1),...
                    deavec(2),monall(1));
                if length(deavec) == 2,
                    pvalues(i) = p;
                    ts(i) = allt(maxplab);
                else
                    pvalall(1)  = p;
                    maxp        = max(pvalall);
                    maxplab     = 0;
                    for k1 = 1:length(pvalall),
                        if maxp == pvalall(k1),
                            maxplab     = k1;
                        end
                    end
                    %%% save ts(i) and pvalues(i) 
                    ts(i) = allt(maxplab);
                    pvalues(i) = pvalall(maxplab);
                    allt(maxplab)    = [];
                    %pvalall;
                    pvalall(maxplab) = [];
                    %monall;
                    monall(maxplab) = [];
                    %%% merge centain items on ttotvec and deavec
                    ttotvec(maxplab) = ttotvec(maxplab) + ttotvec(maxplab+1);
                    ttotvec(maxplab+1) = [];
                    deavec(maxplab)  = deavec(maxplab) + deavec(maxplab+1);
                    deavec(maxplab+1) = [];
                end   
            elseif maxplab > length(allt),
                [a2,p] = exact_pvalue(ttotvec(maxplab-1),ttotvec(maxplab),...
                    deavec(maxplab-1),deavec(maxplab),monall(maxplab-1)); 
                pvalall(length(allt))  = p;
                maxp        = max(pvalall);
                maxplab     = 0;
                for k1 = 1:length(pvalall),
                    if maxp == pvalall(k1),
                        maxplab     = k1;
                    end
                end
                %%% save ts(i) and pvalues(i) 
                ts(i) = allt(maxplab);
                pvalues(i) = pvalall(maxplab);
                allt(maxplab) = [];
                %pvalall;
                pvalall(maxplab) = [];
                %monall;
                monall(maxplab) = [];
                %%% merge centain items on ttotvec and deavec
                ttotvec(maxplab) = ttotvec(maxplab) + ttotvec(maxplab+1);
                ttotvec(maxplab+1) = [];
                deavec(maxplab)  = deavec(maxplab) + deavec(maxplab+1);
                deavec(maxplab+1) = [];
            else 
                [a2,pfront] = exact_pvalue(ttotvec(maxplab-1),ttotvec(maxplab),...
                        deavec(maxplab-1),deavec(maxplab),monall(maxplab-1));
                [a2,pback]  = exact_pvalue(ttotvec(maxplab),ttotvec(maxplab+1),...
                        deavec(maxplab),deavec(maxplab+1),monall(maxplab));                     
                pvalall(maxplab-1) = pfront;
                pvalall(maxplab)   = pback;
                maxp        = max(pvalall);
                maxplab     = 0;
                for k1 = 1:length(pvalall),
                    if maxp == pvalall(k1),
                        maxplab     = k1;
                    end
                end
                %%% save ts(i) and pvalues(i) 
                ts(i) = allt(maxplab);
                pvalues(i) = pvalall(maxplab);
                allt(maxplab) = []; 
                %pvalall;
                pvalall(maxplab) = [];
                %monall;
                monall(maxplab) = [];                
                %%% merge centain items on ttotvec and deavec
                ttotvec(maxplab) = ttotvec(maxplab) + ttotvec(maxplab+1);
                ttotvec(maxplab+1) = [];
                deavec(maxplab)  = deavec(maxplab) + deavec(maxplab+1);
                deavec(maxplab+1) = [];
            end                
        end
    end
end

