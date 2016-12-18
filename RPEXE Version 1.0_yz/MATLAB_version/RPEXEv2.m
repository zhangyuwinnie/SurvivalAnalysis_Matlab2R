function pexeout = RPEXEv2(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% 'Times'     == a sequence of times where the events occur 
% 'TTOT'      == a vector containing the total time on test
% 'Death'     == number of death from the previous time to 
%                the current time 
% 'CutTimes'  == a vector of unique, sorted, possible times to 
%                make the cuts. Default is sorted (from small to 
%                large) event times
%                default == 'EventTime'
% 'Monotone'  == an input having three levels indicating 
%                the monotonic assumption
%                -- 0: there is no monotonic assumption
%                -- 1: the normalized spacing is decreasing over time
%                -- 2: the normalized spacing is increasing over time
%                default == 0
% Output: 
% pexeout.times   ==  times to make the cuts 
% pexeout.pvalues ==  pvalues correspond to the times
% pexeout.powers  ==  power of the tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cuttimes_default = 0;
monotone_default = 0;

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
  case 'times'
    time_die    = value;
    
  case 'ttot'
    ttot        = value;
    
  case 'death'
    deaths      = value;
    
  case 'cuttimes'
    cuttimes    = value;
    cuttimes_default = 1;
    
  case 'monotone'
     monotone           = value;
     monotone_default   = 0;
      
  otherwise
      error(['Unrecognized option ' param '.']);
  end % end switch
end % end k loop on data setup


% Set the default values
if cuttimes_default == 0,
    cuttimes   = time_die(1:(length(time_die)-1));  
end

if monotone_default == 0,
    monotone = 0;
end


if monotone ==0 || monotone == 1 || monotone == 2,
    [ts,pvalues,powers]    = loopcuts_ttot(time_die,ttot,deaths,cuttimes);
elseif monotone == 3,
    [time2, struct] = umbrella(time_die,ttot,deaths,monotone-3);
    cuttimes_trend = time2(1:(length(time2)-1));    
    cuttimes = intersect(cuttimes_trend, cuttimes);
    if struct(1).loglik > struct(2).loglik,
        monotone = 1;
    else 
        monotone = 2;
    end
    [ts,pvalues] = loopcuts(times,censoring,cuttimes,monotone);     
    pexeout.struct = struct;
end

pexeout.times       = ts;
pexeout.pvalues     = pvalues;




