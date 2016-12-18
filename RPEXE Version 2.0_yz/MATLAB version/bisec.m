function a2 = bisec(delta,dea1,dea2,upbd,lowbd);
%%% running bisection algorithm to search for a2, the minimizer of 
%%% (log((a2)^dea1*(1-a2)^dea2-delta))^2
%%% upbd and lowbd are the upper and lower bounds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a2init      = 0.5*lowbd+0.5*upbd;
lowvalue    = dea1*log(lowbd) + dea2*log(1-lowbd)-delta;
a2initvalue = dea1*log(a2init)   + dea2*log(1-a2init)-delta;

if lowvalue < a2initvalue, 
   %%% monotone increasing,
  while abs(upbd-lowbd) > 0.00000001, 
    if 0 < a2initvalue,
        upbd = a2init;
        a2init = 0.5*(lowbd+upbd);
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta;
         %lowbd
         %upbd
    elseif 0 > a2initvalue,
        lowbd = a2init;
        lowvalue = a2initvalue;
        a2init = 0.5*(lowbd+upbd);
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta;
         %lowbd
         %upbd
    else
        lowbd = upbd;
        a2 = a2init;
    end
  end
  a2  = 0.5*lowbd+0.5*upbd;
elseif lowvalue > a2initvalue,
   %%% monotone decreasing
  while abs(upbd-lowbd) > 0.00000001, 
    if 0 > a2initvalue,
        upbd = a2init;
        a2init = 0.5*(lowbd+upbd);
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta;
         %lowbd
         %upbd
    elseif 0 < a2initvalue,
        lowbd = a2init;
        lowvalue = a2initvalue;
        a2init = 0.5*(lowbd+upbd);
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta;
         %lowbd
         %upbd
    else
        lowbd = upbd;
        a2 = a2init;
    end
  end
  a2  = 0.5*lowbd+0.5*upbd;
else 
  %%% no change, means lowvalue == a2initvalue == 0;
  %%%                  a2 = 1 if upbd = 1; 
  %%%                  a2 = 0 if upbd < 1. 
  if upbd == 1,
      a2 = upbd;
  else 
      a2 = lowbd;
  end
end
