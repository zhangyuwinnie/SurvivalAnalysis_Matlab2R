############################################################
#Input:
# 'EventTime' == A sequence of times where the events occur
# 'Censor'    == A sequence of dichotomous values indicating
#                censored or not (0=censored and 1=not censored)
# 'CutTimes'  == A vector of unique, sorted, possible times to
#                make the cuts. Default is sorted (from small to
#                large) event times
#                Default == 'EventTime'
# 'Trend'     == An input having indicating the monotonicity assumption
#                -- 0: no monotonic assumption
#                -- 1: failure rate is decreasing over time
#                -- 2: failure rate is increasing over time
#                -- 3: monotonic failure rate
#                -- 4: failure rate is increasing and then decreasing
#                -- 5: failure rate is decreasing and then increasing
#                -- 6: failure rate is increasing and then decreasing with
#                      the peak removed first
#                -- 7: failure rate is decreasing and then increasing with
#                      the peak removed first
#                Default == 0
#Output:
# pexeout.times   ==  times to make the cuts
# pexeout.pvalues ==  pvalues correspond to the times
# pexeout.trend   ==  trend information
# pexeout.struct  ==  structure information for multiple order restrictions
# pexeout.changet ==  change point in time for umbrella alternatives.
#############################################################################

RPEXEv1 <- function(EventTime=NA,eventtime,Censor=NA,censor,CutTimes=NA,cuttime=1,Trend=NA,trend=0)
{


  cuttimes_default = 0
  monotone_default = 0
  pexeout=list(times=as.null(),pvalues=as.null(),trend=as.null(),struct=as.null(),changet=as.null())



  #reset parameter inputs
    if (!is.na(EventTime))
       times= eventtime
    if (!is.na(Censor))
       censoring= censor
    if (!is.na(CutTimes))
      {
       cuttimes    = cuttime
       cuttimes_default = 1
       }
    if (!is.na(Trend))
      {
       monotone= trend
       monotone_default   = 1
       }
    if (is.na(EventTime)||is.na(Censor))
       cat("Unrecognized option")
  # Compute the time of the death(in increasing order), the total time on
  # test, and the number of deaths corresponding to ttot.
  # These quantities will be used later.
  returnv=totaltest(times,censoring) #returnv=list(time_die,ttot,deaths),the variables owns the same length
  m=dim(returnv)[2]/3
  time_die=returnv[,1:m]
  ttot=returnv[,(m+1):(2*m)]
  deaths=returnv[,(2*m+1):3*m]

  # Set the default values
  if (cuttimes_default == 0)
     cuttimes   = time_die[1:(length(time_die)-1)]
  if (monotone_default == 0)
     monotone = 0

  if (monotone== 0)
     {
     returnva_l=loopcuts(times,censoring,cuttimes,monotone)
     n=dim(returnva_l)[2]/2
     ts=returnva_l[,1:n]
     pvalues=returnva_l[,(n+1):2*n]
     pexeout$trend="No order restriction"
     }
  if (monotone == 1)
     {
     returnva_p=pava_dfr(time_die,ttot,deaths)
     n=dim(returnva_p)[2]/3
     time2=returnva_p[,1:n]
     ttot2=returnva_p[,(n+1):(2*n)]
     deaths2=returnva_p[,(2*n+1):3*n]
     cuttimes_trend= time2[1:(length(time2)-1)]
     cuttimes= sort(intersect(cuttimes_trend, cuttimes))
     returnva_l=loopcuts(times,censoring,cuttimes,monotone)
     n=dim(returnva_l)[2]/2
     ts=returnva_l[,1:n]
     pvalues=returnva_l[,(n+1):2*n]
     pexeout$trend="Decreasing falilure rate"
     }
  if (monotone == 2)
     {
     returnva_p= pava_ifr(time_die,ttot,deaths)
     n=dim(returnva_p)[2]/3
     time2=returnva_p[,1:n]
     ttot2=returnva_p[,(n+1):(2*n)]
     deaths2=returnva_p[,(2*n+1):3*n]
     cuttimes_trend=time2[1:(length(time2)-1)]
     cuttimes= sort(intersect(cuttimes_trend, cuttimes))
     returnva_l=loopcuts(times,censoring,cuttimes,monotone)
     n=dim(returnva_l)/2
     ts=returnva_l[,1:n]
     pvalues=returnva_l[,(n+1):2*n]
     pexeout$trend="Increasing falilure rate"
     }
  if (monotone == 3)
     {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes =sort(intersect(cuttimes_trend, cuttimes))
      if (struct[[1]][[4]][1]>struct[[1]][[4]][2])
         monotone = 1
      else
         monotone = 2
      returnva_l=loopcuts(times,censoring,cuttimes,monotone)
      n=dim(returnva_l)[2]/2
      ts=returnva_l[,1:n]
      pvalues=returnva_l[,(n+1):(2*n)]
      pexeout$struct= struct
      pexeout$trend="Monotone failure rate"
      }
  if (monotone==4) #increasing then decreasing failure rate
      {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      #indx
      #label
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes        = sort(intersect(cuttimes_trend,cuttimes))
      indx=indx[[1]][1]
      changetime      = time_die[indx]
      cuttimes        = c(cuttimes[cuttimes< changetime],changetime,cuttimes[cuttimes>changetime])
      mono            = c(2*matrix(1,nrow=length(cuttimes[cuttimes<changetime]),ncol=1),0,matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_l=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_l)[2]/2
      ts=returnva_l[,1:n]
      pvalues=returnva_l[,(n+1):(2*n)]
      pexeout$struct= struct
      pexeout$changet = changetime
      pexeout$trend   = "Increasing-decreasing failure rate"
     }
  if  (monotone==5) #decreasing then increasing failure rate
     {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes =sort(intersect(cuttimes_trend, cuttimes))
      indx=indx[[1]][1]
      changetime      = time_die[indx]
      cuttimes        = c(cuttimes[cuttimes< changetime],changetime,cuttimes[cuttimes>changetime])
      mono            = c(matrix(1,nrow=length(cuttimes[cuttimes< changetime]),ncol=1),0,2*matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_lu=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_lu)[2]/2
      ts=returnva_lu[,1:n]
      pvalues=returnva_lu[,(n+1):2*n]
  #     [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
  #     cuttimes_trend  = time2(1:(length(time2)-1));
  #     cuttimes        = intersect(cuttimes_trend, cuttimes);
  #     changetime      = time_die(indx);
  #     indsmall        = find(cuttimes< changetime);
  #     indlarge        = find(cuttimes>=changetime);
  #     [ts1,pvalues1]  = loopcuts(times,censoring,cuttimes(indsmall),1);
  #                         % decreasing Failure Rate
  #     [ts2,pvalues2]  = loopcuts(times,censoring,cuttimes(indlarge),2);
  #                         % increasing Failure Rate
  #     ts              = [ts1;ts2];
  #     pvalues         = [pvalues1;pvalues2];

      pexeout$struct  = struct
      pexeout$changet = changetime
      pexeout$trend   = "Decreasing-increasing failure rate"
      }
  if (monotone==6) # increasing then decreasing failure rate, no peak
     {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2] #struct=data.frame(time,ttot,deaths,loglik)
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes =sort(intersect(cuttimes_trend, cuttimes))
      indx=indx[[1]][1]
      changetime      = time_die[indx]
      cuttimes        = c(cuttimes[cuttimes< changetime],cuttimes[cuttimes>changetime])
      mono            = c(2*matrix(1,nrow=length(cuttimes[cuttimes< changetime]),ncol=1),matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_lu=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_lu)[2]/2
      ts=returnva_lu[,1:n]
      pvalues=returnva_lu[,(n+1):2*n]
      pexeout$struct  = struct
      pexeout$changet = changetime
      pexeout$trend   ="Increasing-decreasing failure rate"
     }

  #     [time2,struct,label,indx] = umbrella(time_die,ttot,deaths,monotone-3);
  #     cuttimes_trend  = time2(1:(length(time2)-1));
  #     cuttimes        = intersect(cuttimes_trend, cuttimes);
  #     changetime      = time_die(indx);
  #     indsmall        = find(cuttimes< changetime);
  #     indlarge        = find(cuttimes>changetime);
  #     cuttimes        = [cuttimes(indsmall);changetime;cuttimes(indlarge)];
  #     mono            = [2*ones(length(indsmall),1);0;ones(length(indlarge),1)];
  #     [ts,pvalues]    = loopcuts_umbrella(times,censoring,cuttimes,mono);
  #     pexeout.struct  = struct;
  #     pexeout.changet = changetime;
  #     pexeout.trend   = 'Increasing-decreasing failure rate';

  if (monotone == 7) #decreasing then increasing failure rate, no peak
      {
      returnva_u=umbrella(time_die,ttot,deaths,monotone-3)
      time2=returnva_u[1]
      struct=returnva_u[2]
      label=returnva_u[3]
      indx=returnva_u[4]
      time2=time2[[1]][]
      cuttimes_trend  = time2[1:(length(time2)-1)]
      cuttimes        = sort(intersect(cuttimes_trend,cuttimes))
  #   cuttimes_trend
  #   cuttimes
  #   changetime
      indx=indx[[1]][1]
      changetime = time_die[indx]
  #   cuttimes
      cuttimes        = c(cuttimes[cuttimes<changetime],cuttimes[cuttimes>changetime])
  #   mono
      mono=c(matrix(1,nrow=length(cuttimes[cuttimes<changetime]),ncol=1),2*matrix(1,nrow=length(cuttimes[cuttimes>changetime]),ncol=1))
      returnva_lu=loopcuts_umbrella(times,censoring,cuttimes,mono)
      n=dim(returnva_lu)[2]/2
      ts=returnva_lu[,1:n]
      pvalues=returnva_lu[,(n+1):2*n]
      pexeout$struct  = struct
      pexeout$changet = changetime
      pexeout$trend   ="Decreasing-increasing failure rate"
      }
  pexeout$times       = unique(ts[pvalues<1])
  pexeout$pvalues     = unique(pvalues[pvalues<1])
  return(pexeout)
}