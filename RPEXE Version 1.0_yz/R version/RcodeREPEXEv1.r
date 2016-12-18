bisec=function(delta,dea1,dea2,upbd,lowbd)
{
#Running bisection algorithm to search for a2, the minimizer of (log((a2)^dea1*(1-a2)^dea2-delta))^2
#Upbd and lowbd are the upper and lower bounds.
#########################################################

a2init=0.5*(lowbd+upbd)
lowvalue=dea1*log(lowbd)+dea2*log(1-lowbd)-delta
a2initvalue=dea1*log(a2init)+dea2*log(1-a2init)-delta
a=a2initvalue
l=lowvalue
if (l<a)
 {
  #monotone increasing
  while(abs(upbd-lowbd)> 0.00000001)
   {
    n=a2initvalue
    if (0<n)
        {
        upbd = a2init
        a2init = 0.5*(lowbd+upbd)
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
        #lowbd
        #upbd
        }
     if (0>n)
        {
        lowbd = a2init
        lowvalue = a2initvalue
        a2init = 0.5*(lowbd+upbd)
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
        #lowbd
        #upbd
        }
      if (0==n)
        {
        lowbd = upbd
        a2 = a2init
        }
     }
  a2=0.5*(lowbd+upbd)
  }
if (l>a)
  {
   #monotone decreasing
   while (abs(upbd-lowbd)>0.00000001)
   {
     n=a2initvalue
    if (0>n)
       {
        upbd = a2init
        a2init = 0.5*(lowbd+upbd)
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
         #lowbd
         #upbd
        }
    if (0<n)
        {
        lowbd = a2init
        lowvalue = a2initvalue
        a2init = 0.5*(lowbd+upbd)
        a2initvalue = log(a2init)*dea1+log(1-a2init)*dea2-delta
         #lowbd
         #upbd
        }
    if (0==n)
        {
        lowbd = upbd
        a2 = a2init
        }
    }
  a2=0.5*lowbd+0.5*upbd
  }
if (l==a)
  {
#no change, means lowvalue == a2initvalue == 0;
#                 a2 = 1 if upbd = 1; 
#                 a2 = 0 if upbd < 1. 
    if (upbd == 1)
      a2 = upbd
    else 
      a2 = lowbd
   }
  return(a2)
}

exact_pvalue=function(ttot1,ttot2,dea1,dea2,mono)
{
#Compute the exact P value from the likelihood ratio test
#Input:
#       ttot1, ttot2 = total time on test 1 and 2
#       dea1, dea2   = number of death 1 and 2
#       mono: indicate the type of the test
#           mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
#                == 1: 1-sided hypothesis: H0:lam1>=lam2; H1:lam1 < lam2
#                == 2: 1-sided hypothesis: H0:lam1<=lam2; H1:lam1 > lam2
# Output: 
#       pval         = the exact Pvalue for testing H0: lambda1 = lambda2
# Assumption: 
#       ttot1 and ttot2 are from Gamma(dea1,lambda1) and
#       Gamma(dea2,lambda2)
# Compute the test statistic, delta = (ttot1/(ttot1+ttot2))^dea1*...
#                                     (ttot2/(ttot1+ttot2))^dea2;

a1=ttot1/(ttot1+ttot2)
delta=(log(a1)*dea1)+(log(1-a1)*dea2)

#compute the p-value
if(dea1==1&dea2==1)
 {
  # if dea1 = dea2 = 1
  if (mono==0)
    {
     if(a1<0.5)
        pval = 2*pbeta(a1,dea1,dea2)
     else
        pval = 2*(1-pbeta(a1,dea1,dea2))
     }
   if (mono==1)
    {
     pval = pbeta(a1,dea1,dea2)
     pval = pval/pbeta(dea1/(dea1+dea2),dea1,dea2)
     }
   if (mono==2)
    {
     pval = 1-pbeta(a1,dea1,dea2)
     pval = pval/(1-pbeta(dea1/(dea1+dea2),dea1,dea2))
    }
   a2=1-a1
  }
if (dea1!=1||dea2!=1)
 {
  #search for a2
  if (a1<(dea1/(dea1+dea2)))
    {
      upbd = 1
      lowbd = dea1/(dea1+dea2)
      a2=bisec(delta,dea1,dea2,upbd,lowbd)
     }
   if (a1>=(dea1/(dea1+dea2)))
     {
      upbd = dea1/(dea1+dea2)
      lowbd =0.0000000000000000001
      a2=bisec(delta,dea1,dea2,upbd,lowbd)
      }
    #compute the probability in the beta distribution: 
    #P(bigger than max(a1,a2),smaller than min(a1,a2)) 
   if (mono==0)
        pval = pbeta(min(a1,a2),dea1,dea2)+(1-pbeta(max(a1,a2),dea1,dea2))
   if (mono==1)
      {
        pval = pbeta(a1,dea1,dea2)
        pval = pval/pbeta(dea1/(dea1+dea2),dea1,dea2)
       }
   if (mono==2)
       {
        pval = 1-pbeta(a1,dea1,dea2)
        pval = pval/(1-pbeta(dea1/(dea1+dea2),dea1,dea2))
       }
  }
returnval=c(a2,pval)
return(returnval)
}


gamllik=function(structtime,structttot,structdeaths,time_die,ttot,deaths)
{
##############################################################################
#A function computing the log likelihood from the gamma distribution under
#       an order restriction reduction 
# Inputs
#       structtime   == times under restriction
#       structttot   == ttots between each time point and the previous 
#                           time point (or 0) under restriction
#       structdeaths == number of deaths corresponding to struct.ttot
#       time_die      == all possible times to make the cuts
#       ttot          == ttots corresponding to time_die
#       deaths        == dealthes corresponding to ttot
# Outputs      
#       loglik        == log of the likelihood
###############################################################################
#compute the Gamma parameter
structgamindi=array(0,c(length(structtime),1))
for(j in 1:length(structtime))
    structgamindi[j]= structttot[j]/structdeaths[j]
#the likelihood
#get the indices of the cut time:
structindi=array(0,c(length(structtime),1))
for(j in 1:length(structtime))
  for (jj in 1:length(time_die))
     if (time_die[jj]==structtime[j])
         structindi[j] =jj

#set the scale parameter for Gamma distribution of the ttot
 structgampar=array(0,c(length(time_die),1))
for (ii in 1:length(time_die))
    for (j in 1:length(structtime))
        if (ii<=structindi[j])
          {
            if (j==1)
                structgampar[ii]= structgamindi[j]
            else 
                if (ii>structindi[j-1])
                    structgampar[ii] = structgamindi[j]
           }
loglik = 0
for (ii in 1:length(structgampar))
    loglik = loglik + log(dgamma(ttot[ii]/deaths[ii],shape=deaths[ii],scale=structgampar[ii],log=FALSE))
return(loglik)
}


totaltest=function(time,censor)
{
#Function 'totaltest' computes total time on test
#Inputs: times and censor (0 = censored; 1 = uncensored)
#Return: time_die: times that events occur (in ascending order)
#        ttot: total time on test at each time point in time_die
#        deaths: number of death at each time point in time_die


#save the data in the orginal structure
sew=array(0,c(length(time),1))
for(i in 1:length(time))
   sew[i]=i
tmpdata=array(cbind(sew,censor,time),c(length(sew),10))

#sort the tmp data
tmp2 =array(0,c(dim(tmpdata)))
for (i in 1:dim(tmp2)[1])
   {
    if (i!=dim(tmp2)[1])
    {
    tmp2[i,3]=min(tmpdata[,3])
    j=which.min(tmpdata[,3])
    tmp2[i,1:2]=tmpdata[j,1:2]
    tmpdata=tmpdata[-j,]
     }
   else
     tmp2[i,]=tmpdata
   }
#Compute alpha's for the sequence
for (i in 1:dim(tmp2)[1])
   {
   if (tmp2[i,2]==0)
      tmp2[i,4]=0
   else
      tmp2[i,4]=1
    }
#Deal with alpha > 1
for (i in 1:(dim(tmp2)[1]-1))
    if (tmp2[dim(tmp2)[1]+1-i,3]== tmp2[dim(tmp2)[1]-i,3])
      {
       tmp2[dim(tmp2)[1]-i,4]= tmp2[dim(tmp2)[1]-i,4] + tmp2[dim(tmp2)[1]+1-i,4]
       tmp2[dim(tmp2)[1]+1-i,4]=0
      }
# Delete the repeats
k=as.null()
for (i in 1:dim(tmp2)[1])
    if (tmp2[i,2] == 1&tmp2[i,4]==0)
        k[length(k)+1]=i
tmp3 = tmp2
if (length(k)!=0)
   tmp3=tmp3[-k,]


#Compute the number of patients in the study
for(i in 1:dim(tmp3)[1])
  {
   if (tmp3[i,4]==0)
       tmp3[i,5]= 1
   else
       tmp3[i,5]= tmp3[i,4]
   }
for(i in 1:dim(tmp3)[1]-1)
  {
   tmp3[i,6]= sum(tmp3[,5])-sum(tmp3[1:i,5])
   tmp3[dim(tmp3)[1],6]= 0
  }


#Compute the survival time of this cohort
for(i in 1:dim(tmp3)[1])
   {
   if (i==1)
      tmp3[i,7]= sum(tmp3[,5])*tmp3[i,3]
   else
      ###survival time == patient number * time difference
      tmp3[i,7] = tmp3[i-1,6]* (tmp3[i,3]-tmp3[i-1,3])
   }
tmp3[,8] = tmp3[,7]
for (i in 1:dim(tmp3)[1])
   if (tmp3[i,2]==0)
     {
       if (t(tmp3[i:dim(tmp3)[1],2])%*%tmp3[i:dim(tmp3)[1],2]>0)
           tmp3[i+1,8] = tmp3[i,8]+tmp3[i+1,8]
       if (t(tmp3[i:dim(tmp3)[1],2])%*%tmp3[i:dim(tmp3)[1],2]==0&tmp3[i-1,2]!=0)
             {
            ### put all the credit to the last noncensered data
             k = length(tmp3[i:dim(tmp3)[1],2])
             for (j in 1:k)
                tmp3[i-1,8] = tmp3[i-1,8]+tmp3[i-1+j,8]
              }
       } 


#Build the survival reaction
tmp3[,9] = tmp3[,6]
for (i in 2:length(tmp3[,9]))
   tmp3[i,9]= tmp3[i-1,9]-tmp3[i,4]


#plot (tmp3[,3],tmp3[,9])

###delete all the censered items
k=as.null()
for (i in 1:length(tmp3[,1]))
   if (tmp3[i,2]== 0)
      k[length(k)+1]=i
tmp4 = tmp3
if (length(k)!=0)
   tmp4=tmp4[-k,]


time_die=tmp4[,3]
ttot=  tmp4[,8]
deaths=  tmp4[,5]
returnv=cbind(time_die,ttot,deaths)
return(returnv)
}


pava_dfr=function(time_die,ttot,deaths)
{
#########################################################################
#merge certain entries to make the sequence of ttot to be non decreasing
#Input: 
#      time_die == a sequence of times where deaths happened.
#      ttot     == the total time on test at each time point.
#      deaths   == the number of deaths at each time point.
#Output:
#      time2    == the merged time_die
#      ttot2    == .......... ttot
#      deaths2  == .......... deaths
##########################################################################

len=length(ttot)

if (len == 1)
 {
   time2   = time_die
   ttot2   = ttot
   deaths2 = deaths
   }
if (len!=1)
 {  
 for (j in 1:len)
  {
    #indicate deleted items 
   n=length(ttot)
   if(n!=1)
   {
    for (i in 1:(n-1))
       #indicate items to be deleted
       if ((ttot[i]/deaths[i])>(ttot[i+1]/deaths[i+1]))
          {
           deaths[i+1] = deaths[i]+deaths[i+1]
           deaths[i]= 0
           ttot[i+1]=ttot[i]+ttot[i+1]
           ttot[i]=0
           }
     #delete the indicated items
     k=as.null()
     for (i in 1:length(ttot))
        if(ttot[i]==0)
          k[length(k)+1]=i
    if (!is.null(k))
      {
       ttot=ttot[-k]
       deaths=deaths[-k]
       time_die=time_die[-k]
       }
    }
   }
 }
 time2   = time_die
 ttot2   = ttot
 deaths2 = deaths
 returnval=cbind(time2,ttot2,deaths2)
 return(returnval)
}



pava_ifr=function(time_die,ttot,deaths)
{
########################################################
# merge certain entries to make the sequence of ttot/deaths to be non increasing 
# (Note that the pava function makes it non decreasing.This function directly uses function pava().)
#
#Input: 
#      time_die == a sequence of times where deaths happened.
#      ttot     == the total time on test at each time point.
#      deaths   == the number of deaths at each time point.
#Output:
#      time2    == the merged time_die
#      ttot2    == .......... ttot
#      deaths2  == .......... deaths
############################################################

ttotrev= (-1)*ttot
returnval=pava_dfr(time_die,ttotrev,deaths)
m=dim(returnval)[2]/3
time2=returnval[,1:m]
ttot3=returnval[,(m+1):2*m]
deaths2=returnval[,(2*m+1):3*m]
#[time2, ttot3, deaths2] = retuenval 
ttot2=-ttot3
returnval_t=cbind(time2,ttot2,deaths2)
return(returnval_t)
}


loopcuts=function(time,censor,cuttimes,mono)
{
###################################################################(have problem)
#using a loop format to find out the times to make the cuts
#from large pvalues to small pvalues and the lambdas in between
#
#Input: 
#        time: A sequence of time 
#        censor: and censor (0=censored and 1=not censored)
#        cuttimes: unique, sorted, possible times to make the cuts, 
#             including 0 and the ending time
#       mono: indicate the type of the test
#            mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
#                 == 1: 1-sides hypothesis: H0:lam1>=lam2; H1:lam1 < lam2
#                         decreasing failure rate constraint
#                 == 2: 1-sides hypothesis: H0:lam1<=lam2; H1:lam1 > lam2
#                         increasing failure rate constraint
#Output: 
#        ts: the times where the cuts shall be made 
#        pvalues: the p values for deleting each cutting times
####################################################################

#Sort cuttimes and find cutlen
cuttimes=unique(sort(cuttimes))
cutlen=length(cuttimes)

#prepare the death time, ttot, and the number of deaths.returnval_1 is the four returned valure from function totaltest
returnval_1=totaltest(time,censor)
m=dim(returnval_1)[2]/3
time_die=returnval_1[,1:m]
ttot=returnval_1[,(m+1):(2*m)]
deaths=returnval_1[,(2*m+1):3*m]
#Based on cuttimes, find out the corresponding totaltime on test and the number of deaths.

#Set up ts and pvalues
ts=array(0,c(cutlen,1))
pvalues=array(0,c(cutlen,1))
if (cutlen==1)
{
  for(j in 1:length(time_die))
     if (time_die[j]==cuttimes)
      {
        ttot1 = sum(ttot[1:j])
        ttot2 = sum(ttot[-(1:j)])
        death1 = sum(deaths[1:j])
        death2 = sum(deaths[-(1:j)])
        returnval_2= exact_pvalue(ttot1,ttot2,death1,death2,mono)
        a2=returnval_2[1]
        p=returnval_2[2]
       }
  ts= cuttimes
  pvalues = p
  }
if (cutlen!=1)
   {
   for(i in 1:cutlen)
    {
     if (i==1)
       {
        ###initialize the matrix storing ttots and deaths and p values
        ttotvec= array(0,c(cutlen+1,1))
        deavec= array(0,c(cutlen+1,1))
        allt= cuttimes
        pvalall=array(0,c(cutlen,1))
        for(ii in 1:cutlen)
           for(j in 1:length(time_die))
               if (time_die[j]==cuttimes[ii])
                {  
                 if (ii==1)
                    {
                     ttotvec[ii] = sum(ttot[1:j])
                     deavec[ii]  = sum(deaths[1:j])
                     }
                  if(ii!=1)               
                     {
                     ttotvec[ii] = sum(ttot[1:j])-sum(ttotvec[1:(ii-1)])
                     deavec[ii]=sum(deaths[1:j])-sum(deavec[1:(ii-1)])
                     }
                  }
       #do the last item in ttotvec and deavec
       ttotvec[cutlen+1] = sum(ttot)-sum(ttotvec[1:cutlen])
       deavec[cutlen+1]  = sum(deaths)-sum(deavec[1:cutlen])
       #Compute ts(1) and pvalues(1)
       for (k in 1:length(allt)) 
         {
         ttot1 = ttotvec[k]
         ttot2 = ttotvec[k+1]
         dea1  = deavec[k]
         dea2  = deavec[k+1]
         returnval_3=exact_pvalue(ttot1,ttot2,dea1,dea2,mono)
         a2=returnval_3[1]
         pvalall[k]=returnval_3[2]
          }
       maxp = max(pvalall)
       maxplab = 0
      for (k1 in 1:length(pvalall))
        if (maxp == pvalall[k1])
           maxplab = k1
      #pvalall
      #save ts(i) and pvalues(i)
      ts[i] = allt[maxplab]
      pvalues[i]= pvalall[maxplab]
      allt=allt[-maxplab]
      pvalall=pvalall[-maxplab]
      #merge centain items on ttotvec and deavec
      ttotvec[maxplab] = ttotvec[maxplab]+ ttotvec[maxplab+1]
      ttotvec=ttotvec[-(maxplab+1)]
      deavec[maxplab]= deavec[maxplab]+ deavec[maxplab+1]
      deavec=deavec[-(maxplab+1)]
     }
    if (i!=1) 
     {
      #Use position maxplab, two vectors:(allt,pvalall), and ttotvec, deavec to compute ts(i) and pvalues(i)
      maxplab1=maxplab
      n=length(allt)
      if(maxplab1==1)
       {
         #p = lrtpvalue(ttotvec(1),ttotvec(2),deavec(1),deavec(2),n)
         returnval_4= exact_pvalue(ttotvec[1],ttotvec[2],deavec[1],deavec[2],mono)
         a2=returnval_4[1]
         p=returnval_4[2]
         if (length(deavec)==2)
           {
             pvalues[i]=p
             ts[i]=allt[maxplab]
            }
         else 
            {
             pvalall[1]= p
             maxp=max(pvalall)
             maxplab= 0
             for (k1 in 1:length(pvalall))
                 if (maxp == pvalall[k1])
                      maxplab= k1
              #save ts(i) and pvalues(i) 
              ts[i] = allt[maxplab]
              pvalues[i] = pvalall[maxplab]
              allt=allt[-maxplab]
              #pvalall
              pvalall=pvalall[-maxplab]
              #merge centain items on ttotvec and deavec
              ttotvec[maxplab] = ttotvec[maxplab] + ttotvec[maxplab+1]
              ttotvec=ttotvec[-(maxplab+1)]
              deavec[maxplab]  = deavec[maxplab] + deavec[maxplab+1]
              deavec=deavec[-(maxplab+1)]
              }
         }
       if((maxplab1!=1)&(maxplab1>n))
        {
          returnval_5=exact_pvalue(ttotvec[maxplab-1],ttotvec[maxplab],deavec[maxplab-1],deavec[maxplab],mono)
          a2=returnval_5[1]
          p=returnval_5[2]
          pvalall[length(allt)]= p
          maxp= max(pvalall)
          maxplab=0
          for (k1 in 1:length(pvalall))
            if (maxp == pvalall[k1])
                maxplab= k1
          #save ts(i) and pvalues(i)
          ts[i] = allt[maxplab]
          pvalues[i] = pvalall[maxplab]
          allt=allt[-maxplab]
          #pvalall
          pvalall=pvalall[-maxplab]
          #merge centain items on ttotvec and deavec
          ttotvec[maxplab] = ttotvec[maxplab]+ ttotvec[maxplab+1]
          ttotvec=ttotvec[-(maxplab+1)]
          deavec[maxplab] = deavec[maxplab] + deavec[maxplab+1]
          deavec=deavec[-(maxplab+1)]
          }
       if ((maxplab1!=1)&&(maxplab1<=n))
        {
        returnval_6=exact_pvalue(ttotvec[maxplab-1],ttotvec[maxplab],deavec[maxplab-1],deavec[maxplab],mono)
        a2=returnval_6[1]
        pfront=returnval_6[2]
        returnval_7=exact_pvalue(ttotvec[maxplab],ttotvec[maxplab+1],deavec[maxplab],deavec[maxplab+1],mono)
        a2=returnval_7[1]
        pback=returnval_7[2]                 
        pvalall[maxplab-1]= pfront
        pvalall[maxplab]= pback
        maxp= max(pvalall)
        maxplab= 0
        for (k1 in 1:length(pvalall))
           if (maxp==pvalall[k1])
              maxplab=k1
        #save ts(i) and pvalues(i)
        ts[i]= allt[maxplab]
        pvalues[i]= pvalall[maxplab]
        allt=allt[-maxplab]
        #pvalall
        pvalall=pvalall[-maxplab]
        #merge centain items on ttotvec and deavec
        ttotvec[maxplab]= ttotvec[maxplab]+ ttotvec[maxplab+1]
        ttotvec=ttotvec[-(maxplab+1)]
        deavec[maxplab]= deavec[maxplab]+ deavec[maxplab+1]
        deavec=deavec[-(maxplab+1)]
       }
     }
    }
   }
  returnv=cbind(ts,pvalues)
  return(returnv)
}
           
        
    
umbrella=function(time_die,ttot,deaths,indi)
{
#########################################################################
#Using the umbrella alternative to merge certain entries to make 
#    the sequence of ttot/deaths to increase then decrease 
#    or to decrease then increase 
#   (Note that the pava function makes it non decreasing. 
#    This function directly uses function pava().)
#
#Input: 
#       time_die == a sequence of times where deaths happened.
#       ttot     == the total time on test between each time point 
#                     and the previous time point (or 0).
#      deaths   == the number of deaths at each time point.
#       indi     == an indicator 
#                indi == 0: monotonic failure rate (either decrease ...
#                           or increase)
#                indi == 1: denoting the failure rate increase then
#                           decrease
#                indi == 2: denoting the failure rate decrease then
#                           increase
# Output:
#       time2    == the merged time_die
#       struct   == a structure saves the partition information
#       label    == a note about how the failure rate varies
#       indx     == the position where the change point value is
#########################################################################
n=length(time_die)
time_all = matrix(0,nrow=n+1,ncol=1)
time_all[1]=0
time_r=array(0,c(n,n))
ttot_r=array(0,c(n,n))
deaths_r=array(0,c(n,n))
loglik_r=array(0,c(1,n))


for (i in 2:(n+1))
   time_all[i]= time_die[i-1]
if (indi==0) # monotonic failure rate.
  {
   # in the structure, the first element is decreasing rate; the second element is increasing rate
   returnva_p1=pava_dfr(time_die,ttot,deaths)#returnva_p=list(struct$time[1],struct$ttot[1],struct$deaths[1])
   m=dim(returnva_p1)[2]/3
   w1=length(returnva_p1)/3
   time_r[1:w1,1]=returnva_p1[,1:m]
   ttot_r[1:w1,1]=returnva_p1[,(m+1):2*m]
   deaths_r[1:w1,1]=returnva_p1[,(2*m+1):3*m]
   returnva_p2=pava_ifr(time_die,ttot,deaths)#returnva_p2=list(struct$time[2],struct$ttot[2],struct$deaths[2])
   m=dim(returnva_p2)[2]/3
   w2=length(returnva_p2)/3
   time_r[1:w2,2]=returnva_p2[,1:m]
   ttot_r[1:w2,2]=returnva_p2[,(m+1):2*m]
   deaths_r[1:w2,2]=returnva_p2[,(2*m+1):3*m]
   #computing the likelihood
   for (i in 1:2)
     {
       t=time_r[,i]
       t=t[t!=0]
       tt=ttot_r[,i]
       tt=tt[tt!=0]
       d=deaths_r[1:length(tt),i]
     loglik_r[1,i]=gamllik(t,tt,d,time_die,ttot,deaths) 
     }
   #select the times using the larger loglikelihood
   loglik=loglik_r[1,1]
   time2=time_r[,1]
   i=2
   if (loglik<loglik_r[1,i])
     {
      time2  = time_r[,i]
      loglik = loglik_r[1,i]
      label  = "Decreasing failure rate."
      indx   = 1
      }
   else
      {
      label  ="Increasing failure rate."
      indx   = length(time2)
      }
  }
if (indi==1)#failure rate increases then decreases. ttot/death U shape 
  {
  for (j in 1:length(time_die))
    {
      if (j==1)
        {
        returnva_p3=pava_dfr(time_die,ttot,deaths)#returnva_p3=cbind(struct$time[j],struct$ttot[j],struct$deaths[j])
         m=dim(returnva_p3)[2]/3
         w1=length(returnva_p3)/3
         time_r[1:w1,j]=returnva_p3[,1:m]
         ttot_r[1:w1,j]=returnva_p3[,(m+1):2*m]
         deaths_r[1:w1,j]=returnva_p3[,(2*m+1):3*m]
        }
      if (j==length(time_die))
        {
         returnva_p4=pava_ifr(time_die,ttot,deaths)#returnva_p4=list(struct$time[j],struct$ttot[j],struct$deaths[j])
         m=dim(returnva_p4)[2]/3
         w2=length(returnva_p4)/3
         time_r[1:w2,j]=returnva_p4[,1:m]
         ttot_r[1:w2,j]=returnva_p4[,(m+1):2*m]
         deaths_r[1:w2,j]=returnva_p4[,(2*m+1):3*m]
         }
      if ((j!=1)&&(j<length(time_die))) 
         {
         returnva_p5=pava_ifr(time_die[1:j],ttot[1:j],deaths[1:j])
         m=dim(returnva_p5)[2]/3
         w3=length(returnva_p5)/3
         tmp_time1=returnva_p5[,1:m]
         tmp_ttot1=returnva_p5[,(m+1):2*m]
         tmp_deaths1=returnva_p5[,(2*m+1):3*m]
         returnva_p6=pava_dfr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
         m=dim(returnva_p6)[2]/3
         w4=length(returnva_p6)/3
         tmp_time2=returnva_p6[,1:m]
         tmp_ttot2=returnva_p6[,(m+1):2*m]
         tmp_deaths2=returnva_p6[,(2*m+1):3*m]
         time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
         ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
         deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)
         }
       #compute the loglikelihood for struct(j)
       t=time_r[,j]
       t=t[t!=0]
       tt=ttot_r[,j]
       tt=tt[tt!=0]
       d=deaths_r[1:length(tt),j]
       loglik_r[1,j]= gamllik(t,tt,d,time_die,ttot,deaths)
     }
   label="increasing then decreasing failure rate"
   #select the times using the largest loglikelihood
   loglik=loglik_r[1,1]
   time2=time_r[,1]
   indx    = 1
   for (j in 1:length(time_die))
      if (loglik < loglik_r[1,j])
         {
          time2  = time_r[,j]  
          loglik = loglik_r[1,j]
          indx   = j
          }
   }
if (indi==2)#failure rate decreases then increases. ttot/death arc shape
  {
   for (j in  1:length(time_die))
     {
     if (j==1)
       {
         returnva_p7=pava_ifr(time_die,ttot,deaths)#returnva_p4=list(struct$time[j],struct$ttot[j],struct$deaths[j])
         m=dim(returnva_p7)[2]/3
         w1=length(returnva_p7)/3
         time_r[1:w1,j]=returnva_p7[,1:m]
         ttot_r[1:w1,j]=returnva_p7[,(m+1):2*m]
         deaths_r[1:w1,j]=returnva_p7[,(2*m+1):3*m]
        }
     if (j == length(time_die))
        {
        returnva_p=pava_dfr(time_die,ttot,deaths)#returnva_p4=list(struct$time[j],struct$ttot[j],struct$deaths[j])
        m=dim(returnva_p)[2]/3
        w2=length(returnva_p)/3
        time_r[1:w2,j]=returnva_p[,1:m]
        ttot_r[1:w2,j]=returnva_p[,(m+1):2*m]
        deaths_r[1:w2,j]=returnva_p[,(2*m+1):3*m]
        }
     if ((j!=1)&&j < length(time_die))
        {

         returnva_p8=pava_dfr(time_die[1:j],ttot[1:j],deaths[1:j])
         m=dim(returnva_p8)[2]/3
         w3=length(returnva_p8)/3
         tmp_time1=returnva_p8[,1:m]
         tmp_ttot1=returnva_p8[,(m+1):2*m]
         tmp_deaths1=returnva_p8[,(2*m+1):3*m]
         returnva_p9=pava_ifr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
         m=dim(returnva_p9)[2]/3
         w4=length(returnva_p9)/3
         tmp_time2=returnva_p9[,1:m]
         tmp_ttot2=returnva_p9[,(m+1):2*m]
         tmp_deaths2=returnva_p9[,(2*m+1):3*m]
         time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
         ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
         deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)
         }
      # compute the loglikelihood for struct(j);
       t=time_r[,j]
       t=t[t!=0]
       tt=ttot_r[,j]
       tt=tt[tt!=0]
       d=deaths_r[1:length(tt),j]
       loglik_r[1,j]= gamllik(t,tt,d,time_die, ttot, deaths)
    }
    label = "decreasing then increasing failure rate"
    # select the times using the largest loglikelihood
    loglik  = loglik_r[1,1]
    time2   = time_r[,1]
    indx    = 1
    for (j in 1:length(time_die))
       if (loglik < loglik_r[1,j])
         {
          time2  = time_r[,j]   
          loglik = loglik_r[1,j]
          indx   = j
          }
  }
#  time2
#  loglik
#  indx

if (indi == 3)#failure rate increases then decreases.ttot/death U shape, peak point removed
  {
  for (j in 1:length(time_die))
    {
    if (j==1)
       {
        returnva_p7=pava_dfr(time_die,ttot,deaths)
        m=dim(returnva_p7)[2]/3
        w1=length(returnva_p7)
        time_r[1:w1,j]=returnva_p7[,1:m]
        ttot_r[1:w1,j]=returnva_p7[,(m+1):2*m]
        deaths_r[1:w1,j]=returnva_p7[,(2*m+1):3*m]
        }
     if (j == length(time_die))
        {
        returnva_p=pava_ifr(time_die,ttot,deaths)#returnva_p=list(struct$time[j],struct$ttot[j],struct$deaths[j])
        m=dim(returnva_p)[2]/3
        w2=length(returnva_p)
        time_r[1:w2,j]=returnva_p[,1:m]
        ttot_r[1:w2,j]=returnva_p[,(m+1):2*m]
        deaths_r[1:w2,j]=returnva_p[,(2*m+1):3*m]
        }
     if ((j!=1)&&(j < length(time_die)))
        {
         returnva_p8=pava_ifr(time_die[1:j],ttot[1:j],deaths[1:j])
         m=dim(returnva_p8)[2]/3
         w3=length(returnva_p8)/3
         tmp_time1=returnva_p8[,1:m]
         tmp_ttot1=returnva_p8[,(m+1):2*m]
         tmp_deaths1=returnva_p8[,(2*m+1):3*m]
         returnva_p9=pava_dfr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
         m=dim(returnva_p9)[2]/3
         w4=length(returnva_p9)/3
         tmp_time2=returnva_p9[,1:m]
         tmp_ttot2=returnva_p9[,(m+1):2*m]
         tmp_deaths2=returnva_p9[,(2*m+1):3*m]
         time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
         ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
         deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)       
         }
     # compute the loglikelihood for struct(j)
     t=time_r[,j]
     t=t[t!=0]
     tt=ttot_r[,j]
     tt=tt[tt!=0]
     d=deaths_r[1:length(tt),j]
     loglik_r[1,j]= gamllik(t,tt,d,time_die, ttot, deaths)
 
     if ((j < length(time_die))&& (j>1))
       {
         time_r[length(tmp_time1),j]=NA
         ttot_r[-length(tmp_time1),j]=NA
         deaths_r[length(tmp_time1),j]=NA
        }
   }
    label = "increasing then decreasing failure rate"
    # select the times using the largest loglikelihood
    loglik  = loglik_r[1,1]
    time2   = time_r[,1]
    indx    = 1
    for (j in 1:length(time_die))
       if (loglik < loglik_r[1,j])
          {
          time2  = time_r[,j]    
          loglik = loglik_r[1,j]
          indx   = j
          }
   }
if (indi == 4)# failure rate decreases then increases ttot/death arc shape; peak point removed.
  {
  for (j in 1:length(time_die))
    {
    if (j==1)
       {
        returnva_p7=pava_ifr(time_die,ttot,deaths)
        m=dim(returnva_p7)[2]/3
        w1=length(returnva_p7)/3
        time_r[1:w1,j]=returnva_p7[,1:m]
        ttot_r[1:w1,j]=returnva_p7[,(m+1):2*m]
        deaths_r[1:w1,j]=returnva_p7[,(2*m+1):3*m]
        }
     if (j == length(time_die))
        {
        returnva_p=pava_dfr(time_die,ttot,deaths)#returnva_p=list(struct$time[j],struct$ttot[j],struct$deaths[j])
        m=dim(returnva_p)[2]/3
        w2=length(returnva_p)/3
        time_r[1:w2,j]=returnva_p[,1:m]
        ttot_r[1:w2,j]=returnva_p[,(m+1):2*m]
        deaths_r[1:w2,j]=returnva_p[,(2*m+1):3*m]
        }
     if (j<length(time_die)&&j!=1)
        {
         returnva_p8=pava_dfr(time_die[1:j],ttot[1:j],deaths[1:j])
         m=dim(returnva_p8)[2]/3
         w3=length(returnva_p8)/3
         tmp_time1=returnva_p8[,1:m]
         tmp_ttot1=returnva_p8[,(m+1):2*m]
         tmp_deaths1=returnva_p8[,(2*m+1):3*m]
         returnva_p9=pava_ifr(time_die[-(1:j)],ttot[-(1:j)],deaths[-(1:j)])
         m=dim(returnva_p9)[2]/3
         w4=length(returnva_p9)/3
         tmp_time2=returnva_p9[,1:m]
         tmp_ttot2=returnva_p9[,(m+1):2*m]
         tmp_deaths2=returnva_p9[,(2*m+1):3*m]
         time_r[1:(w3+w4),j] = c(tmp_time1,tmp_time2)
         ttot_r[1:(w3+w4),j]= c(tmp_ttot1,tmp_ttot2)
         deaths_r[1:(w3+w4),j] = c(tmp_deaths1,tmp_deaths2)            
         }
     # compute the loglikelihood for struct(j)
     t=time_r[,j]
     t=t[t!=0]
     tt=ttot_r[,j]
     tt=tt[tt!=0]
     d=deaths_r[1:length(tt),j]
     loglik_r[1,j]= gamllik(t,tt,d,time_die, ttot, deaths)
     if (j < length(time_die) && j > 1)
       {
         time_r[length(tmp_time1),j]=NA
         ttot_r[-length(tmp_time1),j]=NA
         deaths_r[length(tmp_time1),j]=NA
        }
    }
    label = "decreasing then increasing failure rate"
    # select the times using the largest loglikelihood
    loglik  = loglik_r[1,1]
    time2   = time_r[,1]
    indx    = 1
    for (j in 1:length(time_die))
       if (loglik < loglik_r[1,j])
          {
           time2  = time_r[,j]
           loglik = loglik_r[1,j]
           indx   = j
           }
 }
#  loglik
#  indx
time_r=time_r[!is.na(time_r)&time_r!=0]
ttot_r=ttot_r[!is.na(ttot_r)&ttot_r!=0]
deaths_r=deaths_r[!is.na(deaths_r)&deaths_r!=0]
time2=time2[time2!=0&!is.na(time2)]
struct=list(time_r,ttot_r,deaths_r,loglik_r)
returnc=list(time2, struct, label, indx)
return(returnc)
}


loopcuts_umbrella=function(time,censor,cuttimes,mono)
{
####################################################
#using a loop format to find out the times to make the cuts from large pvalues to small pvalues and the lambdas in between
#
#Input: 
#      time: A sequence of time 
#      censor: and censor (0=censored and 1=not censored)
#      cuttimes: unique, sorted, possible times to make the cuts, including 0 and the ending time
#      mono: indicate the type of the test, a vector of the size of "cuttimes"
#      mono == 0: 2-sided hypothesis: H0:lam1=lam2; H1:lam1 \ne lam2
#           == 1: 1-sides hypothesis: H0:lam1>=lam2; H1:lam1 < lam2(decreasing failure rate constraint)
#           == 2: 1-sides hypothesis: H0:lam1<=lam2; H1:lam1 > lam2(increasing failure rate constraint)
#Output: 
#      ts: the times where the cuts shall be made 
#      pvalues: the p values for deleting each cutting times
######################################################

#Sort cuttimes and find cutlen
cuttimes=unique(sort(cuttimes))
cutlen=length(cuttimes)

#prepare the death time, ttot, and the number of deaths
returnval_1=totaltest(time,censor)
m=dim(returnval_1)[2]/3
time_die=returnval_1[,1:m]
ttot=returnval_1[,(m+1):(2*m)]
deaths=returnval_1[,(2*m+1):(3*m)]

#Based on cuttimes, find out the corresponding totaltime on test and the number of deaths.



#Set up ts and pvalues
ts=array(0,c(cutlen,1))
pvalues=array(0,c(cutlen,1))
if (cutlen==1)
 {
    for(j in 1:length(time_die))
      if (time_die[j]==cuttimes)
        {
        ttot1 = sum(ttot[1:j])
        ttot2 = sum(ttot[-(1:j)])
        death1 = sum(deaths[1:j])
        death2 = sum(deaths[-(1:j)])
        returnval_2= exact_pvalue(ttot1,ttot2,death1,death2,mono)
        a2=returnval_2[1]
        p=returnval_2[2]
         }
  ts= cuttimes
  pvalues = p
  }
if(cutlen!=1)
{
   for(i in 1:cutlen)
    {
     if(i==1)
      {
      ###initialize the matrix storing ttots and deaths p values
      ttotvec = array(0,c(cutlen+1,1))
      deavec  = array(0,c(cutlen+1,1))
      allt    = cuttimes
      pvalall = array(0,c(cutlen,1))
      monall  = mono
      for(ii in 1:cutlen)
        for(j in 1:length(time_die))
            if (time_die[j] == cuttimes[ii])
               {
                if (ii==1)
                    {
                     ttotvec[ii] = sum(ttot[1:j])
                     deavec[ii]  = sum(deaths[1:j])
                     }
                 if(ii!=1)              
                     {
                     ttotvec[ii] = sum(ttot[1:j]) - sum(ttotvec[1:(ii-1)])
                     deavec[ii]  = sum(deaths[1:j]) - sum(deavec[1:(ii-1)])
                     }
                 }
       #do the last item in ttotvec and deavec
       ttotvec[cutlen+1] = sum(ttot)-sum(ttotvec[1:cutlen])
       deavec[cutlen+1]  = sum(deaths)-sum(deavec[1:cutlen])
       #Compute ts(1) and pvalues(1)
       for (k in 1:length(allt)) 
         {
         ttot1 = ttotvec[k]
         ttot2 = ttotvec[k+1]
         dea1  = deavec[k]
         dea2  = deavec[k+1]
         returnval_3=exact_pvalue(ttot1,ttot2,dea1,dea2,monall[k])
         a2=returnval_3[1]
         pvalall[k]=returnval_3[2]
          }
       maxp = max(pvalall)
       maxplab = 0
       for (k1 in 1:length(pvalall))
        if (maxp == pvalall[k1])
           maxplab = k1
      #pvalall
      #save ts(i) and pvalues(i)
      ts[i] = allt[maxplab]
      pvalues[i]= pvalall[maxplab]
      allt= allt[-maxplab]
      pvalall=pvalall[-maxplab]
      monall=monall[-maxplab]
      #merge centain items on ttotvec and deavec
      ttotvec[maxplab] = ttotvec[maxplab]+ ttotvec[maxplab+1]
      ttotvec=ttotvec[-(maxplab+1)]
      deavec[maxplab]= deavec[maxplab]+ deavec[maxplab+1]
      deavec=deavec[-(maxplab+1)]
    }
 if (i!=1)
   {
     #Use position maxplab, two vectors:(allt,pvalall), and ttotvec,deavec to compute ts(i) and pvalues
     maxplab1=maxplab
     n=length(allt)
    if(maxplab1==1)
      {
        #p = lrtpvalue(ttotvec(1),ttotvec(2),deavec(1),deavec(2),n)
        returnval_4= exact_pvalue(ttotvec[1],ttotvec[2],deavec[1],deavec[2],monall[1])
        a2=returnval_4[1]
        p=returnval_4[2]
        if (length(deavec)== 2)
           {
             pvalues[i] = p
             ts[i] = allt[maxplab]
            }
         if (length(deavec)!=2)
            {
             pvalall[1]= p
             maxp = max(pvalall)
             maxplab= 0
             for (k1 in 1:length(pvalall))
                 if (maxp == pvalall[k1])
                      maxplab= k1
              #save ts(i) and pvalues(i) 
              ts[i] = allt[maxplab]
              pvalues[i] = pvalall[maxplab]
              allt=allt[-maxplab]
              #pvalall
              pvalall=pvalall[-maxplab]
              #monall
              monall=monall[-maxplab]
              #merge centain items on ttotvec and deavec
              ttotvec[maxplab] = ttotvec[maxplab] + ttotvec[maxplab+1]
              ttotvec=ttotvec[-(maxplab+1)]
              deavec[maxplab]  = deavec[maxplab] + deavec[maxplab+1]
              deavec=deavec[-(maxplab+1)]
             }
         }
 
      if((maxplab1!=1)&(maxplab1>n))
        {
         returnval_5=exact_pvalue(ttotvec[maxplab-1],ttotvec[maxplab],deavec[maxplab-1],deavec[maxplab],monall[maxplab-1])
         a2=returnval_5[1]
         p=returnval_5[2]
         pvalall[length(allt)]= p
         maxp= max(pvalall)
         maxplab= 0
         for (k1 in 1:length(pvalall))
            if (maxp == pvalall[k1])
                maxplab= k1
         #save ts(i) and pvalues(i)
         ts[i] = allt[maxplab]
         pvalues[i] = pvalall[maxplab]
         allt=allt[-maxplab]
         #pvalall
         pvalall=pvalall[-maxplab]
         #monall
         monall=monall[-maxplab]
         #merge centain items on ttotvec and deavec
          ttotvec[maxplab] = ttotvec[maxplab]+ ttotvec[maxplab+1]
          ttotvec=ttotvec[-(maxplab+1)]
          deavec[maxplab] = deavec[maxplab] + deavec[maxplab+1]
          deavec=deavec[-(maxplab+1)]
          }
   
       if((maxplab1!=1)&(maxplab1<=n))
        {
        returnval_6=exact_pvalue(ttotvec[maxplab-1],ttotvec[maxplab],deavec[maxplab-1],deavec[maxplab],monall[maxplab-1])
        a2=returnval_6[1]
        pfront=returnval_6[2]
        returnval_7=exact_pvalue(ttotvec[maxplab],ttotvec[maxplab+1],deavec[maxplab],deavec[maxplab+1],monall[maxplab])
        a2=returnval_7[1]
        pback=returnval_7[2]                 
        pvalall[maxplab-1]= pfront
        pvalall[maxplab]= pback
        maxp= max(pvalall)
        maxplab= 0
        for (k1 in 1:length(pvalall))
           if (maxp == pvalall[k1])
              maxplab= k1
        #save ts(i) and pvalues(i)
        ts[i]= allt[maxplab]
        pvalues[i]= pvalall[maxplab]
        allt=allt[-maxplab]
        #pvalall
        pvalall=pvalall[-maxplab]
        #monall
        monall=monall[-maxplab]
        #merge centain items on ttotvec and deavec
        ttotvec[maxplab]= ttotvec[maxplab]+ ttotvec[maxplab+1]
        ttotvec=ttotvec[-(maxplab+1)]
        deavec[maxplab]=deavec[maxplab]+ deavec[maxplab+1]
        deavec=deavec[-(maxplab+1)]
        }
     }
   }
 }
  returnv=cbind(ts,pvalues)
  return(returnv)
}
           


RPEXEv1=function(EventTime=NA,eventtime,Censor=NA,censor,CutTimes=NA,cuttime=1,Trend=NA,trend=0)
{
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

   
