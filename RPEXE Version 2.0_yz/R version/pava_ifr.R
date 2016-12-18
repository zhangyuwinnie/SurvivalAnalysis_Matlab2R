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

pava_ifr <- function(time_die,ttot,deaths)
{
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
