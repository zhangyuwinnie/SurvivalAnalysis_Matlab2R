function [time2, ttot2, deaths2] = pava_ifr(time_die,ttot,deaths);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% merge certain entries to make the sequence of ttot/deaths to be 
%%% non increasing 
%%% (Note that the pava function makes it non decreasing. 
%%%  This function directly uses function pava().)
%%%
%%% Input: 
%%%       time_die == a sequence of times where deaths happened.
%%%       ttot     == the total time on test at each time point.
%%%       deaths   == the number of deaths at each time point.
%%% Output:
%%%       time2    == the merged time_die
%%%       ttot2    == .......... ttot
%%%       deaths2  == .......... deaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ttotrev                 = -ttot;
[time2, ttot3, deaths2] = pava_dfr(time_die,ttotrev,deaths);
ttot2                   = ttot3*(-1);