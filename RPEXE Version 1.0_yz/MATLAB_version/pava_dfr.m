function [time2, ttot2, deaths2] = pava_dfr(time_die,ttot,deaths)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% merge certain entries to make the sequence of ttot to be non decreasing
%%% Input: 
%%%       time_die == a sequence of times where deaths happened.
%%%       ttot     == the total time on test at each time point.
%%%       deaths   == the number of deaths at each time point.
%%% Output:
%%%       time2    == the merged time_die
%%%       ttot2    == .......... ttot
%%%       deaths2  == .......... deaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len  = length(ttot);
if len == 1,
    time2   = time_die;
    ttot2   = ttot;
    deaths2 = deaths;
else
    for j = 1:len,
        %%% indicate deleted items 
        for i = 1:length(ttot)-1,
            %%% indicate items to be deleted
            if ttot(i)/deaths(i) > ttot(i+1)/deaths(i+1),
                deaths(i+1) = deaths(i)+deaths(i+1);
                deaths(i)   = 0;
                ttot(i+1)= ttot(i) + ttot(i+1);
                ttot(i)  = 0;
            end
        end

        %%% delete the indicated items
        k = 0;
        for i = 1:length(ttot)-1,
            if ttot(i) == 0,
                k(length(k)+1)=i;
            end
        end
        k(1)        = [];
        ttot(k)     = [];
        deaths(k)   = [];
        time_die(k) = [];
    end

    time2   = time_die;
    ttot2   = ttot;
    deaths2 = deaths;
   
end




