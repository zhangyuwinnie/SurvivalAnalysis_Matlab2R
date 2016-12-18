function [time2, struct, label, indx] = umbrella(time_die,ttot,deaths,indi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Using the umbrella alternative to merge certain entries to make 
%%%    the sequence of ttot/deaths to increase then decrease 
%%%    or to decrease then increase 
%%%   (Note that the pava function makes it non decreasing. 
%%%    This function directly uses function pava().)
%%%
%%% Input: 
%%%       time_die == a sequence of times where deaths happened.
%%%       ttot     == the total time on test between each time point 
%%%                     and the previous time point (or 0).
%%%       deaths   == the number of deaths at each time point.
%%%       indi     == an indicator 
%%%                indi == 0: monotonic failure rate (either decrease ...
%%%                           or increase)
%%%                indi == 1: denoting the failure rate increase then
%%%                           decrease
%%%                indi == 2: denoting the failure rate decrease then
%%%                           increase
%%% Output:
%%%       time2    == the merged time_die
%%%       struct   == a structure saves the partition information
%%%       label    == a note about how the failure rate varies
%%%       indx     == the position where the change point value is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = length(time_die);
time_all = zeros(n+1,1);
time_all(1) = 0;
for i = 2:n+1;
    time_all(i) = time_die(i-1);
end


if indi == 0, % monotonic failure rate.
    % in the structure, the first element is decreasing rate; 
    %                   the second element is increasing rate.
    [struct(1).time, struct(1).ttot, struct(1).deaths] = ...
        pava_dfr(time_die,ttot,deaths);
    [struct(2).time, struct(2).ttot, struct(2).deaths] = ...
        pava_ifr(time_die,ttot,deaths);
    % computing the likelihood
    for i = 1:2, 
        struct(i).loglik = gamllik(struct(i).time, struct(i).ttot, ...
            struct(i).deaths, time_die, ttot, deaths);
    end 
    
    % select the times using the larger loglikelihood
    loglik  = struct(1).loglik;
    time2   = struct(1).time;
    i       = 2;
    if loglik < struct(i).loglik,
        time2  = struct(i).time;
        loglik = struct(i).loglik;
        label  = 'Decreasing failure rate.';
        indx   = 1;
    else
        label  = 'Increasing failure rate.';
        indx   = length(time2);
    end   
elseif indi == 1, % failure rate increases then decreases. ttot/death U shape 
    for j = 1:length(time_die),
        if j == 1,
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_dfr(time_die,ttot,deaths); 
        elseif j == length(time_die),
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_ifr(time_die,ttot,deaths);    
        elseif j < length(time_die),
            [tmp_time1, tmp_ttot1, tmp_deaths1] = ...
                pava_ifr(time_die(1:j),ttot(1:j),deaths(1:j));
            [tmp_time2, tmp_ttot2, tmp_deaths2] = ...
                pava_dfr(time_die(1+j:end),ttot(1+j:end),deaths(1+j:end));
            struct(j).time = [tmp_time1;tmp_time2];
            struct(j).ttot = [tmp_ttot1;tmp_ttot2];
            struct(j).deaths = [tmp_deaths1;tmp_deaths2];
        end
        % compute the loglikelihood for struct(j);
        struct(j).loglik = gamllik(struct(j).time, struct(j).ttot, ...
            struct(j).deaths, time_die, ttot, deaths);
    end
    label = 'increasing then decreasing failure rate';
    % select the times using the largest loglikelihood
    loglik  = struct(1).loglik;
    time2   = struct(1).time;
    indx    = 1;
    for j = 1:length(time_die),
        if loglik < struct(j).loglik,
            time2  = struct(j).time;
            loglik = struct(j).loglik;
            indx   = j;
        end  
    end
elseif indi == 2; % failure rate decreases then increases. ttot/death arc shape
    for j = 1:length(time_die),
        if j == 1,
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_ifr(time_die,ttot,deaths); 
        elseif j == length(time_die),
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_dfr(time_die,ttot,deaths);    
        elseif j < length(time_die),
            [tmp_time1, tmp_ttot1, tmp_deaths1] = ...
                pava_dfr(time_die(1:j),ttot(1:j),deaths(1:j));
            [tmp_time2, tmp_ttot2, tmp_deaths2] = ...
                pava_ifr(time_die(1+j:end),ttot(1+j:end),deaths(1+j:end));
            struct(j).time = [tmp_time1;tmp_time2];
            struct(j).ttot = [tmp_ttot1;tmp_ttot2];
            struct(j).deaths = [tmp_deaths1;tmp_deaths2];
        end
        % compute the loglikelihood for struct(j);
        struct(j).loglik = gamllik(struct(j).time, struct(j).ttot, ...
            struct(j).deaths, time_die, ttot, deaths);
    end
    label = 'decreasing then increasing failure rate';
    % select the times using the largest loglikelihood
    loglik  = struct(1).loglik;
    time2   = struct(1).time;
    indx    = 1;
    for j = 1:length(time_die),
        if loglik < struct(j).loglik,
            time2  = struct(j).time;
            loglik = struct(j).loglik;
            indx   = j;
        end  
    end   
%     time2
%     loglik
%     indx
    
elseif indi == 3, % failure rate increases then decreases. 
                  % ttot/death U shape, peak point removed
    for j = 1:length(time_die),
        if j == 1,
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_dfr(time_die,ttot,deaths); 
        elseif j == length(time_die),
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_ifr(time_die,ttot,deaths);    
        elseif j < length(time_die),
            [tmp_time1, tmp_ttot1, tmp_deaths1] = ...
                pava_ifr(time_die(1:j),ttot(1:j),deaths(1:j));
            [tmp_time2, tmp_ttot2, tmp_deaths2] = ...
                pava_dfr(time_die(1+j:end),ttot(1+j:end),deaths(1+j:end));
            struct(j).time = [tmp_time1;tmp_time2];
            struct(j).ttot = [tmp_ttot1;tmp_ttot2];
            struct(j).deaths = [tmp_deaths1;tmp_deaths2];
%             struct(j).time(length(tmp_time1)) = [];
%             struct(j).ttot(length(tmp_time1)) = [];
%             struct(j).deaths(length(tmp_time1)) = [];            
        end
        % compute the loglikelihood for struct(j);
        struct(j).loglik = gamllik(struct(j).time, struct(j).ttot, ...
            struct(j).deaths, time_die, ttot, deaths);
        if j < length(time_die) && j > 1,
            struct(j).time(length(tmp_time1)) = [];
            struct(j).ttot(length(tmp_time1)) = [];
            struct(j).deaths(length(tmp_time1)) = [];              
        end
    end
    label = 'increasing then decreasing failure rate';
    % select the times using the largest loglikelihood
    loglik  = struct(1).loglik;
    time2   = struct(1).time;
    indx    = 1;
    for j = 1:length(time_die),
        if loglik < struct(j).loglik,
            time2  = struct(j).time;
            loglik = struct(j).loglik;
            indx   = j;
        end  
    end     
elseif indi == 4; % failure rate decreases then increases
                  % ttot/death arc shape; peak point removed.
    for j = 1:length(time_die),
        if j == 1,
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_ifr(time_die,ttot,deaths); 
        elseif j == length(time_die),
            [struct(j).time, struct(j).ttot, struct(j).deaths] = ...
                pava_dfr(time_die,ttot,deaths);    
        elseif j < length(time_die),
            [tmp_time1, tmp_ttot1, tmp_deaths1] = ...
                pava_dfr(time_die(1:j),ttot(1:j),deaths(1:j));
            [tmp_time2, tmp_ttot2, tmp_deaths2] = ...
                pava_ifr(time_die(1+j:end),ttot(1+j:end),deaths(1+j:end));
            struct(j).time = [tmp_time1;tmp_time2];
            struct(j).ttot = [tmp_ttot1;tmp_ttot2];
            struct(j).deaths = [tmp_deaths1;tmp_deaths2];
%             struct(j).time(length(tmp_time1)) = [];
%             struct(j).ttot(length(tmp_time1)) = [];
%             struct(j).deaths(length(tmp_time1)) = [];              
        end
        % compute the loglikelihood for struct(j);
        struct(j).loglik = gamllik(struct(j).time, struct(j).ttot, ...
            struct(j).deaths, time_die, ttot, deaths);
        if j < length(time_die) && j > 1,
            struct(j).time(length(tmp_time1)) = [];
            struct(j).ttot(length(tmp_time1)) = [];
            struct(j).deaths(length(tmp_time1)) = [];              
        end        
    end
    label = 'decreasing then increasing failure rate';
    % select the times using the largest loglikelihood
    loglik  = struct(1).loglik;
    time2   = struct(1).time;
    indx    = 1;
    for j = 1:length(time_die),
        if loglik < struct(j).loglik,
            time2  = struct(j).time;
            loglik = struct(j).loglik;
            indx   = j;
        end  
    end   
%     time2
%     loglik
%     indx

end 

       
    


    
    
    
    
% Can be extended to umbrella alternative approach.    
%     
% if indi == 1, % the failure rate is increasing then decreasing.
%     
% % make the cuts and save the cuts in a structure
% for i = 1:length(time_all),
%     if i == 1, % then it is pava with decreasing failure rate
%         [struct(i).time, struct(i).ttot, struct(i).deaths] = ...
%             pava(time_die,ttot,deaths);
%     elseif i == length(time_all),
%         [struct(i).time, struct(i).ttot, struct(i).deaths] = ...
%             pava_ifr(time_die,ttot,deaths);
%     elseif i == 2,
%         struct(i).time1 = time_die(i-1);
%         struct(i).ttot1 = ttot(i-1);
%         struct(i).deaths1 = deaths(i-1);
%         [struct(i).time2 struct(i).ttot2 struct(i).deaths2] = ...
%             pava(time_die(i:end),ttot(i:end),deaths(i:end));
%         struct(i).time  = [struct(i).time1; struct(i).time2];  
%         struct(i).ttot  = [struct(i).ttot1; struct(i).ttot2];
%         struct(i).deaths= [struct(i).deaths1; struct(i).deaths2];
%     elseif i == length(time_all)-1,
%         [struct(i).time1 struct(i).ttot1 struct(i).deaths1] = ...
%             pava_ifr(time_die(1:i-1),ttot(1:i-1),deaths(1:i-1));
%         struct(i).time2 = time_die(end);
%         struct(i).ttot2 = ttot(end);
%         struct(i).deaths2 = deaths(end);
%         struct(i).time  = [struct(i).time1; struct(i).time2];  
%         struct(i).ttot  = [struct(i).ttot1; struct(i).ttot2];
%         struct(i).deaths= [struct(i).deaths1; struct(i).deaths2];
%     else
%         [struct(i).time1 struct(i).ttot1 struct(i).deaths1] = ...
%             pava_ifr(time_die(1:i-1),ttot(1:i-1),deaths(1:i-1));
%         [struct(i).time2 struct(i).ttot2 struct(i).deaths2] = ...
%             pava(time_die(i:end),ttot(i:end),deaths(i:end));
%         struct(i).time  = [struct(i).time1; struct(i).time2];  
%         struct(i).ttot  = [struct(i).ttot1; struct(i).ttot2];
%         struct(i).deaths= [struct(i).deaths1; struct(i).deaths2];
%     end
%     % compute the Gamma parameter
%     for j = 1: length(struct(i).time),
%         struct(i).gamindi(j) = struct(i).ttot(j)/struct(i).deaths(j);
%     end
%     % compute the likelihood
%     % get the indices of the cut time:
%     struct(i).indi = zeros(length(struct(i).time),1);
%     for j = 1:length(struct(i).time),
%         for jj = 1:length(time_die),
%             if time_die(jj) == struct(i).time(j);
%                 struct(i).indi(j) = jj;
%             end
%         end
%     end
%     % set the scale parameter for Gamma distribution of the ttot
%     struct(i).gampar = zeros(length(time_die),1);
%     for ii = 1:length(time_die),
%         for j = 1:length(struct(i).time),
%             if ii <= struct(i).indi(j),
%                 if j ==1, 
%                     struct(i).gampar(ii) = struct(i).gamindi(j);
%                 else 
%                     if ii > struct(i).indi(j-1),
%                         struct(i).gampar(ii) = struct(i).gamindi(j);
%                     end
%                 end
%             end
%         end
%     end
%     loglik = 0;
%     for ii = 1 : length(struct(i).gampar),
%         loglik = loglik + log(gampdf(ttot(ii)/deaths(ii),deaths(ii),...
%             struct(i).gampar(ii)));
%     end  
%     struct(i).loglik = loglik;
% end
% 
% 
% end
% 
% 
% 
% 
% if indi == 2, % the failure rate is decreasing then increasing.
% 
% % make the cuts and save the cuts in a structure
% for i = 1:length(time_all),
%     if i == 1, % then it is pava with increasing failure rate
%         [struct(i).time, struct(i).ttot, struct(i).deaths] = ...
%             pava_ifr(time_die,ttot,deaths);
%     elseif i == length(time_all),
%         [struct(i).time, struct(i).ttot, struct(i).deaths] = ...
%             pava(time_die,ttot,deaths);
%     elseif i == 2,
%         struct(i).time1 = time_die(i-1);
%         struct(i).ttot1 = ttot(i-1);
%         struct(i).deaths1 = deaths(i-1);
%         [struct(i).time2 struct(i).ttot2 struct(i).deaths2] = ...
%             pava_ifr(time_die(i:end),ttot(i:end),deaths(i:end));
%         struct(i).time  = [struct(i).time1; struct(i).time2];  
%         struct(i).ttot  = [struct(i).ttot1; struct(i).ttot2];
%         struct(i).deaths= [struct(i).deaths1; struct(i).deaths2];
%     elseif i == length(time_all)-1,
%         [struct(i).time1 struct(i).ttot1 struct(i).deaths1] = ...
%             pava(time_die(1:i-1),ttot(1:i-1),deaths(1:i-1));
%         struct(i).time2 = time_die(end);
%         struct(i).ttot2 = ttot(end);
%         struct(i).deaths2 = deaths(end);
%         struct(i).time  = [struct(i).time1; struct(i).time2];  
%         struct(i).ttot  = [struct(i).ttot1; struct(i).ttot2];
%         struct(i).deaths= [struct(i).deaths1; struct(i).deaths2];
%     else
%         [struct(i).time1 struct(i).ttot1 struct(i).deaths1] = ...
%             pava(time_die(1:i-1),ttot(1:i-1),deaths(1:i-1));
%         [struct(i).time2 struct(i).ttot2 struct(i).deaths2] = ...
%             pava_ifr(time_die(i:end),ttot(i:end),deaths(i:end));
%         struct(i).time  = [struct(i).time1; struct(i).time2];  
%         struct(i).ttot  = [struct(i).ttot1; struct(i).ttot2];
%         struct(i).deaths= [struct(i).deaths1; struct(i).deaths2];
%     end
%     % compute the Gamma parameter
%     for j = 1: length(struct(i).time),
%         struct(i).gamindi(j) = struct(i).ttot(j)/struct(i).deaths(j);
%     end
%     % compute the likelihood
%     % get the indices of the cut time:
%     struct(i).indi = zeros(length(struct(i).time),1);
%     for j = 1:length(struct(i).time),
%         for jj = 1:length(time_die),
%             if time_die(jj) == struct(i).time(j);
%                 struct(i).indi(j) = jj;
%             end
%         end
%     end
%     % set the scale parameter for Gamma distribution of the ttot
%     struct(i).gampar = zeros(length(time_die),1);
%     for ii = 1:length(time_die),
%         for j = 1:length(struct(i).time),
%             if ii <= struct(i).indi(j),
%                 if j ==1, 
%                     struct(i).gampar(ii) = struct(i).gamindi(j);
%                 else 
%                     if ii > struct(i).indi(j-1),
%                         struct(i).gampar(ii) = struct(i).gamindi(j);
%                     end
%                 end
%             end
%         end
%     end
%     loglik = 0;
%     for ii = 1 : length(struct(i).gampar),
%         loglik = loglik + log(gampdf(ttot(ii)/deaths(ii),deaths(ii),...
%             struct(i).gampar(ii)));
%     end  
%     struct(i).loglik = loglik;
% end
% % select the times using the smallest loglikelihood
% loglik = struct(1).loglik;
% time2  = struct(1).time;
% for i = 2:length(time_all),
%     if loglik < struct(i).loglik,
%         time2  = struct(i).time;
%         loglik = struct(i).loglik;
%     end    
% end
% 
% % remove the last time point from the cuts

