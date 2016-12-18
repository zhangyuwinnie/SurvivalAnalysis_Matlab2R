function [time_die,ttot,deaths] = totaltest(time,censor);
% Function 'totaltest' computes total time on test
% Inputs: times and censor (0 = censored; 1 = uncensored)
% Return: time_die: times that events occur (in ascending order)
%         ttot: total time on test at each time point in time_die
%         deaths: number of death at each time point in time_die

% save the data in the orginal structure
seq     = zeros(length(time),1);
for i = 1:length(time),
    seq(i) = i;
end
tmpdata = [seq censor time];

%%% sort the tmp data
tmp2     = zeros(size(tmpdata));
for i = 1:size(tmp2,1),
    tmp2(i,3) = min(tmpdata(:,3));
    j = 1;
    while  tmpdata(j,3) ~= min(tmpdata(:,3)),
        j= j+1;
    end
    tmp2(i,1:2) = tmpdata(j,1:2);
    tmpdata(j,:) = [];
end

%%% Compute alpha's for the sequence
for i = 1:size(tmp2,1),
    if tmp2(i,2) == 0,
       tmp2(i,4) = 0;
    else
       tmp2(i,4) = 1;
    end
end

%%% Deal with alpha > 1
for i = 1:(size(tmp2,1)-1),
    if tmp2(size(tmp2,1)+1-i,3)== tmp2(size(tmp2,1)-i,3),
        tmp2(size(tmp2,1)-i,4) = tmp2(size(tmp2,1)-i,4) + tmp2(size(tmp2,1)+1-i,4);
        tmp2(size(tmp2,1)+1-i,4) = 0;
    end
end

%%% Delete the repeats
k = 0;
for i = 1:size(tmp2,1),
    if tmp2(i,2) == 1 && tmp2(i,4)==0,
        k(length(k)+1)=i;
    end
end
k(1)=[];
tmp3 = tmp2;
tmp3(k,:) = [];


%%% Compute the number of patients in the study
for i = 1:size(tmp3,1),
    if tmp3(i,4)==0,
       tmp3(i,5)= 1;
    else
       tmp3(i,5)= tmp3(i,4);
    end
end
for i = 1:size(tmp3,1)-1,
    tmp3(i,6)   = sum(tmp3(:,5))-sum(tmp3(1:i,5));
    tmp3(size(tmp3,1),6)   = 0;
end

%%% Compute the survival time of this cohort
for i = 1:size(tmp3,1),
    if i ==1,
       tmp3(i,7) = sum(tmp3(:,5))*tmp3(i,3);
    else
        %%% survival time == patient number * time difference
       tmp3(i,7) = tmp3(i-1,6) * (tmp3(i,3)-tmp3(i-1,3));
    end
end

tmp3(:,8) = tmp3(:,7);
for i = 1:size(tmp3,1),
    if tmp3(i,2)==0,
        if tmp3(i:end,2)'*tmp3(i:end,2)>0,
            tmp3(i+1,8) = tmp3(i,8)+tmp3(i+1,8);
        elseif  tmp3(i:end,2)'*tmp3(i:end,2)==0 && tmp3(i-1,2)~=0;
            %%% put all the credit to the last noncensered data
            k = length(tmp3(i:end,2));
            for j = 1:k,
                tmp3(i-1,8) = tmp3(i-1,8)+tmp3(i-1+j,8);
            end
        end
    end
end

%%% Build the survival reaction
tmp3(:,9) = tmp3(:,6);
for i=2:length(tmp3(:,9)),
    tmp3(i,9) = tmp3(i-1,9)-tmp3(i,4);
end

% plot(tmp3(:,3),tmp3(:,9))
        

%%% delete all the censered items
k = 0;
for i = 1:length(tmp3(:,1)),
    if tmp3(i,2) == 0,
        k(length(k)+1)=i;
    end
end
k(1)=[];
tmp4 = tmp3;
tmp4(k,:) = [];

time_die    =  tmp4(:,3);
ttot        =  tmp4(:,8);
deaths      =  tmp4(:,5);


