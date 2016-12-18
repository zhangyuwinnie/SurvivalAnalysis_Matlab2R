%% Breast cancer example;
% Load Christo's data, keep the time, censor, ER status, 
% and predicted treatment sensitivity;
% Notes fron Christos:

% Hi Gang,
% 
% These are the datasets used in our JAMA paper.  You can use the development dataset (Dev) for training and the validation (Val) for testing.
% The relevant columns are drfs (indicator for distant relapse events) and drfs.time (event time in years).
% There are several covariates that can be used to stratify the cohorts - er.status, her2.status, nodal status (pre.N), response group (pCR.RD).
% 
% Let me know if you have any questions.
% 
% Thanks.
% Christos

% For grade: exclude NA, indeterminate or 0 and 
% categorize the rest as 3 vs 1 or 2.
% 
% T-stage (tumor size) as T3/T4 vs T1/T2.  
% T-stage had a stronger effect in multivariate Cox regression, 
% so it might be a better stratifier.
% 
% You could also use either the pCR/RD or 
% the predicted pCR/RD within each ER subset.  Both would be relevant.
% 



% % significance cut off for monotonic order assumption at alpha = 0.05: 
% N   cutoff
% 20  .0048
% 50  .0029
% 200 .0017
% 800 .0011
% % significance cut off for monotonic order assumption at alpha = 0.1:
% N   cutoff
% 20  .0108
% 50  .0067
% 200 .0038
% 800 .0026


addpath('..\RPEXE Version 1.0 (03032014)\MATLAB version');
addpath('..\data');

[bstnum, bsttxt, alldata] = xlsread('JAMABreast.xlsx');
% Three columns in the number: indicator of the validation, censor
%   indicator, time to drfs; 
% Two columns in the text: er status ("P" or "N"), 
%   chemo prediction: ("Rx Sensitive" or "Rx Insensitive");
vali_indi = bstnum(:,1);
bstcens   = bstnum(:,2);
bsttime   = bstnum(:,3);
ertxt     = bsttxt(2:end,4);
chemotxt  = bsttxt(2:end,5);
nodaltxt  = bsttxt(2:end,6);
resptxt   = bsttxt(2:end,7);
groupnum  = bstnum(:,8);
t_stage   = bsttxt(2:end,9);
predresptxt= bsttxt(2:end,10); 
% Get the indices for er, nodal, group, tumor size,
%    pathological response, predicted pathological response(dlda 30 predictor), 
er_p_indi   = strcmp(ertxt,'P');
er_n_indi   = strcmp(ertxt,'N');
chem_indi = strcmp(chemotxt,'Rx Sensitive');
% chem_s_indi = strcmp(chemotxt,'Rx Sensitive');
% chem_i_indi = strcmp(chemotxt,'Rx Insensitive');
nodal_p_indi= [strcmp(nodaltxt,'N1'); strcmp(nodaltxt,'N2'); ...
    strcmp(nodaltxt,'N3')];
nodal_n_indi= strcmp(nodaltxt,'N0');
group_p_indi= (groupnum==3);
group_n_indi= [(groupnum==1); (groupnum==2)];
tsize_n_indi= [strcmp(t_stage,'T1'); strcmp(t_stage,'T2')];
tsize_p_indi= [strcmp(t_stage,'T3'); strcmp(t_stage,'T4')];
res_p_indi  = strcmp(resptxt,'pCR');
res_n_indi  = strcmp(resptxt,'RD');
pres_p_indi = strcmp(predresptxt,'pCR');
pres_n_indi = strcmp(predresptxt,'RD');

% calculate the sample size
tabulate(vali_indi)
tabulate(er_p_indi(find(vali_indi==0)))
tabulate(er_p_indi(find(vali_indi==1)))
% tabulate(nodal_p_indi(find(vali_indi==0)))
% tabulate(nodal_p_indi(find(vali_indi==1)))
% tabulate(group_p_indi(find(vali_indi==0)))
% tabulate(group_p_indi(find(vali_indi==1)))
% tabulate(tsize_p_indi(find(vali_indi==0)))
% tabulate(tsize_p_indi(find(vali_indi==1)))


res_tab= crosstab(resptxt, predresptxt)
tabulate(resptxt)
tabulate(predresptxt)
%zy: put the following into comment
%tabulate([resptxt predresptxt])

%           pred pCR     pred RD   Total
%true pCR   67           32        99          
%true RD    123          266       389
%Missing    6            14        20
%Total      196          312
for i = 1:508,
    rest{i} = [resptxt{i} predresptxt{i}];
end
tabulate(rest)

%% Create all index
er_P_train = intersect(find(vali_indi==0),find(er_p_indi));
er_N_train = intersect(find(vali_indi==0),find(er_n_indi));
er_P_valid = intersect(find(vali_indi==1),find(er_p_indi));
er_N_valid = intersect(find(vali_indi==1),find(er_n_indi));
% nodal
er_P_train_nodal_p = intersect(er_P_train, find(nodal_p_indi));
er_P_train_nodal_n = intersect(er_P_train, find(nodal_n_indi));
er_N_train_nodal_p = intersect(er_N_train, find(nodal_p_indi));
er_N_train_nodal_n = intersect(er_N_train, find(nodal_n_indi));
er_P_valid_nodal_p = intersect(er_P_valid, find(nodal_p_indi));
er_P_valid_nodal_n = intersect(er_P_valid, find(nodal_n_indi));
er_N_valid_nodal_p = intersect(er_N_valid, find(nodal_p_indi));
er_N_valid_nodal_n = intersect(er_N_valid, find(nodal_n_indi));
length(er_P_train_nodal_p) + length(er_N_train_nodal_p)
length(er_P_train_nodal_n) + length(er_N_train_nodal_n)
length(er_P_valid_nodal_p) + length(er_N_valid_nodal_p)
length(er_P_valid_nodal_n) + length(er_N_valid_nodal_n)


% group
er_P_train_group_p = intersect(er_P_train, find(group_p_indi));
er_P_train_group_n = intersect(er_P_train, find(group_n_indi));
er_N_train_group_p = intersect(er_N_train, find(group_p_indi));
er_N_train_group_n = intersect(er_N_train, find(group_n_indi));
er_P_valid_group_p = intersect(er_P_valid, find(group_p_indi));
er_P_valid_group_n = intersect(er_P_valid, find(group_n_indi));
er_N_valid_group_p = intersect(er_N_valid, find(group_p_indi));
er_N_valid_group_n = intersect(er_N_valid, find(group_n_indi));
length(er_P_train_group_p) + length(er_N_train_group_p)
length(er_P_train_group_n) + length(er_N_train_group_n)
length(er_P_valid_group_p) + length(er_N_valid_group_p)
length(er_P_valid_group_n) + length(er_N_valid_group_n)

% tumor size
er_P_train_tsize_p = intersect(er_P_train, find(tsize_p_indi));
er_P_train_tsize_n = intersect(er_P_train, find(tsize_n_indi));
er_N_train_tsize_p = intersect(er_N_train, find(tsize_p_indi));
er_N_train_tsize_n = intersect(er_N_train, find(tsize_n_indi));
er_P_valid_tsize_p = intersect(er_P_valid, find(tsize_p_indi));
er_P_valid_tsize_n = intersect(er_P_valid, find(tsize_n_indi));
er_N_valid_tsize_p = intersect(er_N_valid, find(tsize_p_indi));
er_N_valid_tsize_n = intersect(er_N_valid, find(tsize_n_indi));
length(er_P_train_tsize_p) + length(er_N_train_tsize_p)
length(er_P_train_tsize_n) + length(er_N_train_tsize_n)
length(er_P_valid_tsize_p) + length(er_N_valid_tsize_p)
length(er_P_valid_tsize_n) + length(er_N_valid_tsize_n)

% true pathological response
er_P_train_res_p = intersect(er_P_train, find(res_p_indi));
er_P_train_res_n = intersect(er_P_train, find(res_n_indi));
er_N_train_res_p = intersect(er_N_train, find(res_p_indi));
er_N_train_res_n = intersect(er_N_train, find(res_n_indi));
er_P_valid_res_p = intersect(er_P_valid, find(res_p_indi));
er_P_valid_res_n = intersect(er_P_valid, find(res_n_indi));
er_N_valid_res_p = intersect(er_N_valid, find(res_p_indi));
er_N_valid_res_n = intersect(er_N_valid, find(res_n_indi));
length(er_P_train_res_p) + length(er_N_train_res_p)
length(er_P_train_res_n) + length(er_N_train_res_n)
length(er_P_valid_res_p) + length(er_N_valid_res_p)
length(er_P_valid_res_n) + length(er_N_valid_res_n)



% predicted pathological response from dlda 30 predictor
er_P_train_pres_p = intersect(er_P_train, find(pres_p_indi));
er_P_train_pres_n = intersect(er_P_train, find(pres_n_indi));
er_N_train_pres_p = intersect(er_N_train, find(pres_p_indi));
er_N_train_pres_n = intersect(er_N_train, find(pres_n_indi));
er_P_valid_pres_p = intersect(er_P_valid, find(pres_p_indi));
er_P_valid_pres_n = intersect(er_P_valid, find(pres_n_indi));
er_N_valid_pres_p = intersect(er_N_valid, find(pres_p_indi));
er_N_valid_pres_n = intersect(er_N_valid, find(pres_n_indi));
length(er_P_train_pres_p) + length(er_N_train_pres_p)
length(er_P_train_pres_n) + length(er_N_train_pres_n)
length(er_P_valid_pres_p) + length(er_N_valid_pres_p)
length(er_P_valid_pres_n) + length(er_N_valid_pres_n)

%% Part 1 Survival for all 
% All training;
figure(1);
[xpart1 ypart1] = km(bsttime(find(vali_indi==0)),bstcens(find(vali_indi==0)), 1);
xlabel('Years');
ylabel('Survival probability');
train.all     = RPEXEv1('EventTime',bsttime(find(vali_indi==0)),...
    'Censor',bstcens(find(vali_indi==0)),'Trend',3);
% chose decreasing;
[train.all.times train.all.pvalues]
% ans =
%     2.3956    0.9836
%     3.4524    0.7891
%     0.9884    0.5783
%     0.1916    0.3970
%     2.8775    0.1503
%     3.7563    0.0005
[length(bstcens(find(vali_indi==0))) sum(bstcens(find(vali_indi==0)))]
%    310    66
% Interpretation: There is one significant changepoint at 3.76

figure(2);
[xpart1 ypart1] = km_log(bsttime(find(vali_indi==0)),bstcens(find(vali_indi==0)), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([3.7563 3.7563], [-.4 0], 'r');
hold off;


% Survival for all validation;
figure(3);
[xpart1 ypart1] = km(bsttime(find(vali_indi==1)),bstcens(find(vali_indi==1)), 1);
xlabel('time');
ylabel('Survival probability');
valid.all     = RPEXEv1('EventTime',bsttime(find(vali_indi==1)),...
    'Censor',bstcens(find(vali_indi==1)),'Trend',3);
% chose decreasing;
[valid.all.times valid.all.pvalues]
%     2.2998    0.9934
%     1.3388    0.8943
%     4.1725    0.5863
%     2.1793    0.2928
%     3.4935    0.0591
% Interpretation: no significant changepiont
[length(bstcens(find(vali_indi==1))) sum(bstcens(find(vali_indi==1)))]
%    198    45
figure(4);
[xpart2 ypart2] = km_log(bsttime(find(vali_indi==1)),bstcens(find(vali_indi==1)), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([3.4935 3.4935], [-.35 0], 'r');
hold off;

% er_p_indi   = strcmp(ertxt,'P');
% er_n_indi   = strcmp(ertxt,'N');
% chem_s_indi = strcmp(chemotxt,'Rx Sensitive');
% chem_i_indi = strcmp(chemotxt,'Rx Insensitive');
% nodal_p_indi= [strcmp(nodaltxt,'N1'); strcmp(nodaltxt,'N2'); ...
%     strcmp(nodaltxt,'N3')];
% nodal_n_indi= strcmp(nodaltxt,'N0');
% res_n_indi  = strcmp(resptxt,'RD');
% res_p_indi  = strcmp(resptxt,'pCR');



%% Part 2. ER

% Survival for ER training;
figure(1);
[xpart1 ypart1] = km(bsttime(er_P_train),bstcens(er_P_train), 1);
xlabel('Years');
ylabel('Survival probability');
[length(er_P_train) sum(bstcens(er_P_train))]
% ans =
%    174    24
train.er     = RPEXEv1('EventTime',bsttime(er_P_train),...
    'Censor',bstcens(er_P_train),'Trend',3);
% chose decreasing;
[train.er.times train.er.pvalues]
%     0.1150    0.9048
%     0.1916    0.5089
%     0.0520    0.1214
%     3.4524    0.0489
% Interpretation: no significant changepiont. 
%   3.453 is the most plausible changepoint.
figure(2);
[xpart1 ypart1] = km_log(bsttime(er_P_train),bstcens(er_P_train), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([3.4524 3.4524], [-.25 0], 'r');
hold off;


% Survival for ER+ validation;
figure(3);
[xpart1 ypart1] = km(bsttime(er_P_valid),bstcens(er_P_valid), 1);
xlabel('time');
ylabel('Survival probability');
[length(er_P_valid) sum(bstcens(er_P_valid))]
%    123    18
valid.er     = RPEXEv1('EventTime',bsttime(er_P_valid),...
    'Censor',bstcens(er_P_valid),'Trend',3);
% chose decreasing;
[valid.er.times valid.er.pvalues]
%     2.2998    0.9934
%     1.3388    0.8943
%     4.1725    0.5863
%     2.1793    0.2928
%     3.4935    0.0591
% Interpretation: no significant changepiont. 
figure(4);
[xpart2 ypart2] = km_log(bsttime(er_P_valid),bstcens(er_P_valid), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([3.4935 3.4935], [-.25 0], 'r');
hold off;



% Survival for ER Negative training;
figure(1);
[xpart1 ypart1] = km(bsttime(er_N_train),bstcens(er_N_train), 1);
xlabel('Years');
ylabel('Survival probability');
[length(er_N_train) sum(bstcens(er_N_train))]
%    131    41
valid.erN     = RPEXEv1('EventTime',bsttime(er_N_train),...
    'Censor',bstcens(er_N_train),'Trend',3);
% chose decreasing;
[valid.erN.times valid.erN.pvalues]
% ans =
%     1.3771    0.8425
%     0.9884    0.2142
%     2.8720    0.0002
figure(2);
[xpart1 ypart1] = km_log(bsttime(er_N_train),bstcens(er_N_train), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.8720 2.8720], [-.7 0],'r');
hold off;


% Survival for ER negative validation;
figure(3);
[xpart1 ypart1] = km(bsttime(er_N_valid),bstcens(er_N_valid), 1);
xlabel('time');
ylabel('Survival probability');
[length(er_N_valid) sum(bstcens(er_N_valid))]
%     74    27
valid.erN     = RPEXEv1('EventTime',bsttime(er_N_valid),...
    'Censor',bstcens(er_N_valid),'Trend',3);
% chose decreasing;
[valid.erN.times valid.erN.pvalues]
%     1.7029    0.8251
%     1.2594    0.6022
%     2.9514    0.2999
%     2.2998    0.0041
% Interpretation: one significant changepiont at 2.30. 
figure(4);
[xpart2 ypart2] = km_log(bsttime(er_N_valid),bstcens(er_N_valid), 1);
xlabel('Years');
ylabel('Log of survival probability');
xlim([0 8]);
ylim([-.55 0]);
hold on;
plot([2.2998 2.2998], [-.5 0],'r');
hold off;
% sample size:
sum(bstcens(er_N_valid))
% 27
% Note: in ER- validation, the number of event is 27. 
% The p-value is right on the borderline of significance level 0.05. 
% Maybe this is the case where the proposed approach finds the cut off 
% while the alternative approach (Hollandar and Prochan) can not.


%% Part 3. Check and compare the nodal status (pre.N)
% nodal status
% er_P_train_nodal_p = intersect(er_P_train, find(nodal_p_indi));
% er_P_train_nodal_n = intersect(er_P_train, find(nodal_n_indi));
% er_N_train_nodal_p = intersect(er_N_train, find(nodal_p_indi));
% er_N_train_nodal_n = intersect(er_N_train, find(nodal_n_indi));
% er_P_valid_nodal_p = intersect(er_P_valid, find(nodal_p_indi));
% er_P_valid_nodal_n = intersect(er_P_valid, find(nodal_n_indi));
% er_N_valid_nodal_p = intersect(er_N_valid, find(nodal_p_indi));
% er_N_valid_nodal_n = intersect(er_N_valid, find(nodal_n_indi));

% ER+ training, nodal+
[length(er_P_train_nodal_p) sum(bstcens(er_P_train_nodal_p))]
%    22     4
train.er_p_nodal_p     = RPEXEv1('EventTime',bsttime(er_P_train_nodal_p),...
    'Censor',bstcens(er_P_train_nodal_p),'Trend',3);
% [train.er_p_nodal_p.times train.er_p_nodal_p.pvalues]
%     1.5797    0.6093
%     0.1150    0.2126

% ER+ training, nodal-
[length(er_P_train_nodal_n) sum(bstcens(er_P_train_nodal_n))]
%     16     0: no event
%  conclusion: 

% ER- training, nodal+
[length(er_N_train_nodal_p) sum(bstcens(er_N_train_nodal_p))]
%    24     9
train.er_n_nodal_p     = RPEXEv1('EventTime',bsttime(er_N_train_nodal_p),...
    'Censor',bstcens(er_N_train_nodal_p),'Trend',3);
[train.er_n_nodal_p.times train.er_n_nodal_p.pvalues]
%     0.3833    0.8988
%     1.8316    0.7158
%     0.5448    0.6572
%     1.2485    0.0927

% ER- training, nodal-
[length(er_N_train_nodal_n) sum(bstcens(er_N_train_nodal_n))]
%     [8     1]
train.er_n_nodal_n     = RPEXEv1('EventTime',bsttime(er_N_train_nodal_n),...
    'Censor',bstcens(er_N_train_nodal_n),'Trend',3);
[train.er_n_nodal_n.times train.er_n_nodal_n.pvalues]
% only one event



% ER+ validation, nodal+
[length(er_P_valid_nodal_p) sum(bstcens(er_P_valid_nodal_p))]
%     55     10
valid.er_p_nodal_p     = RPEXEv1('EventTime',bsttime(er_P_valid_nodal_p),...
    'Censor',bstcens(er_P_valid_nodal_p),'Trend',3);
[valid.er_p_nodal_p.times valid.er_p_nodal_p.pvalues]
%     2.7652    0.8702
%     1.6318    0.6495
%     0.8186    0.2381

% ER+ validation, nodal-
[length(er_P_valid_nodal_n) sum(bstcens(er_P_valid_nodal_n))]
%     46     4
valid.er_p_nodal_n     = RPEXEv1('EventTime',bsttime(er_P_valid_nodal_n),...
    'Censor',bstcens(er_P_valid_nodal_n),'Trend',3);
[valid.er_p_nodal_n.times valid.er_p_nodal_n.pvalues]
%     0.1369    0.3562
%     1.3388    0.1204


% ER- validation, nodal+
[length(er_N_valid_nodal_p) sum(bstcens(er_N_valid_nodal_p))]
%     36    10
valid.er_n_nodal_p     = RPEXEv1('EventTime',bsttime(er_N_valid_nodal_p),...
    'Censor',bstcens(er_N_valid_nodal_p),'Trend',3);
[valid.er_n_nodal_p.times valid.er_n_nodal_p.pvalues]
%     0.8241    0.8337
%     2.2998    0.0369
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_valid_nodal_p),...
    bstcens(er_N_valid_nodal_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.2998 2.2998], [-.4 0], 'r');
hold off;   


% ER- validation, nodal-
[length(er_N_valid_nodal_n) sum(bstcens(er_N_valid_nodal_n))]
%        23     7
valid.er_n_nodal_n     = RPEXEv1('EventTime',bsttime(er_N_valid_nodal_n),...
    'Censor',bstcens(er_N_valid_nodal_n),'Trend',3);
[valid.er_n_nodal_n.times valid.er_n_nodal_n.pvalues]
%     1.2594    0.8098
%     1.6509    0.0404
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_valid_nodal_n),...
    bstcens(er_N_valid_nodal_n), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([1.6509 1.6509], [-.4 0], 'r');
hold off;   

% Conclusion: for nodal status, the training set does not have enough 
% data value but the validation data set shows that there might be a 
% changepoint for ER- patients. The current significance is 
% between level 0.1 and 0.2

%% Part 4. Grade of breast cancer before the treatment: 1,2 vs. 3 
% group
% er_P_train_group_p = intersect(er_P_train, find(group_p_indi));
% er_P_train_group_n = intersect(er_P_train, find(group_n_indi));
% er_N_train_group_p = intersect(er_N_train, find(group_p_indi));
% er_N_train_group_n = intersect(er_N_train, find(group_n_indi));
% er_P_valid_group_p = intersect(er_P_valid, find(group_p_indi));
% er_P_valid_group_n = intersect(er_P_valid, find(group_n_indi));
% er_N_valid_group_p = intersect(er_N_valid, find(group_p_indi));
% er_N_valid_group_n = intersect(er_N_valid, find(group_n_indi));


% ER+ training, group+
[length(er_P_train_group_p) sum(bstcens(er_P_train_group_p))]
%    46     8
train.er_p_group_p     = RPEXEv1('EventTime',bsttime(er_P_train_group_p),...
    'Censor',bstcens(er_P_train_group_p),'Trend',3);
[train.er_p_group_p.times train.er_p_group_p.pvalues]
%         0.1916    0.1190

% ER+ training, group-
[length(er_P_train_group_n) sum(bstcens(er_P_train_group_n))]
%    19     1
train.er_p_group_n     = RPEXEv1('EventTime',bsttime(er_P_train_group_n),...
    'Censor',bstcens(er_P_train_group_n),'Trend',3);
[train.er_p_group_n.times train.er_p_group_n.pvalues]
% Only 1 event!
    

% ER- training, group+
[length(er_N_train_group_p) sum(bstcens(er_N_train_group_p))]
%   103    29
train.er_n_group_p     = RPEXEv1('EventTime',bsttime(er_N_train_group_p),...
    'Censor',bstcens(er_N_train_group_p),'Trend',3);
[train.er_n_group_p.times train.er_n_group_p.pvalues]
%     1.3771    0.5211
%     0.9637    0.0776
%     2.8720    0.0018
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_train_group_p),...
    bstcens(er_N_train_group_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.8720 2.8720], [-.6 0], 'r');
ylim([-0.6 0])
hold off;   



% ER- training, group-
[length(er_N_train_group_n) sum(bstcens(er_N_train_group_n))]
%     0   0
% No event!!



% ER+ validation, group+
[length(er_P_valid_group_p) sum(bstcens(er_P_valid_group_p))]
%     54    8
valid.er_p_group_p     = RPEXEv1('EventTime',bsttime(er_P_valid_group_p),...
    'Censor',bstcens(er_P_valid_group_p),'Trend',3);
[valid.er_p_group_p.times valid.er_p_group_p.pvalues]
%     1.6318    0.8826
%     1.2731    0.0911

% ER+ validation, group-
[length(er_P_valid_group_n) sum(bstcens(er_P_valid_group_n))]
%     10     0
% All censored, no event


% ER- validation, group+
[length(er_N_valid_group_p) sum(bstcens(er_N_valid_group_p))]
%    54    18
valid.er_n_group_p     = RPEXEv1('EventTime',bsttime(er_N_valid_group_p),...
    'Censor',bstcens(er_N_valid_group_p),'Trend',3);
[valid.er_n_group_p.times valid.er_n_group_p.pvalues]
%     1.7029    0.8270
%     2.9514    0.6877
%     2.2998    0.0144
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_valid_group_p),...
    bstcens(er_N_valid_group_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.2998 2.2998], [-.5 0], 'r');
hold off;   


% ER- validation, group-
[length(er_N_valid_group_n) sum(bstcens(er_N_valid_group_n))]
%        2      0
% All censored.

% Conclusion: Only violation is at 
% ER-, grade+ in both training and testing 
% When grade is nagative the number of events could be either 0 or 1.

%% Part 5. Tumor size
% er_P_train_tsize_p = intersect(er_P_train, find(tsize_p_indi));
% er_P_train_tsize_n = intersect(er_P_train, find(tsize_n_indi));
% er_N_train_tsize_p = intersect(er_N_train, find(tsize_p_indi));
% er_N_train_tsize_n = intersect(er_N_train, find(tsize_n_indi));
% er_P_valid_tsize_p = intersect(er_P_valid, find(tsize_p_indi));
% er_P_valid_tsize_n = intersect(er_P_valid, find(tsize_n_indi));
% er_N_valid_tsize_p = intersect(er_N_valid, find(tsize_p_indi));
% er_N_valid_tsize_n = intersect(er_N_valid, find(tsize_n_indi));



% ER+ training, tsize+
[length(er_P_train_tsize_p) sum(bstcens(er_P_train_tsize_p))]
%    20     2
train.er_p_tsize_p     = RPEXEv1('EventTime',bsttime(er_P_train_tsize_p),...
    'Censor',bstcens(er_P_train_tsize_p),'Trend',3);
[train.er_p_tsize_p.times train.er_p_tsize_p.pvalues]
%    0.1150    0.0821
% this is the actual level; p=0.0821 
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_P_train_tsize_p),...
    bstcens(er_P_train_tsize_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([0.1150 0.1150], [-.12 0], 'r');
hold off;       


% ER+ training, tsize-
[length(er_P_train_tsize_n) sum(bstcens(er_P_train_tsize_n))]
%    1  0
% no event


% ER- training, tsize+
[length(er_N_train_tsize_p) sum(bstcens(er_N_train_tsize_p))]
%    16     7
train.er_n_tsize_p     = RPEXEv1('EventTime',bsttime(er_N_train_tsize_p),...
    'Censor',bstcens(er_N_train_tsize_p),'Trend',3);
[train.er_n_tsize_p.times train.er_n_tsize_p.pvalues]
%       1.2485    0.0550
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_train_tsize_p),...
    bstcens(er_N_train_tsize_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([1.2485 1.2485], [-.8 0], 'r');
hold off;       

% ER- training, tsize+
[length(er_N_train_tsize_n) sum(bstcens(er_N_train_tsize_n))]
%   0   0
% no event



% ER+ validation, tsize+
[length(er_P_valid_tsize_p) sum(bstcens(er_P_valid_tsize_p))]
%    41     6
valid.er_p_tsize_p     = RPEXEv1('EventTime',bsttime(er_P_valid_tsize_p),...
    'Censor',bstcens(er_P_valid_tsize_p),'Trend',3);
[valid.er_p_tsize_p.times valid.er_p_tsize_p.pvalues]
%     1.8973    0.8604
%     2.7652    0.5296
%     0.1369    0.3066

% ER+ validation, tsize-
[length(er_P_valid_tsize_n) sum(bstcens(er_P_valid_tsize_n))]
%     8     2
valid.er_p_tsize_n     = RPEXEv1('EventTime',bsttime(er_P_valid_tsize_n),...
    'Censor',bstcens(er_P_valid_tsize_n),'Trend',3);
[valid.er_p_tsize_n.times valid.er_p_tsize_n.pvalues]
%    3.0910    0.5271






% ER- validation, tsize+
[length(er_N_valid_tsize_p) sum(bstcens(er_N_valid_tsize_p))]
%    30    11
valid.er_n_tsize_p     = RPEXEv1('EventTime',bsttime(er_N_valid_tsize_p),...
    'Censor',bstcens(er_N_valid_tsize_p),'Trend',3);
[valid.er_n_tsize_p.times valid.er_n_tsize_p.pvalues]
%     1.6509    0.7964
%     2.2998    0.0367
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_valid_tsize_p),...
    bstcens(er_N_valid_tsize_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.2998 2.2998], [-.6 0], 'r');
ylim([-.6 0]);
hold off;       


% ER- validation, tsize-
[length(er_N_valid_tsize_n) sum(bstcens(er_N_valid_tsize_n))]
%	2   2
valid.er_n_tsize_n     = RPEXEv1('EventTime',bsttime(er_N_valid_tsize_n),...
    'Censor',bstcens(er_N_valid_tsize_n),'Trend',3);
[valid.er_n_tsize_n.times valid.er_n_tsize_n.pvalues]
%    0.4682    0.4800


% Conclusion: will discuss in a table

%% Part 6. pathological response: positive=pCR, 
% er_P_train_res_p = intersect(er_P_train, find(res_p_indi));
% er_P_train_res_n = intersect(er_P_train, find(res_n_indi));
% er_N_train_res_p = intersect(er_N_train, find(res_p_indi));
% er_N_train_res_n = intersect(er_N_train, find(res_n_indi));
% er_P_valid_res_p = intersect(er_P_valid, find(res_p_indi));
% er_P_valid_res_n = intersect(er_P_valid, find(res_n_indi));
% er_N_valid_res_p = intersect(er_N_valid, find(res_p_indi));
% er_N_valid_res_n = intersect(er_N_valid, find(res_n_indi));



% ER+ training, res+
[length(er_P_train_res_p) sum(bstcens(er_P_train_res_p))]
%    11     0
% no event

% ER+ training, res-
[length(er_P_train_res_n) sum(bstcens(er_P_train_res_n))]
%   161    23
train.er_p_res_n     = RPEXEv1('EventTime',bsttime(er_P_train_res_n),...
    'Censor',bstcens(er_P_train_res_n),'Trend',3);
[train.er_p_res_n.times train.er_p_res_n.pvalues]
%     0.0520    0.9078
%     0.1150    0.8159
%     0.1916    0.3468
%     3.4524    0.0618
% p value is bigger than .1 but smaller than .2. 
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_P_train_res_n),...
    bstcens(er_P_train_res_n), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([3.4524 3.4524], [-.25 0], 'r');
hold off;       

% ER- training, res+
[length(er_N_train_res_p) sum(bstcens(er_N_train_res_p))]
%    45     3
train.er_n_res_p     = RPEXEv1('EventTime',bsttime(er_N_train_res_p),...
    'Censor',bstcens(er_N_train_res_p),'Trend',3);
[train.er_n_res_p.times train.er_n_res_p.pvalues]
%     1.3717    0.4478

% ER- training, res-
[length(er_N_train_res_n) sum(bstcens(er_N_train_res_n))]
%    84    37
train.er_n_res_n     = RPEXEv1('EventTime',bsttime(er_N_train_res_n),...
    'Censor',bstcens(er_N_train_res_n),'Trend',3);
[train.er_n_res_n.times train.er_n_res_n.pvalues]
%    1.3771    0.5116
%    2.7324    0.0003
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_train_res_n),...
    bstcens(er_N_train_res_n), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.7324 2.7324], [-1 0], 'r');
hold off;       


% ER+ validation, res+
[length(er_P_valid_res_p) sum(bstcens(er_P_valid_res_p))]
%    19     0
% no event.

% ER+ validation, res-
[length(er_P_valid_res_n) sum(bstcens(er_P_valid_res_n))]
%    94    17
valid.er_p_res_n     = RPEXEv1('EventTime',bsttime(er_P_valid_res_n),...
    'Censor',bstcens(er_P_valid_res_n),'Trend',3);
[valid.er_p_res_n.times valid.er_p_res_n.pvalues]
%     0.1369    0.8511
%     3.4935    0.1599

% ER- validation, res+
[length(er_N_valid_res_p) sum(bstcens(er_N_valid_res_p))]
%	23   3
valid.er_n_res_p     = RPEXEv1('EventTime',bsttime(er_N_valid_res_p),...
    'Censor',bstcens(er_N_valid_res_p),'Trend',3);
[valid.er_n_res_p.times valid.er_n_res_p.pvalues]
%    1.4428    0.4029



% ER- validation, res-
[length(er_N_valid_res_n) sum(bstcens(er_N_valid_res_n))]
%    45    22
valid.er_n_res_n     = RPEXEv1('EventTime',bsttime(er_N_valid_res_n),...
    'Censor',bstcens(er_N_valid_res_n),'Trend',3);
[valid.er_n_res_n.times valid.er_n_res_n.pvalues]
%     2.2998    0.8517
%     1.2594    0.2580
%     2.9514    0.0426
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_valid_res_n),...
    bstcens(er_N_valid_res_n), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.9514 2.9514], [-.9 0], 'r');
hold off;       

% Conclusion: will discuss in a table


%% Part 7. Predicted pathological response: positive=pCR, 
% er_P_train_pres_p = intersect(er_P_train, find(pres_p_indi));
% er_P_train_pres_n = intersect(er_P_train, find(pres_n_indi));
% er_N_train_pres_p = intersect(er_N_train, find(pres_p_indi));
% er_N_train_pres_n = intersect(er_N_train, find(pres_n_indi));
% er_P_valid_pres_p = intersect(er_P_valid, find(pres_p_indi));
% er_P_valid_pres_n = intersect(er_P_valid, find(pres_n_indi));
% er_N_valid_pres_p = intersect(er_N_valid, find(pres_p_indi));
% er_N_valid_pres_n = intersect(er_N_valid, find(pres_n_indi));

% ER+ training, res+
[length(er_P_train_pres_p) sum(bstcens(er_P_train_pres_p))]
%    10     4
train.er_p_pres_p     = RPEXEv1('EventTime',bsttime(er_P_train_pres_p),...
    'Censor',bstcens(er_P_train_pres_p),'Trend',3);
[train.er_p_pres_p.times train.er_p_pres_p.pvalues]
%     0.7940    0.9017
%     0.1150    0.2577
    
    
% ER+ training, res-
[length(er_P_train_pres_n) sum(bstcens(er_P_train_pres_n))]
%   164    20
train.er_p_pres_n     = RPEXEv1('EventTime',bsttime(er_P_train_pres_n),...
    'Censor',bstcens(er_P_train_pres_n),'Trend',3);
[train.er_p_pres_n.times train.er_p_pres_n.pvalues]
%     0.1916    0.9448
%     3.4524    0.1148
%     0.0520    0.0664
% p value is bigger than .1 but probably smaller than .2. 
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_P_train_pres_n),...
    bstcens(er_P_train_pres_n), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([0.0520 0.0520], [-.2 0], 'r');
hold off;       

% ER- training, res+
[length(er_N_train_pres_p) sum(bstcens(er_N_train_pres_p))]
%   109    35
train.er_n_pres_p     = RPEXEv1('EventTime',bsttime(er_N_train_pres_p),...
    'Censor',bstcens(er_N_train_pres_p),'Trend',3);
[train.er_n_pres_p.times train.er_n_pres_p.pvalues]
%     1.3771    0.1106
%     2.7324    0.0007
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_train_pres_p),...
    bstcens(er_N_train_pres_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.7324 2.7324], [-.7 0], 'r');
hold off;  


% ER- training, res-
[length(er_N_train_pres_n) sum(bstcens(er_N_train_pres_n))]
%    22     6
train.er_n_pres_n     = RPEXEv1('EventTime',bsttime(er_N_train_pres_n),...
    'Censor',bstcens(er_N_train_pres_n),'Trend',3);
[train.er_n_pres_n.times train.er_n_pres_n.pvalues]
%    2.3409    0.2676



% ER+ validation, res+
[length(er_P_valid_pres_p) sum(bstcens(er_P_valid_pres_p))]
%    17     5
valid.er_p_pres_p     = RPEXEv1('EventTime',bsttime(er_P_valid_pres_p),...
    'Censor',bstcens(er_P_valid_pres_p),'Trend',3);
[valid.er_p_pres_p.times valid.er_p_pres_p.pvalues]
%     1.2731    0.9488
%     3.3539    0.6803

% ER+ validation, res-
[length(er_P_valid_pres_n) sum(bstcens(er_P_valid_pres_n))]
%   106    13
valid.er_p_pres_n     = RPEXEv1('EventTime',bsttime(er_P_valid_pres_n),...
    'Censor',bstcens(er_P_valid_pres_n),'Trend',3);
[valid.er_p_pres_n.times valid.er_p_pres_n.pvalues]
%     0.1369    0.8511
%     3.4935    0.1599

% ER- validation, res+
[length(er_N_valid_pres_p) sum(bstcens(er_N_valid_pres_p))]
%    57    25
valid.er_n_pres_p     = RPEXEv1('EventTime',bsttime(er_N_valid_pres_p),...
    'Censor',bstcens(er_N_valid_pres_p),'Trend',3);
[valid.er_n_pres_p.times valid.er_n_pres_p.pvalues]
%     1.2594    0.9539
%     1.7029    0.6004
%     2.9514    0.3890
%     2.1793    0.0057
figure(1);
[xpart1 ypart1] = km_log(bsttime(er_N_valid_pres_p),...
    bstcens(er_N_valid_pres_p), 1);
xlabel('Years');
ylabel('Log of survival probability');
hold on;
plot([2.1793 2.1793], [-.7 0], 'r');
hold off;       



% ER- validation, res-
[length(er_N_valid_pres_n) sum(bstcens(er_N_valid_pres_n))]
%    17     2
valid.er_n_pres_n     = RPEXEv1('EventTime',bsttime(er_N_valid_pres_n),...
    'Censor',bstcens(er_N_valid_pres_n),'Trend',3);
[valid.er_n_pres_n.times valid.er_n_pres_n.pvalues]
%    0.7830    0.4825



%% Other codes

% Figure 2a in the paper;
figure(3);
subplot(1,2,1);
[xpart31 ypart31] = km(bsttime(intersect(find(vali_indi==0),find(chem_indi==1))),...
    bstcens(intersect(find(vali_indi==0),find(chem_indi==1))));
ylim([0 1]);
xlim([0 5]);
hold on;
[xpart32 ypart32] = km2(bsttime(intersect(find(vali_indi==0),find(chem_indi==0))),...
    bstcens(intersect(find(vali_indi==0),find(chem_indi==0))));
legend('sensitive', 'insensitive');
hold off;
subplot(1,2,2);
[xpart33 ypart33] = km_log(bsttime(intersect(find(vali_indi==0),find(chem_indi==1))),...
    bstcens(intersect(find(vali_indi==0),find(chem_indi==1))), 0);
xlim([0 5]);
hold on;
[xpart34 ypart34] = km_log2(bsttime(intersect(find(vali_indi==0),find(chem_indi==0))),...
    bstcens(intersect(find(vali_indi==0),find(chem_indi==0))), 0);
legend('sensitive', 'insensitive');
hold off;



% Figure 3b in the paper;
figure(4);
subplot(1,2,1);
[xpart31 ypart31] = km(bsttime(intersect(intersect(find(vali_indi==1),find(chem_indi==1)), find(er_indi==0))),...
    bstcens(intersect(intersect(find(vali_indi==1),find(chem_indi==1)), find(er_indi==0))));
ylim([0 1]);
xlim([0 5]);
hold on;
[xpart32 ypart32] = km2(bsttime(intersect(intersect(find(vali_indi==1),find(chem_indi==0)), find(er_indi==0))),...
    bstcens(intersect(intersect(find(vali_indi==1),find(chem_indi==0)), find(er_indi==0))));
legend('sensitive ER-', 'insensitive ER-');
hold off;
subplot(1,2,2);
[xpart33 ypart33] = km_log(bsttime(intersect(intersect(find(vali_indi==1),find(chem_indi==1)), find(er_indi==0))),...
    bstcens(intersect(intersect(find(vali_indi==1),find(chem_indi==1)), find(er_indi==0))), 0);
xlim([0 5]);
hold on;
[xpart34 ypart34] = km_log2(bsttime(intersect(intersect(find(vali_indi==1),find(chem_indi==0)), find(er_indi==0))),...
    bstcens(intersect(intersect(find(vali_indi==1),find(chem_indi==0)), find(er_indi==0))), 0);
legend('sensitive ER-', 'insensitive ER-');
hold off;

