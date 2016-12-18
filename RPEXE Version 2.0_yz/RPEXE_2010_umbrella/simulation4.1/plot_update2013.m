%cd D:\work\RPEXE\Coding\program\matlab_functions\RPEXE_2010_umbrella\simulation4.1;

load res_gam_fourpanel; 

lenrow = 9;
prop   = [1:9]/10; 

struct=res_gam_fourpanel;

lencol = size(struct(1,1).r1error,1);

% for colomn i in the structure
for i = 1:4,
    setting(i).error1 = zeros(20,9);
    setting(i).error2 = zeros(20,9);
    setting(i).errorI = zeros(20,9);
    setting(i).bias2   = zeros(20,9);
    setting(i).variance = zeros(20,9);
    setting(i).scale  = struct(i,1).scale;
    setting(i).shape  = struct(i,1).shape;
    for j = 1:lenrow,
        setting(i).error1(:,j) = struct(i,j).r1error; % R1 error
        setting(i).error2(:,j) = struct(i,j).r2error; % R2 error
        setting(i).errorI(:,j) = struct(i,j).maxerror; % RInfinite error
        setting(i).bias2(:,j)  = struct(i,j).bias*struct(i,j).bias; % Bias
        setting(i).variance(:,j) = struct(i,j).variance;
    end
    hFig = figure(2*i-1);
    set(hFig, 'Position', [100 100 200 260]);
    lambda = setting(i).shape;
    sigma  = setting(i).scale; 
%    subplot(1,3,1); 
    plot(prop,setting(i).error2(1,:),'b',... %prop,setting(i).error2(3,:),'bo-',...
    prop,setting(i).error2(4,:),'k',prop,setting(i).error2(6,:),'ko-',...
    prop,setting(i).error2(7,:),'m--',prop,setting(i).error2(8,:),'m-.',...
    prop,setting(i).error2(9,:),'m:',prop,setting(i).error2(10,:),'bx-');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    xlabel('Survival Probability');
    ylabel(['RMSPE, shape \lambda = ' num2str(struct(i,1).shape) ...
    '; scale \sigma = ' num2str(struct(i,1).scale)]);
%     legend('kme', 'exp.', 'exp q0.5',...
%                     'gam','Wei','logn','pexe');
%%    subplot(1,3,2);
    hFig = figure(2*i);
    set(hFig, 'Position', [100 100 200 260]);
    plot(prop,setting(i).error2(11,:),'r:',prop,setting(i).error2(12,:),'r-.',...
    prop,setting(i).error2(13,:),'r-o',prop,setting(i).error2(14,:),'r-*',...
    prop,setting(i).error2(15,:),'r',prop,setting(i).error2(16,:),'g:',...
    prop,setting(i).error2(17,:),'g-.',prop,setting(i).error2(18,:),'g-o',...
    prop,setting(i).error2(19,:),'g-*',prop,setting(i).error2(20,:),'g');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    xlabel('Survival Probability');
    ylabel(['RMSPE, shape \lambda = ' num2str(struct(i,1).shape) ...
    '; scale \sigma = ' num2str(struct(i,1).scale)]);
%     legend('umrpexe10','umrpexe25','umrpexe50','umrpexe75','umrpexe100',...
%     'monorpexe10','monorpexe25','monorpexe50','monorpexe75','monorpexe100',...
%     'monorpexe40','umrpexe60','berpexe10','berpexe25','berpexe50',...
%     'berpexe75','berpexe100');
end




%%% plot the legend: removed all labels

i=1;
    setting(i).error1 = zeros(20,9);
    setting(i).error2 = zeros(20,9);
    setting(i).errorI = zeros(20,9);
    setting(i).scale  = struct(i,1).scale;
    setting(i).shape  = struct(i,1).shape;
    for j = 1:lenrow,
        setting(i).error1(:,j) = struct(i,j).r1error; % R1 error
        setting(i).error2(:,j) = struct(i,j).r2error; % R2 error
        setting(i).errorI(:,j) = struct(i,j).maxerror; % RInfinite error
    end
    hFig = figure(2*i-1);
    set(hFig, 'Position', [100 100 200 260]);
    lambda = setting(i).shape;
    sigma  = setting(i).scale; 
%    subplot(1,3,1); 
    plot(prop,setting(i).error2(1,:),'b',... %prop,setting(i).error2(3,:),'bo-',...
    prop,setting(i).error2(4,:),'k',prop,setting(i).error2(6,:),'ko-',...
    prop,setting(i).error2(7,:),'m--',prop,setting(i).error2(8,:),'m-.',...
    prop,setting(i).error2(9,:),'m:',prop,setting(i).error2(10,:),'bx-');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    legend('kme', 'exp.', 'exp q0.5',...
                    'gam','Wei','logn','pexe');
%%    subplot(1,3,2);
    hFig = figure(2*i);
    set(hFig, 'Position', [100 100 200 260]);
    plot(prop,setting(i).error2(11,:),'r:',prop,setting(i).error2(12,:),'r-.',...
    prop,setting(i).error2(13,:),'r-o',prop,setting(i).error2(14,:),'r-*',...
    prop,setting(i).error2(15,:),'r',prop,setting(i).error2(16,:),'g:',...
    prop,setting(i).error2(17,:),'g-.',prop,setting(i).error2(18,:),'g-o',...
    prop,setting(i).error2(19,:),'g-*',prop,setting(i).error2(20,:),'g');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    legend('umrpexe10','umrpexe25','umrpexe50','umrpexe75','umrpexe100',...
    'monorpexe10','monorpexe25','monorpexe50','monorpexe75','monorpexe100',...
    'monorpexe40','umrpexe60','berpexe10','berpexe25','berpexe50',...
    'berpexe75','berpexe100');





%% Compute bias, variance, and minimum predictive error;
struct=res_gam_fourpanel2;

lencol = size(struct(1,1).r1error,1);

% for colomn i in the structure
for i = 1:4,
    setting(i).error1 = zeros(lencol,9);
    setting(i).error2 = zeros(lencol,9);
    setting(i).errorI = zeros(lencol,9);
    setting(i).bias2   = zeros(lencol,9);
    setting(i).variance = zeros(lencol,9);
    setting(i).scale  = struct(i,1).scale;
    setting(i).shape  = struct(i,1).shape;
    for j = 1:lenrow,
        setting(i).error1(:,j) = struct(i,j).r1error; % R1 error
        setting(i).error2(:,j) = struct(i,j).r2error; % R2 error
        setting(i).errorI(:,j) = struct(i,j).maxerror; % RInfinite error
        setting(i).bias(:,j)  = abs(struct(i,j).bias); % absolute bias
        setting(i).std(:,j) = sqrt(struct(i,j).variance);
    end
end


%%% calculate the table that sum up the predictive error 
% Do the table;
% Column 1
maxerrors = zeros(lencol-3-2,4);
for i = 1:4,
    a = max([setting(i).error2(1,:) 
        %setting(i).error2(3,:) 
        setting(i).error2(4,:)
        setting(i).error2(6,:) 
        setting(i).error2(7,:)
        setting(i).error2(8,:) 
        setting(i).error2(9,:) 
        setting(i).error2(10,:)
        setting(i).error2(11,:) 
        setting(i).error2(12,:) 
        setting(i).error2(13,:) 
        setting(i).error2(14,:) 
        setting(i).error2(15,:)
        setting(i).error2(16,:) 
        setting(i).error2(17,:) 
        setting(i).error2(18,:) 
        setting(i).error2(19,:) 
        setting(i).error2(20,:)
%         setting(i).error2(21,:) 
%         setting(i).error2(22,:) 
        setting(i).error2(23,:) 
        setting(i).error2(24,:) 
        setting(i).error2(25,:)
        setting(i).error2(26,:) 
        setting(i).error2(27,:)
        ]');
        maxerrors(:,i) = a'/min(a);
end
%max(maxerrors,[],2)
maxe = [maxerrors max(maxerrors,[],2) rankavec(max(maxerrors,[],2))];
meanmax = zeros(size(maxerrors,1),2);
for i = 1:size(maxerrors,1),
    meanmax(i,1) = sqrt((maxe(i,1)^2+maxe(i,2)^2+maxe(i,3)^2+maxe(i,4)^2)/4);
end
meanmax(:,2) = rankavec(meanmax(:,1));
maxe = [maxe meanmax];



%%% calculate the table that sum up the predictive biase 
% Do the table;
% Column 1
maxbiases = zeros(lencol-3-2,4);
for i = 1:4,
    a = max([setting(i).bias(1,:) 
        %setting(i).bias(3,:) 
        setting(i).bias(4,:)
        setting(i).bias(6,:) 
        setting(i).bias(7,:)
        setting(i).bias(8,:) 
        setting(i).bias(9,:) 
        setting(i).bias(10,:)
        setting(i).bias(11,:) 
        setting(i).bias(12,:) 
        setting(i).bias(13,:) 
        setting(i).bias(14,:) 
        setting(i).bias(15,:)
        setting(i).bias(16,:) 
        setting(i).bias(17,:) 
        setting(i).bias(18,:) 
        setting(i).bias(19,:) 
        setting(i).bias(20,:)
%         setting(i).bias(21,:) 
%         setting(i).bias(22,:) 
        setting(i).bias(23,:) 
        setting(i).bias(24,:) 
        setting(i).bias(25,:)
        setting(i).bias(26,:) 
        setting(i).bias(27,:)
        ]');
        maxbiases(:,i) = a'/min(a);
end
maxe = [maxbiases max(maxbiases,[],2) rankavec(max(maxbiases,[],2))];
meanmax = zeros(size(maxbiases,1),2);
for i = 1:size(maxbiases,1),
    meanmax(i,1) = sqrt((maxe(i,1)^2+maxe(i,2)^2+maxe(i,3)^2+maxe(i,4)^2)/4);
end
meanmax(:,2) = rankavec(meanmax(:,1));
maxe = [maxe meanmax];





%%% calculate the table that sum up the predictive stde 
% Do the table;
% Column 1
maxstdes = zeros(lencol-3-2,4);
for i = 1:4,
    a = max([setting(i).std(1,:) 
        %setting(i).std(3,:) 
        setting(i).std(4,:)
        setting(i).std(6,:) 
        setting(i).std(7,:)
        setting(i).std(8,:) 
        setting(i).std(9,:) 
        setting(i).std(10,:)
        setting(i).std(11,:) 
        setting(i).std(12,:) 
        setting(i).std(13,:) 
        setting(i).std(14,:) 
        setting(i).std(15,:)
        setting(i).std(16,:) 
        setting(i).std(17,:) 
        setting(i).std(18,:) 
        setting(i).std(19,:) 
        setting(i).std(20,:)
%         setting(i).std(21,:) 
%         setting(i).std(22,:) 
        setting(i).std(23,:) 
        setting(i).std(24,:) 
        setting(i).std(25,:)
        setting(i).std(26,:) 
        setting(i).std(27,:)
        ]');
        maxstdes(:,i) = a'/min(a);
end
maxe = [maxstdes max(maxstdes,[],2) rankavec(max(maxstdes,[],2))];
meanmax = zeros(size(maxstdes,1),2);
for i = 1:size(maxstdes,1),
    meanmax(i,1) = sqrt((maxe(i,1)^2+maxe(i,2)^2+maxe(i,3)^2+maxe(i,4)^2)/4);
end
meanmax(:,2) = rankavec(meanmax(:,1));
maxe = [maxe meanmax];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % last figure 
% alpha1 = [0.1 0.25 0.5 0.75 1]';
% maxratio = [6.95 1.31 1.70 1.52 1.31
%             6.95 1.21 1.46 1.46 1.31
%             6.95 1.17 1.24 1.26 1.31
%             6.95 1.20 1.20 1.20 1.31
%             6.95 1.41 1.45 1.29 1.31];
% rmsrratio= [3.91 1.17 1.45 1.34 1.17
%             3.91 1.12 1.28 1.29 1.17
%             3.91 1.10 1.13 1.14 1.17
%             3.91 1.10 1.14 1.12 1.17
%             3.91 1.14 1.24 1.15 1.17];
% figure1 = {'monorpexe','umrpexe','berpexe','kme'};
% figure(5);
% plot(alpha1,maxratio(:,2),'k-o',alpha1,maxratio(:,3),'k-^',...
%     alpha1,maxratio(:,4),'k-*',alpha1,maxratio(:,5),'k--');
% ylim([1.0 1.8]);
% xlabel('\alpha');
% ylabel('Max Ratio');
% legend(figure1)
% figure(6);
% plot(alpha1,rmsrratio(:,2),'k-o',alpha1,rmsrratio(:,3),'k-^',...
%     alpha1,rmsrratio(:,4),'k-*',alpha1,rmsrratio(:,5),'k--');
% ylim([1.0 1.8]);
% xlabel('\alpha');
% ylabel('Root Mean Squared Ratio');
% legend(figure1)







