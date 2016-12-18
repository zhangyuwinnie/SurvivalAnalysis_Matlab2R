function maxe = plotprederror1_inpaper2(struct)
% plot the predictive errors 


lenrow = 9;
prop   = [1:9]/10; 
lencol = size(struct(1,1).r1error,1);

% for colomn i in the structure
for i = 1:4,
    setting(i).error1 = zeros(lencol,9);
    setting(i).error2 = zeros(lencol,9);
    setting(i).errorI = zeros(lencol,9);
    setting(i).scale  = struct(i,1).scale;
    setting(i).shape  = struct(i,1).shape;
    for j = 1:lenrow,
        setting(i).error1(:,j) = struct(i,j).r1error; % R1 error
        setting(i).error2(:,j) = struct(i,j).r2error'; % R2 error
        setting(i).errorI(:,j) = struct(i,j).maxerror; % RInfinite error
    end
    figure(i);
    lambda = setting(i).shape;
    sigma  = setting(i).scale; 
    subplot(1,3,1); 
    plot(prop,setting(i).error2(1,:),'b',... %prop,setting(i).error2(3,:),'bo-',...
    prop,setting(i).error2(4,:),'k',prop,setting(i).error2(6,:),'ko-',...
    prop,setting(i).error2(7,:),'m--',prop,setting(i).error2(8,:),'m-.',...
    prop,setting(i).error2(9,:),'m:',prop,setting(i).error2(10,:),'bx-');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    xlabel('Survival Probability');
    ylabel(['RMSPE, shape \lambda = ' num2str(struct(i,1).shape) ...
    '; scale \sigma = ' num2str(struct(i,1).scale)]);
    legend('kme', 'exp.', 'exp q0.5',...
                    'gam','Wei','logn','pexe');
    subplot(1,3,2);
    plot(prop,setting(i).error2(11,:),'r:',prop,setting(i).error2(12,:),'r-.',...
    prop,setting(i).error2(13,:),'r-o',prop,setting(i).error2(14,:),'r-*',...
    prop,setting(i).error2(15,:),'r',prop,setting(i).error2(16,:),'g:',...
    prop,setting(i).error2(17,:),'g-.',prop,setting(i).error2(18,:),'g-o',...
    prop,setting(i).error2(19,:),'g-*',prop,setting(i).error2(20,:),'g',... %prop,setting(i).error2(21,:),'g->',prop,setting(i).error2(22,:),'r-<',...
    prop,setting(i).error2(23,:),'b:',prop,setting(i).error2(24,:),'b-.',...
    prop,setting(i).error2(25,:),'b-o',prop,setting(i).error2(26,:),'b-*',...
    prop,setting(i).error2(27,:),'b');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    xlabel('Survival Probability');
    legend('umrpexe10','umrpexe25','umrpexe50','umrpexe75','umrpexe100',...
    'monorpexe10','monorpexe25','monorpexe50','monorpexe75','monorpexe100',...     %'monorpexe40','umrpexe60',
    'berpexe10','berpexe25','berpexe50',...
    'berpexe75','berpexe100');
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






alpha1 = [0.1 0.25 0.5 0.75 1]';
maxratio = [6.95 1.31 1.70 1.52 1.31
            6.95 1.21 1.46 1.46 1.31
            6.95 1.17 1.24 1.26 1.31
            6.95 1.20 1.20 1.20 1.31
            6.95 1.41 1.45 1.29 1.31];
rmsrratio= [3.91 1.17 1.45 1.34 1.17
            3.91 1.12 1.28 1.29 1.17
            3.91 1.10 1.13 1.14 1.17
            3.91 1.10 1.14 1.12 1.17
            3.91 1.14 1.24 1.15 1.17];
figure1 = {'monorpexe','umrpexe','berpexe','kme'};
figure(5);
plot(alpha1,maxratio(:,2),'k-o',alpha1,maxratio(:,3),'k-^',...
    alpha1,maxratio(:,4),'k-*',alpha1,maxratio(:,5),'k--');
ylim([1.0 1.8]);
xlabel('\alpha');
ylabel('Max Ratio');
legend(figure1)
figure(6);
plot(alpha1,rmsrratio(:,2),'k-o',alpha1,rmsrratio(:,3),'k-^',...
    alpha1,rmsrratio(:,4),'k-*',alpha1,rmsrratio(:,5),'k--');
ylim([1.0 1.8]);
xlabel('\alpha');
ylabel('Root Mean Squared Ratio');
legend(figure1)

