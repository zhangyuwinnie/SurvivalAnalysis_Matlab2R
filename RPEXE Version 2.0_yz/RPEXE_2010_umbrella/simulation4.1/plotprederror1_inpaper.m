function maxe = plotprederror1_inpaper(struct)
% plot the predictive errors 


lenrow = 9;
prop   = [1:9]/10; 

% for colomn i in the structure
for i = 1:4,
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
    prop,setting(i).error2(19,:),'g-*',prop,setting(i).error2(20,:),'g');
    xlim([0 1]);
    ylim([0 0.01+max(max(setting(i).error2([1 4 6 7 8 9 10],:)))]);
    xlabel('Survival Probability');
    legend('umrpexe10','umrpexe25','umrpexe50','umrpexe75','umrpexe100',...
    'monorpexe10','monorpexe25','monorpexe50','monorpexe75','monorpexe100',...
    'monorpexe40','umrpexe60','berpexe10','berpexe25','berpexe50',...
    'berpexe75','berpexe100');
end


%%% calculate the table that sum up the predictive error 
% Do the table;
% Column 1
maxerrors = zeros(17,4);
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
        setting(i).error2(20,:)]');
        maxerrors(:,i) = a'/min(a);
end
%max(maxerrors,[],2)
maxe = [maxerrors max(maxerrors,[],2)];














