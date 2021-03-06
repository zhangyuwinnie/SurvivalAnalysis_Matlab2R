function plotprederror1_2(struct)
% plot the predictive errors 

%struct = res_gam_melanomas;

lencol = 9;
lenrow = 8;
ts     = [5 30 60 90 120 150 180 210];
prop   = zeros(1,8);
beta   = 4.7;

% for colomn i in the structure
for i = 3:3:lencol,
    setting(i).error1 = zeros(16,7);
    setting(i).error2 = zeros(16,7);
    setting(i).errorI = zeros(16,7);
    setting(i).scale  = struct(i,1).scale;
    setting(i).shape  = struct(i,1).shape;

    for j = 1:lenrow,
        setting(i).error1(:,j) = struct(i,j).r1error; % R1 error
        setting(i).error2(:,j) = struct(i,j).r2error; % R2 error
        setting(i).errorI(:,j) = struct(i,j).maxerror; % RInfinite error
    end    
end

for i = 3:3:lencol,

figure(i);
size(setting(i).error2)
lambda = setting(i).shape;
sigma  = setting(i).scale; 
for its = 1:length(ts),
    if setting(i).shape >= 0, 
        prop(its) = 1-gamcdf(lambda^(-2)*(exp(-beta)*ts(its))^(lambda/sigma),lambda^(-2),1);
    else
        prop(its) = gamcdf(lambda^(-2)*(exp(-beta)*ts(its))^(lambda/sigma),lambda^(-2),1);
    end
end
% subplot(1,3,1);
% plot(ts,setting(i).error1(1,:),'b',ts,setting(i).error1(2,:),'bo-',ts,setting(i).error1(3,:),'b-.',...
%     ts,setting(i).error1(4,:),'gd-',ts,setting(i).error1(5,:),'go-',ts,setting(i).error1(6,:),'g-.',...
%     ts,setting(i).error1(7,:),'k--',ts,setting(i).error1(8,:),'k-.',...
%     ts,setting(i).error1(9,:),'k',ts,setting(i).error1(10,:),'cx-',...
%     ts,setting(i).error1(11,:),'r-.',ts,setting(i).error1(12,:),'r--',ts,setting(i).error1(13,:),'r',...
%     ts,setting(i).error1(14,:),'ro-',ts,setting(i).error1(15,:),'rs-',ts,setting(i).error1(16,:),'r^-');
% xlim([0 220]);
% xlabel('month');
% ylabel(['R 1 penalty, shape(lambda)=' num2str(struct(i,1).shape) '; scale(sigma)=' num2str(struct(i,1).scale)]);
% 
% subplot(1,2,1);
plot(prop,setting(i).error2(1,:),'b',prop,setting(i).error2(3,:),'bo-',...
    prop,setting(i).error2(4,:),'k',prop,setting(i).error2(6,:),'ko-',...
    prop,setting(i).error2(7,:),'k--',prop,setting(i).error2(8,:),'k-.',...
    prop,setting(i).error2(9,:),'k:',prop,setting(i).error2(10,:),'cx-',...
    prop,setting(i).error2(12,:),'g-.',prop,setting(i).error2(11,:),'g^-',...
    prop,setting(i).error2(15,:),'r-.',prop,setting(i).error2(14,:),'r^-');
xlim([0 1]);
% ylim([0 0.01+max(max(setting(i).errorI))]);
xlabel('Survival Probability');
ylabel(['RMSPE, shape \lambda = ' num2str(struct(i,1).shape) '; scale \sigma = ' num2str(struct(i,1).scale)]);
% subplot(1,2,2);
% plot(prop,setting(i).errorI(1,:),'b',prop,setting(i).errorI(3,:),'bo-',...
%     prop,setting(i).errorI(4,:),'m',prop,setting(i).errorI(6,:),'mo-',...
%     prop,setting(i).errorI(7,:),'k--',prop,setting(i).errorI(8,:),'k-.',...
%     prop,setting(i).errorI(9,:),'k:',prop,setting(i).errorI(10,:),'cx-',...
%     prop,setting(i).errorI(12,:),'g-.',prop,setting(i).errorI(11,:),'g^-',...
%     prop,setting(i).errorI(15,:),'r-.',prop,setting(i).errorI(14,:),'r^-');
% xlim([0 1]);
% ylim([0 0.01+max(max(setting(i).errorI))]);
% xlabel('Survival Probability');
% ylabel('Maximum predictive error');
legend('KME raw','KME window 30',...
    'exponential','exp. & censored at 150',...
    'gamma','Weibull','lognormal','PEXE',...
    'UMRPEXE \alpha = .1','UMRPEXE \alpha = .4',...
    'MONORPEXE \alpha =.1','MONORPEXE \alpha =.4');
end




















