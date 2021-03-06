function plotprederror2(struct)


% plot the predictive errors to compare multiple predictors
lencol = 9;
lenrow = 8;
ts     = [5 30 60 90 120 150 180 210];
prop   = zeros(1,8);
beta   = 4.7;


% for colomn i in the structure
for i = 3:3:lencol,
    setting(i).error1 = zeros(15,8);
    setting(i).error2 = zeros(15,8);
    setting(i).errorI = zeros(15,8);
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
lambda = setting(i).shape;
sigma  = setting(i).scale; 
for its = 1:length(ts),
    if setting(i).shape >= 0, 
        prop(its) = 1-gamcdf(lambda^(-2)*(exp(-beta)*ts(its))^(lambda/sigma),lambda^(-2),1);
    else
        prop(its) = gamcdf(lambda^(-2)*(exp(-beta)*ts(its))^(lambda/sigma),lambda^(-2),1);
    end
end
subplot(1,2,1);
plot(prop,setting(i).error2(1,:),'b',prop,setting(i).error2(2,:),'m',prop,setting(i).error2(3,:),'cx-',...
    prop,setting(i).error2(4,:),'gd-',prop,setting(i).error2(5,:),'go-',prop,setting(i).error2(6,:),'g^-',...
    prop,setting(i).error2(7,:),'g--',prop,setting(i).error2(8,:),'g-.',prop,setting(i).error2(9,:),'g',...
    prop,setting(i).error2(10,:),'rd-',prop,setting(i).error2(11,:),'ro-',prop,setting(i).error2(12,:),'r^-',...
    prop,setting(i).error2(13,:),'r--',prop,setting(i).error2(14,:),'r-.',prop,setting(i).error2(15,:),'r');
xlim([0 1]);
ylim([0 0.01+max(max(setting(i).errorI))]);
xlabel('Survival Probability');
ylabel(['RMSPE, Shape(lambda)=' num2str(struct(i,1).shape) '; Scale(sigma)=' num2str(struct(i,1).scale)]);

subplot(1,2,2);
plot(prop,setting(i).errorI(1,:),'b',prop,setting(i).errorI(2,:),'m',prop,setting(i).errorI(3,:),'cx-',...
    prop,setting(i).errorI(4,:),'gd-',prop,setting(i).errorI(5,:),'go-',prop,setting(i).errorI(6,:),'g^-',...
    prop,setting(i).errorI(7,:),'g--',prop,setting(i).errorI(8,:),'g-.',prop,setting(i).errorI(9,:),'g',...
    prop,setting(i).errorI(10,:),'rd-',prop,setting(i).errorI(11,:),'ro-',prop,setting(i).errorI(12,:),'r^-',...
    prop,setting(i).errorI(13,:),'r--',prop,setting(i).errorI(14,:),'r-.',prop,setting(i).errorI(15,:),'r');
xlim([0 1]);
ylim([0 0.01+max(max(setting(i).errorI))]);
xlabel('Survival Probability');
ylabel('Maximum Predictive Error');
legend('kme raw','exponential','pexe',...
    'umrpexe alpha= 1','umrpexe alpha= .7','umrpexe alpha=.4',...
    'umrpexe alpha= .2','umrpexe alpha= .1','umrpexe alpha=.05',...
    'monorpexe alpha=1','monorpexe alpha=.7','monorpexe alpha=.4',...
    'monorpexe alpha=.2','monorpexe alpha=.1','monorpexe alpha=.05');
end
