%%%% Fitting generalized gamma distribution to the data 


Ns    = 30;


textseq = {'1,  Kaplan-Meier',...
'2,  Kaplan-Meier with window of 30',...
'3,  Kaplan-Meier with window of 10',...
'4,  Gamma',...
'5,  Weibull',...
'6,  Lognormal',...
'7,  Exponential with censor at 0',...
'8,  Exponential with censor at 60',...
'9,  Exponential with censor at 90',...
'10,  Exponential with censor at 120',...
'11,  Exponential with censor at 150',...
'12,  Exponential with censor at 180',...
'13,  Exponential with censor at 210',...
'14,  PEXE (Kim and Proschan)',...
'15,  RPEXE'};


beta = 4.7;
comb1 = [0.92 0.8];
comb2 = [1.08 0.8];
comb3 = [0.75 0.5];
comb4 = [1.25 0.5];
comb5 = [0.85 1.05];
comb6 = [0.85 0.95];
comb7 = [0.5 0.75];
comb8 = [0.5 1.25];
comb9 = [0.92 1.2];
comb10 = [1.08 1.2];
comb11 = [0.75 1.5];
comb12 = [1.25 1.5];
comb13 = [1.15 1.05];
comb14 = [1.15 0.95];
comb15 = [1.5 0.75];
comb16 = [1.5 1.25];

% plot([0.92 1.08 0.75 1.25], [0.8 0.8 0.5 0.5],'o'); % increase-decrease
% plot([0.85 0.85 0.5 0.5], [1.05 0.95 0.75 1.25],'^'); % increase
% plot([0.92 1.08 0.75 1.25], [1.2 1.2 1.5 1.5],'p'); % decrease-increase
% plot([1.15 1.15 1.5 1.5], [1.05 0.95 0.75 1.25],'s'); % decrease




addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
res_gam1 = gen_gam_test(beta,comb1(1),comb1(2),Ns);
res_gam2 = gen_gam_test(beta,comb2(1),comb2(2),Ns);
res_gam3 = gen_gam_test(beta,comb3(1),comb3(2),Ns);
res_gam4 = gen_gam_test(beta,comb4(1),comb4(2),Ns);
res_gam5 = gen_gam_test(beta,comb5(1),comb5(2),Ns);
res_gam6 = gen_gam_test(beta,comb6(1),comb6(2),Ns);
res_gam7 = gen_gam_test(beta,comb7(1),comb7(2),Ns);
res_gam8 = gen_gam_test(beta,comb8(1),comb8(2),Ns);
res_gam9 = gen_gam_test(beta,comb9(1),comb9(2),Ns);
res_gam10 = gen_gam_test(beta,comb10(1),comb10(2),Ns);
res_gam11 = gen_gam_test(beta,comb11(1),comb11(2),Ns);
res_gam12 = gen_gam_test(beta,comb12(1),comb12(2),Ns);
res_gam13 = gen_gam_test(beta,comb13(1),comb13(2),Ns);
res_gam14 = gen_gam_test(beta,comb14(1),comb14(2),Ns);
res_gam15 = gen_gam_test(beta,comb15(1),comb15(2),Ns);
res_gam16 = gen_gam_test(beta,comb16(1),comb16(2),Ns);

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.92 0.8)'},'A1:A1');
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam1,'A2:C16');
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam1(:,3)),'D2:D16');

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.08 0.8)'},'A17:A17')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam2,'A18:C32')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam2(:,3)),'D18:D32')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(.75 0.5)'},'A33:A33')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam3,'A34:C48')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam3(:,3)),'D34:D48')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.25 0.5)'},'A49:A49')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam4,'A50:C64')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam4(:,3)),'D50:D64')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.85 1.05)'},'A65:A65')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam5,'A66:C80')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam5(:,3)),'D66:D80')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.85 0.95)'},'A81:A81')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam6,'A82:C96')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam6(:,3)),'D82:D96')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.5 0.75)'},'A97:A97')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam7,'A98:C112')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam7(:,3)),'D98:D112')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.5 1.25)'},'A113:A113')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam8,'A114:C128')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam8(:,3)),'D114:D128')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.92 1.2)'},'A129:A129')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam9,'A130:C144')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam9(:,3)),'D130:D144')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.08 1.2)'},'A145:A145')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam10,'A146:C160')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam10(:,3)),'D146:D160')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(0.75 1.5)'},'A161:A161')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam11,'A162:C176')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam11(:,3)),'D162:D176')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.25 1.5)'},'A177:A177')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam12,'A178:C192')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam12(:,3)),'D178:D192')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.15 1.05)'},'A193:A193')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam13,'A194:C208')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam13(:,3)),'D194:D208')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.15 0.95)'},'A209:A209')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam14,'A210:C224')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam14(:,3)),'D210:D224')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.5 0.75)'},'A225:A225')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam15,'A226:C240')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam15(:,3)),'D226:D240')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',{'gen gamma (lam sig)=(1.5 1.25)'},'A241:A241')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',res_gam16,'A242:C256')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_30.xls',rankavec(res_gam16(:,3)),'D242:D256')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
Ns    = 100;


textseq = {'1,  Kaplan-Meier',...
'2,  Kaplan-Meier with window of 30',...
'3,  Kaplan-Meier with window of 10',...
'4,  Gamma',...
'5,  Weibull',...
'6,  Lognormal',...
'7,  Exponential with censor at 0',...
'8,  Exponential with censor at 60',...
'9,  Exponential with censor at 90',...
'10,  Exponential with censor at 120',...
'11,  Exponential with censor at 150',...
'12,  Exponential with censor at 180',...
'13,  Exponential with censor at 210',...
'14,  PEXE (Kim and Proschan)',...
'15,  RPEXE'};


beta = 4.7;
comb1 = [0.92 0.8];
comb2 = [1.08 0.8];
comb3 = [0.75 0.5];
comb4 = [1.25 0.5];
comb5 = [0.85 1.05];
comb6 = [0.85 0.95];
comb7 = [0.5 0.75];
comb8 = [0.5 1.25];
comb9 = [0.92 1.2];
comb10 = [1.08 1.2];
comb11 = [0.75 1.5];
comb12 = [1.25 1.5];
comb13 = [1.15 1.05];
comb14 = [1.15 0.95];
comb15 = [1.5 0.75];
comb16 = [1.5 1.25];

% plot([0.92 1.08 0.75 1.25], [0.8 0.8 0.5 0.5],'o'); % increase-decrease
% plot([0.85 0.85 0.5 0.5], [1.05 0.95 0.75 1.25],'^'); % increase
% plot([0.92 1.08 0.75 1.25], [1.2 1.2 1.5 1.5],'p'); % decrease-increase
% plot([1.15 1.15 1.5 1.5], [1.05 0.95 0.75 1.25],'s'); % decrease




addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
res_gam1 = gen_gam_test(beta,comb1(1),comb1(2),Ns);
res_gam2 = gen_gam_test(beta,comb2(1),comb2(2),Ns);
res_gam3 = gen_gam_test(beta,comb3(1),comb3(2),Ns);
res_gam4 = gen_gam_test(beta,comb4(1),comb4(2),Ns);
res_gam5 = gen_gam_test(beta,comb5(1),comb5(2),Ns);
res_gam6 = gen_gam_test(beta,comb6(1),comb6(2),Ns);
res_gam7 = gen_gam_test(beta,comb7(1),comb7(2),Ns);
res_gam8 = gen_gam_test(beta,comb8(1),comb8(2),Ns);
res_gam9 = gen_gam_test(beta,comb9(1),comb9(2),Ns);
res_gam10 = gen_gam_test(beta,comb10(1),comb10(2),Ns);
res_gam11 = gen_gam_test(beta,comb11(1),comb11(2),Ns);
res_gam12 = gen_gam_test(beta,comb12(1),comb12(2),Ns);
res_gam13 = gen_gam_test(beta,comb13(1),comb13(2),Ns);
res_gam14 = gen_gam_test(beta,comb14(1),comb14(2),Ns);
res_gam15 = gen_gam_test(beta,comb15(1),comb15(2),Ns);
res_gam16 = gen_gam_test(beta,comb16(1),comb16(2),Ns);

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.92 0.8)'},'A1:A1');
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam1,'A2:C16');
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam1(:,3)),'D2:D16');

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.08 0.8)'},'A17:A17')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam2,'A18:C32')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam2(:,3)),'D18:D32')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(.75 0.5)'},'A33:A33')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam3,'A34:C48')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam3(:,3)),'D34:D48')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.25 0.5)'},'A49:A49')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam4,'A50:C64')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam4(:,3)),'D50:D64')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.85 1.05)'},'A65:A65')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam5,'A66:C80')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam5(:,3)),'D66:D80')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.85 0.95)'},'A81:A81')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam6,'A82:C96')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam6(:,3)),'D82:D96')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.5 0.75)'},'A97:A97')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam7,'A98:C112')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam7(:,3)),'D98:D112')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.5 1.25)'},'A113:A113')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam8,'A114:C128')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam8(:,3)),'D114:D128')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.92 1.2)'},'A129:A129')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam9,'A130:C144')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam9(:,3)),'D130:D144')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.08 1.2)'},'A145:A145')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam10,'A146:C160')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam10(:,3)),'D146:D160')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(0.75 1.5)'},'A161:A161')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam11,'A162:C176')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam11(:,3)),'D162:D176')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.25 1.5)'},'A177:A177')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam12,'A178:C192')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam12(:,3)),'D178:D192')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.15 1.05)'},'A193:A193')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam13,'A194:C208')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam13(:,3)),'D194:D208')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.15 0.95)'},'A209:A209')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam14,'A210:C224')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam14(:,3)),'D210:D224')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.5 0.75)'},'A225:A225')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam15,'A226:C240')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam15(:,3)),'D226:D240')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',{'gen gamma (lam sig)=(1.5 1.25)'},'A241:A241')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',res_gam16,'A242:C256')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest1_100.xls',rankavec(res_gam16(:,3)),'D242:D256')
