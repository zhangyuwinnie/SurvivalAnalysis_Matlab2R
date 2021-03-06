%%%% Fitting generalized gamma distribution to the data 

Ns    = 50;

beta = 4.7;
seerest = [1.73	0.27	1.24	0.54	1.79	0.36	1.93	0.03	1.31	0.57	2.04	-0.1
1.66	0.29	0.43	2.29	1.67	0.37	1.81	-0.22	1.12	0.79	0.55	2.14
1.56	0.42	0.43	2.07	1.7	0.45	1.91	-0.13	1.03	0.93	0.37	2.82
1.63	0.07	1.34	0.37	1.63	0.11	1.49	0.11	1.61	0.23	1.93	-0.34
1.66	0.08	1.43	0.29	1.65	0.09	1.51	0.12	1.63	0.24	0.37	2.56
1.58	0.24	0.46	2.24	1.65	0.13	1.48	0.25	1.75	0.24	0.27	3.34
1.63	0.33	1.31	0.52	1.49	0.23	1.53	0.14	1.45	0.26	1.9	-0.15
1.42	0.26	1.36	0.47	1.47	0.2	1.59	0.11	1.48	0.29	2	-0.06
1.4	0.34	1.29	0.63	1.47	0.26	1.48	0.25	1.41	0.41	1.77	-0.04];

lungscale   = seerest(:,1);
lungshape   = seerest(:,2);

[lungscale lungshape];


plot(lungscale,lungshape,'o');

% plot([0.92 1.08 0.75 1.25], [0.8 0.8 0.5 0.5],'o'); % increase-decrease
% plot([0.85 0.85 0.5 0.5], [1.05 0.95 0.75 1.25],'^'); % increase
% plot([0.92 1.08 0.75 1.25], [1.2 1.2 1.5 1.5],'p'); % decrease-increase
% plot([1.15 1.15 1.5 1.5], [1.05 0.95 0.75 1.25],'s'); % decrease




addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
for i = 1:length(lungscale);
    res_gam(i).res = gen_gam_test(beta,lungshape(i),lungscale(i),Ns);
end

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.73 0.27)'},'A1:A1');
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(1).res,'A2:C19');
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(1).res(:,3)),'D2:D19');

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.66 0.29)'},'A20:A20')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(2).res,'A21:C38')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(2).res(:,3)),'D21:D38')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.56 0.42)'},'A39:A39')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(3).res,'A40:C57')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(3).res(:,3)),'D40:D57')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.63 0.07)'},'A58:A58')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(4).res,'A59:C76')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(3).res(:,3)),'D59:D76')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.66 0.08)'},'A77:A77')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(5).res,'A78:C95')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(5).res(:,3)),'D78:D95')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.58 0.24)'},'A96:A96')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(6).res,'A97:C114')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(6).res(:,3)),'D97:D114')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.63 0.33)'},'A115:A115')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(7).res,'A116:C133')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(7).res(:,3)),'D116:D133')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.42 0.26)'},'A134:A134')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(8).res,'A135:C152')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(8).res(:,3)),'D135:D152')

xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',{'gen gamma (lam sig)=(1.40 0.34)'},'A153:A153')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',res_gam(9).res,'A154:C171')
xlswrite('G:\Projects\HanGang\Piece_wise_Expo\Exponential\program\gengamest2_50.xls',rankavec(res_gam(9).res(:,3)),'D154:D171')

