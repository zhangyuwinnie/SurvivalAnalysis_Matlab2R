% %%%% Fitting generalized gamma distribution to the data 

clear;
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');

% melanomas and pancreatic cancer
%%%% Fitting generalized gamma distribution to the data 
Ns    = 50;
beta  = 1;
prop = [1:9]/10;
N     = 1000;

% shape and scale parameters for regional melanomas period 3
%   localized colon cancer period 3, distant colon cancer period 3, 
%   and regional melanomas period 1.
shape = [3.34 0.93 0.41 -0.34];
scale = [0.27 1.03 1.41  1.93];

reading = {'kme'; 'kme 10'; 'kme 30'; 'exponential';...
    'exp censored at 210';'exp censored at 150';'gam';'wei';'logn';'pexe';...
    'umrpexe10';'umrpexe25';'umrpexe50';'umrpexe75';'umrpexe100';...
    'monorpexe10';'monorpexe25';'monorpexe50';'monorpexe75';'monorpexe100'};



for i = 1:4,
    lambda = shape(i);
    sigma  = scale(i);
    if lambda > 0,
        for k = 1:length(prop),
            ts(k) = (lambda^(2)*gaminv(1-prop(k),lambda^(-2),1))...
                ^(sigma/lambda)*exp(beta);
        end
    else
        for k = 1:length(prop),
            ts(k) = (lambda^(2)*gaminv(prop(k),lambda^(-2),1))...
                ^(sigma/lambda)*exp(beta);
        end
    end  
    for j = 1:length(ts),
        res_gam(i,j).res       = gen_gam_test_inpaper(beta,shape(i),...
            scale(i),Ns,ts(j),N);
        res_gam(i,j).variance  = res_gam(i,j).res(:,1);
        res_gam(i,j).bias      = res_gam(i,j).res(:,2);
        res_gam(i,j).r1error   = res_gam(i,j).res(:,5);
        res_gam(i,j).r2error   = res_gam(i,j).res(:,3);
        res_gam(i,j).maxerror  = res_gam(i,j).res(:,4);
        res_gam(i,j).r1penrank = rankavec(res_gam(i,j).res(:,5));
        res_gam(i,j).r2penrank = rankavec(res_gam(i,j).res(:,3));
        res_gam(i,j).rIpenrank = rankavec(res_gam(i,j).res(:,4));
        res_gam(i,j).label     = reading;
        res_gam(i,j).shape     = shape(i);
        res_gam(i,j).scale     = scale(i);
        res_gam(i,j).beta      = beta;
        res_gam(i,j).t         = ts(j);
        res_gam(i,j).samplesize= Ns;
        res_gam(i,j).simuN     = N;
    end
end
res_gam_fourpanel = res_gam;
save res_gam_fourpanel;


load res_gam_fourpanel;
maxe = plotprederror1_inpaper(res_gam_fourpanel)

