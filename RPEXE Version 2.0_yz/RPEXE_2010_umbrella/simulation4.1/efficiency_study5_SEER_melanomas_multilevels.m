% melanomas and pancreatic cancer
%%%% Fitting generalized gamma distribution to the data 
Ns     = 50;
beta   = 4.7;
N      = 1000;
lencol = 9;
lenrow = 8;
ts     = [5 30 60 90 120 150 180 210];

seerest = [1.73	0.27	1.24	0.54	1.79	0.36	1.93	0.03	1.31	0.57	2.04	-0.1
1.66	0.29	0.43	2.29	1.67	0.37	1.81	-0.22	1.12	0.79	0.55	2.14
1.56	0.42	0.43	2.07	1.7	0.45	1.91	-0.13	1.03	0.93	0.37	2.82
1.63	0.07	1.34	0.37	1.63	0.11	1.49	0.11	1.61	0.23	1.93	-0.34
1.66	0.08	1.43	0.29	1.65	0.09	1.51	0.12	1.63	0.24	0.37	2.56
1.58	0.24	0.46	2.24	1.65	0.13	1.48	0.25	1.75	0.24	0.27	3.34
1.63	0.33	1.31	0.52	1.49	0.23	1.53	0.14	1.45	0.26	1.9	-0.15
1.42	0.26	1.36	0.47	1.47	0.2	1.59	0.11	1.48	0.29	2	-0.06
1.4	0.34	1.29	0.63	1.47	0.26	1.48	0.25	1.41	0.41	1.77	-0.04];

addpath('G:\Projects\HanGang\Piece_wise_Expo\program\matlab_functions\RPEXE_2010_umbrella');
addpath('G:\Projects\HanGang\Piece_wise_Expo\program\commonly used functions');
melanomasscale = seerest(:,11);
melanomasshape = seerest(:,12);

reading = {'kme raw'; 'exp_raw'; 'pexe'; 'umrpexe100'; 'umrpexe70'; ...
    'umrpexe40'; 'umrpexe20'; 'umrpexe10'; 'umrpexe5';...
    'monorpexe100';'monorpexe70';'monorpexe40';'monorpexe20';...
    'monorpexe10';'monorpexe5'};


for i = 3:3:lencol,
    for j = 1:length(ts),
        res_gam_melanomas(i,j).res       = gen_gam_test_multilevels...
            (beta,melanomasshape(i),melanomasscale(i),Ns,ts(j),N);
        res_gam_melanomas(i,j).variance  = res_gam_melanomas(i,j).res(:,1);
        res_gam_melanomas(i,j).bias      = res_gam_melanomas(i,j).res(:,2);
        res_gam_melanomas(i,j).r1error   = res_gam_melanomas(i,j).res(:,5);
        res_gam_melanomas(i,j).r2error   = res_gam_melanomas(i,j).res(:,3);
        res_gam_melanomas(i,j).maxerror  = res_gam_melanomas(i,j).res(:,4);
        res_gam_melanomas(i,j).r1penrank = rankavec(res_gam_melanomas(i,j).res(:,5));
        res_gam_melanomas(i,j).r2penrank = rankavec(res_gam_melanomas(i,j).res(:,3));
        res_gam_melanomas(i,j).rIpenrank = rankavec(res_gam_melanomas(i,j).res(:,4));
        res_gam_melanomas(i,j).label   = reading;
        res_gam_melanomas(i,j).shape     = melanomasshape(i);
        res_gam_melanomas(i,j).scale     = melanomasscale(i);
        res_gam_melanomas(i,j).beta      = beta;
        res_gam_melanomas(i,j).t         = ts(j);
        res_gam_melanomas(i,j).samplesize= Ns;
        res_gam_melanomas(i,j).simuN     = N;
    end
end
res_gam_melanomas_m =  res_gam_melanomas;
struct = res_gam_melanomas_m;

plotprederror2(struct);
save resgam_melanomas_m;


% just do the plot
load resgam_melanomas_m;
plotprederror2_2(struct);
% figure (a) in the paper is figure 6 in the plot
















