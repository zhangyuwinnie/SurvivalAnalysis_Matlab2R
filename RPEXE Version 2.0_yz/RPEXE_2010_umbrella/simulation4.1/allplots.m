


% three panels in figure1
%
load res_gam_melanomas_c;
plotprederror1_2(res_gam_melanomas_c);
% figure 6 is regional melanomas
load res_gam_colon_c;
plotprederror1_2(res_gam_colon_c);
% figures 3 and 9 are local colon and distant colon










% three panels in figure2
%
load resgam_melanomas_m;
plotprederror2_2(struct);
% figure 6 is regional melanomas
load resgam_colon_m;
plotprederror2_2(struct);
% figures 3 and 9 are local colon and distant colon

