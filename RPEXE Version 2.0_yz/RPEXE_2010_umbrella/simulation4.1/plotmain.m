load res_gam_melanomas;
plotprederror(res_gam_melanomas)
load res_gam_pancreas;
plotprederror(res_gam_pancreas)
load res_gam_colon;
plotprederror(res_gam_colon)
load res_gam_lung;
plotprederror(res_gam_lung)




%check the alpha value
Ns = [20 30 50 80 100 120 150];

for i = 1:7,
    crivalmono5(i)  = log(exp(-4.2331-0.3938*log(Ns(i))));
    crivalmono10(i) = log(exp(-3.4481-0.3847*log(Ns(i))));
    crivalmono20(i) = crivalmono10(i) + crivalmono10(i) - crivalmono5(i);
    crivalum5(i)  = log(exp(-1.2234-2.7044*log(Ns(i))));
    crivalum10(i) = log(exp(-1.1951-2.0238*log(Ns(i))));
    crivalum20(i)   = crivalum10(i) + crivalum10(i) - crivalum5(i);
end


plot(Ns, crivalmono5, 'b--', Ns, crivalmono10, 'b.-', Ns, crivalmono20, 'bo-',...
    Ns, crivalum5, 'r--', Ns, crivalum10, 'r.-', Ns ,crivalum20, 'ro-');
title('Plot of log(\alpha^\star) vs sample size');
legend('montone .05','montone .1','montone .2',...
    'umbrella .05','umbrella .1','umbrella .2');
