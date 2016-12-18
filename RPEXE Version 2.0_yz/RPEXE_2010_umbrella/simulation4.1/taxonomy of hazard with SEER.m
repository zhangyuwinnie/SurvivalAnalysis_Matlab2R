%%% Plot the generalized gamma parameters 

x1 = zeros(201,1);
x2 = zeros(201,1);
for i = 1:201,
    x1(i) = 1+(i-1)/100;
    x2(i) = 1/x1(i);
end

figure(1);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
xlabel('\sigma');
ylabel('\lambda');
plot([0.92 1.08 0.75 1.25], [0.8 0.8 0.5 0.5],'o'); % increase-decrease
plot([0.85 0.85 0.5 0.5], [1.05 0.95 0.75 1.25],'^'); % increase
plot([0.92 1.08 0.75 1.25], [1.2 1.2 1.5 1.5],'p'); % decrease-increase
plot([1.15 1.15 1.5 1.5], [1.05 0.95 0.75 1.25],'s'); % decrease
hold off;



%%% Input and Plot the SEER Data
%%% All vectors are of the form: 
% [LocalDecade1 
% LocalDecade2
% LocalDecade3
% RegionalDecade1
% RegionalDecade2
% RegionalDecade3
% DistantDecade1
% DistantDecade2
% DistantDecade3]

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
breastscale = seerest(:,3);
breastshape = seerest(:,4);
stomachscale= seerest(:,5);
stomachshape= seerest(:,6);
pancreasscale = seerest(:,7);
pancreasshape = seerest(:,8);
colonscale  = seerest(:,9);
colonshape  = seerest(:,10);
melanomasscale = seerest(:,11);
melanomasshape = seerest(:,12);

% plot the four figures
x1 = zeros(251,1);
x2 = zeros(251,1);
x3 = zeros(251,1);
for i = 1:251,
    x1(i) = (i-1)/100;
    x2(i) = 1/x1(i);
    x3(i) = 1/(2*x1(i));
end

figure(1);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
xlabel('\sigma');
ylabel('\lambda');
xlim([0 2.5]);
ylim([-0.5 2.5]);
title('Scale (\sigma) vs. shape (\lambda) plot for generalized Gamma distribution','Fontsize',14);
hold off;

figure(2);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(melanomasscale,melanomasshape,'+');
plot(breastscale,breastshape,'*');
xlabel('\sigma');
ylabel('\lambda');
xlim([0 2.5]);
ylim([-0.5 2.5]);
title('Scale (\sigma) vs. shape (\lambda) plot with melanomas (+) and breast cancer (*) estimates','Fontsize',11);
hold off;

figure(3);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(colonscale,colonshape,'^');
xlabel('\sigma');
ylabel('\lambda');
xlim([1 2.5]);
ylim([-0.5 1]);
title('Scale (\sigma) vs. shape (\lambda) plot with colon cancer (triangles) estimate','Fontsize',11);
hold off;

figure(4);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(stomachscale,stomachshape,'d');
plot(pancreasscale,pancreasshape,'s');
plot(lungscale,lungshape,'x');
xlabel('\sigma');
ylabel('\lambda');
xlim([1 2.5]);
ylim([-0.5 1]);
title('Lung cancer (x), stomach cancer (diamond), and pancreas (square) estimates','Fontsize',10);
hold off;




figure(5);
subplot(2,2,1);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
xlabel('\sigma');
ylabel('\lambda');
xlim([0 2.5]);
ylim([-0.5 2.5]);
title('Scale (\sigma) vs. shape (\lambda) plot for generalized Gamma distribution','Fontsize',12);
hold off;
subplot(2,2,2);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(melanomasscale,melanomasshape,'+');
plot(breastscale,breastshape,'*');
xlabel('\sigma');
ylabel('\lambda');
xlim([0 2.5]);
ylim([-0.5 2.5]);
title('Scale (\sigma) vs. shape (\lambda) plot with melanomas (+) and breast cancer (*) estimates','Fontsize',12);
hold off;
subplot(2,2,3);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(colonscale,colonshape,'^');
xlabel('\sigma');
ylabel('\lambda');
xlim([1 2.5]);
ylim([-0.5 1]);
title('Scale (\sigma) vs. shape (\lambda) plot with colon cancer (triangles) estimate','Fontsize',12);
hold off;
subplot(2,2,4);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,x3,'--');
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(stomachscale,stomachshape,'d');
plot(pancreasscale,pancreasshape,'s');
plot(lungscale,lungshape,'x');
xlabel('\sigma');
ylabel('\lambda');
xlim([1 2.5]);
ylim([-0.5 1]);
title('Lung cancer (x), stomach cancer (diamond), and pancreas (square) estimates','Fontsize',12);
hold off;





% ENAR meeting graph
figure(6);
subplot(3,1,1);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,zeros(251,1),'-.');
plot(x1,ones(251,1),':');
xlabel('Scale \sigma','Fontsize',14);
ylabel('Shape \lambda','Fontsize',14);
xlim([0 2.5]);
ylim([-0.5 2.5]);
title('Scale (\sigma) vs. shape (\lambda) for generalized Gamma distribution','Fontsize',14);
hold off;
subplot(3,1,2);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,zeros(251,1),'-.');
plot(x1,ones(251,1),':');
plot(melanomasscale(1:3),melanomasshape(1:3),'g+'); % green: local
plot(melanomasscale(4:6),melanomasshape(4:6),'m+'); % pink: regional
plot(melanomasscale(7:9),melanomasshape(7:9),'k+'); % black: distant
plot(breastscale(1:3),breastshape(1:3),'g*');
plot(breastscale(4:6),breastshape(4:6),'m*');
plot(breastscale(7:9),breastshape(7:9),'k*');
xlabel('Scale \sigma','Fontsize',14);
ylabel('Shape \lambda','Fontsize',14);
xlim([0 2.5]);
ylim([-0.5 3]);
title('Melanomas(+) and breast(*)','Fontsize',14);
hold off;
subplot(3,1,3);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,zeros(251,1),'-.');
plot(x1,ones(251,1),':');
plot(stomachscale(1:3),stomachshape(1:3),'dg');
plot(stomachscale(4:6),stomachshape(4:6),'dm');
plot(stomachscale(7:9),stomachshape(7:9),'kd');
plot(pancreasscale(1:3),pancreasshape(1:3),'gs');
plot(pancreasscale(4:6),pancreasshape(4:6),'ms');
plot(pancreasscale(7:9),pancreasshape(7:9),'ks');
plot(lungscale(1:3),lungshape(1:3),'gx');
plot(lungscale(4:6),lungshape(4:6),'mx');
plot(lungscale(7:9),lungshape(7:9),'kx');
plot(colonscale(1:3),colonshape(1:3),'g^');
plot(colonscale(4:6),colonshape(4:6),'m^');
plot(colonscale(7:9),colonshape(7:9),'k^');
xlabel('Scale \sigma','Fontsize',14);
ylabel('Shape \lambda','Fontsize',14);
xlim([1 2.5]);
ylim([-0.5 1]);
title('Lung(x), stomach(diamond), pancreas(square), and colon(triangle)','Fontsize',14);
hold off;



% separate gragh on in the manuscript;
figure(7);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
xlabel('Scale \sigma','Fontsize',14);
ylabel('Shape \lambda','Fontsize',14);
xlim([0 2.5]);
ylim([-0.5 2.5]);
%title('Scale (\sigma) vs. shape (\lambda) for generalized Gamma distribution','Fontsize',14);
hold off;

figure(8);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(melanomasscale(1:3),melanomasshape(1:3),'g+'); % green: local
plot(melanomasscale(4:6),melanomasshape(4:6),'m+'); % pink: regional
plot(melanomasscale(7:9),melanomasshape(7:9),'k+'); % black: distant
plot(breastscale(1:3),breastshape(1:3),'g*');
plot(breastscale(4:6),breastshape(4:6),'m*');
plot(breastscale(7:9),breastshape(7:9),'k*');
xlabel('Scale \sigma','Fontsize',14);
ylabel('Shape \lambda','Fontsize',14);
xlim([0 2.5]);
ylim([-0.5 3]);
%title('Melanomas(+) and breast(*)','Fontsize',14);
hold off;

figure(9);
plot(x1,x2);
hold on;
plot(x2,x1);
plot(x1,x1);
plot(x2,x2);
plot(x1,zeros(251,1),':');
plot(x1,ones(251,1),'-.');
plot(stomachscale(1:3),stomachshape(1:3),'dg');
plot(stomachscale(4:6),stomachshape(4:6),'dm');
plot(stomachscale(7:9),stomachshape(7:9),'kd');
plot(pancreasscale(1:3),pancreasshape(1:3),'gs');
plot(pancreasscale(4:6),pancreasshape(4:6),'ms');
plot(pancreasscale(7:9),pancreasshape(7:9),'ks');
plot(lungscale(1:3),lungshape(1:3),'gx');
plot(lungscale(4:6),lungshape(4:6),'mx');
plot(lungscale(7:9),lungshape(7:9),'kx');
plot(colonscale(1:3),colonshape(1:3),'g^');
plot(colonscale(4:6),colonshape(4:6),'m^');
plot(colonscale(7:9),colonshape(7:9),'k^');
xlabel('Scale \sigma','Fontsize',14);
ylabel('Shape \lambda','Fontsize',14);
xlim([1 2.5]);
ylim([-0.5 1]);
%title('Lung(x), stomach(diamond), pancreas(square), and colon(triangle)','Fontsize',14);
hold off;



