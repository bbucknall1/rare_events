N = 24;
smth = 1001;
%xs = 0:50000/(2*pi):(2e8)+10000;
xs = (0:25134)*(50000/(2*pi));
L = length(xs);

eccs = zeros(N+1, L);
smoothed = zeros(N+1, L);

eccs0 = readmatrix("sim_0_ecc_dmc.csv");
eccs3 = readmatrix("sim_3_ecc_dmc.csv");
eccs6 = readmatrix("sim_6_ecc_dmc.csv");
eccs9 = readmatrix("sim_9_ecc_dmc.csv");
eccs12 = readmatrix("sim_12_ecc_dmc.csv");
eccs15 = readmatrix("sim_15_ecc_dmc.csv");
eccs18 = readmatrix("sim_18_ecc_dmc.csv");
eccs21 = readmatrix("sim_21_ecc_dmc.csv");

smoothed0 = movmean(eccs0, smth);
smoothed3 = movmean(eccs3, smth);
smoothed6 = movmean(eccs6, smth);
smoothed9 = movmean(eccs9, smth);
smoothed12 = movmean(eccs12, smth);
smoothed15 = movmean(eccs15, smth);
smoothed18 = movmean(eccs18, smth);
smoothed21 = movmean(eccs21, smth);

hold on
plot(xs(1:length(smoothed0)), smoothed0);
plot(xs(1:length(smoothed3)), smoothed3);
plot(xs(1:length(smoothed6)), smoothed6);
plot(xs(1:length(smoothed9)), smoothed9);
plot(xs(1:length(smoothed12)), smoothed12);
plot(xs(1:length(smoothed15)), smoothed15);
plot(xs(1:length(smoothed18)), smoothed18);
plot(xs(1:length(smoothed21)), smoothed21);

ylim([0 0.5]);


% hold on;
% for i = 0:3:N
%     eccs(i+1,:) = readmatrix("sim_"+num2str(i)+"_ecc_dmc.csv");
%     smoothed(i+1,:) = movmean(eccs(i+1,:), smth);
%     plot(xs, smoothed(i+1,:));
% end
% ylim([0 0.5]);