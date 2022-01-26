smth = 101;
xs = 0:50000/(2*pi):(1e8)+10000;

eccs0 = readmatrix("sim_0_ecc_dmc.csv");
eccs1 = readmatrix("sim_1_ecc_dmc.csv");
eccs2 = readmatrix("sim_2_ecc_dmc.csv");
eccs3 = readmatrix("sim_3_ecc_dmc.csv");
eccs4 = readmatrix("sim_4_ecc_dmc.csv");
eccs5 = readmatrix("sim_5_ecc_dmc.csv");
% eccs6 = readmatrix("sim_6_ecc_dmc.csv");
% eccs7 = readmatrix("sim_7_ecc_dmc.csv");
% eccs8 = readmatrix("sim_8_ecc_dmc.csv");
% eccs9 = readmatrix("sim_9_ecc_dmc.csv");

% smoothed0 = movmean(eccs0, smth);
% smoothed1 = movmean(eccs1, smth);
% smoothed2 = movmean(eccs2, smth);
% smoothed3 = movmean(eccs3, smth);
% smoothed4 = movmean(eccs4, smth);
% smoothed5 = movmean(eccs5, smth);
% smoothed6 = movmean(eccs6, smth);
% smoothed7 = movmean(eccs7, smth);
% smoothed8 = movmean(eccs8, smth);
% smoothed9 = movmean(eccs9, smth);

hold on
plot(xs(1:length(eccs0)), eccs0);
plot(xs(1:length(eccs1)), eccs1);
plot(xs(1:length(eccs2)), eccs2);
plot(xs(1:length(eccs3)), eccs3);
plot(xs(1:length(eccs4)), eccs4);
plot(xs(1:length(eccs5)), eccs5);
% plot(xs(1:length(eccs6)), eccs6);
% plot(xs(1:length(eccs7)), eccs7);
% plot(xs(1:length(eccs8)), eccs8);
% plot(xs(1:length(eccs9)), eccs9);
%legend("Sim 0", "Sim 1", "Sim 2", "Sim 3", "location", "northwest");