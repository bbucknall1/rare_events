N = 24;
smth = 1001;
xs = 0:50000/(2*pi):(2e8)+10000;
%xs = (0:105625)*(50000/(2*pi));
L = length(xs);

eccs = zeros(N+1, L);
smoothed = zeros(N+1, L);

%eccs0 = readmatrix("sim_0_ecc_dmc.csv");
eccs3 = readmatrix("sim_3_ecc_dmc.csv");
eccs6 = readmatrix("sim_6_ecc_dmc.csv");
eccs9 = readmatrix("sim_9_ecc_dmc.csv");
eccs12 = readmatrix("sim_12_ecc_dmc.csv");
eccs15 = readmatrix("sim_15_ecc_dmc.csv");
eccs18 = readmatrix("sim_18_ecc_dmc.csv");
eccs21 = readmatrix("sim_21_ecc_dmc.csv");

eccs1 = readmatrix("sim_1_ecc_dmc.csv");
eccs4 = readmatrix("sim_4_ecc_dmc.csv");
eccs7 = readmatrix("sim_7_ecc_dmc.csv");
eccs10 = readmatrix("sim_10_ecc_dmc.csv");
eccs13 = readmatrix("sim_13_ecc_dmc.csv");
eccs16 = readmatrix("sim_16_ecc_dmc.csv");
eccs19 = readmatrix("sim_19_ecc_dmc.csv");
eccs22 = readmatrix("sim_22_ecc_dmc.csv");

eccs2 = readmatrix("sim_2_ecc_dmc.csv");
eccs5 = readmatrix("sim_5_ecc_dmc.csv");
eccs8 = readmatrix("sim_8_ecc_dmc.csv");
eccs11 = readmatrix("sim_11_ecc_dmc.csv");
eccs14 = readmatrix("sim_14_ecc_dmc.csv");
eccs17 = readmatrix("sim_17_ecc_dmc.csv");
eccs20 = readmatrix("sim_20_ecc_dmc.csv");
eccs23 = readmatrix("sim_23_ecc_dmc.csv");

%smoothed0 = movmean(eccs0, smth);
smoothed3 = movmean(eccs3, smth);
smoothed6 = movmean(eccs6, smth);
smoothed9 = movmean(eccs9, smth);
smoothed12 = movmean(eccs12, smth);
smoothed15 = movmean(eccs15, smth);
smoothed18 = movmean(eccs18, smth);
smoothed21 = movmean(eccs21, smth);

smoothed1 = movmean(eccs1, smth);
smoothed4 = movmean(eccs4, smth);
smoothed7 = movmean(eccs7, smth);
smoothed10 = movmean(eccs10, smth);
smoothed13 = movmean(eccs13, smth);
smoothed16 = movmean(eccs16, smth);
smoothed19 = movmean(eccs19, smth);
smoothed22 = movmean(eccs22, smth);

smoothed2 = movmean(eccs2, smth);
smoothed5 = movmean(eccs5, smth);
smoothed8 = movmean(eccs8, smth);
smoothed11 = movmean(eccs11, smth);
smoothed14 = movmean(eccs14, smth);
smoothed17 = movmean(eccs17, smth);
smoothed20 = movmean(eccs20, smth);
smoothed23 = movmean(eccs23, smth);

hold on
%plot(xs(1:length(smoothed0)), smoothed0);
plot(xs(1:length(smoothed3)), smoothed3);
plot(xs(1:length(smoothed6)), smoothed6);
plot(xs(1:length(smoothed9)), smoothed9);
plot(xs(1:length(smoothed12)), smoothed12);
plot(xs(1:length(smoothed15)), smoothed15);
plot(xs(1:length(smoothed18)), smoothed18);
plot(xs(1:length(smoothed21)), smoothed21);

plot(xs(1:length(smoothed1)), smoothed1);
plot(xs(1:length(smoothed4)), smoothed4);
plot(xs(1:length(smoothed7)), smoothed7);
plot(xs(1:length(smoothed10)), smoothed10);
plot(xs(1:length(smoothed13)), smoothed13);
plot(xs(1:length(smoothed16)), smoothed16);
plot(xs(1:length(smoothed19)), smoothed19);
plot(xs(1:length(smoothed22)), smoothed22);

plot(xs(1:length(smoothed2)), smoothed2);
plot(xs(1:length(smoothed5)), smoothed5);
plot(xs(1:length(smoothed8)), smoothed8);
plot(xs(1:length(smoothed11)), smoothed11);
plot(xs(1:length(smoothed14)), smoothed14);
plot(xs(1:length(smoothed17)), smoothed17);
plot(xs(1:length(smoothed20)), smoothed20);
plot(xs(1:length(smoothed23)), smoothed23);

ylim([0 0.5]);