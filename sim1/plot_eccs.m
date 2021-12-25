eccs0 = readmatrix("sim_0_ecc_dmc.csv");
eccs1 = readmatrix("sim_1_ecc_dmc.csv");
eccs2 = readmatrix("sim_2_ecc_dmc.csv");
eccs3 = readmatrix("sim_3_ecc_dmc.csv");

smth = 1001;

smoothed0 = movmean(eccs0, smth);
smoothed1 = movmean(eccs1, smth);
smoothed2 = movmean(eccs2, smth);
smoothed3 = movmean(eccs3, smth);

xs = 0:50000/(2*pi):(2e8)+10000;

plot(xs, smoothed0, xs, smoothed1, xs, smoothed2, xs, smoothed3);
ylim([0 0.5])
legend(["sim0" "sim1" "sim2" "sim3"]);