smth = 101;
xs = 0:50000/(2*pi):(1e7)+10000;

eccs0 = readmatrix("sim_0_ecc_dmc.csv");
eccs1 = readmatrix("sim_1_ecc_dmc.csv");

smoothed0 = movmean(eccs0, smth);
smoothed1 = movmean(eccs1, smth);

hold on
plot(xs(1:length(eccs0)), eccs0);
%plot(xs(1:length(smoothed0)), smoothed0);
plot(xs(1:length(eccs1)), eccs1);
%plot(xs(1:length(smoothed1)), smoothed1);