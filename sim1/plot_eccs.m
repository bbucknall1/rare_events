N = 5;
smth = 1001;
xs = 0:50000/(2*pi):(2e8)+10000;
L = length(xs);

eccs = zeros(N+1, L);
smoothed = zeros(N+1, L);

hold on;
for i = 0:N
    eccs(i+1,:) = readmatrix("sim_"+num2str(i)+"_ecc_dmc.csv");
    smoothed(i+1,:) = movmean(eccs(i+1,:), smth);
    plot(xs, smoothed(i+1,:));
end
ylim([0 0.5]);