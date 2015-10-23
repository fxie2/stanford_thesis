function fdplotPhy(fig,CnList,betaList,k,tstr)
% Plotting routine
% ---
% fig - figure number start-index
% CnList - complex numerical wavespeed for the eigenmodes for optimal, RK,
% Semi-discrete schemes
% beta List- relative weights for the eigenmodes for the three schemes
% k - wavenumber array
% tstr - plot title
% -------------------------------------------------------------------------

%% Extract data sets and define plot styles
plotStyle = {'r-','g-','b-'};
P = size(CnList,1) - 1;
ksize = size(k,2);
Cn = CnList(1,1:ksize);
CnRK = CnList(1,ksize+1:2*ksize);
Cnsd = CnList(1,2*ksize+1:3*ksize);
beta = betaList(1,1:ksize);
betaRK = betaList(1,ksize+1:2*ksize);
betasd = betaList(1,2*ksize+1:3*ksize);

%% Dispersion
figure(fig+1);
plot(k/(P+1),real(Cnsd(1,:)),plotStyle{1},'MarkerSize',2);
hold on;
plot(k/(P+1),real(CnRK(1,:)),plotStyle{2},'MarkerSize',2);
plot(k/(P+1),real(Cn(1,:)),plotStyle{3},'MarkerSize',2);
title(tstr)
xlabel('k/(P+1)'); ylabel('Real{C_n}');
ylim([0 1.5]);
legend('Semi-discrete','RK','ORK')

%% Dissipation
figure(fig+2);
plot(k/(P+1),imag(Cnsd(1,:)),plotStyle{1},'MarkerSize',2);
hold on;
plot(k/(P+1),imag(CnRK(1,:)),plotStyle{2},'MarkerSize',2);
plot(k/(P+1),imag(Cn(1,:)),plotStyle{3},'MarkerSize',2);
title(tstr)
xlabel('k/(P+1)'); ylabel('Imag{C_n}');
ylim([-2 0.1]);
legend('Semi-discrete','RK','ORK')

%% Energy distribution
% figure(fig+3);
% plot(k/(P+1),betasd(1,:),plotStyle{1},'MarkerSize',2);
% hold on;
% plot(k/(P+1),betaRK(1,:),plotStyle{2},'MarkerSize',2);
% plot(k/(P+1),beta(1,:),plotStyle{3},'MarkerSize',2);
% title(tstr)
% xlabel('k/(P+1)'); ylabel('\beta');
% ylim([0 1]);
% legend('Semi-discrete','RK','ORK')
