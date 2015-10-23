function sdplot(fig,Cn,beta,k,tstr)
% Plotting routine
% ---
% fig - figure number start-index
% Cn - complex numerical wavespeed for the eigenmodes
% beta - relative weights for the eigenmodes
% k - wavenumber array
% tstr - plot title
% -------------------------------------------------------------------------

%% Plot styles
plotStyle = {'k-','r--','g--','b--.','r:','g:','b:','r-.','g-.','b-.'};
P = size(Cn,1) - 1;
legendCell = cell(P+1,1);

%% Dispersion
figure(fig+1);
for p = 1:P+1
    plot(k/(P+1),real(Cn(p,:)),plotStyle{p},'MarkerSize',2);
    hold on;
    legendCell{p} = ['Mode ',num2str(p)];
end
title(['Semi-discrete numerical dispersion for ',tstr])
xlabel('k/(P+1)'); ylabel('Real{C_n}');
legend(legendCell{1:P+1},1);
ylim([-2 3]);

%% Dissipation
figure(fig+2);
for p = 1:P+1
    plot(k/(P+1),imag(Cn(p,:)),plotStyle{p},'MarkerSize',2);
    hold on;
    legendCell{p} = ['Mode ',num2str(p)];
end
title(['Semi-discrete numerical dissipation for ',tstr])
xlabel('k/(P+1)'); ylabel('Imag{C_n}');
legend(legendCell{1:P+1},4);
ylim([-2 0.1]);

%% Energy distribution
figure(fig+3);
for p = 1:P+1
    plot(k/(P+1),beta(p,:),plotStyle{p},'MarkerSize',2);
    hold on;
    legendCell{p} = ['Mode ',num2str(p)];
end
title(['Semi-discrete energy distribution for ',tstr])
xlabel('k/(P+1)'); ylabel('\beta');
legend(legendCell{1:P+1},2);
ylim([0 1]);

