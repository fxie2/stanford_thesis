function sdplotPhy(fig,Cn,beta,k,tstr,col)
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

%% Dispersion
figure(fig+1);
p = 1;
if (col==1); plotStyle{p} = '-r'; end
if (col==2); plotStyle{p} = '-g'; end
if (col==3); plotStyle{p} = '-b'; end
plot(k/(P+1),real(Cn(p,:)),plotStyle{p},'MarkerSize',2);
hold on;

title(['Semi-discrete numerical dispersion for ',tstr])
xlabel('k/(P+1)'); ylabel('Real{C_n}');
ylim([0 1.5]);

%% Dissipation
figure(fig+2);
p = 1;
if (col==1); plotStyle{p} = '-r'; end
if (col==2); plotStyle{p} = '-g'; end
if (col==3); plotStyle{p} = '-b'; end
plot(k/(P+1),imag(Cn(p,:)),plotStyle{p},'MarkerSize',2);
hold on;
title(['Semi-discrete numerical dissipation for ',tstr])
xlabel('k/(P+1)'); ylabel('Imag{C_n}');
ylim([-2 0]);

%% Energy distribution
figure(fig+3);
p = 1;
if (col==1); plotStyle{p} = '-r'; end
if (col==2); plotStyle{p} = '-g'; end
if (col==3); plotStyle{p} = '-b'; end
plot(k/(P+1),beta(p,:),plotStyle{p},'MarkerSize',2);
hold on;

title(['Semi-discrete energy distribution for ',tstr])
xlabel('k/(P+1)'); ylabel('\beta');
ylim([0 1]);

