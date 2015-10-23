function plotFR(fig, tstr)
% Plotting routine
% ---
% fig - figure number start-index
% tstr - plot title
% -------------------------------------------------------------------------

global x u xbdry xedge h nelem xi P xex uex

npUnit = 100;                   % No. of plot points per unit spacing

% Exact solution
figure(fig);
pl1 = plot(xex,uex,'k--');
hold on;

% Numerical solution
xlist = reshape(x,(P+1)*nelem,1);
ulist = reshape(u,(P+1)*nelem,1);
scatter(xlist,ulist,'*r');
hold on;
for j = 1: nelem
    nelp = round(npUnit*h(j));
    xel = linspace(xedge(j),xedge(j+1),nelp);
    uel = zeros(nelp,1);
    for k = 1: nelp
        xik = -1 + (2/h(j))*(xel(k)-xedge(j));
        for p = 1: P+1
            uel(k) = uel(k) + u(p,j)*Lagrange(xi,p,xik);
        end
    end
    pl2 = plot(xel,uel,'-r');
end

xlim([xbdry(1), xbdry(2)]); %ylim([-0.2 1.2]);
title(tstr);
legend([pl1,pl2],{'u_{ex}','u'});
