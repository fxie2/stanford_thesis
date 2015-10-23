function plotFR2(fig, tstr, t, uex, ind)
% Routine for plotting
% ---
% fig - figure number start-index
% tstr - plot title
% uex  - Exact solution
% ind - plot index
% -------------------------------------------------------------------------

global x u xbdry xedge h nelem xi P xex nex plotStyle

figure(fig);

% Exact solution
if (ind == 1 || ind ==5)
    xlist = reshape(xex,nex*nelem,1);
    ulist = reshape(uex,nex*nelem,1);
    pl1 = plot(xlist,ulist,plotStyle{1});
    hold on;
end

% Numerical solution
xlist = reshape(x,(P+1)*nelem,1);
ulist = reshape(u,(P+1)*nelem,1);
scatter(xlist,ulist,'*',plotStyle{1+ind}(1));
hold on;
for j = 1: nelem
    uel = zeros(nex,1);
    for k = 1: nex
        xik = -1 + (2/h(j))*(xex(k,j)-xedge(j));
        for p = 1: P+1
            uel(k) = uel(k) + u(p,j)*Lagrange(xi,p,xik);
        end
    end
    pl2 = plot(xex(:,j),uel,plotStyle{1+ind});
end

xlim([xbdry(1), xbdry(2)]);
ylim([-1.1,1.1]);
xlabel('x'); ylabel('u'); title(tstr);
