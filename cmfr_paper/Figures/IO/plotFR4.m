function plotFR4(fig, tstr, t, ind)
% Routine for plotting
% ---
% fig - figure number start-index
% tstr - plot title
% ind - plot index
% -------------------------------------------------------------------------

global u xbdry xedge h nelem xi P xex nex plotStyle

figure(fig);
hold on;

% Numerical solution
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
%ylim([-1,1]);
xlabel('x'); ylabel('u'); title(tstr);

switch ind
    case 1
        t2 = t;
        save('Buffer/pLegend.mat','pl2','t2');
    case 2
        t3 = t;
        pl3 = pl2;
        save('Buffer/pLegend.mat','pl3','t3','-append');
    case 3
        t4 = t;
        pl4 = pl2;
        save('Buffer/pLegend.mat','pl4','t4','-append');
    case 4
        t5 = t;
        pl5 = pl2;
        load('Buffer/pLegend.mat');
        legend([pl2,pl3,pl4,pl5],...
            ['DG    (t=',num2str(t2),')'],['OESFR (t=',num2str(t3),')'],...
            ['OFR   (t=',num2str(t4),')'],['c+    (t=',num2str(t5),')'])
end
