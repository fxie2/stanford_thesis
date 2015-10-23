function sdOptPlot1(fig,eta,c,P)
% Plotting routine
% ---
% fig - figure number
% eta - optimization function
% c - c-VCJH range
% tstr - plot title
% -------------------------------------------------------------------------

%% Linear Search plot
figure(fig);
semilogy(c,eta,'k-');
title(['Variation of optimization function \eta with c for P = ',num2str(P)])
xlabel('c'); ylabel('\eta(c)');

