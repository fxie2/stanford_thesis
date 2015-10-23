% load figures, merge them into a single figure
% and modify the legend


f1 = hgload('eN0_high_k.fig');
set(f1, 'Visible', 'off')
tmpaxis1 = findobj(f1, 'Type','axes');

f2 = hgload('eN1_high_k.fig');
set(f2, 'Visible', 'off');
tmpaxis2 = findobj(f2, 'Type','axes');

h = figs2subplots([f1 f2], [1 2]);

set(h, 'Position',[1200, 0, 1200, 350])

set(gcf,'PaperPositionMode','auto')

saveFlag = 0;
if saveFlag; print -loose -depsc high_k.eps;end