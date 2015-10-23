function tAE = postProc(teN, eNormNum, eNormEx, xPlot, uPlot, uEx, pSch)


%% Plot Style

pS = {'k-','r-.','b--'};
legVals = {'Exact','C1FR','DG'};

%% Solution plots

% Reshape arrays
np = size(xPlot,1) * size(xPlot,2);
xL = reshape(xPlot,np,1);
uL = reshape(uPlot,np,5);
uExL = reshape(uEx,np,5);

%{
figure(11)

if (strcmpi(pSch,'C1FR'))
    
    % Time tf/4
    subplot(2,2,1)
    plot(xL, uExL(:,2), pS{1}); hold on;
    plot(xL, uL(:,2), pS{2});
    xlabel('x'); ylabel('u(x)'); title('$t = t_f / 4$', 'interpreter', 'latex')
    xlim([-1,1]); ylim(1.5*[min(uExL(:,2)),max(uExL(:,2))]);
    
    % Time tf/2
    subplot(2,2,2)
    plot(xL, uExL(:,3), pS{1}); hold on;
    plot(xL, uL(:,3), pS{2});
    xlabel('x'); ylabel('u(x)'); title('$t = t_f / 2$', 'interpreter', 'latex')
    xlim([-1,1]); ylim(1.5*[min(uExL(:,3)),max(uExL(:,3))]);
    
    % Time 3tf/4
    subplot(2,2,3)
    plot(xL, uExL(:,4), pS{1}); hold on;
    plot(xL, uL(:,4), pS{2});
    xlabel('x'); ylabel('u(x)'); title('$t = 3 t_f / 4$', 'interpreter', 'latex')
    xlim([-1,1]); ylim(1.5*[min(uExL(:,4)),max(uExL(:,4))]);
    
    % Time tf
    subplot(2,2,4)
    plot(xL, uExL(:,5), pS{1}); hold on;
    plot(xL, uL(:,5), pS{2});
    xlabel('x'); ylabel('u(x)'); title('$t = t_f$', 'interpreter', 'latex')
    xlim([-1,1]); ylim(1.5*[min(uExL(:,5)),max(uExL(:,5))]);
    
else
    
    % Time tf/4
    subplot(2,2,1)
    plot(xL, uL(:,2), pS{3}); 
    
    % Time tf/2
    subplot(2,2,2)
    plot(xL, uL(:,3), pS{3});

    % Time 3tf/4
    subplot(2,2,3)
    plot(xL, uL(:,4), pS{3});
    
    % Time tf
    subplot(2,2,4)
    plot(xL, uL(:,5), pS{3});
    
    legend(leg,'Location',[0.91,0.6,0.1,0.1]);
    legend('boxoff')
    
end
%}

%% Plots of energies
height = 0.2;
legHeight = 0.01;
% Energy in 0th derivative ||u||_{1,0}
figure(20);
set(gcf,'Position',[1200, 0, 1200, 350])
set(gcf,'PaperPositionMode','auto')
subplot(1,2,1)

if (strcmpi(pSch,'C1FR'))
    
    semilogy(teN, eNormEx(1,:),pS{1},'LineWidth',2); hold on;
    semilogy(teN, eNormNum(1,:),pS{2},'LineWidth',2);
    xlabel('t'); ylabel('$\|u\;  \|^2_{(2,0)}$','interpreter','latex')
    %
    %hold on
else
    semilogy(teN, eNormNum(1,:),pS{3},'LineWidth',2);
    set(gca,'Units','normalized','Position',[0.1 height 0.35 0.75])
    %legend(leg)
end

% Energy in 1st derivative ||u||_{1,1}
subplot(1,2,2)

if (strcmpi(pSch,'C1FR'))
    %
    semilogy(teN, eNormEx(2,:),pS{1},'LineWidth',2); hold on;
    semilogy(teN, eNormNum(2,:),pS{2},'LineWidth',2);
    xlabel('t'); ylabel('$\|u\; \|^2_{(2,1)}$','interpreter','latex')
    %
    %hold on
else
    semilogy(teN, eNormNum(2,:),pS{3},'LineWidth',2);
    set(gca,'Units','normalized','Position',[0.55 height 0.35 0.75])
    leg = legend(legVals);
    set(leg,'Orientation','horizontal','box','off',...
        'interpreter','latex','Units','normalized',...
        'Position',[0.15 legHeight 0.8 0.1])
end


%% Error in energy norms

% Error terms
tAE = tAvgError(teN, eNormEx, eNormNum);

disp(['Error in energy norms for ',pSch,': 0th | 1st'])
disp(tAE)


