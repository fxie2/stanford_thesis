% Order of accuracy for C1FR
% ---
% - Kartikey Asthana
% ---

%% Initialization
path(path,'../Polynomials')
path(path,'../IO')
path(path,'../Error')
clear all;
clc;
close all

global f1u f2u su a1u a2u alp10 alp11 alp20 alp21 alpu nu ...
    x u xbdry xedge  nelem xi   ...
    xex nex plotStyle ...
    phi0L phi0R phi1L phi1R D gL0pv gR0pv gL1pv gR1pv ...
    k0 kN uref plotFlag upref
%% Parameters

% Upwinding coefficients - purely upwinded
alpArrC1FR = [0,0;0,0];     % coefficients for C1FR
alpArrDG = [0;0];           % coefficients for DG

% Polynomial order array                  
P = 1:3;
nP = length(P);               % No. of P values

% Grid spacing

xbdry = [0, 1.3e2];
xdom = xbdry(2)-xbdry(1);

nh = 4;                     % No. of h values
maxNel = 34;
minNel = 17;
nel = floor((2).^(linspace(log(minNel)/log(2),log(maxNel)/log(2),nh)));      % Number of elements
h = xdom./nel;                % Grid spacing

% Time step
CFL = 0.01;
tf = 2*xdom;%1e5*dt;            % Final time

% Frequencies
efreq = 1;                  % Frequency of norm evaluation
plotFreq = 1e7;             % Frequency of plotting

% Flags

plotFlag = 0;

%% Optimal c1byJsq values

% Define stability constraints
c1byJsq_lim = -2./(P.*(P+1).*(2*P+1));

% Set c1 to zero
c1byJsq = zeros(1,nP);

%% Error tables

% Declare error arrays
err_eN0_DG = zeros(nP,nh);
err_eN1_DG = zeros(nP,nh);
err_eN0_C1FR = zeros(nP,nh);
err_eN1_C1FR = zeros(nP,nh);
%%{
% Generate table
simCount = 1;
nSims = nP*nh;
for i = 1:nP
    for j = 1: nh
        if plotFlag; figure; end
        % Status
        
        %dt_temp = CFL*h(j)/(P(i)+1);
        dt_temp = CFL*h(j)*DGCFL(P(i),1);
        n_dt = ceil(xdom/dt_temp);
        dt = xdom/(n_dt);
        
        display(sprintf(['\n\n-*-*-*-*-*-*-*-*-*-*-*-'...
            '*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n\n'...
            'Simulation %i of %i'],simCount,nSims))
        simCount = simCount + 1;
        display(sprintf(['C1FR:\n P = %i, nel = %i,\n'...
            'dt = %.4e, n_dt = %i'],P(i),nel(j),dt,n_dt))

        % C1FR case
        outerLoopTime = tic;
        c1frTime = tic;
        if plotFlag; subplot(1,2,1); end
        [errNorm01] = CmFluxRecons(P(i), h(j), c1byJsq(i), alpArrC1FR, dt, tf, efreq, plotFreq);
        errNorm01
        display(sprintf('C1FR took %.4f',toc(c1frTime)))
        %tAE = tAvgError(teN, errNorm01, normConst01);
        tAE = (errNorm01);%/sqrt(nel(j)*(P(i)+1));
        err_eN0_C1FR(i,j) = tAE(1);
        err_eN1_C1FR(i,j) = tAE(2);
        %figure(100*i+10*j+1)
        if plotFlag
            hold on
            plot(reshape(x,[1 numel(x)]),reshape(uref,[1 numel(uref)]),'r-o')
            title(['CMFR, P=' num2str(P(i)) ', h = ' num2str(h(j))])
        end
        % DG baseline
        if plotFlag; subplot(1,2,2); end
        display('DG')
        dgTime = tic;
        [errNorm01] = FluxRecons(1, P(i), h(j), alpArrDG, dt, tf, efreq, plotFreq);
        errNorm01
        display(sprintf('DG took %.4f',toc(dgTime)))
        %tAE = tAvgError(teN, errNorm01, normConst01);
        tAE = (errNorm01);%/sqrt(nel(j)*(P(i)+1));
        err_eN0_DG(i,j) = tAE(1);
        err_eN1_DG(i,j) = tAE(2);
        endTime = toc(outerLoopTime);
        %figure(100*i+10*j+2)
        display(sprintf(['Total took ' num2str(endTime) ' s\n\n']))
        if plotFlag
            hold on
            plot(reshape(x,[1 numel(x)]),reshape(uref,[1 numel(uref)]),'r-o')
            title(['DG, P=' num2str(P(i)) ', h = ' num2str(h(j))])
        end
    end
end
%}
%% Save data
caseName = sprintf('case_P_%i_nh_%i_C1_adv_CFL_%.2e_dx_%.1e.mat',nP,nh,CFL,xbdry(2));
save(caseName)

%% Plots and order of accuracy
%%{
% Load data
%load 'large_case.mat'
%load('case_P_5_nh_4_C1_adv_fine_dt')
% Plot styles
pS = {'r-.','b--'};
mS = {'r^','bs'};
saveFlag = 0;
% Linear regression matrix
A = [log(nel'), ones(nh,1)];

% Order of accuracy vector
O_eN0_C1FR = zeros(nP,1);
O_eN0_DG = zeros(nP,1);
O_eN1_C1FR = zeros(nP,1);
O_eN1_DG = zeros(nP,1);

allErrors = [err_eN0_C1FR;err_eN0_DG;err_eN1_C1FR;err_eN1_DG];
miny = min(min(allErrors))*0.5;
maxy = max(max(allErrors))*2;

height = 0.2;
legHeight = 0.01;

for i = 1: nP
    
    figure(i)
    set(gcf,'Position',[1200, 0, 1200, 350])
    set(gcf,'PaperPositionMode','auto')
    % (1) Error in 0th derivative
    subplot(1,2,1);
    
    set(gca,'Units','normalized','Position',[0.1 height 0.35 0.75])
    
    % Linear fits
    fitC1FR = A \ log(err_eN0_C1FR(i,:)');
    fitDG = A \ log(err_eN0_DG(i,:)');
    nelL = linspace(nel(1),nel(end),100);
    fit1 = loglog(nelL, exp(fitC1FR(2)) * nelL.^fitC1FR(1), pS{1}); hold on;
    fit2 = loglog(nelL, exp(fitDG(2)) * nelL.^fitDG(1), pS{2});
    set([fit1,fit2],'LineWidth',2)
    grid on
    % Store rate of convergence
    O_eN0_C1FR(i) = -fitC1FR(1);
    O_eN0_DG(i) = -fitDG(1);
    
    % Scatter plot of error values
    pl1 = scatter(nel, err_eN0_C1FR(i,:), mS{1});
    pl2 = scatter(nel, err_eN0_DG(i,:), mS{2});
    set([pl1,pl2],'Linewidth', 2)
    % Annotate
    xlabel('$n_{el}$','interpreter','latex');
    ylabel('$e_{(2,0)}$','interpreter','latex');
    xlim([nel(1)/2,2*nel(end)])
    ylim([miny maxy])
    %leg1 = legend([pl1,pl2],'C1FR','DG');
    %set(leg1,'Interpreter','latex')
    set(gca,'Xtick',nel,'XtickLabel',nel)
    % (2) Error in 1st derivative
    subplot(1,2,2);
    
    set(gca,'Units','normalized','Position',[0.55 height 0.35 0.75])
    
    % Linear fits
    fitC1FR = A \ log(err_eN1_C1FR(i,:)');
    fitDG = A \ log(err_eN1_DG(i,:)');
    nelL = linspace(nel(1),nel(end),100);
    fit1 = loglog(nelL, exp(fitC1FR(2)) * nelL.^fitC1FR(1), pS{1}); hold on;
    fit2 = loglog(nelL, exp(fitDG(2)) * nelL.^fitDG(1), pS{2});
    set([fit1,fit2],'LineWidth', 2)
    % Store rate of convergence
    O_eN1_C1FR(i) = -fitC1FR(1);
    O_eN1_DG(i) = -fitDG(1);
    
    % Scatter plot of error values
    pl1 = scatter(nel, err_eN1_C1FR(i,:), mS{1});
    pl2 = scatter(nel, err_eN1_DG(i,:), mS{2});
    set([pl1,pl2],'Linewidth', 2)
    grid on
    
    
    % Annotate
    xlabel('$n_{el}$','interpreter','latex');
    ylabel('$e_{(2,1)}$','interpreter','latex');
    xlim([nel(1)/2,2*nel(end)])
    ylim([miny maxy])
    
    %set(leg2,'Interpreter','latex')
    set(gca,'Xtick',nel,'XtickLabel',nel)
    
    set(gcf,'PaperPositionMode','auto')
    
    
    
    leg = legend([pl1,pl2,fit1,fit2],'C1FR','DG','C1FR fit','DG fit');%,'Orientation','Horizontal',-1);
    %legendlinestyles(leg,{'^','s'},{'-.','--'},{'k','b'})
    set(leg,'Orientation','horizontal','box','off','interpreter','latex','Units','normalized','Position',[0.1 legHeight 0.8 0.1])
    if saveFlag; print('-loose','-depsc',sprintf('P_%i.eps',i));end
end


O_eN0_C1FR
O_eN0_DG

O_eN1_C1FR
O_eN1_DG

%}

