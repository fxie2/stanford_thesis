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
    k0 kN
%% Parameters

% Upwinding coefficients - purely upwinded
alpArrC1FR = [0,0;0,0];     % coefficients for C1FR
alpArrDG = [0;0];           % coefficients for DG

% Polynomial order array
nP = 3;                     % No. of P values
P = 1:nP;

% Grid spacing
nh = 1;                     % No. of h values
nel = 20* 2.^(0:nh-1);      % Number of elements
h = 20./nel;                % Grid spacing

% Time step
dt = 1e-2;                 % Time step
tf = 1e4*dt;                  % Final time

% Frequencies
efreq = 1;                  % Frequency of norm evaluation
plotFreq = 1e7;             % Frequency of plotting

% Flags
saveFlag = 0;

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
for i = 1: nP
    for j = 1: nh
        
        % Status
        disp(['P = ',num2str(P(i)),' ; nel = ',num2str(nel(j))]);
        
        % C1FR case
        [teN, errNorm01,normConst01] = CmFluxRecons(P(i), h(j), c1byJsq(i), alpArrC1FR, dt, tf, efreq, plotFreq);
        %tAE = tAvgError(teN, errNorm01, normConst01);
        tAE = errNorm01;
        err_eN0_C1FR(i,j) = tAE(1);
        err_eN1_C1FR(i,j) = tAE(2);
        
        % DG baseline
        %[teN, errNorm01] = FluxRecons(1, P(i), h(j), alpArrDG, dt, tf, efreq, plotFreq);
        %tAE = tAvgError(teN, errNorm01, normConst01);
        tAE = errNorm01;
        err_eN0_DG(i,j) = tAE(1);
        err_eN1_DG(i,j) = tAE(2);
    end
end
%}
%% Save data
%save 'data_order_accuracy_C1_0_adv.mat'

%% Plots and order of accuracy
%%{
% Load data
%load 'data_order_accuracy_C1_0_adv.mat'

% Plot styles
pS = {'r-.','b--'};
mS = {'r^','bs'};

% Linear regression matrix
A = [log(nel'), ones(nh,1)];

% Order of accuracy vector
O_eN0_C1FR = zeros(nP,1);
O_eN0_DG = zeros(nP,1);
O_eN1_C1FR = zeros(nP,1);
O_eN1_DG = zeros(nP,1);

for i = 1: nP
    
    figure(i)
    set(gcf,'Position',[0, 0, 1200, 350])
    set(gcf,'PaperPositionMode','auto')
    % (1) Error in 0th derivative
    subplot(1,2,1);
    
    % Linear fits
    fitC1FR = A \ log(err_eN0_C1FR(i,:)');
    fitDG = A \ log(err_eN0_DG(i,:)');
    nelL = linspace(nel(1),nel(end),100);
    loglog(nelL, exp(fitC1FR(2)) * nelL.^fitC1FR(1), pS{1}); hold on;
    loglog(nelL, exp(fitDG(2)) * nelL.^fitDG(1), pS{2});
    
    % Store rate of convergence
    O_eN0_C1FR(i) = -fitC1FR(1);
    O_eN0_DG(i) = -fitDG(1);
    
    % Scatter plot of error values
    pl1 = scatter(nel, err_eN0_C1FR(i,:), mS{1});
    pl2 = scatter(nel, err_eN0_DG(i,:), mS{2}); 
        
    % Annotate
    xlabel('$n_{el}$','interpreter','latex');
    ylabel('$e_{(1,0)}$','interpreter','latex');
    xlim([nel(1)/2,2*nel(end)])
    leg1 = legend([pl1,pl2],'C1FR','DG');
    set(leg1,'Interpreter','latex')
    set(gca,'Xtick',nel,'XtickLabel',nel)
    % (2) Error in 1st derivative
    subplot(1,2,2);
    
    % Linear fits
    fitC1FR = A \ log(err_eN1_C1FR(i,:)');
    fitDG = A \ log(err_eN1_DG(i,:)');
    nelL = linspace(nel(1),nel(end),100);
    loglog(nelL, exp(fitC1FR(2)) * nelL.^fitC1FR(1), pS{1}); hold on;
    loglog(nelL, exp(fitDG(2)) * nelL.^fitDG(1), pS{2});
    
    % Store rate of convergence
    O_eN1_C1FR(i) = -fitC1FR(1);
    O_eN1_DG(i) = -fitDG(1);
    
    % Scatter plot of error values
    pl1 = scatter(nel, err_eN1_C1FR(i,:), mS{1});
    pl2 = scatter(nel, err_eN1_DG(i,:), mS{2}); 
        
    % Annotate
    xlabel('$n_{el}$','interpreter','latex');
    ylabel('$e_{(1,1)}$','interpreter','latex');
    xlim([nel(1)/2,2*nel(end)])
    leg2 = legend([pl1,pl2],'C1FR','DG');
    set(leg2,'Interpreter','latex')
    set(gca,'Xtick',nel,'XtickLabel',nel)
    if saveFlag; print('-loose','-depsc',sprintf('P_%i.eps',i));end
end

O_eN0_C1FR
O_eN0_DG

%}

%save 'order_order_accuracy_C1_0.mat'
