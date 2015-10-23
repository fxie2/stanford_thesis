% - Kartikey Asthana
% ---

%% Initialization
path(path,'../Polynomials')
path(path,'../IO')
clear all;
clc;

%% Parameters

% FR scheme
P = 2;                                                                      % Polynomial order
s = P+3;                                                                    % Order of correction polynomial
c1Arr = [-1/15; -5e-2; -1e-2; 0; 1e-2; 5e-2; 1e-1];                         % Scheme parameter

% Computational domain
h = 1;                                                                      % Width of element
nex = 1e3;                                                                  % No. of plot points
xpts = linspace(-1,1,nex)';                                                 % Location of plot points

% Aux
saveFlag = 1;                                                               % 1 to overwrite figure

%% Correction functions

% Declare arrays
nc1 = length(c1Arr);
gL0v = zeros(nex,nc1); gR0v = zeros(nex,nc1);
gL1v = zeros(nex,nc1); gR1v = zeros(nex,nc1);

for i = 1: nc1
    
    % c1
    c1 = c1Arr(i);
    
    % g0 functions
    theta = gLC1FR(P,h/2,c1,0);
    for k = 0: s
        gL0v(:,i) = gL0v(:,i) + theta(k+1)*(xpts.^k);
        gR0v(:,i) = gR0v(:,i) + ((-1)^k)* theta(k+1)*(xpts.^k);
    end
    
    % g1 functions
    theta = gLC1FR(P,h/2,c1,1);
    for k = 0: s
        gL1v(:,i) = gL1v(:,i) + theta(k+1)*(xpts.^k);
        gR1v(:,i) = gR1v(:,i) - ((-1)^k)* theta(k+1)*(xpts.^k);
    end
    
end

save 'cfs_C1FR.mat'

%% Post-processing

% Plot and marker styles
plotStyle = {'k-','r--','g-.','b:','m-','k-','r-'};
markerStyle = {'','','','','ms','ko','r^'};
%legendStyle = {'k-','r--','g-.','b:','ms-','k-o','r-^'};
mSkip = 50;
pH = zeros(nc1,1);

f1 = figure(1);
set(gcf,'Position',[0, 0, 1200, 700])
% 1,1
subplot(2,2,1)
for i = 1: nc1
    plot(xpts,gL0v(:,i),plotStyle{i},'LineWidth',2);
    set(gca,'Units','normalized','Position',[0.1 0.6 0.35 0.35])
    hold on;
    if i>4; scatter(xpts(1:mSkip:end),gL0v(1:mSkip:end,i),markerStyle{i}); end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{L_0}$','interpreter','latex');

% 1,2
subplot(2,2,2)
for i = 1: nc1
    plot(xpts,gR0v(:,i),plotStyle{i},'LineWidth',2);
    set(gca,'Units','normalized','Position',[0.55 0.6 0.35 0.35])
    hold on;
    if i>4; scatter(xpts(1:mSkip:end),gR0v(1:mSkip:end,i),markerStyle{i}); end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{R_0}$','interpreter','latex');

% 2,1
subplot(2,2,3)
for i = 1: nc1
    plot(xpts,gL1v(:,i),plotStyle{i},'LineWidth',2);
    set(gca,'Units','normalized','Position',[0.1 0.15 0.35 0.35])
    hold on;
    if i>4; scatter(xpts(1:mSkip:end),gL1v(1:mSkip:end,i),markerStyle{i}); end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{L_1}$','interpreter','latex');

%% Hacked figure to create legend handles
%{
subplot(2,3,6)
for i = 1: nc1
    pH(i) = plot(xpts,gR1v(:,i),legendStyle{i},'LineWidth',2);
    %hold off
    set(gca,'Units','normalized','Position',[1.55 0.15 0.35 0.35])
    hold on;
    %if i>4; scatter(xpts(1:mSkip:end),gR1v(1:mSkip:end,i),markerStyle{i}); end
end
%}
%%
% 2,2
subplot(2,2,4)
for i = 1: nc1
    pH(i) = plot(xpts,gR1v(:,i),plotStyle{i},'LineWidth',2);
    set(gca,'Units','normalized','Position',[0.55 0.15 0.35 0.35])
    hold on;
    if i>4; scatter(xpts(1:mSkip:end),gR1v(1:mSkip:end,i),markerStyle{i}); end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{R_1}$','interpreter','latex');


figure(1)
%%
leg = legend(pH,'$c_1 = -1/15$','$c_1 = -5\cdot10^{-2}$','$c_1 = -10^{-2}$','$c_1 = 0$','$c_1 = 10^{-2}$','$c_1 = 5\cdot10^{-2}$'...
    ,'$c_1 = 10^{-1}$','Orientation','horizontal');
legendlinestyles(leg,{'none','none','none','none','s','o','^'},{'-','--','-.',':','-','-','-'},{'k','r','g','b','m','k','r'})

%legend boxoff
set(leg,'box','off','interpreter','latex','Units','normalized','Position',[0.1 0.01 0.8 0.1])
set(gcf,'PaperPositionMode','auto')
if saveFlag; print -loose -depsc C1FR.eps;end
    %saveas(f1,'C1FR.eps','epsc'); end
 %set(gcf, 'Units','normalized', 'Position',[0.15 0.3 0.75 0.7])
