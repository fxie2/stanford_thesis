% - Kartikey Asthana
% ---

%% Initialization
path(path,'../Polynomials')
path(path,'../IO')
clear all;
clc;

%% Parameters

% FR scheme
P = 3;                                                                      % Polynomial order
s = P+5;                                                                    % Order of correction polynomial
c1Arr = [0; 1e-2];                                                          % Scheme parameter
c2Arr = [0; 1e-3];                                                          % Scheme parameter

% Computational domain
h = 1;                                                                      % Width of element
nex = 1e3;                                                                  % No. of plot points
xpts = linspace(-1,1,nex)';                                                 % Location of plot points

% Aux
saveFlag = 1;                                                               % 1 to overwrite figure

%% Correction functions

% Declare arrays
nc1 = length(c1Arr);
nc2 = length(c2Arr);
gL0v = zeros(nex,nc1,nc2); gR0v = zeros(nex,nc1,nc2);
gL1v = zeros(nex,nc1,nc2); gR1v = zeros(nex,nc1,nc2);
gL2v = zeros(nex,nc1,nc2); gR2v = zeros(nex,nc1,nc2);


for j = 1: nc2
    for i = 1: nc1
        
        % c1, c2
        c2 = c2Arr(j);
        c1 = c1Arr(i);
        
        % g0 functions
        theta = gLC2FR(P,h/2,[c1;c2],0);
        for k = 0: s
            gL0v(:,i,j) = gL0v(:,i,j) + theta(k+1)*(xpts.^k);
            gR0v(:,i,j) = gR0v(:,i,j) + ((-1)^k)* theta(k+1)*(xpts.^k);
        end
        
        % g1 functions
        theta = gLC2FR(P,h/2,[c1;c2],1);
        for k = 0: s
            gL1v(:,i,j) = gL1v(:,i,j) + theta(k+1)*(xpts.^k);
            gR1v(:,i,j) = gR1v(:,i,j) - ((-1)^k)* theta(k+1)*(xpts.^k);
        end
        
        % g2 functions
        theta = gLC2FR(P,h/2,[c1;c2],2);
        for k = 0: s
            gL2v(:,i,j) = gL2v(:,i,j) + theta(k+1)*(xpts.^k);
            gR2v(:,i,j) = gR2v(:,i,j) + ((-1)^k)* theta(k+1)*(xpts.^k);
        end
    end
end

save 'cfs_C2FR.mat'

%% Post-processing

% Plot and marker styles
plotStyle = {'k-','g:','r--','b-.'};
markerStyle = {'','','rs','b^'};
mSkip = 50;
pH = zeros(nc1*nc2,1);

f1 = figure(1);
set(gcf,'Position',[1200, 0, 1200, 800])

% 1,1
subplot(3,2,1)
for j = 1: nc2
    for i = 1: nc1
        pH(nc1*(j-1)+i) = plot(xpts,gL0v(:,i,j),plotStyle{nc1*(j-1)+i},'LineWidth',2);
        hold on;
        if nc1*(j-1)+i>2; scatter(xpts(1:mSkip:end),gL0v(1:mSkip:end,i,j),markerStyle{nc1*(j-1)+i}); end
    end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{L_0}$','interpreter','latex');

% 1,2
subplot(3,2,2)
for j = 1: nc2
    for i = 1: nc1
        pH(nc1*(j-1)+i) = plot(xpts,gR0v(:,i,j),plotStyle{nc1*(j-1)+i},'LineWidth',2);
        hold on;
        if nc1*(j-1)+i>2; scatter(xpts(1:mSkip:end),gR0v(1:mSkip:end,i,j),markerStyle{nc1*(j-1)+i}); end
    end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{R_0}$','interpreter','latex');

% 2,1
subplot(3,2,3)
for j = 1: nc2
    for i = 1: nc1
        pH(nc1*(j-1)+i) = plot(xpts,gL1v(:,i,j),plotStyle{nc1*(j-1)+i},'LineWidth',2);
        hold on;
        if nc1*(j-1)+i>2; scatter(xpts(1:mSkip:end),gL1v(1:mSkip:end,i,j),markerStyle{nc1*(j-1)+i}); end
    end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{L_1}$','interpreter','latex');


% 2,2
subplot(3,2,4)
for j = 1: nc2
    for i = 1: nc1
        pH(nc1*(j-1)+i) = plot(xpts,gR1v(:,i,j),plotStyle{nc1*(j-1)+i},'LineWidth',2);
        hold on;
        if nc1*(j-1)+i>2; scatter(xpts(1:mSkip:end),gR1v(1:mSkip:end,i,j),markerStyle{nc1*(j-1)+i}); end
    end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{R_1}$','interpreter','latex');

% 3,1
subplot(3,2,5)
for j = 1: nc2
    for i = 1: nc1
        pH(nc1*(j-1)+i) = plot(xpts,gL2v(:,i,j),plotStyle{nc1*(j-1)+i},'LineWidth',2);
        hold on;
        if nc1*(j-1)+i>2; scatter(xpts(1:mSkip:end),gL2v(1:mSkip:end,i,j),markerStyle{nc1*(j-1)+i}); end
    end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{L_2}$','interpreter','latex');

% 3,2
subplot(3,2,6)
for j = 1: nc2
    for i = 1: nc1
        pH(nc1*(j-1)+i) = plot(xpts,gR2v(:,i,j),plotStyle{nc1*(j-1)+i},'LineWidth',2);
        hold on;
        if nc1*(j-1)+i>2; scatter(xpts(1:mSkip:end),gR2v(1:mSkip:end,i,j),markerStyle{nc1*(j-1)+i}); end
    end
end
xlabel('$\xi$','interpreter','latex');
ylabel('$g_{R_2}$','interpreter','latex');

leg = legend(pH,['$c_1 = 0$',10,'$c_2 = 0$'],['$c_1 = 10^{-2}$',10,'$c_2 = 0$'],...
    ['$c_1 = 0$',10,'$c_2 = 10^{-3}$'],['$c_1 = 10^{-2}$',10,'$c_2 = 10^{-3}$'],'interpreter','latex','Orientation','horizontal');
set(leg,'box','off','interpreter','latex','Units','normalized','Position',[0.15 -0.02 0.8 0.1])

legendlinestyles(leg,{'none','none','s','^'},{'-',':','--','-.'},{'k','g','r','b'})
% plotStyle = {'k-','g:','r--','b-.'};
% markerStyle = {'','','rs','b^'};

%if saveFlag; saveas(f1,'C2FR','fig'); end
set(gcf,'PaperPositionMode','auto')
if saveFlag; print -loose -depsc C2FR.eps;end
