% FR1D_driver - Flux Reconstruction in 1D - driver function
% ---
% Time marches the hyperbolic scalar conservation law dudt + dfdx = s
% where f = f1 + f2,
% f1 is a possibly non-linear function of u
% f2 is a possibly non-linear function of du/dx
% s is a possibly non-linear function of u
%
% - Periodic boundary conditions
% - Correction functions: (1) ESFR, (2) g with known zeros
% - Time stepping: RK4
%
% - Kartikey Asthana
% ---
% function [xL, tQ, uQv, teN, eNormNum, eNormEx] = FluxRecons(alpArr,k0F)
function [teN, errNorm01] = FluxRecons(FRsch, P_in, h_in, alpArr, dt, tf, efreq, plotFreq)
%
% Input:
% FRsch     - 1 for DG, 2 for OESFR, 3 for OFR
% P_in      - polynomial order
% h_in      - grid spacing
% alpArr    - [alp10, alp11; alp20, alp21] upwinding coefficients where:
%             alp10 corresponds to 0th order corrections for 1st order flux
%             alp11 corresponds to 1st order corrections for 1st order flux
%             alp20 corresponds to 0th order corrections for 2nd order flux
%             alp21 corresponds to 1st order corrections for 2nd order flux
% dt        - Time step
% tf        - Final time
% efreq     - frequency of norm evaluation
% plotFreq  - Frequency of plotting
%
% Output:
% teN       - Array of times at which energy norms were evaluated
% errNorm01 - Error norms in 0th and 1st derivative
%
%   errNorm01 = [err0; err1] where
%       err0   - Error norm for the 0th derivative
%                   = \sum_{n=1}^{n_{el}} \int_{x_n}^{x_{n=1}} (u - u^\delta)^2/2 dx
%       err1   - Error norm for the 1st derivative
%                   = \sum_{n=1}^{n_{el}} \int_{x_n}^{x_{n=1}} (dudx - du^\deltadx)^2/2 dx
%
% efreq     - Frequency of energy norm evaluation

%% Global variables and parameters

global f1u f2u su a1u a2u alp1 alp2 alpu nu ...
    x u xbdry xedge h nelem xi P  ...
    xex nex plotStyle ...
    phiL phiR D gpLv gpRv ...
    k0 kN

% Conservation law
f1u = @(u) u;                                                               % Convective flux
a1u = @(u) 1;                                                               % Analytical wavespeed for 1st order flux
nu = 1e-2;                                                                  % Diffusion coefficient
f2u = @(dudx) -nu*dudx;                                                     % Diffusive flux
a2u = @(dudx) -nu;                                                          % Analytical wavespeed for 2nd order flux
su = @(u) 0*u;                                                              % Source term
alp1 = alpArr(1);                                                           % Upwinding coefficient for 1st order flux
alpu = 1;                                                                   % Upwinding coefficient for u
alp2 = alpArr(2);                                                           % Upwinding coefficient for 2nd order flux

% FR scheme
P = P_in;                                                                   % Polynomial order
switch (FRsch)
    case 1
        gpL = @(xi) VCJH_p(P,0,xi,-1);                                      % Left bdry correction
        gpR = @(xi) VCJH_p(P,0,xi,1);                                       % Right bdry correction
    case 2
        load(['~/ACL/DRP_FR/Routines/Semi_discrete_DRP_lt/cVCJH_P',num2str(P),'.mat'],'cOpt');
        gpL = @(xi) VCJH_p(P,cOpt,xi,-1);
        gpR = @(xi) VCJH_p(P,cOpt,xi,1);
    case 3
        load(['~/ACL/DRP_FR/Routines/Semi_discrete_DRP_lt/zgLB_P',num2str(P),'.mat'],'zgLB');
        gLGrid = [-1;zgLB;1];
        gRGrid = sort([1;-zgLB;-1],'ascend');
        gpL = @(xi) Lagrange_p(gLGrid,1,xi)/Lagrange(gLGrid,1,-1);
        gpR = @(xi) Lagrange_p(gRGrid,P+2,xi)/Lagrange(gRGrid,P+2,1);
end

% Computational domain
xbdry = [-10, 10];                                                          % Domain
nelem = round((xbdry(2)-xbdry(1)) / h_in);                                  % No. of elements
xedge = linspace(xbdry(1),xbdry(2),nelem+1);                                % Element edges
xi = GLpoints(P);                                                           % Soln. points

% Initial condition
u0 = @(x) init_cond(x,0);                                                   % Initial condition

% Auxiliary
k0 = 2*pi/(xbdry(2)-xbdry(1));                                              % Base wavenumber for the domain
kN = (P+1)*pi;                                                              % Nyquist frequency for P=P_in, h=1
nex = round(10 * h_in);                                                     % No. of pts. for plotting in each element
tDelay = 0.01;                                                              % Delay between plots

%% Grid setup

% Grid spacing and global coordinates
h = zeros(1,nelem); x = zeros(P+1,nelem);
for j = 1: nelem
    h(j) = xedge(j+1) - xedge(j);
    x(:,j) = 0.5*(1-xi)*xedge(j) + 0.5*(1+xi)*xedge(j+1);
end

% Grid for interpolated solution
xex = zeros(nex,nelem);
for j = 1: nelem
    xex(:,j) = linspace(xedge(j),xedge(j+1),nex);
end

%% Pre-processing

% Lagrange function values at edges
phiL = zeros(P+1,1); phiR = zeros(P+1,1);
for p = 1: P+1
    phiL(p) = Lagrange(xi,p,-1);
    phiR(p) = Lagrange(xi,p,1);
end

% Pth order Lagrange function derivatives at solution points
D = zeros(P+1,P+1);
for k = 1: P+1
    for p = 1: P+1
        D(k,p) = Lagrange_p(xi,p,xi(k));
    end
end

% Correction function derivatives at solution points
gpLv = zeros(P+1,1); gpRv = zeros(P+1,1);
for p = 1: P+1
    gpLv(p) = gpL(xi(p));
    gpRv(p) = gpR(xi(p));
end

% Exact solution
% uex = @(x,t) exact_sol(x, t, 0);
% uex_d = @(x,t) exact_sol(x, t, 1);
uexV = @(t) exact_sol(x, t, 0);
uex_dV = @(t) (2/h_in)*D*exact_sol(x, t, 0);
uex = @(x,t) eval_poly_u(uexV(t), x);
uex_d = @(x,t) eval_poly_u(uex_dV(t), x);

% Style
plotStyle = {'k-','r-','g-','b-','m-','g-'};

%% FR update

% Setup run
u = u0(x); t = zeros(1,1e5);
cnt = 1; iter = 1;

% Setup norm arrays
errNorm01 = zeros(2,1e5);
teN = zeros(1,1e5);

% Update
while 1
    
    % Error norm
%     if (mod(iter,efreq)==0)
%         teN(cnt) = t(iter);
%         [errNorm01(1,cnt), errNorm01(2,cnt)] = errNorm(uex, uex_d, t(iter));
%         cnt = cnt+1;
%     end
   
    % Stage 1
    ui = u;
    Fp0 = FR(u);
           
    % Plots of solution
    if (mod(iter,plotFreq) == 0)
        plotFR2p(1,'DG',t(iter),uex(xex,t(iter)),1);
        hold off;
        pause(tDelay);
    end
    
    u = ui + (dt/2)*Fp0;
    
    % Stage 2
    Fp1 = FR(u);
    u = ui + (dt/2)*Fp1;
    
    % Stage 3
    Fp2 = FR(u);
    u = ui + dt*Fp2;
    
    % Stage 4
    Fp3 = FR(u);
    u = ui + (dt/6)*(Fp0 + 2*Fp1 + 2*Fp2 + Fp3);
    
    % Update
    iter = iter + 1;
    t(iter) = t(iter-1) + min(dt);
    
    % Check for completion
    if (t(iter)>tf); break; end
    
end

%% Truncate zero entries

[errNorm01(1,cnt), errNorm01(2,cnt)] = errNorm(uex, uex_d, t(iter-1));
cnt = 1;

teN = teN(1:cnt);
errNorm01 = errNorm01(:,1:cnt);
