% CmFR1D - Cm Flux Reconstruction in 1D
% ---
% Time marches the hyperbolic scalar conservation law dudt + df1(u)dx + df2(dudx)dx = s
% Time stepping: RK4
%
% - Kartikey Asthana
% ---
function [teN, eNormNum, eNormEx, xPlot, tPlot, uPlot, uEx] = CmFluxRecons(P_in, h_in, c1byJsq, alpArr, k0F, efreq, plotFreq)
%
% Input:
% P_in      - polynomial order
% h_in      - grid spacing
% c1byJsq   - normalized C1FR parameter
% alpArr    - [alp10, alp11; alp20, alp21] upwinding coefficients where:
%             alp10 corresponds to 0th order corrections for 1st order flux
%             alp11 corresponds to 1st order corrections for 1st order flux
%             alp20 corresponds to 0th order corrections for 2nd order flux
%             alp21 corresponds to 1st order corrections for 2nd order flux
% k0F       - Fraction defined as
%             wavenumber of initial condition / Nyquist wavenumber
% efreq     - frequency of norm evaluation
% plotFreq  - Frequency of plotting
%
% Output:
% xPlot     - linear grid for plotting
% tPlot     = [0; t_f/4; t_f/2; 3t_f/4; tf]: Array of times at which a
%             snapshot of solution was recorded
% uPlot     = [u(0); u(t_f/4); u(t_f/2); u(3t_f/4); u(tf)]: Array of
%               solution snapshots
% uEx       = [uex(0); uex(t_f/4); uex(t_f/2); uex(3t_f/4); uex(tf)]: Array of
%               exact solution at every quarter
% teN       - Array of times at which energy norms were evaluated
% eNormNum  - Energy norms of numerical solution
% eNormEx   - Energy norms of exact solution
%
%   eNorm = [eNorm10; eNorm11; eNorm1] where
%       eNorm10   - Energy associated with sol.
%                   = \sum_{n=1}^{n_{el}} \int_{x_n}^{x_{n=1}} u^2/2 dx
%       eNorm11   - Energy associated with derivative
%                   = \sum_{n=1}^{n_{el}} \int_{x_n}^{x_{n=1}} (dudx)^2/2 dx
%       eNorm1    - Energy norm
%                   = eNorm10 + c1 eNorm11
%
% efreq     - Frequency of energy norm evaluation

%% Global variables and parameters

global f1u f2u su a1u a2u alp10 alp11 alp20 alp21 alpu ...
    x u xbdry xedge h nelem xi P  ...
    xex nex plotStyle ...
    phi0L phi0R phi1L phi1R D gL0pv gR0pv gL1pv gR1pv

% Conservation law
f1u = @(u) u;                                                               % Convective flux
a1u = @(u) 1;                                                               % Analytical wavespeed for 1st order flux
nu = 1e-2;                                                                  % Diffusion coefficient
f2u = @(dudx) -nu*dudx;                                                     % Diffusive flux
a2u = @(dudx) -nu;                                                          % Analytical wavespeed for 2nd order flux
su = @(u) 0*u;                                                              % Source term
alp10 = alpArr(1,1);                                                        % Upwinding coefficient for 0th order corrections for 1st order flux
alp11 = alpArr(1,2);                                                        % Upwinding coefficient for 1st order corrections for 1st order flux
alpu = 1;                                                                   % Upwinding coefficient for u
alp20 = alpArr(2,1);                                                        % Upwinding coefficient for 0th order corrections for 2nd order flux
alp21 = alpArr(2,2);                                                        % Upwinding coefficient for 1st order corrections for 2nd order flux

% FR scheme
P = P_in;                                                                   % Polynomial order
s = P+3;                                                                    % Order of correction polynomial
c1 = c1byJsq * (h_in/2)^2;                                                  % c1 parameter

% Computational domain
xbdry = [-10, 10];                                                          % Domain
nelem = round((xbdry(2)-xbdry(1)) / h_in);                                  % No. of elements
xedge = linspace(xbdry(1),xbdry(2),nelem+1);                                % Element edges
xi = GLpoints(P);                                                           % Soln. points
CFL = DGCFL(P,1)* 0.20;                                                     % CFL no.

% Initial condition
hAv = (xbdry(2)-xbdry(1)) / nelem;
k0 = k0F* ((P+1)*pi/hAv);                                                   % Wavenumber
u0 = @(x) sin(k0 *x);                                                       % Initial condition
u0_d = @(x) k0*cos(k0 *x);                                                  % Derivative of initial condition

% Auxiliary
tf = -log(0.01)/(nu*k0^2);                                                  % Final time
nex = 1e2;                                                                  % No. of pts. for ex. soln. per element
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

% Pth order Lagrange function values and derivatives at edges in ref. domain
phi0L = zeros(P+1,1); phi0R = zeros(P+1,1);
phi1L = zeros(P+1,1); phi1R = zeros(P+1,1);
for p = 1: P+1
    phi0L(p) = Lagrange(xi,p,-1);
    phi0R(p) = Lagrange(xi,p,1);
    phi1L(p) = Lagrange_p(xi,p,-1);
    phi1R(p) = Lagrange_p(xi,p,1);
end

% Pth order Lagrange function derivatives at solution points
D = zeros(P+1,P+1);
for k = 1: P+1
    for p = 1: P+1
        D(k,p) = Lagrange_p(xi,p,xi(k));
    end
end

% Derivative of correction function
gL0pv = zeros(P+1,nelem); gR0pv = zeros(P+1,nelem);
gL1pv = zeros(P+1,nelem); gR1pv = zeros(P+1,nelem);
for j = 1: nelem
    theta = gLC1FR(P,h(j)/2,c1,0);
    for k = 1: s
        gL0pv(:,j) = gL0pv(:,j) + k*theta(k+1)*(xi.^(k-1));
        gR0pv(:,j) = gR0pv(:,j) + ((-1)^k)* k*theta(k+1)*(xi.^(k-1));
    end
    theta = gLC1FR(P,h(j)/2,c1,1);
    for k = 1: s
        gL1pv(:,j) = gL1pv(:,j) + k*theta(k+1)*(xi.^(k-1));
        gR1pv(:,j) = gR1pv(:,j) - ((-1)^k)* k*theta(k+1)*(xi.^(k-1));
    end
end

% Exact solution
uex = @(t) exp(-nu*k0^2*t)* u0((xbdry(1)+ mod((xex-t)-xbdry(1),xbdry(2)-xbdry(1)) ));
uex_d = @(t) exp(-nu*k0^2*t)* u0_d((xbdry(1)+ mod((xex-t)-xbdry(1),xbdry(2)-xbdry(1)) ));

% Style
plotStyle = {'k--','r-','g-','b-','m-','g-'};

%% FR update

% Setup run
u = u0(x); t = zeros(1,1e5);
cnt = 1; iter = 1;

% Setup norm arrays
eNormEx = zeros(3,1e5); eNormNum = zeros(3,1e5);
teN = zeros(1,1e5);

% Setup solution arrays
tPlot = [0,1,2,3,4]*(tf/4); xPlot = xex;
uPlot = zeros(nex,nelem,5); uEx = zeros(nex,nelem,5);
uCnt = 1;

% Update
while 1
    
    % Error norm
    if (mod(iter,efreq)==0)
        teN(cnt) = t(iter);
        [eNormEx(3,cnt),eNormEx(1,cnt),eNormEx(2,cnt)] = energyNormEx(uex(t(iter)), uex_d(t(iter)), c1);
        [eNormNum(3,cnt),eNormNum(1,cnt),eNormNum(2,cnt)] = energyNormNum(c1);
        cnt = cnt+1;
    end
    
    % Stage 1
    ui = u;
    [Fp0,amax] = CmFR_adv_diff(u);
    
    % Time step
    dt = CFL*(min(h)/amax);
    
    % Store solution
    logT = t(iter)<= tPlot & t(iter)+dt>=tPlot;
    if (max(logT)==1)
        uEx(:,:,uCnt) = uex(t(iter));
        for j = 1: nelem
            for k = 1: nex
                xik = -1 + (2/h(j))*(xex(k,j)-xedge(j));
                for p = 1: P+1
                    uPlot(k,j,uCnt) = uPlot(k,j,uCnt) + u(p,j)*Lagrange(xi,p,xik);
                end
            end
        end
        uCnt = uCnt + 1;
    end
    
    % Plots of solution
    if (mod(iter,plotFreq) == 0)
        plotFR2(1,'C1 continuous FR',t(iter),uex(t(iter)),1);
        hold off;
        pause(tDelay);
    end
    
    u = ui + (dt/2)*Fp0;
    
    % Stage 2
    Fp1 = CmFR_adv_diff(u);
    u = ui + (dt/2)*Fp1;
    
    % Stage 3
    Fp2 = CmFR_adv_diff(u);
    u = ui + dt*Fp2;
    
    % Stage 4
    Fp3 = CmFR_adv_diff(u);
    u = ui + (dt/6)*(Fp0 + 2*Fp1 + 2*Fp2 + Fp3);
    
    % Update
    iter = iter + 1;
    t(iter) = t(iter-1) + min(dt);
    
    % Check for completion
    if (t(iter)>tf); break; end
    
end

%% Truncate zero entries

teN = teN(1:cnt-1);
eNormEx = eNormEx(:,1:cnt-1);
eNormNum = eNormNum(:,1:cnt-1);
