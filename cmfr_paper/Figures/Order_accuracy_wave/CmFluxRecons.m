% CmFR1D - Cm Flux Reconstruction in 1D
% ---
% Time marches the hyperbolic scalar conservation law dudt + df1(u)dx + df2(dudx)dx = s
% Time stepping: RK4
%
% - Kartikey Asthana
% ---
function [teN, errNorm01,normConst01] = CmFluxRecons(P_in, h_in, c1byJsq, alpArr, dt, tf, efreq, plotFreq)
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
% normConst01 - normalization constants for energy

%% Global variables and parameters

global f1u f2u su a1u a2u alp10 alp11 alp20 alp21 alpu nu ...
    x u xbdry xedge h nelem xi P  ...
    xex nex plotStyle ...
    phi0L phi0R phi1L phi1R D gL0pv gR0pv gL1pv gR1pv ...
    k0 kN

% Conservation law
f1u = @(u) u;                                                               % Convective flux
a1u = @(u) 1;                                                               % Analytical wavespeed for 1st order flux
nu = 0;1e-2;                                                                  % Diffusion coefficient
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

% Initial condition
u0 = @(x) zeros(size(x));%init_cond(x,0);                                                   % Initial condition

% Auxiliary
k0 = 2*pi/(xbdry(2)-xbdry(1));                                              % Base wavenumber for the domain
%kN = (P+1)*pi;                                                              % Nyquist frequency for P=P_in, h=1
k0 = k0 ;
kN = k0;
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
% uex = @(x,t) exact_sol(x, t, 0);
% uex_d = @(x,t) exact_sol(x, t, 1);
uexV = @(t) exact_sol(x, t, 0);
uex_dV = @(t) (2/h_in)*D*exact_sol(x, t, 0);
uex = @(x,t) eval_poly_u(uexV(t), x);
uex_d = @(x,t) eval_poly_u(uex_dV(t), x);

% Style
plotStyle = {'k-','r-','g-','b-','m-','g-'};

%% FR update

% Normalization constants
normConst01 = zeros(2,1);
u = zeros(size(x));
[normConst01(1), normConst01(2)] = errNorm(uex, uex_d, 0);

% Setup run
u = u0(x); t = zeros(1,1e5);
cnt = 1; iter = 1;

% Setup norm arrays
errNorm01 = zeros(2,1e5);
teN = zeros(1,1e5);

 figure(1)
% Update
while 1
    
    % Error norm
%     if (mod(iter,efreq)==0)
%         teN(cnt) = t(iter);
%         [errNorm01(1,cnt), errNorm01(2,cnt)] = errNorm(uex, uex_d, t(iter));
%         cnt = cnt+1;
%     end
    
    % Stage 1
    %u(1,1) = sin((t(iter))*(3));
    %u(2,1) = u(1,1);
    ui = u;
    Fp0 = CmFR_adv_diff(u);
         
    % Plots of solution
    if (mod(iter,plotFreq) == 0)
        plotFR2p(1,'C1 continuous FR',t(iter),uex(xex,t(iter)),1);
        hold off;
        pause(tDelay);
    end
    
    u = ui + (dt/2)*Fp0;
    
    % Stage 2
    u(1,1) = sin((t(iter)+dt/2)*(3));
    Fp1 = CmFR_adv_diff(u);
    u = ui + (dt/2)*Fp1;
    
    % Stage 3
    u(1,1) = sin((t(iter)+dt/2)*(3));
    Fp2 = CmFR_adv_diff(u);
    u = ui + dt*Fp2;
    
    % Stage 4
    u(1,1) = sin((t(iter)+dt)*(3));
    Fp3 = CmFR_adv_diff(u);
    u = ui + (dt/6)*(Fp0 + 2*Fp1 + 2*Fp2 + Fp3);
    
    % Update
    iter = iter + 1;
    t(iter) = t(iter-1) + min(dt);
    
    % Check for completion
    if (t(iter)>tf); break; end
    
     plot(reshape(x,[1 numel(x)]),reshape(u,[1 numel(u)]))
     drawnow
    
end

%% Truncate zero entries

cnt = 1;
[errNorm01(1,cnt), errNorm01(2,cnt)] = errNorm(uex, uex_d, t(iter-1));

teN = teN(1:cnt);
errNorm01 = errNorm01(:,1:cnt);
