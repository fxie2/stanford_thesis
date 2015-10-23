% Computes derivative of C1 continuous flux at solution points
function [Fp,amax] = CmFR_O1(u,fu,au,alp)
% ---
% Input arguments
% u   - discrete discontinuous solution (Pth order polynomial)
% fu  - function handle for the flux
% au  - function handle for the wavespeed
% alp - upwinding coefficient (in [0,1])
% ---
% Output arguments
% Fp   - discrete derivative of C1 continuous flux
% amax - max. wavespeed in the domain
% -----

global h nelem P phi0L phi0R phi1L phi1R D gL0pv gR0pv gL1pv gR1pv

% Edge values of discontinuous solution polynomial and derivative
uE = [phi0L';phi0R']*u;
upE = [phi1L';phi1R']*u;

% Edge values of discontinuous flux polynomial and derivative
fE = [phi0L';phi0R']*fu(u);
fpE = [phi1L';phi1R']*fu(u);

% Upwind fluxes and wavespeeds
fInt0 = zeros(1,nelem+1); fInt1 = zeros(1,nelem+1);
amax = eps;
for j = 1: nelem
    jm1 = j-1;
    if (j==1); jm1=nelem; end;
    
    a0 = au(uE(1,j));
    a1 = au(upE(1,j));
    amax = max(max(abs(a0),abs(a1)),amax);
    
    fInt0(j) = 0.5*(fu(uE(2,jm1)) + fu(uE(1,j)))   - 0.5*(1-alp(1))* abs(a0) * (uE(1,j)-uE(2,jm1));
    fInt1(j) = 0.5*(fu(upE(2,jm1)) + fu(upE(1,j))) - 0.5*(1-alp(2))* abs(a1) * (upE(1,j)-upE(2,jm1));
end
fInt0(nelem+1) = fInt0(1);
fInt1(nelem+1) = fInt1(1);

% Derivative of discontinuous flux
fp = D*fu(u);

% Derivative of C1 continuous flux
Fp = zeros(P+1,nelem);
for j = 1: nelem
    Fp(:,j) = (2/h(j)) * ( fp(:,j) + ...
        gL0pv(:,j)*(fInt0(j)-fE(1,j))   + gR0pv(:,j)*(fInt0(j+1)-fE(2,j)) +...
        (2/h(j)) * (gL1pv(:,j)*(fInt1(j)-fpE(1,j))  + gR1pv(:,j)*(fInt1(j+1)-fpE(2,j))) );
end

end