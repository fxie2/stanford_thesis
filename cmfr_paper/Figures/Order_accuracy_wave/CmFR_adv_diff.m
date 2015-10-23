function [rhs,amax] = CmFR_adv_diff(u)

global f1u a1u f2u a2u alp10 alp11 alp20 alp21 alpu su h nelem P phi0L phi0R phi1L phi1R D gL0pv gR0pv gL1pv gR1pv
 
%% 1st order flux 
[f1Cp,a1max] = CmFR_O1(u,f1u,a1u,[alp10;alp11]);

%% Derivative of Cm continuous solution

% Edge values of discontinuous solution polynomial and derivative
uE = [phi0L';phi0R']*u;
upE = [phi1L';phi1R']*u;

% Discontinuous derivative of solution
Up = D*u;

% Upwind fluxes and wavespeeds
UInt0 = zeros(1,nelem+1); UInt1 = zeros(1,nelem+1);

% Specify left-boundary conditions
%uE(2,1) = 1;
%upE(2,1) = 0.5;
%UInt1(1) = 0.5;

for j = 1: nelem
    jm1 = j-1;
    if (j==1); jm1 = 1;end%jm1=nelem; end;
    
    a0 = 1;
    a1 = 1;
    
    UInt0(j) = 0.5*(uE(2,jm1)+uE(1,j)) - 0.5*abs(a0)*(1-alpu)*(uE(1,j)-uE(2,jm1));
    UInt1(j) = 0.5*(upE(2,jm1)+upE(1,j)) - 0.5*abs(a1)*(1-alpu)*(upE(1,j)-upE(2,jm1));
end
UInt0(nelem+1) = UInt0(1);
UInt1(nelem+1) = UInt1(1);

% Derivative of C1 continuous solution
UCp = zeros(P+1,nelem);
for j = 1: nelem
    UCp(:,j) = (2/h(j)) * ( Up(:,j) + ...
        gL0pv(:,j)*(UInt0(j)-uE(1,j))   + gR0pv(:,j)*(UInt0(j+1)-uE(2,j)) +...
        (2/h(j)) * (gL1pv(:,j)*(UInt1(j)-upE(1,j))  + gR1pv(:,j)*(UInt1(j+1)-upE(2,j))) );
end

%% 2nd order flux
[f2Cp,a2max] = CmFR_O1(UCp,f2u,a2u,[alp20;alp21]);

%% Update equation RHS term

% Source terms
sv = su(u);

% RHS term
rhs = zeros(P+1,nelem);
for j = 1: nelem
    rhs(:,j) = -(f1Cp(:,j) + f2Cp(:,j)) ...
               + sv(:,j);
end

% Maximum wavespeed
amax = max(a1max,a2max/min(h));

end