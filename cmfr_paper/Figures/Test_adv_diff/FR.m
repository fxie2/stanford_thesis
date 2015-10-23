function [rhs,amax] = FR(u)

global f1u a1u f2u a2u alp1 alp2 alpu su h nelem P phiL phiR D gpLv gpRv

%% 1st order flux 

% Discontinuous derivative of discontinuous 1st order flux
f1p = D*f1u(u);

% Edge values of discontinuous solution
uE = [phiL';phiR']*u;

% Upwind fluxes and wavespeeds
f1upw = zeros(1,nelem+1); a1max = eps;
for j = 1: nelem
    jm1 = j-1;
    if (j==1); jm1=nelem; end;
    
    if (uE(1,j)==uE(2,jm1))
        a = a1u(uE(1,j));
    else
        a = (f1u(uE(1,j))-f1u(uE(2,jm1)))/(uE(1,j)-uE(2,jm1));
    end
    a1max = max(abs(a),a1max);
    f1upw(j) = 0.5*(f1u(uE(2,jm1))+f1u(uE(1,j))) - 0.5*abs(a)*(1-alp1)*(uE(1,j)-uE(2,jm1));
end
f1upw(nelem+1) = f1upw(1);

% Edge values of discontinuous 1st order flux
f1E = [phiL';phiR']*f1u(u);

% Discontinuous derivative of continuous 1st order flux
f1Cp = zeros(P+1,nelem);
for j = 1: nelem
    f1Cp(:,j) = (2/h(j)) * (f1p(:,j) + ...
        gpLv*(f1upw(j)-f1E(1,j)) + ...
        gpRv*(f1upw(j+1)-f1E(2,j)));
end

%% Continuous dudx required for 2nd order flux

% Upwind fluxes and wavespeeds
Uupw = zeros(1,nelem+1);
for j = 1: nelem
    jm1 = j-1;
    if (j==1); jm1=nelem; end;
    
    a = 1;
    Uupw(j) = 0.5*(uE(2,jm1)+uE(1,j)) - 0.5*abs(a)*(1-alpu)*(uE(1,j)-uE(2,jm1));
end
Uupw(nelem+1) = Uupw(1);

% Discontinuous derivative of solution
Up = D*u;

% Discontinuous derivative of continuous solution
UCp = zeros(P+1,nelem);
for j = 1: nelem
    UCp(:,j) = (2/h(j)) * ( Up(:,j) + ...
        gpLv*(Uupw(j)-uE(1,j)) + ...
        gpRv*(Uupw(j+1)-uE(2,j)) );
end


%% 2nd order flux

% Discontinuous derivative of discontinuous 2nd order flux
f2p = D*f2u(UCp);

% Edge values of discontinuous derivative of continuous solution
UCpE = [phiL';phiR']*UCp;

% Upwind fluxes and wavespeeds
f2upw = zeros(1,nelem+1); a2max = eps;
for j = 1: nelem
    jm1 = j-1;
    if (j==1); jm1=nelem; end;
    
    if (UCpE(1,j)==UCpE(2,jm1))
        a = a2u(UCpE(1,j));
    else
        a = (f2u(UCpE(1,j))-f2u(UCpE(2,jm1)))/(UCpE(1,j)-UCpE(2,jm1));
    end
    a2max = max(abs(a),a2max);
    f2upw(j) = 0.5*(f2u(UCpE(2,jm1))+f2u(UCpE(1,j))) - 0.5*abs(a)*(1-alp2)*(UCpE(1,j)-UCpE(2,jm1));
end
f2upw(nelem+1) = f2upw(1);

% Edge values of discontinuous 2nd order flux
f2E = [phiL';phiR']*f2u(UCp);

% Discontinuous derivative of continuous 2nd order flux
f2Cp = zeros(P+1,nelem);
for j = 1: nelem
    f2Cp(:,j) = (2/h(j)) * ( f2p(:,j) + ...
        gpLv*(f2upw(j)-f2E(1,j)) + ...
        gpRv*(f2upw(j+1)-f2E(2,j)) );
end


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
amax = max(a1max,a2max);

end