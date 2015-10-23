function [eNormTot,eNorm0,eNorm1] = energyNormEx(uex, uex_d, c1)
% Routine for calculating components of energy norm
% ---
% uex   - Exact solution
% uex_d - Derivative of exact solution
% c1    - CmFR parameter
% -------------------------------------------------------------------------

global nelem xex

eNorm0 = 0;
eNorm1 = 0;

for j = 1: nelem
    
    % 1,0 component
    e0 = uex(:,j).*uex(:,j) / 2;
    eNorm0 = eNorm0 + trapz(xex(:,j),e0);
    
    % 1,1 component
    e1 = uex_d(:,j).*uex_d(:,j) / 2;
    eNorm1 = eNorm1 + trapz(xex(:,j),e1);
       
end

% Total 1 component
eNormTot = eNorm0 + c1*eNorm1;
