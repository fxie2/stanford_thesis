function [eNormTot,eNorm0,eNorm1] = energyNormNum(c1)
% Routine for calculating components of energy norm
% ---
% c1    - CmFR parameter
% -------------------------------------------------------------------------

global u nelem h xedge xex nex P xi

%% Solution at xex grid

uSamp0 = zeros(nex,nelem);  % numerical soln.
uSamp1 = zeros(nex,nelem);  % derivative of numerical soln.

for j = 1: nelem
    
    % Solution value   
    for k = 1: nex
        xik = -1 + (2/h(j))*(xex(k,j)-xedge(j));
        for p = 1: P+1
            uSamp0(k,j) = uSamp0(k,j) + u(p,j)*Lagrange(xi,p,xik);
        end
    end
    
    % Solution derivative
    for k = 1: nex
        xik = -1 + (2/h(j))*(xex(k,j)-xedge(j));
        for p = 1: P+1
            uSamp1(k,j) = uSamp1(k,j) + u(p,j)*Lagrange_p(xi,p,xik);
        end
    end
    uSamp1(:,j) = uSamp1(:,j) / (h(j)/2);
    
end

%% Energy terms

eNorm0 = 0;
eNorm1 = 0;

for j = 1: nelem
    
    % 1,0 component
    e0 = uSamp0(:,j).*uSamp0(:,j) / 2;
    eNorm0 = eNorm0 + trapz(xex(:,j),e0);
    
    % 1,1 component
    e1 = uSamp1(:,j).*uSamp1(:,j) / 2;
    eNorm1 = eNorm1 + trapz(xex(:,j),e1);
       
end

% Total 1 component
eNormTot = eNorm0 + c1*eNorm1;
