%% Lagrange functions
function phi_j = Lagrange(xGrid,j,x)
% ---
% Evaluates the j'th Lagrange function on xGrid at x
% xGrid - real array
% j - index for the function
% x - location of evaluation
% ---
n = numel(xGrid);
phi_j = 1;
for i = 1: n
    if(i~=j); phi_j = phi_j.*(x-xGrid(i))./(xGrid(j)-xGrid(i)); end
end
end