%% Derivative of Lagrange functions
function phip_j = Lagrange_p(xGrid,j,x)
% ---
% Evaluates the derivative of the j'th Lagrange function on xGrid at x
% xGrid - real array
% j - index for the function
% x - location of evaluation
% ---
n = numel(xGrid);%,1);
phip_j = zeros(size(x));
den = 1;
for i = 1: n
    if(i~=j); den = den*(xGrid(j)-xGrid(i)); end;
end
for k = 1: n
    if (k==j); continue; end;
    termk = 1;
    for i = 1: n
        if(i~=j && i~=k); termk = termk.*(x-xGrid(i)); end;
    end
    phip_j = phip_j + termk;
end
phip_j = phip_j / den;
end