% Evaluate the piecewise polynomial from the discrete soln
function us = eval_poly_u(u, xp)
% ---
% u - discrete soln.
% xp - query array
% ---

global x P xedge

us = zeros(size(xp));

for k = 1: length(xp)
    
    % Locate query point element
    j = find(xedge >= xp(k), 1, 'first') - 1;
    
    % Evaluate polynomial
    for i = 1: P+1
        us(k) = us(k) + u(i,j)*Lagrange(x(:,j),i,xp(k));
    end

end