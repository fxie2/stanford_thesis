function tAE = tAvgError(teN, errNorm01, normConst01)
% --- 
% time averaged error in the 0th and 1st derivative
% ---
tAE(1) = trapz(teN, errNorm01(1,:)); % 0th derivative
tAE(2) = trapz(teN, errNorm01(2,:)); % 1st derivative

tAE(1) = tAE(1) ./ ( (teN(end) - teN(1))*normConst01(1) );
tAE(2) = tAE(2) ./ ( (teN(end) - teN(1))*normConst01(2) );
