function tAE = tAvgError(teN, eNormEx, eNormNum)
% --- 
% time averaged error (absolute) in the energies of 0th and 1st derivative
% ---
err = abs(eNormEx - eNormNum) ./ (eNormNum(:,1)*ones(1,length(teN)));
tAE(1) = trapz(teN, err(1,:)) / (teN(end) - teN(1)); % 0th derivative
tAE(2) = trapz(teN, err(2,:)) / (teN(end) - teN(1)); % 1st derivative

