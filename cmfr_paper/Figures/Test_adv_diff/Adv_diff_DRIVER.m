% Numerical test cases for C1FR
% ---
% - Kartikey Asthana
% ---

%% Initialization
path(path,'../Polynomials')
path(path,'../IO')
clear all;
clc;
saveFlag = 1;

for cvals = [-12:1:3]*1e-3

%% Parameters

% Initial sinusoid: wavenumber fraction of the Nyquist limit
k0F_low = 0.25;              % Low wavenumber fraction
k0F_med = 0.50;              % Medium wavenumber fraction
k0F_high = 0.75;             % High wavenumber fraction

% Upwinding coefficients - purely upwinded
alpArrC1FR = [0,0;0,0];     % coefficients for C1FR
alpArrDG = [0;0];           % coefficients for DG

% Polynomial order
P = 5;

% Normalized c1 value - reasonable guess close to the optimal
c1byJsq = cvals;%-12e-3;            % Optimal for 0th derivative = -0.0046

% Grid spacing
h = 1;                      % Unit spacing

% Frequencies
efreq = 1;                  % Frequency of norm evaluation
plotFreq = 0;             % Frequency of plotting

%% Test - 1: Low wavenumber case
%{
% C1FR case
[teN_C1FR, eNormNum_C1FR, eNormEx_C1FR, xPlot, tPlot, uPlot_C1FR, uEx] = CmFluxRecons(P, h, c1byJsq, alpArrC1FR, k0F_low, efreq, plotFreq);
postProc(teN_C1FR, eNormNum_C1FR, eNormEx_C1FR, xPlot, uPlot_C1FR, uEx, 'C1FR');

% DG baseline
[teN_DG, eNormNum_DG, eNormEx_DG, ~, ~, uPlot_DG, ~] = FluxRecons(1, P, h, c1byJsq, alpArrDG, k0F_low, efreq, 0);
postProc(teN_DG, eNormNum_DG, eNormEx_DG, xPlot, uPlot_DG, uEx, 'DG');

saveFlag = 1;
if saveFlag; print('-loose', '-depsc', ['low_k_c_'...
        num2str(ceil(1e4*c1byJsq)) '.eps']);end
%%
% Save figures and data
f1 = figure(11); saveas(f1,['sol_low_k_c_' num2str(ceil(1e4*c1byJsq))],'fig');
f20 = figure(20); saveas(f20,['eN0_low_k_c_' num2str(ceil(1e4*c1byJsq))],'fig');
f21 = figure(21); saveas(f21,['eN1_low_k_c_' num2str(ceil(1e4*c1byJsq))],'fig');
if saveFlag
save( ['data_low_k_c_'...
        num2str(ceil(1e4*c1byJsq)) '.mat'])
end

% Close figures
close all;
%}
%% Test - 2: Med wavenumber case

% C1FR case
[teN_C1FR, eNormNum_C1FR, eNormEx_C1FR, xPlot, tPlot, uPlot_C1FR, uEx] = CmFluxRecons(P, h, c1byJsq, alpArrC1FR, k0F_med, efreq, 0);
postProc(teN_C1FR, eNormNum_C1FR, eNormEx_C1FR, xPlot, uPlot_C1FR, uEx, 'C1FR');

% DG baseline
[teN_DG, eNormNum_DG, eNormEx_DG, ~, ~, uPlot_DG, ~] = FluxRecons(1, P, h, c1byJsq, alpArrDG, k0F_med, efreq, 0);
postProc(teN_DG, eNormNum_DG, eNormEx_DG, xPlot, uPlot_DG, uEx, 'DG');

saveFlag = 1;
if saveFlag; print('-loose', '-depsc', ['med_k_c_'...
        num2str(ceil(1e4*c1byJsq)) '.eps']);end
%%
% Save figures and data
f1 = figure(11); 
if saveFlag
saveas(f1,['sol_med_k_c_' num2str(ceil(1e4*c1byJsq))],'fig'); end
f20 = figure(20); if saveFlag 
    saveas(f20,['eN0_med_k_c_' num2str(ceil(1e4*c1byJsq))],'fig'); end
f21 = figure(21); if saveFlag 
    saveas(f21,['eN1_med_k_c_' num2str(ceil(1e4*c1byJsq))],'fig'); end
if saveFlag
save( ['data_med_k_c_'...
        num2str(ceil(1e4*c1byJsq)) '.mat'])
end

% Close figures
close all;

%% Test - 3: High wavenumber case
close all
clc
% C1FR case
[teN_C1FR, eNormNum_C1FR, eNormEx_C1FR, xPlot, tPlot, uPlot_C1FR, uEx] = CmFluxRecons(P, h, c1byJsq, alpArrC1FR, k0F_high, efreq, 0);
postProc(teN_C1FR, eNormNum_C1FR, eNormEx_C1FR, xPlot, uPlot_C1FR, uEx, 'C1FR');

% DG baseline
[teN_DG, eNormNum_DG, eNormEx_DG, ~, ~, uPlot_DG, ~] = FluxRecons(1, P, h, c1byJsq, alpArrDG, k0F_high, efreq, 0);
postProc(teN_DG, eNormNum_DG, eNormEx_DG, xPlot, uPlot_DG, uEx, 'DG');

saveFlag = 1;
if saveFlag; print('-loose', '-depsc', ['high_k_c_'...
        num2str(ceil(1e4*c1byJsq)) '.eps']);end

%%
% Save figures and data
f1 = figure(11); if saveFlag 
    saveas(f1,['sol_high_k_c_' num2str(ceil(1e4*c1byJsq))],'fig'); end
f20 = figure(20); if saveFlag
    saveas(f20,['eN0_high_k_c_' num2str(ceil(1e4*c1byJsq))],'fig'); end
f21 = figure(21); if saveFlag
    saveas(f21,['eN1_high_k_c_' num2str(ceil(1e4*c1byJsq))],'fig'); end
if saveFlag
save( ['data_high_k_c_'...
        num2str(ceil(1e4*c1byJsq)) '.mat'])
end
% Close figures
close all;

end