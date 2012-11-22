% Variation of the speckle contrast with the ratio correlation time to
% exposure time \tau_c/T

% Ratio correlation time to exposure time t = tau_c/T
t = logspace(-4, 4, 2^10);
% Speckle contrast C
C = sqrt((t/2).*(1-exp(-2./t)));

%% Plot speckle contrast
figure; set(gcf,'color','w')
semilogx(t, C, 'k-', 'LineWidth', 3)
% label fonts
labelFontSize   = 18;
axisFontSize    = 16;
% xlabel('$\displaystyle\frac{\tau_c}{T}$', 'interpreter', 'latex', 'FontSize', labelFontSize, 'FontWeight', 'bold')
xlabel('\tau_c/T', 'FontSize', labelFontSize, 'FontWeight', 'bold')
ylabel('C', 'FontSize', labelFontSize, 'FontWeight', 'bold')
% ylabel('$\displaystyle\frac{dW}{dN}*$','interpreter','latex')
set(gca, 'FontSize', axisFontSize)

%% Save figure
addpath(genpath('D:\Edgar\ssoct\Matlab'))
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures','speckle'),'-png',gcf)
