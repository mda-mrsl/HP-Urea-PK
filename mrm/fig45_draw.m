clear variables
close all

%%
outDir = mfilename;
if exist(outDir, 'dir')
    rmdir(outDir, 's')
end
mkdir(outDir)

cm = lines(5);

%%
load('fig45_snrPerf_Results-W10L-38VYN53-20230720_125633.mat', ...
    'resVb', 'resVei', 'resVe', 'resKve')
% rvk: Results for VIF Known
rvkVb  = resVb;
rvkVei = resVei;
rvkVe  = resVe;
rvkKve = resKve;
clearvars res*

rvkVbMn = mean(rvkVb, 3);
rvkVbSd = std(rvkVb, 0, 3);
rvkVbCv = rvkVbSd ./ rvkVbMn;

rvkVeiMn = mean(rvkVei, 3);
rvkVeiSd = std(rvkVei, 0, 3);
rvkVeiCv = rvkVeiSd ./ rvkVeiMn;

rvkVeMn = mean(rvkVe, 3);
rvkVeSd = std(rvkVe, 0, 3);
rvkVeCv = rvkVeSd ./ rvkVeMn;

rvkKveMn = mean(rvkKve, 3);
rvkKveSd = std(rvkKve, 0, 3);
rvkKveCv = rvkKveSd ./ rvkKveMn;

rvkKveVeMn = mean(rvkKve./rvkVe, 3);
rvkKveVeSd = std(rvkKve./rvkVe, 0, 3);
rvkKveVeCv = rvkKveVeSd ./ rvkKveVeMn;

rvkKveVbMn = mean(rvkKve./rvkVb, 3);
rvkKveVbSd = std(rvkKve./rvkVb, 0, 3);
rvkKveVbCv = rvkKveVbSd ./ rvkKveVbMn;

load('fig45_snrPerf_vscFit_Results-W10L-38VYN53-20230720_133454.mat')
% rvf: Results for VIF scale Fit
rvfVb  = resVb;
rvfVei = resVei;
rvfVe  = resVe;
rvfKve = resKve;
rvfVsc = resVsc;
clearvars res*

rvfVbMn = mean(rvfVb, 3);
rvfVbSd = std(rvfVb, 0, 3);
rvfVbCv = rvfVbSd ./ rvfVbMn;

rvfVeiMn = mean(rvfVei, 3);
rvfVeiSd = std(rvfVei, 0, 3);
rvfVeiCv = rvfVeiSd ./ rvfVeiMn;

rvfVeMn = mean(rvfVe, 3);
rvfVeSd = std(rvfVe, 0, 3);
rvfVeCv = rvfVeSd ./ rvfVeMn;

rvfKveMn = mean(rvfKve, 3);
rvfKveSd = std(rvfKve, 0, 3);
rvfKveCv = rvfKveSd ./ rvfKveMn;

rvfVscMn = mean(rvfVsc, 3);
rvfVscSd = std(rvfVsc, 0, 3);
rvfVscCv = rvfVscSd ./ rvfVscMn;

rvfKveVeMn = mean(rvfKve./rvfVe, 3);
rvfKveVeSd = std(rvfKve./rvfVe, 0, 3);
rvfKveVeCv = rvfKveVeSd ./ rvfKveVeMn;

rvfKveVbMn = mean(rvfKve./rvfVb, 3);
rvfKveVbSd = std(rvfKve./rvfVb, 0, 3);
rvfKveVbCv = rvfKveVbSd ./ rvfKveVbMn;


%% Fig4
fig4 = figure('Position', [10 10 1300 1300]);
tiledlayout(2, 2);

% UL - vb vs. SNR, VIF known
nexttile;
errorbar(tSnr-0.8, rvkVbMn(1,:), rvkVbSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr+0.8, rvkVbMn(3,:), rvkVbSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.27])
xlabel('Reference SNR')
ylabel({'v_b (v/v)', 'Accurate VIF'})
grid on
legend('Model I', 'Model III', 'True v_b', ...
    'location', 'northeast')
set(gca, 'fontsize', 16)

% UR - vb vs. SNR, VIF scale fit
nexttile;
errorbar(tSnr, rvfVbMn(1,:), rvfVbSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvfVbMn(3,:), rvfVbSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.27])
xlabel('Reference SNR')
ylabel({'v_b (v/v)', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model III', 'True v_b', ...
    'location', 'east')
set(gca, 'fontsize', 16)

% LL - ve vs. SNR, VIF known
nexttile;
errorbar(tSnr-0.8, rvkVeMn(1,:), rvkVeSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr+0.8, rvkVeMn(2,:), rvkVeSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr, rvkVeMn(3,:), rvkVeSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.92])
xlabel('Reference SNR')
ylabel({'v_e (v/v)', 'Accurate VIF'})
grid on
lObj = legend('Model I (v_e_e)', 'Model II (v_e_c)', 'Model III (v_e_v)', 'True v_e_e', ...
    'location', 'east');
% lObj = legend('Model I (v_e_e)', 'Model II (v_e_c)', 'Model III (v_e_v)', 'True v_e_e', ...
%     'location', 'best');
set(gca, 'fontsize', 16)

% LR - ve vs. SNR, VIF scale fit
nexttile;
errorbar(tSnr-0.8, rvfVeMn(1,:), rvfVeSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr+0.8, rvfVeMn(2,:), rvfVeSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr, rvfVeMn(3,:), rvfVeSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.92])
xlabel('Reference SNR')
ylabel({'v_e (v/v)', 'VIF Scale Fit'})
grid on
legend('Model I (v_e_e)', 'Model II (v_e_c)', 'Model III (v_e_v)', 'True v_e_e', ...
    'location', 'east')
set(gca, 'fontsize', 16)

% subfig labels
annotation('textbox', [0.0435 0.9253 0.0683 0.0538], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4967 0.9243 0.0683 0.0538], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.0438 0.4622 0.0683 0.0538], 'String', 'C', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4968 0.4623 0.0683 0.0538], 'String', 'D', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

lObj.Position = lObj.Position + [-0.03 0.05 0 0];
drawnow

export_fig(fullfile(outDir, 'MRM_fig4_vb_ve_snrPerf'), '-tif', '-a1', '-nocrop', fig4);
savefig(fig4, fullfile(outDir, 'MRM_fig4_vb_ve_snrPerf'));


%% Fig5
fig5 = figure('Position', [100 100 1300 1300]);
tiledlayout(2, 2);

% UL - kve vs. SNR, VIF known
nexttile;
errorbar(tSnr-0.8, rvkKveMn(1,:), rvkKveSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvkKveMn(2,:), rvkKveSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr+0.8, rvkKveMn(3,:), rvkKveSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.17])
xlabel('Reference SNR')
ylabel({'k_v_e (s^-^1)', 'Accurate VIF'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', ...
    'location', 'northeast')
set(gca, 'fontsize', 16)

% UR - kve vs. SNR, VIF scale fit
nexttile;
errorbar(tSnr-0.8, rvfKveMn(1,:), rvfKveSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvfKveMn(2,:), rvfKveSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr+0.8, rvfKveMn(3,:), rvfKveSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.17])
xlabel('Reference SNR')
ylabel({'k_v_e (s^-^1)', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', ...
    'location', 'best');
set(gca, 'fontsize', 16)

% LL - kve CV vs. SNR, VIF known & scale fit
nexttile;
semilogy(tSnr, rvkKveCv(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
semilogy(tSnr, rvkKveCv(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvkKveCv(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
semilogy(tSnr, rvfKveCv(1,:), 'o:', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
semilogy(tSnr, rvfKveCv(2,:), 's:', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvfKveCv(3,:), 'd:', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
axis([8 52 0 1])
xlabel('Reference SNR')
ylabel('k_v_e Coefficient of Variation')
grid on
lObj = legend('Model I, Accurate VIF', 'Model II, Accurate VIF', 'Model III, Accurate VIF', ...
    'Model I, VIF Scale Fit', 'Model II, VIF Scale Fit', 'Model III, VIF Scale Fit', ...
    'location', 'northeast');
set(gca, 'fontsize', 16)

% LR - kve/ve ratio vs. SNR, VIF scale fit
nexttile;
errorbar(tSnr-0.8, rvfKveVeMn(1,:), rvfKveVeSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvfKveVeMn(2,:), rvfKveVeSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr+0.8, rvfKveVeMn(3,:), rvfKveVeSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tKve./tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.34])
xlabel('Reference SNR')
ylabel({'k_v_e/v_e Ratio', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e/v_e_e', ...
    'location', 'east')
set(gca, 'fontsize', 16)

% subfig labels
annotation('textbox', [0.0481 0.9245 0.0683 0.0538], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.5021 0.9243 0.0683 0.0538], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.0492 0.4622 0.0683 0.0538], 'String', 'C', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.5022 0.4623 0.0683 0.0538], 'String', 'D', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

lObj.Position = lObj.Position + [0.003 0.003 0 0];
drawnow

export_fig(fullfile(outDir, 'MRM_fig5_kve_snrPerf'), '-tif', '-a1', '-nocrop', fig5);
savefig(fig5, fullfile(outDir, 'MRM_fig5_kve_snrPerf'));

%% Fig5v2
fig5v2 = figure('Position', [100 100 1300 1300]);
tiledlayout(2, 2);

% UL - kve vs. SNR, VIF known
nexttile;
errorbar(tSnr-0.8, rvkKveMn(1,:), rvkKveSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvkKveMn(2,:), rvkKveSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr+0.8, rvkKveMn(3,:), rvkKveSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.17])
xlabel('Reference SNR')
ylabel({'k_v_e (s^-^1)', 'Accurate VIF'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', ...
    'location', 'northeast')
set(gca, 'fontsize', 16)

% UR - kve vs. SNR, VIF scale fit
nexttile;
errorbar(tSnr-0.8, rvfKveMn(1,:), rvfKveSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvfKveMn(2,:), rvfKveSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr+0.8, rvfKveMn(3,:), rvfKveSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.17])
xlabel('Reference SNR')
ylabel({'k_v_e (s^-^1)', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', ...
    'location', 'best');
set(gca, 'fontsize', 16)

% LL - kve CV vs. SNR, VIF known & scale fit
nexttile;
semilogy(tSnr, rvkKveCv(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
semilogy(tSnr, rvkKveCv(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvkKveCv(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
semilogy(tSnr, rvfKveCv(1,:), 'o:', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
semilogy(tSnr, rvfKveCv(2,:), 's:', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvfKveCv(3,:), 'd:', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
axis([8 52 0 1])
xlabel('Reference SNR')
ylabel('k_v_e Coefficient of Variation')
grid on
lObj = legend('Model I, Accurate VIF', 'Model II, Accurate VIF', 'Model III, Accurate VIF', ...
    'Model I, VIF Scale Fit', 'Model II, VIF Scale Fit', 'Model III, VIF Scale Fit', ...
    'location', 'northeast');
set(gca, 'fontsize', 16)

% LR - kve/vb ratio vs. SNR, VIF scale fit
nexttile;
errorbar(tSnr-0.8, rvfKveVbMn(1,:), rvfKveVbSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
% errorbar(tSnr, rvfKveVeMn(2,:), rvfKveVeSd(2,:), 's-', ...
%     'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
errorbar(tSnr+0.8, rvfKveVbMn(3,:), rvfKveVbSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], tKve./tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
axis([8 52 0 0.5])
xlabel('Reference SNR')
ylabel({'k_v_e/v_b Ratio', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model III', 'True k_v_e/v_b', ...
    'location', 'northeast')
set(gca, 'fontsize', 16)

% subfig labels
annotation('textbox', [0.0481 0.9245 0.0683 0.0538], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.5021 0.9243 0.0683 0.0538], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.0492 0.4622 0.0683 0.0538], 'String', 'C', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.5022 0.4623 0.0683 0.0538], 'String', 'D', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

lObj.Position = lObj.Position + [0.003 0.003 0 0];
drawnow

export_fig(fullfile(outDir, 'MRM_fig5v2_kve_snrPerf'), '-tif', '-a1', '-nocrop', fig5v2);
savefig(fig5v2, fullfile(outDir, 'MRM_fig5v2_kve_snrPerf'));


%% FigS3
figs3 = figure('Position', [10 10 1900 650]);
tiledlayout(1, 3);

% L - vb CV vs. SNR, VIF known & scale fit
nexttile;
semilogy(tSnr, rvkVbCv(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
semilogy(tSnr, rvkVbCv(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvkVbCv(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
semilogy(tSnr, rvfVbCv(1,:), 'o:', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
semilogy(tSnr, rvkVbCv(2,:), 's:', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvfVbCv(3,:), 'd:', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
axis([8 52 0 0.8])
xlabel('Reference SNR')
ylabel('v_b Coefficient of Variation')
grid on
lObj = legend('Model I, Accurate VIF', 'Model II, Accurate VIF', 'Model III, Accurate VIF', ...
    'Model I, VIF Scale Fit', 'Model II, VIF Scale Fit', 'Model III, VIF Scale Fit', ...
    'location', 'northeast');
set(gca, 'fontsize', 16)

% C - ve CV vs. SNR, VIF known & scale fit
nexttile;
semilogy(tSnr, rvkVeCv(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
semilogy(tSnr, rvkVeCv(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvkVeCv(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
semilogy(tSnr, rvfVeCv(1,:), 'o:', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
semilogy(tSnr, rvfVeCv(2,:), 's:', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
semilogy(tSnr, rvfVeCv(3,:), 'd:', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
axis([8 52 1e-3 1])
xlabel('Reference SNR')
ylabel('v_e Coefficient of Variation')
grid on
% lObj = legend('Model I, Accurate VIF', 'Model II, Accurate VIF', 'Model III, Accurate VIF', ...
%     'Model I, VIF Scale Fit', 'Model II, VIF Scale Fit', 'Model III, VIF Scale Fit', ...
%     'location', 'southeast');
set(gca, 'fontsize', 16)

% R - VIF scale vs. SNR, VIF scale fit
nexttile;
minVal = 0.01;
negErr = [rvfVscMn(1,1) - minVal, rvfVscSd(1,2:end)];
errorbar(tSnr-0.8, rvfVscMn(1,:), negErr, rvfVscSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:), 'MarkerSize', 10);
hold on
errorbar(tSnr, rvfVscMn(2,:), rvfVscSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 12);
negErr = [rvfVscMn(3,1) - minVal, rvfVscSd(3,2:end)];
errorbar(tSnr+0.8, rvfVscMn(3,:), negErr, rvfVscSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:), 'MarkerSize', 10);
line([8 52], [1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('Reference SNR')
ylabel('VIF Scale Factor')
grid on
set(gca, 'fontsize', 16, 'YScale', 'log')
axis([8 52 0.2 30])
legend('Model I', 'Model II', 'Model III', 'True k_v_e/v_e_e', ...
    'location', 'northeast')

% subfig labels
annotation('textbox', [0.0252 0.9196 0.0467 0.1076], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.3367 0.9174 0.0467 0.1076], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.6476 0.9183 0.0467 0.1076], 'String', 'C', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

lObj.Position = lObj.Position + [0.01 0.01 0 0];
drawnow

export_fig(fullfile(outDir, 'MRM_figs3_vb_ve_cv_vscAcc'), '-tif', '-a1', '-nocrop', figs3);
savefig(figs3, fullfile(outDir, 'MRM_figs3_vb_ve_cv_vscAcc'))

close all
return


%% F-tests
% ftp: F-test p-value

% kve, accurate vif, models I vs II/III
[~,ftpKveRvkMd2] = vartest2(squeeze(rvkKve(1,:,:)).', squeeze(rvkKve(2,:,:)).');
[~,ftpKveRvkMd3] = vartest2(squeeze(rvkKve(1,:,:)).', squeeze(rvkKve(3,:,:)).');

figf1 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
plot(tSnr, rvkKveSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
plot(tSnr, rvkKveSd(2,:).^2, 's-', 'Color', cm(1,:));
plot(tSnr, rvkKveSd(3,:).^2, 'd-', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
legend('Model I', 'Model II', 'Model III', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpKveRvkMd2, 's-', 'Color', cm(1,:));
hold on
semilogy(tSnr, ftpKveRvkMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model II', 'Model III', 'p = 0.01', 'location', 'east')
grid on
set(gca, 'fontsize', 16)
suptitle('k_v_e: F-test vs. Model I, Accurate VIF');
export_fig(fullfile(outDir, 'figf1_kve-pval_accVsc_ModComp'), '-tif', '-nocrop', figf1);
savefig(figf1, fullfile(outDir, 'figf1_kve-pval_accVsc_ModComp'));


% vb, accurate vif, models I vs II/III
[~,ftpVbRvkMd3] = vartest2(squeeze(rvkVb(1,:,:)).', squeeze(rvkVb(3,:,:)).');

figf2 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
plot(tSnr, rvkVbSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
plot(tSnr, rvkVbSd(2,:).^2, 's-', 'Color', cm(1,:));
plot(tSnr, rvkVbSd(3,:).^2, 'd-', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
legend('Model I', 'Model II', 'Model III', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpVbRvkMd3, 'd-', 'Color', cm(4,:));
hold on
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('v_b: F-test vs. Model I, Accurate VIF');
export_fig(fullfile(outDir, 'figf2_vb-pval_accVsc_ModComp'), '-tif', '-nocrop', figf2);
savefig(figf2, fullfile(outDir, 'figf2_vb-pval_accVsc_ModComp'));

% ve, accurate vif, models I vs II/III
[~,ftpVeRvkMd2] = vartest2(squeeze(rvkVe(1,:,:)).', squeeze(rvkVe(2,:,:)).');
[~,ftpVeRvkMd3] = vartest2(squeeze(rvkVe(1,:,:)).', squeeze(rvkVe(3,:,:)).');

figf3 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
plot(tSnr, rvkVeSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
plot(tSnr, rvkVeSd(2,:).^2, 's-', 'Color', cm(1,:));
plot(tSnr, rvkVeSd(3,:).^2, 'd-', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
legend('Model I', 'Model II', 'Model III', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpVeRvkMd2, 's-', 'Color', cm(1,:));
hold on
semilogy(tSnr, ftpVeRvkMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model II', 'Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('v_e: F-test vs. Model I, Accurate VIF')
export_fig(fullfile(outDir, 'figf3_ve-pval_accVsc_ModComp'), '-tif', '-nocrop', figf3);
savefig(figf3, fullfile(outDir, 'figf3_ve-pval_accVsc_ModComp'));


% kve, fit vif scale, models I vs II/III
[~,ftpKveRvfMd2] = vartest2(squeeze(rvfKve(1,:,:)).', squeeze(rvfKve(2,:,:)).');
[~,ftpKveRvfMd3] = vartest2(squeeze(rvfKve(1,:,:)).', squeeze(rvfKve(3,:,:)).');

figf4 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
plot(tSnr, rvfKveSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
plot(tSnr, rvfKveSd(2,:).^2, 's-', 'Color', cm(1,:));
plot(tSnr, rvfKveSd(3,:).^2, 'd-', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
legend('Model I', 'Model II', 'Model III', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpKveRvfMd2, 's-', 'Color', cm(1,:));
hold on
semilogy(tSnr, ftpKveRvfMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model II', 'Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('k_v_e: F-test vs. Model I, VIF Scale Fit')
export_fig(fullfile(outDir, 'figf4_kve-pval_fitVsc_ModComp'), '-tif', '-nocrop', figf4);
savefig(figf4, fullfile(outDir, 'figf4_kve-pval_fitVsc_ModComp'));


% vb, fit vif scale, models I vs II/III
[~,ftpVbRvfMd3] = vartest2(squeeze(rvfVb(1,:,:)).', squeeze(rvfVb(3,:,:)).');

figf5 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
plot(tSnr, rvfVbSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
plot(tSnr, rvfVbSd(2,:).^2, 's-', 'Color', cm(1,:));
plot(tSnr, rvfVbSd(3,:).^2, 'd-', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
legend('Model I', 'Model II', 'Model III', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpVbRvfMd3, 'd-', 'Color', cm(4,:));
hold on
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('v_b: F-test vs. Model I, VIF Scale Fit')
export_fig(fullfile(outDir, 'figf5_vb-pval_fitVsc_ModComp'), '-tif', '-nocrop', figf5);
savefig(figf5, fullfile(outDir, 'figf5_vb-pval_fitVsc_ModComp'));


% ve, fit vif scale, models I vs II/III
[~,ftpVeRvfMd2] = vartest2(squeeze(rvfVe(1,:,:)).', squeeze(rvfVe(2,:,:)).');
[~,ftpVeRvfMd3] = vartest2(squeeze(rvfVe(1,:,:)).', squeeze(rvfVe(3,:,:)).');

figf6 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
plot(tSnr, rvfVeSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
plot(tSnr, rvfVeSd(2,:).^2, 's-', 'Color', cm(1,:));
plot(tSnr, rvfVeSd(3,:).^2, 'd-', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
legend('Model I', 'Model II', 'Model III', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpVeRvfMd2, 's-', 'Color', cm(1,:));
hold on
semilogy(tSnr, ftpVeRvfMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model II', 'Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('v_e: F-test vs. Model I, VIF Scale Fit')
export_fig(fullfile(outDir, 'figf6_ve-pval_fitVsc_ModComp'), '-tif', '-nocrop', figf6);
savefig(figf6, fullfile(outDir, 'figf6_ve-pval_fitVsc_ModComp'));


%kve, accurate vs vif scale fit
[~,ftpKveRvcMd1] = vartest2(squeeze(rvkKve(1,:,:)).', squeeze(rvfKve(1,:,:)).');
[~,ftpKveRvcMd2] = vartest2(squeeze(rvkKve(2,:,:)).', squeeze(rvfKve(2,:,:)).');
[~,ftpKveRvcMd3] = vartest2(squeeze(rvkKve(3,:,:)).', squeeze(rvfKve(3,:,:)).');

figf7 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
semilogy(tSnr, rvkKveSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
semilogy(tSnr, rvfKveSd(1,:).^2, 'o:', 'Color', cm(5,:));
semilogy(tSnr, rvkKveSd(2,:).^2, 's-', 'Color', cm(1,:));
semilogy(tSnr, rvfKveSd(2,:).^2, 's:', 'Color', cm(1,:));
semilogy(tSnr, rvkKveSd(3,:).^2, 'd-', 'Color', cm(4,:));
semilogy(tSnr, rvfKveSd(3,:).^2, 'd:', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpKveRvcMd1, 'o-', 'Color', cm(5,:));
hold on
semilogy(tSnr, ftpKveRvcMd2, 's-', 'Color', cm(1,:));
semilogy(tSnr, ftpKveRvcMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model I', 'Model II', 'Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('k_v_e: F-test, VIF Scale Accurate vs Fit')
export_fig(fullfile(outDir, 'figf7_kve-pval_VscFitComp'), '-tif', '-nocrop', figf7);
savefig(figf7, fullfile(outDir, 'figf7_kve-pval_VscFitComp'));


%vb, accurate vs vif scale fit
[~,ftpVbRvcMd1] = vartest2(squeeze(rvkVb(1,:,:)).', squeeze(rvfVb(1,:,:)).');
[~,ftpVbRvcMd3] = vartest2(squeeze(rvkVb(3,:,:)).', squeeze(rvfVb(3,:,:)).');

figf8 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
semilogy(tSnr, rvkVbSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
semilogy(tSnr, rvfVbSd(1,:).^2, 'o:', 'Color', cm(5,:));
semilogy(tSnr, rvkVbSd(3,:).^2, 'd-', 'Color', cm(4,:));
semilogy(tSnr, rvfVbSd(3,:).^2, 'd:', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpVbRvcMd1, 'o-', 'Color', cm(5,:));
hold on
semilogy(tSnr, ftpVbRvcMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model I', 'Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
title('v_b: F-test, VIF Scale Accurate vs Fit')
export_fig(fullfile(outDir, 'figf8_vb-pval_VscFitComp'), '-tif', '-nocrop', figf8);
savefig(figf8, fullfile(outDir, 'figf8_vb-pval_VscFitComp'));


%ve, accurate vs vif scale fit
[~,ftpVeRvcMd1] = vartest2(squeeze(rvkVe(1,:,:)).', squeeze(rvfVe(1,:,:)).');
[~,ftpVeRvcMd2] = vartest2(squeeze(rvkVe(2,:,:)).', squeeze(rvfVe(2,:,:)).');
[~,ftpVeRvcMd3] = vartest2(squeeze(rvkVe(3,:,:)).', squeeze(rvfVe(3,:,:)).');

figf9 = figure('Position', [10 10 1100 420]);
tiledlayout(1, 2);
nexttile;
semilogy(tSnr, rvkVeSd(1,:).^2, 'o-', 'Color', cm(5,:));
hold on
semilogy(tSnr, rvfVeSd(1,:).^2, 'o:', 'Color', cm(5,:));
semilogy(tSnr, rvkVeSd(2,:).^2, 's-', 'Color', cm(1,:));
semilogy(tSnr, rvfVeSd(2,:).^2, 's:', 'Color', cm(1,:));
semilogy(tSnr, rvkVeSd(3,:).^2, 'd-', 'Color', cm(4,:));
semilogy(tSnr, rvfVeSd(3,:).^2, 'd:', 'Color', cm(4,:));
xlabel('Reference SNR')
ylabel('Variance')
grid on
set(gca, 'fontsize', 16)
nexttile;
semilogy(tSnr, ftpVeRvcMd1, 'o-', 'Color', cm(5,:));
hold on
semilogy(tSnr, ftpVeRvcMd2, 's-', 'Color', cm(1,:));
semilogy(tSnr, ftpVeRvcMd3, 'd-', 'Color', cm(4,:));
line([8 52], 0.01*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlim([8 52]);
xlabel('Reference SNR')
ylabel('p-value')
legend('Model I', 'Model II', 'Model III', 'p = 0.01', 'location', 'best')
grid on
set(gca, 'fontsize', 16)
suptitle('v_e: F-test, VIF Scale Accurate vs Fit')
export_fig(fullfile(outDir, 'figf9_ve-pval_VscFitComp'), '-tif', '-nocrop', figf9);
savefig(figf9, fullfile(outDir, 'figf9_ve-pval_VscFitComp'));



%% vb,ve,kve scatter figs
figx1 = figure;
plot(tSnr, squeeze(rvkVb(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvkVb(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0.02 0.36])
xlabel('Reference SNR')
ylabel({'v_b (v/v)', 'Accurate VIF Scale'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx1_vbScatter_accVsc'), '-tif', '-nocrop', figx1);
savefig(figx1, fullfile(outDir, 'figx1_vbScatter_accVsc'));

figx2 = figure;
plot(tSnr, squeeze(rvfVb(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvfVb(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0.02 0.36])
xlabel('Reference SNR')
ylabel({'v_b (v/v)', 'VIF Scale Fit'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx2_vbScatter_fitVsc'), '-tif', '-nocrop', figx2);
savefig(figx2, fullfile(outDir, 'figx2_vbScatter_fitVsc'));

figx3 = figure;
plot(tSnr, squeeze(rvkVe(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvkVe(2,:,:)), 's', 'Color', cm(1,:), 'MarkerSize', 7);
plot(tSnr, squeeze(rvkVe(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 1])
xlabel('Reference SNR')
ylabel({'v_e (v/v)', 'Accurate VIF Scale'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx3_veScatter_accVsc'), '-tif', '-nocrop', figx3);
savefig(figx3, fullfile(outDir, 'figx3_veScatter_accVsc'));

figx4 = figure;
plot(tSnr, squeeze(rvfVe(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvfVe(2,:,:)), 's', 'Color', cm(1,:), 'MarkerSize', 7);
plot(tSnr, squeeze(rvfVe(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 1])
xlabel('Reference SNR')
ylabel({'v_e (v/v)', 'VIF Scale Fit'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx4_veScatter_fitVsc'), '-tif', '-nocrop', figx4);
savefig(figx4, fullfile(outDir, 'figx4_veScatter_fitVsc'));

figx5 = figure;
plot(tSnr, squeeze(rvkKve(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvkKve(2,:,:)), 's', 'Color', cm(1,:), 'MarkerSize', 7);
plot(tSnr, squeeze(rvkKve(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.23])
xlabel('Reference SNR')
ylabel({'k_v_e (s^-^1)', 'Accurate VIF Scale'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx5_kveScatter_accVsc'), '-tif', '-nocrop', figx5);
savefig(figx5, fullfile(outDir, 'figx5_kveScatter_accVsc'));

figx6 = figure;
plot(tSnr, squeeze(rvfKve(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvfKve(2,:,:)), 's', 'Color', cm(1,:), 'MarkerSize', 7);
plot(tSnr, squeeze(rvfKve(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.23])
xlabel('Reference SNR')
ylabel({'k_v_e (s^-^1)', 'VIF Scale Fit'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx6_kveScatter_fitVsc'), '-tif', '-nocrop', figx6);
savefig(figx6, fullfile(outDir, 'figx6_kveScatter_fitVsc'));

%% VIF scale figs
figx7 = figure;
errorbar(tSnr, rvfVscMn(1,:), rvfVscSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:));
hold on
errorbar(tSnr, rvfVscMn(2,:), rvfVscSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
errorbar(tSnr, rvfVscMn(3,:), rvfVscSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:));
line([8 52], [1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 3.5])
xlabel('Reference SNR')
ylabel('VIF Scale Factor')
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx7_vscAcc_fitVsc'), '-tif', '-nocrop', figx7);
savefig(figx7, fullfile(outDir, 'figx7_vscAcc_fitVsc'));

figx8 = figure;
plot(tSnr, squeeze(rvfVsc(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvfVsc(2,:,:)), 's', 'Color', cm(1,:), 'MarkerSize', 7);
plot(tSnr, squeeze(rvfVsc(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], [1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 3.5])
xlabel('Reference SNR')
ylabel('VIF Scale Factor')
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figx8_vscScatter_fitVsc'), '-tif', '-nocrop', figx8);
savefig(figx8, fullfile(outDir, 'figx8_vscScatter_fitVsc'));

%% kve/ve figs
figy1 = figure;
errorbar(tSnr-0.5, rvkKveVeMn(1,:), rvkKveVeSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:));
hold on
errorbar(tSnr, rvkKveVeMn(2,:), rvkKveVeSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
errorbar(tSnr+0.5, rvkKveVeMn(3,:), rvkKveVeSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:));
line([8 52], tKve./tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.38])
xlabel('Reference SNR')
ylabel({'k_v_e/v_e Ratio', 'Accurate VIF Scale'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e/v_e_e', ...
    'location', 'east')
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy1_kveveratioAcc_accVsc'), '-tif', '-nocrop', figy1);
savefig(figy1, fullfile(outDir, 'figy1_kveveratioAcc_accVsc'));

figy2 = figure;
plot(tSnr-0.5, squeeze(rvkKve(1,:,:)./rvkVe(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvkKve(2,:,:)./rvkVe(2,:,:)), 's', 'Color', cm(1,:));
plot(tSnr+0.5, squeeze(rvkKve(3,:,:)./rvkVe(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tKve./tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.5])
xlabel('Reference SNR')
ylabel({'k_v_e/v_e Ratio', 'Accurate VIF Scale'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy2_kveveratioScatter_accVsc'), '-tif', '-nocrop', figy2);
savefig(figy2, fullfile(outDir, 'figy2_kveveratioScatter_accVsc'));

figy3 = figure;
errorbar(tSnr-0.5, rvfKveVeMn(1,:), rvfKveVeSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:));
hold on
errorbar(tSnr, rvfKveVeMn(2,:), rvfKveVeSd(2,:), 's-', ...
    'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
errorbar(tSnr+0.5, rvfKveVeMn(3,:), rvfKveVeSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:));
line([8 52], tKve./tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 50 0 0.38])
xlabel('Reference SNR')
ylabel({'k_v_e/v_e Ratio', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e/v_e_e', ...
    'location', 'east')
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy3_kveveratioAcc_fitVsc'), '-tif', '-nocrop', figy3);
savefig(figy3, fullfile(outDir, 'figy3_kveveratioAcc_fitVsc'));

figy4 = figure;
plot(tSnr-0.5, squeeze(rvfKve(1,:,:)./rvfVe(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr, squeeze(rvfKve(2,:,:)./rvfVe(2,:,:)), 's', 'Color', cm(1,:));
plot(tSnr+0.5, squeeze(rvfKve(3,:,:)./rvfVe(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tKve./tVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.5])
xlabel('Reference SNR')
ylabel({'k_v_e/v_e Ratio', 'VIF Scale Fit'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy4_kveveratioScatter_fitVsc'), '-tif', '-nocrop', figy4);
savefig(figy4, fullfile(outDir, 'figy4_kveveratioScatter_fitVsc'));

%% kve/vb figs
figy5 = figure;
errorbar(tSnr-0.5, rvkKveVbMn(1,:), rvkKveVbSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:));
hold on
errorbar(tSnr+0.5, rvkKveVbMn(3,:), rvkKveVbSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:));
line([8 52], tKve./tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.5])
xlabel('Reference SNR')
ylabel({'k_v_e/v_b Ratio', 'Accurate VIF Scale'})
grid on
legend('Model I', 'Model III', 'True k_v_e/v_b', ...
    'location', 'northeast')
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy5_kvevbratioAcc_accVsc'), '-tif', '-nocrop', figy5);
savefig(figy5, fullfile(outDir, 'figy5_kvevbratioAcc_accVsc'));

figy6 = figure;
plot(tSnr-0.5, squeeze(rvkKve(1,:,:)./rvkVb(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr+0.5, squeeze(rvkKve(3,:,:)./rvkVb(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tKve./tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 1.2])
xlabel('Reference SNR')
ylabel({'k_v_e/v_b Ratio', 'Accurate VIF Scale'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy6_kvevbratioScatter_accVsc'), '-tif', '-nocrop', figy6);
savefig(figy6, fullfile(outDir, 'figy6_kvevbratioScatter_accVsc'));

figy7 = figure;
errorbar(tSnr-0.5, rvfKveVbMn(1,:), rvfKveVbSd(1,:), 'o-', ...
    'LineWidth', 1.3, 'Color', cm(5,:));
hold on
errorbar(tSnr+0.5, rvfKveVbMn(3,:), rvfKveVbSd(3,:), 'd-', ...
    'LineWidth', 1.3, 'Color', cm(4,:));
line([8 52], tKve./tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 0.5])
xlabel('Reference SNR')
ylabel({'k_v_e/v_b Ratio', 'VIF Scale Fit'})
grid on
legend('Model I', 'Model III', 'True k_v_e/v_b', ...
    'location', 'northeast')
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy7_kvevbratioAcc_fitVsc'), '-tif', '-nocrop', figy7);
savefig(figy7, fullfile(outDir, 'figy7_kvevbratioAcc_fitVsc'));

figy8 = figure;
plot(tSnr-0.5, squeeze(rvfKve(1,:,:)./rvfVb(1,:,:)), 'o', 'Color', cm(5,:));
hold on
plot(tSnr+0.5, squeeze(rvfKve(3,:,:)./rvfVb(3,:,:)), 'd', 'Color', cm(4,:));
line([8 52], tKve./tVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
% axis([8 52 0 1.2])
xlabel('Reference SNR')
ylabel({'k_v_e/v_b Ratio', 'VIF Scale Fit'})
grid on
set(gca, 'fontsize', 16)
export_fig(fullfile(outDir, 'figy8_kvevbratioScatter_fitVsc'), '-tif', '-nocrop', figy8);
savefig(figy8, fullfile(outDir, 'figy8_kvevbratioScatter_fitVsc'));



close all
return
