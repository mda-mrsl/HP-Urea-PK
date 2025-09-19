clear variables
close all

%% Prep workspace
outDir   = mfilename;
if exist(outDir, 'dir')
    rmdir(outDir, 's')
end
mkdir(outDir)

c1 = brewermap(320, 'BrBG');
c1 = c1(33:288,:);
c2 = brewermap(320, 'RdBu');
c2 = c2(33:288,:);
c3 = brewermap(320, 'PRGn');
c3 = c3(33:288,:);
cMaps = {c1, c2, c3, brewermap(128, 'Greys'), ...
         brewermap(7, 'Set1'),   brewermap(128, 'PuBuGn')};

%%
load('fig2_acqPerf_Results-W10L-38VYN53-20230720_172310.mat');

tVee = tVei * (1 - tVb);
[meshFa, meshTR] = meshgrid(tFa, tTR);
kern = ones(3)/9;

mnKve = mean(resKve, 3);
sdKve = std(resKve, 0, 3);
cvKve = 100 * (sdKve ./ mnKve);
erKve = 100 * (mnKve - tKve) ./ tKve;

mnVb = mean(resVb, 3);
sdVb = std(resVb, 0, 3);
cvVb = 100 * (sdVb ./ mnVb);
erVb = 100 * (mnVb - tVb) ./ tVb;

mnVee = mean(resVee, 3);
sdVee = std(resVee, 0, 3);
cvVee = 100 * (sdVee ./ mnVee);
erVee = 100 * (mnVee - tVee) ./ tVee;

mnResnm = mean(resnms, 3);

mdEflag = mode(eflags, 3);

%% Fig2
fig2 = figure('Position', [10 10 1300 1300]);
tiledlayout(2, 2);

% UL - kve error
nexttile;
imagesc(tFa, tTR, erKve);
set(gca, 'ydir', 'normal');
hold on
eLvl = [-50, -20, -10, 10, 20, 50, 200];
[cMat,cObj] = contour(meshFa, meshTR, imfilter(erKve, kern), eLvl, ...
    'LineColor', 'k', 'LineWidth', 2);
clabel(cMat, cObj, 'fontsize', 12);
eLim = [-50, 50];
caxis(eLim)
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('k_v_e  Mean Error (%)', 'FontSize', 20)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 18)
colormap(gca, cMaps{1})
colorbar(gca, 'Ticks', -50:25:50);

% UR - kve cv
nexttile;
imagesc(tFa, tTR, cvKve);
set(gca, 'ydir', 'normal');
hold on
cLvl = [10, 20, 50, 200];
[cMat,cObj] = contour(meshFa, meshTR, imfilter(cvKve, kern), cLvl, ...
    'LineColor', 'k', 'LineWidth', 2);
clabel(cMat, cObj, 'fontsize', 12);
cLim = [0, 50];
caxis(cLim)
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('k_v_e  Coefficient of Variation (%)', 'FontSize', 20)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 18)
colormap(gca, cMaps{1}(129:end,:))
colorbar(gca, 'Ticks', 0:10:50);

% LL - vb error
nexttile;
imagesc(tFa, tTR, erVb);
set(gca, 'ydir', 'normal');
hold on
eLvl = [-50, -20, -10, 10, 20, 50];
[cMat,cObj] = contour(meshFa, meshTR, imfilter(erVb, kern), eLvl, ...
    'LineColor', 'k', 'LineWidth', 2);
clabel(cMat, cObj, 'fontsize', 12);
eLim = [-50, 50];
caxis(eLim)
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('v_b  Mean Error (%)', 'FontSize', 20)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 18)
colormap(gca, cMaps{2})
colorbar(gca, 'Ticks', -50:25:50);

% LR - vb cv
nexttile;
imagesc(tFa, tTR, cvVb);
set(gca, 'ydir', 'normal');
hold on
cLvl = [10, 20, 50];
[cMat,cObj] = contour(meshFa, meshTR, imfilter(cvVb, kern), cLvl, ...
    'LineColor', 'k', 'LineWidth', 2);
clabel(cMat, cObj, 'fontsize', 12);
cLim = [0, 50];
caxis(cLim)
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('v_b  Coefficient of Variation (%)', 'FontSize', 20)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 18)
colormap(gca, flipud(cMaps{2}(1:128,:)))
colorbar(gca, 'Ticks', 0:10:50);

% subfig labels
annotation('textbox', [0.0350 0.9245 0.0683 0.0538], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4991 0.9245 0.0683 0.0538], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.0350 0.4561 0.0683 0.0538], 'String', 'C', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4991 0.4561 0.0683 0.0538], 'String', 'D', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

export_fig(fullfile(outDir, 'MRM_fig2_kve_vb_acqPerf'), '-tif', '-a1', '-nocrop', fig2);
savefig(fig2, fullfile(outDir, 'MRM_fig2_kve_vb_acqPerf'));


%% FigS1
figs1 = figure('Position', [10 10 1300 650]);
tiledlayout(1, 2);

% L - vee error
nexttile;
imagesc(tFa, tTR, erVee);
set(gca, 'ydir', 'normal');
hold on
eLvl = [-50, -20, 10, 20, 50];
[cMat,cObj] = contour(meshFa, meshTR, imfilter(erVee, kern), eLvl, ...
    'LineColor', 'k', 'LineWidth', 2);
clabel(cMat, cObj, 'fontsize', 16);
eLim = [-50, 50];
caxis(eLim)
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('v_e_e  Mean Error (%)', 'FontSize', 20)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 18)
colormap(gca, cMaps{3})
colorbar(gca, 'Ticks', -50:25:50);

% R - vee cv
nexttile;
imagesc(tFa, tTR, cvVee);
set(gca, 'ydir', 'normal');
hold on
cLvl = [10, 20, 50, 100];
[cMat,cObj] = contour(meshFa, meshTR, imfilter(cvVee, kern), cLvl, ...
    'LineColor', 'k', 'LineWidth', 2);
clabel(cMat, cObj, 'fontsize', 16);
cLim = [0, 50];
caxis(cLim)
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('v_e_e  Coefficient of Variation (%)', 'FontSize', 20)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 18)
colormap(gca, flipud(cMaps{3}(1:128,:)))
colorbar(gca, 'Ticks', 0:10:50);

% subfig labels
annotation('textbox', [0.0350 0.9245 0.0683 0.0538], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4991 0.9245 0.0683 0.0538], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

export_fig(fullfile(outDir, 'MRM_figs1_vee_acqPerf'), '-tif', '-a1', '-nocrop', figs1);
savefig(figs1, fullfile(outDir, 'MRM_figs1_vee_acqPerf'));


close all
return


%% resnms, eflags, avgtime
return

figz1 = figure;
imagesc(tFa, tTR, mnResnm);
set(gca, 'ydir', 'normal');
hold on
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('Mean Norm of Residuals per Pt', 'FontSize', 16)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 14)
colormap(gca, cMaps{4})
colorbar
export_fig(fullfile(outDir, 'figz1_resnm_acqPerf'), '-tif', '-nocrop', figz1);
savefig(figz1, fullfile(outDir, 'figz1_resnm_acqPerf'));

figz2 = figure;
imagesc(tFa, tTR, mdEflag);
set(gca, 'ydir', 'normal');
hold on
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('Mode of lsqcurvefit Exit Flag', 'FontSize', 16)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 14)
colormap(gca, cMaps{5})
caxis([-2 4])
colorbar
export_fig(fullfile(outDir, 'figz2_eflag_acqPerf'), '-tif', '-nocrop', figz2);
savefig(figz2, fullfile(outDir, 'figz2_eflag_acqPerf'));

figz3 = figure;
imagesc(tFa, tTR, avgTimesMsec);
set(gca, 'ydir', 'normal');
hold on
plot(20, 1, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'k');
title('Mean Fitting Time (ms)', 'FontSize', 16)
xlabel('Excitation Angle (deg)')
ylabel('Repetition Time (s)')
set(gca, 'fontsize', 14)
colormap(gca, cMaps{6})
colorbar
export_fig(fullfile(outDir, 'figz3_avgtime_acqPerf'), '-tif', '-nocrop', figz3);
savefig(figz3, fullfile(outDir, 'figz3_avgtime_acqPerf'));

