% Sim script for comparing accuracy performance:
%   across models (3-param, 2-param) and model param values (kve, ve, vb)
%   Runtime is <1min on laptop with i7-10610U processor
clear variables
close all

%% Simulation parameters
% Acquisition parameters
tAcq = 'conventional';
tTR  = 1;
tFa  = 20;

% Model parameters
%  Test one param at a time, others kept fixed
%  fixed parameter values:
fVb  = 0.09;
fVee = 0.3;
fKve = 0.02;
%  test parameter ranges: (25 each, must be consistent)
tvb  = 0.02:0.02:0.5;
tvee = 0.02:0.02:0.5;
tkve = 0.002:0.002:0.05;

% Coarse ranges of parameters to be fit, for initial guesses
gVb     = [0.001, 0.01, 0.1];
gVei    = [0.005, 0.05, 0.5];
gKve    = [0.001, 0.01, 0.1];
gSize1  = [numel(gVb), numel(gVei), numel(gKve)];
nGuess1 = prod(gSize1);
gSize2  = [numel(gVei), numel(gKve)];
nGuess2 = prod(gSize2);
gSize3  = [numel(gVb), numel(gKve)];
nGuess3 = prod(gSize3);

% Other simulation parameters
opts       = optimoptions(@lsqcurvefit, 'Display', 'off', ...
    'MaxFunctionEvaluations', 3e4, 'MaxIterations', 3e4, ...
    'FunctionTolerance', 1e-9, 'OptimalityTolerance', 1e-9, ...
    'StepTolerance', 1e-9);
nIter      = numel(tvb);
startTime  = now;
timeString = datestr(startTime, 'yyyymmdd_HHMMSS');

% Initialize arrays for results
% row 1: 3-param, row 2: 2-param (kve, vec), row 3: 2-param (kve, vb)
% col 1: test vb, col 2: test vee, col 3: test kve
resVb	= zeros(3, 3, nIter);
resVei	= zeros(3, 3, nIter);
resVe 	= zeros(3, 3, nIter);
resKve	= zeros(3, 3, nIter);
resnms	= zeros(3, 3, nIter);
eflags	= zeros(3, 3, nIter);

%% Run simulations
% Assemble acquisition parameter structure
acqParams = struct( ...
    'acqScheme',    tAcq, ...
    'tEnd',         60,	...         % [s]
    'TR',           tTR,	...     % [s]
    'fa',           tFa,	...     % [deg]
    'errFa',        1,	...
    'tSegment',     10,	...         % [s]
    'faSegment',    90,	...         % [deg]
    'bval',         30, ...         % [s/mm^2]
    'TE',           0);             % [ms]
t    = 0 : acqParams.TR : acqParams.tEnd;
nPts = numel(t);
vif  = segPolyfit(t, @(t) vifFunction(t, 100, 0.8, 8.1), 4);

params = {'vb', 'vei', 'kve'};

fprintf('Starting %d iterations at %s. \n', nIter*3, datestr(now))
for jj = 1:3 % vb, vee, kve
    
    if jj == 1
        tvei = fVee ./ (1-tvb);
    elseif jj == 2
        tvei = tvee ./ (1-fVb);
    else
        tvei = ones(1, nIter) * fVee / (1-fVb); 
    end
    
    modelParams = struct( ...
        'vb',       fVb,  ... % [v/v]
        'vei',      NaN,  ... % [v/v]
        'kve',      fKve, ... % [1/s]
        'T1',       20,      ... % [s]
        'vif',      vif,     ...
        'M0',       1e-6,    ...
        'De',       5e-4,    ... % [mm^2/s]
        'DvFactor', 30,      ...
        'T2',       inf);        % [ms]
        
    for ii = 1:nIter
        
        eval(sprintf('modelParams.%s = t%s(ii);', params{jj}, params{jj}))
        modelParams.vei = tvei(ii);
        
        trueMx = runPerfModel(acqParams, modelParams);
        
        % 3-params, full model
        gMx1 = zeros([gSize1, nPts]);
        for kk = 1:nGuess1
            [x, y, z]       = ind2sub(gSize1, kk);
            guessParams     = modelParams;
            guessParams.vb  = gVb(x);
            guessParams.vei = gVei(y);
            guessParams.kve = gKve(z);
            gMx1(x, y, z, :) = runPerfModel(acqParams, guessParams);
        end
        residMx   = gMx1 - repmat(reshape(trueMx, 1, 1, 1, []), ...
            [gSize1, 1]);
        ssrMx     = sum(residMx.^2, 4);
        [~, idx]  = min(ssrMx(:));
        [x, y, z] = ind2sub(gSize1, idx);
        fitPars   = struct( ...
            'vb',  struct('guess', gVb(x),  'lower', 1e-4, 'upper', 0.9999), ... % [v/v]
            'vei', struct('guess', gVei(y), 'lower', 1e-4, 'upper', 1),      ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10));    ... % [1/s]
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel(trueMx, fitPars, ...
                acqParams, modelParams, opts);
            resVb(1, jj, ii)  = fits.vb;
            resVei(1, jj, ii) = fits.vei;
            resVe(1, jj, ii)  = fits.vei*(1-fits.vb);
            resKve(1, jj, ii) = fits.kve;
            resnms(1, jj, ii) = curResnm / nPts;
            eflags(1, jj, ii) = curEflag;
        catch
            resVb(1, jj, ii)  = NaN;
            resVei(1, jj, ii) = NaN;
            resVe(1, jj, ii)  = NaN;
            resKve(1, jj, ii) = NaN;
            resnms(1, jj, ii) = NaN;
            eflags(1, jj, ii) = NaN;
        end

        % 2-params, (kve, vei), vb = 0
        gMx2 = zeros([gSize2, nPts]);
        for kk = 1:nGuess2
            [y, z]          = ind2sub(gSize2, kk);
            guessParams     = modelParams;
            guessParams.vb  = 0;
            guessParams.vei = gVei(y);
            guessParams.kve = gKve(z);
            gMx2(y, z, :)   = runPerfModel(acqParams, guessParams);
        end
        residMx   = gMx2 - repmat(reshape(trueMx, 1, 1, []), ...
            [gSize2, 1]);
        ssrMx     = sum(residMx.^2, 3);
        [~, idx]  = min(ssrMx(:));
        [y, z] = ind2sub(gSize2, idx);
        fitPars   = struct( ...
            'vei', struct('guess', gVei(y), 'lower', 1e-4, 'upper', 1),      ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10));    ... % [1/s]
        fixPars    = modelParams;
        fixPars.vb = 0;
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel(trueMx, fitPars, ...
                acqParams, fixPars, opts);
            resVb(2, jj, ii)  = 0;
            resVei(2, jj, ii) = fits.vei;
            resVe(2, jj, ii)  = fits.vei;
            resKve(2, jj, ii) = fits.kve;
            resnms(2, jj, ii) = curResnm / nPts;
            eflags(2, jj, ii) = curEflag;
        catch
            resVb(2, jj, ii)  = NaN;
            resVei(2, jj, ii) = NaN;
            resVe(2, jj, ii)  = NaN;
            resKve(2, jj, ii) = NaN;
            resnms(2, jj, ii) = NaN;
            eflags(2, jj, ii) = NaN;
        end

        % 2-params, (kve, vb), vei = 1
        gMx3 = zeros([gSize3, nPts]);
        for kk = 1:nGuess3
            [x, z]          = ind2sub(gSize3, kk);
            guessParams     = modelParams;
            guessParams.vb  = gVb(x);
            guessParams.vei = 1;
            guessParams.kve = gKve(z);
            gMx3(x, z, :)   = runPerfModel(acqParams, guessParams);
        end
        residMx   = gMx3 - repmat(reshape(trueMx, 1, 1, []), ...
            [gSize3, 1]);
        ssrMx     = sum(residMx.^2, 3);
        [~, idx]  = min(ssrMx(:));
        [x, z] = ind2sub(gSize3, idx);
        fitPars   = struct( ...
            'vb',  struct('guess', gVb(x),  'lower', 1e-4, 'upper', 0.9999), ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10));    ... % [1/s]
        fixPars     = modelParams;
        fixPars.vei = 1;
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel(trueMx, fitPars, ...
                acqParams, fixPars, opts);
            resVb(3, jj, ii)  = fits.vb;
            resVei(3, jj, ii) = 1;
            resVe(3, jj, ii)  = 1-fits.vb;
            resKve(3, jj, ii) = fits.kve;
            resnms(3, jj, ii) = curResnm / nPts;
            eflags(3, jj, ii) = curEflag;
        catch
            resVb(3, jj, ii)  = NaN;
            resVei(3, jj, ii) = NaN;
            resVe(3, jj, ii)  = NaN;
            resKve(3, jj, ii) = NaN;
            resnms(3, jj, ii) = NaN;
            eflags(3, jj, ii) = NaN;
        end
        
    end
    clear tvei
end
endTime = now;
fprintf('Finished %d iterations at %s. \n', nIter*3, datestr(endTime))

% return

% % Save workspace
% if ispc
%     compString = getenv('COMPUTERNAME');
% else
%     compString = getenv('HOSTNAME');
% end
% save(sprintf('%s_Results-%s-%s.mat', mfilename, compString, timeString))
% 
% return

%% Plot results
outDir = mfilename;
if exist(outDir, 'dir')
    rmdir(outDir, 's')
end
mkdir(outDir)

cm = lines(5);
close all

%% Fig3
fig3 = figure('Position', [10 10 1900 650]);
tiledlayout(1, 3);

% True kve on x-axis, Fit kve on y-axis
nexttile;
plot(tkve, squeeze(resKve(1,3,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tkve, squeeze(resKve(2,3,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
plot(tkve, squeeze(resKve(3,3,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
plot(tkve, tkve, 'k--', 'linewidth', 1.1);
xlabel('True k_v_e (s^-^1)')
ylabel('Fit k_v_e (s^-^1)')
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', 'location', 'northeast')
set(gca, 'fontsize', 16);
axis([0 0.05 0 0.23])

% True vb on x-axis, Fit kve on y-axis
nexttile;
semilogy(tvb, squeeze(resKve(1,1,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
semilogy(tvb, squeeze(resKve(2,1,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
semilogy(tvb, squeeze(resKve(3,1,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.5], fKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_b (v/v)')
ylabel('Fit k_v_e (s^-^1)')
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', 'location', 'northwest')
set(gca, 'fontsize', 16);
xlim([0 0.5])

% True ve on x-axis, Fit kve on y-axis
nexttile;
semilogy(tvee, squeeze(resKve(1,2,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
semilogy(tvee, squeeze(resKve(2,2,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
semilogy(tvee, squeeze(resKve(3,2,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.5], fKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_e_e (v/v)')
ylabel('Fit k_v_e (s^-^1)')
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', 'location', 'southeast')
set(gca, 'fontsize', 16);
xlim([0 0.5])

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

export_fig(fullfile(outDir, 'MRM_fig3_paramSensitivity'), '-tif', '-a1', '-nocrop', fig3);
savefig(fig3, fullfile(outDir, 'MRM_fig3_paramSensitivity'));


%% FigS2
figs2 = figure('Position', [10 10 1900 1300]);
tiledlayout(2, 3);

% UL- True kve on x-axis, Fit vb on y-axis
nexttile;
plot(tkve, squeeze(resVb(1,3,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tkve, squeeze(resVb(3,3,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.05], fVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True k_v_e (s^-^1)')
ylabel('Fit v_b (v/v)')
grid on
legend('Model I', 'Model III', 'True v_b', 'location', 'northwest')
set(gca, 'fontsize', 16);
axis([0 0.05 0.08 0.17])

% UC- True vb on x-axis, Fit vb on y-axis
nexttile;
plot(tvb, squeeze(resVb(1,1,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tvb, squeeze(resVb(3,1,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
plot(tvb, tvb, 'k--', 'linewidth', 1.1);
xlabel('True v_b (v/v)')
ylabel('Fit v_b (v/v)')
grid on
legend('Model I', 'Model III', 'True v_b', 'location', 'northwest')
set(gca, 'fontsize', 16);
axis([0 0.5 0 0.55])

% UR- True vee on x-axis, Fit vb on y-axis
nexttile;
plot(tvee, squeeze(resVb(1,2,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tvee, squeeze(resVb(3,2,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.5], fVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_e_e (v/v)')
ylabel('Fit v_b (v/v)')
grid on
legend('Model I', 'Model III', 'True v_b', 'location', 'northeast')
set(gca, 'fontsize', 16);
axis([0 0.5 0.08 0.13])

% LL- True kve on x-axis, Fit ve on y-axis
nexttile;
plot(tkve, squeeze(resVe(1,3,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tkve, squeeze(resVe(2,3,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
plot(tkve, squeeze(resVe(3,3,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.05], fVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True k_v_e (s^-^1)')
ylabel('Fit v_e (v/v)')
grid on
legend('Model I (v_e_e)', 'Model II (v_e_c)', 'Model III (v_e_v)', ...
    'True v_e_e', 'location', 'west')
set(gca, 'fontsize', 16);
axis([0 0.05 0 1])

% LC- True vb on x-axis, Fit ve on y-axis
nexttile;
plot(tvb, squeeze(resVe(1,1,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tvb, squeeze(resVe(2,1,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
plot(tvb, squeeze(resVe(3,1,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.5], fVee*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_b (v/v)')
ylabel('Fit v_e (v/v)')
grid on
legend('Model I (v_e_e)', 'Model II (v_e_c)', 'Model III (v_e_v)', ...
    'True v_e_e', 'location', 'northeast')
set(gca, 'fontsize', 16);
axis([0 0.5 0 1])

% LR- True vee on x-axis, Fit ve on y-axis
nexttile;
plot(tvee, squeeze(resVe(1,2,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tvee, squeeze(resVe(2,2,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
plot(tvee, squeeze(resVe(3,2,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
plot(tvee, tvee, 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_e_e (v/v)')
ylabel('Fit v_e (v/v)')
grid on
legend('Model I (v_e_e)', 'Model II (v_e_c)', 'Model III (v_e_v)', ...
    'True v_e_e', 'location', 'west')
set(gca, 'fontsize', 16);
axis([0 0.5 0 1])

% subfig labels
% annotation('textbox', [0.0252 0.9196 0.0467 0.1076], 'String', 'A', ...
%     'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
%     'LineStyle', 'none');
% annotation('textbox', [0.3367 0.9174 0.0467 0.1076], 'String', 'B', ...
%     'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
%     'LineStyle', 'none');
% annotation('textbox', [0.6476 0.9183 0.0467 0.1076], 'String', 'C', ...
%     'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
%     'LineStyle', 'none');

export_fig(fullfile(outDir, 'MRM_figS2_paramSensitivity'), '-tif', '-a1', '-nocrop', figs2);
savefig(figs2, fullfile(outDir, 'MRM_figS2_paramSensitivity'));


close all
return


%% Fig3v1
fig3v1 = figure('Position', [10 10 1300 1300]);
tiledlayout(2, 2);

% True kve on x-axis, Fit kve on y-axis
nexttile;
plot(tkve, squeeze(resKve(1,3,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tkve, squeeze(resKve(2,3,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
plot(tkve, squeeze(resKve(3,3,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
plot(tkve, tkve, 'k--', 'linewidth', 1.1);
xlabel('True k_v_e (s^-^1)')
ylabel('Fit k_v_e (s^-^1)')
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', 'location', 'northeast')
set(gca, 'fontsize', 16);
axis([0 0.05 0 0.25])

% True kve on x-axis, Fit vb on y-axis
nexttile;
plot(tkve, squeeze(resVb(1,3,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
plot(tkve, squeeze(resVb(3,3,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.05], fVb*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True k_v_e (s^-^1)')
ylabel('Fit v_b (v/v)')
grid on
legend('Model I', 'Model III', 'True v_b', 'location', 'northwest')
set(gca, 'fontsize', 16);

% True vb on x-axis, Fit kve on y-axis
nexttile;
semilogy(tvb, squeeze(resKve(1,1,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
semilogy(tvb, squeeze(resKve(2,1,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
semilogy(tvb, squeeze(resKve(3,1,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.5], fKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_b (v/v)')
ylabel('Fit k_v_e (s^-^1)')
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', 'location', 'northwest')
set(gca, 'fontsize', 16);

% True ve on x-axis, Fit kve on y-axis
nexttile;
semilogy(tvee, squeeze(resKve(1,2,:)), 'o-', 'LineWidth', 1.3, 'Color', cm(5,:));
hold on
semilogy(tvee, squeeze(resKve(2,2,:)), 's-', 'LineWidth', 1.3, 'Color', cm(1,:), 'MarkerSize', 7);
semilogy(tvee, squeeze(resKve(3,2,:)), 'd-', 'LineWidth', 1.3, 'Color', cm(4,:));
line([0 0.5], fKve*[1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 1.1);
xlabel('True v_e_e (v/v)')
ylabel('Fit k_v_e (s^-^1)')
grid on
legend('Model I', 'Model II', 'Model III', 'True k_v_e', 'location', 'southeast')
set(gca, 'fontsize', 16);
xlim([0 0.5])

% subfig labels
annotation('textbox', [0.0392 0.9251 0.0640 0.0538], 'String', 'A', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4947 0.9245 0.0640 0.0538], 'String', 'B', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.0375 0.4626 0.0658 0.0538], 'String', 'C', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');
annotation('textbox', [0.4924 0.4622 0.0658 0.0538], 'String', 'D', ...
    'FitBoxToText', 'on', 'FontSize', 36, 'FontWeight', 'bold', ...
    'LineStyle', 'none');

export_fig(fullfile(outDir, 'MRM_fig3v1_paramSensitivity'), '-tif', '-nocrop', fig3v1);
savefig(fig3v1, fullfile(outDir, 'MRM_fig3v1_paramSensitivity'));

