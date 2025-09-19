% Sim script for comparing SNR performance across models (3-param, 2-param)
% VIF shape known and fixed, VIF scale factor is fit
%   Runtime is ~10mins on laptop with i7-10610U processor
clear variables
close all

rng(1234, 'twister'); % for reproducibility

%% Simulation parameters
% Acquisition parameters
tAcq = 'conventional';
tTR  = 1;
tFa  = 20;

% Model parameters
tVb  = 0.09;
tVee = 0.3;
tKve = 0.02;

tVei = tVee / (1-tVb);

% Coarse ranges of parameters to be fit, for initial guesses
gVb     = [0.001, 0.01, 0.1];
gVei    = [0.005, 0.05, 0.5];
gKve    = [0.001, 0.01, 0.1];
gVsc    = [0.4, 0.9, 1.4];
gSize1  = [numel(gVb), numel(gVei), numel(gKve), numel(gVsc)];
nGuess1 = prod(gSize1);
gSize2  = [numel(gVei), numel(gKve), numel(gVsc)];
nGuess2 = prod(gSize2);
gSize3  = [numel(gVb), numel(gKve), numel(gVsc)];
nGuess3 = prod(gSize3);

% Noise parameters
refSig = max(runPerfModel(struct('TR', 1, 'fa', 20)));
tSnr = 10:5:50;

% Other simulation parameters
nReps = 100;
opts       = optimoptions(@lsqcurvefit, 'Display', 'off', ...
    'MaxFunctionEvaluations', 3e4, 'MaxIterations', 3e4, ...
    'FunctionTolerance', 1e-9, 'OptimalityTolerance', 1e-9, ...
    'StepTolerance', 1e-9);
resSize    = numel(tSnr);
nIter      = prod(resSize);
startTime  = now;
timeString = datestr(startTime, 'yyyymmdd_HHMMSS');

% Initialize arrays for results
% row 1: 3-param, row 2: 2-param (kve, ve), row 3: 2-param (kve, vb)
resVb	= zeros(3, nIter, nReps);
resVei	= zeros(3, nIter, nReps);
resVe 	= zeros(3, nIter, nReps);
resKve	= zeros(3, nIter, nReps);
resnms	= zeros(3, nIter, nReps);
eflags	= zeros(3, nIter, nReps);
resVsc  = zeros(3, nIter, nReps);

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
modelParams = struct( ...
    'vb',       tVb,  ... % [v/v]
    'vei',      tVei, ... % [v/v]
    'kve',      tKve, ... % [1/s]
    'T1',       20,      ... % [s]
    'vif',      vif,     ...
    'M0',       1e-6,    ...
    'De',       5e-4,    ... % [mm^2/s]
    'DvFactor', 30,      ...
    'T2',       inf);        % [ms]
trueMx = runPerfModel(acqParams, modelParams);

% Compute signal timecourses for coarse ranges of fit parameters
gMx1 = zeros([gSize1, nPts]);
for jj = 1:nGuess1
    [x, y, z, v]    = ind2sub(gSize1, jj);
    guessParams     = modelParams;
    guessParams.vb  = gVb(x);
    guessParams.vei = gVei(y);
    guessParams.kve = gKve(z);
    tmp             = runPerfModel(acqParams, guessParams);
    gMx1(x, y, z, v, :) = gVsc(v) * tmp;
end
gMx2 = zeros([gSize2, nPts]);
for jj = 1:nGuess2
    [y, z, v]          = ind2sub(gSize2, jj);
    guessParams     = modelParams;
    guessParams.vb  = 0;
    guessParams.vei = gVei(y);
    guessParams.kve = gKve(z);
    tmp             = runPerfModel(acqParams, guessParams);
    gMx2(y, z, v, :) = gVsc(v) * tmp;
end
gMx3 = zeros([gSize3, nPts]);
for jj = 1:nGuess3
    [x, z, v]          = ind2sub(gSize3, jj);
    guessParams     = modelParams;
    guessParams.vb  = gVb(x);
    guessParams.vei = 1;
    guessParams.kve = gKve(z);
    tmp             = runPerfModel(acqParams, guessParams);
    gMx3(x, z, v, :) = gVsc(v) * tmp;
end

fprintf('Starting %d iterations at %s. \n', nIter, datestr(now))
for ii = 1:nIter   % SNR vals
    
    % Fit noisy signal timecourses
    for kk = 1:nReps
        
        noisyMx = trueMx + randn(nPts, 1) * refSig / tSnr(ii);
        
        % 3-params, full model
        residMx   = gMx1 - repmat(reshape(noisyMx, 1, 1, 1, 1, []), ...
            [gSize1, 1]);
        ssrMx     = sum(residMx.^2, 5);
        [~, idx]  = min(ssrMx(:));
        [x, y, z, v] = ind2sub(gSize1, idx);
        fitPars   = struct( ...
            'vb',  struct('guess', gVb(x),  'lower', 1e-4, 'upper', 0.9999), ... % [v/v]
            'vei', struct('guess', gVei(y), 'lower', 1e-4, 'upper', 1),      ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10),     ... % [1/s]
            'vsc', struct('guess', gVsc(v), 'lower', 1e-4, 'upper', inf));
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel2(noisyMx, fitPars, ...
                acqParams, modelParams, opts);
            resVb(1,ii,kk)  = fits.vb;
            resVei(1,ii,kk) = fits.vei;
            resVe(1,ii,kk)  = fits.vei * (1-fits.vb);
            resKve(1,ii,kk) = fits.kve;
            resnms(1,ii,kk) = curResnm / nPts;
            eflags(1,ii,kk) = curEflag;
            resVsc(1,ii,kk) = fits.vsc;
        catch
            resVb(1,ii,kk)  = NaN;
            resVei(1,ii,kk) = NaN;
            resVe(1,ii,kk)  = NaN;
            resKve(1,ii,kk) = NaN;
            resnms(1,ii,kk) = NaN;
            eflags(1,ii,kk) = NaN;
            resVsc(1,ii,kk) = NaN;
        end

        noisyMx = trueMx + randn(nPts, 1) * refSig / tSnr(ii);
        
        % 2-params, (kve, vei), vb = 0
        residMx   = gMx2 - repmat(reshape(noisyMx, 1, 1, 1, []), ...
            [gSize2, 1]);
        ssrMx     = sum(residMx.^2, 4);
        [~, idx]  = min(ssrMx(:));
        [y, z, v] = ind2sub(gSize2, idx);
        fitPars   = struct( ...
            'vei', struct('guess', gVei(y), 'lower', 1e-4, 'upper', 1),      ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10),     ... % [1/s]
            'vsc', struct('guess', gVsc(v), 'lower', 1e-4, 'upper', inf));
            fixPars    = modelParams;
        fixPars.vb = 0;
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel2(noisyMx, fitPars, ...
                acqParams, fixPars, opts);
            resVb(2,ii,kk)  = 0;
            resVei(2,ii,kk) = fits.vei;
            resVe(2,ii,kk)  = fits.vei;
            resKve(2,ii,kk) = fits.kve;
            resnms(2,ii,kk) = curResnm / nPts;
            eflags(2,ii,kk) = curEflag;
            resVsc(2,ii,kk) = fits.vsc;
        catch
            resVb(2,ii,kk)  = NaN;
            resVei(2,ii,kk) = NaN;
            resVe(2,ii,kk)  = NaN;
            resKve(2,ii,kk) = NaN;
            resnms(2,ii,kk) = NaN;
            eflags(2,ii,kk) = NaN;
        end

        noisyMx = trueMx + randn(nPts, 1) * refSig / tSnr(ii);
        
        % 2-params, (kve, vb), vei = 1
        residMx   = gMx3 - repmat(reshape(noisyMx, 1, 1, 1, []), ...
            [gSize3, 1]);
        ssrMx     = sum(residMx.^2, 4);
        [~, idx]  = min(ssrMx(:));
        [x, z, v] = ind2sub(gSize3, idx);
        fitPars   = struct( ...
            'vb',  struct('guess', gVb(x),  'lower', 1e-4, 'upper', 0.9999), ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10),     ... % [1/s]
            'vsc', struct('guess', gVsc(v), 'lower', 1e-4, 'upper', inf));
            fixPars     = modelParams;
        fixPars.vei = 1;
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel2(noisyMx, fitPars, ...
                acqParams, fixPars, opts);
            resVb(3,ii,kk)  = fits.vb;
            resVei(3,ii,kk) = 1;
            resVe(3,ii,kk)  = 1-fits.vb;
            resKve(3,ii,kk) = fits.kve;
            resnms(3,ii,kk) = curResnm / nPts;
            eflags(3,ii,kk) = curEflag;
            resVsc(3,ii,kk) = fits.vsc;
        catch
            resVb(3,ii,kk)  = NaN;
            resVei(3,ii,kk) = NaN;
            resVe(3,ii,kk)  = NaN;
            resKve(3,ii,kk) = NaN;
            resnms(3,ii,kk) = NaN;
            eflags(3,ii,kk) = NaN;
        end
    
    end

end
endTime = now;
fprintf('Finished %d iterations at %s. \n', nIter, datestr(endTime))

% Save workspace
if ispc
    compString = getenv('COMPUTERNAME');
else
    compString = getenv('HOSTNAME');
end
save(sprintf('%s_Results-%s-%s.mat', mfilename, compString, timeString))

return
