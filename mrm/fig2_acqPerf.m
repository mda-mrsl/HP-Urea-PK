% Simulation script for comparing acquisition methods with noise
%   Runtime is ~8hrs on laptop with i7-10610U processor (using parfor)
clear variables
close all

rng(1119, 'twister'); % for reproducibility

%% Simulation parameters
% Acquisition parameters
tAcq = 'conventional';
tTR  = 0.1:0.1:3;
tFa  = 1:1:90;

% Model parameters
tVb  = 0.09;
tVei = 0.3 / (1-tVb);
tKve = 0.02;

% Coarse ranges of parameters to be fit, for initial guesses
gVb    = [0.001, 0.01, 0.1];
gVei   = [0.005, 0.05, 0.5];
gKve   = [0.001, 0.01, 0.1];
gSize  = [numel(gVb), numel(gVei), numel(gKve)];
nGuess = prod(gSize);

% Noise parameters
tSnr     = 25;
refSig   = max(runPerfModel(struct('TR', 1, 'fa', 20)));

% Other simulation parameters
nReps      = 100;
opts       = optimoptions(@lsqcurvefit, 'Display', 'off', ...
    'MaxFunctionEvaluations', 3e4, 'MaxIterations', 3e4, ...
    'FunctionTolerance', 1e-9, 'OptimalityTolerance', 1e-9, ...
    'StepTolerance', 1e-9);
resSize    = [numel(tTR), numel(tFa)];
nIter      = prod(resSize);
startTime  = now;
timeString = datestr(startTime, 'yyyymmdd_HHMMSS');

% Initialize arrays for results
resVb        = zeros(nIter, nReps);
resVei       = zeros(nIter, nReps);
resVee       = zeros(nIter, nReps);
resKve       = zeros(nIter, nReps);
resnms       = zeros(nIter, nReps);
eflags       = zeros(nIter, nReps);
avgTimesMsec = zeros(nIter, 1);

% Randomize iteration order, helps with progress monitoring
iterOrder    = randperm(nIter);
[~, reOrder] = sort(iterOrder);

%% Run simulations
% Setup parallel pool, profiling and progress monitoring
fprintf('Starting %d iterations at %s. \n', nIter, datestr(now))
nCores  = feature('numCores');
poolObj = parpool(nCores);
p       = Par(nIter);
parfor_progress(nIter);
parfor ii = 1:nIter
    Par.tic;
    
    % Assemble acquisition and model parameter structures
    [a, b] = ind2sub(resSize, iterOrder(ii)); % tr, fa
    acqParams = struct( ...
        'acqScheme',    tAcq, ...
        'tEnd',         60,	...         % [s]
        'TR',           tTR(a),	...     % [s]
        'fa',           tFa(b),	...     % [deg]
        'errFa',        1,	...
        'tSegment',     10,	...         % [s]
        'faSegment',    90,	...         % [deg]
        'bval',         30, ...         % [s/mm^2]
        'TE',           0);             % [ms]
    t = 0 : acqParams.TR : acqParams.tEnd;
    nPts = numel(t);
    vif = segPolyfit(t, @(t) vifFunction(t, 100, 0.8, 8.1), 4);
    modelParams = struct( ...
        'vb',       tVb,     ... % [v/v]
        'vei',      tVei,    ... % [v/v]
        'kve',      tKve,    ... % [1/s]
        'T1',       20,      ... % [s]
        'vif',      vif,     ...
        'M0',       1e-6,    ...
        'De',       5e-4,    ... % [mm^2/s]
        'DvFactor', 30,      ...
        'T2',       inf);        % [ms]
    
    % Compute signal timecourses for coarse ranges of fit parameters
    gMx = zeros([gSize, nPts]);
    for jj = 1:nGuess
        [x, y, z]       = ind2sub(gSize, jj);
        guessParams     = modelParams;
        guessParams.vb  = gVb(x);
        guessParams.vei = gVei(y);
        guessParams.kve = gKve(z);
        gMx(x, y, z, :) = runPerfModel(acqParams, guessParams);
    end
    
    % Fit noisy signal timecourses
    trueMx = runPerfModel(acqParams, modelParams);
    tic
    for kk=1:nReps
        noisyMx   = trueMx + randn(nPts, 1) * refSig / tSnr;
        residMx   = gMx - repmat(reshape(noisyMx, 1, 1, 1, []), ...
                                 [gSize, 1]);
        ssrMx     = sum(residMx.^2, 4);
        [~, idx]  = min(ssrMx(:));
        [x, y, z] = ind2sub(gSize, idx);
        fitPars   = struct( ...
            'vb',  struct('guess', gVb(x),  'lower', 1e-4, 'upper', 0.9999), ... % [v/v]
            'vei', struct('guess', gVei(y), 'lower', 1e-4, 'upper', 1),      ... % [v/v]
            'kve', struct('guess', gKve(z), 'lower', 1e-4, 'upper', 10));    ... % [1/s]
            
        try
            [fits, curResnm, ~, curEflag] = fitPerfModel(noisyMx, fitPars, ...
                acqParams, modelParams, opts);
            resVb(ii, kk)  = fits.vb;
            resVei(ii, kk) = fits.vei;
            resVee(ii, kk) = fits.vei*(1-fits.vb);
            resKve(ii, kk) = fits.kve;
            resnms(ii, kk) = curResnm / nPts;
            eflags(ii, kk) = curEflag;
        catch
            resVb(ii, kk)  = NaN;
            resVei(ii, kk) = NaN;
            resVee(ii, kk) = NaN;
            resKve(ii, kk) = NaN;
            resnms(ii, kk) = NaN;
            eflags(ii, kk) = NaN;
        end
    end
    avgTimesMsec(ii) = 1e3 * toc / nReps;
    
    parfor_progress;
    p(ii) = Par.toc;
end
stop(p);
parfor_progress(0);
delete(poolObj);
endTime = now;
fprintf('Finished %d iterations at %s. \n', nIter, datestr(endTime))

% Reshape and de-randomize result arrays
resVb        = reshape(resVb(reOrder,:),  [resSize, nReps]);
resVei       = reshape(resVei(reOrder,:), [resSize, nReps]);
resVee       = reshape(resVee(reOrder,:), [resSize, nReps]);
resKve       = reshape(resKve(reOrder,:), [resSize, nReps]);
resnms       = reshape(resnms(reOrder,:), [resSize, nReps]);
eflags       = reshape(eflags(reOrder,:), [resSize, nReps]);
avgTimesMsec = reshape(avgTimesMsec(reOrder), resSize);

% return

% Save workspace
if ispc
    compString = getenv('COMPUTERNAME');
else
    compString = getenv('HOSTNAME');
end
save(sprintf('%s_Results-%s-%s.mat', mfilename, compString, timeString))

