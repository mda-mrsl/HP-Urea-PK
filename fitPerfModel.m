function [modelParams, resnorm, residual, exitflag, output] = ...
            fitPerfModel(Mx, fitParams, acqParams, fixParams, options)

% Fit hyperpolarized perfusion model to transverse signal timecourse
%##########################################################################
%#  [modelParams, resnorm, residual, exitflag, output] = ...
%#     fitPerfModel(Mx, fitParams, acqParams, fixParams, options)
%#
%#  Fit perfusion model to given transverse signal timecourse obtained
%#  with given acquisition parameters. Physical and physiological model
%#  parameters are fit within specified bounds or fixed during fitting.
%#
%#  REQUIRED INPUTS:
%#      Mx        - Vector of transverse magnetization over time
%#      fitParams - Structure specifying parameters to be fit
%#        > Each field should be a substructure with fields:
%#          guess - Initial guess for parameter value
%#          lower - Lower bound value imposed during fitting
%#          upper - Upper bound value imposed during fitting
%#            > If one or more of these substructure fields are not given,
%#              defaults for the specified parameter are used (see below).
%#        > Currently, only the following parameters should be fit:
%#          vb  - Vascular tissue volume fraction
%#            Defaults: guess = 0.09, lower = 1e-4, upper = 0.9999
%#          vei - Distribution volume fraction of extravascular space
%#            Defaults: guess = 0.3,  lower = 1e-4, upper = 1
%#          kve - Exchange rate btwn vascular and extrascular spaces  [1/s]
%#            Defaults: guess = 0.02, lower = 1e-4, upper = 10
%#          T1  - Time constant of losses other than RF and exchange  [s]
%#            Defaults: guess = 20,   lower = 1,    upper = 1e4
%#
%#  OPTIONAL INPUTS:
%#      acqParams - Structure containing acquisition parameters
%#        > Acquisition parameters not given in acqParams are set to
%#          default values (see parseAcq).
%#      fixParams - Structure containing parameters fixed during fitting
%#        > Model parameters not given in fitParams or fixParams are fixed
%#          during fitting at default values (see parseModel).
%#        > If the same parameter is given in fitParams and fixParams it
%#          will be fit.
%#      options   - Optimization options structure, passed to lsqcurvefit
%#
%#  OUTPUTS:
%#      modelParams - Structure containing all model parameters
%#        > Fixed parameters are passed directly into modelParams.
%#      resnorm     - Squared 2-norm of the residuals for the fit model
%#      residual    - Residuals for the fit model
%#      exitflag    - Value describing lsqcurvefit exit condition
%#      output      - Structure containing lsqcurvefit outputs
%#
%#  LITERATURE:
%#      Bankson JA, et al. Cancer Res 2015.
%#        DOI 10.1158/0008-5472.CAN-15-0171
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 04/2018
%##########################################################################

if nargin<2, help(mfilename); error('fitPerfModel:inputs', ...
                                    'Not enough inputs.'); end

%% Parse inputs
if nargin<5,            options = []; end
if isempty(options),	options = optimoptions(@lsqcurvefit); end
if nargin<4,            fixParams = []; end
if nargin<3,            acqParams = []; end
parseAcq(acqParams)
assert(numel(Mx)==numel(t), 'fitPerfModel:MxSize', ...
    'Mx input does not match given acquisition parameters')

defaultFitParams = struct( ...
    'vb',  struct('guess', 0.09, 'lower', 1e-4, 'upper', 0.9999), ... % [v/v]
    'vei', struct('guess', 0.3,  'lower', 1e-4, 'upper', 1),      ... % [v/v]
    'kve', struct('guess', 0.02, 'lower', 1e-4, 'upper', 10),     ... % [1/s]
    'T1',  struct('guess', 20,   'lower', 1,    'upper', 1e4));       % [s]
defaultFields = fieldnames(defaultFitParams);
fitFields   = fieldnames(fitParams);
for field = 1:numel(fitFields)
    curField = fitFields{field};
    if ~any(strcmp(curField, defaultFields))
        if ~any(strcmpi(curField, defaultFields))
            warning('fitPerfModel:fieldName', ...
                'fitParams field "%s" not recognized', curField)
        else
            match = defaultFields{strcmpi(curField, defaultFields)};
            warning('fitPerfModel:fieldNameCase', ['Case insensitive ', ...
                'match for fitParams field "%s"'], match)
            fitParams.(match) = fitParams.(curField);
            fitParams = rmfield(fitParams, curField);
        end
    end
end

defaultSubFields = fieldnames(defaultFitParams.(defaultFields{1}));
fitFields = fieldnames(fitParams);
guess     = zeros(numel(fitFields), 1);
lower     = zeros(numel(fitFields), 1);
upper     = zeros(numel(fitFields), 1);
for field = 1:numel(fitFields)               % vb, vei, kve and/or T1
    curField = fitFields{field};
    for subField = 1:numel(defaultSubFields) % guess, lower and upper
        curSubField = defaultSubFields{subField};
        if ~isstruct(fitParams.(curField))
            warning('fitPerfModel:fitParamsSubstruct', ['fitParams field ', ...
                '"%s" is not a substruct. Using default initial guess ', ...
                'and bounds for this parameter'], curField)
            fitParams.(curField) = defaultFitParams.(curField);
        else
            if ~any(strcmpi(curSubField, fieldnames(fitParams.(curField))))
                warning('fitPerfModel:fitParamsSubfield', ...
                    ['fitParams substruct "%s" subfield "%s" not ', ...
                    'given. Using default for this parameter'], ...
                    curField, curSubField)
                fitParams.(curField).(curSubField) = ...
                    defaultFitParams.(curField).(curSubField);
            end
        end
        evalString = sprintf( ...
            '%s(%d) = fitParams.(curField).(curSubField);', ...
            curSubField, field);
        eval(evalString)
    end
end

%% Fit model
fcn = @(X, T) fitFcn(acqParams, fixParams, fitFields, X);
[fitValues, resnorm, residual, exitflag, output] = ...
    lsqcurvefit(fcn, guess, t(:), Mx(:), lower, upper, options);

%% Assemble modelParams output
modelParams = fixParams;
for field = 1:numel(fitFields)
    curField = fitFields{field};
    modelParams.(curField) = fitValues(field);
end


%=========================================================================%

function Mx = fitFcn(acqParams, fixParams, fitFields, fitValues)
%   Combine fixed and fit parameters and evaluates perfModel

modelParams = fixParams;
for field = 1:numel(fitFields)
        curField = fitFields{field};
        modelParams.(curField) = fitValues(field);
end
acqParams.errFa = 1;    % Nominal flip angles
Mx = runPerfModel(acqParams, modelParams);
