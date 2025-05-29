function parseModel(modelParams)

% Parse structure containing physical and physiological model parameters
%##########################################################################
%#  parseModel(modelParams)
%#
%#  Unpack structure containing physical and physiological model parameter
%#  values, replacing any missing parameters with defaults. No outputs
%#  should be assigned when using this function, variables are assigned by
%#  name into the workspace of the caller.
%#
%#  OPTIONAL INPUTS:
%#      modelParams - Structure containing model parameters
%#        > If empty or not provided, defaults are used for all parameters.
%#
%#  OUTPUTS: (assigned into workspace of caller)
%#      vb    - Vascular tissue volume fraction
%#        Default = 0.09
%#      vei   - Distribution volume fraction of extravascular space
%#        Default = 0.3297
%#      kve   - Exchange rate between vascular and extrascular spaces [1/s]
%#        Default = 0.02
%#      T1    - Time constant of losses other than RF and exchange    [s]
%#        Default = 20
%#      vif   - One of the following vascular input functions:
%#        > Substructure containing gamma-variate curve parameters passed
%#          to vifFunction, containing fields:
%#          vifMax   - Maximum amplitude of gamma-variate curve
%#          vifShape - Shape parameter of gamma-variate curve
%#          vifTMax  - Time at which maximum amplitude occurs
%#          vifDelay - Delay time to first point on gamma-variate curve
%#        > Handle to function evaluating Mz in blood over time.
%#        > Matrix of size [numel(t)-1, N], where 2<=N<=5. Rows contain
%#          coefficients for the (N-1)th order polynomial fit to the
%#          segments of the VIF in the corresponding TR intervals.
%#        > Vector of same length as t providing Mz in blood at each
%#          excitation time.
%#        Default = struct('vifMax', 100, 'vifShape', 0.8, ...
%#                         'vifTmax', 8.1, 'vifDelay', 0);
%#      M0    - Initial Mz in extravascular space at first excitation
%#        Default = 1e-6
%#      De    - Extravascular apparent diffusion coefficient  [mm^2/s]
%#        Default = 5e-4
%#      DvFactor  - Factor applied to De to obtain the pseudodiffusion
%#                  coefficient characterizing microvascular flow
%#        Default = 30
%#      T2        - Spin-spin relaxation time constant        [ms]
%#        Default = inf
%#
%#  LITERATURE:
%#      Bankson JA, et al. Cancer Res 2015.
%#        DOI 10.1158/0008-5472.CAN-15-0171
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 04/2018
%##########################################################################

if ~nargin, modelParams = []; end
if isempty(modelParams)
    modelParams = struct();
else
    assert(isa(modelParams, 'struct'), 'parseModel:input', ...
        'Input must be a structure')
end

defaultVif = struct( ...
    'vifMax',   100, ...
    'vifShape', 0.8, ...
    'vifTMax',  8.1, ... % [s]
    'vifDelay', 0);      % [s]

if isfield(modelParams, 'vif') && isstruct(modelParams.vif)
    defaultVifFields = fieldnames(defaultVif);
    inputVifFields   = fieldnames(modelParams.vif);
    for field = 1:numel(defaultVifFields)
        curField = defaultVifFields{field};
        if ~any(strcmp(curField, inputVifFields))
            if any(strcmpi(curField, inputVifFields))
                match = inputVifFields{strcmpi(curField, inputVifFields)};
                warning('parseModel:vifFieldNameCase', ...
                    'Case insensitive match for input VIF subfield "%s"', match)
                modelParams.vif.(curField) = modelParams.vif.(match);
                modelParams.vif = rmfield(modelParams.vif, match);
            else
                modelParams.vif.(curField) = defaultVif.(curField);
            end
        end
    end
    inputVifFields   = fieldnames(modelParams.vif);
    for field = 1:numel(inputVifFields)
        curField = inputVifFields{field};
        if ~any(strcmp(curField, defaultVifFields))
            warning('parseModel:vifFieldName', ...
                'Input VIF subfield "%s" not recognized', curField)
            modelParams.vif = rmfield(modelParams.vif, curField);
        end
    end
end

defaultParams = struct( ...
    'vb',       0.09,   ... % [v/v]
    'vei',      0.3297, ... % [v/v]
    'kve',      0.02,   ... % [1/s]
    'T1',       20,     ... % [s]
    'vif',      defaultVif, ...
    'M0',       1e-6,   ...
    'De',       5e-4,   ... % [mm^2/s]
    'DvFactor', 30,     ...
    'T2',       inf);       % [ms]

defaultFields = fieldnames(defaultParams);
inputFields   = fieldnames(modelParams);
for field = 1:numel(inputFields)
    curField = inputFields{field};
    if ~any(strcmp(curField, defaultFields))
        if ~any(strcmpi(curField, defaultFields))
            warning('parseModel:fieldName', ...
                'Input field "%s" not recognized', curField)
        else
            match = defaultFields{strcmpi(curField, defaultFields)};
            warning('parseModel:fieldNameCase', ...
                'Case insensitive match for input field "%s"', match)
            modelParams.(match) = modelParams.(curField);
        end
    end
end

for field = 1:numel(defaultFields)
    curField = defaultFields{field};
    if ~isfield(modelParams, curField)
        modelParams.(curField) = defaultParams.(curField);
    end
    assignin('caller', curField, modelParams.(curField))
end
