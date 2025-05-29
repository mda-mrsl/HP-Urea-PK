function parseAcq(acqParams)

% Parse structure containing acquisition parameters
%##########################################################################
%#  parseAcq(acqParams)
%#
%#  Interpret structure containing acquisition parameter values, replacing
%#  any missing parameters with defaults. No outputs should be assigned
%#  when using this function, variables are assigned by name into the
%#  workspace of the caller.
%#
%#  OPTIONAL INPUTS:
%#      acqParams - Structure containing acquisition parameters
%#        > If empty or not provided, defaults are used for all parameters.
%#        > The following fields may be provided:
%#          acqScheme - One of the following acquisition schemes:
%#            'conventional' - Fixed TR and flip angle spoiled imaging
%#            'periodicRF'   - Alternate flip angle at end of each segment
%#            'periodicDP'   - Diffusion preparation at end of each segment
%#            Default = 'conventional'
%#          tEnd      - Time of acquisition end                      [s]
%#            Default = 60
%#          TR        - Repetition time                              [s]
%#            Default = 1
%#            > Time of first excitation is zero.
%#            > Time of last excitation is TR * floor(tEnd/TR).
%#            > Can be a monotonically increasing vector for uneven spacing
%#                (passed directly as first input to perfModel)
%#          fa        - Flip angle                                   [deg]
%#            Default = 20
%#          errFa     - Error factor in flip angle calibration
%#            Default = 1
%#          tSegment  - Segment duration                             [s]
%#            Default = 10
%#          faSegment - Alternate flip angle at end of each segment  [deg]
%#            Default = 90
%#          bval      - Gradient b-value at end of each segment    [s/mm^2]
%#            Default = 30
%#          TE        - Echo time                                    [ms]
%#            Default = 0
%#
%#  OUTPUTS: (assigned into workspace of caller)
%#      t     - Vector of times at which excitations occur         [s]
%#      fa    - Vector of corresponding flip angles                [deg]
%#      bval  - Vector of b-values for each excitation             [s/mm^2]
%#      TE    - Acquisition echo time                              [ms]
%#
%#  LITERATURE:
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 04/2018
%##########################################################################

if ~nargin, acqParams = []; end
if isempty(acqParams)
    acqParams = struct();
else
    assert(isa(acqParams, 'struct'), 'parseAcq:input', ...
        'Input must be a structure')
end

defaultParams = struct( ...
    'acqScheme',    'conventional', ...
    'tEnd',         60,	... % [s]
    'TR',           1,	... % [s]
    'fa',           20,	... % [deg]
    'errFa',        1,	...
    'tSegment',     10,	... % [s]
    'faSegment',    90,	... % [deg]
    'bval',         30,	... % [s/mm^2]
    'TE',           0);     % [ms]

defaultFields = fieldnames(defaultParams);
inputFields   = fieldnames(acqParams);
for field = 1:numel(inputFields)
    curField = inputFields{field};
    if ~any(strcmp(curField, defaultFields))
        if ~any(strcmpi(curField, defaultFields))
            warning('parseAcq:fieldName', ...
                'Input field "%s" not recognized', curField)
        else
            match = defaultFields{strcmpi(curField, defaultFields)};
            warning('parseAcq:fieldNameCase', ...
                'Case insensitive match for input field "%s"', match)
            acqParams.(match) = acqParams.(curField);
        end
    end
end

for field = 1:numel(defaultFields)
    curField = defaultFields{field};
    if ~isfield(acqParams, curField)
        acqParams.(curField) = defaultParams.(curField);
    end
end

if isscalar(acqParams.TR)
    t = 0 : acqParams.TR : acqParams.tEnd;
else
    assert(all(diff(acqParams.TR) > 0), 'parseAcq:TRvector', ...
        'TR vector must be monotonically increasing')
    t = acqParams.TR;
end
fa   = acqParams.fa * ones(size(t)) * acqParams.errFa;
bval = zeros(size(t));
TE   = acqParams.TE;
switch lower(acqParams.acqScheme)
    case 'conventional'
    case 'periodicrf'
        for segment = 1:floor(t(end)/acqParams.tSegment)
            ts = find(t >= segment*acqParams.tSegment, 1);
            fa(ts) = acqParams.faSegment * acqParams.errFa;
        end
    case 'periodicdp'
        for segment = 1:floor(t(end)/acqParams.tSegment)
            ts = find(t >= segment*acqParams.tSegment, 1);
            bval(ts) = acqParams.bval;
        end
    otherwise
        warning('parseAcq:acqScheme', ['Acquisition scheme not ', ...
            'recognized. Defaulting to "conventional".'])
end

assignParams = {'t', 'fa', 'bval', 'TE'};
for param = 1:numel(assignParams)
    assignin('caller', assignParams{param}, eval(assignParams{param}))
end
