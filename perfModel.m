function [Mx, Mz, Mex, Mez] = perfModel(t, fa, vb, vei, kve, T1, vif, M0, ...
                              De, DvFactor, bval, T2, TE)

% Evaluate two-site model for hyperpolarized perfusion agent
%##########################################################################
%# [Mx, Mz, Mex, Mez] = perfModel(t, fa, vei, vb, kve, T1, vif, M0, ...
%#                      De, DvFactor, bval, T2, TE)
%#
%#  Solve for transverse and longitudinal magnetizations of a
%#  hyperpolarized perfusion agent using a two-site model with
%#  bidirectional exchange
%#
%#  REQUIRED INPUTS:
%#      t     - Vector of times at which excitations occur            [s]
%#      fa    - Vector of corresponding flip angles                   [deg]
%#      vb    - Vascular tissue volume fraction
%#      vei   - Distribution volume fraction of extravascular space
%#      kve   - Exchange rate between vascular and extrascular spaces [1/s]
%#      T1    - Time constant of losses other than RF and exchange    [s]
%#      vif   - One of the following vascular input functions:
%#        > Handle to function evaluating Mz in blood over time. Integrals
%#          are calculated using the MATLAB function integral(). This is
%#          the most numerically accurate method but also the slowest.
%#        > Structure containing gamma-variate curve parameters passed to
%#          vifFunction (see parseModel). A function handle is generated
%#          with these parameters and the VIF is integrated as above.
%#        > Matrix of size [numel(t)-1, N], where 2<=N<=5. Rows contain
%#          coefficients for the (N-1)th order polynomial fit to the
%#          segments of the VIF in the corresponding TR intervals.
%#          Integrals are calculated using the analytical solution for a
%#          4th order polynomial VIF. This method provides the best balance
%#          of speed and accuracy.
%#        > Vector of same length as t providing Mz in blood at each
%#          excitation time. Integrals are calculated using a linear
%#          approximation of the VIF in each TR interval. This method is
%#          the least numerically accurate.
%#
%#  OPTIONAL INPUTS:
%#      M0        - Initial Mz in extravascular space at first excitation
%#        > If M0 is not provided, a default value of 1e-6 is used
%#      De        - Extravascular apparent diffusion coefficient  [mm^2/s]
%#      DvFactor  - Factor applied to De to obtain the pseudodiffusion
%#                  coefficient characterizing microvascular flow
%#      bval      - Vector of b-values for each excitation        [s/mm^2]
%#        > If any of De, DvFactor and bval are not provided, no diffusion
%#          effects are simulated
%#      T2        - Spin-spin relaxation time constant            [ms]
%#      TE        - Acquisition echo time                         [ms]
%#        > If any of T2 and TE are not provided, no T2 effects are
%#          simulated
%#
%#  OUTPUTS:
%#      Mx  - Vector of transverse magnetizations
%#      Mz  - Vector of longitudinal magnetizations
%#      Mex - Vector of transverse extravascular magnetizations
%#      Mez - Vector of longitudinal extravascular magnetizations
%#
%#  LITERATURE:
%#      Bankson JA, et al. Cancer Res 2015.
%#        DOI 10.1158/0008-5472.CAN-15-0171
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 04/2018
%##########################################################################

if nargin<7, help(mfilename); error('perfModel:inputs', ...
                                    'Not enough inputs'); end

%% Parse inputs
nPts = numel(t);
assert(nPts==numel(fa), 'perfModel:faSize', ...
    'Number of elements in t and fa must match')
ve = vei * (1-vb);
A  = 1/T1 + kve/ve;
if isstruct(vif)
    vif = @(t) vifFunction(t, vif.vifMax, vif.vifShape, ...
                              vif.vifTMax, vif.vifDelay);
end
if isa(vif, 'function_handle')
    % Function below calculates Mz in extravascular space
    calMez = @(ii, Mez, t) ...
        ... % Extravascular Mz from previous TR interval w/ non-RF losses:
        Mez(ii-1) * exp(-A*(t(ii)-t(ii-1))) + ...
        ... % Mz extravasating from blood within current TR interval:
        kve * integral( @(tt) exp(-A*(tt-t(ii-1))).*vif(tt), t(ii-1), t(ii) );
    vifVals = vif(t);
elseif isvector(vif)
    assert(nPts==numel(vif), 'perfModel:vifVector', ['Number of ', ...
        'elements in vif vector must match t and fa'])
    vifVals = vif;
    vif     = zeros(nPts-1, 5);
    for ii=1:nPts-1
        vif(ii,:) = [0, 0, 0, ...
            polyfit([0, t(ii+1)-t(ii)], vifVals(ii:ii+1), 1)];
    end
else
    [nCoeffs, polyOrder] = size(vif);
    assert(nCoeffs==nPts-1, 'perfModel:vifCoeffsRows', ['Number of ', ...
        'rows in vifCoeffs must be one fewer than the number of ', ...
        'elements in t and fa'])
    assert(polyOrder<=5, 'perfModel:vifCoeffsCols', ['Number of ', ...
        'colums in vifCoeffs must be <=5 (max 4th order polynomial)'])
    if polyOrder<5
        vif(end,5) = 0;
        vif        = circshift(vif, 5-polyOrder, 2);
    end
    vifVals    = zeros(nPts,1);
    vifVals(1) = polyval(vif(1,:), t(1));
    for ii=2:nPts
        vifVals(ii) = polyval(vif(ii-1,:), t(ii)-t(ii-1));
    end
end
if nargin<8,     M0 = []; end
if isempty(M0),  M0 = 1e-6; end
if nargin<9,     De       = []; end
if nargin<10,    DvFactor = []; end
if nargin<11,    bval     = []; end
if isempty(De) || isempty(DvFactor) || isempty(bval)
    De       = 0;
    DvFactor = 0;
    bval     = zeros(nPts, 1);
end
assert(nPts==numel(bval), 'perfModel:bvalSize', ...
    'Number of elements in bval and t must match')
if nargin<12,   T2 = []; end
if nargin<13,   TE = []; end
if isempty(T2) || isempty(TE)
    T2 = inf;
    TE = 0;
end

%% Evaluate Model
Dv               = De * DvFactor;
attenDv          = exp(-bval*Dv);
attenDv(bval==0) = 1;               % Removes NaNs when DvFactor = Inf
attenDe          = exp(-bval*De);
attenT2          = exp(-TE/T2);

Mx     = zeros(nPts, 1);
Mex    = zeros(nPts, 1);
Mex(1) = M0 * attenDe(1) * sind(fa(1)) * attenT2;
Mx(1)  = Mex(1) + ...
         vb*vifVals(1) * attenDv(1) * sind(fa(1)) * attenT2;
Mz     = zeros(nPts, 1);
Mez    = zeros(nPts, 1);
Mez(1) = M0 * cosd(fa(1));
Mz(1)  = Mez(1) + vb*vifVals(1) * cosd(fa(1));

if isa(vif, 'function_handle')
    for ii=2:nPts
        Mez(ii) = calMez(ii, Mez, t);      % Extravascular Mz
        Mex(ii) = Mez(ii) * attenDe(ii) * sind(fa(ii)) * attenT2;
        Mez(ii) = Mez(ii) * cosd(fa(ii));
        Mx(ii)  = Mex(ii) + ...
                  vb*vifVals(ii) * attenDv(ii) * sind(fa(ii)) * attenT2;
        Mz(ii)  = Mez(ii) + vb*vifVals(ii) * cosd(fa(ii));
    end
else
    tSub = 40;
    for ii=2:nPts
        ta     = linspace(0, t(ii)-t(ii-1), tSub);
        tb     = ta.*ta;
        tc     = tb.*ta;
        td     = tc.*ta;
        c      = vif(ii-1,:);
        kernel = exp(A.*ta).*(((((               24*c(1)     )/A - ...
            (                         6*c(2)    +24*c(1).*ta))/A + ...
            (              2*c(3)    +6*c(2).*ta+12*c(1).*tb))/A - ...
            (     c(4)    +2*c(3).*ta+3*c(2).*tb+ 4*c(1).*tc))/A + ...
            (c(5)+c(4).*ta+  c(3).*tb+  c(2).*tc+   c(1).*td))/A - ...
            ((((24*c(1)/A - 6*c(2))/A + 2*c(3))/A - c(4))/A + c(5))/A;
        % Approximation of extravasation integral at tSub points in TR interval:
        intSub  = ( Mez(ii-1) + kve.*kernel ) .* exp(-A*ta);
        Mez(ii) = intSub(end);
        Mex(ii) = Mez(ii) * attenDe(ii) * sind(fa(ii)) * attenT2;
        Mez(ii) = Mez(ii) * cosd(fa(ii));
        Mx(ii)  = Mex(ii) + ...
                  vb*vifVals(ii) * attenDv(ii) * sind(fa(ii)) * attenT2;
        Mz(ii)  = Mez(ii) + vb*vifVals(ii) * cosd(fa(ii));
    end
end

