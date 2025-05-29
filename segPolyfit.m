function coeffs = segPolyfit(t, fcn, order, tSub)

% Fit Nth order polynomial to segments of a single-variable function
%##########################################################################
%# coeffs = segPolyfit(t, vif, order, tSub)
%#
%#  Fit Nth order polynomials to segments of a single-variable function.
%#
%#  REQUIRED INPUTS:
%#      t     - Vector of points between which to fit function
%#      fcn   - Handle to single-variable function to be fit
%#      order - Integer order of fitting polynomials
%#
%#  OPTIONAL INPUTS:
%#      tSub  - Integer by which to subdivide segments for fitting
%#            default  = 40
%#
%#  OUTPUT:
%#      coeffs - Matrix of size [numel(t)-1, 5] containing 4th order 
%#               polynomial coefficients fit to function segments
%#
%#  LITERATURE:
%#      
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 04/2018
%##########################################################################

if nargin<3, help(mfilename); return; end
if nargin<4, tSub = 40; end
assert(isa(fcn, 'function_handle'), 'vifPolyfit:vifHandle', ...
    'fcn must be a function handle')
nSeg = numel(t)-1;
coeffs = zeros(nSeg, order+1);
for ii=1:nSeg
    tt = linspace(t(ii), t(ii+1), tSub);
    coeffs(ii,:) = polyfit(tt-tt(1), fcn(tt), order);
end