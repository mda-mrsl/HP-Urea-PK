function y = vifFunction(t, ymax, a, tmax, t0)

% Gamma-variate function parameterized by ymax, a, tmax and t0
%##########################################################################
%# y = vifFunction(t, ymax, a, tmax, [t0])
%#
%#  Calculate a gamma-variate vascular input function in a way that 
%#  is more amenable to least-squares curve fitting. Unlike the default
%#  gamma-variate function, coupling between input parameters is 
%#  eliminated and it is more straightforward to estimate parameters
%#  from a given curve.
%#
%#  REQUIRED INPUTS:
%#      t     - Vector of times at which to evaluate function
%#      ymax  - Maximum amplitude of gamma-variate curve
%#      a     - Shape parameter of gamma-variate curve
%#      tmax  - Time at which maximum amplitude occurs
%#
%#  OPTIONAL INPUTS:
%#      t0    - Delay time to first point on gamma-variate curve
%#            default  = 0
%#
%#  OUTPUT:
%#      y     - Gamma variate curve
%#
%#  LITERATURE:
%#      Madsen MT. Phys Med Biol 1992.
%#
%# Keith Michel
%# kamichel at mdanderson.org
%# 04/2018
%##########################################################################

if nargin<4, help(mfilename); return; end
if nargin<5, t0 = 0; end
% The following may need to be handled more gracefully for robust fitting:
assert(tmax > t0, 'vifFunction:delayTime', ...
    'Time of maximum must be greater than onset delay!')
tp = (t - t0) ./ (tmax - t0);
tp(t-t0 < 0) = 0;
y = ymax * (tp.^a) .* exp( a*(1-tp));