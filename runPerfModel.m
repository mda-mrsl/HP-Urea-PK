function [Mx, Mz, Mex, Mez] = runPerfModel(acqParams, modelParams)

% Evaluate two-site model for hyperpolarized perfusion agent
%##########################################################################
%# [Mx, Mz, Mex, Mez] = runPerfModel(acqParams, modelParams)
%#
%#  Solve for transverse and longitudinal magnetizations of a 
%#  hyperpolarized perfusion agent using a two-site model with 
%#  bidirectional exchange. Unlike perfModel() this function accepts
%#  structures of acquisition and model parameters, which allows you to 
%#  keep the caller workspace clean.
%#
%#  OPTIONAL INPUTS:
%#      acqParams - Structure containing acquisition parameters
%#        > Acquisition parameters not given in acqParams are set to
%#          default values (see parseAcq).
%#      modelParams - Structure containing model parameters
%#        > Model parameters not given in modelParams are set to default
%#          values (see parseModel).
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

%% Parse inputs
if ~nargin,   acqParams   = []; end
parseAcq(acqParams)
if nargin<2,  modelParams = []; end
parseModel(modelParams)

%% Evaluate Model
[Mx, Mz, Mex, Mez] = perfModel(t, fa, vb, vei, kve, T1, vif, ...
                               M0, De, DvFactor, bval, T2, TE);
