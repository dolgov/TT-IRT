% normalised Inverse Effective Sample Size, an analog of IACT for
% Importance Weighting. Also an estimate to 1 + chi^2divergence.
% Inputs:
%   lFex: log-exact density
%   lFapp: log-sampling density
% Outputs:
%   tau: N/ESS
%
% See also: iw_prune, mcmc_prune, tt_irt_debias, UWerr

function [tau] = essinv(lFex,lFapp)
dF = lFex - lFapp;
dF = dF - max(dF);
tau = numel(lFapp)*sum(exp(dF*2))/sum(exp(dF))^2;
end
