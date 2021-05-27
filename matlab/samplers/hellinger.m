% Hellinger distance approximation from samples
% 2*H^2 = E_{Fapp}[ \sqrt{(Fex/Zex)/Fapp} - 1 ]^2
% Inputs:
%   lFex: log(exact density Fex), can be unnormalised
%   lFapp: log(sampling density Fapp), must be normalised
% Outputs:
%   H
%
% See also: iw_prune, mcmc_prune, tt_irt_debias, UWerr, essinv

function [H] = hellinger(lFex,lFapp)
dF = lFex - lFapp;
dF = dF - max(dF);
lZex = log(mean(exp(dF))); % up to  +max(dF) which cancels below anyway
H = mean( (exp(0.5*(dF-lZex)) - 1).^2 );
H = sqrt(H/2);
end
