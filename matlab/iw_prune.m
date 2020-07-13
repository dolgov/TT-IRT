% Importance Weighting cleanup
% Inputs:
%   lFex: Samples of rows [log(exact density), QoI]
%   lFapp: Samples of log(proposal density)
%
% Outputs:
%   lFex_: corrected data, lFex_ = lFex.*(Fex/Fapp)
%
%   isstd: relative standard deviation of the ratio 
%          Fex/Fapp = (1/Z)exp(lFex-lFapp), where Z is normalising constant
%   max_ratio: maximum of the ratio Fex/Fapp
%   err1: empirical L1-error, err1 = <|Fex - Fapp|>

function [lFex,isstd,max_ratio,err1]=iw_prune(lFex,lFapp)

% Importance Sampling: exact QoI = QoI*L/Lapp
isstd = exp(lFex(:,1) - lFapp);
renorm = mean(isstd); % IS normalisation constant
isstd = isstd/renorm;

% Extra error checks
max_ratio = max(isstd);
renorm = log(renorm);
err1 = mean(abs(exp(lFex(:,1)-renorm)-exp(lFapp))./exp(lFapp));

lFex = lFex.*repmat(isstd, 1, size(lFex,2));
isstd = sqrt(mean((isstd - 1).^2));

end
