% Importance Weighting cleanup
% Inputs:
%   Fex: [exact density, QoI] samples
%   Fapp: proposal density samples
% Outputs:
%   Fex: corrected data
%   isstd: relative standard deviation of the ratio Lex/Lapp

function [Fex,isstd]=iw_prune(Fex,Fapp)

% Importance Sampling: exact QoI = QoI*L/Lapp
isstd = Fex(:,1)./Fapp;
renorm = mean(isstd); % IS normalisation constant
isstd = isstd/renorm;
Fex = Fex.*repmat(isstd, 1, size(Fex,2));
isstd = sqrt(mean((isstd - 1).^2));

end
