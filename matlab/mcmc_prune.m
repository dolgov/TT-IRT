% simple MCMC rejection loop with independent proposals
% Inputs:
%   y: proposal samples
%   Fex: [exact density, QoI] evaluated at y
%   Fapp: proposal density at y
% Outputs:
%   y: pruned samples
%   Fex: exact data at new y
%   Fapp: proposal density at new y
%   num_of_rejects

function [y,Fex,Fapp,num_of_rejects]=mcmc_prune(y,Fex,Fapp)

M = numel(Fapp);
num_of_rejects = 0;
for i=1:M-1
    alpha = log(Fex(i+1,1)) - log(Fex(i,1)) - log(abs(Fapp(i+1))) + log(abs(Fapp(i)));
    alpha = exp(alpha);
    if (alpha<rand)
        % reject, i.e. copy the i-th data to (i+1)
        y(i+1,:) = y(i,:);
        Fapp(i+1) = Fapp(i);
        Fex(i+1,:) = Fex(i,:);
        num_of_rejects = num_of_rejects+1;
    end
end
fprintf('mcmc_prune completed with %d rejections (%g%%)\n', num_of_rejects, num_of_rejects/M*100);

end
