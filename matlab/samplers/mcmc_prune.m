% simple MCMC rejection loop with independent proposals
% Inputs:
%   y: proposal samples
%   lFex: [log(exact density), QoI] evaluated at y
%   lFapp: log(proposal density at y)
%
% Outputs:
%   y: pruned samples
%   lFex: exact data at new y
%   lFapp: log(proposal density at new y)
%   num_of_rejects: total number of rejections
%   rej_distribution: unnormalised distribution function of the number of
%       consequitive rejections, i.e. rej_distribution(L) ~ Prob(lag==L)
%
% See also: iw_prune, essinv, tt_irt_debias

function [y,lFex,lFapp,num_of_rejects, rej_distribution]=mcmc_prune(y,lFex,lFapp)

rej_distribution = zeros(1,1);

M = numel(lFapp);
num_of_rejects = 0;
rej_seq = 0;
for i=1:M-1
    alpha = lFex(i+1,1) - lFex(i,1) - lFapp(i+1) + lFapp(i);
    alpha = exp(alpha);
    if (alpha<rand)
        % reject, i.e. copy the i-th data to (i+1)
        y(i+1,:) = y(i,:);
        lFapp(i+1) = lFapp(i);
        lFex(i+1,:) = lFex(i,:);
        num_of_rejects = num_of_rejects+1;
        rej_seq = rej_seq + 1;
    elseif (rej_seq>0)
        % accept, save the previous rej_seq and reset it
        if (length(rej_distribution)>=rej_seq)
            rej_distribution(rej_seq,1) = rej_distribution(rej_seq,1)+1;
        else
            rej_distribution(rej_seq,1) = 1;
        end
        rej_seq = 0;
    end
end
fprintf('mcmc_prune completed with %d rejections (%g%%)\n', num_of_rejects, num_of_rejects/M*100);

end
