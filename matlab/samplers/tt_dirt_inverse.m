function [q,lFapp] = tt_dirt_inverse(IRTstruct, x)
% Inverse of the Deep Inverse Rosenblatt Transform
% Inputs:
%   IRTstruct: a structure from tt_dirt_approx with density ratios in TT
%   x:  a M x d array of samples on the original space
%
% Outputs:
%   q: transformed samples (q is on the reference space)
%   lFapp: log(pushforward density function) sampled at x
%
% See also: tt_dirt_sample, tt_dirt_approx, tt_rt_sqr


if (strcmp(IRTstruct.crossmethod, 'build_ftt'))
    error('dirt_inverse is implemented only for spline interpolation of tt_tensors, but "ftt" detected.');
end
if (IRTstruct.interpolation(1)~='s')
    warning('dirt_inverse is implemented only for spline interpolation of tt_tensors, but "%s" detected. Inversion may be inaccurate.', IRTstruct.interpolation);
end

nlvl = numel(IRTstruct.beta)-1;
lFapp = 0;
q = x;

if (IRTstruct.reference(1)~='u')
    % Parse the domain size
    sigma = double(IRTstruct.reference);
    sigma = sigma((sigma==46) | ((sigma>=48) & (sigma<=57))); % extract numbers and dot
    sigma = str2double(char(sigma));
    if (isnan(sigma))
        sigma = 4;
    end
    % Multiply erf by this to get the truncated normal CDF between 0 and 1
    cdf_factor = 0.5/erf(sigma/sqrt(2));
end


% Sample in increasing order
% Level 0 first
[q,dlFapp] = tt_rt_sqr(IRTstruct.x0, IRTstruct.F0, q);
% R^{-1}, back to reference measure
if (IRTstruct.reference(1)~='u')
    q = erfinv((q-0.5)/cdf_factor)*sqrt(2);
end
lFapp = lFapp + dlFapp;

% Other levels from the front
for j=1:nlvl
    if (IRTstruct.reference(1)~='u')
        % Subtract log(reference density). Ignore additive const log(cdf_factor)
        lFapp = lFapp + sum(q.^2,2)/2;
    end
    [q,dlFapp] = tt_rt_sqr(IRTstruct.x, IRTstruct.F{j}, q);
    % R^{-1}, back to reference measure
    if (IRTstruct.reference(1)~='u')
        q = erfinv((q-0.5)/cdf_factor)*sqrt(2);
    end
    lFapp = lFapp + dlFapp; % Add log(Jacobian)
end
end
