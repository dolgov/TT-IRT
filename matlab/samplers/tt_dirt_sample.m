function [z,lFapp,lFex] = tt_dirt_sample(IRTstruct, q, logpostfun, vec)
% Sampling of a density represented by a Deep Inverse Rosenblatt Transform
% Inputs:
%   IRTstruct: a structure from tt_dirt_approx with density ratios in TT
%   q:  a M x d array of seed points on [0,1]^d for uniform reference, or
%                                    on [-S,S]^d for trunc normal reference
%   logpostfun: a @(x)function of exact target log-density (optional)
%   vec: whether logpostfun can process vectorised x (optional, default true)
%
% Outputs:
%   z: transformed samples from q
%   lFapp: log(pushforward density)  (inv. Jacobian of z)
%   lFex: samples of log(exact density) (if logpostfun was given)
%
% See also: tt_dirt_approx, randref, tt_irt_sqr, tt_irt_fourier

nlvl = numel(IRTstruct.beta)-1;
lFapp = 0;
z = q;

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

% Sample in reversed order
for j=nlvl:-1:1
    if (IRTstruct.reference(1)~='u')
        z = erf(z/sqrt(2))*cdf_factor + 0.5; % Transform TNormal -> uniform
    end
    % Transform the current level
    if (strcmp(IRTstruct.crossmethod, 'build_ftt'))
        [z,Fapp] = eval_irt(IRTstruct.F{j}, z');
        z = z';
        Fapp = Fapp';
        dlFapp = log(Fapp);
    else
        if (IRTstruct.interpolation(1)=='s')
            [z,dlFapp] = tt_irt_sqr(IRTstruct.x, IRTstruct.F{j}, z);
        else
            [z,dlFapp] = tt_irt_fourier(IRTstruct.x, IRTstruct.F{j}, z);
        end
    end
    lFapp = lFapp + dlFapp; % Add log(Jacobian)
    if (IRTstruct.reference(1)~='u')
        % Add log(reference density). Ignore additive const log(cdf_factor)
        lFapp = lFapp + sum(z.^2,2)/2;
    end
end

if (IRTstruct.reference(1)~='u')
    z = erf(z/sqrt(2))*cdf_factor + 0.5;
end
% Sample Level 0
if (strcmp(IRTstruct.crossmethod, 'build_ftt'))
    [z,Fapp] = eval_irt(IRTstruct.F0, z');
    z = z';
    Fapp = Fapp';
    dlFapp = log(Fapp);
else
    % Unlikely the original tempered density is band-limited, use splines
    % regardless of IRTstruct.interpolation
    [z,dlFapp] = tt_irt_sqr(IRTstruct.x0, IRTstruct.F0, z);
end
lFapp = lFapp + dlFapp;

% Exact density
if (nargin>2)&&(~isempty(logpostfun))
    if (nargin<4)||(isempty(vec))
        vec = true;
    end
    lFex = logpostfun_vec(z, logpostfun, vec);
end

end



% A wrapper for vectorised (or not) user fun
function [y]=logpostfun_vec(x, logpostfun, vec)
if (vec)
    y = logpostfun(x);
else
    y = zeros(size(x,1), 1);
    if isempty(gcp('nocreate'))
        for i=1:size(x,1)
            y(i) = logpostfun(x(i,:));
        end
    else % Use parallel pool if possible
        parfor i=1:size(x,1)
            y(i) = logpostfun(x(i,:));
        end
    end
end
end

