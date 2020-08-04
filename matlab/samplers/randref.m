% Samples uniform or truncated normal pseudorandom numbers
%   function [y] = RANDREF(reference, sizes)
% Inputs:
%   reference: reference density (['UNIform'], 'Normal' or 'Normal S',
%              where S is the number of sigmas defining the [-S,S] support
%              of the reference variables. S=4 if absent)
%   sizes: arbitrary combination of array sizes compatible with rand.
%          Alternatively, you can pass an array of numbers in [0,1] 
%          (e.g. QMC lattice points) that will be transformed to reference.
% Outputs:
%   y: array of samples
function [y] = randref(reference, varargin)

if (isscalar(varargin{1}))
    % varargin contains sizes, sample corresponding Uniform
    y = rand(varargin{:});
else
    % These are some samples on [0,1] -- could be QMC
    y = varargin{1};
end
    
if (lower(reference(1))~='u')
    % Truncated normal
    sigma = double(lower(reference));
    sigma = sigma((sigma==46) | ((sigma>=48) & (sigma<=57))); % extract numbers and dot
    sigma = str2double(char(sigma));
    if (isnan(sigma))
        sigma = 4;
    end
    % Multiply erf by this to get the truncated CDF between 0 and 1
    cdf_ifactor = erf(sigma/sqrt(2))/0.5;
    
    y = erfinv((y-0.5)*cdf_ifactor)*sqrt(2);
end
end
