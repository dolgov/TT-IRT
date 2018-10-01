function [params] = parse_diffusion_inputs(varargin)
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end
if (~isfield(params, 'sigma'))
    params.sigma = input('Prior model variance sigma = ? (default 1): ');
    if (isempty(params.sigma))
        params.sigma = 1;
    end
end
if (~isfield(params, 'corr_length'))
    params.corr_length = input('Prior model corr_length = ? (default 1): ');
    if (isempty(params.corr_length))
        params.corr_length = 1;
    end
end
if (~isfield(params, 'nu'))
    params.nu = input('Decay rate nu = ? (default 2): ');
    if (isempty(params.nu))
        params.nu = 2;
    end
end

% Only uniform distribution here
check_lgwt; 

if (~isfield(params, 'runs'))
    params.runs = input('Number of test runs = ? (default 16): ');
    if (isempty(params.runs))
        params.runs = 16;
    end
end
if (~isfield(params, 'meshlevel'))
    params.meshlevel = input('Spatial meshlevel = ? ("1" <=> h=1/32) (default 2): ');
    if (isempty(params.meshlevel))
        params.meshlevel = 2;
    end
end

if (~isfield(params, 'sigma_n'))
    params.sigma_n = input('Noise variance sigma_n = ? (default 1e-2): ');
    if (isempty(params.sigma_n))
        params.sigma_n = 1e-2;
    end
end
if (~isfield(params, 'm0'))
    params.m0 = input('Number of measurements in each direction m0 = ? (default 3): ');
    if (isempty(params.m0))
        params.m0 = 3;
    end
end
if (~isfield(params, 'y0'))
    params.y0 = input('Truth value of parameters for synthetic data y0 = ? (default 1.5): ');
    if (isempty(params.y0))
        params.y0 = 1.5;
    end
end

if (~isfield(params, 'log2N'))
    params.log2N = input('log2-number of samples in the chain log2N = ? (default 14): ');
    if (isempty(params.log2N))
        params.log2N = 14;
    end
end

end


function check_lgwt()
if (exist('lgwt', 'file')==0)
    try
        fprintf('Gauss-Legendre rule is not found. Downloading...\n');
        urlwrite('https://uk.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/4540/versions/1/download/zip', 'lgwt.zip');
        unzip('lgwt.zip');
        rehash;
    catch ME
        error('%s. Please download lgwt.m (from e.g. MatlabCentral)', ME.message);
    end
end

end

