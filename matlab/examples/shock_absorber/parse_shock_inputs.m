function [params] = parse_shock_inputs(varargin)
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end

if (~isfield(params, 'D'))
    params.D = input('Number of covariates D = ? (default 6): ');
    if (isempty(params.D))
        params.D = 6;
    end
end
if (~isfield(params, 'x'))
    params.x = input('Covariates data x = ? (default empty, means random generation): ');
end
if (~isfield(params, 'log2N'))
    params.log2N = input('log2(number of samples in the chain) log2N = ? (default 14): ');
    if (isempty(params.log2N))
        params.log2N = 14;
    end
end
if (~isfield(params, 'runs'))
    params.runs = input('Number of test runs = ? (default 8): ');
    if (isempty(params.runs))
        params.runs = 8;
    end
end

end