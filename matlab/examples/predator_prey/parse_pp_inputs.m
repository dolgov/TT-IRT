function [params] = parse_pp_inputs(varargin)
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end

if (~isfield(params, 'sigma_n'))
    params.sigma_n = input('Noise variance sigma_n = ? (default 2): ');
    if (isempty(params.sigma_n))
        params.sigma_n = 2;
    end
end
if (~isfield(params, 'xtrue'))
    params.xtrue = input('True parameters xtrue = ? (default [50, 5, 0.6, 100, 1.2, 25, 0.5, 0.3]): ');
    if (isempty(params.xtrue))
        params.xtrue = [50, 5, 0.6, 100, 1.2, 25, 0.5, 0.3];
    end
end
if (~isfield(params, 'obs_times'))
    params.obs_times = input('Observed time points obs_times = ? (default linspace(0,50,13)): ');
    if (isempty(params.obs_times))
        params.obs_times = linspace(0,50,13);
    end
end

if (~isfield(params, 'data'))
    params.data = input('Observation data = ? (default empty, means random generation): ');
end

if (~isfield(params, 'Nsamples'))
    params.Nsamples = input('Number of samples in the chain Nsamples = ? (default 1e4): ');
    if (isempty(params.Nsamples))
        params.Nsamples = 1e4;
    end
end
if (~isfield(params, 'runs'))
    params.runs = input('Number of test runs = ? (default 2): ');
    if (isempty(params.runs))
        params.runs = 2;
    end
end

end