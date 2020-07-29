% A script to test DIRT on the Lorenz-40 example
% The following parameters can be passed in the form 
% 'param_name1', param_value1, 'param_name2', param_value2 and so on.
%   Model parameters:
%       d: problem dimension
%       sigma_n: standard deviation of the observation noise (must be >0)
%       x0true: synthetic true initial state (vector of d values)
%       sigma_truth: standard deviation of the initial state perturbation
%   Test parameters:
%       Nsamples: length of MCMC chain to produce
%   Approximation (DIRT) parameters:
%       n: number of grid points in each variable
%       a: size of the domain [-a,a]^d in the original space (must be >0)
%       R0: TT rank (all TT ranks are set to R0)
%       beta: a vector of tempering powers (must increase and end at 1)
%
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   tau_iact: IACT
%   tau_ess: N/ESS
%   tracecov: Trace of posterior covariance
%   IRT: DIRT structure (contains all TT decompositions)
%   z: DIRT samples (params.Nsamples x d matrix) (before rejection)
%   lFapp: log(proposal density) values (before rejection)
%   lFex: log(exact density) values (before rejection)
%
function test_lorenz(varargin)
% Parse model parameters
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end

if (~isfield(params, 'd'))
    params.d = input('Dimension d = ? (default 10): ');
    if (isempty(params.d))
        params.d = 10;
    end
end
if (~isfield(params, 'sigma_n'))
    params.sigma_n = input('Standard deviation of the observation noise sigma_n = ? (default 0.1): ');
    if (isempty(params.sigma_n))
        params.sigma_n = 0.1; % std
    end
end
if (~isfield(params, 'x0true'))
    params.x0true = input('Baseline initial state x0true = ? (default ones(1,d)): ');
    if (isempty(params.x0true))
        params.x0true = ones(1,params.d);
    end
end
if (~isfield(params, 'sigma_truth'))
    params.sigma_truth = input('Standard deviation of the initial perturbation sigma_truth = ? (default 0.01): ');
    if (isempty(params.sigma_truth))
        params.sigma_truth = 1e-2; % std
    end
end

% Log(prior) function
lpriorfun = @(x)-0.5*sum((x-params.x0true).^2, 2);

% Synth some data
params.x0true = params.x0true + params.sigma_truth*randn(1,params.d);
[t,xd] = ode45(@(t,x)lorenz_rhs(t,x,params.d), [0 0.1], params.x0true);
data = xd(end, 2:2:end) + params.sigma_n*randn(1,params.d/2);


% Grid for DIRT Level 0
if (~isfield(params, 'n'))
    params.n = input('Number of tensor grid points n = ? (default 17): ');
    if (isempty(params.n))
        params.n = 17;
    end
end
if (~isfield(params, 'a'))
    params.a = input('Level 0 domain size a = ? (default 10): ');
    if (isempty(params.a))
        params.a = 10;
    end
end
x = linspace(-params.a,params.a,params.n)';
xsf = repmat({x},params.d,1);

% TT rank
if (~isfield(params, 'R0'))
    params.R0 = input('TT rank R0 = ? (default 15): ');
    if (isempty(params.R0))
        params.R0 = 10;
    end
end

% Tempering powers
if (~isfield(params, 'beta'))
    params.beta = input('Tempering powers beta = ? (default 10.^(-4:1/2:0)): ');
    if (isempty(params.beta))
        params.beta = 10.^(-4:1/2:0);
    end
end

% Number of samples
if (~isfield(params, 'Nsamples'))
    params.Nsamples = input('Number of MCMC samples Nsamples = ? (default 1e4): ');
    if (isempty(params.Nsamples))
        params.Nsamples = 1e4;
    end
end

% Build DIRT
IRT = tt_dirt_approx(xsf, @(x,b1,b2)(lorenz_ll(x,data,params.sigma_n)*(b2-b1) + lpriorfun(x)*(b2^0.25-b1^0.25)), ...
                     params.beta, 'nswp', 1, 'kickrank', 0, 'y0', params.R0, 'interpolation', 's', ...
                     'boundary', true, 'stoptol', 0.1, 'reference', 'n3', 'IRTdenom', false);

% Sample MC points
[z,lFapp,lFex] = tt_dirt_sample(IRT, randref(IRT.reference, params.Nsamples,params.d), @(x)(lorenz_ll(x,data,params.sigma_n) + lpriorfun(x)));

% MCMC rejection and its properties
z2 = mcmc_prune(z,lFex,lFapp);
[~,~,~,tau_iact,~,~] = UWerr(z2,1.5,size(z2,1),0);  tau_iact = tau_iact*2;
tau_ess = essinv(lFex,lFapp);

% Trace of posterior covariance
w = exp(lFex-lFapp);
w = w'/sum(w);
tracecov = sum(w*(z-w*z).^2);

fprintf('\nIACT: %g\n', tau_iact);
fprintf('N/ESS: %g\n', tau_ess);
fprintf('Trace(Cov): %g\n', tracecov);

% Inferred initial state
hold off
ha = area((1:params.d)', (w*z)'*[1 0]+sqrt(w*(z-w*z).^2)'*[-2 4],  min(w*z-2*sqrt(w*(z-w*z).^2)), 'EdgeColor', 'none' );
set(ha(1), 'FaceColor', 'none');
set(ha(2), 'FaceColor', [0.9,0.9,1]);
hold on
plot(1:params.d, w*z, 1:params.d, params.x0true);
legend(' ', ' ', 'Mean', 'Truth');
hold off
title('Inferred initial state, mean +- 2 std. dev.');


% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end


