% Predator-prey parameter calibration example
% Using Stein Variational Newton
% The following parameters can be passed in the form 
% 'param_name1', param_value1, 'param_name2', param_value2 and so on.
%   Model parameters:
%       sigma_n: observation noise variance (>0)
%       xtrue: synthetic true parameters (vector of 8 values)
%       obs_times: vector of time points where states are observed
%       data: P,Q observed data at obs_times (numel(obs_times) x 2 matrix)
%   Test parameters:
%       Nsamples: number of samples to produce
%       runs: number of runs of the test (must be >1 for error estimation)
%   Approximation (SVN) parameters:
%       stepsize: damped Newton step size
%       itermax: maximum number of Newton iterations
%       initial_std: standard deviation of the initial cloud of points
%
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   ttimes_sample: CPU time of SVN sampling
%   Fdist: average Forstner-Moonen distance
% Each of these variables is a vector of params.runs values produced in the
% corresponding runs. The script computes expectations:
%   meanZ: mean variables (1 x d x params.runs tensor)
%   covZ:  covariances (d x d x params.runs tensor)
% If the standard sequential 'for' is used, additional variables created are
%   z: transformed samples (params.Nsamples x d matrix)

function test_predator_prey_svn(varargin)
% Check prerequisites
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    cd(mydir); cd('..'); cd('..'); cd('utils'); check_ttirt;
end
cd(mydir); check_svn; cd(mydir);

params = parse_pp_inputs(varargin{:});

% SVN parameters
if (~isfield(params, 'stepsize'))
    params.stepsize = input('Newton step size = ? (default 2e-2): ');
    if (isempty(params.stepsize))
        params.stepsize = 2e-2;
    end
end
if (~isfield(params, 'itermax'))
    params.itermax = input('Max number of iterations = ? (default 23): ');
    if (isempty(params.itermax))
        params.itermax = 23;
    end
end
if (~isfield(params, 'initial_std'))
    params.initial_std = input('Standard deviation of initial samples = ? (default 2e-2): ');
    if (isempty(params.initial_std))
        params.initial_std = 2e-2;
    end
end

ind = 8:-1:1; % which parameters we infer
d = numel(ind);

if (nargin<1)||(isempty(params.data))
    % Synthesize some data
    opts = odeset('RelTol', 1e-6);
    [~,params.data] = ode45(@(t,y)PP_RHS(t,y,params.xtrue), params.obs_times, params.xtrue(1:2), opts);
    params.data = params.data + sqrt(params.sigma_n)*randn(size(params.data));
end

for irun=1:params.runs    
    % Model structure for SVN
    model = struct();
    model.n = d;
    model.data = params.data;
    model.obs_times = params.obs_times;
    model.sigma_n = params.sigma_n;
    model.xtrue = params.xtrue;
    model.ind = ind;
    
    obs = struct();
    obs.std2 = params.sigma_n;
    
    prior = struct();
    prior.C0i = zeros(d,d);
    
    tttic = tic;
    [z, stepsize, timeave] = SVN_H(1+params.initial_std*randn(d, params.Nsamples), params.stepsize, params.itermax, model, prior, obs);
    ttimes_sample(irun) = toc(tttic);
    z = z.';   
    
    plot(z(:,ind));
    title('SVN samples');
    legend('P_0', 'Q_0', 'r', 'K', 'a', 's', 'u', 'v');    
    
    meanZ(:,:,irun) = mean(z);
    covZ(:,:,irun) = squeeze(mean((z-mean(z)).*reshape(z-mean(z), [], 1, d)));        
end

% Forstner distance
covZmean = mean(covZ,3);
for irun=1:size(covZ,3)
    Fdist(irun) = sum(log(eig(covZ(:,:,irun), covZmean)).^2);
end

% Print statistical data
fprintf('SVN time:\t %g+-%g\n', mean(ttimes_sample), std(ttimes_sample));
fprintf('FM distance:\t %g+-%g\n', mean(Fdist), std(Fdist));


% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
