% Predator-prey parameter calibration example
% Using DIRT-MCMC sampling
% The following parameters can be passed in the form 
% 'param_name1', param_value1, 'param_name2', param_value2 and so on.
%   Model parameters:
%       sigma_n: observation noise variance (>0)
%       xtrue: synthetic true parameters (vector of 8 values)
%       obs_times: vector of time points where states are observed
%       data: P,Q observed data at obs_times (numel(obs_times) x 2 matrix)
%       domain: support interval of each variable
%   Test parameters:
%       Nsamples: length of MCMC chain to produce
%       runs: number of runs of the test (must be >1 for error estimation)
%   Approximation (DIRT) parameters:
%       n: number of grid points in each variable
%       R0: TT rank (all TT ranks are set to R0)
%       beta: a vector of tempering powers
%
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   num_of_rejects: numbers of rejections in MH method
%   tau: IACT
%   tau_ess: N/ESS
%   evalcnt: total number of function evaluations in DIRT
%   ttimes_approx: CPU time of DIRT construction
%   ttimes_sample: CPU time of DIRT sampling
%   Fdist: average Forstner-Moonen distance
% Each of these variables is a vector of params.runs values produced in the
% corresponding runs. The script computes expectations:
%   meanZ: mean variables (1 x d x params.runs tensor)
%   covZ:  covariances (d x d x params.runs tensor)
% If the standard sequential 'for' is used, additional variables created are
%   IRT: DIRT structure (contains all TT decompositions)
%   z: MCMC samples (params.Nsamples x d matrix)
%   lFapp: log(proposal density) values (after rejection)
%   lFex: log(exact density) values (after rejection)

function test_predator_prey_dirt(varargin)
% Check prerequisites
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    cd(mydir); cd('..'); cd('..'); cd('utils'); check_ttirt;
end
check_tt;
check_mcmc;
cd(mydir);

params = parse_pp_inputs(varargin{:});
% extra parameters for TT
if (~isfield(params, 'n'))
    params.n = input('Uniform grid size for TT n = ? (default 18): ');
    if (isempty(params.n))
        params.n = 18;
    end
end
if (~isfield(params, 'R0'))
    params.R0 = input('TT rank R0 = ? (default 13): ');
    if (isempty(params.R0))
        params.R0 = 13;
    end
end
if (~isfield(params, 'domain'))
    params.domain = input('Normalised domain of all variables = ? (default [0.6 1.6]): ');
    if (isempty(params.domain))
        params.domain = [0.6 1.6];
    end
end
if (~isfield(params, 'beta'))
    params.beta = input('Tempering powers beta = ? (default 10.^(-4:1/2:0)): ');
    if (isempty(params.beta))
        params.beta = 10.^(-4:1/2:0);
    end
end

% Initial grid for DIRT (Level 0)
xsf = arrayfun(@(i)linspace(params.domain(1),params.domain(2),params.n)', (1:numel(params.xtrue))', 'UniformOutput', 0);

ind = 8:-1:1; % which parameters we infer
d = numel(ind);

if (nargin<1)||(isempty(params.data))
    % Synthesize some data
    opts = odeset('RelTol', 1e-6);
    [~,params.data] = ode45(@(t,y)PP_RHS(t,y,params.xtrue), params.obs_times, params.xtrue(1:2), opts);
    params.data = params.data + sqrt(params.sigma_n)*randn(size(params.data));
end

for irun=1:params.runs
    tic;
    % Bad things may happen if you change the remaining parameters, so
    % those are not a part of the params structure
    IRT = tt_dirt_approx(xsf(ind), @(x,beta1,beta2)PP_loglikelihood(x,params.data,params.obs_times,params.sigma_n,params.xtrue,ind)*(beta2-beta1), ...
                         params.beta, 'IRTdenom', false, 'crossmethod', 'amen_cross_s', ...
                         'nswp', 1, 'kickrank', 0, 'y0', params.R0, 'interpolation', 's', ...
                         'boundary', true, 'reference', 'n4', 'plotdiag', false, 'testsamples', 100);
    ttimes_approx(irun) = toc;
    evalcnt(irun) = sum(IRT.evalcnt);
    
    tic;
    [z,lFapp,lFex] = tt_dirt_sample(IRT, randref(IRT.reference, params.Nsamples, d), @(x)PP_loglikelihood(x,params.data,params.obs_times,params.sigma_n,params.xtrue,ind));
    ttimes_sample(irun) = toc;
    
    plot(z(:,ind));
    title('DIRT-MCMC Chain');
    legend('P_0', 'Q_0', 'r', 'K', 'a', 's', 'u', 'v');
    
    % N/ESS
    tau_ess(irun) = essinv(lFex,lFapp);
    % number of rejections
    [z,lFex,lFapp,num_of_rejects(irun)] = mcmc_prune(z,lFex,lFapp);
    % autocorr time of MCMC
    [~,~,~,tau(irun),~,~] = UWerr(z,1.5,size(z,1),0);  tau(irun) = tau(irun)*2;
    
    % posterior moments
    meanZ(:,:,irun) = mean(z);
    covZ(:,:,irun) = squeeze(mean((z-mean(z)).*reshape(z-mean(z), [], 1, d)));        
end

% Forstner distance
covZmean = mean(covZ,3);
for irun=1:size(covZ,3)
    Fdist(irun) = sum(log(eig(covZ(:,:,irun), covZmean)).^2);
end

% Print statistical data
fprintf('%%Rejected:\t %g+-%g\n', mean(num_of_rejects)/params.Nsamples*100, std(num_of_rejects)/params.Nsamples*100)
fprintf('IACT:\t\t %g+-%g\n', mean(tau), std(tau))
fprintf('N/ESS:\t\t %g+-%g\n', mean(tau_ess), std(tau_ess))
fprintf('Eval/lvl:\t %g\n', mean(evalcnt)/numel(params.beta)-100);
fprintf('DIRT time:\t %g+-%g\n', mean(ttimes_approx), std(ttimes_approx));
fprintf('Sampling time:\t %g+-%g\n', mean(ttimes_sample), std(ttimes_sample));
fprintf('FM distance:\t %g+-%g\n', mean(Fdist), std(Fdist));


% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
