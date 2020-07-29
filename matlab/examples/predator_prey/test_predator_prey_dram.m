% Predator-prey parameter calibration example
% Using DRAM sampling
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
%
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   num_of_rejects: numbers of rejections in MH method
%   tau: IACT
%   ttimes_sample: CPU time of DRAM
%   Fdist: average Forstner-Moonen distance
% Each of these variables is a vector of params.runs values produced in the
% corresponding runs. The script computes expectations:
%   meanZ: mean variables (1 x d x params.runs tensor)
%   covZ:  covariances (d x d x params.runs tensor)
% If the standard sequential 'for' is used, additional variables created are
%   results: DRAM results structure
%   z: DRAM samples (params.Nsamples x d matrix)

function test_predator_prey_dram(varargin)
% Check prerequisites
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    cd(mydir); cd('..'); cd('..'); cd('utils'); check_ttirt;
end
check_mcmc;
cd(mydir);

params = parse_pp_inputs(varargin{:});
if (~isfield(params, 'domain'))
    params.domain = input('Normalised domain of all variables = ? (default [0.6 1.6]): ');
    if (isempty(params.domain))
        params.domain = [0.6 1.6];
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
    % DRAM parameters
    options = struct();
    options.nsimu    = params.Nsamples; % number of samples
    options.adaptint = 10;
    options.drscale  = 2;
    options.adascale = 2.4/sqrt(d); % scale for adaptation
    options.qcov     = eye(d)*5; % initial covariance
    
    % create input arguments for the dramrun function
    model = struct();
    model.ssfun    = @(x,ddd)-2*PP_loglikelihood(x,params.data,params.obs_times,params.sigma_n,params.xtrue,ind);
    domain = struct();
    domain.par0    = ones(1,d); % initial value
    domain.bounds = repmat([params.domain(1); params.domain(2)], 1, d);
    
    tic;
    [results,z] = dramrun(model,struct(),domain,options);
    ttimes_sample(irun) = toc;
    
    num_of_rejects(irun) = (1-results.accepted)*size(z,1);
    % autocorr time of MCMC
    [~,~,~,tau(irun),~,~] = UWerr(z,1.5,size(z,1),0);  tau(irun) = tau(irun)*2;
    
    plot(z(:,ind));
    title('DRAM Chain');
    legend('P_0', 'Q_0', 'r', 'K', 'a', 's', 'u', 'v');    
    
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
fprintf('DRAM time:\t %g+-%g\n', mean(ttimes_sample), std(ttimes_sample));
fprintf('FM distance:\t %g+-%g\n', mean(Fdist), std(Fdist));            

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
