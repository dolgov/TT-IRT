% Shock absorber tests with DRAM
% The following parameters can be passed in the form 
% 'param_name1', param_value1, 'param_name2', param_value2 and so on:
%   x: a matrix of covariates (can be empty, in this case drawn from normal
%      distribution)
%   log2N: log2(number of samples in the chain)
%   D: number of covariates
%   runs: number of runs of the test (must be >1 for error estimation)
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   ttimes_dram: CPU times (for each run)
%   Q_dram: quantities of interest
%   N_dram: total number of pi evaluations
%   tauint_dram = zeros(nruns, 1): IACT
function test_shock_absorber_dram(varargin)
% Check prerequisites
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    cd(mydir); cd('..'); cd('..'); cd('utils'); check_ttirt;
end
check_mcmc;
cd(mydir);

params = parse_shock_inputs(varargin{:});

% Data table from the paper
y = [6700 6950 7820 8790 9120 9660 9820 11310 11690 11850 11880 12140 12200 ...
     12870 13150 13330 13470 14040 14300 17520 17540 17890 18420 18960  ...
     18980 19410 20100 20100 20150 20320 20900 22700 23490 26510 27410  ...
     27490 27890 28100];
censind = [0 1 1 1 0 1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 0 1 0 1 0 1 1];

fprintf('Shock absorber DRAM test:\n');
% check for covariates data
if (isempty(params.x))
    % Simulate artificial covariates
    params.x = randn(params.D, numel(y))*1/params.D;
    fprintf('\tcovariates: randomly synthesized, mean(x)=%g\n', mean(params.x(:)));
else
    fprintf('\tcovariates: given, mean(x)=%g\n', mean(params.x(:)));
end
fprintf('\tNumber of samples: N=2^%d\n', params.log2N);
fprintf('\tNumber of covariates: D=%d\n', params.D);
fprintf('\tNumber of runs: nruns=%d\n', params.runs);

% scales for beta
beta_mean = zeros(1,params.D+1);
beta_mean(1) = log(30796);
beta_var = ones(1,params.D+1);
beta_var(1) = 0.1563;

% Ranges of domain
a = beta_mean - 3*sqrt(beta_var);
b = 2*beta_mean - a;
a = [a 0]';
b = [b 13]';


% Allocate storage
ttimes_dram = zeros(params.runs, 1); % CPU times
Q_dram = zeros(params.runs, 2); % quantities of interest
num_of_rejects = zeros(params.runs, 1); % number of rejections in DRAM
tauint_dram = zeros(params.runs, 1); % IACT

for irun=1:params.runs
    % DRAM parameters
    options = struct();
    options.nsimu    = 2^params.log2N; % number of samples
    options.adaptint = 10;
    options.drscale  = 2;
    options.adascale = 2.4/sqrt(params.D+2); % scale for adaptation
    options.qcov     = eye(params.D+2)*5; % initial covariance
    options.printint = 2^(params.log2N-4);
    
    % create input arguments for the dramrun function
    model = struct();
    model.ssfun    = @(theta,ddd)(-2*(shock_log_weibull(theta, params.x, y, censind)+shock_log_prior(theta, beta_mean, beta_var)));
    data = struct();
    domain = struct();
    domain.par0    = [beta_mean 2]; % initial value
    domain.bounds = [a, b]';
    
    tic;
    [results,Z_dram] = dramrun(model,data,domain,options);
    ttimes_dram(irun) = toc;
    num_of_rejects(irun) = (1-results.accepted)*(2^params.log2N);

    % Remove burn-in
    Z_dram = Z_dram(((2^params.log2N)/4+1):(2^params.log2N), :);

    Q_dram(irun,:) = shock_quantiles(Z_dram, 0*ones(params.D,1));
        
    % Estimate autocorrr time
    [~,~,~,tauint_dram(irun),~,~] = UWerr(Z_dram,1.5,length(Z_dram),0);
    tauint_dram(irun) = tauint_dram(irun)*2;        
end

% Print some statsy information
fprintf('DRAM Shock absorber completed. Some average values:\n');
fprintf('\tttimes_dram: %g\n', mean(ttimes_dram));
fprintf('\tnum_of_rejects: %g (out of N=%d)\n', mean(num_of_rejects), 2^params.log2N);
fprintf('\tIACT: %g\n', mean(tauint_dram));
fprintf('\tQ_dram: [%g %g]+-[%g %g]\n', mean(Q_dram(:,1)), mean(Q_dram(:,2)), sqrt(sum((Q_dram(:,1)-mean(Q_dram(:,1))).^2)/(params.runs-1)),  sqrt(sum((Q_dram(:,2)-mean(Q_dram(:,2))).^2)/(params.runs-1)));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
