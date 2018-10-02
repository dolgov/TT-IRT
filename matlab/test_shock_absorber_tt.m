% Shock absorber tests with TT decomposition
% The following parameters can be passed in the form 
% 'param_name1', param_value1, 'param_name2', param_value2 and so on:
%   x: a matrix of covariates (can be empty, in this case drawn from normal
%      distribution)
%   log2N: log2(number of samples in the chain)
%   D: number of covariates
%   delta: TT approximation/stopping threshold
%   n: number of grid points in each variable
%   runs: number of runs of the test (must be >1 for error estimation)
% The test script does not require output parameters. Instead, it copies
% all variables to the main Matlab workspace. Among those are:
%   N_cross:  number of pi evaluations in different runs of TT Cross
%   ttimes_cross: CPU times of TT cross runs
%   ttimes_invcdf: CPU times of IRT runs
%   num_of_rejects: numbers of rejections in MH method
%   Q_mh: quantiles computed by MH method
%   Q_iw: quantiles computed by IW method with QMC
%   tauint_tt: integrated autocorrelation times in each run
%   err_TT: estimate of the TT approximation error (empirical std)
function test_shock_absorber_tt(varargin)
% Check prerequisites
check_tt;
check_qmc;
check_mcmc;

params = parse_shock_inputs(varargin{:});
% extra parameters for TT
if (~isfield(params, 'n'))
    params.n = input('Uniform grid size for the posterior n = ? (default 16): ');
    if (isempty(params.n))
        params.n = 16;
    end
end
if (~isfield(params, 'delta'))
    params.delta = input('TT approx. threshold for the posterior delta = ? (default 0.05): ');
    if (isempty(params.delta))
        params.delta = 0.05;
    end
end

% Data table from the paper
y = [6700 6950 7820 8790 9120 9660 9820 11310 11690 11850 11880 12140 12200 ...
     12870 13150 13330 13470 14040 14300 17520 17540 17890 18420 18960  ...
     18980 19410 20100 20100 20150 20320 20900 22700 23490 26510 27410  ...
     27490 27890 28100];
censind = [0 1 1 1 0 1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 0 1 0 1 0 1 1];

fprintf('Shock absorber TT test:\n');
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
fprintf('\tTT approx. tolerance: delta=%g\n', params.delta);
fprintf('\tNumber of grid points: n=%d\n', params.n);
fprintf('\tNumber of runs: runs=%d\n', params.runs);

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

% Vectors of grid points
theta = cell(params.D+2,1);
for i=1:params.D+2
    h = (b(i)-a(i))/params.n;
    theta{i} = (a(i)+h:h:b(i))';
    theta{i} = tt_tensor(theta{i});
end
% Replicate those in all dimensions as TT tensors
% These must be without boundary points
theta_tt = tt_meshgrid_vert(theta);
% Add boundary points for the Rosenblatt transform
for i=1:params.D+2
    theta{i} = full(theta{i});
    theta{i} = [a(i); theta{i}];
end

% Exact posterior function
pifun = @(theta)exp(shock_log_weibull(theta, params.x, y, censind) + shock_log_prior(theta, beta_mean, beta_var));

n_err_TT = round(params.runs/4); % number of runs used to estimate error in pi
Pi = cell(1,n_err_TT); % storage for pi approximations from several runs

% Allocate storage
N_cross = zeros(params.runs,2);  % number of pi evaluations in TT Cross
ttimes_cross = zeros(params.runs, 1); % CPU time of TT cross
ttimes_invcdf = zeros(params.runs, 1); % times of IRT
num_of_rejects = zeros(params.runs, 1); % number of rejections in MH
Q_mh = zeros(params.runs, 2); % quantiles computed by MH method
Q_iw = zeros(params.runs, 2); % quantiles computed by IS method with QMC
tauint_tt = zeros(params.runs, 1); % IACT (after MH)


for irun=1:params.runs
    tic;
    [pi,~,~,~,N_cross(irun,:)] = amen_cross_s(theta_tt, pifun, 0, 'kickrank', 2, 'y0', 8, 'tol_exit', params.delta);
    ttimes_cross(irun) = toc;
    pi = round(pi, 1e-6); % truncate some ranks we can
    if (irun<=n_err_TT)
        Pi{irun} = pi;
    end

    [Z_mh,~,num_of_rejects(irun),ttimes_invcdf(irun)] = tt_irt_debias(2^params.log2N, pifun, pi, theta, 'mcmc');
    [Q_mh(irun,:)] = shock_quantiles(Z_mh, 0*ones(params.D,1));

    % Importance weighting with QMC
    Z_iw = qmcnodes(params.D+2, params.log2N)';
    [Z_iw,pi_app] = tt_irt(theta,pi,Z_iw);
    pi_ex = pifun(Z_iw);
    Q_iw(irun,:) = shock_quantiles(Z_iw, 0*ones(params.D,1), pi_ex./pi_app);

    % Estimate autocorr time of Q1
    [~,~,~,tauint_tt(irun),~,~] = UWerr(Z_mh,1.5,length(Z_mh),0);
    tauint_tt(irun) = tauint_tt(irun)*2;    
end
N_cross = N_cross(:,2);

% Estimate the actual TT approximation error
n_err_TT = min(n_err_TT,params.runs);
err_TT = nan;
if (n_err_TT>1)
    pi = amen_sum(Pi, (1/n_err_TT)*ones(n_err_TT,1), 1e-6, 'fkick', true, 'y0', Pi{1}, 'kickrank', 16);
    err_TT = norm(pi)*ones(n_err_TT,1);
    for irun=1:n_err_TT
        err_TT(irun) = norm(Pi{irun} - pi)/err_TT(irun);
    end
    err_TT = sqrt(sum(err_TT.^2)/(n_err_TT-1)); % empirical standard deviation
end


% Print some statsy information
fprintf('TT Shock absorber completed. Some average values:\n');
fprintf('\tN_cross: %g\n', mean(N_cross));
fprintf('\tttimes_cross: %g\n', mean(ttimes_cross));
fprintf('\tttimes_invcdf: %g\n', mean(ttimes_invcdf));
fprintf('\tnum_of_rejects: %g (out of N=%d)\n', mean(num_of_rejects), 2^params.log2N);
fprintf('\tIACT: %g\n', mean(tauint_tt));
fprintf('\tTT error: %g\n', err_TT);
fprintf('\tQ_mh: [%g %g]+-[%g %g]\n', mean(Q_mh(:,1)), mean(Q_mh(:,2)), sqrt(sum((Q_mh(:,1)-mean(Q_mh(:,1))).^2)/(params.runs-1)),  sqrt(sum((Q_mh(:,2)-mean(Q_mh(:,2))).^2)/(params.runs-1)));
fprintf('\tQ_iw: [%g %g]+-[%g %g]\n', mean(Q_iw(:,1)), mean(Q_iw(:,2)), sqrt(sum((Q_iw(:,1)-mean(Q_iw(:,1))).^2)/(params.runs-1)),  sqrt(sum((Q_iw(:,2)-mean(Q_iw(:,2))).^2)/(params.runs-1)));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
