% DRAM Inverse Diffusion test
function test_diffusion_dram(varargin)
% Download DRAM if necessary
check_mcmc;

% Parse parameters or ask a user for them
params = parse_model_inputs(varargin{:});

% A priori fitted function to map spatial meshlevel into space discr. error
htolfun = @(x)(7.6742e-03*4^(-x-1)); % For the flux in the inverse problem


% PDE params
tol = htolfun(params.meshlevel);
fprintf('Solving for lvl=%d, tol=%3.3e\n', params.meshlevel, tol);

% Build the discretization and KLE
tol_kle = tol*3;
[~,bound,W1g,W1m,spind,Pua,phi,lambda,Mass] = build_grid_and_kle(params.meshlevel, 'DN', params.nu, params.corr_length, tol_kle, params.m0);
% weighted KLE components
L = numel(lambda);
phil = full(phi*spdiags(sqrt(lambda), 0, L, L));

% Simulate some observations
if (exist(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0), 'file')>0)
    % Load the same KLE for all experiments
    fprintf('Found Obs file for nu=%g, ml=%d, sn=%g, m0=%d, y0=%g\nRemove it to regenerate the observations\n', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0);
    load(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0));
    Q_obs = Q_obs + 0; % Parallel toolbox still sucks
else
    error('Please generate observation data by running test_diffusion_tt');
end

Q_dram = zeros(params.runs, 2);
ttimes_dram = zeros(params.runs, 1);
tau_dram = zeros(params.runs, 3);
num_of_rejects = zeros(params.runs, 1);

parfor irun=1:params.runs
    % DRAM parameters
    options = struct();
    options.nsimu    = 2^params.log2N; % number of samples
    options.adaptint = 10;
    options.drscale  = 2;
    options.adascale = 2.4/sqrt(L); % scale for adaptation
    options.qcov     = eye(L)*5; % initial covariance
    
    % create input arguments for the dramrun function
    model = struct();    
    model.ssfun    = @(y,d)(-2*log(diffusion_likelihood(y, phil, params.sigma, bound, W1g, W1m, spind, Mass, Q_obs, params.sigma_n)*[1;zeros(2,1)]));
    
    data = struct();
    domain = struct();
    domain.par0    = zeros(1,L); % initial value
    domain.bounds = sqrt(3)*[-ones(1,L); ones(1,L)];
    
    tic;
    [results,y] = dramrun(model,data,domain,options);
    ttimes_dram(irun) = toc;
    
    num_of_rejects(irun) = (1-results.accepted)*(2^params.log2N);
    
    N_burn_in = 6000; % hardcode it for now
    y = y((N_burn_in+1):(2^params.log2N), :);  % eliminate the burn-in
    
    Fex = diffusion_likelihood(y, phil, params.sigma, bound, W1g, W1m, spind, Mass, Q_obs, params.sigma_n);
    Q_dram(irun,:) = mean(Fex(:,2:3));
    
    % Estimate autocorrr time
    [~,~,~,tauint_y,~,~] = UWerr(y,1.5,length(y),0);
    [~,~,~,tauint_f,~,~] = UWerr(Fex(:,2),1.5,length(Fex),0);
    [~,~,~,tauint_p,~,~] = UWerr(Fex(:,3),1.5,length(Fex),0);
    tau_dram(irun, :) = [tauint_y  tauint_f  tauint_p]*2;
end

% Print some statsy information
fprintf('DRAM Diffusion completed. Some average values:\n');
fprintf('\tCPU time of DRAM: %g\n', mean(ttimes_dram));
fprintf('\tnum_of_rejects: %g (out of N=%d)\n', mean(num_of_rejects), 2^params.log2N);
fprintf('\tIACT for [y, F, P]: [%g %g %g]\n', mean(tau_dram(:,1)), mean(tau_dram(:,2)), mean(tau_dram(:,3)));
fprintf('\tQ_dram: [%g %g]+-[%g %g]\n', mean(Q_dram(:,1)), mean(Q_dram(:,2)), sqrt(sum((Q_dram(:,1)-mean(Q_dram(:,1))).^2)/(params.runs-1)),  sqrt(sum((Q_dram(:,2)-mean(Q_dram(:,2))).^2)/(params.runs-1)));

end

