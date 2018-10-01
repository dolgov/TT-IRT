% TT-IRT Inverse Diffusion test
function test_diffusion_tt(varargin)
% Check for and download TT-Toolbox
check_tt;

% Parse parameters or ask a user for them
params = parse_diffusion_inputs(varargin{:});
% Extra parameters (only for TT)
if (~isfield(params, 'ny'))
    params.ny = input('Gauss grid size for the forward map ny = ? (default 7): ');
    if (isempty(params.ny))
        params.ny = 7;
    end
end
if (~isfield(params, 'npi'))
    params.npi = input('Uniform grid size for the posterior npi = ? (default 32): ');
    if (isempty(params.npi))
        params.npi = 32;
    end
end
if (~isfield(params, 'delta'))
    params.delta = input('TT approx. threshold for the posterior delta = ? (default 0.1): ');
    if (isempty(params.delta))
        params.delta = 0.1;
    end
end
if (~isfield(params, 'rmax'))
    params.rmax = input('Max TT rank rmax = ? (default 800): ');
    if (isempty(params.rmax))
        params.rmax = 800;
    end
end
if (~isfield(params, 'correction'))
    params.correction = input('Algorithm of debiasing correction = ? (''mcmc'' or ''iw'') (default ''mcmc''): ');
    if (isempty(params.correction))
        params.correction = 'mcmc';
    end
end

if (strcmpi(params.correction, 'mcmc'))
    % We need DRAM and UWerr
    check_mcmc;
else
    % We need QMC generating vector
    check_qmc;
end

% A priori fitted function to map spatial meshlevel into space discr. error
htolfun = @(x)(7.6742e-03*4^(-x-1)); % For the flux in the inverse problem


tol = htolfun(params.meshlevel);
fprintf('Solving for lvl=%d, tol=%3.3e, ny=%d\n', params.meshlevel, tol, params.ny);

% Build the discretization and KLE
tol_kle = tol*3;
[~,bound,W1g,W1m,spind,Pua,phi,lambda,Mass] = build_grid_and_kle(params.meshlevel, 'DN', params.nu, params.corr_length, tol_kle, params.m0);


% weighted KLE components
L = numel(lambda);
phil = full(phi*spdiags(sqrt(lambda), 0, L, L));

% Anisotripic grid in parameters
ni = log(lambda);
ni = round(params.ny + (1-params.ny)*(ni/ni(L)));

Y = cell(L,1);
for i=1:L
    % Parametric grid points and weights
    [Y{i},~]=lgwt(ni(i),-sqrt(3),sqrt(3));
end

u_av = cell(1, params.runs);  % Storage for prior observations
ttimes_forward = zeros(params.runs, 1);
nsolves_forward = zeros(params.runs, 1);
% Run the forward model, compute prior observables
for irun=1:params.runs
    tic;
    
    % Create the affine expansion.
    % Within the inner loop to estimate the cpu time
    log_a = [];
    for i=1:L
        af = tt_tensor(phil(:,i)*sqrt(params.sigma));
        for j=1:L
            if (j==i)
                af = tkron(af, tt_tensor(Y{i}));
            else
                af = tkron(af, tt_ones(ni(j)));
            end
        end
        log_a = af+log_a;
    end
    
    % Determine tolerance for a based on its range
    log_a_bound = tt_stat(log_a, 'sr','lr');
    log_a_bound = exp(log_a_bound(2)-log_a_bound(1));
    tol_a = min(1/log_a_bound, tol);
    % Create log-(normal or uniform) coefficient via TT-Cross
    af = amen_cross_s({log_a}, @(x)exp(x), tol_a, 'y0', params.rmax, 'nswp', 1, 'kickrank', 0);
    
    % ALS-Cross solver is here
    [u,~,nsolves_forward(irun)] = als_cross_parametric(af, @(C)diffusion_assem_solve(C,bound,W1g,W1m,spind), tol, 'Pua', Pua, 'random_init', params.rmax, 'nswp', 1, 'kickrank', 0);
    %       ^^ number of deterministic pde solves
    
    % Compute observables
    u1 = u{1};
    u1 = reshape(u1, size(u1,2), size(u1,3));
    u_av_1 = zeros(params.m0^2, size(u1,2));
    for j=1:params.m0
        for i=1:params.m0
            u_av_1(i+(j-1)*params.m0, :) = sum(Mass{i,j}*u1, 1);
        end
    end
    u_av{irun} = u_av_1*chunk(u,2,L+1);
    
    ttimes_forward(irun) = toc; % time to compute forward map
end % irun

% Simulate some observations
if (exist(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0), 'file')>0)
    % Load the same KLE for all experiments
    fprintf('Found Q_obs file for nu=%g, ml=%d, sn=%g, m0=%d, y0=%g\nRemove it to regenerate the observations\n', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0);
    load(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0));
else
    fprintf('Generating Q_obs from the current solution\n');
    Q_obs = tt_sample_lagr(u_av{1}, Y, params.y0*ones(1,L)) + randn(1,params.m0^2)*sqrt(params.sigma_n);
    save(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0), 'Q_obs');
end

% Interpolate from lgwt to uniform
ys = 2*sqrt(3)/params.npi;
ys = (-sqrt(3)+ys:ys:sqrt(3))';
ys = repmat({ys}, L, 1);
P_lgwt_uni = cell(L,1);
for i=1:L
    P_lgwt_uni{i} = reshape(lagrange_interpolant(Y{i}, ys{i}), 1, params.npi, ni(i));
    ys{i} = [-sqrt(3); ys{i}]; % left boundary is needed for IRT
end
P_lgwt_uni = cell2core(tt_matrix, P_lgwt_uni);


err_Pi = zeros(params.runs, 1);
ttimes_pi = zeros(params.runs, 1);
ttimes_invcdf = zeros(params.runs, 1);
ttimes_debias = zeros(params.runs, 1);
bias = zeros(params.runs, 1);
Q_tt = zeros(params.runs, 2);
tau_tt = zeros(params.runs, 3);

% Fire up Bayesian runs
Pi = cell(1,params.runs);
for irun=1:params.runs
    u_av{irun} = P_lgwt_uni*u_av{irun};
    u_av{irun} = tt_reshape(u_av{irun}, repmat(factor(params.npi)', L, 1), tol*1e-2, params.m0^2, 1);
    
    % Compute approx posterior density
    err_Pi(irun) = params.delta;
    tic;
    pi = amen_cross_s(u_av(irun), @(x)exp(-sum((x-repmat(Q_obs,size(x,1),1)).^2, 2)/(2*params.sigma_n)), err_Pi(irun), 'y0', params.rmax, 'kickrank', 0, 'nswp', 1);
    pi = tt_reshape(pi, params.npi*ones(L,1));
    ttimes_pi(irun) = toc;
    Pi{irun} = pi;
    
    % TT-MH or TT-IW algorithms
    if (strcmpi(params.correction, 'mcmc'))
        Z = rand(2^params.log2N, L);
    else
        Z = qmcnodes(L, params.log2N)';
    end
    tic;
    [y,Fex,bias(irun),ttimes_invcdf(irun)] = tt_irt_debias(Z, @(y)diffusion_likelihood(y, phil, params.sigma, bound, W1g, W1m, spind, Mass, Q_obs, params.sigma_n), pi, ys, params.correction);
    ttimes_debias(irun) = toc;
    
    Q_tt(irun,:) = mean(Fex(:,2:3));
    
    % Estimate IACT
    if (strcmpi(params.correction, 'mcmc'))
        [~,~,~,tauint_y,~,~] = UWerr(y,1.5,length(y),0);
        [~,~,~,tauint_f,~,~] = UWerr(Fex(:,2),1.5,length(Fex),0);
        [~,~,~,tauint_p,~,~] = UWerr(Fex(:,3),1.5,length(Fex),0);
        tau_tt(irun, :) = [tauint_y  tauint_f  tauint_p]*2;
    end
end % irun

% Estimate error in Pi
pi = amen_sum(Pi, (1/params.runs)*ones(params.runs,1), 1e-4, 'y0', Pi{1}, 'fkick', true, 'kickrank', 64);
for irun=1:params.runs
    err_Pi(irun) = norm(Pi{irun} - pi);
end
err_Pi = err_Pi/norm(pi);


% Print some statsy information
fprintf('TT Diffusion completed. Some average values:\n');
fprintf('\tCPU time of forward model: %g\n', mean(ttimes_forward));
fprintf('\tnsolves in forward model: %g\n', mean(nsolves_forward));
fprintf('\tCPU time of approx. pi: %g\n', mean(ttimes_pi));
fprintf('\tCPU time of IRT: %g\n', mean(ttimes_invcdf));
fprintf('\tCPU time of exact pi: %g\n', mean(ttimes_debias));
if (strcmpi(params.correction, 'mcmc'))
    fprintf('\tnum_of_rejects: %g (out of N=%d)\n', mean(bias), 2^params.log2N);
    fprintf('\tIACT for [y, F, P]: [%g %g %g]\n', mean(tau_tt(:,1)), mean(tau_tt(:,2)), mean(tau_tt(:,3)));
else
    fprintf('\tstd of importance weights: %g\n', mean(bias));
end
fprintf('\tTT error for pi: %g\n', mean(err_Pi));
fprintf('\tQ_tt: [%g %g]+-[%g %g]\n', mean(Q_tt(:,1)), mean(Q_tt(:,2)), sqrt(sum((Q_tt(:,1)-mean(Q_tt(:,1))).^2)/(params.runs-1)),  sqrt(sum((Q_tt(:,2)-mean(Q_tt(:,2))).^2)/(params.runs-1)));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
