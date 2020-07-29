% DIRT Inverse Diffusion test
function test_diffusion_dirt(varargin)
% Check for and download TT-Toolbox
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    cd(mydir); cd('..'); cd('..'); cd('utils'); check_ttirt;
end
check_tt;
cd(mydir);

% Parse parameters or ask a user for them
params = parse_diffusion_inputs(varargin{:});
% Extra parameters (only for TT)
if (~isfield(params, 'ny'))
    params.ny = input('Gauss grid size for the forward map ny = ? (default 7): ');
    if (isempty(params.ny))
        params.ny = 7;
    end
end
if (~isfield(params, 'rmax'))
    params.rmax = input('Max TT rank rmax = ? (default 800): ');
    if (isempty(params.rmax))
        params.rmax = 800;
    end
end
if (~isfield(params, 'npi'))
    params.npi = input('Uniform grid size for the posterior npi = ? (default 17): ');
    if (isempty(params.npi))
        params.npi = 17;
    end
end
if (~isfield(params, 'rpi'))
    params.rpi = input('TT rank of ratios R = ? (default 8): ');
    if (isempty(params.rpi))
        params.rpi = 8;
    end
end
if (~isfield(params, 'beta'))
    params.beta = input('Tempering powers (default 10.^(-1:0.5:0)): ');
    if (isempty(params.beta))
        params.beta = 10.^(-1:0.5:0);
    end
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
ni = round(params.ny + (2-params.ny)*(ni/ni(L)));

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
ys = 2*sqrt(3)/(params.npi-1);
ys = (-sqrt(3):ys:sqrt(3))';
ys = repmat({ys}, L, 1);

ttimes_dirt = zeros(params.runs, 1);
ttimes_debias = zeros(params.runs, 1);
evalcnt = zeros(params.runs, 1);
num_of_rejects = zeros(params.runs, 1);
tau_ess = zeros(params.runs, 1);
tau_tt = zeros(params.runs, 1);

% Fire up Bayesian runs
for irun=1:params.runs
    % Tempered Log(posterior) function   
    lpfun = @(theta,beta0,beta)-sum((tt_sample_lagr(u_av{irun},Y,theta) - Q_obs).^2, 2)*(beta-beta0)/(2*params.sigma_n);
    
    % Construct DIRT
    tic;
    IRTstruct = tt_dirt_approx(ys, lpfun, params.beta, 'testsamples', 1e2, ...
                               'nswp', 1, 'y0', params.rpi, 'kickrank', 0, 'boundary', true, ...
                               'reference', 'n4', 'interpolation', 'f');
    ttimes_dirt(irun) = toc;
    % Total number of evals
    evalcnt(irun) = sum(IRTstruct.evalcnt-1e2);    

    % Sample
    tic;
    q = randref(IRTstruct.reference, 2^params.log2N, L);                           
    [z,lFapp,lFex] = tt_dirt_sample(IRTstruct, q, @(x)lpfun(x,0,1));
    ttimes_debias(irun) = toc;
    
    % Reject (a.k.a. check error)
    [z2,~,~,num_of_rejects(irun)] = mcmc_prune(z,lFex,lFapp);
    % IACT and ESS
    [~,~,~,tau_tt(irun),~,~] = UWerr(z2,1.5,length(z2),0);
    tau_tt(irun) = tau_tt(irun)*2;
    tau_ess(irun) = essinv(lFex,lFapp);    
end % irun

% Print some statsy information
fprintf('DIRT Diffusion completed. Some average values:\n');
fprintf('\tCPU time of forward model: %g\n', mean(ttimes_forward));
fprintf('\tnsolves in forward model: %g\n', mean(nsolves_forward));
fprintf('\tCPU time of DIRT: %g\n', mean(ttimes_dirt));
fprintf('\tNumber of evaluations in DIRT: %g\n', mean(evalcnt));
fprintf('\tCPU time of exact pi: %g\n', mean(ttimes_debias));
fprintf('\tnum_of_rejects: %g (out of N=%d)\n', mean(num_of_rejects), 2^params.log2N);
fprintf('\tIACT: %g\n', mean(tau_tt));
fprintf('\tN/ESS: %g\n', mean(tau_ess));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
