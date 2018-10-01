% Bayesian ratio + QMC Inverse Diffusion test
function test_diffusion_qmcrat(varargin)
% Download prerequisites if necessary
check_mcmc;
check_qmc;

% Parse parameters or ask a user for them
params = parse_diffusion_inputs(varargin{:});

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
    fprintf('Found Q_obs file for nu=%g, ml=%d, sn=%g, m0=%d, y0=%g\nRemove it to regenerate the observations\n', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0);
    load(sprintf('Q_obs_nu%g_ml%d_sigman%g_m0%d_ytrue%g.mat', params.nu, params.meshlevel, params.sigma_n, params.m0, params.y0));
    Q_obs = Q_obs + 0; % Parallel toolbox still sucks
else
    error('Please generate observation data by running test_diffusion_tt');
end

Q_qmcrat = zeros(params.runs, 2);
ttimes_qmcrat = zeros(params.runs, 1);

for irun=1:params.runs
    tic;
    Z = qmcnodes(L, params.log2N); % size L x I
    Z = (Z-0.5)*2*sqrt(3); % map to [-sqrt(3),sqrt(3)]
        
    Fsum = zeros(1,3); % total quadrature
    
    % Implement a loop here to save memory
    for k=1:(2^params.log2N)        
        Fex = diffusion_likelihood(Z(:,k)', phil, params.sigma, bound, W1g, W1m, spind, Mass, Q_obs, params.sigma_n);
        Fsum(1) = Fsum(1) + Fex(1);
        Fsum(2) = Fsum(2) + Fex(2)*Fex(1); % Fex(1)/1 is the ratio pi/uniform
        Fsum(3) = Fsum(3) + Fex(3)*Fex(1);                
        
        if (mod(k,100)==0)
            fprintf('QMC-rat solved problem %d out of %d\n', k, 2^params.log2N);
        end
    end
    ttimes_qmcrat(irun) = toc;
    Q_qmcrat(irun,:) = Fsum(2:3)/Fsum(1);
end

% Print some statsy information
fprintf('QMC-rat Diffusion completed. Some average values:\n');
fprintf('\tCPU time: %g\n', mean(ttimes_qmcrat));
fprintf('\tQ_qmcrat: [%g %g]+-[%g %g]\n', mean(Q_qmcrat(:,1)), mean(Q_qmcrat(:,2)), sqrt(sum((Q_qmcrat(:,1)-mean(Q_qmcrat(:,1))).^2)/(params.runs-1)),  sqrt(sum((Q_qmcrat(:,2)-mean(Q_qmcrat(:,2))).^2)/(params.runs-1)));


% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end

