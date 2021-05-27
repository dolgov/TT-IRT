function [IRTstruct] = tt_dirt_approx(x0, logpostfun, beta, varargin)
% The constructor of the Deep Inverse Rosenblatt Transform
%  function [IRTstruct] = TT_DIRT_APPROX(x0, logpostfun, beta, varargin)
% Inputs:
%   x0: d x 1 cell array {x1...xd} of initial grid vectors on level 0
%   logpostfun: @(x,\beta_min,\beta_max) function for producing
%               log(density ratios) lpi_{\beta_max}(x) - lpi_{\beta_min}(x)
%   beta: array of bridging parameters (e.g. reciprocal temperatures)
%         in increasing order
%
% Optional inputs given in form 'ParamName', ParamValue, ...
%   nq: number of grid points in levels 1 to L (default same as in x0)
%   stoptol: stopping tolerance for TT cross (default 0.4)
%   trunctol: truncation tolerance for TT cross (default 0)
%   crossmethod: TT cross algorithm: 'amen_cross_s' (default supplied), or
%                                    'greedy2_cross' (from TT-Toolbox), or
%                                    'build_ftt'    (from ftt.m)
%   y0: Initial TT rank in TT cross (default 1)
%   kickrank: rank enrichment step in TT cross (default 4)
%   nswp: number of TT cross sweeps (default 4)
%   vec: whether logpostfun is vectorized (default true)
%   boundary: should densities be evaluated on boundary points [false]
%   testsamples: test approximations via MCMC using min(N,testsamples) 
%                extra samples per level [1e4]
%   recompute: if N/ESS > recompute, recompute the current level  [50]
%   reference: reference density (['UNIform'], 'Normal' or 'Normal S',
%              where S is the number of sigmas defining the (-S,S] support
%              of the truncated normal reference variables. S=4 if absent)
%   IRTdenom: if true, divide by IRT density in the next importance ratio,
%             otherwise take a ratio of two exact densities    [false]
%   plotdiag: plot diagnostic marginals and samples            [true]
%   interpolation:  'spline' (linear) or 'fourier'             ['spline']
%
% Outputs:
%   IRTstruct: a structure for sampling via tt_dirt_sample 
%              or further transformation
%
% See also: amen_cross_s, tt_dirt_sample


% Parse parameters
stoptol = 0.4;
trunctol = 0;
nq = [];
y0 = 1;
kickrank = 4;
nswp = 4;
vec = true;
boundary = false;
testsamples = 1e4;
recompute = 50;
reference = 'uni';
crossmethod = 'amen_cross_s';
IRTdenom = false;
plotdiag = true;
interpolation = 'spline';
IRTstruct = [];
for i=1:2:numel(varargin)
    switch (lower(varargin{i}))
        case 'stoptol'
            stoptol = varargin{i+1};
        case 'trunctol'
            trunctol = varargin{i+1};            
        case 'nq'
            nq = varargin{i+1};
        case 'y0'
            y0 = varargin{i+1};
        case 'kickrank'
            kickrank = varargin{i+1};
        case 'nswp'
            nswp = varargin{i+1};
        case 'vec'
            vec = varargin{i+1};
        case 'boundary'
            boundary = varargin{i+1};
        case 'testsamples'
            testsamples = varargin{i+1};
        case 'recompute'
            recompute = varargin{i+1};
        case 'irtdenom'
            IRTdenom = varargin{i+1};            
        case 'reference'
            reference = lower(varargin{i+1});
        case 'crossmethod'
            crossmethod = lower(varargin{i+1});
        case 'plotdiag'
            plotdiag = varargin{i+1};
        case 'interpolation'
            interpolation = lower(varargin{i+1});
        case 'irtstruct'
            IRTstruct = varargin{i+1};
        otherwise
            error('unknown optional parameter %s', varargin{i});
    end
end

nlvl = numel(beta)-1; % number of increment levels
d = numel(x0); % number of variables

% Allow vector-valued cross parameters to fine-tune different levels
if (numel(nswp)==1)
    nswp = repmat(nswp, 1, nlvl+1);
end
if (numel(kickrank)==1)
    kickrank = repmat(kickrank, 1, nlvl+1);
end
if (numel(stoptol)==1)
    stoptol = repmat(stoptol, 1, nlvl+1);
end
if (numel(trunctol)==1)
    trunctol = repmat(trunctol, 1, nlvl+1);
end
if (numel(IRTdenom)==1)
    IRTdenom = repmat(IRTdenom, 1, nlvl+1);
end
if (size(y0,1)==1)
    y0 = repmat(y0, d+1, 1);    % we may need to fine-tune ranks
end                             % over dimensions and levels
if (size(y0,2)==1)
    y0 = repmat(y0, 1, nlvl+1);
end

if (interpolation(1)~='s')&&(~boundary)
    boundary = true;
    warning('Overriding boundary->true for Fourier interpolation');
end

% Zero level needs special treatment, as the box for X can be arbitrary
if (isa(x0{1}, 'struct'))&&(strcmp(crossmethod, 'build_ftt'))
    % This is for build_ftt, x0 is a cell of oned polys
    if (isempty(nq))
        nq = cellfun(@(p)p.order, x0);
    end
else
    % This is for gridpoints and TT-Toolbox
    if (isempty(nq))
        nq = cellfun(@numel, x0);
    end
    if (boundary)
        X = tt_meshgrid_vert(cellfun(@(x)tt_tensor(x), x0, 'UniformOutput', false));
    else
        X = tt_meshgrid_vert(cellfun(@(x)tt_tensor(x(2:end-1)), x0, 'UniformOutput', false));
    end
end
if (numel(nq)==1)
    nq = repmat(nq, d, 1);
end


% Storage for densities and functions
F = cell(nlvl, 1);
if (isempty(IRTstruct))
    % Initialize empty IRT structure
    ilvl = 0;
    IRTstruct = struct();
    IRTstruct.x0 = x0;
    IRTstruct.beta = beta(1);
    IRTstruct.reference = reference;
    IRTstruct.crossmethod = crossmethod;
    IRTstruct.interpolation = interpolation;
    % Initial guess
    Fprev = max(y0(:, min(2,size(y0,2))));
else
    ilvl = numel(IRTstruct.beta);
    if (ilvl>1)
        F(1:ilvl-1) = IRTstruct.F(1:ilvl-1);
    end
    beta(1:ilvl) = IRTstruct.beta;
    lFshift = IRTstruct.lFshift;
    Fprev = IRTstruct.Fprev;
end

if (ilvl==0)
    fprintf('Approximating level 0, for beta=%g\n', beta(1));
    switch (crossmethod)
        case 'amen_cross_s'
            [F0,~,~,~,evalcnt1] = amen_cross_s(X, @(x)exp(logpostfun_vec(x, 0,beta(1), logpostfun, vec)*0.5), ...
                trunctol(1), 'tol_exit', stoptol(1), 'y0', max(y0(:,1)), 'kickrank', kickrank(1), 'nswp', nswp(1), 'verb', 1);
        case 'greedy2_cross'
            nx0 = cellfun(@numel, x0);
            y0mid = round((nx0-1)/2);
            if (~boundary)      % greedy2_cross is defined by indices,
                nx0 = nx0 - 2;  % give it correct sizes accordingly
            end
            [F0,~,~,~,~,evalcnt1]=greedy2_cross(nx0, @(i)zeros(size(i,1),1), 1e-12, 'nswp', nswp(1), 'y0', y0mid, 'verb', true, ...
                'aux', X, 'tol_exit', stoptol(1), 'auxfun', @(x)exp(logpostfun_vec(x, 0,beta(1), logpostfun, vec)*0.5));
            
        case 'build_ftt'
            debug_size = max(y0(:,1)) + kickrank(1)*nswp(1);
            debug_x = zeros(d, debug_size);
            for k = 1:d
                debug_x(k,:) = sample_oned_domain(x0{k}, debug_size);
            end
            opts = ftt_options('method', 'AMEN', 'ng_flag', false, 'oned_ref', x0, ...
                'err_tol', stoptol(1), 'loc_err_tol', trunctol(1), 'max_als', nswp(1), ...
                'kick_rank', kickrank(1), 'max_rank', 50, 'init_rank', max(y0(:,1)));
            F0 = build_ftt(@(x)exp(logpostfun_vec(x', 0,beta(1), logpostfun, vec)*0.5)', d, [], opts, 'sample_x', debug_x, 'debug_x', debug_x);
            F0.ng_flag = true; % set this manually after constructing the sqrt(function) explicitly
            evalcnt1 = nan;
            F0 = build_irt(F0);
    end
    
    IRTstruct.evalcnt = [sum(evalcnt1); zeros(nlvl,1)];
    
    if (isa(F0, 'tt_tensor'))
        if (plotdiag)
            % draw 1D marginals
            Fdiag = zeros(max(F0.n), d);
            Fdiag(1:F0.n(1),1) = full(dot(tt_ones(F0.n(2:d)), F0, 2, d));
            for i=2:d-1
                Fdiag(1:F0.n(i),i) = full(dot(tt_ones(F0.n(i+1:d)), dot(tt_ones(F0.n(1:i-1)), F0, 1, i-1), 2, d-i+1));
            end
            Fdiag(1:F0.n(d),d) = full(dot(tt_ones(F0.n(1:d-1)), F0, 1, d-1));
            figure(1);
            plot(Fdiag);
            legend toggle;
            title('1D marginal sqrt(densities)');
            
            % Try drawing the x_1 x_2 marginal PDF
            figure(2);
            if (d==2)
                surf(full(F0, F0.n(1:2)'), 'EdgeColor', 'none');
            else
                surf(full(dot(tt_ones(F0.n(3:end))/prod(F0.n(3:end)), F0, 3, d), F0.n(1:2)'), 'EdgeColor', 'none');
            end
            shading interp;
            title('2D x_1 x_2 marginal')
            drawnow;
        end
        % Desintegrate tt_tensor into cell array for faster referencing
        F0 = core2cell(F0);
    end
    
    
    % Populate the DIRT structure with the zeroth level
    IRTstruct.F0 = F0;
    IRTstruct.Fprev = max(y0(:, min(2,size(y0,2))));
    
    lFshift = 0;
    if (testsamples>0)
        % Test approximation
        y = randref(reference,  min(sum(evalcnt1),testsamples),d);
        [y,lFapp,lFex] = tt_dirt_sample(IRTstruct, y, @(x)logpostfun(x,0,beta(1)), vec);
        [y2,~,~,num_of_rejects] = mcmc_prune(y, lFex, lFapp);
        num_of_rejects = num_of_rejects*100/size(y,1);
        tau = essinv(lFex, lFapp);
        fprintf('N/ESS = %g\n\n', tau);        
        if (plotdiag)
            figure(3);
            plot(y2);
            title(sprintf('Chain: #rejects = %g%%, N/ESS = %g', num_of_rejects, tau));
            drawnow;
        end
        IRTstruct.evalcnt(1) = IRTstruct.evalcnt(1) + size(y,1);
        lFshift = max(lFex); % we will subtract this to prevent overflows
        if (IRTdenom(1))
            lFshift = lFshift - max(lFapp);
        end
        IRTstruct.lFshift = lFshift;
    end
    
    ilvl = ilvl+1;
end

if (reference(1)~='u')
    % Parse the domain size of truncated normal
    sigma = double(reference);
    sigma = sigma((sigma==46) | ((sigma>=48) & (sigma<=57))); % extract numbers and dot
    sigma = str2double(char(sigma));
    if (isnan(sigma))
        sigma = 4;
    end
    fprintf('Using normal reference on [%g,%g]\n', -sigma, sigma);
end

% Set up 1D ansatz for other levels
if (strcmp(crossmethod, 'build_ftt'))
    if (reference(1)=='u')
        x = arrayfun(@(n)setup_oned(n, 'type', 'Chebyshev1st', 'domain', [0,1], 'bc', 'Neumann'), nq, 'UniformOutput', false);
    else
        x = arrayfun(@(n)setup_oned(n, 'type', 'Fourier', 'domain', [-sigma,sigma]), nq, 'UniformOutput', false);
    end
else
    if (reference(1)=='u')
        x = arrayfun(@(n)0.5*(cos(pi*((n-1):-1:0)'/(n-1))+1), nq, 'UniformOutput', false);
    else
        if (interpolation(1)=='s')
            x = arrayfun(@(n)(0:1/(n-1):1)'*(2*sigma)-sigma, nq, 'UniformOutput', false);
        else
            x = arrayfun(@(n)(1:n)'*(2*sigma/n)-sigma, round(nq/2)*2, 'UniformOutput', false);
        end
    end
    if (boundary)
        X = tt_meshgrid_vert(cellfun(@(x)tt_tensor(x), x, 'UniformOutput', false));
    else
        X = tt_meshgrid_vert(cellfun(@(x)tt_tensor(x(2:end-1)), x, 'UniformOutput', false));
    end
    x = cell2mat(x);
end
IRTstruct.x = x;

recompute_count = 0;  % Count the number of unsuccessful TT-Cross'es

while (ilvl<=nlvl)
    fprintf('Approximating level %d, for beta=%g\n', ilvl, beta(ilvl+1));    
    switch (crossmethod)
        case 'amen_cross_s'
            [F{ilvl},~,~,~,evalcnt1] = amen_cross_s(X, @(x)dualbetafun(x,IRTstruct,logpostfun,beta(ilvl+1),beta(ilvl),lFshift, vec, reference, IRTdenom(ilvl+1)), ...
                                       trunctol(ilvl+1), 'tol_exit', stoptol(ilvl+1), 'y0', Fprev, 'kickrank', kickrank(ilvl+1), 'nswp', nswp(ilvl+1), 'verb', 1);
            
        case 'greedy2_cross'
            nx0 = nq;
            y0mid = round((nx0-1)/2);
            if (~boundary)
                nx0 = nx0 - 2;
            end
            [F{ilvl},~,~,~,~,evalcnt1]=greedy2_cross(nx0, @(i)zeros(size(i,1),1), 1e-12, 'nswp', nswp(ilvl+1), 'y0', y0mid, 'verb', true, ...
                'aux', X, 'tol_exit', stoptol(ilvl+1), 'auxfun', @(x)dualbetafun(x,IRTstruct,logpostfun,beta(ilvl+1),beta(ilvl),lFshift, vec, reference, IRTdenom(ilvl+1)));
            
        case 'build_ftt'
            debug_size = max(y0(:,ilvl+1)) + kickrank(ilvl+1)*nswp(ilvl+1);
            debug_x = zeros(d, debug_size);
            for k = 1:d
                debug_x(k,:) = sample_oned_domain(x{k}, debug_size);
            end
            opts = ftt_options('method', 'AMEN', 'ng_flag', false, 'oned_ref', x, ...
                'err_tol', stoptol(ilvl+1), 'loc_err_tol', trunctol(ilvl+1), 'max_als', nswp(ilvl+1), ...
                'kick_rank', kickrank(ilvl+1), 'max_rank', 50, 'init_rank', max(y0(:,ilvl+1)));
            F{ilvl} = build_ftt(@(x)dualbetafun(x',IRTstruct,logpostfun,beta(ilvl+1),beta(ilvl),lFshift, vec, reference, IRTdenom(ilvl+1))', d, [], opts, 'sample_x', debug_x, 'debug_x', debug_x);
            F{ilvl}.ng_flag = true; % set this manually after constructing the sqrt(function) explicitly
            evalcnt1 = nan;
            F{ilvl} = build_irt(F{ilvl});
    end
    
    IRTstruct.evalcnt(ilvl+1) = IRTstruct.evalcnt(ilvl+1) + sum(evalcnt1); % Record the number of evaluations
    
    if (isa(F{ilvl}, 'tt_tensor'))
        if (plotdiag)
            % draw 1D marginals
            Fdiag = zeros(max(F{ilvl}.n), d);
            Fdiag(1:F{ilvl}.n(1),1) = full(dot(tt_ones(F{ilvl}.n(2:d)), F{ilvl}, 2, d));
            for i=2:d-1
                Fdiag(1:F{ilvl}.n(i),i) = full(dot(tt_ones(F{ilvl}.n(i+1:d)), dot(tt_ones(F{ilvl}.n(1:i-1)), F{ilvl}, 1, i-1), 2, d-i+1));
            end
            Fdiag(1:F{ilvl}.n(d),d) = full(dot(tt_ones(F{ilvl}.n(1:d-1)), F{ilvl}, 1, d-1));
            figure(1);
            plot(Fdiag);
            legend toggle;
            title('1D marginal sqrt(densities)');
            
            % draw 2D marginal
            figure(2);
            if (d==2)
                surf(full(F{ilvl}, F{ilvl}.n(1:2)'), 'EdgeColor', 'none');
            else
                surf(full(dot(tt_ones(F{ilvl}.n(3:end))/prod(F{ilvl}.n(3:end)), F{ilvl}, 3, d), F{ilvl}.n(1:2)'), 'EdgeColor', 'none');
            end
            shading interp;
            title('2D u_1 u_2 marginal')
            drawnow;
        end
        
        % Initial guess for the next step
        if (size(y0,2)<(ilvl+2)) % Continue with the last prescribed TT rank gracefully
            y0(:,ilvl+2) = y0(:,end);                                      %#ok
        end
        Fprev = round(F{ilvl}, 0, y0(:,ilvl+2)); % Initial guess with rank y0
        IRTstruct.Fprev = Fprev;
        
        % disintegrate into cells for faster sampling
        F{ilvl} = core2cell(F{ilvl});
    end
    
    % Record the current DIRT stack
    IRTstruct.F = F(1:ilvl);
    IRTstruct.beta = beta(1:ilvl+1);
    
    if (testsamples>0)
        % Test approximation
        y = randref(reference,  min(sum(evalcnt1),testsamples),d);
        [y,lFapp,lFex] = tt_dirt_sample(IRTstruct, y, @(x)logpostfun(x,0,beta(ilvl+1)), vec);
        [y2,~,~,num_of_rejects] = mcmc_prune(y, lFex, lFapp);
        num_of_rejects = num_of_rejects*100/size(y,1);
        tau = essinv(lFex, lFapp);
        hl = hellinger(lFex, lFapp);
        fprintf('N/ESS = %g, Hellinger = %3.3e\n\n', tau, hl);
        if (plotdiag)
            figure(3);
            plot(y2);
            title(sprintf('Chain: #rejects = %g%%, N/ESS = %g, H = %3.3e', num_of_rejects, tau, hl));
            drawnow;
        end
        IRTstruct.evalcnt(ilvl+1) = IRTstruct.evalcnt(ilvl+1) + size(y,1);
        if (tau > recompute)
            % Recompute the current level
            ilvl = ilvl-1;
            recompute_count = recompute_count + 1;
            if (recompute_count>4)
                error('Too poor approximation at beta=%g after 5 attempts, giving up', beta(ilvl+2));
            end
        else
            % Update the baseline of log-density
            if (ilvl<nlvl)
                if (IRTdenom(ilvl+1))
                    lFshift = max(lFex)*beta(ilvl+2)/beta(ilvl+1) - max(lFapp);
                else
                    lFshift = max(lFex)*(beta(ilvl+2)-beta(ilvl+1))/beta(ilvl+1);
                end
                IRTstruct.lFshift = lFshift;
            end
            recompute_count = 0;
        end
    end
    
    ilvl = ilvl+1;
end
end




% Helper function for faster evaluation of importance ratios
function [F]=dualbetafun(x,IRTstruct,logpostfun,beta,beta_prev,lFshift, vec, reference, IRTdenom)
% Sample existing DIRT
[z,lFapp] = tt_dirt_sample(IRTstruct, x);
if (IRTdenom) % Compute ratio with the IRT density in the denominator
    beta_prev = 0;
end
% Exact log-ratio
F = logpostfun_vec(z,beta_prev, beta, logpostfun, vec);
F = F-lFshift; % Remove baseline to prevent overflow
if (IRTdenom)
    F = F - lFapp;
end
if (reference(1)~='u')
    F = F - sum(x.^2,2)/2; % Remove the reference log-density
end
F = exp(F*0.5); % SQRT(Density)
end


% A wrapper for vectorised (or not) user function
function [y]=logpostfun_vec(x, beta_min, beta_max, logpostfun, vec)
if (vec)
    y = logpostfun(x,beta_min,beta_max);
else
    y = zeros(size(x,1), 1);
    if isempty(gcp('nocreate'))
        for i=1:size(x,1)
            y(i) = logpostfun(x(i,:),beta_min, beta_max);
        end
    else  % if there is a parallel pool, use it
        parfor i=1:size(x,1)
            y(i) = logpostfun(x(i,:),beta_min, beta_max);
            % no idea what's wrong about broadcasting logpostfun...
        end
    end
end
end

