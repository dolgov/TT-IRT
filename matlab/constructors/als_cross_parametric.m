function [u,time_extern,funevals]=als_cross_parametric(coeff, assem_solve_fun, tol, varargin)
% TT ALS-Cross algorithm.
%   [u,time_extern,funevals]=ALS_CROSS_PARAMETRIC(coeff, assem_solve_fun, tol, varargin)
%
% Inputs:
%   coeff: (d+1)-dimensional block TT format storing coefficients, with the
%       first TT rank Mc corresponding to different coefficients, and the 
%       entire first TT block corresponding to the deterministic part with 
%       Nxc degrees of freedom in the deterministic variable. The other 
%       TT blocks correspond to the d (stochastic) parameters.
%   assem_solve_fun: a handle to the function which should take the 
%       coefficients Ci of size [Mc, Nxc, r] and return cell arrays U,A,F, 
%       respectively solutions, matrices and RHS at the given Ci. 
%       Here r is the TT rank of the solution in the first TT block.
%       U,A,F must be cell arrays of size 1 x r, each snapshot U{i} must be
%       a column vector of size Nxu x 1, each snapshot F{i} must be a 
%       column vector of size Nxa x 1, and each matrix A{i} must be a
%       Nxa x Nxa matrix.
%       Alternatively, if the 'funarg' option (see below) is switched to 
%       'indices', assem_solve_fun should take indices of parameters where
%       the systems must be solved, in the form of r x d integer matrix,
%       with elements in the k-th column ranging from 1 to n_k, the number
%       of grid points in the k-th parameter. The output format is the same
%       !!! In both cases, A and F should depend near linearly on coeff !!!
%   tol: cross truncation and stopping tolerance
%
% Optional inputs given in varargin:
%   Pua: a matrix to project spatial block of solution to spat. block of matrix
%       For good performance, Pua should be a full rank projector, with
%       size(Pua)==[Nxa,Nxu] with Nxu>=Nxa.
%       Default empty Pua assumes Nxu==Nxa.
%   nswp: max number of iterations (default 5)
%   kickrank: max TT rank of the residual/enrichment (default 10)
%   random_init: if greater than 0, take random_init random indices at
%       start; if 0 (default), take maxvol indices of coeff
%   funarg: selects the type of input for assem_solve_fun:
%           'values' (default) assumes that the function takes values of
%                    the coefficients,
%           'indices' assumes that the function takes indices of parameters
%
% Outputs:
%   u: the solution in the TT format
%   time_extern: for profiling: a vector of cpu times for:
%       [solving determ. systems, projecting determ. systems]
%   funevals: total number of deterministic solves

nswp = 5;
kickrank = 10; % The base for enrichment rank. The actual ranks are scaled to the coeff. ranks
Pua = []; % A matrix that maps spatial DOFS of the solution to the spatial DOFs of the matrix/rhs
random_init = 0; % If >0, start with random indices, instead of the coefficient
funarg = 'values'; % whether assem_solve_fun takes values or indices

% Parse parameters
i = 1;
vars = varargin;
while (i<length(vars))
    switch lower(vars{i})
        case 'nswp'
            nswp=vars{i+1};
        case 'kickrank'
            kickrank=vars{i+1};
        case 'pua'
            Pua = vars{i+1};
        case 'random_init'
            random_init = vars{i+1};
        case 'funarg'
            funarg = lower(vars{i+1});
        otherwise
            warning('Option %s was not recognized', vars{i});
    end
    i=i+2;
end

% All grid sizes
ny = coeff.n;
d = numel(ny);
Nxc = ny(1); % Spatial grid size of the coefficient
ny = ny(2:d); % Parametric grid sizes
d = d-1; % Consider d to be param. dimension only
rc = coeff.r; % TT ranks of the coefficient
Mc = rc(1); % number of coefficient components (can be [coeff rhs], etc.)
rc = rc(2:d+2);
ru = rc; % these will be TT ranks of solution
coeff = core2cell(coeff);
C0 = coeff{1}; % spatial coeff block
coeff = coeff(2:d+1); % parametric TT part

% Prepare storage for reduction/sampling matrices
if (nswp>1)
    UAU = cell(d+1,1); % left Galerkin reductions, matrix
    UF = cell(d+1,1); % left Galerkin reductions, RHS
    % Otherwise we don't need to store all left reductions
end
UC = cell(d+1,1); % right cross samples of C on U-indices
UC{d+1} = 1;
% Prepare storage for the residual
if (kickrank>0)
    ZZ = cell(d+1,1);
    ZZ{d+1} = 1;
    ZU = ZZ; % ZAU from the left, ZU from the right
    ZC = ZZ;
    rz = round(kickrank.*rc/max(rc));
    rz(rz<1)=1;
    rz(d+1) = 1;
end

xi = ones(1,random_init);
if (strcmp(funarg, 'indices'))
    % Initialise global indices if the user function works with them
    Ju = [];
end

% First, orthogonalize the coefficient.
% We can derive its optimal indices (as an initial guess), or use random
% ones
v = 1;
for i=d:-1:1
    crc = reshape(coeff{i}, rc(i)*ny(i), []);
    crc = crc*v.';
    crc = reshape(crc, rc(i), ny(i)*rc(i+1));
    crc = crc.';
    [crc,v]=qr(crc, 0);
    v = v(:,1:rc(i));
    rc(i) = size(crc,2);
    crc = crc.';
    ind = maxvol2(crc.');
    CC = crc(:,ind);
    crc = CC\crc;
    v = CC.'*v;
    coeff{i} = reshape(crc, rc(i), ny(i), rc(i+1));
    
    if (strcmp(funarg, 'indices'))
        % Collect global parameter indices for the deterministic solver
        Ju = [repmat((1:ny(i))', rc(i+1),1), reshape(repmat(Ju(:)', ny(i), 1), ny(i)*rc(i+1), d-i)];
        Ju = Ju(ind, :);
    end
    
    if (random_init)&&(i>1)
        % Random-sample coeff from the right
        ind = rand(1,random_init);
        ind = ind*(ny(i)-1);
        ind = round(ind)+1; % should be >=1 and <=ny now
        for k=1:random_init
            xi(1:rc(i),k) = reshape(coeff{i}(:,ind(k),:), rc(i), rc(i+1))*xi(1:rc(i+1),k);
        end
        xi = xi(1:rc(i),:);
        UC{i} = xi; % this is of size rc(i) x nq always
        ru(i) = random_init;
    else
        UC{i} = eye(rc(i));
        ru(i) = rc(i);
    end
    
    if (kickrank>0)
        % Initialize the residual
        crz = randn(ny(i)*rz(i+1), rz(i));
        [crz,~]=qr(crz,0);
        crz = crz.';
        rz(i) = size(crz,1);
        ind = maxvol2(crz.');
        % Sample the coefficient and solution at Z-indices
        ZC{i} = reshape(coeff{i}, rc(i)*ny(i), rc(i+1))*ZC{i+1};
        ZC{i} = reshape(ZC{i}, rc(i), ny(i)*rz(i+1));
        ZC{i} = ZC{i}(:,ind);
        ZU{i} = ZC{i}; % in the first iteration this is the same
    end
end
% Init empty solution storage
u = cell(d,1);
U0 = [];

% The coefficient+rhs at sampled indices
C0 = reshape(C0, Mc*Nxc, []); % size Mc*Nxc, rc1
C0 = C0*v.';
% This is the spatial block of the coefficient, in the representation when
% all parametric blocks contain identity matrices.
% The coeff will not change anymore
C0 = reshape(C0, Mc, Nxc, rc(1));

% Initialise cost profilers
time_extern = [0, 0];  % 1 - solve, 2 - project
funevals = 0; % init number of det. solves

% Now start iterating
i = 0;
max_dx = 0;
dir = 1;
swp = 1;
while (swp<=nswp)
    if (i==0)
        %%%%%%%%%%%%%%%%%% Work on the spatial block
        % Previous guess (if we have one)
        Uprev = U0;
        % Solve deterministic problems at the U-indices
        if (strcmp(funarg, 'indices'))
            Ci = Ju;
        else
            % Construct the coeff there
            Ci = reshape(C0, Mc*Nxc, [])*UC{1}; % size Mc*Nxc, rc1
            Ci = reshape(Ci, Mc, Nxc, ru(1));
        end
        t1__uc = tic;
        if (swp==1)
            [U0,A0s,F0] = assem_solve_fun(Ci);
            Nxa = size(A0s{1},1); % number of spatial DOFS in A (can differ from Nxc!)
            F0 = cell2mat(F0);
            % In the first sweep, Ci==C0, and we need the corresponding
            % matrices (A0s) and RHS (F0), since we'll use them in
            % subsequent iterations ...
        else
            % ... where it's enough to compute the solution only
            U0 = assem_solve_fun(Ci);
        end
        time_extern(1) = time_extern(1) + toc(t1__uc);
        clear Ci; % save some mem
        funevals = funevals + ru(1);
        U0 = cell2mat(U0);
        Nxu = size(U0,1); % this, again, can differ from Nxa or Nxc
        if (Nxu~=Nxa)&&(isempty(Pua))
            error('Numbers of spatial DOFs in u and A differ, and no transformation matrix is given. Unable to reduce model');
        end
        
        % Check the error
        dx = 1;
        if (~isempty(Uprev))
            dx = norm(U0-Uprev,'fro')/norm(U0,'fro');
        end
        max_dx = max(max_dx,dx);
        
        fprintf('=als_cross_parametric= 0 swp=%d, max_dx=%3.3e, max_rank=%d\n', swp, max_dx, max(ru));
        
        % Unusual ALS exit condition: after solving for the spatial block,
        % to have the non-orth center here.
        if ((max_dx<tol)||(swp>nswp))
            break;
        end
        max_dx = 0;
        
        % Truncate U0 via full-pivoted cross approx
        [U0,v]=localcross(U0,tol/sqrt(d));
        if (swp>1)
            u{1} = reshape(u{1}, ru(1), ny(1)*ru(2));
            u{1} = v*u{1};
        end
        ru(1) = size(U0,2);
        
        if (kickrank>0)
            % Compute the residual
            % Compute A*U at Z indices
            cru = U0*v*ZU{1};
            if (Nxa~=Nxu)
                cru = Pua*cru; % size Nxa x rz1
            end
            Z0 = zeros(Nxa, rz(1));
            for j=1:rz(1)
                crA = A0s{1}*ZC{1}(1,j);
                for k=2:rc(1)
                    crA = crA + A0s{k}*ZC{1}(k,j);
                end
                Z0(:,j) = crA*cru(:,j);
            end
            Z0 = Z0 - F0*ZC{1};
            % QR it
            [Z0,~]=qr(Z0,0);
            rz(1) = size(Z0,2);
            if (Nxa~=Nxu)
                cru = [U0, Pua'*Z0]; % make it Nxu x rz(1)
            else
                cru = [U0, Z0];
            end
            % QR U0
            [U0,v]=qr(cru,0);
            v = v(:,1:ru(1));
            if (swp>1)
                u{1} = reshape(u{1}, ru(1), ny(1)*ru(2));
                u{1} = v*u{1};
            end
            ru(1) = size(U0,2);
        end
        
        % Project the model onto the solution basis U0
        % UAU is very large here
        % Need to run some loops to save mem
        t1__uc = tic;
        UAU_new = cell(1,rc(1));
        Uprev = U0;
        if (Nxa~=Nxu)
            Uprev = Pua*U0; % Project U0 to the size of A
        end
        for j=1:rc(1)
            UAU_new{j} = A0s{j}*Uprev;
            UAU_new{j} = Uprev'*UAU_new{j};
            UAU_new{j} = reshape(UAU_new{j}, [], 1);
        end
        if (nswp==1)
            % we don't need to save UAU projections for all blocks, if we
            % will never iterate back
            UAU = cell2mat(UAU_new);
            UF = Uprev'*F0; % size ru', rc
        else
            UAU{1} = cell2mat(UAU_new);
            % Project the RHS block
            UF{1} = Uprev'*F0; % size ru', rc
        end
        time_extern(2) = time_extern(2) + toc(t1__uc);
        clear UAU_new
        % Project onto the residual
        if (kickrank>0)
            % Project onto resid basis
            ZU_new = cell(1,rc(1));
            for j=1:rc(1)
                ZU_new{j} = Z0'*A0s{j};
                ZU_new{j} = ZU_new{j}*Uprev;
                ZU_new{j} = reshape(ZU_new{j}, [], 1);
            end
            ZU{1} = cell2mat(ZU_new);
            ZC{1} = Z0'*F0;
            clear ZU_new
        end
        % Save some mem if we only have 1 iteration
        if (nswp==1)
            clear A0s;
            clear F0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%% Done with i==0
    end
    
    % Loop for reduced system
    if (i>0)
        % Solve block-diag reduced system
        crC = reshape(coeff{i}, rc(i)*ny(i), rc(i+1));
        crC = crC*UC{i+1};  % now these are indices from the right, and Galerkin from the left
        crC = reshape(crC, rc(i), ny(i)*ru(i+1)); % complexity rc1*n x rc2 x ru2
        if (nswp==1)
            UAUi = UAU; % cell-free local reduction from the left
            UFi = UF;
        else
            UAUi = UAU{i};
            UFi = UF{i};
        end
        crF = UFi*crC;
        try
            % Use faster MEX solver if it's available
            crC = reshape(crC, rc(i), ny(i), ru(i+1));
            UAUi = reshape(UAUi, ru(i), ru(i), rc(i));
            cru = solve_blockdiag_mex(UAUi, crC, crF);
        catch ME
            fprintf(ME.message);  fprintf('\n'); % no MEX file
            % Solve the system via Matlab loop
            crA = reshape(UAUi, ru(i)*ru(i), rc(i));
            cru = cell(1, ny(i)*ru(i+1));
            crC = reshape(crC, rc(i), ny(i)*ru(i+1));
            crC = num2cell(crC,1);
            crF = num2cell(crF,1);
            for j=1:ny(i)*ru(i+1)
                Ai = crA*crC{j}; % ru1*ru1 x rc1 x 1
                Ai = reshape(Ai, ru(i), ru(i));
                cru{j} = Ai\crF{j};
            end
            cru = cell2mat(cru);
        end
        
        cru = reshape(cru, ru(i)*ny(i), ru(i+1));
        % Error check
        dx = 1;
        if (~isempty(u{i}))
            dx = norm(cru - reshape(u{i}, ru(i)*ny(i), ru(i+1)), 'fro')/norm(cru,'fro');
        end
        max_dx = max(max_dx, dx);
        u{i} = reshape(cru, ru(i), ny(i), ru(i+1)); % if it's d-th block, we don't do anything else
        
        if (i<d)&&(dir>0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Left-Right sweep
            % Truncate cru with cross
            [cru,v]=localcross(cru,tol/sqrt(d));
            rv = v;
            if (kickrank>0)
                % Update the residual and enrichment
                crC = reshape(coeff{i}, rc(i)*ny(i), rc(i+1));
                crC = crC*ZC{i+1};  % now these are indices from the right, and Galerkin from the left
                crC = reshape(crC, rc(i), ny(i)*rz(i+1));
                crC = num2cell(crC,1);
                Uprev = cru*rv*ZU{i+1};
                Uprev = reshape(Uprev, ru(i), ny(i)*rz(i+1));
                Uprev = num2cell(Uprev,1);
                % First, enrichment
                crA = reshape(UAUi, ru(i)*ru(i), rc(i));
                crz = cell(1, ny(i)*rz(i+1));
                for j=1:ny(i)*rz(i+1)
                    Ai = crA*crC{j};
                    Ai = reshape(Ai, ru(i), ru(i));
                    crz{j} = Ai*Uprev{j};
                end
                crz = cell2mat(crz); % size ru1 x n*rz2
                crz = crz - UFi*cell2mat(crC);
                crz = reshape(crz, ru(i)*ny(i), rz(i+1));
                v = [cru, crz];
                cru = v;
                % QR u
                [cru,v]=qr(cru,0);
                v = v(:,1:size(rv,1));
                rv = v*rv;
                % Now the residual itself
                crA = reshape(ZU{i}, rz(i)*ru(i), rc(i));
                crz = cell(1, ny(i)*rz(i+1));
                for j=1:ny(i)*rz(i+1)
                    Ai = crA*crC{j};
                    Ai = reshape(Ai, rz(i), ru(i));
                    crz{j} = Ai*Uprev{j};
                end
                crz = cell2mat(crz); % size rz1 x n*rz2
                crz = crz - ZC{i}*cell2mat(crC);
                crz = reshape(crz, rz(i)*ny(i), rz(i+1));
            end
            % cast the non-orth factor to the next block
            if (swp>1)
                u{i+1} = reshape(u{i+1}, ru(i+1), ny(i+1)*ru(i+2));
                u{i+1} = rv*u{i+1};
                u{i+1} = reshape(u{i+1}, [], ny(i+1), ru(i+2));
            end
            ru(i+1) = size(cru, 2);
            u{i} = reshape(cru, ru(i), ny(i), ru(i+1));
            
            % Projection from the left -- Galerkin
            % With matrix
            crC = reshape(coeff{i}, rc(i), ny(i), rc(i+1));
            cru = reshape(cru, ru(i), ny(i), ru(i+1));
            try
                % Faster MEX projection
                UAUi = reshape(UAUi, ru(i), ru(i), rc(i));
                UAU_new = project_blockdiag_mex(UAUi,crC,cru);
            catch ME
                fprintf(ME.message); fprintf('\n'); % MEX is not available or damaged
                UAU_new = zeros(ru(i+1), ru(i+1)*rc(i+1));
                UAUi = reshape(UAUi, ru(i), ru(i)*rc(i));
                crC = permute(crC, [1,3,2]);
                cru = permute(cru, [1,3,2]);
                for j=1:ny(i)
                    v = cru(:,:,j);
                    crA = v'*UAUi; % ru2 x ru1 x ru1*rc1
                    crA = reshape(crA, ru(i+1)*ru(i), rc(i)); % ru2',ru1,rc1
                    crA = crA*crC(:,:,j); % ru2*ru1 x rc1 x rc2
                    crA = reshape(crA, ru(i+1), ru(i)*rc(i+1)); % ru2',ru1,rc2
                    crA = crA.';
                    crA = reshape(crA, ru(i), rc(i+1)*ru(i+1)); % ru1,rc2,ru2'
                    crA = v.'*crA; % ru2 x ru1 x rc2*ru2
                    crA = reshape(crA, ru(i+1)*rc(i+1), ru(i+1)); % ru2,rc2,ru2'
                    crA = crA.';
                    UAU_new = UAU_new+crA;
                end
                cru = permute(cru, [1,3,2]);
            end
            % Save some mem
            if (nswp==1)
                UAU = UAU_new;
            else
                UAU{i+1} = UAU_new;
            end
            clear UAU_new
            clear UAUi
            
            % With RHS
            crC = reshape(coeff{i}, rc(i), ny(i)*rc(i+1));
            UFi = UFi*crC;
            UFi = reshape(UFi, ru(i)*ny(i), rc(i+1));
            cru = reshape(cru, ru(i)*ny(i), ru(i+1));
            UFi = cru'*UFi;
            if (nswp==1)
                UF = UFi;
            else
                UF{i+1} = UFi;
            end
            clear UFi;
            
            % Projections with Z
            if (kickrank>0)
                [crz,~]=qr(crz,0);
                rz(i+1) = size(crz,2);
                % With matrix
                crC = reshape(crC, rc(i), ny(i), rc(i+1));
                cru = reshape(cru, ru(i), ny(i), ru(i+1));
                crz = reshape(crz, rz(i), ny(i), rz(i+1));
                ZU{i+1} = zeros(rz(i+1), ru(i+1)*rc(i+1));
                ZU{i} = reshape(ZU{i}, rz(i), ru(i)*rc(i));
                crC = permute(crC, [1,3,2]);
                cru = permute(cru, [1,3,2]);
                crz = permute(crz, [1,3,2]); % size rz1,rz2,n
                for j=1:ny(i)
                    v = crz(:,:,j);
                    crA = v'*ZU{i}; % rz2 x rz1 x ru1*rc1
                    crA = reshape(crA, rz(i+1)*ru(i), rc(i)); % ru2',ru1,rc1
                    crA = crA*crC(:,:,j); % ru2*ru1 x rc1 x rc2
                    crA = reshape(crA, rz(i+1), ru(i)*rc(i+1)); % ru2',ru1,rc2
                    crA = crA.';
                    crA = reshape(crA, ru(i), rc(i+1)*rz(i+1)); % ru1,rc2,ru2'
                    v = cru(:,:,j);
                    crA = v.'*crA; % ru2 x ru1 x rc2*ru2
                    crA = reshape(crA, ru(i+1)*rc(i+1), rz(i+1)); % ru2,rc2,ru2'
                    crA = crA.';
                    ZU{i+1} = ZU{i+1}+crA;
                end
                crz = permute(crz, [1,3,2]);
                crz = reshape(crz, rz(i)*ny(i), rz(i+1));
                % With RHS
                crC = reshape(coeff{i}, rc(i), ny(i)*rc(i+1));
                ZC{i+1} = ZC{i}*crC;
                ZC{i+1} = reshape(ZC{i+1}, rz(i)*ny(i), rc(i+1));
                ZC{i+1} = crz'*ZC{i+1};
                if (nswp==1)
                    ZC{i} = []; % clear old reductions unused anymore
                    ZU{i} = [];
                end
            end
            
        elseif (dir<0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right-Left sweep
            cru = reshape(cru, ru(i), ny(i)*ru(i+1));
            % Truncate cru with cross
            [v,cru]=localcross(cru,tol/sqrt(d));
            % now cru is not orthogonal
            rv = v;
            if (kickrank>0)
                % Update the residual and enrichment
                % First, enrichment
                crC = reshape(coeff{i}, rc(i)*ny(i), rc(i+1));
                crC = crC*UC{i+1};  % now these are indices from the right, and Galerkin from the left
                crC = reshape(crC, rc(i), ny(i)*ru(i+1));
                crC = num2cell(crC,1);
                Uprev = reshape(rv*cru, ru(i)*ny(i), ru(i+1));
                Uprev = reshape(Uprev, ru(i), ny(i)*ru(i+1));
                Uprev = num2cell(Uprev,1);
                crA = reshape(ZU{i}, rz(i)*ru(i), rc(i));
                crz = cell(1, ny(i)*ru(i+1));
                for j=1:ny(i)*ru(i+1)
                    Ai = crA*crC{j};
                    Ai = reshape(Ai, rz(i), ru(i));
                    crz{j} = Ai*Uprev{j};
                end
                crz = cell2mat(crz); % size rz1 x n*ru2
                crz = crz - ZC{i}*cell2mat(crC);
                crz = reshape(crz, rz(i)*ny(i), ru(i+1));
                crz = reshape(crz, rz(i), ny(i)*ru(i+1));
                v = [cru; crz];
                % Now the residual itself
                crC = reshape(coeff{i}, rc(i)*ny(i), rc(i+1));
                crC = crC*ZC{i+1};  % now these are indices from the right, and Galerkin from the left
                crC = reshape(crC, rc(i), ny(i)*rz(i+1));
                crC = num2cell(crC,1);
                Uprev = reshape(rv*cru, ru(i)*ny(i), ru(i+1));
                Uprev = Uprev*ZU{i+1};
                Uprev = reshape(Uprev, ru(i), ny(i)*rz(i+1));
                Uprev = num2cell(Uprev,1);
                crz = cell(1, ny(i)*rz(i+1));
                for j=1:ny(i)*rz(i+1)
                    Ai = crA*crC{j};
                    Ai = reshape(Ai, rz(i), ru(i));
                    crz{j} = Ai*Uprev{j};
                end
                crz = cell2mat(crz); % size rz1 x n*rz2
                crz = crz - ZC{i}*cell2mat(crC);
                crz = reshape(crz, rz(i)*ny(i), rz(i+1));
                cru = v;
            end
            ru(i) = size(rv,2);
            % QR u
            [cru,v]=qr(cru.', 0);
            v = v(:,1:ru(i));
            rv = rv*v.';
            cru = cru.';
            % Maxvol to determine new indices
            ind = maxvol2(cru.');
            UU = cru(:,ind);
            cru = UU\cru;
            rv = rv*UU;
            ru(i) = size(rv,2);
            % Cast non-orth factor to the next block
            if (i>1)
                u{i-1} = reshape(u{i-1}, ru(i-1)*ny(i-1), []);
                u{i-1} = u{i-1}*rv; % it's now Mc+1 x r x n x r
                u{i-1} = reshape(u{i-1}, ru(i-1), ny(i-1), []);
            else
                U0 = U0*rv;
            end
            ru(i) = size(cru,1);
            u{i} = reshape(cru, ru(i), ny(i), ru(i+1));
            
            if (strcmp(funarg, 'indices'))
                Ju = [repmat((1:ny(i))', ru(i+1),1), reshape(repmat(Ju(:)', ny(i), 1), ny(i)*ru(i+1), d-i)];
                Ju = Ju(ind, :);
            end            
            
            % Projection from the right -- sample C on U indices
            UC{i} = reshape(coeff{i}, rc(i)*ny(i), rc(i+1));
            UC{i} = UC{i}*UC{i+1};
            UC{i} = reshape(UC{i}, rc(i), ny(i)*ru(i+1));
            UC{i} = UC{i}(:, ind);
            % Reductions with Z
            if (kickrank>0)
                % QR and maxvol Z
                crz = reshape(crz, rz(i), ny(i)*rz(i+1));
                [crz,~]=qr(crz.', 0);
                rz(i) = size(crz,2);
                ind = maxvol2(crz);
                % Sample C and U
                ZC{i} = reshape(coeff{i}, rc(i)*ny(i), rc(i+1));
                ZC{i} = ZC{i}*ZC{i+1};
                ZC{i} = reshape(ZC{i}, rc(i), ny(i)*rz(i+1));
                ZC{i} = ZC{i}(:, ind);
                ZU{i} = reshape(u{i}, ru(i)*ny(i), ru(i+1));
                ZU{i} = ZU{i}*ZU{i+1};
                ZU{i} = reshape(ZU{i}, ru(i), ny(i)*rz(i+1));
                ZU{i} = ZU{i}(:, ind);
            end
        end
        
        fprintf('=als_cross_parametric= swp=%d (%d), i=%d, dx=%3.3e, rank=[%d,%d]\n', swp, dir, i, dx, ru(i),ru(i+1));
    end
    
    i = i+dir;
    
    if (dir>0)&&(i==d+1)&&(swp==nswp)
        break; % Last block & last sweep
    end
    if (dir>0)&&(i==d)&&(swp<nswp)
        % Turn at the right end
        fprintf('=als_cross_parametric= fwd swp=%d, max_dx=%3.3e, max_rank=%d\n', swp, max_dx, max(ru));
        dir = -dir;
        swp = swp+1;
        max_dx = 0;
        if (strcmp(funarg, 'indices'))
            Ju = [];
        end
    end
    if (i==0)&&(dir<0)
        % turn at the left end
        dir = -dir;
        swp = swp+1;
    end
end

U0 = reshape(U0, 1, Nxu, ru(1));
u = [{U0}; u];
u = cell2core(tt_tensor, u);
end
