function [y,statvals,statind,Jy,evalcnt]=amen_cross_s(inp, fun, tol, varargin)
%   [y,statvals,statind,Jy,evalcnt]=AMEN_CROSS_S(inp, fun, tol, varargin)
% Block cross with error-based enrichment ("S": stabilized, statistics).
% Tries to interpolate the function(s) via the error-enriched maxvol-cross.
% Allows both index-defined (traditional) functions
% and elementwise functions depending on other tt_tensors,
% computation of min,max values, and returning of maxvol indices.
%
% If inp==n is a column vector of mode sizes:
%       fun = @(ind)fun(ind) is a sought function of index ind,
%       ind comes as an array of sizes M x d,
%       fun should return the block M x B (vectorized variant).
%       M=1 if vec=false, hence the return may be either 1 x B or B x 1.
% If inp=={x1,x2,...,xN} a cell array of tt_tensors x1,x2,...xN:
%       fun = @(x)fun(x) is a sought function of elements of x=[x1,x2,...].
%       it should receive a 2d array V of sizes M x N, where the
%       first dimension stands for the reduced set of spatial indices, and the
%       second dimension is the enumerator of vectors in X.
%       The returned sizes should be M x B, where B is the number of
%       components in FUNS.
% In addition to these obligatory three inputs, the second function of another
% type may be added via optional parameters 'auxinp', 'auxfun' (see below).
%
% Optional arguments are provided in the form
% 'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so on.
% The list of option names and default values:
%       o y0 - initial approximation. Can be either a tt_tensor, or an
%              integer value. In the latter case, y0 random right indices
%              are used in the first sweep. [4]
%       o nswp - maximal number of sweeps [20]
%       o stop_sweep - number of extra iterations carried out after err<tol [0]
%       o kickrank - error/enrichment rank. If integer, prescribes the
%              exact rank of the error. If between 0 and 1, prescribes the
%              error ranks as a fraction of the solution ranks, i.e.
%              rz = kickrank*ry.  [4]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o vec - whether funs can accept and return vectorized values [true]
%       o tol_exit - stopping tolerance [tol]
%       o exitdir - if 1, return after the forward sweep, if -1, return the
%                   backward sweep, if 0, after any [0]
%       o auxinp - secondary input data
%       o auxfun - secondary input function
%
% Moreover, it's possible to estimate minimal and maximal values in
% modulus, real or imaginary part. Pass an arbitrary number of varargins
% from the set {'sr', 'lr', 'sm', 'lm', 'si', 'li'}, which is consistent
% with eigs notation. The output statvals contains the corresponding values
% of the EXACT y (NOT the TT decomposition of y, contrary to tt_stat),
% and statind contains the corresponding indices (optimizers).
%
% Jy contains the final indices. Their direction depends on the
% direction of the last sweep, use exitdir parameter to control it.
%
% evalcnt is the total number of evaluations, returned as a 1 x 2 array of
% numbers of @(ind)fun(ind) and @(x)fun(x) evaluations, respectively.
%
%********
%   Alternating optimization with enrichment (AMEn):
%   S. Dolgov, D. Savostyanov. SIAM J. Sci. Comput., 36 (2014), p. A2248.
%
%   Please send feedback to: {sergey.v.dolgov,dmitry.savostyanov}@gmail.com
%
%---------------------------
%
% Example 1: compare with amen_cross on 1/r
% (1) y = amen_cross(192*ones(50,1), @(i)sqrt(1./sum(i.^2,2)), tol, 'vec', true, 'kickrank', 4, 'zrank', 4, 'zrank2', 0, 'exitdir', -1);
% (2) y2 = amen_cross_s(192*ones(50,1), @(i)sqrt(1./sum(i.^2,2)), tol, 'exitdir', -1);
% tol    Time(1)      error(1)     Time(2)      error(2)
% 1e-2   2.436728     2.9794e-01   2.829848     1.0947e-02
% 1e-4   5.586820     1.8110e-03   4.442022     7.2383e-05
% 1e-6   11.721951    1.0534e-05   6.052327     4.9593e-06
% 1e-8   21.167288    6.5322e-08   13.641011    1.0377e-09
% 1e-10  31.672621    4.1914e-10   29.067563    2.0816e-11
%
% Example 2: compute the Fresnel integral int_0^1000 cos(t^2) via QTT
% t = tt_x(2,30)/2^30*1000;
% y = amen_cross_s({t}, @(t)cos(t.^2), 1e-4, 'kickrank', 0.3, 'nswp', 40);
% fprintf('error: %e\n', dot(y, tt_ones(y.n)/2^30*1000) - sqrt(pi/8));
%

y = 4; % means 4 random indices at start
nswp = 20;
kickrank = 4;
verb = 1;
vec = true;
exitdir=0;
tol_exit = tol;
stop_sweep = 0;

auxinp = [];
auxfun = [];

i = 1;
vars = varargin;
% Distinguish stat QoI from parameters
soughts = cell(1,numel(vars));
while (i<length(vars))
    switch lower(vars{i})
        case 'y0'
            y=vars{i+1};
        case 'nswp'
            nswp=vars{i+1};
        case 'stop_sweep'
            stop_sweep=vars{i+1};
        case 'kickrank'
            kickrank=vars{i+1};
        case 'verb'
            verb = vars{i+1};
        case 'vec'
            vec = vars{i+1};
        case 'auxinp'
            auxinp = vars{i+1};
        case 'auxfun'
            auxfun = vars{i+1};
        case 'exitdir'
            exitdir=vars{i+1};
        case 'tol_exit'
            tol_exit=vars{i+1};
        case {'sr', 'lr', 'sm', 'lm', 'si', 'li'} % Stat params
            soughts{i}=vars{i};
            i=i-1;
        otherwise
            warning('Option %s was not recognized', vars{i});
    end
    i=i+2;
end
% Process the last entry -- it's missed in the loop
if (~isempty(vars))&&(isa(vars{length(vars)},'char'))
    switch lower(vars{length(vars)})
        case {'sr', 'lr', 'sm', 'lm', 'si', 'li'} % Stat params
            soughts{length(vars)}=vars{length(vars)};
    end
end

% Delete unused cells
soughts(cellfun('isempty', soughts))=[];

% We need all those empty guys to call subroutines uniformly. Cf. NULL in C
X = [];
rx = [];
ifun = [];
ffun = [];
n = ones(1,1);

% Distinguish general and TT-fun inputs
if (isa(inp, 'cell'))
    % First input is an array of TT tensors
    X = inp;
    ffun = fun;
elseif (isa(inp, 'double'))
    % First input is mode sizes
    ifun = fun;
    n = inp;
end
if (~isempty(auxinp))&&(~isempty(auxfun))
    if (isa(auxinp, 'cell'))
        % Second input is an array of TT tensors
        if (isempty(ffun))
            X = auxinp;
            ffun = auxfun;
        else
            error('Cannot use ffun on both inputs');
        end
    else
        if (isempty(ifun))
            ifun = auxfun;
            n = auxinp;
        else
            error('Cannot use ifun on both inputs');
        end
    end
end

YX = cell(1);
ZX = cell(1);

% If there is a TT-fun part, prepare it for computations
if (~isempty(X))
    nx = numel(X);
    d = X{1}.d;
    n = X{1}.n;
    rx = zeros(d+1,nx);
    X = reshape(X, 1, nx);
    X = [X; cell(d-1,nx)];
    for i=1:nx
        rx(:,i) = X{1,i}.r;
        X(:,i) = core2cell(X{1,i});
    end
    % Interface matrices
    YX = cell(d+1,nx);
    % these are for stuffing z indices into x
    ZX = cell(d+1,nx);
    for i=1:nx
        YX{1,i}=eye(rx(1,i)); YX{d+1,i}=eye(rx(d+1,i));
        ZX{1,i}=eye(rx(1,i)); ZX{d+1,i}=eye(rx(d+1,i));
    end
end

d = numel(n);

tol_local = tol/sqrt(d);

% Some place to store global indices
Jy = cell(d+1,1);
Jz = cell(d+1,1);

% rank of z -- will be updated
% rz = [1; kickrank*ones(d-1,1); 1];
% projection
ZY = cell(d+1,1);
ZY{1} = 1;
ZY{d+1} = 1;

% Find out what our initial guess is
if (isa(y,'double'))
    % Random samples for the right indices
    if (isscalar(y))
        nq = y;
        Z = rand(nq,d);
        Z = Z.*repmat(n.'-1,nq,1);
        Z = round(Z)+1; % should be >=1 and <=n now
    else
        Z = y;
        nq = size(Z,1);
    end
    
    if (~isempty(X))
        % Interface matrices
        % Sample x from the right
        xi = cell(1,nx);
        % The last block is processed separately
        for j=1:nx
            xi{j} = X{d,j}(:,Z(:,d),:); % size rxd x nq x rxr
            YX{d,j} = reshape(xi{j}, rx(d,j), nq*rx(d+1,j));
        end
        for i=d-1:-1:2
            for j=1:nx
                for k=1:nq
                    xi{j}(1:rx(i,j),k,:) = reshape(X{i,j}(:,Z(k,i),:), rx(i,j), rx(i+1,j))*reshape(xi{j}(1:rx(i+1,j),k,:), rx(i+1,j), rx(d+1,j));
                end
                xi{j} = xi{j}(1:rx(i,j),:,:);
                YX{i,j} = reshape(xi{j}, rx(i,j), nq*rx(d+1,j));
            end
        end
    end
    
    if (~isempty(ifun))||(~isempty(soughts))
        % Store the indices
        for i=d:-1:2
            Jy{i} = Z(:, i:d);
        end
    end
    
    % Initialize y
    y = cell(d,1);
    ry = [1; nq*ones(d-1,1); 1];
    
elseif isa(y,'tt_tensor')
    % Initial guess is TT
    ry = y.r;
    y = core2cell(y);
    % Compute QR and initial indices
    for i=d:-1:2
        [y{i-1},y{i},ry(i),ind] = qrmaxvol_block(y{i-1}, y{i}, -1, []);
        % Store indices if we have ifun
        if (~isempty(ifun))||(~isempty(soughts))
            Jy{i} = indexmerge((1:n(i))', Jy{i+1});
            Jy{i} = Jy{i}(ind, :);
        end
        % Restrict X for ffun
        if (~isempty(X))
            for j=1:nx
                YX{i,j} = reshape(X{i,j}, rx(i,j)*n(i), rx(i+1,j));
                YX{i,j} = YX{i,j}*YX{i+1,j};
                YX{i,j} = reshape(YX{i,j}, rx(i,j), n(i)*ry(i+1), rx(d+1,j));
                YX{i,j} = YX{i,j}(:, ind, :);
                YX{i,j} = reshape(YX{i,j}, rx(i,j), ry(i)*rx(d+1,j));
            end
        end
    end
else
    % Initial guess is a set of right indices
    Jy = y;
    % Initialize y
    y = cell(d,1);
    ry = ones(d+1,1);
    for i=d:-1:2
        ry(i) = size(Jy{i},1);
        % Restrict X for ffun
        if (~isempty(X))
            % Find local indices
            J = indexmerge((1:n(i))', Jy{i+1});
            ind = zeros(ry(i),1);
            for j=1:ry(i)
                indj = find(all(repmat(Jy{i}(j,:), n(i)*ry(i+1),1)==J, 2));
                ind(j) = indj(round(rand(1)*(numel(indj)-1))+1);
            end
            for j=1:nx
                YX{i,j} = reshape(X{i,j}, rx(i,j)*n(i), rx(i+1,j));
                YX{i,j} = YX{i,j}*YX{i+1,j};
                YX{i,j} = reshape(YX{i,j}, rx(i,j), n(i)*ry(i+1), rx(d+1,j));
                YX{i,j} = YX{i,j}(:, ind, :);
                YX{i,j} = reshape(YX{i,j}, rx(i,j), ry(i)*rx(d+1,j));
            end
        end
    end
end

% Generate projections to Z
if (kickrank>0)
    if (abs(kickrank-round(kickrank))<1e-8)
        % kickrank is integer, this is the rank of z
        rz = kickrank;
    else
        % kickrank is a fraction of the rank of y
        rz = ceil(kickrank*max(ry));
    end
    Z = rand(rz,d);
    Z = Z.*repmat(n.'-1,rz,1);
    Z = round(Z)+1; % should be >=1 and <=n now
    
    % If y is already here
    if (~isempty(y{1}))
        xi = ones(1,rz);
        for i=d:-1:2
            for k=1:rz
                xi(1:ry(i),k) = reshape(y{i}(:,Z(k,i),:), ry(i), ry(i+1))*xi(1:ry(i+1),k);
            end
            xi = xi(1:ry(i),:);
            ZY{i} = xi;
        end
    else
        % Otherwise populate with some random
        for i=d:-1:2
            ZY{i} = randn(ry(i), rz);
        end
    end
    
    if (~isempty(X))
        % Interface matrices
        % Sample x from the right
        xi = cell(1,nx);
        for j=1:nx
            xi{j} = X{d,j}(:,Z(:,d),:); % size rxd x nq x rxr
            ZX{d,j} = reshape(xi{j}, rx(d,j), rz*rx(d+1,j));
        end
        for i=d-1:-1:2
            for j=1:nx
                for k=1:rz
                    xi{j}(1:rx(i,j),k,:) = reshape(X{i,j}(:,Z(k,i),:), rx(i,j), rx(i+1,j))*reshape(xi{j}(1:rx(i+1,j),k,:), rx(i+1,j), rx(d+1,j));
                end
                xi{j} = xi{j}(1:rx(i,j),:,:);
                ZX{i,j} = reshape(xi{j}, rx(i,j), rz*rx(d+1,j));
            end
        end
    end
    if (~isempty(ifun))||(~isempty(soughts))
        % Store the indices
        for i=d:-1:2
            Jz{i} = Z(:, i:d);
        end
    end
    
    rz = [1; rz*ones(d-1,1); 1];
end


% Start the loop
swp = 1;
i = 1;
dir = 1;
last_swp = 0;
max_dx = 0;
fevalcnt = 0;
ievalcnt = 0;
while (swp<=nswp)
    if (swp==1)||((dir>0)&&(i>1))||((dir<0)&&(i<d))
        % Evaluate the new core
        [cry,ievalcnt,fevalcnt] = evaluate_fun(i,Jy,Jy,n,ifun,ffun,X,rx,YX,YX,vec,ievalcnt,fevalcnt);
    else
        cry = reshape(y{i}, [], b);
    end
    
    if (~isempty(soughts))
        Jy{i,2} = indexmerge(Jy{i}, (1:n(i))', Jy{i+1});
    end
    
    % block size
    if ((swp==1)&&(i==1))
        b = size(cry, 2);
        % Stat quantities
        statvals = zeros(numel(soughts),1,b);
        statind = zeros(numel(soughts),d,b);
    end
    
    % Check if the user function is sane
    if (numel(cry)~=ry(i)*n(i)*ry(i+1)*b)
        error('%d elements requested, but %g values received. Check your function or use vec=false', ry(i)*n(i)*ry(i+1), numel(cry)/b);
    end
    
    % Stat outputs
    for j=1:numel(soughts)
        switch(lower(soughts{j}))
            case 'lm'
                [val,ind]=max(abs(cry));
                to_update = val>abs(statvals(j,:));
            case 'sm'
                [val,ind]=min(abs(cry));
                to_update = val<abs(statvals(j,:));
            case 'lr'
                [val,ind]=max(real(cry));
                to_update = val>real(statvals(j,:));
            case 'sr'
                [val,ind]=min(real(cry));
                to_update = val<real(statvals(j,:));
            case 'li'
                [val,ind]=max(imag(cry));
                to_update = val>imag(statvals(j,:));
            case 'si'
                [val,ind]=min(imag(cry));
                to_update = val<imag(statvals(j,:));
        end
        % If any new values are better, update them
        if ((swp==1)&&(i==1)) % Of course, in the first step update all
            to_update = true(1,b);
        end
        if (any(to_update))
            statvals(j,1,to_update) = diag(cry(ind(to_update),to_update));
            ind = tt_ind2sub([ry(i), n(i), ry(i+1)], ind(:));
            if (i>1)
                statind(j,1:i-1,to_update) = Jy{i}(ind(to_update,1),:).';
            end
            statind(j,i,to_update) = ind(to_update,2);
            if (i<d)
                statind(j,i+1:d,to_update) = Jy{i+1}(ind(to_update,3), :).';
            end
        end
    end
    
    % Estimate the error -- now in C-norm
    if (isempty(y{i}))
        y{i} = zeros(ry(i)*n(i)*ry(i+1)*b, 1);
    end
    y{i} = reshape(y{i}, ry(i)*n(i)*ry(i+1)*b, 1);
    cry = reshape(cry, ry(i)*n(i)*ry(i+1)*b, 1);
    dx = max(abs(cry-y{i}))/max(abs(cry));
    max_dx = max(max_dx, dx);
    
    % Switch to the next block
    cry = reshape(cry, ry(i), n(i), ry(i+1), b);
    if (dir>0)&&(i<d)
        % Truncation
        ry_old = ry(i+1);
        [y{i},y{i+1},ry(i+1),cry] = truncate_block(cry, y{i+1}, tol_local, dir);
        % Enrichment
        crz = [];
        if (kickrank>0)
            cry = reshape(cry, ry(i)*n(i)*ry_old, b);
            cry = cry.';
            cry = reshape(cry, b*ry(i)*n(i), ry_old);
            cry = cry*ZY{i+1};
            cry = reshape(cry, b, ry(i)*n(i)*rz(i+1));
            cry = cry.';
            % Evaluate the Z-core
            [crz,ievalcnt,fevalcnt] = evaluate_fun(i,Jy,Jz,n,ifun,ffun,X,rx,YX,ZX,vec,ievalcnt,fevalcnt);
            crz = crz - cry;
            crz = reshape(crz, ry(i)*n(i), rz(i+1)*b);
            nrmz = norm(crz, 'fro');
            if (nrmz==0)
                crz = randn(ry(i)*n(i), rz(i+1)*b);
            else
                crz = crz/nrmz; % prevent underflows
            end
            % Truncate cry if it's too large
            if (abs(kickrank-round(kickrank))<1e-8)&&(rz(i+1)*b>kickrank)
                [crz,~]=localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), kickrank));
            elseif (rz(i+1)*b>ceil(kickrank*ry(i+1)))
                [crz,~]=localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), ceil(kickrank*ry(i+1))));
            end
        end
        
        % Enrich, QR, maxvol
        [y{i},y{i+1},ry(i+1),ind] = qrmaxvol_block(y{i},y{i+1},dir,crz);
        
        % Restrict left indices/matrices
        if (~isempty(ifun))||(~isempty(soughts))
            Jy{i+1} = indexmerge(Jy{i}, (1:n(i))');
            Jy{i+1} = Jy{i+1}(ind, :);
        end
        if (~isempty(X))
            for j=1:nx
                YX{i+1,j} = X{i,j};
                YX{i+1,j} = reshape(YX{i+1,j}, rx(i,j), n(i)*rx(i+1,j));
                YX{i+1,j} = YX{i,j}*YX{i+1,j};
                YX{i+1,j} = reshape(YX{i+1,j}, rx(1,j), ry(i)*n(i), rx(i+1,j));
                YX{i+1,j} = YX{i+1,j}(:, ind, :);
                YX{i+1,j} = reshape(YX{i+1,j}, rx(1,j)*ry(i+1), rx(i+1,j));
            end
        end
        
        % Update the residual itself
        if (kickrank>0)
            cry = reshape(cry, ry(i), n(i)*rz(i+1)*b);
            cry = ZY{i}*cry;
            cry = reshape(cry, rz(i)*n(i)*rz(i+1), b);
            % Evaluate the Z-core
            [crz,ievalcnt,fevalcnt] = evaluate_fun(i,Jz,Jz,n,ifun,ffun,X,rx,ZX,ZX,vec,ievalcnt,fevalcnt);
            crz = crz - cry;
            % Eliminate b
            crz = reshape(crz, rz(i)*n(i), rz(i+1)*b);
            nrmz = norm(crz, 'fro');
            if (nrmz==0)
                crz = randn(rz(i)*n(i), rz(i+1)*b);
            else
                crz = crz/nrmz; % prevent underflows
            end
            if (abs(kickrank-round(kickrank))<1e-8)&&(rz(i+1)*b>kickrank)
                [crz,~] = localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), kickrank));
            elseif (rz(i+1)*b>ceil(kickrank*ry(i+1)))
                [crz,~] = localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), ceil(kickrank*ry(i+1))));
            else
                % Expand crz up to kickrank*ry if necessary
                if (abs(kickrank-round(kickrank))>=1e-8)&&(rz(i+1)*b<ceil(kickrank*ry(i+1)))
                    crz = [crz, randn(rz(i)*n(i), ceil(kickrank*ry(i+1))-rz(i+1)*b)]; %#ok
                end
                [crz,~] = qr(crz, 0);
            end
            rz(i+1) = size(crz, 2);
            ind = maxvol2(crz);
            % Restrict left indices/matrices
            if (~isempty(ifun))||(~isempty(soughts))
                Jz{i+1} = indexmerge(Jz{i}, (1:n(i))');
                Jz{i+1} = Jz{i+1}(ind, :);
            end
            if (~isempty(X))
                for j=1:nx
                    ZX{i+1,j} = X{i,j};
                    ZX{i+1,j} = reshape(ZX{i+1,j}, rx(i,j), n(i)*rx(i+1,j));
                    ZX{i+1,j} = ZX{i,j}*ZX{i+1,j};
                    ZX{i+1,j} = reshape(ZX{i+1,j}, rx(1,j), rz(i)*n(i), rx(i+1,j));
                    ZX{i+1,j} = ZX{i+1,j}(:, ind, :);
                    ZX{i+1,j} = reshape(ZX{i+1,j}, rx(1,j)*rz(i+1), rx(i+1,j));
                end
            end
            ZY{i+1} = reshape(y{i}, ry(i), n(i)*ry(i+1));
            ZY{i+1} = ZY{i}*ZY{i+1};
            ZY{i+1} = reshape(ZY{i+1}, rz(i)*n(i), ry(i+1));
            ZY{i+1} = ZY{i+1}(ind, :);
        end
        
    elseif (dir<0)&&(i>1)
        % Truncation
        ry_old = ry(i);
        [y{i-1},y{i},ry(i),cry] = truncate_block(y{i-1}, cry, tol_local, dir);
        % Enrichment
        crz = [];
        if (kickrank>0)
            cry = reshape(cry, ry_old, n(i)*ry(i+1)*b);
            cry = ZY{i}*cry;
            cry = reshape(cry, rz(i)*n(i)*ry(i+1), b);
            % Evaluate the Z-core
            [crz,ievalcnt,fevalcnt] = evaluate_fun(i,Jz,Jy,n,ifun,ffun,X,rx,ZX,YX,vec,ievalcnt,fevalcnt);
            crz = crz - cry;
            crz = crz.';
            crz = reshape(crz, b*rz(i), n(i)*ry(i+1));
            crz = crz.';
            nrmz = norm(crz, 'fro');
            if (nrmz==0)
                crz = randn(n(i)*ry(i+1), b*rz(i));
            else
                crz = crz/nrmz; % prevent underflows
            end
            % Truncate cry, as we don't need b
            if (abs(kickrank-round(kickrank))<1e-8)&&(b*rz(i)>kickrank)
                [crz,~]=localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), kickrank));
            elseif (b*rz(i)>ceil(kickrank*ry(i)))
                [crz,~] = localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), ceil(kickrank*ry(i))));
            end
            crz = crz.';
        end
        
        % Enrich, QR, maxvol
        [y{i-1},y{i},ry(i),ind] = qrmaxvol_block(y{i-1},y{i},dir,crz);
        
        % Restrict left indices/matrices
        if (~isempty(ifun))||(~isempty(soughts))
            Jy{i} = indexmerge((1:n(i))', Jy{i+1});
            Jy{i} = Jy{i}(ind, :);
        end
        if (~isempty(X))
            for j=1:nx
                YX{i,j} = X{i,j};
                YX{i,j} = reshape(YX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                YX{i,j} = YX{i,j}*YX{i+1,j};
                YX{i,j} = reshape(YX{i,j}, rx(i,j), n(i)*ry(i+1), rx(d+1,j));
                YX{i,j} = YX{i,j}(:, ind, :);
                YX{i,j} = reshape(YX{i,j}, rx(i,j), ry(i)*rx(d+1,j));
            end
        end
        
        % Update the residual itself
        if (kickrank>0)
            cry = cry.';
            cry = reshape(cry, b*rz(i)*n(i), ry(i+1));
            cry = cry*ZY{i+1};
            cry = reshape(cry, b, rz(i)*n(i)*rz(i+1));
            cry = cry.';
            % Evaluate the Z-core
            [crz,ievalcnt,fevalcnt] = evaluate_fun(i,Jz,Jz,n,ifun,ffun,X,rx,ZX,ZX,vec,ievalcnt,fevalcnt);
            crz = crz - cry;
            % Eliminate b
            crz = crz.';
            crz = reshape(crz, b*rz(i), n(i)*rz(i+1));
            crz = crz.';
            nrmz = norm(crz, 'fro');
            if (nrmz==0)
                crz = randn(n(i)*rz(i+1), b*rz(i));
            else
                crz = crz/nrmz; % prevent underflows
            end
            % Truncate cry, as we don't need b
            if (abs(kickrank-round(kickrank))<1e-8)&&(b*rz(i)>kickrank)
                [crz,~]=localcross(crz, tol_local);
                crz = crz(:, 1:min(size(crz,2), kickrank));
            elseif (b*rz(i)>ceil(kickrank*ry(i)))
                [crz,~] = localcross(crz, tol_local);
                crz = crz(:,1:min(size(crz,2), ceil(kickrank*ry(i))));
            else
                % Expand crz up to kickrank*ry if necessary
                if (abs(kickrank-round(kickrank))>=1e-8)&&(b*rz(i)<ceil(kickrank*ry(i)))
                    crz = [crz, randn(n(i)*rz(i+1), ceil(kickrank*ry(i))-b*rz(i))]; %#ok
                end
                [crz,~]=qr(crz, 0);
            end
            rz(i) = size(crz, 2);
            ind = maxvol2(crz);
            % Restrict left indices/matrices
            if (~isempty(ifun))||(~isempty(soughts))
                Jz{i} = indexmerge((1:n(i))', Jz{i+1});
                Jz{i} = Jz{i}(ind, :);
            end
            if (~isempty(X))
                for j=1:nx
                    ZX{i,j} = X{i,j};
                    ZX{i,j} = reshape(ZX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                    ZX{i,j} = ZX{i,j}*ZX{i+1,j};
                    ZX{i,j} = reshape(ZX{i,j}, rx(i,j), n(i)*rz(i+1), rx(d+1,j));
                    ZX{i,j} = ZX{i,j}(:, ind, :);
                    ZX{i,j} = reshape(ZX{i,j}, rx(i,j), rz(i)*rx(d+1,j));
                end
            end
            ZY{i} = reshape(y{i}, ry(i)*n(i), ry(i+1));
            ZY{i} = ZY{i}*ZY{i+1};
            ZY{i} = reshape(ZY{i}, ry(i), n(i)*rz(i+1));
            ZY{i} = ZY{i}(:, ind);
        end
    else
        y{i} = cry;
    end
    
    if (verb>1)
        fprintf('\t-amen_cross_s- swp=%d, i=%d, dx=%3.3e, ranks=[%d,%d], n=%d\n', swp, i, dx, ry(i), ry(i+1), n(i));
    end
    
    
    i = i+dir;
    % Change direction, check for exit
    if ((dir>0)&&(i==d+1))||((dir<0)&&(i==0))
        if (verb>0)
            fprintf('=amen_cross_s= swp=%d, max_dx=%3.3e, max_rank=%d, max_n=%d, cum#ievals=%d, cum#fevals=%d\n', swp, max_dx, max(ry), max(n), ievalcnt, fevalcnt);
        end
        if (max_dx<tol_exit)
            last_swp = last_swp+1;
        end
        if ((last_swp>stop_sweep)||(swp>=nswp))&&((dir~=exitdir)||(exitdir==0))
            break;
        end
        dir = -dir;
        swp = swp+1;
        max_dx = 0;
        i = i+dir;
    end
end

if (dir<0)
    y{1} = reshape(y{1}, n(1), ry(2), b);
    y{1} = permute(y{1}, [3,1,2]);
else
    y{d} = reshape(y{d}, ry(d), n(d), b);
end
y = cell2core(tt_tensor, y);

evalcnt = [ievalcnt fevalcnt];
end


% Truncate a block via the full cross
function [yl,yr,r1,y_trunc] = truncate_block(yl, yr, tol, dir)
if (dir>0)
    % Reshape it
    [r0,nl,r1,b]=size(yl);
    yl = reshape(yl, r0*nl, []); % remaining dimensions may contain b
    if (tol>0)
        % Full-pivot cross should be more accurate
        [yl,rv]=localcross(yl, tol);
    else
        [yl,rv] = qr(yl, 0);
    end
    y_trunc = reshape(yl*rv, r0, nl, r1, b);
    rv = reshape(rv, [], b);
    rv = rv.';
    rv = reshape(rv, [], r1);
    r1 = size(yl,2);
    yl = reshape(yl, r0, nl, r1);
    if (~isempty(yr))
        [~,nr,r2] = size(yr);
        yr = reshape(yr, [], nr*r2);
        yr = rv*yr;
        yr = reshape(yr, b, r1*nr*r2);
        yr = yr.';
        yr = reshape(yr, r1, nr, r2, b);
    end
else
    % Reshape it
    [r1,nr,r2,b]=size(yr);
    yr = reshape(yr, r1, []);
    yr = yr.';
    yr = reshape(yr, nr*r2, b*r1);
    if (tol>0)
        % Full-pivot cross should be more accurate
        [yr,rv]=localcross(yr, tol);
    else
        [yr,rv] = qr(yr,0);
    end
    y_trunc = reshape(yr*rv, nr*r2*b, r1);
    y_trunc = y_trunc.';
    y_trunc = reshape(y_trunc, r1, nr, r2, b);
    rv = reshape(rv, [], r1);
    r1 = size(yr,2);
    yr = yr.';
    yr = reshape(yr, r1, nr, r2);
    if (~isempty(yl))
        [r0,nl,~] = size(yl);
        yl = reshape(yl, r0*nl, []);
        yl = yl*rv.';
        yl = reshape(yl, r0, nl, r1, b);
    end
end
end

% Move non-orth center between the blocks, computing enrich, QR and maxvol
function [yl,yr,r1,ind] = qrmaxvol_block(yl,yr,dir,z)
[r0,nl,~,bl]=size(yl);
[~,nr,r2,br]=size(yr);
if (dir>0)
    % Reshape all
    yl = reshape(yl, r0*nl, []);
    r1 = size(yl,2);
    if (nargin>3)&&(~isempty(z))
        % Enrich
        z = reshape(z, r0*nl, []);
        yl = [yl, z];
    end
    % QR
    [yl,rv]=qr(yl,0);
    rv = rv(:,1:r1);
    % Maxvol and divide
    ind = maxvol2(yl);
    YY = yl(ind, :);
    yl = yl/YY;
    % Update r
    r1 = size(yl,2);
    yl = reshape(yl, r0, nl, r1, bl);
    % Cast non-orths
    if (~isempty(yr))
        rv = YY*rv;
        yr = reshape(yr, [], nr*r2*br);
        yr = rv*yr;
        yr = reshape(yr, r1, nr, r2, br);
    end
else
    % Reshape all
    yr = reshape(yr, [], nr*r2);
    r1 = size(yr,1);
    if (nargin>3)&&(~isempty(z))
        % Enrich
        z = reshape(z, [], nr*r2);
        yr = [yr; z];
    end
    % QR
    yr = yr.';
    [yr,rv]=qr(yr,0);
    rv = rv(:,1:r1);
    % Maxvol and divide
    ind = maxvol2(yr);
    yr = yr.';
    YY = yr(:, ind);
    yr = YY\yr;
    % Update r
    r1 = size(yr,1);
    yr = reshape(yr, r1, nr, r2, br);
    % Cast non-orths
    if (~isempty(yl))
        rv = rv.'*YY;
        yl = reshape(yl, [], bl);
        yl = yl.';
        yl = reshape(yl, bl*r0*nl, []);
        yl = yl*rv;
        yl = reshape(yl, bl, r0*nl*r1);
        yl = yl.';
        yl = reshape(yl, r0, nl, r1, bl);
    end
end
end


% Evaluate the user function at cross indices
function [cry,ievalcnt,fevalcnt] = evaluate_fun(i,Jl,Jr,n,ifun,ffun,crX,rx,YXl,YXr,vec,ievalcnt,fevalcnt)
if (~isempty(ifun))
    J = indexmerge(Jl{i}, (1:n(i))', Jr{i+1});
    if (vec)
        cry = ifun(J);
    else
        % We need to vectorize the fun
        cry = ifun(J(1,:));
        b = numel(cry);
        cry = reshape(cry, 1, b);
        cry = [cry; zeros(size(J,1)-1, b)];
        for j=2:size(J,1)
            cry(j,:) = ifun(J(j,:));
        end
    end
    ievalcnt = ievalcnt + size(J,1);
end
if (~isempty(ffun))
    ry1 = size(YXl{i,1},1)/rx(1,1);
    ry2 = size(YXr{i+1,1},2)/rx(end,1);
    nx = size(crX,2);
    % Compute the X at Y indices
    fx = zeros(ry1*n(i)*ry2, sum(rx(1,:).*rx(end,:)));
    pos = 1;
    for j=1:nx
        cr = crX{i,j};
        cr = reshape(cr, rx(i,j), n(i)*rx(i+1,j));
        cr = YXl{i,j}*cr;
        cr = reshape(cr, rx(1,j)*ry1*n(i), rx(i+1,j));
        cr = cr*YXr{i+1,j};
        cr = reshape(cr, rx(1,j), ry1*n(i)*ry2*rx(end,j));
        cr = cr.';
        cr = reshape(cr, ry1*n(i)*ry2, rx(end,j)*rx(1,j));
        fx(:,pos:pos+rx(end,j)*rx(1,j)-1) = cr;
        pos = pos + rx(end,j)*rx(1,j);
    end
    % Call the function
    fevalcnt = fevalcnt+ry1*n(i)*ry2;
    if (vec)
        fy = ffun(fx); % sizes: rnr x nx -> rnr x d2
    else
        fy = ffun(fx(1,:));
        b = numel(fy);
        fy = reshape(fy, 1, b);
        fy = [fy; zeros(ry1*n(i)*ry2-1, b)];
        for j=2:(ry1*n(i)*ry2)
            fy(j,:) = ffun(fx(j,:));
        end
    end
    if (isempty(ifun))
        cry = fy;
    else
        cry = cry+fy;
    end
end
end


% Merges two or three indices in the little-endian manner
function [J]=indexmerge(varargin)
sz1 = max(size(varargin{1},1),1);
sz2 = max(size(varargin{2},1),1);
sz3 = 1;
if (nargin>2) % Currently allows only 3
    sz3 = max(size(varargin{3}, 1), 1);
end
% J1 goes to the fastest index, just copy it
J1 = repmat(varargin{1}, sz2*sz3, 1);
% J2 goes to the middle
J2 = reshape(varargin{2}, 1, []);
J2 = repmat(J2, sz1, 1); % now sz1 ones will be the fastest
J2 = reshape(J2, sz1*sz2, []);
J2 = repmat(J2, sz3, 1);
J = [J1,J2];
if (nargin>2)
    % J3 goes to the slowest
    J3 = reshape(varargin{3}, 1, []);
    J3 = repmat(J3, sz1*sz2, 1); % now sz1 ones will be the fastest
    J3 = reshape(J3, sz1*sz2*sz3, []);
    J = [J,J3];
end
end

