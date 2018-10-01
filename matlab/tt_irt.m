function [xq, Fapp]=tt_irt(xsf, f, q)
% Inverse CDF (Rosenblatt) transform through the density F
% This version computes inverses in each variable via quadratic splines
% Inputs (sizes):
%   xsf: (d x 1) cell array of grid points (inc. boundaries) for all dimensions
%   f: TT format of the PDF (tt_tensor), computed on a grid 
%      _without left borders_, i.e. f must be evaluated at xsf{i}(2:end) 
%      in the i-th direction.
%   q: seed points from [0,1]^d, defining the samples (an M x d matrix)
% Outputs (sizes):
%   xq: samples mapped from q by the inverse CDF (M x d matrix)
%   Fapp: approximate PDF (inv. Jacobian of xq) sampled on q


if (isa(f, 'tt_tensor'))
    f = core2cell(f);
end
% Parse TT: extract dimensions and ranks
d = numel(f);
n = cellfun(@(x)size(x,2), f);
rf = [1; cellfun(@(x)size(x,3), f)];

% Check grid points
if (isa(xsf, 'cell'))
    xsf = cell2mat(xsf);
end
if (size(xsf,1)~=sum(n+1))
    error('number of grid points in xf should be n+1');
end
pos = cumsum([1; n+1]); % positions where each x grid vector starts

% Prepare integrated right TT blocks (C array)
C = cell(d+1,1);
C{d+1} = 1;
h = cell(d,1); % array of grid intervals
for k=d:-1:1
    % Grid intervals
    h{k} = xsf(pos(k)+1:pos(k)+n(k)) - xsf(pos(k):pos(k)+n(k)-1);
    % Linear-spline integration of (k..d)-th TT blocks
    C{k} = reshape(f{k}, rf(k)*n(k), rf(k+1));
    C{k} = C{k}*C{k+1};
    C{k} = reshape(C{k}, rf(k), n(k));
    % extrapolate it linearly to the leftmost point, since f^(k) does not
    % contain it
    Aq = (xsf(pos(k)+2)-xsf(pos(k)))/h{k}(2);
    Bq = (xsf(pos(k))-xsf(pos(k)+1))/h{k}(2);
    fk = Aq*C{k}(:,1) + Bq*C{k}(:,2);
    C{k} = [fk, C{k}];
    
    % integrate via trapezoidal rule. It's faster using a sparse matrix
    fk = reshape(h{k}, 1, n(k));
    fk = repmat(fk, rf(k), 1);
    S = spdiags(0.5*ones(n(k)+1,2), -1:0, n(k)+1, n(k));
    C{k} = (C{k}*S).*fk;
    C{k} = [zeros(rf(k),1), C{k}];
    C{k} = sum(C{k}, 2);
end
% Now each C{k} = \int f^{(>k)}(x_k...x_d) dx_{>k}

% Storage for x and approx. density
M = size(q,1);
xq = zeros(M,d);
if (nargout>1)
    Fapp = ones(M,1);
end

% Iterate forward, computing conditional PDFs
% Prevent OOM by blocking
Mb = 2^16;
num_blocks = ceil(M/Mb);
for i_block=1:num_blocks
    % Determine the position of the current block
    start_pos = (i_block-1)*Mb + 1;
    if (i_block*Mb>M)
        Mb = M - start_pos + 1;
    end
    fkm1 = ones(Mb,1); % This will store conditioned left TT blocks f^{(<k)}(x_{<k}^*)
    
    for k=1:d
        %%% Prepare 1D PDF and CDF
        % Prepare semi-marginal f^{(k)}*C{k}, to be integrated into a CDF
        Ck = reshape(f{k}, rf(k)*n(k), rf(k+1));
        Ck = Ck*C{k+1};
        Ck = reshape(Ck, rf(k), n(k));
        % extrapolate it to leftmost point
        Aq = (xsf(pos(k)+2)-xsf(pos(k)))/h{k}(2);
        Bq = (xsf(pos(k))-xsf(pos(k)+1))/h{k}(2);
        fk = Aq*Ck(:,1) + Bq*Ck(:,2);
        fk = [fk, Ck];
        % Condition on left variables -- multiply with sampled left TT blocks
        fk = fkm1*fk;
        % make this nonnegative. ALT: subtract the minimal neg value?
        fk = abs(fk);
%         fk = fk - repmat(min(fk,[],2), 1, n(k)+1); % doesn't work well...
        % Integrate up to whole points using sparse mv
        Ck = reshape(h{k}, 1, n(k));
        Ck = repmat(Ck, Mb, 1);
        S = spdiags(0.5*ones(n(k)+1,2), -1:0, n(k)+1, n(k));
        Ck = (fk*S).*Ck;
        Ck = [zeros(Mb,1), Ck];
        Ck = cumsum(Ck, 2);
        
        % normalize
        Cmax = Ck(:,n(k)+1);
        ind_zero = find(Cmax==0);
        Cmax = repmat(Cmax, 1, n(k)+1);
        Ck = Ck./Cmax;
        fk = fk./Cmax;
        % Eliminate possible exact zeros
        Ck(ind_zero, :) = repmat(0:(1/n(k)):1, numel(ind_zero), 1);
        fk(ind_zero, :) = (1/n(k))*ones(numel(ind_zero), n(k)+1);
        
        %%% Invert the marginal CDF Ck
        % Binary search for the closest to q points
        qk = q(start_pos:(start_pos+Mb-1),k);
        i0 = ones(Mb,1);
        i2 = (n(k)+1)*ones(Mb,1);
        while (any((i2-i0)>1))
            i1 = floor((i0+i2)/2);
            C1 = tracemult(Ck, i1); % MEXed extraction C1(i) = Ck(i,i1(i));
            ind_left = qk>C1; % at these indices i0 becomes i1
            i0(ind_left) = i1(ind_left);
            i2(~ind_left) = i1(~ind_left);
        end
        
        % Sample C and f onto i0:i0+1
        C1 = tracemult(Ck,i0);      % C1(i) = Ck(i,i0(i));
        C2 = tracemult(Ck, i0+1);   % C2(i) = Ck(i,i0(i)+1);
        f1 = tracemult(fk,i0);      % f1(i) = fk(i,i0(i));
        f2 = tracemult(fk, i0+1);   % f1(i) = fk(i,i0(i)+1);
        % We can miss a point, check that
        ind_missed = find((C1>(qk+1e-8)) | (C2<(qk-1e-8)));
        % Solve the _quadratic_ equation, since C is an integrand of linear
        % interpolant f
        x1 = xsf(pos(k)+i0-1);  % grid points at i0, i0+1
        x2 = xsf(pos(k)+i0);
        h3 = x2 - x1; % grid intervals
        Aq = 0.5*(f2-f1)./h3;
        Bq = (f1.*x2 - f2.*x1)./h3;
        Dq = (2*Aq.*x1+Bq).^2 + 4*Aq.*(qk-C1);
        xk = 0.5*(-Bq + sqrt(abs(Dq)))./Aq;
        ind_left = find(Aq==0); % At these indices we need to solve a linear equation
        xk(ind_left) = x1(ind_left) + (qk(ind_left)-C1(ind_left))./Bq(ind_left);
        % wrong bin search can lead to a wrong solution. Replace those xk
        % by midpoints
        xk(ind_missed) = (x1(ind_missed) + x2(ind_missed))*0.5;
        
        xq(start_pos:(start_pos+Mb-1),k) = xk;
        
        % Linear interpolation coefficients into xk
        Aq = (x2-xk)./h3;
        Bq = (xk-x1)./h3;
        if (nargout>1)
            % Interpolate fk into xk:
            % this is our actual conditional probability at the k-th step
            % The whole density is a product of those
            fk = f1.*Aq + f2.*Bq;
            Fapp(start_pos:(start_pos+Mb-1)) = Fapp(start_pos:(start_pos+Mb-1)).*fk;
        end
        
        if (k<d)
            % Sample the k-th block onto xk
            % Build the linear interpolant
            fk = permute(f{k}, [1,3,2]);
            fk = reshape(fk, rf(k)*rf(k+1), n(k));
            % extrapolate fk to the leftmost point
            C1 = (xsf(pos(k)+2)-xsf(pos(k)))/h{k}(2);
            C2 = (xsf(pos(k))-xsf(pos(k)+1))/h{k}(2);
            Ck = C1*fk(:,1) + C2*fk(:,2);
            fk = [Ck, fk];
            fk = reshape(fk, rf(k), rf(k+1), n(k)+1);
            Aq = repmat(Aq, 1, rf(k+1));
            Bq = repmat(Bq, 1, rf(k+1));
            % On-the-fly product
            % fkm1(i,:)*fk(:,i0(i),:)*Aq(i) + fkm1(i,:)*fk(:,i0(i)+1,:)*Bq(i)
            fkm1 = tracemult(fkm1, i0, fk).*Aq + tracemult(fkm1, i0+1, fk).*Bq;
        end
    end
end
end
