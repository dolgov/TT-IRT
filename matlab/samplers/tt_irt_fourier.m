function [xq, lFapp]=tt_irt_fourier(xsf, f, q)
% Inverse CDF (Rosenblatt) transform through the square root of the density F
% This version uses Fourier basis to interpolate and invert the CDFs
% Inputs (sizes):
%   xsf: (d x 1) cell array of grid points for all dimensions
%       !!! must be uniform without left border, for example ((-n+1):n)*h
%   f: TT format of the SQRT(PDF) (tt_tensor), computed on a grid given in xsf, or
%      a structure containing precomputed irt_fourier data
%   q: seed points from [0,1]^d, defining the samples (an M x d matrix)
%      If empty, the first output is the precomputed structure for further
%      sampling
% Outputs (sizes):
%   xq: samples mapped from q by the inverse CDF (M x d matrix)
%       (or an IRT structure if q is missing)
%   lFapp: log(approximate PDF) (inv. Jacobian of xq) sampled on q
%
% See also: tt_irt_lin, tt_irt_sqr


newton_iter_max = 16;       % Those should be OK for most purposes
newton_resid_tol = 1e-7;


if (isa(f, 'tt_tensor'))
    f = core2cell(f);
end

if (isa(f, 'cell'))
    % Parse TT: extract dimensions and ranks, and precompute
    d = numel(f);
    n = cellfun(@(x)size(x,2), f);
    rf = [1; cellfun(@(x)size(x,3), f)];
    
    % Check grid points
    if (isa(xsf, 'cell'))
        xsf = cell2mat(xsf);
    end
    if (size(xsf,1)~=sum(n))
        error('number of grid points in xsf should be sum of mode sizes in f');
    end
    pos = cumsum([1; n]); % positions where each x grid vector starts
    
    if (any(mod(n,2)~=0))
        error('All numbers of grid points should be even');
    end
    n = n/2; % for convenient FFT
    
    % Prepare integrated sqr(right TT blocks)
    P = cell(d+1,1);
    P{d+1} = 1;
    Rprev = 1;
    h = zeros(d,1); % step sizes of (uniform) grids
    S = zeros(d,1); % Width of the period interval (-S,S]
    X0 = zeros(d,1); % midpoint of the original interval (a,b]
    fft_factor = cell(d,1); % these will map discrete and continuous FT 
    ifft_factor = cell(d,1);
    x_fine = cell(d,1); % 2x finer grid needed for squared functions
    for k=d:-1:1
        % Grid intervals
        h(k) = xsf(pos(k)+1) - xsf(pos(k));
        S(k) = xsf(pos(k+1)-1);
        % Check for consistency
        if any(abs(xsf(pos(k)+1:pos(k+1)-1) - xsf(pos(k):pos(k+1)-2) - h(k))>1e-12)
            error('grid is not uniform in dimension %d', k);
        end
        X0(k) = (xsf(pos(k))-h(k) + S(k))*0.5; % Keep the midpoint separately
        S(k) = S(k) - X0(k);
        P{k} = reshape(f{k}, rf(k)*2*n(k), rf(k+1));
        P{k} = P{k}*Rprev; % Rprev is the right TT chunk integrated over x_{>k}
        P{k} = reshape(P{k}, rf(k), 2*n(k), []);
        if (k>1)
            % Integrate over x_k for the next step implicitly using
            % QR decomposition since (diag(w)P)'*(diag(w)P) = R'*R
            % Rectangular quadrature is best for Fourier
            Rprev = reshape(P{k}, rf(k), [])*sqrt(h(k));           
            [~,Rprev] = qr(Rprev.', 0); % size of R is rnew x rf(k)
            Rprev = Rprev.';
        end
        % Prepare pointwise factors for doing continuous FTs via fft
        fft_factor{k} = exp(-1i*(pi/n(k))*(1-n(k))*(0:(2*n(k)-1))) * (exp(-1i*0.5*(pi/n(k))*(1-n(k))^2)*sqrt(h(k)));
        fft_factor{k} = reshape(fft_factor{k}, 1, 1, 2*n(k));
        ifft_factor{k} = exp(1i*0.5*(pi/n(k))*(1-2*n(k))*(0:(4*n(k)-1))) * (exp(1i*0.25*(pi/n(k))*(1-2*n(k))^2)*sqrt(n(k))/S(k));
        x_fine{k} = ((-2*n(k)+1):(2*n(k)))*(h(k)/2);
        if (k<d)
            % Precompute FT of the k-th core for conditioning interpolation
            f{k} = permute(f{k}, [1,3,2]);
            f{k} = fft_factor{k}.*fft(fft_factor{k}.*f{k}, [], 3);
            f{k} = reshape(f{k}, rf(k)*rf(k+1), 2*n(k));
        end
        % Prepare semi-marginal (f^{(k)}*R{k+1})^2, to be integrated into a CDF
        P{k} = permute(P{k}, [1,3,2]);
        P{k} = fft_factor{k}.*fft(fft_factor{k}.*P{k}, [], 3);
        % Get Fourier coefficients of the square, simultaneously integrating
        % \int f^{(>=k)}(x_k...x_d) x (f^{(>=k)}(x_k...x_d))^T dx_{>k}
        % due to right TT cores already collapsed by QR
        fk = zeros(rf(k), rf(k), 4*n(k));
        for i=(-n(k)+1):n(k)
            for j=(-n(k)+1):n(k)
                fk(:,:, i-j+2*n(k)) = fk(:,:, i-j+2*n(k)) + P{k}(:,:,i+n(k))*P{k}(:,:,j+n(k))';
            end
        end
        fk = reshape(fk, rf(k)^2, 4*n(k));
        P{k} = fk;
    end
    
elseif (isa(f, 'struct'))
    % We have all precomputies, ready to sample
    h = f.h;
    S = f.S;
    X0 = f.X0;
    x_fine = f.x_fine;
    fft_factor = f.fft_factor;
    n = cellfun(@numel, fft_factor)/2;
    d = numel(n);
    ifft_factor = f.ifft_factor;
    P = f.P;
    f = f.f;
    rf = cellfun(@(x)size(x,1), P);
    rf = round(sqrt(rf));
else
    error('Unknown type of f. It can be a tt_tensor, a cell of TT cores or a precomputed struct');
end


if (nargin<3)||(isempty(q))
    % Just return the structure
    xq = struct();
    xq.h = h;
    xq.S = S;
    xq.X0 = X0;
    xq.x_fine = x_fine;
    xq.fft_factor = fft_factor;
    xq.ifft_factor = ifft_factor;
    xq.f = f;
    xq.P = P;
    return;
end



M = size(q,1);
% Iterate forward, computing conditional PDFs
% Prevent OOM by blocking
Mb = 2^8;
% We can embarassingly parallelise, but we don't want to create a pool if
% this is not wanted
num_blocks = gcp('nocreate');
if (~isempty(num_blocks))
    num_blocks = num_blocks.NumWorkers;
    Mb = min(Mb, M/num_blocks);
end
num_blocks = ceil(M/Mb);
% Determine the local range
ibeg = round((0:num_blocks-1)'*M/num_blocks)+1;
iend = round((1:num_blocks)'*M/num_blocks);
rerange = arrayfun(@(i1,i2)(i1:i2)', ibeg, iend, 'UniformOutput', false);
% Local seeds
qrange = cellfun(@(r)q(r,:), rerange, 'UniformOutput', false);
% Storage for approx. density
lFapp = cell(numel(rerange), 1);
if (nargout>1)
    lFapp = cellfun(@(r)zeros(numel(r),1), rerange, 'UniformOutput', false);
end
% Storage for x
xq = cell(numel(rerange), 1);
% Each (parallel) block will work with its own qrange{i_block}, etc.

for i_block=1:num_blocks
    Mb = numel(rerange{i_block});
    fkm1 = ones(1, Mb); % This will store conditioned left TT blocks f^{(<k)}(x_{<k}^*)
    
    for k=1:min(d,size(q,2))    % Simple trick to sample marginal if needed
        %%% Prepare 1D PDF and CDF
        fkm1 = reshape(fkm1, rf(k), 1, Mb);
        % Square conditioned \pi^{(<k)}(x_1...x_{k-1})
        fk = tracemult(fkm1, (1:Mb)', permute(conj(fkm1),[2,1,3]));    
        fk = reshape(fk, rf(k)^2, Mb);
        % Condition k-th marginal on x_{<k}
        fk = fk.'*P{k}; % size Mb x 4n(k)
        % Eliminate possible exact zeros
        ind_zero = all(fk==0, 2);
        fk(ind_zero, 2*n(k)) = 1; % Fourier image of a const
        % fk is already the Fourier image of a nonneg function
        % Integrate into CDF in the Fourier space directly
        CkF = fk./(1i*(pi/S(k))*((-2*n(k)+1):(2*n(k)))); % Fourier image for (CDF - linear)
        CkF(:,2*n(k)) = 0;
        Ak = real(fk(:,2*n(k))/(4*S(k)^2)); % Coeff of the linear part
        Bk = Ak*S(k) - real(CkF*exp(-1i*pi*((-2*n(k)+1):(2*n(k)))'))/(4*S(k)^2); % constant shift 
        % Spatial domain
        Ck = Ak.*x_fine{k} + Bk + real(ifft_factor{k}.*ifft(ifft_factor{k}.*CkF, [], 2));

        % normalize
        Cmax = Ck(:,4*n(k));
        Ck = Ck./Cmax;   % Space domain for binary search
        CkF = CkF./Cmax; % still in Fourier for sampling
        Ak = Ak./Cmax;
        Bk = Bk./Cmax;
        fk = fk./Cmax; 
        
        %%% Invert the marginal CDF Ck
        % Binary search for the closest to q points
        qk = qrange{i_block}(:, k);
        i0 = ones(Mb,1);
        i2 = 4*n(k)*ones(Mb,1);
        while (any((i2-i0)>1))
            i1 = floor((i0+i2)/2);
            C1 = tracemult(Ck, i1); % MEXed extraction C1(i) = Ck(i,i1(i));
            ind_left = qk>C1; % at these indices i0 becomes i1
            i0(ind_left) = i1(ind_left);
            i2(~ind_left) = i1(~ind_left);
        end
        
        % Sample C and f onto i0:i0+1
        C1 = tracemult(Ck, i0);     % C1(i) = Ck(i,i0(i));
        x1 = reshape(x_fine{k}(i0), 1, 1, Mb);
        % Reshape for tracemult 
        qk = reshape(qk, 1, 1, Mb);
        C1 = reshape(C1, 1, 1, Mb);
        Ak = reshape(Ak, 1, 1, Mb);
        Bk = reshape(Bk, 1, 1, Mb);        
        CkF = reshape(CkF.', 4*n(k), 1, Mb);
        fk = reshape(fk.', 4*n(k), 1, Mb);
        % Absorb FFT norms up front
        CkF = CkF/(4*S(k)^2);
        fk = fk/(4*S(k)^2);                
        % Frequency grid
        ifreq = 1i*(pi/S(k))*((-2*n(k)+1):(2*n(k)));
        
        % Solve the quadratic interpolant on C for initial guess
        iF = exp(ifreq.*x1);
        f1 = tracemult(iF, 1:Mb, fk);
        iF = iF.*exp(ifreq*h(k)*0.5);
        f2 = tracemult(iF, 1:Mb, fk);
        f1 = real(f1);
        f2 = real(f2);
        Aq = 0.5*(f2-f1)./(0.5*h(k));
        Dq = f1.^2 + 4*Aq.*(qk-C1);
        xk = x1 + (-f1 + sqrt(abs(Dq)))./(2*Aq);
        ind_missed = find(Aq==0); % At these indices we need to solve a linear equation
        xk(ind_missed) = x1(ind_missed) + (qk(ind_missed)-C1(ind_missed))./f1(ind_missed);
        ind_missed = find((f1==0)&(Aq==0));
        xk(ind_missed) = x1(ind_missed);
        % Restrict the points to the current interval exactly
        xk = max(xk, -S(k));
        xk = min(xk,  S(k));

        % Newton search for the exact root
        x_active = reshape(xk, 1, 1, Mb);
        active_pos = 1:Mb;
        Jac_all = ones(Mb,1); % storage for Jacobians
        for iter=1:newton_iter_max
            % Fourier transform into unstructured points must be done explicitly
            iF = exp(ifreq.*x_active);
            Jac = tracemult(iF, 1:numel(x_active), fk);
            Jac = abs(Jac);
            Resid = tracemult(iF, 1:numel(x_active), CkF);
            Resid = real(Resid) + Ak.*x_active + Bk - qk;
            xk(active_pos) = x_active;
            Jac_all(active_pos) = Jac;
            active_ind = abs(Resid)>newton_resid_tol; % Release already converged points
            x_active = x_active(active_ind);
%             fprintf('iter=%d, xq(1)=%g, max(Resid)=%3.3e\n', iter, xk(1), max(abs(Resid)));
            if (isempty(x_active))
                break;
            end            
            Ak = Ak(active_ind);
            Bk = Bk(active_ind);
            qk = qk(active_ind);
            fk = fk(:, :, active_ind);
            CkF = CkF(:, :, active_ind);
            Resid = Resid(active_ind);
            Jac = Jac(active_ind);
            active_pos = active_pos(active_ind);
            % Stabilized Newton since Jac can be <<|Resid| on tails
            x_active = x_active - Resid./(Jac+abs(Resid));
            % Restrict the points to the current interval exactly
            x_active = max(x_active, -S(k));
            x_active = min(x_active,  S(k));
        end
        
        % Record the roots into global storage
        xq{i_block}(:, k) = X0(k) + xk;
        
        if (~isempty(lFapp{i_block}))
            % Jacobian is our actual conditional probability at the k-th step
            % The whole density is a product of those
            lFapp{i_block} = lFapp{i_block} + log(Jac_all);
        end
        
        if (k<d)
            % Sample the k-th block onto xk
            % Full forward FT is precomputed next to P{k} above
            % Inverse transform (original size) into sample points
            xk = reshape(xk, 1, Mb);
            iF = exp(1i*(pi/S(k))*((-n(k)+1):n(k))'.*xk)/(2*S(k));
            fk = f{k}*iF;
            fk = reshape(fk, rf(k), rf(k+1), Mb);
            fkm1 = reshape(fkm1, 1, rf(k), Mb);
            fkm1 = tracemult(fkm1, (1:Mb)', fk); % size 1 x r2 x Mv
            fkm1 = reshape(fkm1, rf(k+1), Mb);
        end
    end
end

xq = cat(1,xq{:});
if (nargout>1)
    lFapp = cat(1,lFapp{:});
end
end
