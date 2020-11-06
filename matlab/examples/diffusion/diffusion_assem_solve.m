function [U,A,F] = diffusion_assem_solve(coeff,bound,W1g,W1m,spind)
% Deterministic solver for the fruitfly
% If Mass is not given, returns a cell of U,A,F at the coefficients coeff
% (used in als_cross).

% Parse coeff
[num_coeffs,nx,I]=size(coeff);
if (num_coeffs>1)
    % Extract RHS separately
    rhs = coeff(2,:,:);
    rhs = reshape(rhs, nx, I);
end
coeff = coeff(1,:,:); % the actual coeff (permeability) is here
% We need to generate separate matrices
n = size(W1g,2);
coeff = reshape(coeff, nx, I);
if (min(coeff(:))<=0)
    warning('non-positive coeff, %g elements', sum(coeff(:)<0)/numel(coeff));
    coeff(coeff<=0) = 1e-8;
end
if (any(isinf(coeff(:))))
    warning('infinite coefficient');
    coeff(isinf(coeff))=1e9;
end

if (nargout==3)
    % Return U,A,F
    A = cell(1,I);
    F = cell(1,I);
    U = cell(1,I);
elseif (nargout==1)&&(nargin<=5)
    % Return U only
    U = cell(1,I);
else
    error('Wrong input/output combination on assem_solve_deterministic');
end
for i=1:I
    % Assemble matrix from W2g and coeff
    B = reshape(coeff(1:n^2,i), n, n);
    B = sparse(B);
    B = W1g*B*W1m' + W1m*B*W1g';
    B = sparse(double(spind(:,1)),double(spind(:,2)),B(spind(:,3)),n^2,n^2); % permute
    if (num_coeffs>1)
        g = reshape(rhs(:,i), n, n);
        localMass = reshape(W1m*sparse(ones(n,1)), n, n);
        g = localMass*g*localMass';
        g = reshape(g, n^2, 1);
    else
        % Eliminate the BC to the RHS
        bound_l = bound(1:numel(bound)/2);
        g = B(:,bound_l)*ones(numel(bound_l),1);
        g = -g;
    end
    g(bound) = [];
    B(bound,:) = [];                                                   %#ok
    B(:,bound) = [];                                                   %#ok
    
    % Prolong matrices
    Pb = speye(n^2);
    Pb(:,bound) = [];                                                  %#ok
    ud = zeros(n^2,1);
    if (num_coeffs==1)
        ud(bound_l) = 1;
    end
    
    % Solve
    u = B\g;
    u = Pb*u;
    u = u+ud;
    
    if (nargout==3)
        % Find the full solution and Return U,A,F
        U{i} = u;
        A{i} = B;
        F{i} = g;
        if (mod(i,100)==0); fprintf('diffusion_assem_solve: assemble/solve %d out of %d problems (full)\n', i, I); end
    elseif (nargout==1)&&(nargin<=5)
        U{i} = u;
        if (mod(i,100)==0); fprintf('diffusion_assem_solve: assemble/solve %d out of %d problems (U)\n', i, I); end
    end
end
end
