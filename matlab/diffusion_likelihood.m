function [LQF]=diffusion_likelihood(y, phil0, sigma, bound, W1g, W1m, spind, Mass, Q_obs, sigma_n)
% This will be a function in the exact posterior
M = size(y,1);
m0 = size(Mass,1);

% Storages for outputs
Q = zeros(M, m0^2);
L = zeros(M,1);
Fl = zeros(M,1);

n = size(W1g,2);

% Matrices for flux
W1m_sum = reshape(W1m, n, n^2); % ind i, jk
W1m_sum = sum(W1m_sum, 1); % sum over i2 here
W1m_sum = reshape(W1m_sum, n, n); % i,k
W1g_sum = reshape(W1g, n, n^2); % ind i, jk
W1g_sum = sum(W1g_sum, 1); % sum over i2 here
W1g_sum = reshape(W1g_sum, n, n); % i,k

% Solve the forward problems
for i=1:M
    C = phil0*y(i,:).'; % size n x M
    C = exp(C*sqrt(sigma));
    
    % Assemble matrix from W1g and coeff
    B = reshape(C(1:n^2,1), n, n);
    B = sparse(B);
    B = W1g*B*W1m' + W1m*B*W1g';
    B = sparse(double(spind(:,1)),double(spind(:,2)),B(spind(:,3)),n^2,n^2); % permute
    % Eliminate the BC to the RHS
    bound_l = bound(1:numel(bound)/2);
    g = B(:,bound_l)*ones(numel(bound_l),1);
    g = -g;
    g(bound) = [];
    B(bound,:) = [];                                                   %#ok
    B(:,bound) = [];                                                   %#ok
    
    % Prolong matrices
    Pb = speye(n^2);
    Pb(:,bound) = [];                                                  %#ok
    ud = zeros(n^2,1);
    ud(bound_l) = 1;
    
    % Solve
    u = B\g;
    u = Pb*u;
    u = u+ud;
    
    if (mod(i,100)==0); fprintf('diffusion_likelihood: computed %d out of %d values\n', i, M); end

    % Observations and residual
    for j=1:m0
        for k=1:m0
            Q(i, k+(j-1)*m0) = sum(Mass{k,j}*u, 1);
            L(i) = L(i) + (Q(i, k+(j-1)*m0) - Q_obs(k+(j-1)*m0)).^2;
        end
    end

    % Flux
    B = reshape(C(1:n^2,1), n, n);
    B = sparse(B);
    B = W1g*B*W1m_sum' + W1m*B*W1g_sum'; % indices i1,j1,j2
    B = reshape(B, n, n^2);
    % 1D flux
    B = B*u; % index i1
    B = full(B);
    Fl(i) = -B(n);
end

% Likelihood
L = exp(-L/(2*sigma_n));

IFL = double(Fl>1.5);

LQF = [L, Fl, IFL]; % size M x 3, first component is L, last is Flux indic
end
