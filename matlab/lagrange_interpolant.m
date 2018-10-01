function [p]=lagrange_interpolant(x, y)
%Lagrange interpolation matrix of size numel(y) x numel(x), interpolates
%from x onto y. Backward interpolation is provided by PINV(P). This is a
%vectorised version

n = numel(x);
M = numel(y);
p = ones(M, n);

% All possible differences of xc and xb
dXY = repmat(y(:), 1, n) - repmat(x(:).', M, 1);
% The same for xb,xb
dXX = repmat(x(:), 1, n) - repmat(x(:).', n, 1);
dXX = dXX+eye(n); % Eliminate zeros on the diag
dXX = 1./dXX;
% Agglomerate Lagrange quotients in p
for i=1:n
    % These lengths are necessary to prevent overflow
    pi = dXY.*repmat(dXX(i,:),M,1);
    pi = pi(:, [1:i-1, i+1:n]);
    signs = prod(sign(pi), 2);
    pi = abs(pi);
    pi = log(pi);
    pi = sum(pi,2);
    pi = exp(pi);
    p(:,i) = signs.*pi;
end

end