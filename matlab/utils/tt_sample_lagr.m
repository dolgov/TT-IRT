function [Q] = tt_sample_lagr(u, x, y)
% Lagrange interpolation of u onto a set of points y.
% Inputs:
%   u: tt_tensor or its core2cell
%   x: (d,1) cell of Cheb2 points in each dimension
%   y: (M,d) double of points
% Output:
%   Q: (M,r) double of values, where r is a max. boundary TT rank of u

if (isa(u, 'tt_tensor'))
    u = core2cell(u);
end

d = numel(u);
P = cell(d,1);
for i=1:d
    P{i} = lagrange_interpolant(x{i}, y(:,i));
end

r = cellfun(@(x)size(x,1), u);
r = [r; size(u{d},3)];

M = size(y,1);
if (r(1)>1)
    Q = zeros(M,r(1));
    for m=1:M
        um = 1;
        for i=d:-1:1
            ui = u{i}; % (:, ind{i}(m,:), :);
            ui = reshape(ui, [], r(i+1));
            um = ui*um; % complexity rn x r 
            um = reshape(um, r(i), []);
            um = um*P{i}(m,:).'; % complexity r x n
        end
        Q(m,:) = um;
    end
else
    Q = zeros(M,r(d+1));
    for m=1:M
        um = 1;
        for i=1:d
            ui = u{i}; % (:, ind{i}(m,:), :);
            ui = reshape(ui, r(i), []);
            um = um*ui;
            um = reshape(um, [], r(i+1));
            um = P{i}(m,:)*um;
        end
        Q(m,:) = um;
    end
end
end
