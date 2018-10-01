function [Q] = tt_sample_lagr(u, x, y, nloc)
% Lagrange interpolation of u onto a set of points y.
% Inputs:
%   u: tt_tensor or its core2cell
%   x: (d,1) cell of Cheb2 points in each dimension
%   y: (M,d) double of points
%   nloc: (1) or (d,1) sizes of local h-intervals in hp-grids
% Output:
%   Q: (M,r) double of values, where r is a max. boundary TT rank of u

if (isa(u, 'tt_tensor'))
    u = core2cell(u);
end

d = numel(u);
P = cell(d,1);
ind = cell(d,1);
if (nargin>3)&&(~isempty(nloc))
    if (numel(nloc)==1)
        nloc = nloc*ones(d,1);
    end
    for i=1:d
        [P{i}, ind{i}] = cheb2_interpolant(x{i}, y(:,i), nloc(i)); % size M x ni
    end
else
    for i=1:d
        P{i} = lagrange_interpolant(x{i}, y(:,i));
    end
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
