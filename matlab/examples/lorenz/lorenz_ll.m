% A function for computing log-likelihood of the Lorenz model
function [L] = lorenz_ll(x0, data, sigma_n)
d = size(x0,2);
x0 = x0(:);
[t,Y] = ode45(@(t,x)lorenz_rhs(t,x,d), [0, 0.1], x0);

Y = Y(end,:); % size nt x (I*d)
Y = reshape(Y, [], d);
Y = Y(:, 2:2:end);

% data is a row vector, auto-broadcast should work
L = -sum((data-Y).^2, 2)*0.5/(sigma_n^2);
end

