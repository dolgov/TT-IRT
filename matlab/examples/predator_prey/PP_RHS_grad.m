function rhs = PP_RHS_grad(t, y, x)
% Evaluation of the [state; gradient] RHS of the predator_prey model
%
% Tiangang Cui, 11 Nov 2014

rhs      = zeros(2, 9);
% for state
rhs(1,1) = x(3)*y(1)*(1-y(1)/x(4)) - x(5)*y(1)*y(2)/(x(6)+y(1));
rhs(2,1) = x(7)*y(1)*y(2)/(x(6)+y(1)) - x(8)*y(2);

% for gradient, d/dt(dx/dtheta) = df/dtheta + (df/dx)*(dx/dtheta)
% df/dtheta first
rhs(1,4) = y(1)*(1-y(1)/x(4)); % f depends only on parameters, not P0,Q0
rhs(1,5) = x(3)*y(1)^2/x(4)^2;
rhs(1,6) = -y(1)*y(2)/(x(6)+y(1));
rhs(1,7) = x(5)*y(1)*y(2)/(x(6)+y(1))^2;
rhs(2,7) = -x(7)*y(1)*y(2)/(x(6)+y(1))^2;
rhs(2,8) = y(1)*y(2)/(x(6)+y(1));
rhs(2,9) = -y(2);

% df/dx
dFdx = zeros(2,2);
g1 = y(2)/(x(6)+y(1)) - y(1)*y(2)/(x(6)+y(1))^2;
g2 = y(1)/(x(6)+y(1));
dFdx(1,1) = x(3)*(1-2*y(1)/x(4)) - x(5)*g1;
dFdx(1,2) = -x(5)*g2;
dFdx(2,1) = x(7)*g1;
dFdx(2,2) = x(7)*g2 - x(8);

rhs(:,2:9) = rhs(:,2:9) + dFdx*reshape(y(3:end), 2, 8);

rhs = rhs(:);
end
