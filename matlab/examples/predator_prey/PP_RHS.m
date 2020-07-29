function rhs = PP_RHS(t, y, x)
% Vectorised evaluation of the RHS of the predator_prey model
%
% Tiangang Cui, 11 Nov 2014

y       = reshape(y, 2, []);
rhs     = zeros(size(y));

P       = y(1,:).';
Q       = y(2,:).';

r   = x(:,3);
K   = x(:,4);
s   = x(:,5);
a   = x(:,6);
u   = x(:,7);
v   = x(:,8);

tmp     = P.*Q./(a + P);
rhs(1,:)  = r.*P.*(1-P./K) - s.*tmp;
rhs(2,:)  = u.*tmp - v.*Q;

rhs = rhs(:);
end
