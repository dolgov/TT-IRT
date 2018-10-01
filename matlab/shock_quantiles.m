% computation of quantiles, with possible importance weights
function [Q]=shock_quantiles(Z,x,ISweight)
theta1 = exp(Z(:,1) + Z(:,2:end-1)*x);
theta2 = Z(:,end);
if (nargin>2)
    Zex = sum(ISweight); % normalization constant
end
% Quantile function
q = 0.05; % right, set to 0.95 for left
fquant = @(theta1,theta2)theta1.*((-log(q)).^(1./theta2));
q_post = fquant(theta1, theta2);
q_post(q_post>1e7) = 0;
if (nargin>2)
    q_post = sum(q_post.*ISweight)/Zex;
else
    q_post = mean(q_post);
end

% Compute quantile from f_post via Newton method for F_post = q
q_post_new = q_post;
for iter=1:20
    % residual
    R = exp(-(q_post_new./theta1).^theta2);
    if (nargin>2)
        R = (sum(R.*ISweight)/Zex)/q-1;
    else
        R = mean(R)/q-1;
    end
    % -Jacobian
    J = exp(-(q_post_new./theta1).^theta2).*(theta2./theta1).*(q_post_new./theta1).^(theta2-1);
    if (nargin>2)
        J = (sum(J.*ISweight)/Zex)/q;
    else
        J = mean(J)/q;
    end
    % update
    q_post_new = q_post_new + J\R;
    fprintf('newton for f_post quantile iter=%d, resid=%3.3e\n', iter, abs(R));
    if (abs(R)<1e-8), break; end
end
Q = [q_post q_post_new];
end

