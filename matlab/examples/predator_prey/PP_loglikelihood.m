function [lF] = PP_loglikelihood(x, data, obs_times, sigma_n, xtrue,ind)
% Vectorised log-likelihood of the predator-prey model

I = size(x,1);
X = repmat(xtrue,I,1);
X(:,ind) = x.*xtrue(ind);
opts = odeset('RelTol', 1e-6);
[~,states] = ode45(@(t,y)PP_RHS(t,y,X), obs_times, reshape(X(:,1:2)', 2*I, 1), opts);
states = reshape(states, numel(obs_times)*2, I);
lF = -sum((states - data(:)).^2)*0.5/sigma_n;
lF = lF.';
end
