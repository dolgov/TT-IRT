% prior distribution for shock absorber
function [F] = shock_log_prior(theta, beta_mean, beta_var)
d = size(theta,2)-2;
I = size(theta,1);
% Normal-Gamma
alpha = 6.8757;
beta = 2.2932;
theta2 = theta(:,end);
F = (alpha-0.5).*log(theta2) -beta*theta2 + sum(-((theta(:,1:d+1) - repmat(beta_mean,I,1)).^2).*0.5.*repmat(theta2,1,d+1)./repmat(beta_var,I,1), 2);
end

