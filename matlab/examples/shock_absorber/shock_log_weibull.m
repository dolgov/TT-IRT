% Weibull likelihood
function [F]=shock_log_weibull(theta, x, y, c)
d = size(theta,2)-2;
beta = theta(:,1:d+1);
lambda = theta(:,d+2);
F = zeros(size(theta,1), 1);
m = numel(y);
if (nargin<4)
    c = zeros(m,1); % censoring
end
for i=1:m
    logeta = beta(:,1) + beta(:,2:d+1)*x(:,i);
    eta = exp(logeta);
    yeta = y(i)./eta;
    if (c(i))
        f = -(yeta.^lambda); % censored PDF = 1-CDF
    else
        f = log(lambda) - logeta + (lambda-1).*(log(y(i)) - logeta) - (yeta.^lambda);        
        f = f + log(3e4); % to prevent underflow
    end
    F = F + f;
end
end
