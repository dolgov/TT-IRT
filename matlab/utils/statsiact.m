% Convenient IACT function returning the "stats" result (IACT(iid chain)=1)
% in contrast to the "physics" result (IACT(iid chain)=0.5) given by UWerr.
function [tau] = statsiact(x)
[~,~,~,tau,~,~] = UWerr(x,1.5,size(x,1),0);  
tau = tau*2;
end
