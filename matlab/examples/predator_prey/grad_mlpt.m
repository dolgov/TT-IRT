function [G] = grad_mlpt(x, model, J, prior, obs)
% Gradient of minus-log-posterior of the predator-prey model in the Stein
% interface

G = (model.states(:) - model.data(:))'*J/model.sigma_n;

end