function [model, J] = forward_solve(x, model)
% Evaluation of the predator-prey forward model and its Jacobian.
% Stein interface.

ind = model.ind;
obs_times = model.obs_times;
xtrue = model.xtrue;

x = reshape(x, 1, model.n);
X = xtrue;
X(:,ind) = x.*xtrue(ind);
opts = odeset('RelTol', 1e-6);
i = 1;
try
    [~,states]  = ode45(@(t,y)PP_RHS_grad(t,y,X(i,:)), obs_times, [X(i,1:2) reshape(eye(2),1,4) zeros(1,12)], opts);
    gstates = reshape(states(:, 3:end), numel(obs_times), 2, 8);
    model.states = states(:,1:2);
    gstates = reshape(gstates, [], 8);
catch ME
    fprintf(sprintf('%s at parameter %g\n', ME.message, X(i,:)));
end

J = gstates(:, ind).*xtrue(ind);

end