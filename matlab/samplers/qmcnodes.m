% Generates 2^l QMC lattice nodes in d dimensions
function [Y]=qmcnodes(d,l)
% N = 2^l
fname = 'lattice-39102-1024-1048576.3600.txt';
z = load(fname);
Y = 0:(2^l-1);
Y = Y/(2^l);
Y = z(1:d,2)*Y;
% Random shift
Delta = rand(d,1);
Delta = repmat(Delta, 1, 2^l);
Y = Y+Delta;
Y = Y-floor(Y);
end
