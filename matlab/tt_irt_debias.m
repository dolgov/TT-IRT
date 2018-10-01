% Debiasing of Inverse Rosenblatt TT using MCMC or Importance Weighting
% Inputs:
%   M: samples (in a form of a M x d matrix) or number of samples
%   Ffun: function of y (in form M x d) returning the exact PDF and
%         optionally any other R functions at y in form of M x (R+1) matrix.
%         The first column should contain the exact PDF values.
%   f: TT format of the PDF (tt_tensor), computed on a grid 
%      _without left borders_, i.e. f must be evaluated at xsf{i}(2:end) 
%      in the i-th direction.
%   xsf: (d x 1) cell array of grid points (inc. boundaries) for all dimensions
%   correction: debiasing algorithm, 'mcmc' or 'iw'
%   [vec]: whether Ffun accepts vectorized input [true]
% Output:
%   y: M x d matrix of samples, distributed according to exact PDF for
%      MCMC, and approximate PDF for IS
%   Fex: correct Ffun sampled on y, in the form M x (R+1)
%   bias: a measure of bias in the proposed samples, that is:
%           the number of rejections for MCMC, and 
%           the relative standard deviation of the ratio Fex/Fapp for IW
%   time_invcdf: CPU time spent in IRT

function [y,Fex,bias,time_invcdf] = tt_irt_debias(M, Ffun, f, xsf, correction, vec)
d = f.d;
if (isscalar(M))
    Z = rand(M,d);
else
    Z = M;
    M = size(Z,1);
end
if (sum(cellfun(@length, xsf))~=sum(f.n+1))
    error('number of grid points in xsf should be n+1');
end

if (nargin<6)
    vec = true;
end

% Propose IRT samples
t_mcmc = tic;
[y,Fapp] = tt_irt_mex(f.n, cell2mat(xsf), f.r, f.core, Z);
time_invcdf = toc(t_mcmc);
fprintf('invcdf time = %g\n', time_invcdf);

t_mcmc = tic;
if (vec)
    Fex = Ffun(y);
else
    % user function is not vectorised
    Fex = Ffun(y(1,:));
    b = numel(Fex);
    Fex = reshape(Fex, 1, b);
    Fex = [Fex; zeros(M-1, b)];
    for i=2:M
        Fex(i,:) = Ffun(y(i,:));
    end
end
fprintf('Lex time=%g\n', toc(t_mcmc));

if (strcmpi(correction, 'mcmc'))
    % Rejection loop here
    [y,Fex,~,bias] = mcmc_prune(y,Fex,Fapp);    
else
    % Importance Weighting: exact QoI = QoI*Fex/Fapp
    [Fex,bias] = iw_prune(Fex,Fapp);
end
end

