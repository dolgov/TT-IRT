% Debiasing of Inverse Rosenblatt TT using MCMC or Importance Weighting
% This function uses tt_irt_lin based on linear spline interpolation of PDF
% Inputs:
%   M: samples (in a form of a M x d matrix) or number of samples
%   lFfun: function of y (in form M x d) returning the exact log(PDF) and
%          optionally any other R functions at y in form of M x (R+1) matrix.
%          The first column should contain the exact log(PDF) values.
%   f: TT format of the PDF (tt_tensor) approximating exp(lFfun)
%   xsf: (d x 1) cell array of grid points (inc. boundaries) for all dimensions
%   correction: debiasing algorithm, 'mcmc' or 'iw'
%   [vec]: whether lFfun accepts vectorized input [true]
%
% Output:
%   y: M x d matrix of samples, distributed according to the exact PDF for
%      MCMC, and approximate PDF for IS
%   lFex: corrected lFfun sampled on y, in the form M x (R+1)
%
%   bias: a measure of bias in the proposed samples, that is:
%           if correction=='mcmc', this is the total number of rejections
%           if correction=='iw', this is the relative standard deviation of
%               the ratio Fex/Fapp for IW
%   time_invcdf: CPU time spent in IRT
%   worst: if correction=='mcmc', this is the distribution function of the
%               number of consequitive rejections
%          if correction=='iw', this is the maximal importance ratio
%   err1: [if correction=='iw' only!] empirical L1-error <|Fex-Fapp|>

function [y,lFex,bias,time_invcdf,worst,err1] = tt_irt_debias(M, lFfun, f, xsf, correction, vec)
d = f.d;
if (isscalar(M))
    Z = rand(M,d);
else
    Z = M;
    M = size(Z,1);
end
if (sum(cellfun(@length, xsf))~=sum(f.n))
    error('number of grid points in xsf should be n');
end

if (nargin<6)
    vec = true;
end

% Propose IRT samples
t_mcmc = tic;
[y,lFapp] = tt_irt(xsf, f, Z);
time_invcdf = toc(t_mcmc);
fprintf('invcdf time = %g\n', time_invcdf);

t_mcmc = tic;
if (vec)
    lFex = lFfun(y);
else
    % user function is not vectorised
    lFex = lFfun(y(1,:));
    b = numel(lFex);
    lFex = reshape(lFex, 1, b);
    lFex = [lFex; zeros(M-1, b)];
    for i=2:M
        lFex(i,:) = lFfun(y(i,:));
    end
end
fprintf('lFex time=%g\n', toc(t_mcmc));

if (strcmpi(correction, 'mcmc'))
    % Rejection loop here
    [y,lFex,~,bias,worst] = mcmc_prune(y,lFex,lFapp);    
    err1 = nan;
else
    % Importance Weighting: exact QoI = QoI*Fex/Fapp
    [lFex,bias,worst,err1] = iw_prune(lFex,lFapp);    
end
end

