% Constructs Q1 FEM discretization of the fruitfly problem and generates
% cosine KLE
function [p,bound,W1g,W1m,spind,Pua,phi,lambda,Mass] = build_grid_and_kle(meshlevel, bc_type, nu, corr_length, tol_kle, m0)
%% Initialize Q1 FEM
h = 2^(-4-meshlevel);
n = 2^(4+meshlevel)+1;
% Gridpoints
p = (0:h:1)';
x = repmat(p,1,n);
y = repmat(p',n,1);
p = [x(:)'; y(:)'];

if (strcmpi(bc_type, 'dn'))
    % Dirichlet-Neumann
    % Indices of Dirichlet boundary nodes. Make it separately left & right
    bound = [find(p(1,:)==0)'; find(p(1,:)==1)'];
else
    % Full Dirichlet
    bound = [find(p(1,:)==0)'; find(p(1,:)==1)'; find(p(2,:)==0)'; find(p(2,:)==1)'];
    bound = unique(bound, 'stable');
end

% Here we want the same representation for a and u, hence p==t (except
% the BC).
% Prepare a sparse 3D tensor for creating a matrix out of the coeff.
%%%%%%%%%%%%%%%% 1D gradient
W1g = sparse(n^2, n);
i0 = (1:n)';
% i,i,i-1
i = i0;
j = i;
k = i-1;
i(k<1)=[];
j(k<1)=[];
k(k<1)=[];
wijk = (0.5/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i-1,i-1
i = i0;
j = i-1;
k = i-1;
i(k<1)=[];
j(k<1)=[];
k(k<1)=[];
wijk = (-0.5/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i
i = i0;
j = i;
k = i;
wijk = (1/h)*ones(numel(i),1);
wijk(1) = 0.5/h; % BC
wijk(n) = 0.5/h; % BC
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i-1,i
i = i0;
j = i-1;
k = i;
i(j<1)=[];
k(j<1)=[];
j(j<1)=[];
wijk = (-0.5/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i+1,i
i = i0;
j = i+1;
k = i;
i(j>n)=[];
k(j>n)=[];
j(j>n)=[];
wijk = (-0.5/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i+1
i = i0;
j = i;
k = i+1;
i(k>n)=[];
j(k>n)=[];
k(k>n)=[];
wijk = (0.5/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i+1,i+1
i = i0;
j = i+1;
k = i+1;
i(k>n)=[];
j(k>n)=[];
k(k>n)=[];
wijk = (-0.5/h)*ones(numel(i),1);
W1g = W1g + sparse(i+(j-1)*n,k,wijk,n^2,n);

%%%%%%%%%%%%%%% 1D Mass
W1m = sparse(n^2, n);
i0 = (1:n)';
% i,i,i-1
i = i0;
j = i;
k = i-1;
i(k<1)=[];
j(k<1)=[];
k(k<1)=[];
wijk = (h/12)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i-1,i-1
i = i0;
j = i-1;
k = i-1;
i(k<1)=[];
j(k<1)=[];
k(k<1)=[];
wijk = (h/12)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i
i = i0;
j = i;
k = i;
wijk = (h/2)*ones(numel(i),1);
wijk(1) = h/4; % BC
wijk(n) = h/4; % BC
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i-1,i
i = i0;
j = i-1;
k = i;
i(j<1)=[];
k(j<1)=[];
j(j<1)=[];
wijk = (h/12)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i+1,i
i = i0;
j = i+1;
k = i;
i(j>n)=[];
k(j>n)=[];
j(j>n)=[];
wijk = (h/12)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i,i+1
i = i0;
j = i;
k = i+1;
i(k>n)=[];
j(k>n)=[];
k(k>n)=[];
wijk = (h/12)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);
% i,i+1,i+1
i = i0;
j = i+1;
k = i+1;
i(k>n)=[];
j(k>n)=[];
k(k>n)=[];
wijk = (h/12)*ones(numel(i),1);
W1m = W1m + sparse(i+(j-1)*n,k,wijk,n^2,n);

% Unfortunately, we have to return W1g and W1m separately, since Matlab
% is so stupid that performs sparse*sparse somehow with the full memory
% We will assemble W2g = W1g*C*W1m' + W1m*C*W1g';
% However, we need to perform [1324] permute afterwards. Prepare
% indices.
A = W1g*sparse(ones(n,n))*W1g'; % We only need indices, any matrix is ok
[ij1,ij2]=find(A);
j1 = floor((ij1-1)/n)+1;
i1 = ij1-(j1-1)*n;
j2 = floor((ij2-1)/n)+1;
i2 = ij2-(j2-1)*n;
i = i1+(i2-1)*n;
i = int64(i);
j = j1+(j2-1)*n;
j = int64(j);
ij1 = int64(ij1); % Otherwise we will overflow double
ij2 = int64(ij2);
% Now the permute is computed as follows:
% A_real = sparse(i,j,A(ij1+(ij2-1)*n^2),n^2,n^2);
spind = [i, j, ij1+(ij2-int64(1))*int64(n)^2]; % keep them in the same storage

% Projection from u to inner u and A
Pua = speye(n^2);
Pua(bound,:)=[];


%% KLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the KLE C0*(k-startpos)^(-nu)*cos(2*pi*f1(k)*x)*cos(2*pi*f2(k)*y),
% where f1(k) = k - s(k)*(s(k)+1)/2, f2(k) = s(k)-f1(k),
% s(k) = floor(-1/2 + sqrt(1/4+2*k)),
% and C0 is chosen such that sum(lambda)=sigma

startpos = 1/corr_length - 1;

fprintf('Computing cosine KLE\n');
% Generate a redundant nu-decay
if (isinf(nu))
    L = ceil(-log2(tol_kle*0.1));
    L = min(L, size(p,2));
    ind = (1:L)';
    lambda = [ones(1,startpos), 2.^(-ind)];
else
    L = ceil(exp(-log(tol_kle*0.1)/(nu+1)));
    L = min(L, size(p,2));
    ind = (1:L)';
    lambda = [ones(1,startpos), ind.^(-nu-1)];
end
% Normalize the variance
lambda = lambda/sum(lambda);

% Actually truncate lambda according to tol_kle
L = min([find(lambda<tol_kle*lambda(1), 1), L]);
lambda = lambda(1:L);
ind = 1:L;

% Map ind to two-dimensional orders
s = floor(-1/2+sqrt(1/4+2*ind));
f1 = ind - s.*(s+1)*0.5;
f2 = s-f1;
phi = cos(2*pi*p(1,:)'*f1).*cos(2*pi*p(2,:)'*f2);



%% Windowed mass matrices for pressure observations
% In 1D first, since the domain is symmetric
Mass1 = cell(m0, 1);
for i=1:m0
    % our interval is [(i-1)H, (i+1)H], where H=1/(m0+1)
    ind = double(((0:(p(1,2)-p(1,1)):1)>=(i-1)/(m0+1))&((0:(p(1,2)-p(1,1)):1)<=(i+1)/(m0+1)))/(0.5/(m0+1));
    Mass1{i,1} = W1m*sparse(ind(:));
    Mass1{i,1} = reshape(Mass1{i,1}, size(W1m,2), []);
    % Correct the border elements wrongly set
    ind_l = find(ind,1,'first');
    if (ind_l>1)
        Mass1{i,1}(ind_l-1,:)=0;
        Mass1{i,1}(:,ind_l-1)=0;
        Mass1{i,1}(ind_l, ind_l) = Mass1{i,1}(ind_l+1, ind_l+1)*0.5;
    end
    ind_r = find(ind,1,'last');
    if (ind_r<size(W1m,2))
        Mass1{i,1}(ind_r+1,:)=0;
        Mass1{i,1}(:,ind_r+1)=0;
        Mass1{i,1}(ind_r, ind_r) = Mass1{i,1}(ind_r-1, ind_r-1)*0.5;
    end
end
% 2D Masses
Mass = cell(m0, m0);
for j=1:m0
    for i=1:m0
        Mass{i,j} = kron(Mass1{j}, Mass1{i});
    end
end
end
