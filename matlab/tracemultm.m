% Computes "trace" matrix product and indexing operations
% C(i,:) = A(i,:)*B(:,:,j(i)) and
% C(i) = A(i,j(i))            (if B is missing)
% Inputs:
%   A,j and (optional) B
% Outputs:
%   C
function [C] = tracemultm(A,j,B)
n = size(A,1);
if (nargin<3)||(isempty(B))    
    C = zeros(n,1);
    for i=1:n
        C(i) = A(i,j(i));
    end
else
    m = size(B,2);
    C = zeros(n,m);
    for i=1:n
        C(i,:) = A(i,:)*B(:,:,j(i));
    end
end
end
