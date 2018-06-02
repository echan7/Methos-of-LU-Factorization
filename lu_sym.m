% LU factor algo for sym matrix A, non singular
function [L, U] = lu_sym(A)

% get size of row
dim = size(A, 1);
% initialize L to identity and fill lower triangle
L = eye(dim);

% each row
for i = 1: dim -1,
    % each row under row i
    for j = i + 1 : dim,
        
        % compute factor and store in L(j, i)
        L(j, i) = A(i,j) / A(i,i);
        % mutiply the upper triangular part of row i by the factor
        % subtract result from upper traingular of row j
        A(j,j:dim) = A(j,j:dim) - A(i,j:dim).*L(j,i);
        
    end
end

% U is upper triangle of A with 0s in lower
U = triu(A);
        
