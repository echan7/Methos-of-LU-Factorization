% lu factor without pivoting
function [L, U] = lu_np(A)
% get size of row
dim = size(A, 1);
% initialize L to identity and fill lower triangle
L = eye(dim);

% each row
for i = 1: dim -1,
    % each row under row i
    for j = i + 1 : dim,
        
        % compute factor and store in L(j, i)
        L(j, i) = A(j,i) / A(i,i);
        % mutiply the nonzero elements of row i by factor. 
        % Subtract this result from nonzero elements of row j
        A(j,i+1:dim) = A(j,i+1:dim) - A(i,i+1:dim).*L(j,i);
        
    end
end

% U is upper triangle of A with 0s in lower
U = triu(A);