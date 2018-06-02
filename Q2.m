function [C, e1, e2, e3] = Q2(A, x) 
C = cond(A, 2)
b = A*x;

%Gauss Elimination
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

% LUx = b; Ly = b; Ux = y
y1(1, 1) = b(1)/L(1,1);
for i = 2:dim
    y1(i,1) = (b(i)-L(i,1:i-1)*y1(1:i-1,1))/L(i,i);
end

x1(dim,1) = y1(dim)/U(dim,dim);
for i = dim-1: -1:1
    x1(i,1) = (y1(i)-U(i, i+1:dim) * x1(i+1:dim,1))/U(i,i);
end
e1 = norm(x-x1,2)

% matlab LU
[L2, U2, P] = lu(A);
b2 = P*b;

% LUx = b; Ly = b; Ux = y
y2(1, 1) = b2(1)/L2(1,1);
for i = 2:dim
    y2(i,1) = (b2(i)-L2(i,1:i-1)*y2(1:i-1,1))/L2(i,i);
end

x2(dim,1) = y2(dim)/U2(dim,dim);
for i = dim-1: -1:1
    x2(i,1) = (y2(i)-U2(i, i+1:dim) * x2(i+1:dim,1))/U2(i,i);
end
e2 = norm(x-x2,2)

%full pvot
% size of A
dim = length(A);
% Initialize P and Q to the identity matrices.
P = eye(dim);
Q = eye(dim);


% For each row in A,
for i=1:dim-1,
    % Find the element with largest magnitude in each
    % submatrix which will be the new pivot.
    pivot = max(max(abs(A([i:dim],[i:dim]))));

 % find the indeces of the new pivot
 [c,m] = find(abs(A([i:dim],[i:dim])) == pivot);
 if i~=1;
    c(1) = c(1) + (i-1);
    m(1) = m(1) + (i-1);
 end;
 
% interchange the rows and columns of the new pivot
% with the old one
A([i,c(1)],:) = A([c(1),i],:);
A(:,[i,m(1)]) = A(:,[m(1),i]);

% store the changes in the matrices P and Q
P([i,c(1)],:) = P([c(1),i],:);
Q(:,[i,m(1)]) = Q(:,[m(1),i]);

% Compute the factor.
A(i+1:dim,i) = A(i+1:dim,i) / A(i,i);

% Multiply the nonzero elements of row i by the factor. Subtract this result from the "nonzero"
% elements of row j.
A(i+1:dim,i+1:dim) = A(i+1:dim,i+1:dim) - A(i,i+1:dim)*A(i+1:dim,i);
end


% The U factor is the upper triangle of A, with 0s in lower
U3 = triu(A);
% The L factor is the lower triangle of A, with 0s in ower
L3 = tril(A,-1) + eye(dim);

% LUQ'x = Pb; Ly =Pb; UQ'x = y; Q'x = z; x = Qz
b3 = P*b;
y3(1,1) = b3(1)/L3(1,1);
for i=2:dim
    y3(i,1) = (b3(i)-L3(i,1:i-1) *y3(1:i-1,1))/L3(i,i);
end
z3(dim,1) = y3(dim)/U3(dim-1,dim-1);

for i = dim-1:-1:1
    z3(i,1) = (y3(i)-U3(i,i+1:dim)*z3(i+1:dim,1))/U3(i,i);
end

x3 = Q*z3;
e3 = norm(x-x3,2)
