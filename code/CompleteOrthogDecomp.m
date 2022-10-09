%% Complete Orthogonal Decomposition

clc; clear; close all; 
n = 100;
A = randn(n);

% Deliberately making singular matrices
[P,S,V] = svd(A);
for i = 1:randi([1,n])
    entry = randi([1,n]);
    S(entry,entry) = 0;
end
A_new = P * S * V';

if rank(A_new) == n 
    fprintf("Full Rank! No need <strong>Complete orthogonal " + ...
        "decomposition</strong>.\n");
    return;
end
% QR by column pivoting
[Q,R] = qr(A_new);

% Householder
K = rank(A_new) + 1;
for i = K : 1 : n
    R = house_elim(R, i);
end

% Householder on first row
R = eliminate(R, K - 1);

% remove small elements
tol = norm(A) * eps;
for i = 1:n
    for j = 1:n
        if abs(R(i,j)) < tol
            R(i,j) = 0;
        end
    end
end

disp(R);

%% Supplimentary functions
function A = house_elim(A,i) 
    x = A(:,i);
    n = length(x); sigma = x(2:n)'*x(2:n); 
    mu = sqrt(x(1)*x(1) + sigma);
    if x(1) >= 0
        x(1) = x(1) + mu;
    else
        x(1) = -sigma/(x(1) + mu);
    end
    b = 2/(x(1)*x(1) + sigma); 
    P = eye(n) - b * (x * x'); 
    A(:,[i,n]) = P * A(:,[i,n]);
end

function [x,b] = house(x)
n = length(x); sigma = x(2:n)'*x(2:n); 
if sigma == 0, b = 2/(x(1)*x(1)); 
else
    mu = sqrt(x(1)*x(1) + sigma);
    if x(1) >= 0, x(1) = x(1) + mu;
    else, x(1) = -sigma/(x(1) + mu); end
    b = 2/(x(1)*x(1) + sigma); 
end
end

function R = eliminate(R, rank1)
    % rank1 is the rank of A
    R = R';
    x = R(rank1:end,1);
    n = length(x);
    [v,b] = house(x);
    P = eye(n) - b * (v * v');
    x = P * x;
    R(rank1:end,1) = x;
    R = R';
end















