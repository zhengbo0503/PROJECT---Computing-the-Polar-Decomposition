function [X, H, its, delta, Rk] = poldec(A, TOL)
%POLDEC Polar decomposition.
% [U, H, ITS, DELTA, RK] = poldec(A) computes the polar decomposition 
% A = U*H of the square, nonsingular matrix A. ITS is the number of
% iterations for convergence. DELTA is the relative error between 
% successive iterations. RK is the norm(X'X - I) at each iteration.
[n,m] = size(A);
if n ~= m
    error('Input matrix is not square!');
end
if nargin == 2
    err = TOL;
else
    err = eps("double");
end
I = eye(n); term_tol = sqrt(2*err)*sqrt(n); maxiter = 1e6;
X = A; switched = false; 
delta = zeros(maxiter,1); Rk = zeros(maxiter,1);
for its = 1:1:maxiter
    % check if switch to Newton Schulz
    Rk(its) = norm(X'*X - I, inf);
    if Rk(its) <= 0.6
        switched = true;
    end
    % apply appropriate iterations
    if switched % Newton Schulz iteration
        tmp = X;
        X = 1.5*X - 0.5*X*(X'*X);
        delta(its+1) = norm(X - tmp, inf)/norm(X, inf); % store error
        % termination condition
        if delta(its+1) < term_tol || (delta(its+1) > delta(its)/(2) && its ~= 1)
            break;
        end
    else % Newton iteration
        tmp = X;
        X = 0.5*(inv(X)' + X);
        delta(its+1) = norm(X - tmp, inf)/norm(X, inf); % store error
    end
end
H = 0.5*(X'*A + A'*X); % construct H
end

