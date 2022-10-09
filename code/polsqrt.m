function H = polsqrt(A)
%POLSQRT Matrix square root using polar decomposition.
% H = poldec(A) computes the matrix square root of the symmetric positive
% definite matrix A, such that norm(H^2 - A) is small.
R = chol(A);
[~,H] = poldec(R);
end
