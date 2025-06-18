function [U, D, Q] = gevd(R_x_k, R_n_k)
% [U, D, Q] = eig(R_x_k, R_n_k);  % D is diagonal matrix of eigenvalues, U contains right eigenvectors
[U, D, Q] = eig(R_n_k \ R_x_k);

% Extract eigenvalues and sort in descending order
[eigVals, idx] = sort(diag(D), 'descend');

% Reorder eigenvalues and eigenvectors
D = diag(eigVals);
U = U(:, idx);
Q = Q(:, idx);

end

