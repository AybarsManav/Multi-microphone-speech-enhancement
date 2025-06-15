function [U, D] = gevd(R_x_k, R_n_k)
[U, D] = eig(R_x_k, R_n_k);  % D is diagonal matrix of eigenvalues, U contains right eigenvectors

% Extract eigenvalues and sort in descending order
[eigVals, idx] = sort(diag(D), 'descend');

% Reorder eigenvalues and eigenvectors
D = diag(eigVals);
U = U(:, idx);
end

