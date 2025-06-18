function [A_hats, R_s_hats] = estimate_rtf_gevd(R_x, R_n)
% estimate_rtf_gevd Estimates the Relative Transfer Function (RTF) using GEVD.
%
% Inputs:
%   R_x     - Mixture covariance matrix (K x M x M)
%   R_n     - Noise covariance matrix (K x M x M)
%
% Outputs:
%   A_hats     - Estimated RTFs (K x M), one per frequency bin
%   R_s_hats   - Estimated signal covariance matrices (K x M x M)

    [K, M, ~] = size(R_x);

    A_hats = zeros(K, M);         % RTF estimates
    R_s_hats = zeros(K, M, M);    % Estimated Rs (signal covariances)

    for k = 1:K
        R_x_k = squeeze(R_x(k, :, :));
        R_n_k = squeeze(R_n(k, :, :));

        [U, D, Q] = gevd(R_x_k, R_n_k);     % GEVD: R_x_k * U = R_n_k * U * D

        % Use dominant eigenvector and normalize
        A_hats(k, :) = Q(:, 1).' / Q(1, 1);  % Ensure 1st mic is reference

        % Estimate signal covariance matrix
        lambda_1 = D(1, 1);
        R_s_hats(k, :, :) = (lambda_1 - 1) * Q(:, 1) * Q(:, 1)' / ...
            (Q(:, 1)' / R_n_k * Q(:, 1));


    end
end
