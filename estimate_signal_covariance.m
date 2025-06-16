function R_x = estimate_signal_covariance(X)
% estimate_signal_covariance Computes the spatial covariance of the mixture signal.
%
% Input:
%   X    - STFT of mixture signal (freq x time x channels)
%
% Output:
%   R_x  - Covariance matrix (freq x channels x channels)

    [K, L, M] = size(X);  % K: frequency bins, L: time frames, M: channels

    % Initialize
    R_x = zeros(K, M, M);

    % Loop over frequency bins
    for k = 1:K
        for l = 1:L
            x_k = squeeze(X(k,l,:));           % Mx1 vector
            R_x(k,:,:) = squeeze(R_x(k,:,:)) + (x_k * x_k');
        end
        R_x(k,:,:) = squeeze(R_x(k,:,:)) / L;
    end
end