function R_n_cheat = estimate_noise_covariance_cheat(X_noise_only_scaled)
% estimate_noise_covariance_cheat Computes the noise spatial covariance matrix
% using all available noise-only STFT frames (ideal / oracle case).
%
% Input:
%   X_noise_only_scaled - STFT of noise-only signal (freq x time x channels)
%
% Output:
%   R_n_cheat            - Noise covariance matrix (freq x channels x channels)

    [K, L, M] = size(X_noise_only_scaled);  % K: frequency bins, L: time frames, M: channels

    % Initialize
    R_n_cheat = zeros(K, M, M);

    % Estimate for each frequency bin
    for k = 1:K
        for l = 1:L
            x_k = squeeze(X_noise_only_scaled(k,l,:));        % Mx1 vector
            R_n_cheat(k,:,:) = squeeze(R_n_cheat(k,:,:)) + (x_k * x_k');
        end
        R_n_cheat(k,:,:) = squeeze(R_n_cheat(k,:,:)) / L;
    end
end
