function R_n = estimate_noise_covariance(X, fs, stft_params, duration)
% estimate_noise_covariance Estimates noise spatial covariance matrix from initial silent frames.
%
% Inputs:
%   X             - STFT data (freq x time x channels)
%   fs            - Sampling frequency
%   stft_params   - Struct with field .window_length (e.g., 512)
%   duration      - Duration (in seconds) to assume as noise-only region
%
% Output:
%   R_n           - Noise covariance matrix (freq x channels x channels)

    [K, L, M] = size(X);  % K = freq bins, L = time frames, M = channels

    n_silent_time_bins = floor(duration * fs / stft_params.overlap);

    % Initialize covariance matrix
    R_n = zeros(K, M, M);

    for k = 1:K
        for l = 1:n_silent_time_bins
            x_k = squeeze(X(k, l, :));              % Mx1 vector
            R_n(k,:,:) = squeeze(R_n(k,:,:)) + (x_k * x_k');  % Outer product
        end
        R_n(k,:,:) = squeeze(R_n(k,:,:)) / n_silent_time_bins;
    end
end
