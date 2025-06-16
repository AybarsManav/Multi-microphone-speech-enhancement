function R_s_hat = estimate_signal_covariance_cheat(X_signal_only, A_s)
% estimate_r_s_cheat Estimates rank-1 spatial covariance of the target signal using clean STFT.
%
% Assumes reference microphone is always mic 1.
%
% Inputs:
%   X_signal_only - Clean target STFT (freq x time x channels)
%   A_s           - Ground-truth RTFs (freq x channels), normalized by mic 1
%
% Output:
%   R_s_hat       - Rank-1 spatial covariances (freq x time x channels x channels)

    [K, L, M] = size(X_signal_only);

    R_s_hat = zeros(K, L, M, M);

    for k = 1:K
        a_k = A_s(k, :).';  % M x 1 RTF for frequency bin k

        for l = 1:L
            % Magnitude squared of the clean STFT at reference mic (mic 1)
            sigma_s_kl = abs(X_signal_only(k, l, 1))^2;

            % Rank-1 covariance estimate: σ_s * a_k * a_k^H
            R_s_hat(k, l, :, :) = sigma_s_kl * (a_k * a_k');
        end
    end
end
