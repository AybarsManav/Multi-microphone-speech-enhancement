function S_hat = single_channel_wiener(R_n, A_hats, R_s_hats, S_hat)
% SINGLE_CHANNEL_WIENER Applies the single-channel Wiener filter
%
% Inputs:
%   R_n       - Noise covariance matrices (K x M x M)
%   A_hats    - RTF corresponding to the target source (K x M)
%   R_s_hats  - Source signal covariance matrices (K x L x M x M)
%   S_hat     - Initial estimated source power spectrogram (K x L)
%
% Output:
%   S_hat     - Updated source power spectrogram after Wiener filtering

[K, M, ~] = size(R_n);
[~, L] = size(S_hat);

for k = 1:K
    Rnk = squeeze(R_n(k,:,:));          % (M x M)
    ak = A_hats(k,:).';                 % (M x 1)

    for l = 1:L
        Rskl = squeeze(R_s_hats(k,l,:,:));     % (M x M)
        sigma2_s_k = real(Rskl(1, 1));       % scalar
        denom = sigma2_s_k + 1 / real(ak' / Rnk * ak);
        S_hat(k,l) = sigma2_s_k / denom * S_hat(k,l);  % scalar
    end
end
end
