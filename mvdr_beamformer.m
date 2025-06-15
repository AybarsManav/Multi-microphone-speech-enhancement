function S_hat = mvdr_beamformer(X, R_n, A_hats)
% MVDR_BEAMFORMER Performs MVDR beamforming on multimicrophone input data.
%
% Inputs:
%   X       - Input signal (K x L x M), where:
%               K = number of frequency bins,
%               L = number of time frames,
%               M = number of microphones
%   R_n     - Noise covariance matrices (K x M x M)
%   A_hats  - RTF vectors (K x M)
%
% Output:
%   S_hat   - Beamformed output signal (K x L)

    [K, L, M] = size(X);
    S_hat = zeros(K, L);

    for k = 1:K
        Rnk = squeeze(R_n(k,:,:));   % (M x M)
        ak = A_hats(k,:).';          % (M x 1)

        % MVDR weights
        w_mvdr = (Rnk \ ak) / (ak' / Rnk * ak);  % (M x 1)

        for l = 1:L
            x_k = squeeze(X(k,l,:));            % (M x 1)
            S_hat(k,l) = w_mvdr' * x_k;         % scalar
        end
    end
end
