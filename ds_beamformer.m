function S_hat = ds_beamformer(X, A_hats)
% DS_BEAMFORMER Performs Delay-and-Sum beamforming on multimicrophone input data.
%
% Inputs:
%   X       - Input signal (K x L x M), where:
%               K = number of frequency bins,
%               L = number of time frames,
%               M = number of microphones
%   A_hats  - RTF vectors / steering vectors (K x M)
%
% Output:
%   S_hat   - Beamformed output signal (K x L)

    [K, L, M] = size(X);
    S_hat = zeros(K, L);

    for k = 1:K
        ak = A_hats(k,:).';             % (M x 1)
        w_ds = ak ./ (ak' * ak);

        for l = 1:L
            x_k = squeeze(X(k,l,:));    % (M x 1)
            S_hat(k,l) = w_ds' * x_k;   % scalar
        end
    end
end
