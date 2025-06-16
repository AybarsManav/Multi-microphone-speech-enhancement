function A_s = compute_rtf_ground_truth(hs, K)
% compute_rtf_ground_truth Computes the ground-truth Relative Transfer Function (RTF).
%
% Inputs:
%   hs  - Cell array of room impulse responses, where hs{5} is for the target source.
%         Each hs{i} is a matrix of size (M x T) for M microphones and T samples.
%   K   - Number of frequency bins (FFT length)
%
% Output:
%   A_s - Relative Transfer Function (RTF), size (K x M)

    M = size(hs{5}, 1);        % Number of microphones
    A_s = zeros(K, M);         % Preallocate RTF matrix

    for m = 1:M
        A_m = fftshift(fft(hs{5}(m, :), K));  % Frequency response for mic m
        A_s(:, m) = A_m;                      % Store full FFT result
    end

    % Normalize with respect to reference mic (mic 1)
    A_s = A_s ./ A_s(:, 1);
end
