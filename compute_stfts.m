function [X, X_noise_only, X_signal_only, f, t] = compute_stfts(x, ...
    x_noise_only, x_signal_only, fs, stft_params)
% compute_stfts Computes the STFT of signal, noise-only, and signal-only
% inputs.
%
% Inputs:
%   x               - Matrix of mixed signals (samples x channels)
%   x_noise_only    - Matrix of noise-only signals (same size as x)
%   x_signal_only   - Matrix of signal-only signals (same size as x)
%   stft_params     - Struct with STFT parameters:
%                       .win_len   - Window length
%                       .overlap   - Overlap length
%                       .nfft      - FFT length
%                       .win       - Window function (vector of length win_len)
%
% Outputs:
%   X               - STFT of mixed signals (freq x time x channels)
%   X_noise_only    - STFT of noise-only signals
%   X_signal_only   - STFT of signal-only signals
%   f               - Frequency vector
%   t               - Time vector


    % Get dimensions
    [L, M] = size(x);

    % Estimate number of frequency bins and time frames
    num_freq_bins = stft_params.nfft;
    num_frames = floor((L - stft_params.overlap) / stft_params.overlap);
    X = zeros(num_freq_bins, num_frames, M);
    X_noise_only = zeros(num_freq_bins, num_frames, M);
    X_signal_only = zeros(num_freq_bins, num_frames, M);

    % Compute STFT for each microphone
    for m = 1:M
        [X(:,:,m), f, t] = stft(x(:,m), fs, ...
            'Window', stft_params.win, 'OverlapLength', stft_params.overlap, ...
            'FFTLength', stft_params.nfft);
        
        X_noise_only(:,:,m) = stft(x_noise_only(:,m), fs, ...
            'Window', stft_params.win, 'OverlapLength', stft_params.overlap, ...
            'FFTLength', stft_params.nfft);
        
        X_signal_only(:,:,m) = stft(x_signal_only(:,m), fs, ...
            'Window', stft_params.win, 'OverlapLength', stft_params.overlap, ...
            'FFTLength', stft_params.nfft);
    end
end
