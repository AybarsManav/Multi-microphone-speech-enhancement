function [x, x_signal_only, x_noise_only_scaled, N] = generate_data(s1, s2, n1, n2, n3, hs, SNR_db, stft_params)
% generate_data Clips signals, simulates microphone recordings, and scales noise.
%
% Inputs:
%   s1, s2         - Clean source signals
%   n1, n2, n3     - Noise signals
%   hs             - Room impulse responses as a 1x5 cell array
%                    (each cell is MxT where M is num. of mics, T is impulse length)
%   SNR_db         - Desired Signal-to-Noise Ratio (in dB)
%   stft_params     - Struct with STFT parameters:
%                       .win_len   - Window length
%                       .overlap   - Overlap length
%                       .nfft      - FFT length
%                       .win       - Window function (vector of length win_len)
%
% Outputs:
%   x              - Mixed microphone signals (signal + scaled noise), size N x M
%   x_signal_only  - Target signal only at each mic, size N x M
%   x_noise_only_scaled - Scaled noise at each mic, size N x M
%   N              - Length of clipped signals

    % --- Clip signals to common length ---
    N = min([length(s1), length(s2), length(n1), length(n2), length(n3)]);
    s1 = s1(1:N); s2 = s2(1:N);
    n1 = n1(1:N); n2 = n2(1:N); n3 = n3(1:N);

    % --- Initialize ---
    M = size(hs{1}, 1);  % Number of microphones (assumes hs{i} is MxT)
    x_signal_only = zeros(N, M);
    x_noise_only = zeros(N, M);
    noise_signals = {n1, n2, n3};

    % --- Convolve signals with room responses ---
    for m = 1:M
        % Add noise components from multiple noise sources
        for noise_index = 1:length(noise_signals)
            noise_signal = noise_signals{noise_index};
            room_response = squeeze(hs{noise_index}(m, :));
            noise_comp = conv(noise_signal, room_response, 'same');
            x_noise_only(:, m) = x_noise_only(:, m) + noise_comp;
        end

        % Add signal component from the target source
        room_response = squeeze(hs{5}(m, :));  % hs{5} is target signal room response
        signal_comp = conv(s1, room_response, 'same');
        x_signal_only(:, m) = signal_comp;
    end

    % Pad signals to be integer multiple of overlap
    pad_length = stft_params.overlap - mod(N, stft_params.overlap);
    if pad_length > 0
        x_signal_only = [x_signal_only; zeros(pad_length, M)];
        x_noise_only = [x_noise_only; zeros(pad_length, M)];
    end

    % --- Scale noise to achieve target SNR ---
    x_noise_only_scaled = scale_noise_snr(x_signal_only, x_noise_only, SNR_db);

    % --- Final mixed signal ---
    x = x_signal_only + x_noise_only_scaled;
end
