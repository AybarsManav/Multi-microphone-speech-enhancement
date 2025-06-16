function scaled_noise = scale_noise_snr(clean_signal, noise_signal, target_snr_db)
% scale_noise_snr Scales noise to achieve a specific SNR with respect to a clean signal.
%
% Usage:
%   scaled_noise = scale_noise_snr(clean_signal, noise_signal, target_snr_db)
%
% Inputs:
%   clean_signal   - The original clean signal (vector)
%   noise_signal   - The noise signal to be scaled (vector, same length as clean_signal)
%   target_snr_db  - Desired SNR in decibels (dB)
%
% Output:
%   scaled_noise   - Scaled noise that, when added to clean_signal, gives the desired SNR

    % Ensure signals are the same length
    if length(clean_signal) ~= length(noise_signal)
        error('clean_signal and noise_signal must be the same length.');
    end

    % Calculate signal and noise powers
    signal_power = mean(clean_signal .^ 2);
    noise_power = mean(noise_signal .^ 2);

    % Compute target noise power based on desired SNR
    target_noise_power = signal_power / (10^(target_snr_db / 10));

    % Compute required scaling factor
    scaling_factor = sqrt(target_noise_power / noise_power);

    % Scale the noise
    scaled_noise = noise_signal * scaling_factor;
end
