function scaled_noise_signal = scale_signal_according_to_SNR(signal, noise, desired_snr_db)
    % Compute energies
    noise_energy = sum(noise.^2);
    signal_energy = sum(signal.^2);
    % Compute target energy
    target_noise_energy = signal_energy * 10^(-desired_snr_db * 0.1);
    % Scale the noise
    scaled_noise_signal = noise_comp / sqrt(noise_energy) * sqrt(target_noise_energy);
    % Checks
%     noise_energy_new = sum(noise_comp_scaled.^2);
%     SNR_db_now = 10 * log10(signal_energy / noise_energy_new)
end

