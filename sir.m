%% Load Data
clc
clear all
close all

[s1, fs] = audioread('inputs/clean_speech.wav');                   % Source 1
[s2, ~]  = audioread('inputs/clean_speech_2.wav');                 % Source 2
[n1, ~]  = audioread('inputs/babble_noise.wav');                   % Babble noise
[n2, ~]  = audioread('inputs/Speech_shaped_noise.wav');            % Speech shaped noise
[n3, ~]  = audioread('inputs/aritificial_nonstat_noise.wav');      % Artificial non-stationary noise
load('inputs/impulse_responses.mat');                              % 5 impulses for 4 microphone each
hs = {h_inter1, h_inter2, h_inter3, h_inter4, h_target};
%%
% Visualize Room Responses
% plot_room_responses(hs, fs);

% Generate signals
M = size(hs{1}, 1);                                         % Number of microphones
desired_SNR_db = -50:10:20;                                 % Desired SNR levels (dB)
% STFT parameters
stft_params.win_len = 512;                                  % Corresponds to 32ms (can be considered stationary)                                                   
stft_params.overlap = stft_params.win_len / 2;              % 50% overlap
stft_params.nfft = 512;                                     % 512 point fft
stft_params.win = sqrt(hann(stft_params.win_len, 'periodic'));  % To avoid edge artifacts
% Noise covariance estimation 
silent_duration = 0.45;                                     % Consider first 0.45 seconds as silent period
snr_scores_ds = zeros(size(desired_SNR_db));
snr_scores_mwdr = zeros(size(desired_SNR_db));
snr_scores_mwdr_Rn_cheat = zeros(size(desired_SNR_db));
snr_scores_wiener = zeros(size(desired_SNR_db));
snr_scores_wiener_Rs_cheat = zeros(size(desired_SNR_db));

for index = 1:length(desired_SNR_db)

    SNR_db = desired_SNR_db(index);                             % Current SNR level
    
    [x, x_signal_only, x_noise_only_scaled, N] = generate_data(s1, s2, n1, n2, n3, hs, SNR_db, stft_params);
    
    % p = play_scaled_sound(x(:, 1), fs); % Uncomment if you want to listen to the received signal at the first microphone
    % stop(p); % If you want to stop playing
    
    [X, X_noise_only_scaled, X_signal_only, f, t] = compute_stfts(x, x_noise_only_scaled, x_signal_only, fs, ...
                                                    stft_params);
    
    [K, L] = size(X(:,:,1));  % frequency bins x time frames
    % plot_spectrograms(X, t, f); sgtitle("Received Signal Spectrograms");
    % plot_spectrograms(X_noise_only_scaled, t, f); sgtitle("Noise Spectrograms");
    % plot_spectrograms(X_signal_only, t, f); sgtitle("Clean Signal Spectrograms");
    
    % Estimate Noise Spatial Covariance Assuming Ergodicity
    R_n = estimate_noise_covariance(X, fs, stft_params, silent_duration);
    
    % Estimate Noise Spatial Covariance Cheating
    R_n_cheat = estimate_noise_covariance_cheat(X_noise_only_scaled);

    % Estimate Received Signal Spatial Covariance Assuming Stationarity in a Band
    R_x = estimate_signal_covariance(X);
    
    % Ground truth RTF
    A_ground_truth = compute_rtf_ground_truth(hs, K);
    
    % Estimate RTF
    R_s = zeros(K, 1, M, M);
    [A_s, R_s(:, 1, :, :)] = estimate_rtf_gevd(R_x, R_n);
    R_s = repmat(R_s, 1, L, 1, 1);

    % Estimate Signal Power Cheating
    R_s_cheat = zeros(K, L, 1, 1);  % Construct it as an Rs
    R_s_cheat(:, :, 1, 1) = cheat_signal_power(x_signal_only(:, 1), fs, stft_params); 
    
    % Delay and Sum Beamformer
    S_hat_delay_sum = ds_beamformer(X, A_s);
    y_delay_sum = istft(S_hat_delay_sum, fs, 'Window', stft_params.win, 'OverlapLength', ...
        stft_params.overlap, 'FFTLength', stft_params.nfft, ...
        'ConjugateSymmetric', true);

    % MVDR
    S_hat_mwdr = mvdr_beamformer(X, R_n, A_s);
    y_mwdr = istft(S_hat_mwdr, fs, 'Window', stft_params.win, 'OverlapLength', ...
        stft_params.overlap, 'FFTLength', stft_params.nfft, ...
        'ConjugateSymmetric', true);

    % MWDR Rn cheat
    S_hat_mwdr_cheat = mvdr_beamformer(X, R_n_cheat, A_s);
    y_mwdr_cheat = istft(S_hat_mwdr_cheat, fs, 'Window', stft_params.win, 'OverlapLength', ...
        stft_params.overlap, 'FFTLength', stft_params.nfft, ...
        'ConjugateSymmetric', true);

    % Wiener 
    S_hat_wiener = mvdr_beamformer(X, R_n, A_s);
    S_hat_wiener = single_channel_wiener(R_n, A_s, R_s, S_hat_wiener);
    y_wiener = istft(S_hat_wiener, fs, 'Window', stft_params.win, 'OverlapLength', ...
    stft_params.overlap, 'FFTLength', stft_params.nfft, ...
    'ConjugateSymmetric', true);
    
    % Wiener Rs Cheat
    S_hat_wiener_cheat = mvdr_beamformer(X, R_n, A_s);
    S_hat_wiener_cheat = single_channel_wiener(R_n, A_s, R_s_cheat, S_hat_wiener_cheat);
    y_wiener_cheat = istft(S_hat_wiener_cheat, fs, 'Window', stft_params.win, 'OverlapLength', ...
    stft_params.overlap, 'FFTLength', stft_params.nfft, ...
    'ConjugateSymmetric', true);
    clean = x_signal_only(:, 1);  % Reference clean signal
    
    snr_scores_mwdr(index) = 10 * log10(sum(clean.^2) / sum((clean - y_mwdr).^2));
    snr_scores_mwdr_Rn_cheat(index) = 10 * log10(sum(clean.^2) / sum((clean - y_mwdr_cheat).^2));
    snr_scores_wiener(index) = 10 * log10(sum(clean.^2) / sum((clean - y_wiener).^2));
    snr_scores_wiener_Rs_cheat(index) = 10 * log10(sum(clean.^2) / sum((clean - y_wiener_cheat).^2));
    snr_scores_ds(index) = 10 * log10(sum(clean.^2) / sum((clean - y_delay_sum).^2));


    % Results
    %audiowrite("outputs/enhanced_ds_snr_" + SNR_db + ".wav", normalize(y_delay_sum, 'range', [-1, 1]), fs);
    %audiowrite("outputs/enhanced_mvdr_snr_" + SNR_db + ".wav", normalize(y_mwdr, 'range', [-1, 1]), fs);
    %audiowrite("outputs/enhanced_mvdr_cheat_snr_" + SNR_db + ".wav", normalize(y_mwdr_cheat, 'range', [-1, 1]), fs);
    %audiowrite("outputs/enhanced_wiener_snr_" + SNR_db + ".wav", normalize(y_wiener, 'range', [-1, 1]), fs);
    %audiowrite("outputs/enhanced_wiener_Rs_cheat_snr_" + SNR_db + ".wav", normalize(y_wiener_cheat, 'range', [-1, 1]), fs);

    % figure;
    % subplot(3,1,1); spectrogram(x_signal_only(:, 1), stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Clean speech');
    % subplot(3,1,2); spectrogram(x(:,1), stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Noisy Mic 1');
    % subplot(3,1,3); spectrogram(y, stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Enhanced output');

end

%% Plots
figure;
hold on;
plot(desired_SNR_db, snr_scores_ds, '-x', 'DisplayName', 'DS');
plot(desired_SNR_db, snr_scores_mwdr, '-o', 'DisplayName', 'MWDR');
plot(desired_SNR_db, snr_scores_mwdr_Rn_cheat, '-s', 'DisplayName', 'MWDR + Rn cheat');
plot(desired_SNR_db, snr_scores_wiener, '-*', 'DisplayName', 'Wiener');
plot(desired_SNR_db, snr_scores_wiener_Rs_cheat, '-d', 'DisplayName', 'Wiener + Rs cheat');
hold off;

xlabel('Input SNR (dB)');
ylabel('Output SNR (dB)');
title('Comparison of Output SNRs');
legend('Location', 'best');
grid on;
