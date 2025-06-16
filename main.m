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
SNR_db = 0;                                                 % Desired SNR levels
% STFT parameters
stft_params.win_len = 512;                                  % Corresponds to 32ms (can be considered stationary)                                                   
stft_params.overlap = stft_params.win_len / 2;              % 50% overlap
stft_params.nfft = 512;                                     % 512 point fft
stft_params.win = sqrt(hann(stft_params.win_len, 'periodic'));  % To avoid edge artifacts
% Noise covariance estimation 
silent_duration = 0.45;                                     % Consider first 0.45 seconds as silent period


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
[A_s, R_s] = estimate_rtf_gevd(R_x, R_n);


%% MVDR
S_hat = mvdr_beamformer(X, R_n, A_s);
y = istft(S_hat, fs, 'Window', stft_params.win, 'OverlapLength', ...
    stft_params.overlap, 'FFTLength', stft_params.nfft, ...
    'ConjugateSymmetric', true);

% Results
p = play_scaled_sound(y, fs);
audiowrite('enhanced_mvdr.wav', normalize(y, 'range', [-1, 1]), fs);
figure;
% subplot(3,1,1); spectrogram(x_signal_only(:, 1), stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Clean speech');
% subplot(3,1,2); spectrogram(x(:,1), stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Noisy Mic 1');
% subplot(3,1,3); spectrogram(y, stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Enhanced output');
score = stoi(y, x_signal_only(:, 1), fs)


%% Multi Channel Wiener
% S_hat = mvdr_beamformer(X, R_n, A_s);
% S_hat = single_channel_wiener(R_n, A_s, R_s, S_hat);
% y = istft(S_hat, fs, 'Window', stft_params.win, 'OverlapLength', ...
%     stft_params.overlap, 'FFTLength', stft_params.nfft, ...
%     'ConjugateSymmetric', true);
% 
% % Results
% p = play_scaled_sound(y, fs);
% audiowrite('enhanced_wiener.wav', normalize(y, 'range', [-1, 1]), fs);
% figure;
% subplot(3,1,1); spectrogram(x_signal_only(:, 1), stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Clean speech');
% subplot(3,1,2); spectrogram(x(:,1), stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Noisy Mic 1');
% subplot(3,1,3); spectrogram(y, stft_params.win, stft_params.overlap, stft_params.nfft, fs, 'yaxis'); title('Enhanced output');
% score = stoi(y, x_signal_only(:, 1), fs)
