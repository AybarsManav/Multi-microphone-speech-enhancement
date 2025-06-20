%% Load Data
clc
clear all
close all

[s1, fs] = audioread('inputs/clean_speech.wav');                   % Source 1
[s2, ~]  = audioread('inputs/clean_speech_2.wav');                 % Source 2
[n2, ~]  = audioread('inputs/babble_noise.wav');                   % Babble noise
[n3, ~]  = audioread('inputs/Speech_shaped_noise.wav');            % Speech shaped noise
[n1, ~]  = audioread('inputs/aritificial_nonstat_noise.wav');      % Artificial non-stationary noise
load('inputs/impulse_responses.mat');                              % 5 impulses for 4 microphone each
hs = {h_inter1, h_inter2, h_inter3, h_inter4, h_target};


%% Clip signals
N = min([length(s1), length(s2), length(n1), length(n2), length(n3)]);      % Number of samples for each signal
s1 = s1(1:N); s2 = s2(1:N);                                                 % Clip source signals
n1 = n1(1:N); n2 = n2(1:N); n3 = n3(1:N);                                   % Clip noise signals

%% Obtain Received Signals at Microphones
% Here we choose a target signal and the remaining (possibly except clean signal 2) signals as interferers.
M = 4;                                                              % Number of microphones                                          
x = zeros(N, M);                                                    % Received signals
x_noise_only = zeros(N, M);                                         % Received signals but only noise
SNR_db = 0;                                                         % If desired, decide on noise levels
noise_signals = {n1, n2, n3};                                       % Chosen noise signals (can choose the other clean signal 2 as well)


for m = 1:M % Obtain received signals at each microphone

    % Noise components
    for noise_index = 1:length(noise_signals) % Loop through each noise component to add
        noise_signal = noise_signals{noise_index};                  % Noise signal
        room_response = squeeze(hs{noise_index}(m, :));             % Room response for current microphone and noise source
        noise_comp = conv(noise_signal, room_response, "same");     % Convolve with the room response
        x(:,m) = x(:, m) + noise_comp;                              % Mic signal
        x_noise_only(:,m) = x_noise_only(:, m) + noise_comp;
    end
    
    % Signal component
    clean_signal = s2;                                              % Clean signal
    room_response = squeeze(hs{5}(m, :));                           % Room response for current micrphone and signal source
    signal_comp = conv(clean_signal, room_response, "same");        % Convolve with the room response
    x(:,m) = x(:, m) + signal_comp;                                 % Mic signal

    clean_signal = s1;                                              % Clean signal
    room_response = squeeze(hs{4}(m, :));                           % Room response for current micrphone and signal source
    signal_comp = conv(clean_signal, room_response, "same");        % Convolve with the room response
    x(:,m) = x(:, m) + signal_comp;                
end
       
% p = play_scaled_sound(x(:, 1), fs); % Uncomment if you want to listen to the received signal at the first microphone
% stop(p); % If you want to stop playing

%% Take STFT
win_len = 512;                                                      % Corresponds to 32ms (can be considered stationary)                                                   
overlap = win_len / 2;                                              % 50% overlap
nfft = 512;                                                         % 512 point fft
win = sqrt(hann(win_len, 'periodic'));                              % To avoid edge artifacts

for m = 1:M
    [X(:,:,m), f, t] = stft(x(:,m), fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);
    [X_noise_only(:,:,m), ~, ~] = stft(x_noise_only(:,m), fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);
end
[K, L] = size(X(:,:,1));  % frequency bins x time frames
% plot_spectrograms(X, t, f);
%% Check clean signal's representation at mic4 to determine silent period
% figure;
% spectrogram(signal_comp, win, overlap, nfft, fs, 'yaxis');
% impixelinfo();

%% Estimate Noise Spatial Covariance Assuming Ergodicity
R_n = zeros(K, M, M);
n_silent_time_bins = floor(0.45 * fs / 512) ; % First 0.45 seconds is silent

% Rn is 512x4x4 (frequency bins(DFT) x channels x channels)
% estimating the noise covariance matrix based on the first n_silent_time_bins 
% frames for each frequency band
for k = 1:K
    for l = 1:n_silent_time_bins                            % First n_silent_time_bins frames
        x_k = squeeze(X(k,l,:));                            % Mx1
        R_n(k,:,:) = squeeze(R_n(k,:,:)) + (x_k * x_k');
    end
    R_n(k,:,:) = squeeze(R_n(k,:,:)) / n_silent_time_bins;
end

%% Estimate Noise Spatial Covariance Cheating
R_n_cheat = zeros(K, M, M);
for k = 1:K
    for l = 1:L
        x_k = squeeze(X_noise_only(k,l,:));                            % Mx1
        R_n_cheat(k,:,:) = squeeze(R_n_cheat(k,:,:)) + (x_k * x_k');
    end
    R_n_cheat(k,:,:) = squeeze(R_n_cheat(k,:,:)) / L;
end

%% Estimate Received Signal Spatial Covariance Assuming Stationarity in a Band
R_x = zeros(K, M, M);
% Rx is 512x4x4 (frequency bins(DFT) x channels x channels)
for k = 1:K
    for l = 1:L                                             
        x_k = squeeze(X(k,l,:));                            % Mx1
        R_x(k,:,:) = squeeze(R_x(k,:,:)) + (x_k * x_k');
    end
    R_x(k,:,:) = squeeze(R_x(k,:,:)) / L;
end

%% Obtain RTF (Ground truth)
A_s = zeros(K, M);
for m = 1:M
    A_m = fftshift(fft(hs{1}(m, :), K)); 
    A_s(:, m) = A_m;                 
end
A_s = A_s ./ A_s(:, 1);

%% Estimate RTF
A_hats = zeros(K, M, 2);                   % Estimated RTF per frequency band
R_s_hats = zeros(K, M, M);              % Estimated Rs per frequency band
for k = 1:K
    R_x_k = squeeze(R_x(k, :, :));
    R_n_k = squeeze(R_n_cheat(k, :, :));
    [U, D, Q] = gevd(R_x_k, R_n_k);        % Generalized eigenvalue decomp, colums of U are the right eigenvectors
    A_hats(k, :, :) = Q(:, 1:2) ./ Q(1, 1:2);   % Normalize with the first element
    R_s_hats(k, :, :) = Q(:, 1:2) * (D(1:2, 1:2) - 1) * Q(:, 1:2)';
end

%%
for k = 1:K
    Rxk = squeeze(R_x(k,:,:));   % (M x M)
    ak = squeeze(A_hats(k,:, :));          % (M x 1)

    % MVDR weights
    w_mvdr = inv(Rxk) * ak * inv(ak' *inv(Rxk) * ak);  % (M x 1)

    for l = 1:L
        x_k = squeeze(X(k,l,:));            % (M x 1)
        S_hat(k,l, :) = w_mvdr' * x_k;         % scalar
    end
end

y1 = istft(S_hat(:, :, 1), fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft, 'ConjugateSymmetric', true);
y2 = istft(S_hat(:, :, 2), fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft, 'ConjugateSymmetric', true);
audiowrite('enhanced_mvdr_src1_.wav', normalize(y1(50:end-50), 'range', [-1, 1]), fs);
audiowrite('enhanced_mvdr_src2_.wav', normalize(y2(50:end-50), 'range', [-1, 1]), fs);

