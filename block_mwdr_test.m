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
hs = { h_inter1, h_inter2, h_inter3, h_inter4, h_target};

%% Parameters
M = size(hs{1}, 1);                                         % Number of microphones
desired_SNR_db = 0;                                 % Desired SNR levels (dB)
% STFT parameters
stft_params.win_len = 512;                                  % Corresponds to 32ms (can be considered stationary)                                                   
stft_params.overlap = stft_params.win_len / 2;              % 50% overlap
stft_params.nfft = 512;                                     % 512 point fft
stft_params.win = sqrt(hann(stft_params.win_len, 'periodic'));  % To avoid edge artifacts

for index = 1:length(desired_SNR_db)
    SNR_db = desired_SNR_db(index);                             % Current SNR level

    % Data Generation Part
    N = min([length(s1), length(s2), length(n1), length(n2), length(n3)]);      % Number of samples for each signal
    s1 = s1(1:N); s2 = s2(1:N);                                                 % Clip source signals
    n1 = n1(1:N); n2 = n2(1:N); n3 = n3(1:N);                                   % Clip noise signals
    
    % Here we choose two target signals and the remaining signals as interferers.                                
    x = zeros(N, M);                                                    % Received signals
    x_noise_only = zeros(N, M);                                         % Received signals but only noise
    x_signal_only = zeros(N, M);                                        % Received signals but only clean signals
    noise_signals = {n1, n2, n3};                                       % Chosen noise signals
    
    for m = 1:M % Obtain received signals at each microphone
    
        % Noise components
        for noise_index = 1:length(noise_signals) % Loop through each noise component to add
            noise_signal = noise_signals{noise_index};                  % Noise signal
            room_response = squeeze(hs{noise_index}(m, :));             % Room response for current microphone and noise source
            noise_comp = conv(noise_signal, room_response, "same");     % Convolve with the room response                           % Mic signal
            x_noise_only(:,m) = x_noise_only(:, m) + noise_comp;
        end
        
        % Signal components
        clean_signal = s2;                                              % Clean signal
        room_response = squeeze(hs{5}(m, :));                           % Room response for current micrphone and signal source
        signal_comp = conv(clean_signal, room_response, "same");        % Convolve with the room response
        x_signal_only(:, m) = x_signal_only(:, m) + signal_comp;
        if m == 1
            s1_at_mic1  = signal_comp;
        end
    
        clean_signal = s1;                                              % Clean signal
        room_response = squeeze(hs{4}(m, :));                           % Room response for current micrphone and signal source
        signal_comp = conv(clean_signal, room_response, "same");        % Convolve with the room response
        x_signal_only(:, m) = x_signal_only(:, m) + signal_comp;
        if m == 1
            s2_at_mic1  = signal_comp;
        end
    end
    
    pad_length = stft_params.overlap - mod(N, stft_params.overlap);
    if pad_length > 0
        x_signal_only = [x_signal_only; zeros(pad_length, M)];
        x_noise_only = [x_noise_only; zeros(pad_length, M)];
        s1_at_mic1 = [s1_at_mic1; zeros(pad_length, 1)];
        s2_at_mic1 = [s2_at_mic1; zeros(pad_length, 1)];
    end
    
    % Scale noise to achieve target SNR
    x_noise_only_scaled = scale_noise_snr(x_signal_only, x_noise_only, SNR_db);
    
    % Final mixed signal
    x = x_signal_only + x_noise_only_scaled;
    
    % Take STFT
    [X, X_noise_only_scaled, X_signal_only, f, t] = compute_stfts(x, x_noise_only_scaled, x_signal_only, fs, ...
                                                    stft_params);
    [K, L] = size(X(:,:,1));  % frequency bins x time frames
    
    % Estimate Noise Spatial Covariance Cheating - Clean signal 2 does not
    % have an initial silent period
    R_n_cheat = estimate_noise_covariance_cheat(X_noise_only_scaled);
    
    % Estimate Received Signal Spatial Covariance Assuming Stationarity in a Band
    R_x = estimate_signal_covariance(X);

    % Estimate RTF for both signals
    A_s = zeros(K, 4, 2);
    for k = 1:K
        R_x_k = squeeze(R_x(k, :, :));
        R_n_k = squeeze(R_n_cheat(k, :, :));
        [U, D, Q] = gevd(R_x_k, R_n_k);    % GEVD: R_x_k * U = R_n_k * U * D
        A_s(k, :, :) = Q(:, 1:2) ./ Q(1, 1:2);  % Ensure 1st mic is reference
        Rs = Q(:, 1:2) * (D(1:2, 1:2) - 1) * Q(:, 1:2)';
    end
    
    % BLOCK MWDR
    for k = 1:K
        Rxk = squeeze(R_x(k,:,:));   % (M x M)
        ak = squeeze(A_s(k,:, :));          % (M x 2)
    
        % MVDR weights
        W_mvdr = inv(Rxk) * ak * inv(ak' *inv(Rxk) * ak);  % (M x 2)
    
        for l = 1:L
            x_k = squeeze(X(k,l,:));            % (M x 1)
            S_hat(k,l, :) = W_mvdr' * x_k;      % (K, L, 2)
        end
    end
    
    y1 = istft(S_hat(:, :, 1), fs, 'Window', stft_params.win, 'OverlapLength', ...
    stft_params.overlap, 'FFTLength', stft_params.nfft, ...
    'ConjugateSymmetric', true);
    y2 = istft(S_hat(:, :, 2), fs, 'Window', stft_params.win, 'OverlapLength', ...
    stft_params.overlap, 'FFTLength', stft_params.nfft, ...
    'ConjugateSymmetric', true);
    audiowrite('enhanced_mvdr_src1_.wav', normalize(y1(2000:end-5000), 'range', [-1, 1]), fs);
    audiowrite('enhanced_mvdr_src2_.wav', normalize(y2(2000:end-5000), 'range', [-1, 1]), fs);
    stoi_s1 = max(stoi(y1, s1_at_mic1, fs), stoi(y1, s2_at_mic1, fs))
    stoi_s2 = max(stoi(y2, s2_at_mic1, fs), stoi(y2, s1_at_mic1, fs))
end
