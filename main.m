clc
clear
close all

[s1, fs] = audioread('clean_speech.wav');  
[s2, ~]  = audioread('aritificial_nonstat_noise.wav');  
h = load('impulse_responses.mat');  

N = min(length(s1), length(s2));
s1 = s1(1:N); s2 = s2(1:N);

M = 4;  
x = zeros(N-1+length(h.h_inter1(1,:)), M);  %convolution lenght

for m = 1:M
    h_s1 = squeeze(h.h_inter1(m,:));  % ir from target to mic m
    h_s2 = squeeze(h.h_inter2(m,:));  % ir from noise to mic m
    
    x(:,m) = conv(s1, h_s1) + conv(s2, h_s2);  % mic signal
end
x = x(1:N,:); 
%soundsc(x(:,1), fs);
win_len = 512;
overlap = win_len / 2;
nfft = 512;
win = sqrt(hann(win_len, 'periodic'));  

for m = 1:M
    [X(:,:,m), f, t] = stft(x(:,m), fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);
   
end
[K, L] = size(X(:,:,1));  % frequency bins x time frames

R_n = zeros(K, M, M);

%Rn is 512x4x4 (frequency bins(DFT) x channels x channels)
%estimating the noise covariance matrix based on the first 20 frames for
%each frequency band
for k = 1:K
    for l = 1:20  % first 20 frames - noise
        x_k = squeeze(X(k,l,:));  % Mx1
        R_n(k,:,:) = squeeze(R_n(k,:,:)) + (x_k * x_k');
    end
    R_n(k,:,:) = squeeze(R_n(k,:,:)) / 20;
end

%d = ones(K, M); no phase shifts

K = 512;
d = zeros(K, M);

for m = 1:M
    H_m = fft(h.h_inter1(m, :), K); 
    d(:, m) = H_m;                 
end

for k = 1:K
    Rnk = squeeze(R_n(k,:,:));
    dk = d(k,:).'; 
    
    % MVDR weights
    w_mvdr = (Rnk \ dk) / (dk' / Rnk * dk); 
    
    for l = 1:L
        x_k = squeeze(X(k,l,:));
        S_hat(k,l) = w_mvdr' * x_k;
    end
end


y = real(istft(S_hat, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft));
soundsc(y, fs);
audiowrite('enhanced_mvdr.wav', y, fs);
figure;
subplot(3,1,1); spectrogram(s1, win, overlap, nfft, fs, 'yaxis'); title('Clean speech');
subplot(3,1,2); spectrogram(x(:,1), win, overlap, nfft, fs, 'yaxis'); title('Noisy Mic 1');
subplot(3,1,3); spectrogram(y, win, overlap, nfft, fs, 'yaxis'); title('Enhanced output'); 
%%
s3  = audioread("aritificial_nonstat_noise.wav");
figure;
spectrogram(s3,win,overlap,nfft,fs)