load('impulse_responses.mat');
[cs1, fs1] = audioread("clean_speech.wav");
[cs2, fs2] = audioread("clean_speech_2.wav");
h = {h_inter1, h_inter2};
s = {cs1, cs2};

%%
[x] = stft(h_inter1(1, :), fs1, "OverlapLength", 50, "FFTLength", 128 );
figure; imagesc(20 * log(1e-12 + abs(x)));
%%
d = 2;                              % Number of sources
M = 4;                              % Number of microphones
x = zeros(M, N);                    % Received signals in the microphones
N = min(length(cs1), length(cs2));  % Number of samples

% Clip the source signals
for i = 1:d
    s{i} = s{i}(1:N);
end

% Convolve with the channels
for j = 1:M
    for i = 1:d
        x(j, :) = x(j, :) + conv(s{i}, h{i}(j, :), "same")';
    end
end



%%
sound(2 * normalize(x{1}, "range") - 1, fs1);











