function plot_spectrograms(S, t, f)
% S is [F x T x M]
% f is frequency vector (F x 1)
% t is time vector (1 x T)
% Plotting M spectrograms

M = size(S, 3);

% Compute global dB range for consistent color scaling
S_dB = 20 * log10(abs(S));
globalMax_dB = max(S_dB(:));
globalMin_dB = min(S_dB(:));  % 60 dB dynamic range

figure;
for m = 1:M
    subplot(ceil(M/2), 2, m);  % Use 4 columns, ceil(M/4) rows
    imagesc(t, f, S_dB(:, :, m));
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram %d', m));
    colorbar;
    colormap jet;
    caxis([globalMin_dB, globalMax_dB]);
end

end

