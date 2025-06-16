function plot_room_responses(hs, fs)
% plot_room_responses(hs, fs)
% Plots impulse responses for multiple locations and microphones.
%
% Parameters:
%   hs - a 1x5 cell array where each cell contains a 4xN matrix (4 mics x N samples)
%   fs - sampling frequency

    % Time vector
    t = 1:length(hs{1}(1, :));
    t = t / fs;

    % Determine global min and max for y-axis
    ymin = inf;
    ymax = -inf;
    for i = 1:5
        for j = 1:4
            data = hs{i}(j, :);
            ymin = min(ymin, min(data));
            ymax = max(ymax, max(data));
        end
    end

    % Plot
    figure;
    for i = 1:5  % Room responses (rows)
        for j = 1:4  % Microphones (columns)
            idx = (i - 1) * 4 + j;  % Convert (i,j) to subplot index
            subplot(5, 4, idx);
            plot(t, hs{i}(j, :));
            title(sprintf('Location %d - Mic %d', i, j));
            xlabel('Time (s)');
            ylabel('Amplitude');
            ylim([ymin ymax]);
        end
    end
end
