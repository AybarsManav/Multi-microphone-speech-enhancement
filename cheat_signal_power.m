function S_pow = cheat_signal_power(x, fs, stft_params)
% cheat_signal_power Computes the sigma_s^2(k, l)
%
% Inputs:
%   x               - Reference microphone's clean and isolated signal (samples x 1)
%   stft_params     - Struct with STFT parameters:
%                       .win_len   - Window length
%                       .overlap   - Overlap length
%                       .nfft      - FFT length
%                       .win       - Window function (vector of length win_len)
%
% Outputs:
%   S_pow               - STFT of reference mic's signal (freq x time)

    % Compute STFT
    S_pow = stft(x, fs, ...
        'Window', stft_params.win, 'OverlapLength', stft_params.overlap, ...
        'FFTLength', stft_params.nfft);

    S_pow = abs(S_pow).^2;
end
