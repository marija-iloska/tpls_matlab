function [psd, frequencies] = powerSD(feature, t)

% Extract the time series data up to time t
hj = sum(feature(1:t, :), 2);

% Compute the FFT
fft_result = fft(hj);

% Compute the Power Spectral Density
psd = abs(fft_result).^2;

% Frequency axis for plotting
fs = 1; % Assuming unit sample rate for simplicity
frequencies = linspace(0, fs/2, length(psd)/2 + 1);

end