clear all; clc;

% Prompt user for file input
filename = input('Enter file name : ', 's');
file_xlsx = strcat(filename, '.xlsx');
file_csv = strcat(filename, '.csv');

% Try to load Excel, fallback to CSV
try
    data = readmatrix(file_xlsx);
    disp(['Loaded Excel file: ', file_xlsx]);
catch
    warning('Could not load Excel file. Trying CSV...');
    try
        data = readmatrix(file_csv);
        disp(['Loaded CSV file: ', file_csv]);
    catch
        error('Error: Neither Excel nor CSV file could be loaded.');
    end
end

% Extract time and amplitude
time = data(:, 1);
amplitude = data(:, 2);

% Parameters
Fs = 54000;        % Sampling frequency
nfft = 54;         % FFT points
overlap = floor(3 * nfft / 4); % 75% overlap
step_size = nfft - overlap;
filter_order = 50;  % FIR filter order

% Design FIR Band-Pass Filters
fir_3kHz  = fir1(filter_order, [2970 3020]   / (Fs / 2), 'bandpass');
fir_9kHz  = fir1(filter_order, [8970 9020]   / (Fs / 2), 'bandpass');
fir_6kHz  = fir1(filter_order, [5970 6020]   / (Fs / 2), 'bandpass');
fir_15kHz = fir1(filter_order, [14970 15020] / (Fs / 2), 'bandpass');

% Filter signals
filtered_3kHz  = filtfilt(fir_3kHz,  1, amplitude);
filtered_9kHz  = filtfilt(fir_9kHz,  1, amplitude);
filtered_6kHz  = filtfilt(fir_6kHz,  1, amplitude);
filtered_15kHz = filtfilt(fir_15kHz, 1, amplitude);

% Initialize result arrays
fft_3kHz_magnitudes  = [];
fft_9kHz_magnitudes  = [];
fft_6kHz_magnitudes  = [];
fft_15kHz_magnitudes = [];
time_segments        = [];

% Frequency vector and index for target bands
frequencies = (0:nfft/2) * (Fs / nfft);
[~, idx_3kHz]  = min(abs(frequencies - 3000));
[~, idx_9kHz]  = min(abs(frequencies - 9000));
[~, idx_6kHz] = min(abs(frequencies - 6000));
[~, idx_15kHz] = min(abs(frequencies - 15000));

% FFT with 75% overlap
for i = 1:step_size:(length(amplitude) - nfft + 1)
    % Extract segments
    seg3  = filtered_3kHz(i:i+nfft-1)  .* hamming(nfft);
    seg9  = filtered_9kHz(i:i+nfft-1)  .* hamming(nfft);
    seg6  = filtered_6kHz(i:i+nfft-1)  .* hamming(nfft);
    seg15 = filtered_15kHz(i:i+nfft-1) .* hamming(nfft);

    % FFT
    fft3  = fft(seg3,  nfft);
    fft9  = fft(seg9,  nfft);
    fft6  = fft(seg6, nfft);
    fft15 = fft(seg15, nfft);

    % Magnitudes (Volts)
    mag3  = abs(fft3(1:nfft/2+1))  * 2 / nfft;
    mag9  = abs(fft9(1:nfft/2+1))  * 2 / nfft;
    mag6 = abs(fft6(1:nfft/2+1)) * 2 / nfft;
    mag15 = abs(fft15(1:nfft/2+1)) * 2 / nfft;

    % Store peak magnitude values at each center frequency
    fft_3kHz_magnitudes  = [fft_3kHz_magnitudes,  mag3(idx_3kHz)];
    fft_9kHz_magnitudes  = [fft_9kHz_magnitudes,  mag9(idx_9kHz)];
    fft_6kHz_magnitudes = [fft_6kHz_magnitudes, mag6(idx_6kHz)];
    fft_15kHz_magnitudes = [fft_15kHz_magnitudes, mag15(idx_15kHz)];

    % Store midpoint time
    time_segments = [time_segments, time(i + floor(nfft / 2))];
end

% Harmonic magnitudes from FFT
V1 = fft_3kHz_magnitudes.*1000;
V3 = fft_9kHz_magnitudes.*1000;
V5 = fft_15kHz_magnitudes.*1000;

% Calculate feature
X = (V3 .^ 0.5) ./ (V1 .^ 0.68);
w = 1 ./ (1 + exp(-16 .* (X - 0.69)));

ratio = ((1 - w) .* (V3 .^ 0.5) + w .* ((V5.*3.2) .^ 0.5)) ./ ((V1 .^ 0.66)) * 10;

% --- Plotting Section ---
figure;

% 3 kHz amplitude
subplot(4, 1, 1);
plot(time_segments, fft_3kHz_magnitudes, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); 
ylabel('3 kHz Amplitude (V)');
title('Amplitude at 3 kHz vs Time');
grid on;

% 9 kHz amplitude
subplot(4, 1, 2);
plot(time_segments, fft_9kHz_magnitudes, 'g', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('9 kHz Amplitude (V)');
title('Amplitude at 9 kHz vs Time');
grid on;

% 15 kHz amplitude
subplot(4, 1, 3);
plot(time_segments, fft_15kHz_magnitudes, 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('15 kHz Amplitude (V)');
title('Amplitude at 15 kHz vs Time');
grid on;

% Ratio plot
subplot(4, 1, 4);
plot(time_segments, ratio, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Ratio (9 kHz / 3kHz^{1.5})');
title('Ratio of Amplitudes vs Time');
grid on;

% --- Add Logo as Overlay (Top-Left Corner) ---
ax_logo = axes('Position', [0.15 0.92 0.08 0.08]); % [left bottom width height]
imshow('Picture1.png');
axis off;
uistack(ax_logo, 'top');  % Ensure it stays on top of plots
