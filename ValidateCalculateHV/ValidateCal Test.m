clear all; clc;

% ===================================================
% Helper Functions: Calculate Ratio and HV (from C++ code)
% ===================================================
function ratio = CalculateRatio(v1, v3)
    EPSILON = 1e-10;
    if abs(v3) > EPSILON
        ratio = 10 * v1 / v3;
    else
        ratio = 0.000;
    end
end

function hvValue = CalculateHV(v1, v3)
    EPSILON = 1e-10;

    % First get the ratio (C in the formula)
    ratio = CalculateRatio(v1, v3);

    % Check for valid ratio to prevent division by zero or invalid calculations
    if ratio < EPSILON
        hvValue = 0.000;
        return;
    end

    % term1 = (1 - (1 / (1 + exp((ratio - 3.4) / 0.4)))) * (ratio / 300.0)^(-4/3)
    term1 = (1 - (1 / (1 + exp((ratio - 3.4) / 0.4)))) * (ratio / 300.0)^(-4/3);

    % term2 = (1 / (1 + exp((ratio - 3.4) / 0.4))) * (ratio / 26.1)^(-20/7)
    term2 = (1 / (1 + exp((ratio - 3.4) / 0.4))) * (ratio / 26.1)^(-20/7);

    % Final HV = term1 + term2
    hvValue = term1 + term2;

    % Ensure non-negative value
    hvValue = max(hvValue, 0.000);
end

% Prompt user for file input (default: Tek717 TestCalculate.xlsx)
filename = input('Enter file name (default: Tek717 TestCalculate.xlsx): ', 's');

% Use default if empty
if isempty(filename)
    filename = 'Tek717 TestCalculate.xlsx';
    disp('Using default filename: Tek717 TestCalculate.xlsx');
end

% Check if filename already has extension
if endsWith(filename, '.xlsx') || endsWith(filename, '.csv')
    [~, name, ext] = fileparts(filename);
    if strcmp(ext, '.xlsx')
        file_xlsx = filename;
        file_csv = strcat(name, '.csv');
    else
        file_csv = filename;
        file_xlsx = strcat(name, '.xlsx');
    end
else
    file_xlsx = strcat(filename, '.xlsx');
    file_csv = strcat(filename, '.csv');
end

% Try to load CSV first (more reliable in Octave)
if exist(file_csv, 'file')
    try
        data = dlmread(file_csv, ',', 1, 0);  % Skip header row
        disp(['Loaded CSV file: ', file_csv]);
    catch
        try
            data = load(file_csv);
            disp(['Loaded CSV file: ', file_csv]);
        catch
            error(['Error: Could not load CSV file: ', file_csv]);
        end
    end
elseif exist(file_xlsx, 'file')
    try
        pkg load io
        [~, ~, data_cell] = xlsread(file_xlsx);
        data = cell2mat(data_cell(2:end, :));  % Skip header row
        disp(['Loaded Excel file: ', file_xlsx]);
    catch
        error(['Error: Could not load Excel file: ', file_xlsx]);
    end
else
    error('Error: Neither Excel nor CSV file exists.');
end

% Extract time and amplitude
time = data(:, 1);
amplitude = data(:, 2);

% Parameters
Fs = 50000;        % Sampling frequency (adjust if different)
nfft = 50;         % FFT points
overlap = floor(3 * nfft / 4); % 75% overlap
step_size = nfft - overlap;

% Constants from C++ code
NUMBER_OF_3KHz = 3000;
NUMBER_OF_9KHz = 9000;
compensation_factor = 2.38;  % 1.0 / 0.42 for Blackman window

% Create Blackman window
blackman_win = blackman(nfft);

% Initialize result arrays
time_segments = [];
magnitude_3k_array = [];
magnitude_9k_array = [];
V1_array = [];
V3_array = [];
Ratio_array = [];
HV_array = [];

% Calculate frequency indices
frequencies = (0:nfft/2) * (Fs / nfft);
[~, idx_3kHz] = min(abs(frequencies - NUMBER_OF_3KHz));
[~, idx_9kHz] = min(abs(frequencies - NUMBER_OF_9KHz));

disp(['3 kHz index: ', num2str(idx_3kHz), ' at frequency: ', num2str(frequencies(idx_3kHz)), ' Hz']);
disp(['9 kHz index: ', num2str(idx_9kHz), ' at frequency: ', num2str(frequencies(idx_9kHz)), ' Hz']);
disp(' ');
disp('Computing FFT with 75% overlap...');

% FFT with 75% overlap
segment_count = 0;
for i = 1:step_size:(length(amplitude) - nfft + 1)
    segment_count = segment_count + 1;

    % Extract segment
    segment = amplitude(i:i+nfft-1);

    % Apply Blackman window
    windowed_segment = segment .* blackman_win;

    % Compute FFT
    fft_result = fft(windowed_segment, nfft);

    % Get positive frequencies only
    fft_positive = fft_result(1:nfft/2+1);

    % Calculate magnitude at 3 kHz
    real_3k = real(fft_positive(idx_3kHz));
    imag_3k = imag(fft_positive(idx_3kHz));
    magnitude_3k = (2.0 / nfft) * sqrt(real_3k^2 + imag_3k^2);

    % Calculate magnitude at 9 kHz
    real_9k = real(fft_positive(idx_9kHz));
    imag_9k = imag(fft_positive(idx_9kHz));
    magnitude_9k = (2.0 / nfft) * sqrt(real_9k^2 + imag_9k^2);

    % Apply compensation factor (for Blackman window)
    V1 = magnitude_3k * compensation_factor;
    V3 = magnitude_9k * compensation_factor;

    % Calculate Ratio and HV (according to C++ code)
    Ratio = CalculateRatio(V1, V3);
    HV = CalculateHV(V1, V3);

    % Store results
    magnitude_3k_array = [magnitude_3k_array, magnitude_3k];
    magnitude_9k_array = [magnitude_9k_array, magnitude_9k];
    V1_array = [V1_array, V1];
    V3_array = [V3_array, V3];
    Ratio_array = [Ratio_array, Ratio];
    HV_array = [HV_array, HV];

    % Store midpoint time
    time_segments = [time_segments, time(i + floor(nfft / 2))];
end

disp(['Total segments processed: ', num2str(segment_count)]);
disp(' ');

% Prepare data for export
disp('Preparing data for export...');
export_data = [time_segments', magnitude_3k_array', magnitude_9k_array', V1_array', V3_array', Ratio_array', HV_array'];

% Create header
header = {'Time (s)', 'Magnitude_3kHz', 'Magnitude_9kHz', 'V1', 'V3', 'Ratio', 'HV'};

% Output filename
[~, base_name, ~] = fileparts(filename);
output_file = strcat(base_name, '_FFT_V1_V3_results.xlsx');

% Export to Excel
disp(['Exporting results to: ', output_file]);
try
    % Write header
    xlswrite(output_file, header, 'Sheet1', 'A1');

    % Write data
    xlswrite(output_file, export_data, 'Sheet1', 'A2');

    disp('Export successful!');
catch
    warning('Excel export failed. Trying CSV export...');
    output_csv = strcat(base_name, '_FFT_V1_V3_results.csv');

    % Combine header and data
    fid = fopen(output_csv, 'w');
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n', header{:});
    fclose(fid);

    % Append data
    dlmwrite(output_csv, export_data, '-append', 'delimiter', ',', 'precision', '%.6f');
    disp(['CSV export successful: ', output_csv]);
end

% Display summary statistics
disp(' ');
disp('=== Summary Statistics ===');
fprintf('Sampling frequency: %.0f Hz\n', Fs);
fprintf('FFT points: %d\n', nfft);
fprintf('Overlap: %.0f%%\n', (overlap/nfft)*100);
fprintf('Total segments: %d\n', segment_count);
fprintf('Compensation factor: %.2f\n', compensation_factor);
disp(' ');
fprintf('Magnitude 3 kHz - Mean: %.6f, Max: %.6f, Min: %.6f\n', ...
        mean(magnitude_3k_array), max(magnitude_3k_array), min(magnitude_3k_array));
fprintf('Magnitude 9 kHz - Mean: %.6f, Max: %.6f, Min: %.6f\n', ...
        mean(magnitude_9k_array), max(magnitude_9k_array), min(magnitude_9k_array));
fprintf('V1 (3 kHz)      - Mean: %.6f, Max: %.6f, Min: %.6f\n', ...
        mean(V1_array), max(V1_array), min(V1_array));
fprintf('V3 (9 kHz)      - Mean: %.6f, Max: %.6f, Min: %.6f\n', ...
        mean(V3_array), max(V3_array), min(V3_array));
fprintf('Ratio (10*V1/V3)- Mean: %.6f, Max: %.6f, Min: %.6f\n', ...
        mean(Ratio_array), max(Ratio_array), min(Ratio_array));
fprintf('HV Value        - Mean: %.6f, Max: %.6f, Min: %.6f\n', ...
        mean(HV_array), max(HV_array), min(HV_array));
disp('==========================');


