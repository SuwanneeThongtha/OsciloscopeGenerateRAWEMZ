% Generate RAWEMZ file from 64 Excel files
% Each Excel file contains time and signal columns
% This script:
% 1. Finds start time >= 0 for each channel
% 2. Truncates all channels to match the shortest length after time 0
% 3. Generates .rawemz file
% 4. Saves processed data as Excel file (64 rows × total_records columns)

clear all;
close all;
clc;

%% Configuration
NUM_CHANNELS = 64;
SAMPLING_RATE = 54000;  % Hz - Adjust according to your data
ODOMETER_DIAMETER = 56.0;  % mm - Adjust according to your setup
TOOL_SPEED = 1.0;  % m/s - Adjust according to your setup

%% Step 1: Select Excel files
fprintf('Please select 64 Excel files...\n');
[filenames, pathname] = uigetfile('*.xlsx', 'Select Excel Files', 'MultiSelect', 'on');

if ~iscell(filenames)
    if filenames == 0
        error('No files selected!');
    else
        filenames = {filenames};  % Convert to cell array if only one file
    end
end

if length(filenames) ~= NUM_CHANNELS
    error('You must select exactly %d Excel files! You selected %d files.', NUM_CHANNELS, length(filenames));
end

%% Step 2: Read all Excel files and find time >= 0
fprintf('Reading Excel files...\n');
all_data = cell(NUM_CHANNELS, 1);
start_indices = zeros(NUM_CHANNELS, 1);
lengths_after_zero = zeros(NUM_CHANNELS, 1);

for ch = 1:NUM_CHANNELS
    fprintf('Reading channel %d/%d: %s\n', ch, NUM_CHANNELS, filenames{ch});

    % Read Excel file (assuming columns: time, signal)
    filepath = fullfile(pathname, filenames{ch});
    data = xlsread(filepath);

    % Assuming first column is time, second column is signal
    time_col = data(:, 1);
    signal_col = data(:, 2);

    % Find first index where time >= 0
    idx_start = find(time_col >= 0, 1, 'first');

    if isempty(idx_start)
        error('Channel %d: No time values >= 0 found!', ch);
    end

    % Store data and indices
    all_data{ch}.time = time_col(idx_start:end);
    all_data{ch}.signal = signal_col(idx_start:end);
    start_indices(ch) = idx_start;
    lengths_after_zero(ch) = length(all_data{ch}.time);

    fprintf('  Start index: %d, Length after time 0: %d\n', idx_start, lengths_after_zero(ch));
end

%% Step 3: Find minimum length (shortest file after time 0)
min_length = min(lengths_after_zero);
fprintf('\nMinimum length after time 0: %d records\n', min_length);
fprintf('Truncating all channels to %d records...\n', min_length);

%% Step 4: Truncate all channels and create data matrix
% Matrix: 64 rows (channels) × min_length columns (records)
signal_matrix = zeros(NUM_CHANNELS, min_length);
time_vector = all_data{1}.time(1:min_length);  % Use time from first channel

for ch = 1:NUM_CHANNELS
    signal_matrix(ch, :) = all_data{ch}.signal(1:min_length);
end

%% Step 5: Normalize signals to int16 range [-32767, 32767]
fprintf('Normalizing signals...\n');
signal_matrix_normalized = signal_matrix;

% Find global min and max for normalization
global_min = min(signal_matrix(:));
global_max = max(signal_matrix(:));

fprintf('Signal range: [%.6f, %.6f]\n', global_min, global_max);

% Normalize to [-32767, 32767]
if global_max ~= global_min
    signal_matrix_normalized = ((signal_matrix - global_min) / (global_max - global_min) * 2 - 1) * 32767;
else
    signal_matrix_normalized = zeros(size(signal_matrix));
end

signal_matrix_int16 = int16(round(signal_matrix_normalized));

%% Step 6: Generate synthetic odometer data
fprintf('Generating odometer data...\n');
odometerCircumference = pi * ODOMETER_DIAMETER;  % mm
odometerSpeed = TOOL_SPEED * 1000;  % mm/s
rotationsPerSecond = odometerSpeed / odometerCircumference;

% Calculate time per sample
dt = 1.0 / SAMPLING_RATE;

odometer_count = zeros(3, min_length, 'uint32');
odometer_phase = zeros(3, min_length, 'uint16');

% Phase shifts for 3 odometers
phaseShifts = [0.0, odometerCircumference/4.0, odometerCircumference/2.0];

for i = 1:min_length
    t = (i-1) * dt;

    for odom = 1:3
        adjustedTime = t + phaseShifts(odom) / odometerSpeed;
        totalRotations = adjustedTime * rotationsPerSecond;
        revolutions = floor(totalRotations);
        rotationAngle = mod(totalRotations * 360.0, 360.0);
        odometerPhaseValue = mod(floor((rotationAngle / 360.0) * 4095), 4096);

        odometer_count(odom, i) = revolutions;
        odometer_phase(odom, i) = odometerPhaseValue;
    end
end

%% Step 7: Create FILEHEADER structure
fprintf('Creating file header...\n');

% Get current date/time
t = clock;

% Initialize header structure (will be written as binary)
header = struct();
header.FirmwareVersion = uint8(1);
header.FirmwareSubVersion = uint8(0);
header.FirmwareRevision = uint8(0);
header.FirmwareDate = uint8(t(3));  % Day
header.FirmwareMonth = uint8(t(2));  % Month
header.FirmwareYear = uint16(t(1));  % Year
header.OperationCode = 'RAW';  % 3 bytes
header.SamplingRateInHz = uint16(SAMPLING_RATE);
header.NumberOfVariableGroups = uint8(4);
header.OdometerDiameterInMM = single(ODOMETER_DIAMETER);
header.BodyDiameterInMM = single(254.0);
header.NumberOfOdometers = uint8(3);
header.SetupTime_SecondMM = int32(0);
header.HeaderSizeInByte = uint32(10240);
header.nOdometerType = uint8(0);
header.InternalPipeDiameterInMM = single(200.0);

%% Step 8: Ask for output filename
[output_file, output_path] = uiputfile('*.RAWEMZ', 'Save RAWEMZ File');
if output_file == 0
    error('No output file selected!');
end

output_filepath = fullfile(output_path, output_file);

%% Step 9: Write RAWEMZ file
fprintf('Writing RAWEMZ file: %s\n', output_filepath);

fid = fopen(output_filepath, 'wb');
if fid == -1
    error('Cannot create output file!');
end

try
    % Write header (simplified - first part only)
    fwrite(fid, header.FirmwareVersion, 'uint8');
    fwrite(fid, header.FirmwareSubVersion, 'uint8');
    fwrite(fid, header.FirmwareDate, 'uint8');
    fwrite(fid, header.FirmwareMonth, 'uint8');
    fwrite(fid, header.FirmwareYear, 'uint16');
    fwrite(fid, header.OperationCode, 'char');  % 3 bytes
    fwrite(fid, header.SamplingRateInHz, 'uint16');
    fwrite(fid, header.NumberOfVariableGroups, 'uint8');
    fwrite(fid, header.OdometerDiameterInMM, 'float32');
    fwrite(fid, header.BodyDiameterInMM, 'float32');
    fwrite(fid, header.NumberOfOdometers, 'uint8');
    fwrite(fid, header.SetupTime_SecondMM, 'int32');
    fwrite(fid, header.HeaderSizeInByte, 'uint32');
    fwrite(fid, header.FirmwareRevision, 'uint8');
    fwrite(fid, header.nOdometerType, 'uint8');
    fwrite(fid, header.InternalPipeDiameterInMM, 'float32');

    % Pad rest of header to 10240 bytes
    current_pos = ftell(fid);
    padding_size = 10240 - current_pos;
    if padding_size > 0
        fwrite(fid, zeros(padding_size, 1), 'uint8');
    end

    % Write records
    fprintf('Writing %d records...\n', min_length);

    for rec = 1:min_length
        % Counter (4 bytes, uint32)
        fwrite(fid, uint32(rec-1), 'uint32');

        % Signals (64 channels × 2 bytes = 128 bytes, int16)
        fwrite(fid, signal_matrix_int16(:, rec), 'int16');

        % Odometer Count (3 × 4 bytes, uint32)
        fwrite(fid, odometer_count(:, rec), 'uint32');

        % Odometer Phase (3 × 2 bytes, uint16)
        fwrite(fid, odometer_phase(:, rec), 'uint16');

        % Progress indicator
        if mod(rec, floor(min_length/10)) == 0
            fprintf('  Progress: %.1f%%\n', (rec/min_length)*100);
        end
    end

    fclose(fid);
    fprintf('RAWEMZ file created successfully!\n');

catch ME
    fclose(fid);
    rethrow(ME);
end

%% Step 10: Save processed data as Excel file
fprintf('Saving processed data to Excel...\n');

% Create Excel filename based on RAWEMZ filename
[~, base_name, ~] = fileparts(output_file);
excel_output = fullfile(output_path, [base_name '_processed_data.xlsx']);

% Prepare data for Excel: 64 rows × total_records columns
% Add channel labels in first column
excel_data = signal_matrix;  % Use original normalized data, not int16

% Create header row with record numbers
header_row = cell(1, min_length + 1);
header_row{1} = 'Channel';
for i = 1:min_length
    header_row{i+1} = sprintf('Record_%d', i);
end

% Create data with channel labels
excel_data_with_labels = cell(NUM_CHANNELS + 1, min_length + 1);
excel_data_with_labels(1, :) = header_row;

for ch = 1:NUM_CHANNELS
    excel_data_with_labels{ch+1, 1} = sprintf('CH%02d', ch);
    for rec = 1:min_length
        excel_data_with_labels{ch+1, rec+1} = signal_matrix(ch, rec);
    end
end

% Write to Excel
fprintf('Writing Excel file: %s\n', excel_output);
fprintf('Note: This may take a while for large datasets...\n');

xlswrite(excel_output, excel_data_with_labels);

fprintf('Excel file created successfully!\n');

%% Summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('Total channels: %d\n', NUM_CHANNELS);
fprintf('Total records per channel: %d\n', min_length);
fprintf('Sampling rate: %d Hz\n', SAMPLING_RATE);
fprintf('Total file size: %.2f MB\n', (10240 + min_length * 150) / 1024 / 1024);
fprintf('RAWEMZ file: %s\n', output_filepath);
fprintf('Excel file: %s\n', excel_output);
fprintf('=============================\n');

msgbox(sprintf('Conversion completed!\n\nRAWEMZ: %s\nExcel: %s', output_filepath, excel_output), 'Success');


