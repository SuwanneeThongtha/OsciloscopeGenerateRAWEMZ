clear all; clc;

function GenerateRAWEMZ(data_folder_or_pattern, output_filename)
% GENERATE_RAWEMZ_FROM_EXCEL - Generate .rawemz file from Excel files
%
% Usage: GenerateRAWEMZ('DataFile', 'output.rawemz')
%        GenerateRAWEMZ('DataFile/', 'output.rawemz')
%        GenerateRAWEMZ('*.xlsx', 'output.rawemz')
%        GenerateRAWEMZ({'file1.xlsx', 'file2.xlsx', ...}, 'output.rawemz')
%
% Parameters:
%   data_folder_or_pattern - Folder path containing Excel files, file pattern, or cell array of filenames
%   output_filename        - Output .rawemz filename
%
% Requirements:
%   - Each Excel file should have columns: Time, Signal
%   - Time should be in seconds
%   - Signal values should be numeric
%   - Expects up to 96 Excel files (will pad with zeros if fewer)
%
% Behavior:
%   - Auto-detects all Excel files in specified folder/pattern
%   - For each file, keep rows where time >= 0 (first non-negative onward)
%   - Use the shortest record count (after time>=0) among all files
%   - NO resampling/interpolation; simply truncate every channel to min length
%   - Creates dummy odometer data
%   - Writes fixed sampling rate of 54000 Hz into header

    % Load io package for xlsx support in Octave
    if exist('OCTAVE_VERSION', 'builtin')
        pkg load io
    end

    % Defaults
    if nargin == 0
        data_folder_or_pattern = 'DataFile/*.xlsx';
        output_filename = 'output.rawemz';
        fprintf('Using default parameters:\n');
        fprintf('  Input: %s\n', data_folder_or_pattern);
        fprintf('  Output: %s\n', output_filename);
    end

    % Constants
    NUMBER_OF_ANALOG_CH = 64;
    NUMBER_OF_ODOMETERS = 3;
    SAMPLING_RATE = 54000; % Hz (fixed value written to header)

    %% Collect file list
    if ischar(data_folder_or_pattern)
        if exist(data_folder_or_pattern, 'dir')
            fprintf('Scanning folder: %s\n', data_folder_or_pattern);
            if data_folder_or_pattern(end) ~= '/' && data_folder_or_pattern(end) ~= '\'
                data_folder_or_pattern = [data_folder_or_pattern '/'];
            end
            excel_patterns = {'*.xlsx', '*.xls', '*.xlsm'};
            file_list = {};
            for pattern_idx = 1:length(excel_patterns)
                pattern = [data_folder_or_pattern excel_patterns{pattern_idx}];
                try
                    found_files = glob(pattern);
                catch
                    d = dir(pattern);
                    if ~isempty(d)
                        found_files = {d.name};
                        for j = 1:length(found_files)
                            found_files{j} = [data_folder_or_pattern found_files{j}];
                        end
                        found_files = found_files';
                    else
                        found_files = {};
                    end
                end
                file_list = [file_list; found_files];
            end
            if ~isempty(file_list), file_list = sort(file_list); end
        else
            try
                file_list = glob(data_folder_or_pattern);
            catch
                d = dir(data_folder_or_pattern);
                if ~isempty(d)
                    file_list = {d.name}';
                else
                    file_list = {};
                end
            end
        end
    else
        file_list = data_folder_or_pattern; % cell array
    end

    if isempty(file_list)
        error('No Excel files found in the specified location: %s', data_folder_or_pattern);
    end

    fprintf('Found %d Excel files\n', length(file_list));
    for i = 1:min(10, length(file_list))
        fprintf('  %d: %s\n', i, file_list{i});
    end
    if length(file_list) > 10
        fprintf('  ... and %d more files\n', length(file_list) - 10);
    end

    % Respect channel count (96)
    if length(file_list) > NUMBER_OF_ANALOG_CH
        fprintf('Warning: Found %d files, but only using first %d files for %d channels\n', ...
            length(file_list), NUMBER_OF_ANALOG_CH, NUMBER_OF_ANALOG_CH);
        file_list = file_list(1:NUMBER_OF_ANALOG_CH);
    elseif length(file_list) < NUMBER_OF_ANALOG_CH
        fprintf('Warning: Found only %d files, will pad remaining %d channels with zeros\n', ...
            length(file_list), NUMBER_OF_ANALOG_CH - length(file_list));
    end

    fprintf('Processing %d Excel files...\n', length(file_list));

    %% Read files, dedupe time, and keep time>=0
    all_data = struct([]);
    num_files_to_process = length(file_list);
    min_records = inf;

    for i = 1:num_files_to_process
        fprintf('Reading file %d/%d: %s\n', i, num_files_to_process, file_list{i});
        try
            [num_data, ~, ~] = xlsread(file_list{i});
            if size(num_data, 2) < 2
                error('File %s must have at least 2 columns (Time, Signal)', file_list{i});
            end
            time_col   = num_data(:,1);
            signal_col = num_data(:,2);

            % Remove NaNs
            valid_idx = ~isnan(time_col) & ~isnan(signal_col);
            time_col   = time_col(valid_idx);
            signal_col = signal_col(valid_idx);

            if isempty(time_col)
                error('No valid data in file %s', file_list{i});
            end

            % Sort and remove duplicate times (avoid interp1 warnings; though we won't use interp1)
            [time_col, sidx] = sort(time_col(:));
            signal_col = signal_col(sidx);
            [time_unique, keep_idx] = unique(time_col, 'stable');
            signal_col = signal_col(keep_idx);
            time_col   = time_unique;

            % Keep rows from first time >= 0
            start_idx = find(time_col >= 0, 1, 'first');
            if isempty(start_idx)
                error('File %s has no data with time >= 0', file_list{i});
            end

            t0 = time_col(start_idx:end);
            y0 = signal_col(start_idx:end);

            all_data(i).filename     = file_list{i};
            all_data(i).time_from0   = t0;
            all_data(i).signal_from0 = y0;

            min_records = min(min_records, numel(t0));
            fprintf('  -> %d records kept (time >= 0). Starts at %.6f s\n', numel(t0), t0(1));

        catch err
            error('Error reading file %s: %s', file_list{i}, err.message);
        end
    end

    if ~isfinite(min_records) || min_records < 1
        error('No records found with time >= 0 across input files.');
    end

    fprintf('Minimum records (time >= 0) across all files: %d\n', min_records);

    %% Build signal matrix by truncation (no interpolation)
    num_samples = min_records;
    trunc_signals = zeros(num_samples, NUMBER_OF_ANALOG_CH);

    % Real channels
    for i = 1:num_files_to_process
        yi = all_data(i).signal_from0;
        trunc_signals(:, i) = yi(1:num_samples);
    end

    % Dummy channels if fewer than 96 files
    if num_files_to_process < NUMBER_OF_ANALOG_CH
        for i = (num_files_to_process+1):NUMBER_OF_ANALOG_CH
            trunc_signals(:, i) = 0; % zeros length = num_samples
        end
    end

    %% Scale to int16 range
    max_signal = max(abs(trunc_signals(:)));
    if max_signal > 0
        scale_factor = 30000 / max_signal; % headroom
        trunc_signals = trunc_signals * scale_factor;
        fprintf('Signals scaled by factor: %.6f\n', scale_factor);
    else
        fprintf('All-zero signals; no scaling applied.\n');
    end
    trunc_signals = int16(round(trunc_signals));

    %% Generate dummy odometer data using sample index / sampling rate
    fprintf('Generating dummy odometer data...\n');
    sync_time = ((0:num_samples-1).') / SAMPLING_RATE; % seconds, consistent with header SR
    odo_counts = zeros(num_samples, NUMBER_OF_ODOMETERS, 'uint32');
    odo_phases = zeros(num_samples, NUMBER_OF_ODOMETERS, 'uint16');

    for k = 1:NUMBER_OF_ODOMETERS
        speed_factor = 0.1 + (k-1) * 0.05; % different speeds
        odo_counts(:, k) = uint32(floor(sync_time * speed_factor * 100));      % counts
        odo_phases(:, k) = uint16(mod(sync_time * speed_factor * 1000, 4096)); % 12-bit phase
    end

    %% Write .rawemz
    write_rawemz_file(output_filename, trunc_signals, odo_counts, odo_phases, SAMPLING_RATE);

    fprintf('Successfully generated %s with %d records\n', output_filename, num_samples);

    %% Summary
    duration = (num_samples - 1) / SAMPLING_RATE;
    fprintf('\n=== SUMMARY ===\n');
    fprintf('Excel files found: %d\n', num_files_to_process);
    fprintf('Channels used: %d (padded with %d dummy channels)\n', NUMBER_OF_ANALOG_CH, NUMBER_OF_ANALOG_CH - num_files_to_process);
    fprintf('Sampling rate (header): %d Hz\n', SAMPLING_RATE);
    fprintf('Total samples: %d\n', num_samples);
    fprintf('Duration: %.6f seconds (based on SR)\n', duration);
    fprintf('Output file: %s\n', output_filename);
    fprintf('===============\n');
end

function write_rawemz_file(filename, signals, odo_counts, odo_phases, sampling_rate)
% WRITE_RAWEMZ_FILE - Write data to .rawemz binary file format

    NUMBER_OF_ANALOG_CH = 64;
    NUMBER_OF_ODOMETERS = 3;

    [num_samples, num_channels] = size(signals);
    if num_channels ~= NUMBER_OF_ANALOG_CH
        error('Expected %d signal channels, got %d', NUMBER_OF_ANALOG_CH, num_channels);
    end

    % Open file (native endianness). If a specific endianness is required, use 'ieee-le' or 'ieee-be'.
    fid = fopen(filename, 'wb');
    if fid == -1
        error('Cannot create file %s', filename);
    end

    try
        % --- HEADERINFO (example header) ---
        fwrite(fid, 1, 'uint8');          % FirmwareVersion
        fwrite(fid, 0, 'uint8');          % FirmwareSubVersion
        fwrite(fid, 2, 'uint8');          % FirmwareDate
        fwrite(fid, 10, 'uint8');         % FirmwareMonth
        fwrite(fid, 2025, 'uint16');      % FirmwareYear
        fwrite(fid, 'RAW', 'char');       % OperationCode[3]
        fwrite(fid, sampling_rate, 'uint16'); % SamplingRateInHz
        fwrite(fid, 4, 'uint8');          % NumberOfVariableGroups
        fwrite(fid, 914.4, 'float32');    % OdometerDiameterInMM
        fwrite(fid, 800.0, 'float32');    % BodyDiameterInMM
        fwrite(fid, NUMBER_OF_ODOMETERS, 'uint8'); % NumberOfOdometers
        fwrite(fid, 0, 'int32');          % SetupTime_SecondMM
        fwrite(fid, 0, 'uint8');          % EnableAlignData (bool)
        header_size = 512 * 20;           % 10240 bytes
        fwrite(fid, header_size, 'uint32'); % HeaderSizeInByte
        fwrite(fid, 1, 'uint8');          % FirmwareRevision
        fwrite(fid, 0, 'uint8');          % nOdometerType
        fwrite(fid, 750.0, 'float32');    % InternalPipeDiameterInMM

        % Pad to header_size
        current_pos = ftell(fid);
        padding_needed = header_size - current_pos;
        if padding_needed > 0
            fwrite(fid, zeros(padding_needed, 1), 'uint8');
        end

        % --- DATA RECORDS ---
        % Each record: UINT32 counter + INT16 signals[96] + UINT32 odometerCounts[3] + UINT16 odometerPhases[3]
        fprintf('Writing %d data records...\n', num_samples);
        for i = 1:num_samples
            fwrite(fid, i-1, 'uint32');         % sample counter (0-based)
            fwrite(fid, signals(i, :), 'int16');% analog signals
            fwrite(fid, odo_counts(i, :), 'uint32');
            fwrite(fid, odo_phases(i, :), 'uint16');

            if mod(i, 10000) == 0
                fprintf('Written %d/%d records (%.1f%%)\n', i, num_samples, i/num_samples*100);
            end
        end
        fprintf('File writing completed successfully.\n');

    catch err
        fclose(fid);
        rethrow(err);
    end

    fclose(fid);
end

