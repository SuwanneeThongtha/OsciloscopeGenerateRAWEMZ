clear all; clc;

function GenerateRAWEMZ(excel_files_pattern, output_filename)
% GENERATE_RAWEMZ - Generate .rawemz file from multiple Excel files
%
% Usage: generate_rawemz('*.xlsx', 'output.rawemz')
%        generate_rawemz({'file1.xlsx', 'file2.xlsx', ...}, 'output.rawemz')
%
% Parameters:
%   excel_files_pattern - Pattern to match Excel files or cell array of filenames
%   output_filename     - Output .rawemz filename
%
% Requirements:
%   - Each Excel file should have columns: Time, Signal1, Signal2, ..., Signal64
%   - Time should be in seconds
%   - Signal values should be numeric

    % Load io package for xlsx support in Octave
    if exist('OCTAVE_VERSION', 'builtin')
        pkg load io
    end

    % Set default values if called without arguments
    if nargin == 0
        excel_files_pattern = 'DataFile/*.xlsx';
        output_filename = 'output.RAWEMZ';
        fprintf('Using default parameters:\n');
        fprintf('  Input: %s\n', excel_files_pattern);
        fprintf('  Output: %s\n', output_filename);
    end

    % Constants from the C++ code
    NUMBER_OF_ANALOG_CH = 64;
    NUMBER_OF_ODOMETERS = 3;

    % Get list of Excel files
    if ischar(excel_files_pattern)
        files_info = dir(excel_files_pattern);
        file_list = cell(length(files_info), 1);
        for idx = 1:length(files_info)
            file_list{idx} = fullfile(files_info(idx).folder, files_info(idx).name);
        end
    else
        file_list = excel_files_pattern;
    end

    if length(file_list) ~= NUMBER_OF_ANALOG_CH
        warning('Expected %d files, got %d. Will pad with zeros or truncate as needed.', ...
                NUMBER_OF_ANALOG_CH, length(file_list));
    end

    fprintf('Processing %d Excel files...\n', length(file_list));

    % Read all Excel files and find common time range
    all_data = {};
    min_start_time = inf;
    max_end_time = -inf;
    min_records = inf;

    for i = 1:min(length(file_list), NUMBER_OF_ANALOG_CH)
        fprintf('Reading file %d/%d: %s\n', i, length(file_list), file_list{i});

        try
            % Read Excel file (assuming first sheet)
            [num_data, txt_data, raw_data] = xlsread(file_list{i});

            % Assume first column is time, rest are signals
            if size(num_data, 2) < 2
                error('File %s must have at least 2 columns (Time, Signal)', file_list{i});
            end

            time_col = num_data(:, 1);
            signal_col = num_data(:, 2); % Use first signal column

            % Remove NaN values
            valid_idx = ~isnan(time_col) & ~isnan(signal_col);
            time_col = time_col(valid_idx);
            signal_col = signal_col(valid_idx);

            if isempty(time_col)
                error('No valid data in file %s', file_list{i});
            end

            % Store data
            all_data{i}.time = time_col;
            all_data{i}.signal = signal_col;
            all_data{i}.filename = file_list{i};

            % Update time range
            min_start_time = min(min_start_time, min(time_col));
            max_end_time = max(max_end_time, max(time_col));
            min_records = min(min_records, length(time_col));

        catch err
            fprintf('Error reading file %s: %s\n', file_list{i}, err.message);
            % Create dummy data for this channel
            all_data{i}.time = [0; 1];
            all_data{i}.signal = [0; 0];
            all_data{i}.filename = file_list{i};
        end
    end

    % Fill remaining channels with dummy data if needed
    for i = (length(file_list) + 1):NUMBER_OF_ANALOG_CH
        all_data{i}.time = [0; 1];
        all_data{i}.signal = [0; 0];
        all_data{i}.filename = sprintf('dummy_channel_%d', i);
    end

    fprintf('Time range: %.3f to %.3f seconds\n', min_start_time, max_end_time);
    fprintf('Minimum records in any file: %d\n', min_records);

    % Determine sampling rate and create common time vector
    % Use the file with most data points to estimate sampling rate
    max_points = 0;
    best_file_idx = 1;
    for i = 1:NUMBER_OF_ANALOG_CH
        if length(all_data{i}.time) > max_points
            max_points = length(all_data{i}.time);
            best_file_idx = i;
        end
    end

##    % Estimate sampling rate from the best file
##    time_diffs = diff(all_data{best_file_idx}.time);
##    median_dt = median(time_diffs);
##    sampling_rate = round(1 / median_dt);

    % Use fixed sampling rate
    sampling_rate = 1;

##    % Ensure reasonable sampling rate
##    if sampling_rate < 1000 || sampling_rate > 100000
##        sampling_rate = 54000; % Default from the C++ code
##        fprintf('Warning: Estimated sampling rate %d Hz seems unreasonable, using default 54000 Hz\n', sampling_rate);
##    else
##        fprintf('Estimated sampling rate: %d Hz\n', sampling_rate);
##    end

    % Create synchronized time vector starting from 0
    duration = max_end_time - min_start_time;
    num_samples = min(min_records, floor(duration * sampling_rate));

    sync_time = linspace(0, duration, num_samples)';
    actual_time = sync_time + min_start_time;

    fprintf('Synchronized data: %d samples over %.3f seconds\n', num_samples, duration);

    % Interpolate all signals to common time vector
    synchronized_signals = zeros(num_samples, NUMBER_OF_ANALOG_CH);

    for i = 1:NUMBER_OF_ANALOG_CH
        if length(all_data{i}.time) > 1
            % Interpolate signal to common time vector
            synchronized_signals(:, i) = interp1(all_data{i}.time, all_data{i}.signal, ...
                                                actual_time, 'linear', 0);
        else
            % Fill with zeros for dummy channels
            synchronized_signals(:, i) = zeros(num_samples, 1);
        end
    end

    % Convert signals to INT16 range (-32768 to 32767)
    % Normalize signals to reasonable range
    max_signal = max(abs(synchronized_signals(:)));
    if max_signal > 0
        scale_factor = 30000 / max_signal; % Leave some headroom
        synchronized_signals = synchronized_signals * scale_factor;
    end

    % Convert to INT16
    synchronized_signals = int16(round(synchronized_signals));

    % Generate dummy odometer data
    odo_counts = zeros(num_samples, NUMBER_OF_ODOMETERS, 'uint32');
    odo_phases = zeros(num_samples, NUMBER_OF_ODOMETERS, 'uint16');

    % Create realistic odometer progression
    for i = 1:NUMBER_OF_ODOMETERS
        % Simulate different odometer speeds
        speed_factor = 0.1 + (i-1) * 0.05; % Different speeds for each odometer
        odo_counts(:, i) = uint32(floor(sync_time * speed_factor * 10)); % Counts
        odo_phases(:, i) = uint16(mod(sync_time * speed_factor * 1000, 4096)); % Phases (12-bit)
    end

    % Save to Excel file before writing .rawemz
    excel_output_filename = strrep(output_filename, '.RAWEMZ', '.xlsx');
    if strcmp(excel_output_filename, output_filename)
        excel_output_filename = [output_filename '.xlsx'];
    end
    save_to_excel(excel_output_filename, sync_time, synchronized_signals, odo_counts, odo_phases, NUMBER_OF_ANALOG_CH);

    % Write .rawemz file
    write_rawemz_file(output_filename, synchronized_signals, odo_counts, odo_phases, sampling_rate);

    fprintf('Successfully generated %s with %d records\n', output_filename, num_samples);
    fprintf('Successfully generated %s (Excel format)\n', excel_output_filename);
end

function write_rawemz_file(filename, signals, odo_counts, odo_phases, sampling_rate)
% WRITE_RAWEMZ_FILE - Write data to .rawemz binary file format

    NUMBER_OF_ANALOG_CH = 64;
    NUMBER_OF_ODOMETERS = 3;

    [num_samples, num_channels] = size(signals);

    if num_channels ~= NUMBER_OF_ANALOG_CH
        error('Expected %d signal channels, got %d', NUMBER_OF_ANALOG_CH, num_channels);
    end

    % Open file for binary writing
    fid = fopen(filename, 'wb');
    if fid == -1
        error('Cannot create file %s', filename);
    end

    try
        % Write HEADERINFO structure (based on C++ structure)
        % typedef struct tagHEADERINFO

        % BYTE FirmwareVersion
        fwrite(fid, 1, 'uint8');

        % BYTE FirmwareSubVersion
        fwrite(fid, 0, 'uint8');

        % BYTE FirmwareDate
        fwrite(fid, 1, 'uint8');

        % BYTE FirmwareMonth
        fwrite(fid, 1, 'uint8');

        % unsigned short FirmwareYear
        fwrite(fid, 2024, 'uint16');

        % BYTE OperationCode[3] - "RAW"
        fwrite(fid, 'RAW', 'char');

        % unsigned short SamplingRateInHz
        fwrite(fid, sampling_rate, 'uint16');

        % BYTE NumberOfVariableGroups
        fwrite(fid, 4, 'uint8'); % Sample counter, analog signals, odo counts, odo phases

        % float OdometerDiameterInMM
        fwrite(fid, 914.4, 'float32'); % Default 36 inch diameter

        % float BodyDiameterInMM
        fwrite(fid, 800.0, 'float32'); % Default body diameter

        % BYTE NumberOfOdometers
        fwrite(fid, NUMBER_OF_ODOMETERS, 'uint8');

        % int SetupTime_SecondMM
        fwrite(fid, 0, 'int32');

        % bool EnableAlignData
        fwrite(fid, 0, 'uint8');

        % unsigned int HeaderSizeInByte
        header_size = 512 * 20; % Default header size
        fwrite(fid, header_size, 'uint32');

        % BYTE FirmwareRevision
        fwrite(fid, 1, 'uint8');

        % BYTE nOdometerType
        fwrite(fid, 0, 'uint8');

        % float InternalPipeDiameterInMM
        fwrite(fid, 750.0, 'float32');

        % Pad header to required size (10240 bytes = 512 * 20)
        current_pos = ftell(fid);
        padding_needed = header_size - current_pos;
        if padding_needed > 0
            fwrite(fid, zeros(padding_needed, 1), 'uint8');
        end

        % Write data records
        % Each record: UINT32 counter + INT16 signals[NUMBER_OF_ANALOG_CH] + UINT32 odometerCounts[3] + UINT16 odometerPhases[3]

        fprintf('Writing %d data records...\n', num_samples);

        for i = 1:num_samples
            % UINT32 counter (sample counter)
            fwrite(fid, i-1, 'uint32'); % 0-based counter

            % INT16 signals[NUMBER_OF_ANALOG_CH]
            fwrite(fid, signals(i, :), 'int16');

            % UINT32 odometerCounts[NUMBER_OF_ODOMETERS]
            fwrite(fid, odo_counts(i, :), 'uint32');

            % UINT16 odometerPhases[NUMBER_OF_ODOMETERS]
            fwrite(fid, odo_phases(i, :), 'uint16');

            % Progress indicator
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

function save_to_excel(filename, time_vector, signals, odo_counts, odo_phases, num_channels)
% SAVE_TO_EXCEL - Save synchronized data to Excel file
%
% Parameters:
%   filename     - Output Excel filename
%   time_vector  - Time vector (seconds)
%   signals      - Matrix of signals [num_samples x num_channels]
%   odo_counts   - Odometer counts [num_samples x 3]
%   odo_phases   - Odometer phases [num_samples x 3]
%   num_channels - Number of signal channels

    fprintf('Saving data to Excel file: %s\n', filename);

    [num_samples, ~] = size(signals);

    % Create column headers
    headers = {'Time'};
    for i = 1:num_channels
        headers{end+1} = sprintf('Signal%d', i);
    end
    for i = 1:size(odo_counts, 2)
        headers{end+1} = sprintf('Odo%d_Count', i);
    end
    for i = 1:size(odo_phases, 2)
        headers{end+1} = sprintf('Odo%d_Phase', i);
    end

    % Prepare data matrix
    % Convert signals back to double for Excel
    data_matrix = [time_vector, double(signals), double(odo_counts), double(odo_phases)];

    % Write to Excel - Compatible with both MATLAB and Octave
    try
        % Method 1: Try using xlswrite (works in both MATLAB and Octave with io package)
        fprintf('Writing Excel file...\n');

        % Write headers
        xlswrite(filename, headers, 'Sheet1', 'A1');

        % Write data
        xlswrite(filename, data_matrix, 'Sheet1', 'A2');

        fprintf('Excel file saved successfully: %d rows x %d columns\n', num_samples, length(headers));

    catch err
        fprintf('Warning: xlswrite failed: %s\n', err.message);
        fprintf('Attempting CSV export as alternative...\n');

        % Fallback: Save as CSV file
        try
            csv_filename = strrep(filename, '.xlsx', '.csv');

            % Open file for writing
            fid_csv = fopen(csv_filename, 'w');

            % Write headers
            fprintf(fid_csv, '%s', headers{1});
            for i = 2:length(headers)
                fprintf(fid_csv, ',%s', headers{i});
            end
            fprintf(fid_csv, '\n');

            % Write data rows
            for row = 1:num_samples
                fprintf(fid_csv, '%.6f', data_matrix(row, 1));
                for col = 2:size(data_matrix, 2)
                    fprintf(fid_csv, ',%.6f', data_matrix(row, col));
                end
                fprintf(fid_csv, '\n');

                % Progress indicator for large files
                if mod(row, 10000) == 0
                    fprintf('Written %d/%d rows (%.1f%%)\n', row, num_samples, row/num_samples*100);
                end
            end

            fclose(fid_csv);
            fprintf('CSV file saved successfully: %s\n', csv_filename);
            fprintf('You can open this CSV file in Excel.\n');

        catch err2
            fprintf('Error: Could not save file: %s\n', err2.message);
        end
    end
end

% Example usage function
##function example_usage()
##    % Example of how to use the function
##
##    % Method 1: Using file pattern
##    % generate_rawemz('data_*.xlsx', 'output.rawemz');
##
##    % Method 2: Using explicit file list
##    % files = {'channel1.xlsx', 'channel2.xlsx', 'channel3.xlsx'};
##    % generate_rawemz(files, 'output.rawemz');
##
##    fprintf('Example usage:\n');
##    fprintf('  generate_rawemz(''*.xlsx'', ''output.rawemz'');\n');
##    fprintf('  generate_rawemz({''file1.xlsx'', ''file2.xlsx''}, ''output.rawemz'');\n');
##end
