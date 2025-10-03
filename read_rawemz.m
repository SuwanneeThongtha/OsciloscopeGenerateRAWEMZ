% Octave script to read the first 64 records of all 64 channels from a RAWEMZ file
% and display the header contents in a meaningful way, along with plotting the signals.
% Octave script to read the first 64 records of all 64 channels from a RAWEMZ file
% and display the header contents in a meaningful way, along with plotting the signals.

function read_rawemz(filename)
  % Open the file
  fid = fopen(filename, 'rb');
  if fid == -1
    error('Cannot open the file: %s', filename);
  end

  % Read the HEADERINFO structure
  header = struct();
  header.FirmwareVersion = fread(fid, 1, 'uint8');
  header.FirmwareSubVersion = fread(fid, 1, 'uint8');
  header.FirmwareDate = fread(fid, 1, 'uint8');
  header.FirmwareMonth = fread(fid, 1, 'uint8');
  header.FirmwareYear = fread(fid, 1, 'uint16');
  header.OperationCode = char(fread(fid, 3, 'char')');
  header.SamplingRateInHz = fread(fid, 1, 'uint16');
  header.NumberOfVariableGroups = fread(fid, 1, 'uint8');
  header.OdometerDiameterInMM = fread(fid, 1, 'float32');
  header.BodyDiameterInMM = fread(fid, 1, 'float32');
  header.NumberOfOdometers = fread(fid, 1, 'uint8');
  header.SetupTime_SecondMM = fread(fid, 1, 'int32');
  header.HeaderSizeInByte = fread(fid, 1, 'uint32');
  header.FirmwareRevision = fread(fid, 1, 'uint8');
  header.nOdometerType = fread(fid, 1, 'uint8');
  header.InternalPipeDiameterInMM = fread(fid, 1, 'float32');

  % Display header information
  fprintf('Firmware Version: %d.%d\n', header.FirmwareVersion, header.FirmwareSubVersion);
  fprintf('Firmware Date: %d/%d/%d\n', header.FirmwareDate, header.FirmwareMonth, header.FirmwareYear);
  fprintf('Operation Code: %s\n', header.OperationCode);
  fprintf('Sampling Rate: %d Hz\n', header.SamplingRateInHz);
  fprintf('Number of Variable Groups: %d\n', header.NumberOfVariableGroups);
  fprintf('Odometer Diameter: %.2f mm\n', header.OdometerDiameterInMM);
  fprintf('Body Diameter: %.2f mm\n', header.BodyDiameterInMM);
  fprintf('Number of Odometers: %d\n', header.NumberOfOdometers);
  fprintf('Setup Time (SecondMM): %d\n', header.SetupTime_SecondMM);
  fprintf('Header Size: %d bytes\n', header.HeaderSizeInByte);
  fprintf('Firmware Revision: %d\n', header.FirmwareRevision);
  fprintf('Odometer Type: %d\n', header.nOdometerType);
  fprintf('Internal Pipe Diameter: %.2f mm\n', header.InternalPipeDiameterInMM);

  % Skip the remaining part of the header
  fseek(fid, header.HeaderSizeInByte, 'bof');

  % Define the structure for the EMZRECORD
  num_signals = 64;
  num_records = 54;

  % Preallocate signal matrix for plotting
  count_matrix = zeros(num_records, 1);
  signal_matrix = zeros(num_records, num_signals);
  odometer_count_matrix = zeros(num_records, 3);
  odometer_phase_matrix = zeros(num_records, 3);

  % Read the first 54 records
  for record_idx = 1:num_records
    counter = fread(fid, 1, 'uint32');
    signals = fread(fid, num_signals, 'int16');
    odometer_count = fread(fid, 3, 'uint32');
    odometer_phase = fread(fid, 3, 'uint16');

    % Store signals and odometer data for plotting
    count_matrix(record_idx, 1) = counter;
    signal_matrix(record_idx, :) = signals / 32767;  % convert to volt
    odometer_count_matrix(record_idx, :) = odometer_count;
    odometer_phase_matrix(record_idx, :) = odometer_phase;

    % Display the data for each record
    fprintf('\nRecord %d:\n', record_idx);
    fprintf('Counter: %d\n', counter);
    fprintf('Signals: ');
    fprintf('%d ', signals);
    fprintf('\n');
    fprintf('Odometer Count: %d %d %d\n', odometer_count);
    fprintf('Odometer Phase: %d %d %d\n', odometer_phase);
  end

  figure; plot(0:num_records-1, count_matrix(:,1)); zoom on; grid on;
  % Plot the signals grouped by 8 channels each
  figure;
  for group_idx = 1:8
    subplot(4, 2, group_idx);
    plot(signal_matrix(:, (group_idx - 1) * 8 + 1 : group_idx * 8) );
    title(sprintf('Channels %d-%d', (group_idx - 1) * 8 + 1, group_idx * 8));
    xlabel('Record Index');
    ylabel('Amplitude'); zoom on; grid on;
  end
  % Use title instead of sgtitle for compatibility
  annotation('textbox', [0.5, 0.95, 0, 0], 'string', 'Signal Groups from First 54 Records of RAWEMZ File', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'EdgeColor', 'none');

  % Plot the odometer data
  figure;
  subplot(4, 1, 1);
  plot(odometer_count_matrix);
  title('Odometer Count (First 64 Records)');
  xlabel('Record Index');
  ylabel('Odometer Count');
  legend('Odometer 1', 'Odometer 2', 'Odometer 3');

  % Plot odometer phases individually
  for i = 1:3
    subplot(4, 1, i + 1);
    plot(odometer_phase_matrix(:, i));
    title(sprintf('Odometer Phase %d (First 54 Records)', i));
    xlabel('Record Index');
    ylabel('Odometer Phase');
    legend(sprintf('Odometer %d', i));
  end

  % Perform FFT on the signals and plot amplitude vs frequency grouped by 8 channels each
  figure;
  fft_length = 54;
  freq_axis = (0:(fft_length/2)-1) * (header.SamplingRateInHz / fft_length);
  for group_idx = 1:8
    subplot(4, 2, group_idx);
    for channel_idx = (group_idx - 1) * 8 + 1 : group_idx * 8
      signal = signal_matrix(:, channel_idx);
      fft_result = fft(signal , fft_length);
      fft_amplitude = (2/fft_length)*(abs(fft_result(1:fft_length/2)));
      plot(freq_axis, fft_amplitude);
      hold on;
    end
    title(sprintf('FFT of Channels %d-%d', (group_idx - 1) * 8 + 1, group_idx * 8));
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    hold off;
    zoom on; grid on;
  end
  annotation('textbox', [0.5, 0.95, 0, 0], 'string', 'FFT of Signals Grouped by 8 Channels from First 54 Records of RAWEMZ File', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'EdgeColor', 'none');

  % Close the file
  fclose(fid);
end

% Example usage
% filename = 'd:\example.RAWEMZ';
% read_rawemz(filename);

