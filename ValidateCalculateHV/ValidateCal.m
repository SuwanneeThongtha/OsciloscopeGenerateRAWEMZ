% ฟังก์ชันสำหรับคำนวณ FFT จากไฟล์ Excel พร้อมคำนวณ V1 (3kHz) และ V3 (9kHz)
% ใช้วิธีการเดียวกับโค้ด C++ (Blackman window + compensation factor)
% คำนวณ FFT แบบ sliding window สำหรับแต่ละจุดข้อมูล
%
% การใช้งาน:
%   ValidateCal(input_file, output_file, window_size)
%
% Input:
%   input_file   - ชื่อไฟล์ Excel ที่มี time และ signal (default: 'Tek717 TestCalculate.xlsx')
%   output_file  - ชื่อไฟล์ Excel สำหรับ export (default: 'FFT_Results_V1_V3.xlsx')
%   window_size  - ขนาดของ FFT window (default: 54, เหมือนใน C++)
%
% Output Excel Columns:
%   - Time (s)
%   - Signal_Original
%   - Magnitude_3kHz
%   - Magnitude_9kHz
%   - V1 (3kHz with compensation)
%   - V3 (9kHz with compensation)
%   - Ratio (C = 10*V1/V3)
%   - HV (calculated from V1 and V3)

function ValidateCal(input_file, output_file, window_size)
    % ตั้งค่า default parameters
    if nargin < 1
        input_file = 'Tek717 TestCalculate.xlsx';
    end
    if nargin < 2
        output_file = 'FFT_Results_V1_V3.xlsx';
    end
    if nargin < 3
        window_size = 54;  % Default window size (เหมือนใน C++)
    end

    % ตรวจสอบ package ที่จำเป็น
    pkg load io

    % ===================================================
    % อ่านข้อมูลจาก Excel
    % ===================================================
    fprintf('Reading data from %s...\n', input_file);
    [num, txt, raw] = xlsread(input_file);

    % สมมติว่าคอลัมน์แรกคือ time และคอลัมน์ที่สองคือ signal
    time = num(:, 1);
    signal_original = num(:, 2);

    % ลบค่า NaN ถ้ามี
    valid_idx = ~isnan(time) & ~isnan(signal_original);
    time = time(valid_idx);
    signal_original = signal_original(valid_idx);

    % จำนวนจุดข้อมูล
    N_total = length(signal_original);

    % คำนวณ sampling frequency
    dt = mean(diff(time));
    Fs = 1 / dt;

    fprintf('Number of data points: %d\n', N_total);
    fprintf('Sampling frequency: %.2f Hz\n', Fs);
    fprintf('Time duration: %.4f s\n', time(end) - time(1));
    fprintf('FFT window size: %d points\n', window_size);

    % ===================================================
    % สร้าง Blackman Window (เหมือนใน C++)
    % ===================================================
    fprintf('Creating Blackman window (size %d)...\n', window_size);
    blackman_window = zeros(window_size, 1);
    for i = 1:window_size
        blackman_window(i) = 0.42 - 0.5 * cos(2.0 * pi * (i-1) / (window_size - 1)) ...
                                  + 0.08 * cos(4.0 * pi * (i-1) / (window_size - 1));
    end

    % ===================================================
    % คำนวณ index สำหรับ 3kHz และ 9kHz
    % ===================================================
    target_freq_3k = 3000;  % Hz
    target_freq_9k = 9000;  % Hz

    % Compensation factor สำหรับ Blackman window (เหมือนใน C++)
    compensation_factor = 2.38;  % 1.0 / 0.42

    % คำนวณ index (ตามสูตรใน C++)
    index_3khz = round((target_freq_3k / Fs) * window_size);
    index_9khz = round((target_freq_9k / Fs) * window_size);

    % ตรวจสอบให้ index อยู่ในช่วงที่ถูกต้อง
    index_3khz = min(max(index_3khz, 1), floor(window_size/2));
    index_9khz = min(max(index_9khz, 1), floor(window_size/2));

    % คำนวณความถี่จริง
    freq_resolution = Fs / window_size;
    actual_freq_3k = index_3khz * freq_resolution;
    actual_freq_9k = index_9khz * freq_resolution;

    fprintf('\nFFT Configuration:\n');
    fprintf('  Frequency resolution: %.2f Hz\n', freq_resolution);
    fprintf('  3kHz: index = %d, actual freq = %.2f Hz\n', index_3khz, actual_freq_3k);
    fprintf('  9kHz: index = %d, actual freq = %.2f Hz\n', index_9khz, actual_freq_9k);
    fprintf('  Compensation factor: %.2f\n', compensation_factor);

    % ===================================================
    % คำนวณ FFT แบบ sliding window สำหรับแต่ละจุดข้อมูล
    % ===================================================
    fprintf('\nCalculating FFT with sliding window...\n');

    % เตรียม arrays สำหรับเก็บผลลัพธ์
    magnitude_3k_array = zeros(N_total, 1);
    magnitude_9k_array = zeros(N_total, 1);
    V1_array = zeros(N_total, 1);
    V3_array = zeros(N_total, 1);
    Ratio_array = zeros(N_total, 1);
    HV_array = zeros(N_total, 1);

    % Progress bar setup
    progress_step = max(1, floor(N_total / 20));

    % คำนวณ FFT สำหรับแต่ละจุด
    for i = 1:N_total
        % แสดง progress
        if mod(i, progress_step) == 0 || i == N_total
            fprintf('  Progress: %d/%d (%.1f%%)\n', i, N_total, (i/N_total)*100);
        end

        % กำหนดช่วงของ window (circular buffer)
        window_indices = mod((i-1):(i-1+window_size-1), N_total) + 1;

        % ดึงข้อมูล signal ในช่วง window
        signal_window = signal_original(window_indices);

        % Apply Blackman window
        signal_windowed = signal_window .* blackman_window;

        % คำนวณ FFT
        Y = fft(signal_windowed);

        % คำนวณ magnitude ที่ 3kHz (ตามสูตรใน C++)
        real_3k = real(Y(index_3khz));
        imag_3k = imag(Y(index_3khz));
        mag_3k = 2.0 / window_size * sqrt(real_3k^2 + imag_3k^2);

        % คำนวณ magnitude ที่ 9kHz
        real_9k = real(Y(index_9khz));
        imag_9k = imag(Y(index_9khz));
        mag_9k = 2.0 / window_size * sqrt(real_9k^2 + imag_9k^2);

        % บันทึกผลลัพธ์
        magnitude_3k_array(i) = mag_3k;
        magnitude_9k_array(i) = mag_9k;
        V1_array(i) = mag_3k * compensation_factor;
        V3_array(i) = mag_9k * compensation_factor;
        
        % คำนวณ Ratio และ HV (ตามสูตรใน C++)
        Ratio_array(i) = CalculateRatio(V1_array(i), V3_array(i));
        HV_array(i) = CalculateHV(V1_array(i), V3_array(i));
    end

    fprintf('✓ FFT calculation complete!\n');

    % แสดงสถิติผลลัพธ์
    fprintf('\n=== Results Statistics ===\n');
    fprintf('V1 (3kHz):\n');
    fprintf('  Mean:   %.6f (x1000 = %.4f)\n', mean(V1_array), mean(V1_array)*1000);
    fprintf('  Min:    %.6f (x1000 = %.4f)\n', min(V1_array), min(V1_array)*1000);
    fprintf('  Max:    %.6f (x1000 = %.4f)\n', max(V1_array), max(V1_array)*1000);
    fprintf('  Std:    %.6f\n', std(V1_array));

    fprintf('\nV3 (9kHz):\n');
    fprintf('  Mean:   %.6f (x1000 = %.4f)\n', mean(V3_array), mean(V3_array)*1000);
    fprintf('  Min:    %.6f (x1000 = %.4f)\n', min(V3_array), min(V3_array)*1000);
    fprintf('  Max:    %.6f (x1000 = %.4f)\n', max(V3_array), max(V3_array)*1000);
    fprintf('  Std:    %.6f\n', std(V3_array));

    fprintf('\nRatio (C = 10*V1/V3):\n');
    fprintf('  Mean:   %.6f\n', mean(Ratio_array));
    fprintf('  Min:    %.6f\n', min(Ratio_array));
    fprintf('  Max:    %.6f\n', max(Ratio_array));
    fprintf('  Std:    %.6f\n', std(Ratio_array));

    fprintf('\nHV Value:\n');
    fprintf('  Mean:   %.6f\n', mean(HV_array));
    fprintf('  Min:    %.6f\n', min(HV_array));
    fprintf('  Max:    %.6f\n', max(HV_array));
    fprintf('  Std:    %.6f\n', std(HV_array));

    % ===================================================
    % Export ไปยัง Excel
    % ===================================================
    fprintf('\nExporting to %s...\n', output_file);

    % Sheet 1: Raw Data with Results
    header1 = {'Time (s)', 'Signal_Original', 'Magnitude_3kHz', 'Magnitude_9kHz', 'V1', 'V3', 'Ratio', 'HV'};
    data1 = [time, signal_original, magnitude_3k_array, magnitude_9k_array, V1_array, V3_array, Ratio_array, HV_array];

    xlswrite(output_file, header1, 'Raw_Data', 'A1');
    xlswrite(output_file, data1, 'Raw_Data', 'A2');

    % Sheet 2: Summary Statistics
    summary_header = {'Parameter', 'Value'};
    summary_data = {
        'Input File', input_file;
        'Number of Data Points', num2str(N_total);
        'Sampling Frequency (Hz)', num2str(Fs);
        'Time Duration (s)', num2str(time(end) - time(1));
        'FFT Window Size', num2str(window_size);
        'Frequency Resolution (Hz)', num2str(freq_resolution);
        'Compensation Factor', num2str(compensation_factor);
        '', '';
        '--- 3kHz (V1) ---', '';
        'Target Frequency (Hz)', num2str(target_freq_3k);
        'Actual Frequency (Hz)', num2str(actual_freq_3k);
        'FFT Index', num2str(index_3khz);
        'V1 Mean', num2str(mean(V1_array), '%.6f');
        'V1 Min', num2str(min(V1_array), '%.6f');
        'V1 Max', num2str(max(V1_array), '%.6f');
        'V1 Std', num2str(std(V1_array), '%.6f');
        'V1 Mean x 1000', num2str(mean(V1_array) * 1000, '%.4f');
        '', '';
        '--- 9kHz (V3) ---', '';
        'Target Frequency (Hz)', num2str(target_freq_9k);
        'Actual Frequency (Hz)', num2str(actual_freq_9k);
        'FFT Index', num2str(index_9khz);
        'V3 Mean', num2str(mean(V3_array), '%.6f');
        'V3 Min', num2str(min(V3_array), '%.6f');
        'V3 Max', num2str(max(V3_array), '%.6f');
        'V3 Std', num2str(std(V3_array), '%.6f');
        'V3 Mean x 1000', num2str(mean(V3_array) * 1000, '%.4f');
        '', '';
        '--- Ratio (C = 10*V1/V3) ---', '';
        'Ratio Mean', num2str(mean(Ratio_array), '%.6f');
        'Ratio Min', num2str(min(Ratio_array), '%.6f');
        'Ratio Max', num2str(max(Ratio_array), '%.6f');
        'Ratio Std', num2str(std(Ratio_array), '%.6f');
        '', '';
        '--- HV Value ---', '';
        'HV Mean', num2str(mean(HV_array), '%.6f');
        'HV Min', num2str(min(HV_array), '%.6f');
        'HV Max', num2str(max(HV_array), '%.6f');
        'HV Std', num2str(std(HV_array), '%.6f')
    };

    xlswrite(output_file, summary_header, 'Summary', 'A1');
    xlswrite(output_file, summary_data, 'Summary', 'A2');

    fprintf('✓ Export successful!\n');
    fprintf('  - Sheet 1: Raw_Data (%d rows)\n', N_total);
    fprintf('  - Sheet 2: Summary (Statistics)\n');

    % ===================================================
    % Plot กราฟ
    % ===================================================
    fprintf('\nCreating plots...\n');

    figure('Name', 'FFT Analysis with V1 and V3 (Sliding Window)', 'Position', [100, 100, 1400, 900]);

    % Plot 1: Signal ต้นฉบับ
    subplot(3, 2, 1);
    plot(time, signal_original, 'b-', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Original Signal');
    grid on;

    % Plot 2: V1 (3kHz) ตามเวลา
    subplot(3, 2, 2);
    plot(time, V1_array, 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('V1 (3kHz)');
    title(sprintf('V1 over Time (Mean = %.6f)', mean(V1_array)));
    grid on;

    % Plot 3: V3 (9kHz) ตามเวลา
    subplot(3, 2, 3);
    plot(time, V3_array, 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('V3 (9kHz)');
    title(sprintf('V3 over Time (Mean = %.6f)', mean(V3_array)));
    grid on;

    % Plot 4: V1 และ V3 เปรียบเทียบ
    subplot(3, 2, 4);
    plot(time, V1_array, 'g-', 'LineWidth', 1.5);
    hold on;
    plot(time, V3_array, 'r-', 'LineWidth', 1.5);
    hold off;
    xlabel('Time (s)');
    ylabel('Magnitude');
    title('V1 and V3 Comparison');
    legend('V1 (3kHz)', 'V3 (9kHz)', 'Location', 'best');
    grid on;

    % Plot 5: Histogram ของ V1
    subplot(3, 2, 5);
    histogram(V1_array, 30, 'FaceColor', 'g');
    xlabel('V1 Value');
    ylabel('Count');
    title(sprintf('V1 Distribution (Mean=%.6f, Std=%.6f)', mean(V1_array), std(V1_array)));
    grid on;

    % Plot 6: Histogram ของ V3
    subplot(3, 2, 6);
    histogram(V3_array, 30, 'FaceColor', 'r');
    xlabel('V3 Value');
    ylabel('Count');
    title(sprintf('V3 Distribution (Mean=%.6f, Std=%.6f)', mean(V3_array), std(V3_array)));
    grid on;

    fprintf('\n✓ Finished!\n');
end

% ===================================================
% ฟังก์ชันช่วย: คำนวณ Ratio (ตามโค้ด C++)
% ===================================================
function ratio = CalculateRatio(v1, v3)
    EPSILON = 1e-10;
    if abs(v3) > EPSILON
        ratio = 10 * v1 / v3;
    else
        ratio = 0.000;
    end
end

% ===================================================
% ฟังก์ชันช่วย: คำนวณ HV (ตามโค้ด C++)
% ===================================================
function hvValue = CalculateHV(v1, v3)
    EPSILON = 1e-10;
    
    % First get the ratio (C in the formula)
    ratio = CalculateRatio(v1, v3);
    
    % Check for valid ratio to prevent division by zero or invalid calculations
    if ratio < EPSILON
        hvValue = 0.000;
        return;
    end
    
    term1 = (1 - (1 / (1 + exp((ratio - 3.4) / 0.4)))) * (ratio / 300.0)^(-4/3);
    term2 = (1 / (1 + exp((ratio - 3.4) / 0.4))) * (ratio / 26.1)^(-20/7);
    
    % Final HV = term1 + term2
    hvValue = term1 + term2;
    
    % Ensure non-negative value
    hvValue = max(hvValue, 0.000);
end
