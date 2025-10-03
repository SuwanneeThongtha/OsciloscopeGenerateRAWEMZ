clear all; clc;

% Load data from Excel
filename = 'HeatX52.xlsx';


% Octave-compatible: Use xlsread instead of readmatrix
try
    data = xlsread(filename);
catch
    % If xlsread fails, try loading with io package
    pkg load io;
    data = xlsread(filename);
end

% Extract and flip distance and hardness values
distance_cm_original = data(1,:);  % First row is distance in cm
hardness_values_raw_original = data(2:end,:); % Remaining rows are hardness values

% Flip data in X direction (columns)Tek
distance_cm = (distance_cm_original);
hardness_values_raw = (hardness_values_raw_original);

% Apply moving average (period 3) along X direction (columns)
window_size = 3;
kernel = ones(1, window_size) / window_size;
hardness_values = conv2(hardness_values_raw, kernel, 'same');

% Define y axis labels
y_labels = 0:5:475; % Distance in mm

% Interpolate data for smoothing only in Y direction
[Xq, Yq] = meshgrid(distance_cm, linspace(min(y_labels), max(y_labels), 384)); %Last number must be 4 times of maximum x+5
hardness_smooth = interp2(distance_cm, y_labels', hardness_values, Xq, Yq, 'cubic');

% Plot heatmap
figure;
imagesc(distance_cm, Yq(:,1), hardness_smooth);
colorbar;
colormap(jet);
caxis([200, 340]); % Adjust color range
xlabel('Distance (mm)');
ylabel('Distance (mm)');
title('Hardness mapping (HV)');
set(gca, 'YDir', 'normal'); % Ensure correct Y-axis direction


% --- Add Logo as Overlay (Top-Left Corner) ---
try
    ax_logo = axes('Position', [0.15 0.92 0.08 0.08]); % [left bottom width height]
    imshow('New folder/Picture1.png');  % Updated path for logo
    axis off;
    uistack(ax_logo, 'top');  % Ensure it stays on top of plots
catch
    warning('Could not load logo image');
end
