%% Hardness heatmap with hotspot size filter (Y>5 mm, X>10 mm)
clear; clc; close all;

% -------- I/O --------
excel_filename = 'HeatX52.xlsx';   % first row: distance (cm); remaining rows: HV
logo_filename  = 'Picture1.png';   % optional

% -------- Parameters --------
window_size   = 3;                 % moving average along X
y_labels_mm   = 0:5:475;
n_rows_interp = 384;               % Y-interpolation for display
caxis_min     = 170;
caxis_max     = 320;

T_fixed      = 260;                % HV threshold
min_y_mm     = 5;                  % minimum vertical span
min_x_mm     = 10;                 % minimum horizontal span

% Morphology parameters for sharper blobs
r_close_mm = 0.8;
r_open_mm  = 0.5;

draw_boxes     = true;
box_linewidth  = 1.5;

%% -------- Load & prepare --------
raw = readmatrix(excel_filename);
distance_cm = raw(1,:); 
HV_raw      = raw(2:end,:);

distance_mm = distance_cm * 10;

% Smooth along X
kernel = ones(1, window_size) / window_size;
HV_smoothX = conv2(HV_raw, kernel, 'same');

% Interpolation
[Xq, Yq]     = meshgrid(distance_mm, linspace(min(y_labels_mm), max(y_labels_mm), n_rows_interp));
HV_map_disp  = interp2(distance_mm, y_labels_mm', HV_smoothX, Xq, Yq, 'cubic');
HV_map_det   = interp2(distance_mm, y_labels_mm', HV_smoothX, Xq, Yq, 'linear');

%% -------- Plot heatmap --------
figure('Color','w');
imagesc(distance_mm, Yq(:,1), HV_map_disp);
set(gca,'YDir','normal');
colormap(jet); caxis([caxis_min, caxis_max]); colorbar;
xlabel('Distance (mm)'); ylabel('Distance (mm)');
title('Hardness mapping (HV)');
axis tight;

%% -------- Hotspot detection --------
T  = T_fixed;
BW = HV_map_det >= T;

dx = mean(diff(distance_mm));
dy = mean(diff(Yq(:,1)));

% Morphological cleanup
BW = imclose(BW, strel('disk', max(1,round(r_close_mm/min(dx,dy)))));
BW = imopen(BW,  strel('disk', max(1,round(r_open_mm/min(dx,dy)))));

% Connected components filtering by X/Y span
CC = bwconncomp(BW);
keep_mask = false(size(BW));
for k = 1:CC.NumObjects
    [rows, cols] = ind2sub(size(BW), CC.PixelIdxList{k});
    x_vals = distance_mm(cols);
    y_vals = Yq(rows,1);
    if (max(y_vals)-min(y_vals) >= min_y_mm) && (max(x_vals)-min(x_vals) >= min_x_mm)
        keep_mask(CC.PixelIdxList{k}) = true;
    end
end
BW = keep_mask;

%% -------- Draw contours and boxes --------
hold on;
[C,h] = contour(distance_mm, Yq(:,1), HV_map_det, [T T], 'LineWidth',2,'LineColor','w');
hold off;

if draw_boxes
    CC = bwconncomp(BW);
    hold on;
    for k = 1:CC.NumObjects
        [rows, cols] = ind2sub(size(BW), CC.PixelIdxList{k});
        x_vals = distance_mm(cols);
        y_vals = Yq(rows,1);
        x_min = min(x_vals); x_max = max(x_vals);
        y_min = min(y_vals); y_max = max(y_vals);

        rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], ...
                  'EdgeColor', 'w', 'LineWidth', box_linewidth, 'LineStyle','--');

        hv_region = HV_map_det(CC.PixelIdxList{k});
        [hv_max, idx_max] = max(hv_region);
        r_max = rows(idx_max); c_max = cols(idx_max);
        text(distance_mm(c_max), Yq(r_max,1), sprintf('%.0f', hv_max), ...
             'Color', 'w', 'FontWeight','bold', ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'BackgroundColor','k','Margin',1);
    end
    hold off;
end

%% -------- Logo overlay --------
try
    if exist(logo_filename,'file')==2
        ax_logo = axes('Position',[0.15 0.92 0.08 0.08]);
        imshow(logo_filename); axis off; uistack(ax_logo,'top');
    end
catch ME
    warning('Logo overlay failed: %s',ME.message);
end

disp('Done: Regions kept only if Y span > 5 mm AND X span > 10 mm, with HV >= 260.');
