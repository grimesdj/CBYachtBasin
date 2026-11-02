%% drifter_post_processing_sheet_standalone.m
% Standalone version of drifter_post_processing_sheet (no cfg needed).
% Rebuilds drifter tracks from the per-drifter sheet files (data_XX.xlsx or
% data_XX_v2.xlsx), trims them to deployment times, plots them, and
% optionally saves the result.
%
% Requires on the MATLAB path (both are in CBYachtBasin\code_v2\FUNCTIONS):
%   convertToDecimalDegrees.m
%   trimDrifterDataByDeployment_old.m

clear; close all;

%% ======================= USER SETTINGS ===================================

% Deployment times file (same format the logged data use)
dep_fin = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\BHI\DRIFTERS\DEP_TIMES\BHI_DEP_V1.xlsx';

% Folder holding the data_XX.xlsx / data_XX_v2.xlsx files
sheet_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\CB_YachtBasin\Archived_Sheet_Data\RAW\';

% Which version of the sheet files to use (ANYTHING BEFORE 06/2026 is v1 (BHI), AFTER (FALLFLOOD) is v2):
%   false -> data_01.xlsx, data_02.xlsx, ...
%   true  -> data_01_v2.xlsx, data_02_v2.xlsx, ...
use_v2 = false;

% Label used in plot titles (e.g. 'IOP4')
plot_label = 'IOP4';

% Loose pre-filter: before the deployment trimming, keep only each drifter's
% data from (its earliest drop - buffer) to (its latest pull + buffer).
% Set to inf to turn the pre-filter off.
prefilter_buffer = hours(2);

% Save the trimmed tracks? Leave save_file empty ('') to skip saving.
save_file = '';   % e.g. 'C:\path\to\output\IOP4_drifters_level_1.mat'

%% =========================================================================
 
if use_v2
    sheet_file_suffix = '_v2';
else
    sheet_file_suffix = '';
end
 
proj = projcrs(32119);
 
%% Load the drifter release times
dep_times = readtable(dep_fin);
 
deployments.drifter_id = table2array(dep_times(:,1));
% Drop times
deployments.Starttime = table2array(dep_times(:,3));
deployments.Starttime = datetime(deployments.Starttime, 'ConvertFrom', 'excel');
deployments.Starttime = timeofday(deployments.Starttime);
deployments.Startdate = table2array(dep_times(:,4));
deployments.Startdattime = deployments.Startdate + deployments.Starttime;
% Pull times
deployments.Pulltime = table2array(dep_times(:,5));
deployments.Pulltime = datetime(deployments.Pulltime, 'ConvertFrom', 'excel');
deployments.Pulltime = timeofday(deployments.Pulltime);
deployments.Pulldate = table2array(dep_times(:,6));
deployments.Pulldattime = deployments.Pulldate + deployments.Pulltime;
 
%% Read each drifter's sheet file
maxID = max(deployments.drifter_id);   % highest drifter number in the deployment log
Sheet_values = struct();
 
for k = 1:maxID
    fname = sprintf('Drifter_%02d', k);
    file  = fullfile(sheet_dir, sprintf('data_%02d%s.xlsx', k, sheet_file_suffix));
 
    if ~ismember(k, deployments.drifter_id)   % not in this deployment: placeholder keeps drifters(k) = drifter k
        Sheet_values.(fname) = nan(1,4);
        continue
    end
    if ~isfile(file)
        warning('%s not found', file);
        Sheet_values.(fname) = nan(1,4);
        continue
    end
 
    fprintf('Reading %s\n', file);
 
    opts = detectImportOptions(file);
    % Skip empty or incomplete files (need at least 6 columns: date, time, ..., lat, ..., lon)
    if numel(opts.VariableNames) < 6
        warning('%s is empty or has only %d column(s), skipping.', file, numel(opts.VariableNames));
        Sheet_values.(fname) = nan(1,4);
        continue
    end
 
    time_col = opts.VariableNames{2};        % 2nd column, whatever it's named (Var2 or a header)
    opts = setvaropts(opts, time_col, 'Type', 'string');
    T = readtable(file, opts);
 
    if height(T) == 0
        warning('%s has no data rows, skipping.', file);
        Sheet_values.(fname) = nan(1,4);
        continue
    end
 
    time_str = T{:,2};
    t2 = str2double(time_str);               % time of day as fraction of a day (Excel), NaN if text
    if any(~isnan(t2))
        secs = round(t2*86400);              % time of day in seconds
    else
        secs = seconds(duration(time_str, 'InputFormat', 'hh:mm:ss.SSS'));
    end
    dnum = floor(datenum(table2array(T(:,1))));   % date of each post as a datenum
    lat  = table2array(T(:,4));                   % latitude as posted (ddmm.mmmm)
    lon  = table2array(T(:,6));                   % longitude as posted
 
    % --- Loose time pre-filter around this drifter's deployments ---
    this_dep = deployments.drifter_id == k;
    win_start = min(deployments.Startdattime(this_dep)) - prefilter_buffer;
    win_end   = max(deployments.Pulldattime(this_dep))  + prefilter_buffer;
 
    t_row = datetime(dnum, 'ConvertFrom', 'datenum') + seconds(secs);
    keep  = t_row >= win_start & t_row <= win_end;   % NaT rows are dropped too
 
    fprintf('  Drifter %02d: kept %d of %d rows (%s to %s)\n', k, nnz(keep), numel(keep), ...
        datestr(win_start, 'mm/dd HH:MM'), datestr(win_end, 'mm/dd HH:MM'));
 
    dnum = dnum(keep);
    secs = secs(keep);
    lat  = lat(keep);
    lon  = lon(keep);
 
    bad = lat == 0 | lon == 0 | isnan(lat) | isnan(lon) | isnan(dnum) | isnan(secs) | dnum > datenum(today);
    [lat, lon] = convertToDecimalDegrees(lat, lon);   % decimal degrees
    lat(bad) = 0;  lon(bad) = 0;              % zeroed points become NaN below
 
    Sheet_values.(fname) = [dnum(:) secs(:) lat(:) lon(:)];
end

 
%% Build the drifter structure
fn = fieldnames(Sheet_values);
 
for i = 1:numel(fn)
 
    drifterData = Sheet_values.(fn{i});
 
    drifterID = split(fn{i}, '_');
    drifterID = str2double(drifterID{end});
 
    dt = datetime(drifterData(:,1), 'ConvertFrom', 'datenum') + seconds(drifterData(:,2));
    dt.Format = 'dd-MMM-yyyy HH:mm:ss.SSSS';
    e_time = drifterData(:,2);
    Lat = drifterData(:,3);
    Lon = drifterData(:,4);
 
    zeroIndices = find(Lon == 0 | Lat == 0);
    Lon(zeroIndices) = nan;
    Lat(zeroIndices) = nan;
    dt(zeroIndices) = NaT;
 
    % Convert to projected coordinates
    [Easting, Northing] = projfwd(proj, Lat, Lon);
 
    % Velocities
    delta_east  = diff(Easting);
    delta_north = diff(Northing);
    t_abs       = drifterData(:,1)*86400 + drifterData(:,2);   % absolute time (s), safe across midnight
    delta_time  = diff(t_abs);
 
    V_e = delta_east  ./ delta_time;
    V_n = delta_north ./ delta_time;
    V   = sqrt(V_e.^2 + V_n.^2);
 
    drifters(i).ID            = drifterID;
    drifters(i).Date_and_time = datenum(dt);
    drifters(i).etime         = e_time;
    drifters(i).Lat           = Lat;
    drifters(i).Lon           = Lon;
    drifters(i).East          = Easting;
    drifters(i).North         = Northing;
    drifters(i).v_e           = V_e;
    drifters(i).v_n           = V_n;
    drifters(i).V             = V;
    drifters(i).Speed         = nan(size(Lat));   % not in sheet data; trimmer expects the field
    drifters(i).Alt           = nan(size(Lat));
    drifters(i).Angle         = nan(size(Lat));
    drifters(i).Volt          = nan(size(Lat));
end

raw_drifters = drifters; %for VeloQAQC (in code_v2\FUNCTIONS)
 
%% Trim the drifters to when they were in the water
trimmed = trimDrifterDataByDeployment_old(drifters, deployments);
 
%% Plot the trimmed drifters
colors = hsv(numel(trimmed(1,:)));
 
for d = 1:numel(trimmed(:,1))
 
    fig = figure;
    gx = geoaxes;
    hold(gx, 'on');
    geobasemap(gx, 'satellite')
 
    title(gx, sprintf('%s D%d', plot_label, d))
 
    trackHandles = [];
    legendEntries = {};
 
    for i = 1:numel(trimmed(d,:))
 
        lat_tr = trimmed(d,i).Lat;
        lon_tr = trimmed(d,i).Lon;
 
        if isempty(lat_tr)
            continue
        end
 
        t_start = min(trimmed(d,i).Date_and_time);
        t_end   = max(trimmed(d,i).Date_and_time);
 
        h = geoplot(gx, lat_tr, lon_tr, 'Color', colors(i,:), 'LineWidth', 1.5);
 
        trackHandles(end+1) = h; %#ok<SAGROW>
        legendEntries{end+1} = sprintf('D%d T%d | %s-%s', d, i, ...
            datestr(t_start, 'mm/dd HH:MM'), datestr(t_end, 'HH:MM')); %#ok<SAGROW>
 
        % Final point
        geoplot(gx, lat_tr(end), lon_tr(end), 'p', 'MarkerSize', 12, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
    end
 
    if ~isempty(trackHandles)
        legend(gx, trackHandles, legendEntries, 'Location', 'bestoutside');
    end
end
 
%% Save
drifters_level_1 = trimmed;

% if ~isempty(save_file)
%     save(save_file, 'drifters_level_1', 'deployments', 'use_v2', 'dep_fin');
%     fprintf('Saved %s\n', save_file);
% end
