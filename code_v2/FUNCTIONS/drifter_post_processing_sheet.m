function [drifters_level_1] = drifter_post_processing_sheet(cfg)
% Rebuilds the drifter tracks from the 5-minute data posted to the Google Sheet.
% Used only when a SHEET observation period is chosen in main_driver.
% Reads cfg.sheet_data_file (sheet_data.mat, variable Sheet_values) and the same deployment
% times file the logged data use (cfg.drifters_dep_times).

% sheet_file = cfg.sheet_data_file;        % full path of sheet_data.mat
dep_fin    = cfg.drifters_dep_times;     % deployment times file for this period

% S = load(sheet_file);                    % the saved sheet data
% Sheet_values = S.Sheet_values;           % one entry per drifter

proj = projcrs(32119);
%% Load the drifter release times
dep_times = readtable(dep_fin);

deployments.drifter_id = table2array(dep_times(:,1));
%Drop times
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

%% Parse the sheet data to be trimmed

%% Read each drifter's sheet file (data_01.xlsx, data_02.xlsx, ...)
sheet_dir = cfg.sheet_data_dir;                 % folder holding the data_XX.xlsx files
maxID = max(deployments.drifter_id);            % highest drifter number in the deployment log
clear Sheet_values
for k = 1:maxID
    fname = sprintf('Drifter_%02d', k);         % field name the loop below turns back into the drifter number
        file  = fullfile(sheet_dir, sprintf('data_%02d%s.xlsx', k, cfg.sheet_file_suffix));
    if ~ismember(k, deployments.drifter_id)     % not in this deployment: empty placeholder keeps drifters(k) = drifter k
        Sheet_values.(fname) = nan(1,4);
        continue
    end
    if ~isfile(file)
        warning('%s not found', file);
        Sheet_values.(fname) = nan(1,4);
        continue
    end

    opts = detectImportOptions(file);
    opts = setvaropts(opts, 'Var2', 'Type', 'string');
    T = readtable(file, opts);

    t2 = str2double(T.Var2);                    % time of day as a fraction of a day (Excel), NaN if text
    if any(~isnan(t2))
        secs = round(t2*86400);                 % time of day in seconds
    else
        secs = seconds(duration(T.Var2, 'InputFormat', 'hh:mm:ss.SSS'));
    end
    dnum = floor(datenum(table2array(T(:,1))));   % date of each post as a datenum
    lat  = table2array(T(:,4));                   % latitude as posted (ddmm.mmmm)
    lon  = table2array(T(:,6));                   % longitude as posted

    bad = lat == 0 | lon == 0 | isnan(lat) | isnan(lon) | isnan(dnum) | isnan(secs) | dnum > datenum(today);
    [lat, lon] = convertToDecimalDegrees(lat, lon);   % decimal degrees
    lat(bad) = 0;  lon(bad) = 0;                % the loop below turns zero lat/lon into NaN, as before

    Sheet_values.(fname) = [dnum(:) secs(:) lat(:) lon(:)];
end

fn = fieldnames(Sheet_values);

for i = 1:numel(fn);
    
    drifterData = Sheet_values.(fn{i});
    
    drifterID = split(fn{i},'_');
    drifterID = str2num(drifterID{end});
    
    dt = datetime(drifterData(:,1),'ConvertFrom','datenum') + seconds(drifterData(:,2));
    dt.Format = 'dd-MMM-yyyy HH:mm:ss.SSSS';
    e_time = drifterData(:,2);
    Lat = drifterData(:,3);
    Lon = drifterData(:,4);
    
    
    zeroIndices = find(Lon == 0 | Lat==0); %| [0; diff(time_mat)]<0.5/86400);
    
    
    Lon(zeroIndices) = nan;
    Lat(zeroIndices) = nan;
    dt(zeroIndices) = NaT;
    
    % Convert the Lat and Lon to decimal degrees
    % [Lat,Lon] = convertToDecimalDegrees(Lat,Lon);

    % Convert to ENU
    [Easting,Northing] = projfwd(proj,Lat,Lon);
    
    % Calculate velocities
    delta_east = diff(Easting);
    delta_north = diff(Northing);
    
    delta_time = diff(e_time);
    
    V_e = delta_east./delta_time;
    V_n = delta_north./delta_time;
    V = sqrt((V_e.^2) + (V_n.^2));
    
    % log the values
    drifters(i).ID = drifterID;
%     drifters(i).Date_and_time = dt;
    drifters(i).Date_and_time = datenum(dt);
    drifters(i).etime = e_time;
    drifters(i).Lat = Lat;
    drifters(i).Lon = Lon;
    drifters(i).East = Easting;
    drifters(i).North = Northing;
    drifters(i).v_e = V_e;
    drifters(i).v_n = V_n;
    drifters(i).V = V;
    drifters(i).Speed = nan(size(Lat));      % not in the sheet data; the trimmer expects the field
    drifters(i).Alt   = nan(size(Lat));
    drifters(i).Angle = nan(size(Lat));
    drifters(i).Volt  = nan(size(Lat));
    

end

%% trim the drifters to when in the water
trimmed = trimDrifterDataByDeployment_old(drifters,deployments);

%% Plot the trimmed drifters recovered from the google sheet

% fixed lander location
lat_N = 34.05633533;
lon_N = -77.88924400;

lat_S = 34.04653800;
lon_S = -77.88883175;

colors = hsv(numel(trimmed(1,:)));

for d = 1:numel(trimmed(:,1))
    % Create a new figure + geoaxes for each deployment
    gcf1 = figure;
    gx = geoaxes;
    hold(gx,'on');
    geobasemap(gx,'satellite')

    % Plot fixed locations (N and S)

%     hN = geoplot(gx, lat_N, lon_N, 's', ...
%     'MarkerSize', 10, ...
%     'MarkerEdgeColor', 'r', ...
%     'MarkerFaceColor', 'r');   % Black square (N)
% 
%     hS = geoplot(gx, lat_S, lon_S, 's', ...
%     'MarkerSize', 10, ...
%     'MarkerEdgeColor', 'y', ...
%     'MarkerFaceColor', 'y');   % Red square (S)
    
    % Title for this deployment
    title(gx, sprintf('IOP4 D%d', d))
    
    % Reset legend storage for this deployment
    trackHandles = [];
    legendEntries = {};
    
    % Loop through all trajectories for this deployment
    for i = 1:numel(trimmed(d,:))
        
        lat_tr = trimmed(d,i).Lat;
        lon_tr = trimmed(d,i).Lon;
        t_start = min(trimmed(d,i).Date_and_time);
        t_end = max(trimmed(d,i).Date_and_time);
        
        dt = t_end - t_start;
        dt_hr = hours(dt);
        
        if isempty(lat_tr)
            continue
        end
        
        % Plot track and store handle
        h = geoplot(gx, lat_tr, lon_tr, 'color', colors(i,:), 'LineWidth', 1.5);
        
        % Store for legend
        trackHandles(end+1) = h; 
        % legendEntries{end+1} = sprintf('d%d_t%d (%.1f h)', d, i, dt_hr); 
        legendEntries{end+1} = sprintf('D%d T%d | %s–%s', d, i, datestr(t_start,'mm/dd HH:MM'), datestr(t_end,'HH:MM'));
        
        % Plot final point as yellow star
        geoplot(gx, lat_tr(end), lon_tr(end), 'p', 'MarkerSize', 12, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
    end
    
    % Add legend
    legend(gx, trackHandles, legendEntries, 'Location', 'bestoutside');
    
    % Export each figure
%     exportgraphics(gcf1, sprintf('%s D%d.png', d), 'Resolution', 600);
end

%% Save drifter information for comparison to ADCP *PAY ATTENTION TO OBSERVATION PERIOD*
%% Return the tracks; main_driver saves them
drifters_level_1 = trimmed;


end 

