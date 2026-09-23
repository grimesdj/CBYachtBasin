function [drifters_level_1, raw_drifters] = SHEET_drifter_post_processing(dep_times,drifters_root_dir, cfg);
proj = projcrs(32119);

%% Load the drifter release times
%dep_times = readtable(dep_fin);

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

%% Parse and convert the onboard drifter data

csvFiles = dir(fullfile(drifters_root_dir, '*.csv'));

if isempty(csvFiles);
    csvFiles = dir(fullfile(drifters_root_dir, '*.xlsx'));
else
    csvFiles = csvFiles;
end

for i = 1:numel(csvFiles);

    filename = fullfile(drifters_root_dir, csvFiles(i).name);
%     opts = detectImportOptions(filename);
%     opts = setvaropts(opts, 'Var2', 'Type', 'string');
% 
%     dataTable = readtable(filename,opts);
    dataTable = readtable(filename);
%     % 4. Convert the text column to duration with milliseconds
%     dataTable.Var2 = duration(dataTable.Var2, 'InputFormat', 'hh:mm:ss.SSS', 'Format', 'hh:mm:ss.SSS'); 
%     
    drifterID = split(csvFiles(i).name,'.');
    drifterID = split(drifterID{1},'_');
    drifterID = str2num(drifterID{end});
    
    drifters(i).ID = drifterID;
    drifters(i).File = filename;
    date = table2array(dataTable(:,1));
    date.Format = 'dd-MMM-yyyy HH:mm:ss.SSSS';
    time = table2array(dataTable(:,2));
    date_and_time = date + time;
    

    time_mat = datenum(date_and_time);

    e_time = table2array(dataTable(:,3));

    Lat = table2array(dataTable(:,4));
    Lon = table2array(dataTable(:,6));
    
    Speed = table2array(dataTable(:,9));
    
    Antenna = table2array(dataTable(:,12));
    

    
    if size(dataTable,2)>12
    trash = table2array(dataTable(:,13));
    trashFlag = ~cellfun(@isempty,trash) | Antenna<1 | Antenna>3 | Lon == 0 | Lat==0 | time_mat>datenum(today) | isnan(time_mat);
    else
    trashFlag = Antenna<1 | Antenna>3 | Lon == 0 | Lat==0 | time_mat>datenum(today) | isnan(time_mat);
    end
    % Take out the Trash!
    Lon(trashFlag)=[];
    Lat(trashFlag)=[];
    time_mat(trashFlag)=[];
    e_time(trashFlag)=[];
    date_and_time(trashFlag)=[];
    Speed(trashFlag)=[];
    
    if isempty(Lat)
        continue
    end
   
    
    % Convert the Lat and Lon to decimal degrees
    [Lat,Lon] = convertToDecimalDegrees(Lat,Lon);
    
    % Sort based on the datetime mat
    [time_mat, sorter] = sort(time_mat);
    
    Lon = Lon(sorter);
    Lat = Lat(sorter);
    e_time = e_time(sorter);
    date_and_time = date_and_time(sorter);
    Speed = Speed(sorter);
    
    % remove repeated times
    [time_mat, uni] = unique(time_mat);
    Lon = Lon(uni);
    Lat = Lat(uni);
    e_time = e_time(uni);
    date_and_time = date_and_time(uni);
    Speed = Speed(uni);
    
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
    drifters(i).Date_and_time = date_and_time;
    drifters(i).etime = e_time;
    drifters(i).Lat = Lat;
    drifters(i).Lon = Lon;
    drifters(i).East = Easting;
    drifters(i).North = Northing;

    drifters(i).Alt       = table2array(dataTable(:,8));
    drifters(i).Speed     = Speed;
    drifters(i).Angle     = table2array(dataTable(:,10));
    drifters(i).Volt      = table2array(dataTable(:,11));
    drifters(i).v_e       = V_e;
    drifters(i).v_n       = V_n;
    drifters(i).V         = V;
    

    
end

% Save the raw drifter data
raw_drifters = drifters;

%% trim the drifters to when in the water
trimmed = trimDrifterDataByDeployment(drifters,deployments);

%% Plot the trimmed drifter trajectories

colors = hsv(numel(trimmed(1,:)));

for d = 1:numel(trimmed(:,1))
    % Create a new figure + geoaxes for each deployment
    gcf1 = figure;
    gx = geoaxes;
    hold(gx,'on');
    geobasemap(gx,'satellite')
    
    % Title for this deployment
    title(gx, sprintf('%s D%d',cfg.name, d))
    
    % Reset legend storage for this deployment
    trackHandles = [];
    legendEntries = {};
    legendStartTimes = datetime.empty;
    
    % Loop through all trajectories for this deployment
    for i = 1:numel(trimmed(d,:))
        
        lat_tr = trimmed(d,i).Lat;
        lon_tr = trimmed(d,i).Lon;
        t_start = trimmed(d,i).Start_time;
        t_end = trimmed(d,i).End_time;
        
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

        % Store actual datetime for chronological sorting
        legendStartTimes(end+1) = t_start;
        
        % Plot final point as yellow star
        geoplot(gx, lat_tr(end), lon_tr(end), 'p', 'MarkerSize', 12, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
    end
    
    [~, sortIdx] = sort(legendStartTimes);

    trackHandles  = trackHandles(sortIdx);
    legendEntries = legendEntries(sortIdx);
    
    % Add legend
    legend(gx, trackHandles, legendEntries, 'Location', 'bestoutside');
    
    % Export each figure
%     exportgraphics(gcf1, fullfile(cfg.out.drifters_figures, [cfg.name, sprintf('_D%d.png', d)]), 'Resolution', 600);
%     savefig(gcf1, fullfile(cfg.out.drifters_figures,[cfg.name, sprintf('_D%d.fig', d)]));
end

%% Save stuff and veloQAQC
drifters_level_1 = trimmed;

save(fullfile(cfg.out.drifters_data, [cfg.name '_raw_drifters.mat']), 'raw_drifters');

end