function [drifters_level_1, raw_drifters] = drifter_post_processing(cfg);
dep_times = readtable(cfg.drifters_dep_times);
drifters_root_dir = cfg.drifters_raw_dir;

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
    opts = detectImportOptions(filename);
    opts = setvaropts(opts, 'Var2', 'Type', 'string');

    dataTable = readtable(filename,opts);
%     dataTable = readtable(filename);
    % 4. Convert the text column to duration with milliseconds
    dataTable.Var2 = duration(dataTable.Var2, 'InputFormat', 'hh:mm:ss.SSS', 'Format', 'hh:mm:ss.SSS'); 
    
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

    

    % 0) Find and remove the zeros from before the GPS fix
%     firstFix = find(Lon ~= 0 & Lat~=0 & time_mat<datenum(today), 1, 'first');
%     
%     Lon(1:firstFix)=[];
%     Lat(1:firstFix)=[];
%     time_mat(1:firstFix)=[];
%     e_time(1:firstFix)=[];
%     date_and_time(1:firstFix)=[];
% 
%     
%     % Now do last fix
%     lastFix = find(Lon ~= 0 & Lat~=0 & time_mat<datenum(today) , 1, 'last');
%     
%     Lon(lastFix:end)=[];
%     Lat(lastFix:end)=[];
%     time_mat(lastFix:end)=[];
%     e_time(lastFix:end)=[];
%     date_and_time(lastFix:end)=[];
    
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
    
    
    
    % 1) Create a 1Hz time vector starting from first fix time...
    e_time = e_time - e_time(1);
%     time_1hz = (1:1:e_time(end-1))';
    time_1hz = (1:1:size(e_time))';
    time_1hz_mat = time_mat(1)+time_1hz/86400;
    time_sec = seconds(date_and_time-date_and_time(1));
    
    % 2) find remaining zero indices
    
    zeroIndices=(Lon == 0 | Lat==0);
    
    % 3) interpolate Lat/Lon/e_time/date_and_time/time_mat to 1Hz excluding
    % the zero indices
   
    Lat = interp1(time_mat(~zeroIndices),Lat(~zeroIndices),time_1hz_mat);
    Lon = interp1(time_mat(~zeroIndices),Lon(~zeroIndices),time_1hz_mat);
    Easting = interp1(time_mat(~zeroIndices),Easting(~zeroIndices),time_1hz_mat);
    Northing = interp1(time_mat(~zeroIndices),Northing(~zeroIndices),time_1hz_mat);
    Speed = interp1(time_mat(~zeroIndices),Speed(~zeroIndices),time_1hz_mat);
    
    e_time = time_1hz;
    time_mat = time_1hz_mat;
    date_and_time = datetime(time_mat,'convertFrom','datenum');
    
    % 4) calculate V and flag
    
    % Calculate velocities
    delta_east = diff(Easting);
    delta_north = diff(Northing);

    delta_time = diff(e_time);

    V_e = delta_east./delta_time;
    V_n = delta_north./delta_time;
    V = sqrt((V_e.^2) + (V_n.^2));
    
    % Flag large jumps in V
    flag = V > 2.5;
    
    % 5) use the flag to interpolate across E/N/Lat/Lon/Speed
    
    Easting(flag) = interp1(time_1hz(~flag),Easting(~flag),time_1hz(flag));
    Northing(flag) = interp1(time_1hz(~flag),Northing(~flag),time_1hz(flag));
    Lat(flag) = interp1(time_1hz(~flag),Lat(~flag),time_1hz(flag));
    Lon(flag) = interp1(time_1hz(~flag),Lon(~flag),time_1hz(flag));
    % keep raw speeds
    %Speed(flag) = interp1(time_1hz(~flag),Speed(~flag),time_1hz(flag));
    
    % QAQC fig
%     figure;
%     geoscatter(Lat, Lon, 10, e_time, 'filled')
%     colorbar
    
    % 6) Smooth the delta_trajectories and re-calculate V
    flt = hanning(61); 
    flt = flt/sum(flt);
    
    
    delta_east = diff(Easting);
    delta_north = diff(Northing);
    
%     delta_east  = conv(delta_east,flt,'same');
%     delta_north = conv(delta_north,flt,'same');
% 
%     tmp = delta_east + sqrt(-1)*delta_north;
%     tmp = conv(tmp,flt,'same');
%     
%     delta_east = real(tmp);
%     delta_north = imag(tmp);
%     
    delta_time = diff(e_time);
% 
    V_e = delta_east./delta_time;
    V_n = delta_north./delta_time;
    V = sqrt((V_e.^2) + (V_n.^2));
%     

    %     
%     flt = hanning(61); 
%     flt = flt/sum(flt);
%     V = conv(V,flt,'same');
%     V_e = conv(V_e,flt,'same');
%     V_n = conv(V_n,flt,'same');
%     
%     % 7) after computing V = East(2:end)-East(1:end-1); must average
%     % Lon/Lat/Time to Lon = 0.5*(Lon(1:end-1)+Lon(2:end));
    Lon = 0.5*(Lon(1:end-1)+Lon(2:end));
    Lat = 0.5*(Lat(1:end-1)+Lat(2:end));
    Easting = 0.5*(Easting(1:end-1)+Easting(2:end));
    Northing = 0.5*(Northing(1:end-1)+Northing(2:end));
    Speed = 0.5*(Speed(1:end-1)+Speed(2:end));
    time_mat = 0.5*(time_mat(1:end-1)+time_mat(2:end));
    e_time = 0.5*(e_time(1:end-1)+e_time(2:end));
    time_1hz = 0.5*(time_1hz(1:end-1)+time_1hz(2:end));
    date_and_time = datetime(time_mat,'convertfrom','datenum');
    
% %     % Make the 5 min res data
% % %     time_1hz_5min = time_1hz(300+1:end - 300);
% %     time_1hz_5min = time_1hz(151:end - 150);
% % 
% %     delta_time_5min = time_1hz(300+1:end) - time_1hz(1:end-300);
% % %     easting_5min = Easting(300+1:end - 300);
% % %     northing_5min = Northing(300+1:end - 300);
% %     easting_5min = Easting(150+1:end - 150);
% %     northing_5min = Northing(150+1:end - 150);
% %     delta_east_5min = Easting(300+1:end) - Easting(1:end-300);
% %     delta_north_5min = Northing(300+1:end) - Northing(1:end-300);
% %     
% % %     delta_east_5min = conv(delta_east_5min,flt,'same');
% % %     delta_north_5min = conv(delta_north_5min,flt,'same');
% %     
% %     V_e_5min = delta_east_5min./delta_time_5min;
% %     V_n_5min = delta_north_5min./delta_time_5min;
% %     V_5min = sqrt((V_e_5min.^2) + (V_n_5min.^2));
% %     
% % %     time_1hz_5min = time_1hz_5min(1:end-1) + delta_time_5min/2;
% %     
% %     
% %     figure;
% %     plot(e_time,V,'k');
% %     hold on;
% %     plot(e_time,Speed,'r');
% %     plot(time_1hz_5min,V_5min,'m');
% %     legend('V 1hz smoothed de/dn','raw doppler speed','V 5min')
%     
%     Lon(zeroIndices) = nan;
%     Lat(zeroIndices) = nan;
%     time_mat(zeroIndices) = nan;
%     date_and_time = datetime(time_mat,'convertfrom','datenum');
% 
%     
%     
%     
%     
%     % recalculate Lat/Lon/E/N/Time/zeroIndices to match V/Ve/Vn
%     zeroIdx = 0.5*(zeroIndices(1:end-1) + zeroIndices(2:end));
%     
%     % Interpolate across these large jumps
%     
%     Lat(flag) = interp1(e_time(~flag),Lat(~flag),e_time(flag));
%     Lon(flag) = interp1(e_time(~flag),Lon(~flag),e_time(flag));
%     V(flag) = interp1(e_time(~flag),V(~flag),e_time(flag));
%     
%     
%     flt = hanning(61); 
%     flt = flt/sum(flt);
%     V = conv(V,flt,'same');
%     V_e = conv(V_e,flt,'same');
%     V_n = conv(V_n,flt,'same');
    
    % QA/QC figure
%     figure;
%     plot(date_and_time, table2array(dataTable(:,9), 'r'));  
%     hold on
%     plot(date_and_time(1:end-1), V,'k');
%     ylim([0 3]);

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

% Make some histograms


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


%% Convert 'trimmed' structure into 1.5 hr time chunks
% time_chunks = sortDrifterDataByTime(trimmed);
% %% Fix the colors so they are consistent for time chunk plots
% allIDs = [];
% 
% for i = 1:numel(time_chunks)
%     allIDs = [allIDs, [time_chunks(i).Drifters.ID]];
% end
% 
% uniqueIDs = unique(allIDs);
% 
% % One color per unique drifter
% cmap = lines(numel(uniqueIDs));
% 
% IDcolor = containers.Map('KeyType','double','ValueType','any');
% 
% for k = 1:numel(uniqueIDs)
%     IDcolor(uniqueIDs(k)) = cmap(k,:);
% end

%% Plot the trimmed trajectories in 1.5 hour chunks

% for i = 1:numel(time_chunks)
%     
%     gcf2 = figure;
%     gx = geoaxes;
%     hold(gx,'on');
%     
%     % Reset legend storage for this deployment
%     trackHandles = [];
%     legendEntries = {};
%     legendStartTimes = datetime.empty;
%     
%     drifters = time_chunks(i).Drifters;
%     
%     for j = 1:numel(drifters)
%         
%         ID = drifters(j).ID;
%         
%         % Find the color for this ID
%         thisColor = IDcolor(ID);
%         
%         h = geoplot(gx, drifters(j).Lat,drifters(j).Lon,'LineWidth', 1.5,'Color',thisColor);
%         
%         % Plot final point as yellow star
%         geoplot(gx, drifters(j).Lat(end),drifters(j).Lon(end), 'p', 'MarkerSize', 12, ...
%             'MarkerEdgeColor', 'k','MarkerFaceColor',thisColor);
% 
%         % Store for legend
%         trackHandles(end+1) = h; 
%         % legendEntries{end+1} = sprintf('d%d_t%d (%.1f h)', d, i, dt_hr); 
%         legendEntries{end+1} = sprintf('D%d | %s–%s', drifters(j).ID, datestr(drifters(j).Start_time,'mm/dd HH:MM'), datestr(drifters(j).End_time,'HH:MM'));
%         % Store actual datetime for chronological sorting
%         legendStartTimes(end+1) = drifters(j).Start_time;
%         
%     end
%     
%     geobasemap(gx,'satellite');
%     
% 
%     title(gx,sprintf('%s - %s', ...
%         string(time_chunks(i).Start_time), ...
%         string(time_chunks(i).End_time)))
%     
%     [~, sortIdx] = sort(legendStartTimes);
% 
%     trackHandles  = trackHandles(sortIdx);
%     legendEntries = legendEntries(sortIdx);
%     
%     
%     % Add legend
%     legend(gx, trackHandles, legendEntries, 'Location', 'bestoutside');
%     
%     % Export each figure
% %     exportgraphics(gcf2, fullfile(cfg.out.drifters_figures, [cfg.name, sprintf('_D%d.png', d)]), 'Resolution', 600);
% %     savefig(gcf2, fullfile(cfg.out.drifters_figures,[cfg.name, sprintf('_D%d.fig', d)]));
%     
%     
%     
%     
% end

%% Save stuff and veloQAQC
drifters_level_1 = trimmed;

save(fullfile(cfg.out.drifters_data, [cfg.name '_raw_drifters.mat']), 'raw_drifters');

% [drifters_QAQC] = veloQAQC(cfg);




end