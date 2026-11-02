%% Regional Velocity Vectors
clear all;
close all;

%% Settings

temporal_res = 8;
% options: [1 2 4 8 16 32 64 128 256]

root_data = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'; % point the script to your CB_YachBasin data folder
save_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES\REGIONAL_ANALYSIS';

obsFiles = { ...
    'CB_YachtBasin\MarshMadness\PROCESSED_DATA\DRIFTERS\MarshMadness_TemporalRes_RMSE.mat'
    'CB_YachtBasin\MarshMayhem\PROCESSED_DATA\DRIFTERS\MarshMayhem_TemporalRes_RMSE.mat'
    'CB_YachtBasin\FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_100825_TemporalRes_RMSE.mat'
    'CB_YachtBasin\FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_101025_TemporalRes_RMSE.mat'
    'CB_YachtBasin\NovDep\PROCESSED_DATA\DRIFTERS\NovDep_TemporalRes_RMSE.mat'
    'CB_YachtBasin\JuneJamboree\PROCESSED_DATA\DRIFTERS\JuneJamboree_TemporalRes_RMSE.mat'
    'CB_YachtBasin\SummerFlood\PROCESSED_DATA\DRIFTERS\SummerFlood_TemporalRes_RMSE.mat'};

names = { ...
    'MarshMadness', ...
    'MarshMayhem', ...
    'FallFrolic_100825', ...
    'FallFrolic_101025', ...
    'NovDep', ...
    'JuneJamboree', ...
    'SummerFlood'};

%% Initialize master struct

RegionalData = struct();
counter = 0;

%% Load all the observation periods

for i = 1:numel(obsFiles)
    
    filename = fullfile(root_data,obsFiles{i});
    
    S = load(filename);
    Results = S.Results;
    
    obs = names{i};
    
    %% Loop through the drifters in this obs period
    
    for d = 1:size(Results,1)
        
        for j = 1:size(Results,2)
            
       if isempty(Results(d,j).ID) 
           continue
       end
       
       %% Obtain desired temporal resolution
       
       if temporal_res == 1
           
           Time = Results(d,j).Time;
           V = Results(d,j).V_1Hz;
           Regions = Results(d,j).Region;
           TidePhase = Results(d,j).TidePhase;
           Lat = Results(d,j).Lat;
           Lon = Results(d,j).Lon;
           Easting = Results(d,j).Easting;
           Northing = Results(d,j).Northing;
           Speed = Results(d,j).Speed;
           V_e = Results(d,j).V_e_1Hz;
           V_n = Results(d,j).V_n_1Hz;
           
           
       else
           
           res_idx = find([Results(d,j).TemporalRes.Seconds] == temporal_res,1);
           
           if isempty(res_idx)
               continue
           end
           
           TR = Results(d,j).TemporalRes(res_idx);
           
           Time = TR.Time;
           V = TR.V;
           Regions = TR.Region;
           TidePhase = TR.TidePhase;
           Lat = TR.Lat;
           Lon = TR.Lon;
           Easting = TR.Easting;
           Northing = TR.Northing;
           Speed = TR.Speed;
           V_e = TR.V_e;
           V_n = TR.V_n;
           
       end
       
       %% Make column vectors
       
       Time =  Time(:);
       V = V(:);
       V_e = V_e(:);
       V_n = V_n(:);
       Regions = Regions(:);
       TidePhase = string(TidePhase(:));
       Lat = Lat(:);
       Lon = Lon(:);
       Easting = Easting(:);
       Northing = Northing(:);
       Speed = Speed(:);
       
       %% Store this trajectory
       
       counter = counter + 1;
       
       RegionalData(counter).Observation = obs;
       RegionalData(counter).ID = Results(d,j).ID;
       RegionalData(counter).Time = Time;
       RegionalData(counter).V = V;
       RegionalData(counter).V_e = V_e;
       RegionalData(counter).V_n = V_n;
       RegionalData(counter).Region = Regions;
       RegionalData(counter).TidePhase = TidePhase;
       RegionalData(counter).Lat = Lat;
       RegionalData(counter).Lon = Lon;

        end
    end
end

%% Calculate the average velocity in each region

load('regions.mat');

nRegions = numel(regions);

RegionStats = struct();

for r = 1:nRegions

    V_region   = [];
    Ve_region  = [];
    Vn_region  = [];

    for k = 1:numel(RegionalData)

        idx = RegionalData(k).Region == r;

        V_region  = [V_region;  RegionalData(k).V(idx)];
        Ve_region = [Ve_region; RegionalData(k).V_e(idx)];
        Vn_region = [Vn_region; RegionalData(k).V_n(idx)];

    end

    RegionStats(r).Region = r;

    % Mean of velocity magnitude
    RegionStats(r).MeanVelocity = mean(V_region,'omitnan');

    % Mean vector components
    RegionStats(r).MeanVe = mean(Ve_region,'omitnan');
    RegionStats(r).MeanVn = mean(Vn_region,'omitnan');

    % Magnitude of the mean vector
    RegionStats(r).MeanVectorMagnitude = hypot( ...
        RegionStats(r).MeanVe, ...
        RegionStats(r).MeanVn);

    RegionStats(r).StdVelocity = std(V_region,'omitnan');

    RegionStats(r).N = sum(~isnan(V_region));

end

%% Plot regional mean velocity vectors

figure

gx = geoaxes;
hold(gx,'on')

geobasemap(gx,'satellite')

legend_handles = gobjects(nRegions,1);
legend_labels = strings(nRegions,1);

% Scale velocity vectors for visualization
velocity_scale = 3000;   
% Increase/decrease this to change arrow lengths

for r = 1:nRegions

    lat = regions(r).Lat;
    lon = regions(r).Lon;

    %% Plot region boundary

    hRegion = geoplot(gx, ...
        lat, ...
        lon, ...
        'LineWidth',2);

    legend_handles(r) = hRegion;

    %% Region center

    lat_center = mean(lat(1:end-1),'omitnan');
    lon_center = mean(lon(1:end-1),'omitnan');

    %% Mean velocity components

    Ve = RegionStats(r).MeanVe;
    Vn = RegionStats(r).MeanVn;

    Vmag = RegionStats(r).MeanVectorMagnitude;

    %% Convert velocity vector to geographic displacement

    meters_per_deg_lat = 111320;
    meters_per_deg_lon = 111320*cosd(lat_center);

    % End point of velocity vector
    lat_end = lat_center + ...
        (Vn * velocity_scale) / meters_per_deg_lat;

    lon_end = lon_center + ...
        (Ve * velocity_scale) / meters_per_deg_lon;

    %% Plot vector shaft

    geoplot(gx, ...
        [lat_center lat_end], ...
        [lon_center lon_end], ...
        'w-', ...
        'LineWidth',2.5);

    %% Make arrow head

    % Vector displacement in meters
    dx = Ve * velocity_scale;
    dy = Vn * velocity_scale;

    L = hypot(dx,dy);

    if L > 0

        % Unit vector in direction of flow
        ux = dx/L;
        uy = dy/L;

        % Unit vector perpendicular to flow
        px = -uy;
        py = ux;

        % Arrowhead dimensions
        head_length = 0.20 * L;
        head_width  = 0.10 * L;

        % Arrowhead points in meters relative to arrow endpoint
        x1 = dx - head_length*ux + head_width*px;
        y1 = dy - head_length*uy + head_width*py;

        x2 = dx - head_length*ux - head_width*px;
        y2 = dy - head_length*uy - head_width*py;

        % Convert arrowhead points back to lat/lon
        lat_head1 = lat_center + y1/meters_per_deg_lat;
        lon_head1 = lon_center + x1/meters_per_deg_lon;

        lat_head2 = lat_center + y2/meters_per_deg_lat;
        lon_head2 = lon_center + x2/meters_per_deg_lon;

        % Plot arrowhead
        geoplot(gx, ...
            [lat_head1 lat_end lat_head2], ...
            [lon_head1 lon_end lon_head2], ...
            'w-', ...
            'LineWidth',2.5);

    end

    %% Region number

    text(gx, ...
        lat_center, ...
        lon_center, ...
        sprintf('%d',r), ...
        'FontSize',14, ...
        'FontWeight','bold', ...
        'Color','w', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');

    %% Legend label

    legend_labels(r) = sprintf( ...
        'Region %d: |V| = %.3f m/s', ...
        r, ...
        Vmag);

end

%% Legend

legend( ...
    legend_handles, ...
    legend_labels, ...
    'Location','bestoutside');

%% Title

title(sprintf( ...
    'Mean Drifter Velocity Vectors by Region'));

%% Save overall regional velocity figure

file_name = sprintf( ...
    'RegionalMeanVelocityVectors_%ds', ...
    temporal_res);

exportgraphics(gcf, ...
    fullfile(save_dir,[file_name '.png']), ...
    'Resolution',300);

savefig(gcf, ...
    fullfile(save_dir,[file_name '.fig']));

%% Calculate average velocity by region and tidal phase

phases = {'High','Low','Rising','Falling'};

RegionPhaseStats = struct();

for p = 1:numel(phases)

    phase = phases{p};

    for r = 1:nRegions

        V_region_phase  = [];
        Ve_region_phase = [];
        Vn_region_phase = [];

        for k = 1:numel(RegionalData)

            idx = ...
                RegionalData(k).Region == r & ...
                RegionalData(k).TidePhase == phase;

            % Velocity magnitude
            V_region_phase = [ ...
                V_region_phase; ...
                RegionalData(k).V(idx)];

            % East velocity component
            Ve_region_phase = [ ...
                Ve_region_phase; ...
                RegionalData(k).V_e(idx)];

            % North velocity component
            Vn_region_phase = [ ...
                Vn_region_phase; ...
                RegionalData(k).V_n(idx)];

        end

        %% Store statistics

        RegionPhaseStats.(phase)(r).Region = r;

        % Mean speed
        RegionPhaseStats.(phase)(r).MeanVelocity = ...
            mean(V_region_phase,'omitnan');

        % Mean vector components
        RegionPhaseStats.(phase)(r).MeanVe = ...
            mean(Ve_region_phase,'omitnan');

        RegionPhaseStats.(phase)(r).MeanVn = ...
            mean(Vn_region_phase,'omitnan');

        % Magnitude of mean velocity vector
        RegionPhaseStats.(phase)(r).MeanVectorMagnitude = ...
            hypot( ...
                RegionPhaseStats.(phase)(r).MeanVe, ...
                RegionPhaseStats.(phase)(r).MeanVn);

        % Standard deviation of speed
        RegionPhaseStats.(phase)(r).StdVelocity = ...
            std(V_region_phase,'omitnan');

        % Number of observations
        RegionPhaseStats.(phase)(r).N = ...
            sum(~isnan(V_region_phase));

    end

end

%% Plot regional mean velocity vectors for each tidal phase

for p = 1:numel(phases)

    phase = phases{p};

    figure

    gx = geoaxes;
    hold(gx,'on')

    geobasemap(gx,'satellite')

    legend_handles = gobjects(nRegions,1);
    legend_labels = strings(nRegions,1);

    for r = 1:nRegions

        lat = regions(r).Lat;
        lon = regions(r).Lon;

        %% Plot region boundary

        hRegion = geoplot(gx, ...
            lat, ...
            lon, ...
            'LineWidth',2);

        legend_handles(r) = hRegion;

        %% Region center

        lat_center = mean(lat(1:end-1),'omitnan');
        lon_center = mean(lon(1:end-1),'omitnan');

        %% Mean velocity components for this tidal phase

        Ve = RegionPhaseStats.(phase)(r).MeanVe;
        Vn = RegionPhaseStats.(phase)(r).MeanVn;

        Vmag = RegionPhaseStats.(phase)(r).MeanVectorMagnitude;

        %% Convert velocity vector to geographic displacement

        meters_per_deg_lat = 111320;
        meters_per_deg_lon = 111320*cosd(lat_center);

        lat_end = lat_center + ...
            (Vn * velocity_scale) / meters_per_deg_lat;

        lon_end = lon_center + ...
            (Ve * velocity_scale) / meters_per_deg_lon;

        %% Plot vector

        if ~isnan(Ve) && ~isnan(Vn) && Vmag > 0

            % Arrow shaft
            geoplot(gx, ...
                [lat_center lat_end], ...
                [lon_center lon_end], ...
                'w-', ...
                'LineWidth',2.5);

            %% Arrowhead

            dx = Ve * velocity_scale;
            dy = Vn * velocity_scale;

            L = hypot(dx,dy);

            ux = dx/L;
            uy = dy/L;

            px = -uy;
            py = ux;

            head_length = 0.20 * L;
            head_width  = 0.10 * L;

            x1 = dx - head_length*ux + head_width*px;
            y1 = dy - head_length*uy + head_width*py;

            x2 = dx - head_length*ux - head_width*px;
            y2 = dy - head_length*uy - head_width*py;

            lat_head1 = lat_center + ...
                y1/meters_per_deg_lat;

            lon_head1 = lon_center + ...
                x1/meters_per_deg_lon;

            lat_head2 = lat_center + ...
                y2/meters_per_deg_lat;

            lon_head2 = lon_center + ...
                x2/meters_per_deg_lon;

            geoplot(gx, ...
                [lat_head1 lat_end lat_head2], ...
                [lon_head1 lon_end lon_head2], ...
                'w-', ...
                'LineWidth',2.5);

        end

        %% Region number

        text(gx, ...
            lat_center, ...
            lon_center, ...
            sprintf('%d',r), ...
            'FontSize',14, ...
            'FontWeight','bold', ...
            'Color','w', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom');

        %% Legend label

        legend_labels(r) = sprintf( ...
            'Region %d: |V| = %.3f m/s', ...
            r, ...
            Vmag);

    end

    %% Legend

    legend( ...
        legend_handles, ...
        legend_labels, ...
        'Location','bestoutside');

    %% Title

    title(sprintf( ...
        '%s Tide Mean Drifter Velocity Vectors by Region', ...
        phase))
    
    %% Save tidal-phase regional velocity figure

file_name = sprintf( ...
    '%sTide_RegionalMeanVelocityVectors_%ds', ...
    phase, ...
    temporal_res);

exportgraphics(gcf, ...
    fullfile(save_dir,[file_name '.png']), ...
    'Resolution',300);

savefig(gcf, ...
    fullfile(save_dir,[file_name '.fig']));

end

