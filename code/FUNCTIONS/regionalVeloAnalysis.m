%% Regional velocity analysis
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
    'JuneJamboree'...
    'SummerFlood'};

%% Intitialize master struct

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
           
       end
       
       %% Make column vectors
       
       Time =  Time(:);
       V = V(:);
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
       RegionalData(counter).Region = Regions;
       RegionalData(counter).TidePhase = TidePhase;
       RegionalData(counter).Lat = Lat;
       RegionalData(counter).Lon = Lon;

        end
    end
end

%% Calculate the average velcity in each region

load('regions.mat');

nRegions = numel(regions);

RegionStats = struct();

for r = 1:nRegions;
    
    V_region = [];
    
    for k = 1:numel(RegionalData)
        
        idx = RegionalData(k).Region == r;
        
        V_region = [V_region; RegionalData(k).V(idx)];
        
    end
    
    RegionStats(r).Region = r;
    RegionStats(r).MeanVelocity = mean(V_region,'omitnan');
    RegionStats(r).StdVelocity = std(V_region, 'omitnan');
    RegionStats(r).N = sum(~isnan(V_region));
    
end

%% Plot regional mean velocities

figure

gx = geoaxes;
hold(gx,'on')

geobasemap(gx,'satellite')

legend_handles = gobjects(nRegions,1);
legend_labels = strings(nRegions,1);

for r = 1:nRegions

    lat = regions(r).Lat;
    lon = regions(r).Lon;

    % Plot region boundary
    hRegion = geoplot(gx, ...
        lat, ...
        lon, ...
        'LineWidth',2);

    % Save handle for legend
    legend_handles(r) = hRegion;

    % Region center
    lat_center = mean(lat(1:end-1),'omitnan');
    lon_center = mean(lon(1:end-1),'omitnan');

    % Display only region ID on map
    text(gx, ...
        lat_center, ...
        lon_center, ...
        sprintf('%d',r), ...
        'FontSize',14, ...
        'FontWeight','bold', ...
        'Color','w', ...
        'HorizontalAlignment','center');

    % Legend label
    legend_labels(r) = sprintf( ...
        'Region %d: %.2f m/s  (N = %d)', ...
        r, ...
        RegionStats(r).MeanVelocity, ...
        RegionStats(r).N);
end

legend( ...
    legend_handles, ...
    legend_labels, ...
    'Location','bestoutside');

title(sprintf( ...
    'Mean Drifter Velocity by Region'))

%% Save overall regional velocity figure

file_name = sprintf( ...
    'RegionalMeanScalarVelo_%ds', ...
    temporal_res);

exportgraphics(gcf, ...
    fullfile(save_dir,[file_name '.png']), ...
    'Resolution',300);

savefig(gcf, ...
    fullfile(save_dir,[file_name '.fig']));

%% Display regional results

Region = [RegionStats.Region]';
MeanVelocity = [RegionStats.MeanVelocity]';
StdVelocity = [RegionStats.StdVelocity]';
N = [RegionStats.N]';

RegionTable = table(Region,MeanVelocity,StdVelocity,N);

disp(RegionTable)

%% Calculate average velocity by region and tidal phase

phases = {'High','Low','Rising','Falling'};

RegionPhaseStats = struct();

for p = 1:numel(phases)

    phase = phases{p};

    for r = 1:nRegions

        V_region_phase = [];

        for k = 1:numel(RegionalData)

            idx = ...
                RegionalData(k).Region == r & ...
                RegionalData(k).TidePhase == phase;

            V_region_phase = [ ...
                V_region_phase; ...
                RegionalData(k).V(idx)];

        end

        RegionPhaseStats.(phase)(r).Region = r;

        RegionPhaseStats.(phase)(r).MeanVelocity = ...
            mean(V_region_phase,'omitnan');

        RegionPhaseStats.(phase)(r).StdVelocity = ...
            std(V_region_phase,'omitnan');

        RegionPhaseStats.(phase)(r).N = ...
            sum(~isnan(V_region_phase));

    end

end

%% Display region and tidal phase results

for p = 1:numel(phases)

    phase = phases{p};

    Region = [RegionPhaseStats.(phase).Region]';
    MeanVelocity = [RegionPhaseStats.(phase).MeanVelocity]';
    StdVelocity = [RegionPhaseStats.(phase).StdVelocity]';
    N = [RegionPhaseStats.(phase).N]';

    T = table(Region,MeanVelocity,StdVelocity,N);

    fprintf('\n%s TIDE\n',upper(phase))
    disp(T)

end

%% Plot regional mean velocities for each tidal phase

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

        % Plot region boundary
        hRegion = geoplot(gx, ...
            lat, ...
            lon, ...
            'LineWidth',2);

        % Save handle for legend
        legend_handles(r) = hRegion;

        % Region center
        lat_center = mean(lat(1:end-1),'omitnan');
        lon_center = mean(lon(1:end-1),'omitnan');

        % Display only region ID on map
        text(gx, ...
            lat_center, ...
            lon_center, ...
            sprintf('%d',r), ...
            'FontSize',14, ...
            'FontWeight','bold', ...
            'Color','w', ...
            'HorizontalAlignment','center');

        % Legend label using region + phase stats
        legend_labels(r) = sprintf( ...
            'Region %d: %.2f m/s  (N = %d)', ...
            r, ...
            RegionPhaseStats.(phase)(r).MeanVelocity, ...
            RegionPhaseStats.(phase)(r).N);

    end

    legend( ...
        legend_handles, ...
        legend_labels, ...
        'Location','bestoutside');

    title(sprintf( ...
        '%s Tide Mean Drifter Velocity by Region', ...
        phase))
    
    %% Save tidal-phase regional velocity figure

    file_name = sprintf( ...
        '%sTide_RegionalMeanScalarVelo_%ds', ...
        phase, ...
        temporal_res);

    exportgraphics(gcf, ...
        fullfile(save_dir,[file_name '.png']), ...
        'Resolution',300);

    savefig(gcf, ...
        fullfile(save_dir,[file_name '.fig']));

end






