clear all;
close all;

%% Load the stationary drifter data

filename = fullfile('C:\Users\bcm3620\OneDrive - UNC-Wilmington\CB_YachtBasin\StationaryDep\PROCESSED_DATA\DRIFTERS\StationaryDep_drifters_level_1.mat');

% filename = fullfile('C:\Users\Ben\OneDrive - UNC-Wilmington\CB_YachtBasin\StationaryDep\PROCESSED_DATA\DRIFTERS\StationaryDep_drifters_level_1.mat');


S = load(filename);

drifters_level_1 = S.drifters_level_1;

%% Temporal Res settings

res_sec = [2 4 8 16 32 64 128 256];

max_halfLag = max(res_sec)/2;

[nDeployments,nDrifters] = size(drifters_level_1);

%% Preallocate output structure

Results = struct(...
    'ID', [], ...
    'Resolution_Sec', [], ...
    'RMSE', [], ...
    'RMSE_v_Speed', [], ...
    'TemporalRes', [], ...
    'Time', [], ...
    'V_1Hz', [], ...
    'Speed', []);

Results = repmat(Results,nDeployments,nDrifters);

%% Loop through the drifters

for d = 1:nDeployments

    for j = 1:nDrifters

        D = drifters_level_1(d,j);

        %% Skip empty deployment/drifter combinations

        if isempty(D.e_time) || ...
           isempty(D.East) || ...
           isempty(D.North)

            continue

        end


        %% Determine drifter ID

        if isfield(D,'ID') && ~isempty(D.ID)

            drifterID = D.ID;

        else

            drifterID = j;

        end


        %% Extract 1 Hz data

        time_1hz = D.e_time;

        Speed = D.Speed;

        % Convert GPS Doppler speed from knots to m/s
        Speed = Speed * 0.51444;

        Easting  = D.East;
        Northing = D.North;
        Lat      = D.Lat;
        Lon      = D.Lon;
%         TideP    = D.TidePhase;
%         Regions  = D.Regions;
        V        = D.V;


        %% Make everything column vectors

        time_1hz = time_1hz(:);
        Speed    = Speed(:);
        Easting  = Easting(:);
        Northing = Northing(:);
        Lat      = Lat(:);
        Lon      = Lon(:);
        V        = V(:);


        %% Check trajectory length

        if length(time_1hz) <= 2*max_halfLag

            fprintf( ...
                'Skipping Deployment %d | Drifter %d: trajectory too short\n', ...
                d,drifterID);

            continue

        end


        %% Common center indices for ALL resolutions

        idx_center = ...
            (1 + max_halfLag): ...
            (length(time_1hz) - max_halfLag);


        time_1hz_centered = time_1hz(idx_center);

        Speed_centered = Speed(idx_center);

        V_centered = V(idx_center);

        Easting_centered = Easting(idx_center);

        Northing_centered = Northing(idx_center);

        Lat_centered = Lat(idx_center);

        Lon_centered = Lon(idx_center);

%         tidalPhase_centered = TideP(idx_center);
% 
%         Regions_centered = Regions(idx_center);


        %% Calculate temporal-resolution velocities

        TemporalRes = struct([]);

        for r = 1:length(res_sec)

            lag = res_sec(r);

            halfLag = lag/2;


            % Indices on either side of center point
            idx_before = idx_center - halfLag;

            idx_after = idx_center + halfLag;


            % Change across full time window
            delta_time = ...
                time_1hz(idx_after) - ...
                time_1hz(idx_before);

            delta_east = ...
                Easting(idx_after) - ...
                Easting(idx_before);

            delta_north = ...
                Northing(idx_after) - ...
                Northing(idx_before);


            % Velocity components
            new_V_e = delta_east ./ delta_time;

            new_V_n = delta_north ./ delta_time;


            % Position Derived Velocity magnitude
            new_V = sqrt( ...
                new_V_e.^2 + ...
                new_V_n.^2 );

            % Smooth Doppler Speed over same temporal resolution

            Speed_smoothed = movmean( ...
                Speed, ...
                lag, ...
                'omitnan');

            Speed_smoothed_centered = ...
                Speed_smoothed(idx_center);


            %% Store temporal-resolution data

            TemporalRes(r).Seconds = lag;

            TemporalRes(r).Time = ...
                time_1hz_centered;

            TemporalRes(r).Lat = ...
                Lat_centered;

            TemporalRes(r).Lon = ...
                Lon_centered;

            TemporalRes(r).Easting = ...
                Easting_centered;

            TemporalRes(r).Northing = ...
                Northing_centered;

            TemporalRes(r).V_e = ...
                new_V_e;

            TemporalRes(r).V_n = ...
                new_V_n;

            TemporalRes(r).V = ...
                new_V;

            TemporalRes(r).Speed = ...
                Speed_smoothed_centered;

        end
        
        %% Calculate RMSE relatice to max temporal res

idx_maxres = find([TemporalRes.Seconds] == max(res_sec),1);

V_maxres = TemporalRes(idx_maxres).V;

% 1Hz vs maximum temporal res
RMSE_1s = sqrt(mean((V_centered - V_maxres).^2, 'omitnan'));

% 1 Hz vs GPS Doppler speed
RMSE_1s_Speed = sqrt(mean((V_centered - Speed_centered).^2, 'omitnan'));

% 1 Hz Doppler vs max res GPS Doppler speed
RMSE_Speed_v_Speed = sqrt(mean((Speed_centered - TemporalRes(end).Speed).^2, 'omitnan'));

%% Preallocate RMSE arrays

RMSE_res = nan(size(res_sec));
RMSE_Speed = nan(size(res_sec));
RMSE_Speed_Smoothed = nan(size(res_sec));

%% RMSE for every other temporal res

for r = 1:length(res_sec)
    
   % Position-derived velocity vs max temporal resolution
   RMSE_res(r) = sqrt(mean( ...
                    (TemporalRes(r).V - V_maxres).^2, ...
                    'omitnan'));
                
   % Position-derived veolicty vs. Doppler Speed at 1Hz
   RMSE_Speed(r) = sqrt(mean( ...
                    (TemporalRes(r).V - Speed_centered).^2, ...
                    'omitnan'));
               
   % Position-derived velocity vs. Doppler Speed at corresponding temporal
   % res
   RMSE_Speed_Smoothed(r) = sqrt(mean( ...
                    (TemporalRes(r).V - TemporalRes(r).Speed).^2, ...
                    'omitnan'));

   %  1Hz Doppler Speed vs. Doppler Speed at corresponding temporal
   % res
   RMSE_Doppler_res(r) = sqrt(mean( ...
       (TemporalRes(r).Speed - TemporalRes(end).Speed).^2, ...
       'omitnan'));
    
    
end

        %% Store results for this deployment/drifter

        Results(d,j).ID = drifterID;

        Results(d,j).Resolution_Sec = ...
            [1 res_sec];

        Results(d,j).RMSE = ...
            [RMSE_1s RMSE_res];

        Results(d,j).RMSE_v_Speed = ...
            [RMSE_1s_Speed RMSE_Speed];
        
        Results(d,j).RMSE_v_Speed_Smoothed = ...
                [RMSE_1s_Speed RMSE_Speed_Smoothed];

        Results(d,j).RMSE_Speed_v_Speed_Smoothed = ...
            [RMSE_Speed_v_Speed RMSE_Doppler_res];

        Results(d,j).TemporalRes = ...
            TemporalRes;

        Results(d,j).Time = ...
            time_1hz_centered;
        
        Results(d,j).Lat = ...
                Lat_centered;
            
        Results(d,j).Lon = ...
                Lon_centered;

        Results(d,j).V_1Hz = ...
            V_centered;

        Results(d,j).Speed = ...
            Speed_centered;

    end

end

%% RMSE convergance across all trajectories

valid = ~arrayfun( ...
        @(x) isempty(x.RMSE), ...
        Results);
    
allRMSE = ...
        vertcat(Results(valid).RMSE);
    
meanRMSE = ...
        mean(allRMSE,1,'omitnan');
    
allRMSE_dopp = ...
        vertcat(Results(valid).RMSE_v_Speed);
    
meanRMSE_dopp = ...
        mean(allRMSE_dopp,1,'omitnan');
    
allRMSE_dopp_smoothed = ...
        vertcat(Results(valid).RMSE_v_Speed_Smoothed);
    
meanRMSE_dopp_smoothed = ...
        mean(allRMSE_dopp_smoothed,1,'omitnan');

allRMSE_dopp_v_dopp_smoothed = ...
    vertcat(Results(valid).RMSE_Speed_v_Speed_Smoothed);

meanRMSE_dopp_v_dopp_smoothed = ...
    mean(allRMSE_dopp_v_dopp_smoothed,1,'omitnan');

Resolutions = [1 res_sec];

%% Plot the RMSEs

figure
hold on

for k = 1:numel(Results)

if isempty(Results(k).RMSE)

continue

end

% RMSE relative to max temporal resolution
scatter( ...
Results(k).Resolution_Sec, ...
Results(k).RMSE, ...
20, ...
'b', ...
'filled', ...
'HandleVisibility','off');

% RMSE relative to Doppler speed
scatter( ...
Results(k).Resolution_Sec, ...
Results(k).RMSE_v_Speed, ...
20, ...
'r', ...
'filled', ...
'HandleVisibility','off');

% RMSE relative to Smoothed Doppler speed
scatter( ...
Results(k).Resolution_Sec, ...
Results(k).RMSE_v_Speed_Smoothed, ...
20, ...
'g', ...
'filled', ...
'HandleVisibility','off');

% RMSE of 1Hz Speed relative to Smoothed Doppler speed
scatter( ...
    Results(k).Resolution_Sec, ...
    Results(k).RMSE_Speed_v_Speed_Smoothed, ...
    20, ...
    'm', ...
    'filled', ...
    'HandleVisibility','off');

end

%% Plot mean RMSE curves

h1 = scatter( ...
nan,nan, ...
20, ...
'b', ...
'filled');

h2 = scatter( ...
nan,nan, ...
20, ...
'r', ...
'filled');

h3 = scatter( ...
nan,nan, ...
20, ...
'g', ...
'filled');

h4 = scatter( ...
    nan,nan, ...
    20, ...
    'm', ...
    'filled');

h5 = plot( ...
Resolutions, ...
meanRMSE, ...
'b-', ...
'LineWidth',2);

h6 = plot( ...
Resolutions, ...
meanRMSE_dopp, ...
'r-', ...
'LineWidth',2);

h7 = plot( ...
Resolutions, ...
meanRMSE_dopp_smoothed, ...
'g-', ...
'LineWidth',2);

h8 = plot( ...
    Resolutions, ...
    meanRMSE_dopp_v_dopp_smoothed, ...
    'm-', ...
    'LineWidth',2);

xlabel('Temporal Resolution (s)')

ylabel('RMSE (m/s)')

legend( ...
[h1 h2 h3 h4 h5 h6 h7 h8], ...
'RMSE vs Max Lag', ...
'RMSE vs Doppler Speed', ...
'RMSE vs Smoothed Doppler Speed', ...
'RMSE 1Hz Dopp vs Smoothed Doppler Speed', ...
'Mean RMSE vs Max Lag', ...
'Mean RMSE vs Doppler Speed', ...
'Mean RMSE vs Smoothed Doppler Speed', ...
'Mean RMSE 1Hz Dopp vs Smoothed Doppler Speed', ...
'Location','best')

title(sprintf( ...
'RMSE Across All Trajectories - Stationary Deployment'))

grid on

%% Doppler Speed Variance Spectra


fs = 1;     % Hz - data are still sampled at 1 Hz

figure;
hold on;

% Raw 1 Hz Doppler speed
Speed_raw = Speed;

% Remove NaNs
valid = ~isnan(Speed_raw);
Speed_temp = Speed_raw(valid);

% Remove mean
Speed_temp = detrend(Speed_temp,0);

% PSD
[Pxx,f] = pwelch(Speed_temp,2^12,[],[],fs);

loglog(f,Pxx,'LineWidth',1.5,...
    'DisplayName','1 s');

% %% Smoothed Doppler speeds
% for r = 1:numel(TemporalRes)
% 
%     Speed_temp = TemporalRes(r).Speed;
% 
%     % Remove NaNs
%     valid = ~isnan(Speed_temp);
%     Speed_temp = Speed_temp(valid);
% 
%     % Remove mean
%     Speed_temp = detrend(Speed_temp,0);
% 
%     % PSD
%     [Pxx,f] = pwelch(Speed_temp,2^12,[],[],fs);
% 
%     loglog(f,Pxx,...
%         'LineWidth',1.5,...
%         'DisplayName',sprintf('%d s',TemporalRes(r).Seconds));
% 
% end

xlabel('Frequency (Hz)')
ylabel('PSD of Doppler Speed [(m/s)^2/Hz]')

legend('Location','best')
grid on
box on

title('Doppler Speed Variance Spectra')

