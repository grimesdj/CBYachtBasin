function [TidalData] = TideandRegionAnalysis(cfg)
%%
obsFiles = {'MarshMadness\PROCESSED_DATA\DRIFTERS\MarshMadness_TidalPhase.mat'
    'MarshMayhem\PROCESSED_DATA\DRIFTERS\MarshMayhem_TidalPhase.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_100825_TidalPhase.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_101025_TidalPhase.mat'
    'NovDep\PROCESSED_DATA\DRIFTERS\NovDep_TidalPhase.mat'
    'JuneJamboree\PROCESSED_DATA\DRIFTERS\JuneJamboree_TidalPhase.mat'};

names = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree'};

%% Choose temporal differencing window

temporal_res = 8;   % seconds
% 1 = use original 1 Hz velocities
% Other options: 2, 4, 8, 16, 32, 64, 128, 256
%% Load in the data fromm each obs period
for i = 1:numel(obsFiles)

    filename = fullfile(cfg.root_data, obsFiles{i});

    S = load(filename);

    data.(names{i}) = S.TidalPhase;

end
%% Now organize by phase
phases = {'High', 'Low', 'Rising', 'Falling'};

for p = 1:numel(phases)
    
    phase = phases{p};
    
    for i = 1:numel(names)
        
        obs = names{i};
        if ~isempty(data.(obs).(phase))
            
            TidalData.(phase).(obs) = data.(obs).(phase);
  
        end
    end
end

%% Recalculate velocity using selected temporal differencing window
% 
% if temporal_res > 1
% 
%     proj = projcrs(32119);
% 
%     halfLag = temporal_res/2;
% 
%     for p = 1:numel(phases)
% 
%         phase = phases{p};
% 
%         obs_names = fieldnames(TidalData.(phase));
% 
%         for i = 1:numel(obs_names)
% 
%             obs = obs_names{i};
% 
%             for j = 1:numel(TidalData.(phase).(obs))
% 
%                 D = TidalData.(phase).(obs)(j);
% 
%                 %% Skip empty trajectories
% 
%                 if isempty(D.Lat) || isempty(D.Lon)
%                     continue
%                 end
% 
%                 N = length(D.Lat);
% 
%                 %% Make sure trajectory is long enough
% 
%                 if N <= temporal_res
%                     continue
%                 end
% 
% 
%                 %% Convert Lat/Lon to Easting/Northing
% 
%                 [Easting,Northing] = projfwd( ...
%                     proj, ...
%                     D.Lat, ...
%                     D.Lon);
% 
% 
%                 %% Centered differencing indices
% 
%                 idx_center = (1 + halfLag):(N - halfLag);
% 
%                 idx_before = idx_center - halfLag;
% 
%                 idx_after = idx_center + halfLag;
% 
% 
%                 %% Position differences
% 
%                 delta_east = ...
%                     Easting(idx_after) - Easting(idx_before);
% 
%                 delta_north = ...
%                     Northing(idx_after) - Northing(idx_before);
% 
% 
%                 %% Calculate new velocities
% 
%                 new_v_e = delta_east ./ temporal_res;
% 
%                 new_v_n = delta_north ./ temporal_res;
% 
%                 new_V = sqrt(new_v_e.^2 + new_v_n.^2);
% 
% 
%                 %% Replace velocity fields
% 
%                 TidalData.(phase).(obs)(j).v_e = new_v_e;
% 
%                 TidalData.(phase).(obs)(j).v_n = new_v_n;
% 
%                 TidalData.(phase).(obs)(j).V = new_V;
% 
% 
%                 %% Center corresponding data
% 
%                 TidalData.(phase).(obs)(j).Lat = ...
%                     D.Lat(idx_center);
% 
%                 TidalData.(phase).(obs)(j).Lon = ...
%                     D.Lon(idx_center);
% 
%                 TidalData.(phase).(obs)(j).Time = ...
%                     D.Time(idx_center);
% 
%                 TidalData.(phase).(obs)(j).Regions = ...
%                     D.Regions(idx_center);
% 
%             end
%         end
%     end
% 
% end

%% Histograms

for p = 1:numel(phases)

    phase = phases{p};
    obs_names = fieldnames(TidalData.(phase));

    V_all = [];

    for i = 1:numel(obs_names)

        obs = obs_names{i};

        V = vertcat(TidalData.(phase).(obs).V);

        V_all = [V_all; V];

    end

    Velocity.(phase) = V_all;

end

edges = 0:0.05:2.5;

for p = 1:numel(phases)

    phase = phases{p};

    figure

    histogram(Velocity.(phase), edges)

    xlabel('Velocity (m/s)')
    ylabel('Count')
    title([phase ' Tide'])

    xlim([0 2.5])
    grid on

end

%% Scatter trajectories with velocities as color

 for p = 1:numel(phases)

    phase = phases{p};

    figure
    gx = geoaxes;
    hold(gx,'on')

    % Observation periods that actually exist for this phase
    obs_names = fieldnames(TidalData.(phase));

    for i = 1:numel(obs_names)

        obs = obs_names{i};

        % Loop through individual drifters
        for j = 1:numel(TidalData.(phase).(obs))

            Lat = TidalData.(phase).(obs)(j).Lat;
            Lon = TidalData.(phase).(obs)(j).Lon;
            V   = TidalData.(phase).(obs)(j).V;

            if isempty(Lat) || isempty(Lon) || isempty(V)
                continue
            end

            geoscatter(gx, Lat, Lon, 15, V, 'filled')

        end
    end

    geobasemap(gx,'satellite')
    gx.CLim = [0 2];

    cb = colorbar;
    cb.Label.String = 'Velocity (m/s)';

    title([phase ' Tide Drifter Velocities'])

end


% TidalData = TidalData;

end

%% QAQC
% V_thresh = 2.5;
% 
% phases = fieldnames(TidalData);
% 
% for p = 1:numel(phases)
% 
%     phase = phases{p};
% 
%     obs_names = fieldnames(TidalData.(phase));
% 
%     for i = 1:numel(obs_names)
% 
%         obs = obs_names{i};
% 
%         for j = 1:numel(TidalData.(phase).(obs))
% 
%             D = TidalData.(phase).(obs)(j);
% 
%             idx = find(D.V > V_thresh);
% 
%             if ~isempty(idx)
% 
%                 DrifterID = repmat(D.ID, numel(idx), 1);
% 
%                 T = table( ...
%                     DrifterID, ...
%                     idx, ...
%                     D.Time(idx), ...
%                     D.Lat(idx), ...
%                     D.Lon(idx), ...
%                     D.V(idx), ...
%                     'VariableNames', ...
%                     {'DrifterID','Index','Time','Lat','Lon','V'});
% 
%                 fprintf('\n%s | %s\n', phase, obs)
%                 disp(T)
% 
%             end
%         end
%     end
% end
% 
% figure;geoscatter(TidalData.Rising.MarshMadness(11).Lat,TidalData.Rising.MarshMadness(11).Lon,10,TidalData.Rising.MarshMadness(11).V,'filled');geobasemap('satellite');
% idx = find(c.Position(1) == TidalData.Rising.MarshMadness(11).Lat & c.Position(2) == TidalData.Rising.MarshMadness(11).Lon)
% idx = idx(1);
% TidalData.Rising.MarshMadness(11).Time(idx)