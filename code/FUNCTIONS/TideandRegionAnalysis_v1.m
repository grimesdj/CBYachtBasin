function [TidalData] = TideandRegionAnalysis_v1(cfg)

%% Settings

temporal_res = 8;

obsFiles = { ...
    'MarshMadness\PROCESSED_DATA\DRIFTERS\MarshMadness_TemporalRes_RMSE.mat'
    'MarshMayhem\PROCESSED_DATA\DRIFTERS\MarshMayhem_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_100825_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_101025_TemporalRes_RMSE.mat'
    'NovDep\PROCESSED_DATA\DRIFTERS\NovDep_TemporalRes_RMSE.mat'
    'JuneJamboree\PROCESSED_DATA\DRIFTERS\JuneJamboree_TemporalRes_RMSE.mat'};

names = { ...
    'MarshMadness', ...
    'MarshMayhem', ...
    'FallFrolic_100825', ...
    'FallFrolic_101025', ...
    'NovDep', ...
    'JuneJamboree'};

phases = {'High','Low','Rising','Falling'};

save_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES';

%% Initialize structures

TidalData = struct();
RegionProb = struct();
TrajectoryTransitions = struct([]);
traj_counter = 0;

for p = 1:numel(phases)
    phase = phases{p};
    RegionProb.(phase).StartRegion = [];
    RegionProb.(phase).NextRegion = [];
end

%% Load and process Results

for i = 1:numel(obsFiles)

    filename = fullfile(cfg.root_data,obsFiles{i});
    S = load(filename);
    Results = S.Results;
    obs = names{i};

    for d = 1:size(Results,1)

        for j = 1:size(Results,2)

            if isempty(Results(d,j).ID)
                continue
            end

            %% Region transition analysis

            regions_transition = Results(d,j).Region(:);
            tides_transition = string(Results(d,j).TidePhase(:));

            valid_transition = ~isnan(regions_transition);

            regions_transition = regions_transition(valid_transition);
            tides_transition = tides_transition(valid_transition);

            if numel(regions_transition) >= 2

                changeIdx = [1; find(diff(regions_transition) ~= 0) + 1];

                regionSequence = regions_transition(changeIdx);
                phaseSequence = tides_transition(changeIdx);

                traj_counter = traj_counter + 1;

                TrajectoryTransitions(traj_counter).Observation = obs;
                TrajectoryTransitions(traj_counter).ID = Results(d,j).ID;
                TrajectoryTransitions(traj_counter).RegionSequence = regionSequence;
                TrajectoryTransitions(traj_counter).PhaseSequence = phaseSequence;
                TrajectoryTransitions(traj_counter).StartRegion = regionSequence(1);
                TrajectoryTransitions(traj_counter).StartPhase = phaseSequence(1);
                TrajectoryTransitions(traj_counter).Exited = ...
                    any(regionSequence(2:end) == 0);

                for kk = 1:numel(regionSequence)-1

                    phase_here = char(phaseSequence(kk));

                    if isfield(RegionProb,phase_here)

                        RegionProb.(phase_here).StartRegion(end+1,1) = ...
                            regionSequence(kk);

                        RegionProb.(phase_here).NextRegion(end+1,1) = ...
                            regionSequence(kk+1);

                    end
                end
            end

            %% Select temporal resolution

            if temporal_res == 1

                Time = Results(d,j).Time;
                Lat = Results(d,j).Lat;
                Lon = Results(d,j).Lon;
                V = Results(d,j).V_1Hz;
                TidePhase = Results(d,j).TidePhase;
                Regions = Results(d,j).Region;

                v_e = [];
                v_n = [];

            else

                res_idx = find( ...
                    [Results(d,j).TemporalRes.Seconds] == temporal_res,1);

                if isempty(res_idx)
                    continue
                end

                TR = Results(d,j).TemporalRes(res_idx);

                Time = TR.Time;
                Lat = TR.Lat;
                Lon = TR.Lon;
                V = TR.V;
                v_e = TR.V_e;
                v_n = TR.V_n;
                TidePhase = TR.TidePhase;
                Regions = TR.Region;

            end

            TidePhase = string(TidePhase);

            %% Separate by tidal phase

            for p = 1:numel(phases)

                phase = phases{p};
                idx = TidePhase == phase;

                if ~any(idx)
                    continue
                end

                if ~isfield(TidalData,phase) || ...
                   ~isfield(TidalData.(phase),obs)

                    k = 1;

                else

                    k = numel(TidalData.(phase).(obs)) + 1;

                end

                TidalData.(phase).(obs)(k).ID = Results(d,j).ID;
                TidalData.(phase).(obs)(k).Lat = Lat(idx);
                TidalData.(phase).(obs)(k).Lon = Lon(idx);
                TidalData.(phase).(obs)(k).Time = Time(idx);
                TidalData.(phase).(obs)(k).Regions = Regions(idx);
                TidalData.(phase).(obs)(k).V = V(idx);

                if temporal_res > 1

                    TidalData.(phase).(obs)(k).v_e = v_e(idx);
                    TidalData.(phase).(obs)(k).v_n = v_n(idx);

                else

                    TidalData.(phase).(obs)(k).v_e = [];
                    TidalData.(phase).(obs)(k).v_n = [];

                end
            end
        end
    end
end

%% Region transition probabilities

for p = 1:numel(phases)

    phase = phases{p};

    StartRegion = RegionProb.(phase).StartRegion;
    NextRegion = RegionProb.(phase).NextRegion;

    if isempty(StartRegion)
        continue
    end

    nRegions = max([StartRegion; NextRegion]);

    counts = zeros(nRegions,nRegions+1);

    for k = 1:numel(StartRegion)

        r1 = StartRegion(k);
        r2 = NextRegion(k);

        if r1 == 0
            continue
        end

        if r2 == 0

            counts(r1,nRegions+1) = ...
                counts(r1,nRegions+1) + 1;

        else

            counts(r1,r2) = counts(r1,r2) + 1;

        end
    end

    probability = zeros(size(counts));

    for r = 1:nRegions

        N = sum(counts(r,:));

        if N > 0
            probability(r,:) = counts(r,:) ./ N;
        end
    end

    RegionProb.(phase).Counts = counts;
    RegionProb.(phase).Probability = probability;

end

%% Exit probabilities

ExitProb = struct();

for p = 1:numel(phases)

    phase = phases{p};

    phase_idx = strcmp( ...
        string({TrajectoryTransitions.StartPhase}),phase);

    start_inside = ...
        [TrajectoryTransitions.StartRegion] > 0;

    use_idx = phase_idx & start_inside;

    if ~any(use_idx)
        continue
    end

    start_regions = ...
        [TrajectoryTransitions(use_idx).StartRegion];

    exited = ...
        [TrajectoryTransitions(use_idx).Exited];

    nRegions = max(start_regions);

    N_Start = zeros(nRegions,1);
    N_Exit = zeros(nRegions,1);
    Probability = NaN(nRegions,1);

    for r = 1:nRegions

        idx = start_regions == r;

        N_Start(r) = sum(idx);
        N_Exit(r) = sum(exited(idx));

        if N_Start(r) > 0
            Probability(r) = N_Exit(r) / N_Start(r);
        end
    end

    ExitProb.(phase).N_Start = N_Start;
    ExitProb.(phase).N_Exit = N_Exit;
    ExitProb.(phase).Probability = Probability;

end

%% Add analyses to output

TidalData.RegionProb = RegionProb;
TidalData.TrajectoryTransitions = TrajectoryTransitions;
TidalData.ExitProb = ExitProb;

%% Print probabilities

for p = 1:numel(phases)

    phase = phases{p};

    if ~isfield(TidalData.RegionProb.(phase),'Probability')
        continue
    end

    fprintf('\n%s TIDE REGION TRANSITIONS\n',upper(phase));

    P = TidalData.RegionProb.(phase).Probability;
    C = TidalData.RegionProb.(phase).Counts;

    for r = 1:size(P,1)

        N = sum(C(r,:));

        if N == 0
            continue
        end

        fprintf('\nStarting Region %d (N = %d transitions)\n',r,N);

        for dest = 1:size(P,2)

            if C(r,dest) == 0
                continue
            end

            if dest == size(P,2)

                fprintf('   Region %d -> EXIT: %5.1f%% (n = %d)\n', ...
                    r,P(r,dest)*100,C(r,dest));

            else

                fprintf('   Region %d -> Region %d: %5.1f%% (n = %d)\n', ...
                    r,dest,P(r,dest)*100,C(r,dest));

            end
        end
    end

    fprintf('\n%s TIDE EXIT PROBABILITIES\n',upper(phase));

    if isfield(TidalData.ExitProb,phase)

        Pexit = TidalData.ExitProb.(phase).Probability;
        Nstart = TidalData.ExitProb.(phase).N_Start;
        Nexit = TidalData.ExitProb.(phase).N_Exit;

        for r = 1:numel(Pexit)

            if Nstart(r) == 0
                continue
            end

            fprintf('Starting Region %d: %5.1f%% exited (%d/%d trajectories)\n', ...
                r,Pexit(r)*100,Nexit(r),Nstart(r));

        end
    end
end

%% Velocity histograms

Velocity = struct();

for p = 1:numel(phases)

    phase = phases{p};

    if ~isfield(TidalData,phase)
        continue
    end

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

    if ~isfield(Velocity,phase)
        continue
    end

    figure

    histogram(Velocity.(phase),edges)

    xlabel('Velocity (m/s)')
    ylabel('Count')

    title(sprintf( ...
        '%s Tide - Temporal differencing window: %d s', ...
        phase,temporal_res))

    xlim([0 2.5])
    grid on

end

%% Velocity geoscatter

for p = 1:numel(phases)

    phase = phases{p};

    if ~isfield(TidalData,phase)
        continue
    end

    figure

    gx = geoaxes;
    hold(gx,'on')

    obs_names = fieldnames(TidalData.(phase));

    for i = 1:numel(obs_names)

        obs = obs_names{i};

        for j = 1:numel(TidalData.(phase).(obs))

            Lat = TidalData.(phase).(obs)(j).Lat;
            Lon = TidalData.(phase).(obs)(j).Lon;
            V = TidalData.(phase).(obs)(j).V;

            if isempty(Lat) || isempty(Lon) || isempty(V)
                continue
            end

            geoscatter(gx,Lat,Lon,15,V,'filled')

        end
    end

    geobasemap(gx,'satellite')
    gx.CLim = [0 0.6];

    cb = colorbar;
    cb.Label.String = 'Velocity (m/s)';

    title([phase ' Tide Drifter Velocities'])

end

%% Load regions

load('regions.mat')

%% Standalone region map

figure

gx = geoaxes;
hold(gx,'on')

for r = 1:numel(regions)

    lat = regions(r).Lat;
    lon = regions(r).Lon;

    geoplot(gx,lat,lon,'LineWidth',2)

    lat_center = mean(lat(1:end-1),'omitnan');
    lon_center = mean(lon(1:end-1),'omitnan');

    text(gx, ...
        lat_center, ...
        lon_center, ...
        sprintf('%d',r), ...
        'FontSize',14, ...
        'Color','w', ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center');

end

geobasemap(gx,'satellite')

title('Yacht Basin Regions')

%% Region transition pathway maps

min_prob = 0.05;

for p = 1:numel(phases)

    phase = phases{p};

    if ~isfield(TidalData.RegionProb.(phase),'Probability')
        continue
    end

    P = TidalData.RegionProb.(phase).Probability;

    if isempty(P)
        continue
    end

    figure

    gx = geoaxes;
    hold(gx,'on')
    geobasemap(gx,'satellite')

    nRegions = numel(regions);

    region_lat = nan(nRegions,1);
    region_lon = nan(nRegions,1);
    region_color = nan(nRegions,3);

    %% Plot regions and save their colors

    for r = 1:nRegions

        lat = regions(r).Lat;
        lon = regions(r).Lon;

        hRegion = geoplot(gx,lat,lon,'LineWidth',2);

        region_color(r,:) = hRegion.Color;

        region_lat(r) = mean(lat(1:end-1),'omitnan');
        region_lon(r) = mean(lon(1:end-1),'omitnan');

        text(gx, ...
            region_lat(r), ...
            region_lon(r), ...
            sprintf('%d',r), ...
            'FontSize',14, ...
            'Color','w', ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');

    end

    %% Plot transition pathways

    for r1 = 1:size(P,1)

        for r2 = 1:size(P,2)

            if P(r1,r2) < min_prob || r1 == r2
                continue
            end

            % Use the color of the starting region
            path_color = region_color(r1,:);

            line_width = 0.5 + 2.5*P(r1,r2);
            marker_size = 4 + 8*P(r1,r2);

            %% Exit pathway

            if r2 == nRegions + 1

                if r1 ~= 1
                    continue
                end

                lat1 = max(regions(1).Lat);

                lon1 = mean( ...
                    regions(1).Lon( ...
                    regions(1).Lat == lat1));

                lat2 = lat1 + 0.001;
                lon2 = lon1;

                geoplot(gx, ...
                    [lat1 lat2], ...
                    [lon1 lon2], ...
                    'LineWidth',line_width, ...
                    'Color',path_color);

                lat_mark = lat1 + 0.75*(lat2-lat1);
                lon_mark = lon1;

                geoplot(gx, ...
                    lat_mark, ...
                    lon_mark, ...
                    'LineStyle','none', ...
                    'Marker','^', ...
                    'MarkerSize',marker_size, ...
                    'MarkerFaceColor',path_color, ...
                    'MarkerEdgeColor',path_color);

                continue
            end

            %% Normal region-to-region pathway

            lat1 = region_lat(r1);
            lon1 = region_lon(r1);

            lat2 = region_lat(r2);
            lon2 = region_lon(r2);

            geoplot(gx, ...
                [lat1 lat2], ...
                [lon1 lon2], ...
                'LineWidth',line_width, ...
                'Color',path_color);

            lat_mark = lat1 + 0.75*(lat2-lat1);
            lon_mark = lon1 + 0.75*(lon2-lon1);

            dlat = lat2-lat1;
            dlon = lon2-lon1;

            if abs(dlon) >= abs(dlat)

                if dlon >= 0
                    mk = '>';
                else
                    mk = '<';
                end

            else

                if dlat >= 0
                    mk = '^';
                else
                    mk = 'v';
                end

            end

            geoplot(gx, ...
                lat_mark, ...
                lon_mark, ...
                'LineStyle','none', ...
                'Marker',mk, ...
                'MarkerSize',marker_size, ...
                'MarkerFaceColor',path_color, ...
                'MarkerEdgeColor',path_color);

        end
    end

    title(sprintf('%s Tide Region Transport Pathways',phase))

    exportgraphics( ...
        gcf, ...
        fullfile( ...
            save_dir, ...
            [phase '_Tide_RegionTransitionPathways.png']), ...
        'Resolution',300);

    savefig( ...
        gcf, ...
        fullfile( ...
            save_dir, ...
            [phase '_Tide_RegionTransitionPathways.fig']));

end

end