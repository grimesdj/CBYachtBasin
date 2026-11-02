%% Exiting trajectory analysis
clear all;
close all;

%% Settings

root_data = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\';
% root_data = 'C:\Users\Ben\OneDrive - UNC-Wilmington';

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

%% Find trajectories that exit the domain through Region 1

ExitTrajectories = struct();
counter = 0;

AllStartRegions = [];
AllStartPhases = strings(0,1);

for i = 1:numel(obsFiles)

    filename = fullfile(root_data,obsFiles{i});

    S = load(filename);
    Results = S.Results;

    obs = names{i};

    for d = 1:size(Results,1)

        for j = 1:size(Results,2)

            if isempty(Results(d,j).ID)
                continue
            end

            Region = Results(d,j).Region(:);

            % Skip trajectories with no region information
            if isempty(Region)
                continue
            end

            % Only include trajectories that START inside the domain
            if Region(1) <= 0
                continue
            end
            
            % Store all start regions and tidal phase for exit percentage
            AllStartRegions(end+1) = Region(1);
            AllStartPhases(end+1,1) = string(Results(d,j).TidePhase(1));

            %% Find first exit from domain

            % Transition from inside (>0) to outside (0)
            exit_idx = find( ...
                Region(1:end-1) > 0 & ...
                Region(2:end) == 0, ...
                1,'first');

            % Skip trajectories that never exit
            if isempty(exit_idx)
                continue
            end

            % Only keep trajectories exiting through northern mouth
            % (Region 1 -> Region 0)
            if Region(exit_idx) ~= 1
                continue
            end

            %% Save ENTIRE trajectory

            counter = counter + 1;

            ExitTrajectories(counter).Observation = obs;
            ExitTrajectories(counter).ID = Results(d,j).ID;
            ExitTrajectories(counter).Time = Results(d,j).Time(:);
            ExitTrajectories(counter).Lat = Results(d,j).Lat(:);
            ExitTrajectories(counter).Lon = Results(d,j).Lon(:);
            ExitTrajectories(counter).Region = Region;
            ExitTrajectories(counter).TidePhase = Results(d,j).TidePhase(:);
            ExitTrajectories(counter).Start_Region = Region(1);
            ExitTrajectories(counter).Start_TidePhase = string(Results(d,j).TidePhase(1));
            ExitTrajectories(counter).Exit_TidePhase = string(Results(d,j).TidePhase(exit_idx+1));
            ExitTrajectories(counter).TimetoExit = (Results(d,j).Time(exit_idx + 1) - Results(d,j).Time(1)) / 60;

            % Useful to save where the exit occurred
            ExitTrajectories(counter).ExitIndex = exit_idx + 1;

            % Which deployment row / drifter column this track came from
            % (needed to look up its absolute time in the level 2 file)
            ExitTrajectories(counter).Deployment = d;
            ExitTrajectories(counter).Track = j;

        end
    end
end

%% Display exiting trajectories

fprintf('\n%d trajectories exited the domain.\n\n', ...
    numel(ExitTrajectories));

for k = 1:numel(ExitTrajectories)

    fprintf('%s | Drifter %d\n', ...
        ExitTrajectories(k).Observation, ...
        ExitTrajectories(k).ID);

end


%% Calculate percent exit for each region

regions = unique(AllStartRegions);

fprintf('\nExit percentage for each region:\n');

for r = regions
    
    total_started = sum(AllStartRegions == r);
    total_exited = sum([ExitTrajectories.Start_Region] == r);
    
    pct_exited = 100 * total_exited / total_started;
    
    fprintf('Region %d: %d of %d exited = %.1f%%\n', ...
        r, total_exited, total_started, pct_exited)
end

%% Calculate percent exit for each starting tidal phase

phases = ["High","Low","Rising","Falling"];

fprintf('\nExit percentage by STARTING tidal phase:\n')

for p = 1:numel(phases)

    phase = phases(p);

    total_started = sum(AllStartPhases == phase);

    total_exited = sum( ...
        string({ExitTrajectories.Start_TidePhase}) == phase);

    pct_exited = 100 * total_exited / total_started;

    fprintf('%s: %d of %d exited = %.1f%%\n', ...
        phase, total_exited, total_started, pct_exited)

end

%% Count tidal phase at actual exit

fprintf('\nTidal phase at EXIT:\n')

exit_phases = string({ExitTrajectories.Exit_TidePhase});

for p = 1:numel(phases)

    phase = phases(p);

    n_exit = sum(exit_phases == phase);

    pct_exit = 100 * n_exit / numel(ExitTrajectories);

    fprintf('%s: %d exits = %.1f%%\n', ...
        phase, n_exit, pct_exit)

end


%% Surface current (U_mapped) at the moment of each exit
% U_mapped = North ADCP principal-axis velocity. Positive = toward the mouth
% (outflow); this was established because every exit occurred on positive U.
% Only the TOP row of U_mapped is used (surface; partly extrapolated),
% smoothed with a 7-sample (~70 min) moving mean. It is NOT a depth average.

smooth_samples = 7;

for k = 1:numel(ExitTrajectories)
    ExitTrajectories(k).Exit_Datenum = NaN;
    ExitTrajectories(k).Exit_U       = NaN;
end

for i = 1:numel(obsFiles)

    these = find(strcmp({ExitTrajectories.Observation}, names{i}));

    if isempty(these)
        continue
    end

    drifter_dir = fileparts(fullfile(root_data,obsFiles{i}));
    adcp_dir    = fullfile(fileparts(drifter_dir),'ADCP');

    A1 = load(fullfile(adcp_dir,[names{i} '_NorthADCP_adcp_level_1.mat'])).adcp_level_1;
    A2 = load(fullfile(adcp_dir,[names{i} '_NorthADCP_adcp_level_2.mat'])).adcp_level_2;
    L2 = load(fullfile(drifter_dir,[names{i} '_drifters_level_2.mat'])).drifters_level_2;

    tADCP = A1.time(:);
    Usurf = movmean(A2.U_mapped(end,:)', smooth_samples, 'omitnan');

    for k = these

        D = L2(ExitTrajectories(k).Deployment, ExitTrajectories(k).Track);

        % Seconds-since-first-fix at the exit sample (same clock as level 2 e_time)
        t_exit = ExitTrajectories(k).Time(ExitTrajectories(k).ExitIndex);
        m = find(abs(D.e_time(:) - t_exit) < 1e-6, 1);

        if isempty(m)
            warning('Could not match exit time for %s drifter %d', ...
                names{i}, ExitTrajectories(k).ID);
            continue
        end

        te = datenum(D.Date_and_time(m));

        ExitTrajectories(k).Exit_Datenum = te;
        ExitTrajectories(k).Exit_U = interp1(tADCP, Usurf, te);

    end
end

exit_U      = [ExitTrajectories.Exit_U];
exit_phases = string({ExitTrajectories.Exit_TidePhase});

fprintf('\nSurface current at each exit (positive = toward mouth):\n')

for k = 1:numel(ExitTrajectories)

    if isnan(ExitTrajectories(k).Exit_Datenum)
        tstr = 'no match';
    else
        tstr = datestr(ExitTrajectories(k).Exit_Datenum,'mm/dd HH:MM');
    end

    fprintf('%-18s | Drifter %2d | %s | %-7s | U = %+.3f m/s\n', ...
        ExitTrajectories(k).Observation, ...
        ExitTrajectories(k).ID, ...
        tstr, ...
        exit_phases(k), ...
        exit_U(k));

end

fprintf('\nSurface current at exit, by tidal phase at exit:\n')

for p = 1:numel(phases)

    sel = exit_phases == phases(p) & ~isnan(exit_U);

    if ~any(sel)
        fprintf('%-8s no exits\n', phases(p))
        continue
    end

    fprintf('%-8s n = %2d | mean %+.3f | min %+.3f | max %+.3f m/s | positive %d of %d\n', ...
        phases(p), sum(sel), mean(exit_U(sel)), min(exit_U(sel)), ...
        max(exit_U(sel)), sum(exit_U(sel) > 0), sum(sel))

end

fprintf('All exits: n = %d | mean %+.3f | min %+.3f | max %+.3f m/s | positive %d of %d\n', ...
    sum(~isnan(exit_U)), mean(exit_U,'omitnan'), min(exit_U), max(exit_U), ...
    sum(exit_U > 0), sum(~isnan(exit_U)))

% Same information as a table in the workspace
ExitU_table = table( ...
    string({ExitTrajectories.Observation})', ...
    [ExitTrajectories.ID]', ...
    exit_phases', ...
    exit_U', ...
    'VariableNames',{'Period','Drifter','ExitPhase','Exit_U'});


%% Plot each exiting trajectory in its own figure

% Plot exit trajectories colored by tidal phase

save_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES\EXITS';

phases = ["High","Low","Rising","Falling"];

for k = 1:numel(ExitTrajectories)

    lat = ExitTrajectories(k).Lat;
    lon = ExitTrajectories(k).Lon;
    phase = string(ExitTrajectories(k).TidePhase);

    figure

    gx = geoaxes;
    hold(gx,'on')

    geobasemap(gx,'satellite')

    %% Plot each tidal phase separately

    legend_handles = gobjects(numel(phases),1);

    for p = 1:numel(phases)

        idx = phase == phases(p);

        % Copy full trajectory
        lat_plot = lat;
        lon_plot = lon;

        % Break line anywhere this phase is not active
        lat_plot(~idx) = NaN;
        lon_plot(~idx) = NaN;

        legend_handles(p) = geoplot(gx, ...
            lat_plot, ...
            lon_plot, ...
            'LineWidth',4);

    end

    %% Legend

    legend(gx, ...
        legend_handles, ...
        phases, ...
        'Location','best')

    %% Title

    title(sprintf('%s | Drifter %d | \\Deltat to Exit = %.1f min', ...
        ExitTrajectories(k).Observation, ...
        ExitTrajectories(k).ID, ...
        ExitTrajectories(k).TimetoExit)); 

    % Create unique trajectory number for this Observation + Drifter ID

current_obs = string(ExitTrajectories(k).Observation);
current_id  = ExitTrajectories(k).ID;

obs_all = string({ExitTrajectories(1:k).Observation});
id_all  = [ExitTrajectories(1:k).ID];

trajectory_number = sum( ...
    obs_all == current_obs & ...
    id_all == current_id);

%% Create unique filename

file_name = sprintf('%s_Drifter_%d_Trajectory_%02d_TidalPhase', ...
    current_obs, ...
    current_id, ...
    trajectory_number);

% Clean filename
file_name = strrep(file_name,' ','_');
file_name = strrep(file_name,'/','-');
file_name = strrep(file_name,'\','-');

%% Save PNG

exportgraphics(gcf, ...
    fullfile(save_dir,[file_name '.png']), ...
    'Resolution',300);

%% Save FIG

savefig(gcf, ...
    fullfile(save_dir,[file_name '.fig']));

end