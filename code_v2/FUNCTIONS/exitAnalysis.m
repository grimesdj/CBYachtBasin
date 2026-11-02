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