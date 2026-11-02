clear all; 
close all;
%% Main driver script

% Wherever you have the data folder on your local machine
root_data = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'; %This is the path to the CB YACHT BASIN shared folder on your machine
% Wherever you have the post processing GIT repo
root_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\CBYachtBasin\'; % Path to the Post Processing folder on your machine
% 
% root_data = 'C:\Users\Ben\OneDrive - UNC-Wilmington\'; %This is the path to the CB YACHT BASIN shared folder on your machine
% root_dir = 'C:\Users\Ben\OneDrive - UNC-Wilmington\THESIS\CBYachtBasin\'; % Path to the Post Processing folder on your machine

% DATES FOR SOPHIE MARCH 10-11 & MAY 27-28



%% Choose the observation period
obsPeriod = 'JuneJamboree';
% Options: 
% 'MarshMadness' *South ADCP was HR (ONLY SAMPLED HALF WATER COLUMN*
% 'MarshMayhem' *South ADCP was HR (ONLY SAMPLED HALF WATER COLUMN*
% 'FallFrolic 100825'
% 'FallFrolic 100825 SHEET'
% 'FallFrolic 101025'
% 'FallFrolic 101025 SHEET'
% 'NovDep' *ONLY HAD NORTH ADCP*
% 'JuneJamboree'
% 'SummerFlood" 
% 'FallFlood 093026' *ADCP data not available'
% 'FallFlood 102826' *ADCP data not available'
% 'FallFlood 102826 SHEET' *ADCP data not available'
% 'ALL' **THIS IS FOR ALL OBS-PERIOD ANALYSIS, NOT INTENDED FOR INTIAL PP**
% 'StationaryDep' ** This is for drifter position error analysis **

%% Choose North or South ADCP for analysis
adcpLoc = 'NorthADCP';
% Options:
% 'NorthADCP'
% 'SouthADCP'

% Get all the necessary configuration for that observation period
cfg = get_obs_config(obsPeriod, adcpLoc, root_dir, root_data);

disp(cfg);

addpath(genpath(fullfile(cfg.root_dir, 'FUNCTIONS')));

%% Stop here for the ALL observation periods
% 'ALL' is only for all-period analysis (functions that take cfg), not for processing.
if strcmp(obsPeriod,'ALL')
    fprintf('obsPeriod = ALL: cfg is ready, stopping before North ADCP level 1.\n');
    return
end

%% Run ADCP level 1

if ~cfg.is_sheet
    
    [adcp_level_1] = load_and_save_adcp(cfg);

    save(fullfile(cfg.out.adcp_data, [cfg.obsTag '_adcp_level_1.mat']), 'adcp_level_1');

end


%% Run Drifter level 1

if cfg.is_sheet

    % SHEET period: rebuild the tracks from the 5-minute data posted to the Google Sheet
    [drifters_level_1] = drifter_post_processing_sheet(cfg);

    save(fullfile([cfg.out.drifters_data], ...
        [cfg.name '_drifters_level_1.mat']), ...
        'drifters_level_1');
%     save(fullfile(cfg.out.drifters_data, ...
%         [cfg.name '_raw_drifters.mat']), 'raw_drifters');

elseif cfg.run_drifter_level_1

    [drifters_level_1,raw_drifters] = drifter_post_processing(cfg);

    save(fullfile(cfg.out.drifters_data, ...
        [cfg.name '_drifters_level_1.mat']), ...
        'drifters_level_1');
    save(fullfile(cfg.out.drifters_data, ...
        [cfg.name '_raw_drifters.mat']), 'raw_drifters');

else

    fprintf('Skipping Drifter Level 1 for %s\n', cfg.name);

end


%% Run Drifter level 2: Temporal Res and Assigning Tidal phase and region

[TidalPhase, drifters_level_2] = findTidalPhaseandRegion(cfg, adcpLoc);
[drifters_level_2, StuckInfo] = flagStuckDrifters(cfg, drifters_level_2);
save(fullfile(cfg.out.drifters_data, [cfg.name '_drifters_level_2.mat']), 'drifters_level_2');
save(fullfile(cfg.out.drifters_data, [cfg.name '_TidalPhase.mat']), 'TidalPhase');


if cfg.run_drifter_level_1

    % Results = veloAnalysis(cfg);
    Results = veloAnalysis_v2(cfg);

    save(fullfile(cfg.out.drifters_data, [cfg.name '_TemporalRes_RMSE.mat']), 'Results');

%     save(fullfile(cfg.out.drifters_data, ...
%         [cfg.name '_raw_drifters.mat']), 'raw_drifters');

else

    fprintf('Skipping VELO ANALYSIS (CREATING DIFF TEMPORAL DIFFERENCING WINDOWS) for %s\n', cfg.name);

end




%% Run ADCP level 2

if ~cfg.is_sheet

    % [adcp_level_2] = adcp_post_processing(cfg);
    [adcp_level_2] = adcp_post_processing(cfg);

    save(fullfile(cfg.out.adcp_data, [cfg.obsTag '_adcp_level_2.mat']), 'adcp_level_2');

end

%% Run Level 3: Drifter and ADCP comparison

[ROI_ADCP, ROI_EXTRAP, ROI_EOF, COMP] = comparison(cfg, adcpLoc);
% [ROI_ADCP, ROI_EXTRAP, ROI_EOF, COMP] = comparison_old(cfg)

save(fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_ADCP.mat']), 'ROI_ADCP');
save(fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_EXTRAP.mat']), 'ROI_EXTRAP');
save(fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_EOF.mat']), 'ROI_EOF');
save(fullfile(cfg.out.comp_data, [cfg.obsTag '_COMP.mat']), 'COMP');

%% Run Level 4: Range of Instrument analysis

ROI(cfg, adcpLoc);


