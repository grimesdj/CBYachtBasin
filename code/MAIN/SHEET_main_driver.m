clear all; 
close all;
%% *SHEET 5-MIN RES DATA* Main driver script

% Wherever you have the data folder on your local machine
root_data = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'; %This is the path to the CB YACHT BASIN shared folder on your machine
% Wherever you have the post processing GIT repo
root_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\CBYachtBasin\'; % Path to the Post Processing folder on your machine

% root_data = 'C:\Users\Ben\OneDrive - UNC-Wilmington\'; %This is the path to the CB YACHT BASIN shared folder on your machine
% root_dir = 'C:\Users\Ben\OneDrive - UNC-Wilmington\THESIS\CBYachtBasin\'; % Path to the Post Processing folder on your machine
%% Choose the observation period
obsPeriod = 'FallFrolic 101025 SHEET';
% Options:
% 'MarshMadness SHEET' *South ADCP was HR (ONLY SAMPLED HALF WATER COLUMN*
% 'MarshMayhem SHEET' *South ADCP was HR (ONLY SAMPLED HALF WATER COLUMN*
% 'FallFrolic 100825 SHEET'
% 'FallFrolic 101025 SHEET'
% 'NovDep SHEET' *ONLY HAD NORTH ADCP*
% 'JuneJamboree SHEET'
% 'SummerFlood SHEET'

%% Choose North or South ADCP for analysis
adcpLoc = 'NorthADCP';
% Options:
% 'NorthADCP'
% 'SouthADCP'

% Get all the necessary configuration for that observation period
cfg = get_obs_config(obsPeriod, adcpLoc, root_dir, root_data);

disp(cfg);

addpath(genpath(fullfile(cfg.root_dir, 'FUNCTIONS')));
%% Run ADCP Level 1

%% Run Drifter level 1

dep_fin = readtable(cfg.drifters_dep_times);
drifters_root_dir = cfg.drifters_raw_dir;

[drifters_level_1,raw_drifters] = SHEET_drifter_post_processing(dep_fin,drifters_root_dir, cfg);

save(fullfile(cfg.out.drifters_data, [cfg.name '_drifters_level_1.mat']), 'drifters_level_1');
save(fullfile(cfg.out.drifters_data, [cfg.name '_raw_drifters.mat']), 'raw_drifters');