%% run_all_periods.m
% Re-run selected processing levels for every period in allObsPeriods.m,
% then make the all-period plots. Use this for quick sensitivity tests
% instead of running main_driver once per period.
%
% Which levels to re-run depends on the setting you changed:
%   adcp_min_valid_frac, roi_radius, drifter_depth,
%   max_time_offset_min, adcp_latlon            -> Level 3 only
%   p_sensor_hab, rho, extrap_method,
%   extrap_nfit, eof_nmodes                     -> Levels 2 and 3
% Level 4 (ROI profile figures) is only needed for the profile plots.

clear; close all

%% Paths (same as main_driver)
root_data = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\';
root_dir  = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\CBYachtBasin\';

adcpLoc = 'NorthADCP';

%% What to run
runLevel2 = false;      % adcp_post_processing
runLevel3 = true;      % comparison
runLevel4 = true;     % ROI profile figures
runAll    = true;      % all-period plots at the end
allPlotsFcn = 'comparisonAllPlots';      % all-period plot function

%% Temporary setting overrides (applied on top of get_obs_config)
% Leave as struct() to use exactly what get_obs_config sets.
overrides = struct();
% overrides.adcp_min_valid_frac = eps;   % test: top bin valid anywhere in transit

%% Run
P = allObsPeriods();

cfg0 = get_obs_config('ALL', adcpLoc, root_dir, root_data);
addpath(genpath(fullfile(cfg0.root_dir, 'FUNCTIONS')));

oldVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');     % faster; figures are still saved
cleanupVis = onCleanup(@() set(0, 'DefaultFigureVisible', oldVis));

for n = 1:numel(P)

    obsPeriod = strrep(P(n).name, '_', ' ');   % 'FallFrolic_100825' -> 'FallFrolic 100825'
    fprintf('\n=============== %s ===============\n', P(n).obsTag);

    cfg = get_obs_config(obsPeriod, adcpLoc, root_dir, root_data);
    ovr = fieldnames(overrides);
    for k = 1:numel(ovr)
        cfg.(ovr{k}) = overrides.(ovr{k});
        fprintf('  override: cfg.%s = %s\n', ovr{k}, mat2str(overrides.(ovr{k})));
    end

    if runLevel2
        adcp_level_2 = adcp_post_processing(cfg);
        save(fullfile(cfg.out.adcp_data, [cfg.obsTag '_adcp_level_2.mat']), 'adcp_level_2');
    end

    if runLevel3
        [ROI_ADCP, ROI_EXTRAP, ROI_EOF, COMP] = comparison_v100(cfg);
        save(fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_ADCP.mat']),   'ROI_ADCP');
        save(fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_EXTRAP.mat']), 'ROI_EXTRAP');
        save(fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_EOF.mat']),    'ROI_EOF');
        save(fullfile(cfg.out.comp_data, [cfg.obsTag '_COMP.mat']),       'COMP');
    end

    if runLevel4
        ROI(cfg);
    end

    close all
end

if runAll
    cfgAll = cfg0;
    ovr = fieldnames(overrides);
    for k = 1:numel(ovr)
        cfgAll.(ovr{k}) = overrides.(ovr{k});
    end
    ALL = feval(allPlotsFcn, cfgAll);
end
