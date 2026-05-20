function cfg = get_obs_config(obsPeriod, adcpLoc, root_dir, root_data)

% Point to post Porcessing folder
cfg.root_dir = [root_dir, 'POST_PROCESSING_V2\'];

% Wherever you hva the ADCP data on your machine
cfg.root_data = [root_data, 'CB_YachtBasin\'];
% cfg.obsTag = strrep(obsPeriod,' ','_');
cfg.obsTag = [strrep(obsPeriod,' ','_') '_' adcpLoc];



%% Set figure output dirs
cfg.out.adcp_figures = fullfile(cfg.root_dir, 'FIGURES', 'ADCP');
cfg.out.drifters_figures = fullfile(cfg.root_dir, 'FIGURES', 'DRIFTERS');
cfg.out.comp_figures = fullfile(cfg.root_dir, 'FIGURES', 'COMP');


% Make sure output folders exist
if ~exist(cfg.out.adcp_figures, 'dir'); mkdir(cfg.out.adcp_figures); end
if ~exist(cfg.out.drifters_figures, 'dir'); mkdir(cfg.out.drifters_figures); end
if ~exist(cfg.out.comp_figures, 'dir'); mkdir(cfg.out.comp_figures); end


switch [obsPeriod, adcpLoc]
    
    case { ['MarshMadness', 'NorthADCP'] }
        
        cfg.name = 'MarshMadness';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Marsh_Madness\');
        
        cfg.root_data = [cfg.root_data, 'MarshMadness\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_032625\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_032625\Dep_times\Dep_01_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N101';
        
     case { ['MarshMadness', 'SouthADCP'] }
        
        cfg.name = 'MarshMadness';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Marsh_Madness\');
        
        cfg.root_data = [cfg.root_data, 'MarshMadness\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_032625\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_032625\Dep_times\Dep_01_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD13862\'];
        cfg.adcp_input_file = 'CB_S202';
        
        
    case { ['MarshMayhem', 'NorthADCP'] }
        
        cfg.name = 'MarshMayhem';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\MarshMayhem\');
        
        cfg.root_data = [cfg.root_data, 'MarshMayhem\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_052725\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_052725\Dep_times\Dep_02_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N201';
        
     case { ['MarshMayhem', 'SouthADCP'] }
        
        cfg.name = 'MarshMayhem';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\MarshMayhem\');
        
        cfg.root_data = [cfg.root_data, 'MarshMayhem\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_052725\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_052725\Dep_times\Dep_02_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD13862\'];
        cfg.adcp_input_file = 'CB_S202';
        
    case { ['FallFrolic 100825', 'NorthADCP'] }
        
        cfg.name = 'FallFrolic';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\Dep_times\Dep_03_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N301';
        
    case { ['FallFrolic 100825', 'SouthADCP'] }
        
        cfg.name = 'FallFrolic';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\Dep_times\Dep_03_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S301';
        
    case { ['FallFrolic 101025', 'NorthADCP'] }
        
        cfg.name = 'FallFrolic';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\Dep_times\Dep_04_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N301';
        
     case { ['FallFrolic 101025', 'SouthADCP'] }
        
        cfg.name = 'FallFrolic';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\Dep_times\Dep_04_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S301';
        
    case { ['NovDep', 'NorthADCP'] }
        
        cfg.name = 'NovDep';
        
        cfg.root_data = [cfg.root_data, 'NovDep\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_110725\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_110725\Dep_times\Dep_05_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N401_01';        
        

        
    
    
end


 %% Set data output dirs
        cfg.out.adcp_data = fullfile(cfg.root_data, 'PROCESSED_DATA', 'ADCP');
        cfg.out.drifters_data = fullfile(cfg.root_data, 'PROCESSED_DATA', 'DRIFTERS');
        cfg.out.comp_data = fullfile(cfg.root_data, 'PROCESSED_DATA', 'COMP');
        %cfg.obsTag = strrep(obsPeriod,' ','_');

        % Make sure output folders exist
        if ~exist(cfg.out.adcp_data, 'dir'); mkdir(cfg.out.adcp_data); end
        if ~exist(cfg.out.drifters_data, 'dir'); mkdir(cfg.out.drifters_data); end
        if ~exist(cfg.out.comp_data, 'dir'); mkdir(cfg.out.comp_data); end


end