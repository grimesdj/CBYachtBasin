function cfg = get_obs_config(obsPeriod, adcpLoc, root_dir, root_data)

% Point to post Porcessing folder
cfg.root_dir = [root_dir, 'code_v2\'];

% Wherever you hva the ADCP data on your machine
cfg.root_data = [root_data, 'CB_YachtBasin\'];
% cfg.obsTag = strrep(obsPeriod,' ','_');
cfg.obsTag = [strrep(obsPeriod,' ','_') '_' adcpLoc];

% Processing flags
cfg.is_sheet = contains(obsPeriod,'SHEET');          % true only for the SHEET observation periods
cfg.run_drifter_level_1 = ~cfg.is_sheet;             % the normal level 1 (and veloAnalysis) are skipped for SHEET periods
cfg.sheet_data_dir = fullfile(root_data, 'CB_YachtBasin', 'Archived_Sheet_Data', 'RAW');   % folder with data_01.xlsx, data_02.xlsx, ... (5 OR 2-minute sheet data)
cfg.sheet_file_suffix = '';     % '' = read data_XX.xlsx (v1, the FallFrolic sheet periods); later sheet periods set '_v2' in their case




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
        
        cfg.name = 'FallFrolic_100825';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\Dep_times\Dep_03_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N301';
        
    case { ['FallFrolic 100825', 'SouthADCP'] }
        
        cfg.name = 'FallFrolic_100825';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\Dep_times\Dep_03_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S301';
        
        case { ['FallFrolic 100825 SHEET', 'NorthADCP'] }
        
        cfg.name = 'FallFrolic_100825_sheet';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'Archived_Sheet_Data\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\Dep_times\Dep_03_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N301';
        
    case { ['FallFrolic 100825 SHEET', 'SouthADCP'] }
        
        cfg.name = 'FallFrolic_100825_sheet';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'Archived_Sheet_Data\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_100825\Dep_times\Dep_03_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S301';
        
    case { ['FallFrolic 101025', 'NorthADCP'] }
        
        cfg.name = 'FallFrolic_101025';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\Dep_times\Dep_04_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N301';
        
     case { ['FallFrolic 101025', 'SouthADCP'] }
        
        cfg.name = 'FallFrolic_101025';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\Dep_times\Dep_04_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S301';
        
        case { ['FallFrolic 101025 SHEET', 'NorthADCP'] }
        
        cfg.name = 'FallFrolic_101025_sheet';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'Archived_Sheet_Data\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_101025\Dep_times\Dep_04_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N301';
        
     case { ['FallFrolic 101025 SHEET', 'SouthADCP'] }
        
        cfg.name = 'FallFrolic_101025_sheet';
        % cfg.main_dir = ('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\Fall_Frolic\');
        
        cfg.root_data = [cfg.root_data, 'FallFrolic\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'Archived_Sheet_Data\RAW'];
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
        
    case { ['JuneJamboree', 'NorthADCP'] }
        
        cfg.name = 'JuneJamboree';
        
        cfg.root_data = [cfg.root_data, 'JuneJamboree\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\RAW'];
%        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\Dep_times\Dep_06_times.xlsx'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\Dep_times\Dep_06_times_01.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N502'; 
        
     case { ['JuneJamboree', 'SouthADCP'] }
        
        cfg.name = 'JuneJamboree';
        
        cfg.root_data = [cfg.root_data, 'JuneJamboree\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\RAW'];
%        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\Dep_times\Dep_06_times.xlsx'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\Dep_times\Dep_06_times_01.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S501'; 

     case { ['SummerFlood', 'NorthADCP'] }
        
        cfg.name = 'SummerFlood';
        
        cfg.root_data = [cfg.root_data, 'SummerFlood\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_081126\RAW'];

        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_081126\Dep_times\Dep_07_times_01.xlsx'];        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N601'; 
        
     case { ['SummerFlood', 'SouthADCP'] }
        
        cfg.name = 'SummerFlood';
        
        cfg.root_data = [cfg.root_data, 'SummerFlood\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_081125\RAW'];
%        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\Dep_times\Dep_06_times.xlsx'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_061426\Dep_times\Dep_07_times_01.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S601'; 
        
        
     case { ['FallFlood 093026', 'NorthADCP'] }
        
        cfg.name = 'FallFlood_093026';
        
        cfg.root_data = [cfg.root_data, 'FallFlood\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_093026\RAW'];

        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_093026\Dep_times\Dep_08_times_v2.xlsx'];        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N701'; 
        
     case { ['FallFlood 093026', 'SouthADCP'] }
        
        cfg.name = 'FallFlood_093026';
        
        cfg.root_data = [cfg.root_data, 'FallFlood\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_093026\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_093026\Dep_times\Dep_08_times_v2.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S701'; 
        
     
    case { ['FallFlood 102826', 'NorthADCP'] }
        
        cfg.name = 'FallFlood_102826';
        
        cfg.root_data = [cfg.root_data, 'FallFlood\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_102826\RAW'];

        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_102826\Dep_times\Dep_09_times.xlsx'];        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N701'; 
        
     case { ['FallFlood 102826', 'SouthADCP'] }
        
        cfg.name = 'FallFlood_102826';
        
        cfg.root_data = [cfg.root_data, 'FallFlood\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Deployment_102826\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_102826\Dep_times\Dep_09_times.xlsx'];
        
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S701'; 
        
     
    case { ['FallFlood 102826 SHEET', 'NorthADCP'] }
        
        cfg.name = 'FallFlood_102826_sheet';
        cfg.sheet_file_suffix = '_v2';      % later sheet periods read the data_XX_v2.xlsx files
        
        cfg.root_data = [cfg.root_data, 'FallFlood\'];
        
        cfg.drifters_raw_dir = cfg.sheet_data_dir;    % not read on the sheet route; sheet files come from cfg.sheet_data_dir
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_102826\Dep_times\Dep_09_times.xlsx'];
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD5459\'];
        cfg.adcp_input_file = 'CB_N701';
        
     case { ['FallFlood 102826 SHEET', 'SouthADCP'] }
        
        cfg.name = 'FallFlood_102826_sheet';
        cfg.sheet_file_suffix = '_v2';
        
        cfg.root_data = [cfg.root_data, 'FallFlood\'];
        
        cfg.drifters_raw_dir = cfg.sheet_data_dir;
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Deployment_102826\Dep_times\Dep_09_times.xlsx'];
        
        cfg.adcp_input_dir = [cfg.root_data, 'AQD2710\'];
        cfg.adcp_input_file = 'CB_S701';
        
        
        case { ['ALL', 'NorthADCP'] }
        
        cfg.name = 'ALL_OBS_PERIODS';
        
        cfg.root_data = [cfg.root_data];
        
        
        case { ['ALL', 'SouthADCP'] }
        
        cfg.name = 'ALL_OBS_PERIODS';
        
        cfg.root_data = [cfg.root_data];
        
        
        case { ['StationaryDep', 'NorthADCP'] }
        
        cfg.name = 'StationaryDep';
        
        cfg.root_data = [cfg.root_data 'StationaryDep\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Stationary_Deployment_093026\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Stationary_Deployment_093026\Dep_times\Stationary_Dep_times.xlsx'];
        
        
        case { ['StationaryDep', 'SouthADCP'] }
        
        cfg.name = 'StationaryDep';
        
        cfg.root_data = [cfg.root_data 'StationaryDep\'];
        
        cfg.drifters_raw_dir = [cfg.root_data ,'DRIFTERS\DATA\Stationary_Deployment_093026\RAW'];
        cfg.drifters_dep_times = [cfg.root_data ,'DRIFTERS\DATA\Stationary_Deployment_093026\Dep_times\Stationary_Dep_times.xlsx'];
        
        
        
        

        
    
    
end

%% Processing settings (read by adcp_post_processing and comparison)

% Height of the ADCP pressure sensor above the bed (m), per instrument.
% Measure from the frame/mount: bed to the sensor on the Aquadopp head.
switch adcpLoc
    case 'NorthADCP'
        cfg.p_sensor_hab = 0.20;   % <-- AQD5459: replace with measured value
    case 'SouthADCP'
        cfg.p_sensor_hab = 0.00;   % <-- AQD2710 / AQD13862: replace
end

% Level 2 surface extrapolation and EOF
cfg.extrap_method = 'linear';     % 'linear' or 'constant'
cfg.extrap_nfit   = 3;            % top valid bins used for the surface shear
cfg.eof_nmodes    = 3;            % modes in the EOF reconstruction

% Level 3 comparison
cfg.adcp_latlon         = [34.05673145 -77.88949171];   % North ADCP
cfg.roi_radius          = 50;     % m
cfg.drifter_depth       = 0.5;    % m below surface
cfg.max_time_offset_min = 15;     % min
cfg.adcp_min_valid_frac = 0.5;    % raw top bin valid in this fraction of transit
cfg.adcp_window_min = 20;   % total ADCP averaging window, centred on each transit (min)


 %% Set data output dirs
        cfg.out.adcp_data = fullfile(cfg.root_data, 'PROCESSED_DATA_V2', 'ADCP');
        cfg.out.drifters_data = fullfile(cfg.root_data, 'PROCESSED_DATA_V2', 'DRIFTERS');
        cfg.out.comp_data = fullfile(cfg.root_data, 'PROCESSED_DATA_V2', 'COMP');
        %cfg.obsTag = strrep(obsPeriod,' ','_');

        % Make sure output folders exist
        if ~exist(cfg.out.adcp_data, 'dir'); mkdir(cfg.out.adcp_data); end
        if ~exist(cfg.out.drifters_data, 'dir'); mkdir(cfg.out.drifters_data); end
        if ~exist(cfg.out.comp_data, 'dir'); mkdir(cfg.out.comp_data); end
        
 %% Set figure output dirs
        cfg.out.adcp_figures = fullfile(cfg.root_data, 'PROCESSED_DATA_V2', 'ADCP','FIGURES');
        cfg.out.drifters_figures = fullfile(cfg.root_data, 'PROCESSED_DATA_V2', 'DRIFTERS','FIGURES');
        cfg.out.comp_figures = fullfile(cfg.root_data,'PROCESSED_DATA_V2', 'COMP', 'FIGURES');


 % Make sure output folders exist
        if ~exist(cfg.out.adcp_figures, 'dir'); mkdir(cfg.out.adcp_figures); end
        if ~exist(cfg.out.drifters_figures, 'dir'); mkdir(cfg.out.drifters_figures); end
        if ~exist(cfg.out.comp_figures, 'dir'); mkdir(cfg.out.comp_figures); end


end