clear all
close all

%% 1) list of deployment directories/files
ADCPlist = {'/Users/derekgrimes/OneDriveUNCW/Data/CB_YachtBasin/MarshMadness/AQD5459/CB_N101_L0.mat';...
            '/Users/derekgrimes/OneDriveUNCW/Data/CB_YachtBasin/MarshMayhem/AQD5459/CB_N201_L0.mat';...
            '/Users/derekgrimes/OneDriveUNCW/Data/CB_YachtBasin/FallFrolic/AQD5459/CB_N301_L0.mat';...
            '/Users/derekgrimes/OneDriveUNCW/Data/CB_YachtBasin/NovDep/AQD5459/CB_N401_02_L0.mat'};

waterlevel_time_lags = [];
velocity_time_lags = [];
for jj=1:length(ADCPlist)
    adcpFile = ADCPlist{jj};
    %% 2) modify estimate_... and plot_... scripts to take above inputs
    tmp_wl = estimate_WrightsvilleBeach_tide_phase_lag(adcpFile);

    [tmp_da, tmp_sa] = estimate_CBYachtBasin_WaterLevel_Velocity_PhaseLag(adcpFile)
    %% 3) compile phase lags, etc., into a table
    waterlevel_time_lags = cat(1,waterlevel_time_lags,tmp_wl); 
    velocity_time_lags = cat(1,velocity_time_lags,[tmp_da, tmp_sa]); 
end

table(waterlevel_time_lags(:,1),waterlevel_time_lags(:,2),waterlevel_time_lags(:,3),'VariableNames',{'Max X-Corr Lag', 'Avg. High Tide Lag', 'Avg. Low Tide Lag'})

table(velocity_time_lags(:,1),velocity_time_lags(:,2),'VariableNames',{'Depth Avg. Lag', 'Surface Avg. Lag'})

