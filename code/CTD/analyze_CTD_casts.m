clear all
close all
% -------------------------------------------------
%% Analyze casts along navigation channel
% -------------------------------------------------
%% 0) list of input L0 files: ADCP (pressure) & CTD
rootDIR ='~/OneDriveUNCW/DATA/CB_YachtBasin/'
CTDFfiles= {fullfile(rootDIR,'JuneJamboree','PROCESSED_DATA','CTD','230442_CTDF_L0.mat');...
            fullfile(rootDIR,'SummerFlood' ,'PROCESSED_DATA','CTD','230442_CTDF_L0.mat');...
            fullfile(rootDIR,'FallFlood' ,'PROCESSED_DATA','RBR','230442_20260930_CTDF_L0.mat');...
            fullfile(rootDIR,'FallFlood' ,'PROCESSED_DATA','RBR','230442_20261028_CTDF_L0.mat')};
ADCPfiles= {fullfile(rootDIR,'JuneJamboree','PROCESSED_DATA','ADCP','JuneJamboree_NorthADCP_adcp_level_1.mat');...
            fullfile(rootDIR,'SummerFlood','PROCESSED_DATA','ADCP','SummerFlood_NorthADCP_adcp_level_1.mat');...
            fullfile(rootDIR,'FallFlood','PROCESSED_DATA','ADCP','FallFlood_NorthADCP_adcp_level_1.mat');...            
            fullfile(rootDIR,'FallFlood','PROCESSED_DATA','ADCP','FallFlood_NorthADCP_adcp_level_1.mat')};
% 
% ------------
%% 1) define ROI for casts:
% 
% 1.1) pre-allocate regular along-channel and depth grid
%
%
% -----------
%% 2) load indivual dates:
% for ii=1:length(CTDFfiles)
    ii=3;
    CTDFfile = CTDFfiles{ii};
    ADCPfile = ADCPfiles{ii};
    %
    % 2.0) load CTDF
    CTD = load(CTDFfile);
    
    % 2.1) get the time/pressure grid for ADCP
    
    %
    % 2.2) get casts in ROI
    %
    % 2.4) separate by tidal phase & interpolate to grid
%
% __________
%% 3) Plot individual casts:
%
% 3.1) Temperature (y,z)
%
% 3.2) Salinity (y,z)
%
% 3.3) Density (y,z)
%
% 3.4) NTU (y,z)
%
% 3.5) CDOM (y,z)
%
% 3.7) Chl-a (y,z)
%
end
% ---------
%% 4) Plot compendium by tidal phase:
