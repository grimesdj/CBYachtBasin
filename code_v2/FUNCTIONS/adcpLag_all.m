clear all
close all

%% Set dir

main_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\CB_YachtBasin\';

save_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES\LAG';

%% Observation periods

obsNames = { ...
    'MarshMadness'
    'MarshMayhem'
    'FallFrolic_100825'
    'NovDep'
    'JuneJamboree'
    'SummerFlood'};

%% ADCP files

inputFiles1 = { ...
    'MarshMadness\PROCESSED_DATA_V2\ADCP\MarshMadness_NORTHADCP_adcp_level_1.mat'
    'MarshMayhem\PROCESSED_DATA_V2\ADCP\MarshMayhem_NORTHADCP_adcp_level_1.mat'
    'FallFrolic\PROCESSED_DATA_V2\ADCP\FallFrolic_100825_NORTHADCP_adcp_level_1.mat'
    'NovDep\PROCESSED_DATA_V2\ADCP\NovDep_NORTHADCP_adcp_level_1.mat'
    'JuneJamboree\PROCESSED_DATA_V2\ADCP\JuneJamboree_NORTHADCP_adcp_level_1.mat'
    'SummerFlood\PROCESSED_DATA_V2\ADCP\SummerFlood_NORTHADCP_adcp_level_1.mat'};

inputFiles2 = { ...
    'MarshMadness\PROCESSED_DATA_V2\ADCP\MarshMadness_NORTHADCP_adcp_level_2.mat'
    'MarshMayhem\PROCESSED_DATA_V2\ADCP\MarshMayhem_NORTHADCP_adcp_level_2.mat'
    'FallFrolic\PROCESSED_DATA_V2\ADCP\FallFrolic_100825_NORTHADCP_adcp_level_2.mat'
    'NovDep\PROCESSED_DATA_V2\ADCP\NovDep_NORTHADCP_adcp_level_2.mat'
    'JuneJamboree\PROCESSED_DATA_V2\ADCP\JuneJamboree_NORTHADCP_adcp_level_2.mat'
    'SummerFlood\PROCESSED_DATA_V2\ADCP\SummerFlood_NORTHADCP_adcp_level_2.mat'};
% 
% inputFiles1 = { ...
%     'MarshMadness\PROCESSED_DATA\ADCP\MarshMadness_NORTHADCP_adcp_level_1.mat'
%     'MarshMayhem\PROCESSED_DATA\ADCP\MarshMayhem_NORTHADCP_adcp_level_1.mat'
%     'FallFrolic\PROCESSED_DATA\ADCP\FallFrolic_100825_NORTHADCP_adcp_level_1.mat'
%     'NovDep\PROCESSED_DATA\ADCP\NovDep_NORTHADCP_adcp_level_1.mat'
%     'JuneJamboree\PROCESSED_DATA\ADCP\JuneJamboree_NORTHADCP_adcp_level_1.mat'
%     'SummerFlood\PROCESSED_DATA\ADCP\SummerFlood_NORTHADCP_adcp_level_1.mat'};
% 
% inputFiles2 = { ...
%     'MarshMadness\PROCESSED_DATA\ADCP\MarshMadness_NORTHADCP_adcp_level_2.mat'
%     'MarshMayhem\PROCESSED_DATA\ADCP\MarshMayhem_NORTHADCP_adcp_level_2.mat'
%     'FallFrolic\PROCESSED_DATA\ADCP\FallFrolic_100825_NORTHADCP_adcp_level_2.mat'
%     'NovDep\PROCESSED_DATA\ADCP\NovDep_NORTHADCP_adcp_level_2.mat'
%     'JuneJamboree\PROCESSED_DATA\ADCP\JuneJamboree_NORTHADCP_adcp_level_2.mat'
%     'SummerFlood\PROCESSED_DATA\ADCP\SummerFlood_NORTHADCP_adcp_level_2.mat'};

%% Preallocate results structure

LagResults = struct();

%% Loop through observation periods

for i = 1:numel(obsNames)

    fprintf('\nProcessing %s...\n',obsNames{i})

    %% Load ADCP

    S1 = load(fullfile(main_dir,inputFiles1{i}));
    adcp_press = S1.adcp_level_1;

    S2 = load(fullfile(main_dir,inputFiles2{i}));
    adcp_curr = S2.adcp_level_2;

    %% Water level

    time = adcp_press.time;

    waterLevel = ...
        adcp_press.pressure - ...
        mean(adcp_press.pressure,'omitnan');

    %% Near-surface current

    nearSurfaceCurrent = ...
        adcp_curr.U_mapped(end,:)';

    %% Format data

    time = time(:);
    waterLevel = waterLevel(:);
    nearSurfaceCurrent = nearSurfaceCurrent(:);

    %% Interpolate NaNs in water level

    flag = isnan(waterLevel);

    if any(flag)

        waterLevel(flag) = interp1( ...
            time(~flag), ...
            waterLevel(~flag), ...
            time(flag), ...
            'linear', ...
            'extrap');

    end

    %% Interpolate NaNs in current

    flag = isnan(nearSurfaceCurrent);

    if any(flag)

        nearSurfaceCurrent(flag) = interp1( ...
            time(~flag), ...
            nearSurfaceCurrent(~flag), ...
            time(flag), ...
            'linear', ...
            'extrap');

    end

    %% Sampling interval

    dt = median(diff(time))*86400;

    %% Filter

    Nf = round(7200/dt);

    if mod(Nf,2) == 0
        Nf = Nf + 1;
    end

    flt = hanning(Nf);
    flt = flt/sum(flt);

    waterLevel_filt = ...
        conv(waterLevel,flt,'same');

    current_filt = ...
        conv(nearSurfaceCurrent,flt,'same');

    %% FIGURE 1
    % Water level and current


    fig1 = figure;

    ax1 = subplot(2,1,1);

%     plot(datetime(time,'ConvertFrom','datenum'), ...
%         waterLevel_filt, ...
%         '-k', ...
%         'LineWidth',1.5)

      plot(datetime(time,'ConvertFrom','datenum'), ...
            waterLevel, ...
            '-k', ...
            'LineWidth',1.5)

    ylabel('Water level')
    grid on

    title(obsNames{i}, ...
        'Interpreter','none')

    ax2 = subplot(2,1,2);

%     plot(datetime(time,'ConvertFrom','datenum'), ...
%         current_filt, ...
%         '-r', ...
%         'LineWidth',1.5)
% 

    plot(datetime(time,'ConvertFrom','datenum'), ...
        nearSurfaceCurrent, ...
        '-r', ...
        'LineWidth',1.5)



    ylabel('Near-surface U [m/s]')
    xlabel('Time')

    grid on

    linkaxes([ax1 ax2],'x')
    
    %% Save Figure 1

file_name = sprintf('%s_Signals',obsNames{i});

% exportgraphics(fig1, ...
%     fullfile(save_dir,[file_name '.png']), ...
%     'Resolution',300);
% 
% savefig(fig1, ...
%     fullfile(save_dir,[file_name '.fig']));

    %% Cross correlation

%     WL = ...
%         waterLevel_filt - ...
%         mean(waterLevel_filt,'omitnan');
% 
%     U = ...
%         current_filt - ...
%         mean(current_filt,'omitnan');

    WL = ...
        waterLevel - ...
        mean(waterLevel,'omitnan');

    U = ...
        nearSurfaceCurrent - ...
        mean(nearSurfaceCurrent,'omitnan');

    maxLag_seconds = 6*60*60;

    maxLag = round(maxLag_seconds/dt);

    [R,lags] = ...
        xcorr(WL,U,maxLag,'coeff');

    %% Find lag

    [maxR,idxMax] = max(R);

    lag_samples = lags(idxMax);

    time_lag_seconds = ...
        lag_samples*dt;

    time_lag_minutes = ...
        time_lag_seconds/60;

    [minR,idxMin] = min(R);

    lag_samples_negative = ...
        lags(idxMin);

    time_lag_negative_seconds = ...
        lag_samples_negative*dt;

    time_lag_negative_minutes = ...
        time_lag_negative_seconds/60;

%% Store lag results

LagResults(i).Observation = obsNames{i};

LagResults(i).R = R;
LagResults(i).lags = lags;
LagResults(i).dt = dt;

LagResults(i).MaxCorrelation = maxR;
LagResults(i).MaxCorrelationLag_min = time_lag_minutes;

LagResults(i).MinCorrelation = minR;
LagResults(i).MinCorrelationLag_min = time_lag_negative_minutes;

    %% Print results

    fprintf('Maximum correlation: %.3f at %.1f min\n', ...
        maxR,time_lag_minutes)

    fprintf('Minimum correlation: %.3f at %.1f min\n', ...
        minR,time_lag_negative_minutes)

    %% FIGURE 2
    % Cross correlation

    fig2 = figure;

    plot(lags*dt/60, ...
        R, ...
        'k', ...
        'LineWidth',1.5)

    hold on

    plot(time_lag_minutes, ...
        maxR, ...
        'ro', ...
        'MarkerFaceColor','r')

    plot(time_lag_negative_minutes, ...
        minR, ...
        'bo', ...
        'MarkerFaceColor','b')

    xline(0,'--')

    xlabel('Lag [min]')
    ylabel('Cross correlation')

    title(obsNames{i}, ...
        'Interpreter','none')

    grid on

    legend( ...
        'Correlation', ...
        'Maximum positive correlation', ...
        'Maximum negative correlation', ...
        'Zero lag', ...
        'Location','best')
    
    %% Save Figure 2

file_name = sprintf('%s_CrossCorrelation',obsNames{i});

% exportgraphics(fig2, ...
%     fullfile(save_dir,[file_name '.png']), ...
%     'Resolution',300);
% 
% savefig(fig2, ...
%     fullfile(save_dir,[file_name '.fig']));

    %% Find high and low water

    [idx_high,idx_low,~,~] = ...
        min_max(waterLevel_filt);

    %% Find slack water

    idx_slack = find( ...
        current_filt(1:end-1).* ...
        current_filt(2:end) <= 0);

    %% FIGURE 3
    % High/low water and slack current
    

    fig3 = figure;

    ax3 = subplot(2,1,1);

    plot(datetime(time,'ConvertFrom','datenum'), ...
        waterLevel_filt, ...
        'k', ...
        'LineWidth',1.5)

    hold on

    plot(datetime(time(idx_high), ...
        'ConvertFrom','datenum'), ...
        waterLevel_filt(idx_high), ...
        'ro')

    plot(datetime(time(idx_low), ...
        'ConvertFrom','datenum'), ...
        waterLevel_filt(idx_low), ...
        'bo')

    yline(0,'--')

    ylabel('Water level')

    title(obsNames{i}, ...
        'Interpreter','none')

    grid on

    hold off

    ax4 = subplot(2,1,2);

    plot(datetime(time,'ConvertFrom','datenum'), ...
        current_filt, ...
        'r', ...
        'LineWidth',1.5)

    hold on

    plot(datetime(time(idx_slack), ...
        'ConvertFrom','datenum'), ...
        current_filt(idx_slack), ...
        'ko', ...
        'MarkerFaceColor','y')

    yline(0,'--')

    ylabel('Near-surface U [m/s]')
    xlabel('Time')

    grid on

    hold off

    linkaxes([ax3 ax4],'x')
    
    %% Save Figure 3

% file_name = sprintf('%s_HighLow_SlackWater',obsNames{i});
% 
% exportgraphics(fig3, ...
%     fullfile(save_dir,[file_name '.png']), ...
%     'Resolution',300);
% 
% savefig(fig3, ...
%     fullfile(save_dir,[file_name '.fig']));

end


%% Combined lag across all observation periods

% Common lag axis in minutes
commonLag_min = (-360:1:360)';

R_all = nan(length(commonLag_min),numel(LagResults));

for i = 1:numel(LagResults)

    lag_min = ...
        LagResults(i).lags .* ...
        LagResults(i).dt ./ 60;

    R_all(:,i) = interp1( ...
        lag_min, ...
        LagResults(i).R, ...
        commonLag_min, ...
        'linear', ...
        NaN);

end

%% Mean correlation across all observation periods

R_mean = mean(R_all,2,'omitnan');

%% Find overall positive correlation lag

[maxR_all,idxMax_all] = max(R_mean);

overallLag_positive_min = ...
    commonLag_min(idxMax_all);

%% Find overall negative correlation lag

[minR_all,idxMin_all] = min(R_mean);

overallLag_negative_min = ...
    commonLag_min(idxMin_all);

%% Print

fprintf('\nCOMBINED RESULTS ACROSS ALL OBSERVATION PERIODS\n')

fprintf( ...
    'Maximum mean correlation: %.3f at %.1f min\n', ...
    maxR_all, ...
    overallLag_positive_min)

fprintf( ...
    'Minimum mean correlation: %.3f at %.1f min\n', ...
    minR_all, ...
    overallLag_negative_min)

%% Plot combined cross correlation

fig4 = figure;
hold on

% Individual observation periods
for i = 1:numel(LagResults)

    plot( ...
        commonLag_min, ...
        R_all(:,i), ...
        'LineWidth',0.75, ...
        'HandleVisibility','off');

end

% Mean across all observation periods
plot( ...
    commonLag_min, ...
    R_mean, ...
    'k', ...
    'LineWidth',2);

% plot( ...
%     overallLag_positive_min, ...
%     maxR_all, ...
%     'ro', ...
%     'MarkerFaceColor','r');

plot( ...
    overallLag_negative_min, ...
    minR_all, ...
    'bo', ...
    'MarkerFaceColor','b');

xline(0,'--')

xline(overallLag_negative_min, ...
    '--b', ...
    sprintf('Minimum at %.0f min',overallLag_negative_min), ...
    'LineWidth',1.5);

xlabel('Lag [min]')
ylabel('Cross correlation')

title('Mean Cross Correlation Across All Observation Periods')

% legend( ...
%     'Mean correlation', ...
%     'Maximum positive correlation', ...
%     'Maximum negative correlation', ...
%     'Zero lag', ...
%     'Location','best')

legend( ...
    'Mean correlation', ...
    'Maximum negative correlation', ...
    'Zero lag', ...
    'Location','best')

grid on

 file_name = sprintf('mean_lag');

% exportgraphics(fig4, ...
%     fullfile(save_dir,[file_name '.png']), ...
%     'Resolution',300);
% 
% savefig(fig4, ...
%     fullfile(save_dir,[file_name '.fig']));

%% Display lag results

Observation = string({LagResults.Observation})';

MaxCorrelation = ...
    [LagResults.MaxCorrelation]';

MaxCorrelationLag_min = ...
    [LagResults.MaxCorrelationLag_min]';

MinCorrelation = ...
    [LagResults.MinCorrelation]';

MinCorrelationLag_min = ...
    [LagResults.MinCorrelationLag_min]';

LagTable = table( ...
    Observation, ...
    MaxCorrelation, ...
    MaxCorrelationLag_min, ...
    MinCorrelation, ...
    MinCorrelationLag_min);

disp(LagTable)

%% Min/max function

function [idx_max,idx_min,cnt_max,cnt_min] = min_max(y)

N = length(y);

idx_max = [];
idx_min = [];

cnt_max = 0;
cnt_min = 0;

for ii = 2:N-1

    if y(ii-1) > y(ii) && ...
            y(ii+1) > y(ii)

        cnt_min = cnt_min + 1;

        idx_min(cnt_min) = ii;

    elseif y(ii-1) < y(ii) && ...
            y(ii+1) < y(ii)

        cnt_max = cnt_max + 1;

        idx_max(cnt_max) = ii;

    end

end

end