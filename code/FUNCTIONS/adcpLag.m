clear all
close all
%% Set dir

main_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\CB_YachtBasin\';

%% Load ADCP

% inputDir  = 'C:\Users\Ben\OneDrive - UNC-Wilmington\CB_YachtBasin\MarshMadness\PROCESSED_DATA\ADCP\';
inputDir  = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\CB_YachtBasin\';
inputFiles1 = {'MarshMadness\PROCESSED_DATA\ADCP\MarshMadness_NORTHADCP_adcp_level_1.mat'
               'MarshMayhem\PROCESSED_DATA\ADCP\MarshMayhem_NORTHADCP_adcp_level_1.mat'
               'FallFrolic\PROCESSED_DATA\ADCP\FallFrolic_100825_NORTHADCP_adcp_level_1.mat'
               'NovDep\PROCESSED_DATA\ADCP\NovDep_NORTHADCP_adcp_level_1.mat'
               'JuneJamboree\PROCESSED_DATA\ADCP\JuneJamboree_NORTHADCP_adcp_level_1.mat'
               'SummerFlood\PROCESSED_DATA\ADCP\SummerFlood_NORTHADCP_adcp_level_1.mat'};

inputFiles2 = {'MarshMadness\PROCESSED_DATA\ADCP\MarshMadness_NORTHADCP_adcp_level_2.mat'
               'MarshMayhem\PROCESSED_DATA\ADCP\MarshMayhem_NORTHADCP_adcp_level_2.mat'
               'FallFrolic\PROCESSED_DATA\ADCP\FallFrolic_100825_NORTHADCP_adcp_level_2.mat'
               'NovDep\PROCESSED_DATA\ADCP\NovDep_NORTHADCP_adcp_level_2.mat'
               'JuneJamboree\PROCESSED_DATA\ADCP\JuneJamboree_NORTHADCP_adcp_level_2.mat'
               'SummerFlood\PROCESSED_DATA\ADCP\SummerFlood_NORTHADCP_adcp_level_2.mat'};


adcp_press = load(fullfile(inputDir, inputFile1));

adcp_press = adcp_press.adcp_level_1;

adcp_curr = load(fullfile(inputDir, inputFile2));

adcp_curr = adcp_curr.adcp_level_2;

%% Water level

time = adcp_press.time;

waterLevel = adcp_press.pressure - mean(adcp_press.pressure,'omitnan');

%% Near-surface current

nearSurfaceCurrent = adcp_curr.U_mapped(end,:)';

%% Format data

time               = time(:);
waterLevel         = waterLevel(:);
nearSurfaceCurrent = nearSurfaceCurrent(:);

flag = isnan(waterLevel);

waterLevel(flag) = interp1( ...
    time(~flag), ...
    waterLevel(~flag), ...
    time(flag), ...
    'linear','extrap');

flag = isnan(nearSurfaceCurrent);

nearSurfaceCurrent(flag) = interp1( ...
    time(~flag), ...
    nearSurfaceCurrent(~flag), ...
    time(flag), ...
    'linear','extrap');

dt = median(diff(time))*86400;

%% Filter

Nf = round(7200/dt);

if mod(Nf,2) == 0
    Nf = Nf + 1;
end

flt = hanning(Nf);
flt = flt/sum(flt);

waterLevel_filt = conv(waterLevel,flt,'same');

current_filt = conv(nearSurfaceCurrent,flt,'same');

%% Plot signals

fig1 = figure;

% Top panel: Water level
subplot(2,1,1)

plot(datetime(time,'ConvertFrom','datenum'), ...
    waterLevel_filt, ...
    '-k', ...
    'LineWidth',1.5)

ylabel('Water level')

grid on

% Bottom panel: Near-surface current
subplot(2,1,2)

plot(datetime(time,'ConvertFrom','datenum'), ...
    current_filt, ...
    '-r', ...
    'LineWidth',1.5)

ylabel('Near-surface U [m/s]')

xlabel('Time')

grid on

%% Cross correlation

WL = waterLevel_filt - mean(waterLevel_filt,'omitnan');
U  = current_filt - mean(current_filt,'omitnan');

maxLag_seconds = 6*60*60;

maxLag = round(maxLag_seconds/dt);

[R,lags] = xcorr(WL,U,maxLag,'coeff');

%% Find lag

[maxR,idxMax] = max(R);

lag_samples = lags(idxMax);

time_lag_seconds = lag_samples*dt;

time_lag_minutes = time_lag_seconds/60

[minR,idxMin] = min(R);

lag_samples_negative = lags(idxMin);

time_lag_negative_seconds = lag_samples_negative*dt;

time_lag_negative_minutes = time_lag_negative_seconds/60

%% Plot correlation

fig2 = figure;

plot(lags*dt/60,R,'k','LineWidth',1.5)

hold on

plot(time_lag_minutes,maxR, ...
    'ro', ...
    'MarkerFaceColor','r')

plot(time_lag_negative_minutes,minR, ...
    'bo', ...
    'MarkerFaceColor','b')

xline(0,'--')

xlabel('Lag [min]')

ylabel('Cross correlation')

grid on

legend( ...
    'Correlation', ...
    'Maximum positive correlation', ...
    'Maximum negative correlation', ...
    'Zero lag')

%% Find high and low water

[idx_high,idx_low,~,~] = min_max(waterLevel_filt);

%% Find slack water

idx_slack = find( ...
    current_filt(1:end-1).*current_filt(2:end) <= 0);

%% Plot water level and slack current

fig3 = figure;

yyaxis left

subplot(2,1,1);

plot(datetime(time,'ConvertFrom','datenum'), ...
    waterLevel_filt, ...
    'k', ...
    'LineWidth',1.5)

hold on

plot(datetime(time(idx_high),'ConvertFrom','datenum'), ...
    waterLevel_filt(idx_high), ...
    'ro')

plot(datetime(time(idx_low),'ConvertFrom','datenum'), ...
    waterLevel_filt(idx_low), ...
    'bo')

yline(0,'--');

grid on

hold off

% yyaxis right

subplot(2,1,2);

plot(datetime(time,'ConvertFrom','datenum'), ...
    current_filt, ...
    'r', ...
    'LineWidth',1.5)

hold on

plot(datetime(time(idx_slack),'ConvertFrom','datenum'), ...
    current_filt(idx_slack), ...
    'ko', ...
    'MarkerFaceColor','y')

yline(0,'--')

xlabel('Time')

grid on
hold off

%% Min/max function

function [idx_max,idx_min,cnt_max,cnt_min] = min_max(y)

N = length(y);

idx_max = [];
idx_min = [];

cnt_max = 0;
cnt_min = 0;

for ii = 2:N-1

    if y(ii-1)>y(ii) && y(ii+1)>y(ii)

        cnt_min = cnt_min+1;

        idx_min(cnt_min) = ii;

    elseif y(ii-1)<y(ii) && y(ii+1)<y(ii)

        cnt_max = cnt_max+1;

        idx_max(cnt_max) = ii;

    end

end

end