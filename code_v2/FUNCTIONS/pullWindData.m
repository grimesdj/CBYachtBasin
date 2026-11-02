function out = pullWindData()
%% pullWindData
% Downloads wind from a NOAA station near Carolina Beach for each observation period
% and saves it, so wind can be compared with the drifter tracks later.
%
% WHAT IT DOES, IN ORDER
%   STEP 1  Settings: which station, and the dates of each observation period.
%   STEP 2  For each period, asks NOAA for the wind record (6-minute values).
%   STEP 3  Adds east/north wind components (direction the wind blows TOWARD,
%           the same way drifter velocity is described).
%   STEP 4  Prints a short summary and saves the results.
%
% WIND CONVENTIONS (read this before using the numbers)
%   speed  : m/s
%   dir    : degrees the wind blows FROM, clockwise from north (90 = from the east).
%   gust   : m/s
%   u_to, v_to : east and north components of the wind blowing TOWARD, m/s.
%   Times are UTC. Drifter times are UTC.
%
% WHAT IT WRITES  (folder CB_YachtBasin\WIND_DATA, created if missing)
%   Wind_AllPeriods.mat     struct Wind, one entry per period
%   Wind_<period>.csv       the same data as a table, one file per period
%
% Needs internet access from MATLAB. It does not use the cfg.

%% STEP 1: settings
stationID   = '8658163';              % NOAA CO-OPS station number
stationName = 'Wrightsville Beach';   % about 20 km north-northeast of Carolina Beach
padHours    = 12;                     % extra hours added before and after each period

% One row per observation period: name, first day, last day, as 'YYYY-MM-DD'.
% Type the dates in. A row with empty dates is skipped.
periods = {
    'MarshMadness',       '2025-03-26', '2025-03-28'
    'MarshMayhem',        '2025-05-27', '2025-05-28'
    'FallFrolic_100825',  '2025-10-08', '2025-10-08'
    'FallFrolic_101025',  '2025-10-10', '2025-10-10'
    'NovDep',             '2025-11-07', '2025-11-07'
    'JuneJamboree',       '2026-06-14', '2026-06-14'
    'SummerFlood',        '2026-08-11', '2026-08-11'
    };

roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', ...
         'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if isempty(k), error('Neither root path contains CB_YachtBasin on this machine.'); end
outDir = fullfile(roots{k}, 'CB_YachtBasin', 'WIND_DATA');
if ~isfolder(outDir), mkdir(outDir); end

%% STEP 2 and 3: download each period
Wind = struct('name',{}, 'station',{}, 'first',{}, 'last',{}, 'data',{});

for p = 1:size(periods,1)
    name = periods{p,1};
    if isempty(periods{p,2}) || isempty(periods{p,3})
        fprintf('%s: dates not filled in, skipped\n', name);
        continue
    end

    % window to request = first day 00:00 minus padding, to last day 24:00 plus padding
    t1 = datetime(periods{p,2}, 'InputFormat','yyyy-MM-dd') - hours(padHours);
    t2 = datetime(periods{p,3}, 'InputFormat','yyyy-MM-dd') + days(1) + hours(padHours);

    % NOAA allows about one month per request, so long periods are done in 30-day pieces
    T = [];
    t = t1;
    while t < t2
        tEnd = min(t + days(30), t2);
        T = [T; getWind(stationID, t, tEnd)]; %#ok<AGROW>
        t = tEnd;
    end

    if isempty(T)
        fprintf('%s: NOAA returned no wind data for this window\n', name);
        continue
    end

    [~, iu] = unique(T.time);          % remove duplicate rows where pieces meet
    T = T(iu,:);

    % wind blows TOWARD the opposite of the direction it comes FROM
    T.u_to = -T.speed .* sind(T.dir);
    T.v_to = -T.speed .* cosd(T.dir);

    Wind(end+1).name    = name; %#ok<AGROW>
    Wind(end).station   = stationID;
    Wind(end).first     = t1;
    Wind(end).last      = t2;
    Wind(end).data      = T;

    writetable(T, fullfile(outDir, ['Wind_' name '.csv']));
end

%% STEP 4: summary and save
fprintf('\n=== Wind from %s (NOAA %s), UTC ===\n', stationName, stationID);
fprintf('%-20s %8s %9s %10s %10s %12s\n', 'period', 'records', 'missing', 'mean m/s', 'max gust', 'mean toward');
for i = 1:numel(Wind)
    T = Wind(i).data;
    nExpected = round(hours(Wind(i).last - Wind(i).first) * 10) + 1;   % 10 records per hour
    pctMissing = 100 * (1 - sum(~isnan(T.speed)) / nExpected);
    meanU = mean(T.u_to, 'omitnan');
    meanV = mean(T.v_to, 'omitnan');
    towardDeg = mod(atan2d(meanU, meanV), 360);       % compass bearing the average wind blows toward
    fprintf('%-20s %8d %8.0f%% %10.1f %10.1f %11.0f deg\n', Wind(i).name, height(T), ...
        pctMissing, mean(T.speed,'omitnan'), max(T.gust), towardDeg);
end

save(fullfile(outDir, 'Wind_AllPeriods.mat'), 'Wind', 'stationID', 'stationName');
fprintf('\nSaved Wind_AllPeriods.mat and one csv per period to %s\n', outDir);

out = Wind;
end


%% ========================= Local function =========================

function T = getWind(stationID, t1, t2)
% One request to NOAA (no more than about 30 days). Returns a table with
% time (UTC), speed (m/s), dir (deg, wind FROM) and gust (m/s). Empty if nothing came back.

url = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter';
raw = webread(url, ...
    'product',    'wind', ...
    'station',    stationID, ...
    'begin_date', char(t1, 'yyyyMMdd HH:mm'), ...
    'end_date',   char(t2, 'yyyyMMdd HH:mm'), ...
    'time_zone',  'gmt', ...
    'units',      'metric', ...
    'format',     'json', ...
    'application','CB_YachtBasin_thesis', ...
    weboptions('Timeout', 60));

T = [];
if isfield(raw, 'error')
    fprintf('  NOAA message: %s\n', raw.error.message);
    return
end
if ~isfield(raw, 'data') || isempty(raw.data), return, end

d = raw.data;
if iscell(d), d = [d{:}]; end            % in case the answer comes back as a cell array

time  = datetime({d.t}, 'InputFormat','yyyy-MM-dd HH:mm', 'TimeZone','UTC')';
speed = str2double({d.s})';
dir   = str2double({d.d})';
gust  = str2double({d.g})';
T = table(time, speed, dir, gust);
end
