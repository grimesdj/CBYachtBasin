function [drifters_level_2, Stuckinfo] = flagStuckDrifters(cfg, drifters_level_2)
%% flagStuckDrifters
% Adds a logical field  Stuck  to every drifter in the level 2 file. Nothing is removed
% and no existing field is changed, so existing scripts keep working.
%
% A 1 Hz sample is flagged Stuck when, in a centered window of win_min minutes:
%   (1) the position stays inside a box whose diagonal is < radius_m, AND
%   (2) the median 1 Hz Doppler speed is < dop_max (m/s),
% and the flagged run lasts at least min_run_min minutes.
% These are the same settings and the same rule as stuckDrifterDiagnostic.m
% (keep the two in sync if you change either).
%
% Samples inside the first/last half window take the flag of the nearest sample that has
% a full window, so a drifter still stuck when it was pulled is flagged to the end.
%
% Use in main_driver, between findTidalPhaseandRegion and the save of level 2:
%     [TidalPhase, drifters_level_2] = findTidalPhaseandRegion(cfg);
%     [drifters_level_2, StuckInfo]  = flagStuckDrifters(cfg, drifters_level_2);
%     save(... '_drifters_level_2.mat' ..., 'drifters_level_2');
% The input struct is optional: with only cfg, the function loads the saved level 2 file
% (cfg.out.drifters_data, [cfg.name '_drifters_level_2.mat']) and you save the result yourself.
%
% Output info: one row per track with the stuck share and trailing stuck time.

%% Settings
win_min     = 10;      % window length, minutes
radius_m    = 10;      % box diagonal, metres
dop_max     = 0.02;    % median Doppler speed, m/s
min_run_min = 10;      % shortest stuck run kept, minutes

%% Get the level 2 struct
if nargin < 2 || isempty(drifters_level_2)
    drifter_filename = fullfile(cfg.out.drifters_data, [cfg.name '_drifters_level_2.mat']);
    S = load(drifter_filename);
    if isfield(S,'drifters_level_2')
        drifters_level_2 = S.drifters_level_2;
    elseif isfield(S,'trimmed')
        drifters_level_2 = S.trimmed;
    else
        error('File does not contain drifters_level_2 or trimmed')
    end
end

[nD,nC] = size(drifters_level_2);
W      = round(win_min*60);
minRun = round(min_run_min*60);

rows = struct('Deployment',{},'Track',{},'ID',{},'Duration_h',{},'StuckPct',{},'TrailingStuck_h',{});
nAll = 0; nStuck = 0; dtWarn = false;

for d = 1:nD
    for c = 1:nC

        D = drifters_level_2(d,c);

        % default: empty flag for tracks that cannot be assessed
        drifters_level_2(d,c).Stuck = [];

        if isempty(D.East) || isempty(D.North) || isempty(D.Date_and_time) || isempty(D.Speed)
            continue
        end

        t   = seconds(D.Date_and_time(:) - D.Date_and_time(1));
        x   = D.East(:);
        y   = D.North(:);
        spd = D.Speed(:) * 0.51444;                 % knots -> m/s

        n = min([numel(t) numel(x) numel(y) numel(spd)]);
        if n < 2, continue, end

        if abs(median(diff(t(1:n))) - 1) > 0.2, dtWarn = true; end

        f = stuckFlag(x(1:n), y(1:n), spd(1:n), W, radius_m, dop_max, minRun);

        % Stuck has the same length as the track's East/North vectors
        flag = false(numel(D.East),1);
        flag(1:n) = f;
        drifters_level_2(d,c).Stuck = flag;

        nAll   = nAll + n;
        nStuck = nStuck + sum(f);

        [s,e] = runs(f);
        trail = 0;
        if ~isempty(e) && e(end) == n
            trail = (t(e(end)) - t(s(end)))/3600;
        end
        if isfield(D,'ID') && ~isempty(D.ID), id = D.ID; else, id = NaN; end
        rows(end+1) = struct('Deployment',d,'Track',c,'ID',id, ...
            'Duration_h',t(n)/3600,'StuckPct',100*sum(f)/n,'TrailingStuck_h',trail); %#ok<AGROW>
    end
end

if dtWarn
    warning('flagStuckDrifters: some tracks are not sampled at about 1 Hz; windows are in samples assuming 1 Hz.');
end

Stuckinfo = struct2table(rows);
fprintf('flagStuckDrifters (%s): %d tracks, %.1f%% of samples flagged stuck (win %d min, box %g m, Doppler < %g m/s)\n', ...
    cfg.name, numel(rows), 100*nStuck/max(nAll,1), win_min, radius_m, dop_max);
end


%% ========================= Local functions =========================

function f = stuckFlag(x, y, spd, W, radius, dopmax, minRun)
n = numel(x);
f = false(n,1);
if n < W + 1, return, end

bx = movmax(x,W,'omitnan','Endpoints','fill') - movmin(x,W,'omitnan','Endpoints','fill');
by = movmax(y,W,'omitnan','Endpoints','fill') - movmin(y,W,'omitnan','Endpoints','fill');
dg = hypot(bx,by);
sm = movmedian(spd,W,'omitnan','Endpoints','fill');

ok = ~isnan(dg) & ~isnan(sm);
g  = double(dg < radius & sm < dopmax);
g(~ok) = NaN;
g = fillmissing(g,'nearest');          % edges take the nearest full-window value
g(isnan(g)) = 0;
f = g > 0;
f = dropShortRuns(f, minRun);
end


function f = dropShortRuns(f, minLen)
[s,e] = runs(f);
for k = 1:numel(s)
    if e(k) - s(k) + 1 < minLen
        f(s(k):e(k)) = false;
    end
end
end


function [s,e] = runs(f)
d = diff([false; f(:); false]);
s = find(d == 1);
e = find(d == -1) - 1;
end
