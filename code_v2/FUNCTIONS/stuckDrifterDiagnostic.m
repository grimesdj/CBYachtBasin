function out = stuckDrifterDiagnostic(cfg)
%% stuckDrifterDiagnostic
% How much of each drifter track is spent stuck or grounded (not moving) before the
% pull time, and how much that lowers regional mean speeds.
%
% A 1 Hz sample is flagged "stuck" when, in a centered window of win_min minutes:
%   (1) the position stays inside a box whose diagonal is < radius_m, AND
%   (2) the median 1 Hz Doppler speed is < dop_max (m/s),
% and the flagged run lasts at least min_run_min minutes.
% Thresholds come from the stationary deployment (95% position radius about own mean
% 4.2 m; Doppler 95th percentile at rest 2.1 cm/s), so a drifter sitting still should
% be flagged and one drifting at a few cm/s should not. A sensitivity table is printed.
%
% Samples inside the first/last half window use the flag of the nearest sample that has
% a full window (so a drifter still stuck at the end of its track is flagged to the end).
%
% Uses the level 2 drifter files (1 Hz: East, North, Speed in knots, Regions, Date_and_time).
% Saves to THESIS\THESIS_DRAFTING\FIGURES\STUCK_DRIFTERS:
%   StuckDrifters_timeline   one row per track, red = stuck
%   StuckDrifters_byRegion   stuck share and mean Doppler speed by region
%   StuckDrifters_byPhase    stuck share by tidal phase, and phase when stuck runs start/end
%   StuckDrifters_tracks.csv one row per track
%   StuckDrifters_flags.mat  per-sample stuck flags, for use in later analysis

close all
if nargin < 1, cfg = []; end %#ok<NASGU>

%% Settings
win_min     = 10;      % window length, minutes
radius_m    = 10;      % box diagonal, metres
dop_max     = 0.02;    % median Doppler speed, m/s
min_run_min = 10;      % shortest stuck run kept, minutes
report_min  = 30;      % report tracks that end with a stuck run at least this long
sens_radius = [5 10 20];
sens_dop    = [0.01 0.02 0.03];

names   = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};
folders = {'MarshMadness','MarshMayhem','FallFrolic','FallFrolic','NovDep','JuneJamboree','SummerFlood'};

roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', ...
    'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if isempty(k), error('Neither root path contains CB_YachtBasin on this machine.'); end
root_data = roots{k};
fprintf('Using root: %s\n', root_data);

save_dir = fullfile(root_data,'THESIS','THESIS_DRAFTING','FIGURES','STUCK_DRIFTERS');
if ~isfolder(save_dir), mkdir(save_dir); end

%% Load every track
TR = struct('Period',{},'Deployment',{},'Track',{},'ID',{},'t',{},'x',{},'y',{},'spd',{},'reg',{},'ph',{},'stk',{});
dtWarn = false;

for i = 1:numel(names)
    
    fn = fullfile(root_data,'CB_YachtBasin',folders{i},'PROCESSED_DATA_V2','DRIFTERS', ...
        [names{i} '_drifters_level_2.mat']);
    if ~isfile(fn)
        warning('%s: file not found, skipped', names{i});
        continue
    end
    
    S = load(fn);
    if isfield(S,'drifters_level_2')
        L2 = S.drifters_level_2;
    elseif isfield(S,'trimmed')
        L2 = S.trimmed;
    else
        warning('%s: no drifters_level_2 or trimmed variable, skipped', names{i});
        continue
    end
    
    [nD,nC] = size(L2);
    for d = 1:nD
        for c = 1:nC
            D = L2(d,c);
            if isempty(D.East) || isempty(D.North) || isempty(D.Date_and_time) || isempty(D.Speed)
                continue
            end
            
            t   = seconds(D.Date_and_time(:) - D.Date_and_time(1));
            x   = D.East(:);
            y   = D.North(:);
            spd = D.Speed(:) * 0.51444;          % knots -> m/s
            if isfield(D,'Regions') && ~isempty(D.Regions)
                reg = double(D.Regions(:));
            else
                reg = nan(size(x));
            end
            
            if ~isfield(D,'Stuck') || numel(D.Stuck) ~= numel(D.East)
                error('%s has no usable Stuck field in level 2. Run flagStuckDrifters first.', names{i});
            end
            stk = logical(D.Stuck(:));           % the flag flagStuckDrifters wrote into level 2
            
            n = min([numel(t) numel(x) numel(y) numel(spd) numel(reg) numel(stk)]);
            if n < 2, continue, end
            t = t(1:n); x = x(1:n); y = y(1:n); spd = spd(1:n); reg = reg(1:n); stk = stk(1:n);
            
            ph = strings(n,1);
            if isfield(D,'TidePhase') && ~isempty(D.TidePhase)
                pp = string(D.TidePhase(:));
                m  = min(n, numel(pp));
                ph(1:m) = pp(1:m);
            end
            
            if abs(median(diff(t)) - 1) > 0.2, dtWarn = true; end
            
            if isfield(D,'ID') && ~isempty(D.ID), id = D.ID; else, id = NaN; end
            
            TR(end+1) = struct('Period',i,'Deployment',d,'Track',c,'ID',id, ...
                't',t,'x',x,'y',y,'spd',spd,'reg',reg,'ph',ph,'stk',stk); %#ok<AGROW>        end
        end
    end
end


nT = numel(TR);
fprintf('\nLoaded %d tracks.\n', nT);
if dtWarn
    warning('Some tracks are not sampled at about 1 Hz; window lengths are in samples assuming 1 Hz.');
end

W      = round(win_min*60);
minRun = round(min_run_min*60);

%% Flag stuck samples and summarise each track
Flag = cell(nT,1);
rows = struct('Period',{},'Deployment',{},'ID',{},'Duration_h',{},'StuckPct',{}, ...
    'TrailingStuck_h',{},'MidTrackStuck_h',{},'NRuns',{},'StuckRegion',{});

for q = 1:nT
    T = TR(q);
    f = T.stk;                           % Stuck flag from level 2 (flagStuckDrifters)
    Flag{q} = f;
    
    n   = numel(f);
    dur = T.t(end)/3600;
    [s,e] = runs(f);
    
    trail = 0;
    if ~isempty(e) && e(end) == n
        trail = (T.t(e(end)) - T.t(s(end)))/3600;
    end
    stuckH = sum(f) * median(diff(T.t)) / 3600;
    
    sr = NaN;
    m = f & ~isnan(T.reg);
    if any(m), sr = mode(T.reg(m)); end
    
    rows(end+1) = struct('Period',string(strrep(names{T.Period},'_',' ')), 'Deployment',T.Deployment, ...
        'ID',T.ID, 'Duration_h',dur, 'StuckPct',100*sum(f)/n, ...
        'TrailingStuck_h',trail, 'MidTrackStuck_h',max(stuckH - trail,0), ...
        'NRuns',numel(s), 'StuckRegion',sr); %#ok<AGROW>
end

TT = struct2table(rows);
writetable(TT, fullfile(save_dir,'StuckDrifters_tracks.csv'));

%% Print results
nAll   = sum(cellfun(@numel,Flag));
nStuck = sum(cellfun(@sum,Flag));
fprintf('\n=== Overall (win %d min, box %g m, Doppler < %g m/s, run >= %d min) ===\n', ...
    win_min, radius_m, dop_max, min_run_min);
fprintf('Samples flagged stuck: %.1f%% of %d samples (%d tracks)\n', 100*nStuck/nAll, nAll, nT);
fprintf('Tracks with any stuck run: %d (%.0f%%)\n', sum(TT.NRuns>0), 100*mean(TT.NRuns>0));
hasTrail = TT.TrailingStuck_h*60 >= report_min;
fprintf('Tracks ending with a stuck run >= %d min: %d (%.0f%%); median trailing stuck time %.1f h\n', ...
    report_min, sum(hasTrail), 100*mean(hasTrail), median(TT.TrailingStuck_h(hasTrail)));
fprintf('Median track duration %.1f h\n', median(TT.Duration_h));

fprintf('\n=== By observation period ===\n');
fprintf('%-20s %7s %12s %22s\n','Period','tracks','% stuck','tracks ending stuck');
for i = 1:numel(names)
    qi = [TR.Period] == i;
    if ~any(qi), continue, end
    ns = sum(cellfun(@sum,Flag(qi))); na = sum(cellfun(@numel,Flag(qi)));
    fprintf('%-20s %7d %11.1f%% %15d (%.0f%%)\n', strrep(names{i},'_',' '), sum(qi), ...
        100*ns/na, sum(hasTrail(qi)), 100*mean(hasTrail(qi)));
end

regs = 0:7;
nReg = numel(regs);
allSp = nan(nReg,1); movSp = nan(nReg,1); pctStuck = nan(nReg,1); nSamp = zeros(nReg,1);
for a = 1:nReg
    sAll = 0; nA = 0; sMov = 0; nM = 0; nS = 0;
    for q = 1:nT
        m  = TR(q).reg == regs(a) & ~isnan(TR(q).spd);
        f  = Flag{q};
        sAll = sAll + sum(TR(q).spd(m));         nA = nA + sum(m);
        sMov = sMov + sum(TR(q).spd(m & ~f));    nM = nM + sum(m & ~f);
        nS = nS + sum(m & f);
    end
    nSamp(a) = nA;
    if nA > 0
        allSp(a) = sAll/nA; pctStuck(a) = 100*nS/nA;
    end
    if nM > 0, movSp(a) = sMov/nM; end
end

fprintf('\n=== By region (1 Hz Doppler speed) ===\n');
fprintf('%-7s %10s %9s %14s %16s\n','Region','samples','% stuck','mean all m/s','mean moving m/s');
for a = 1:nReg
    fprintf('%-7d %10d %8.1f%% %14.3f %16.3f\n', regs(a), nSamp(a), pctStuck(a), allSp(a), movSp(a));
end


%% Tidal phase: do stuck samples cluster at a particular phase?
% If drifters ground on falling/low water and float off on rising water, stuck samples
% should cluster at Low (and Falling) and run starts/ends should line up with the phases.
phases  = ["High","Falling","Low","Rising"];
isMarsh = ismember([TR.Period], [1 2]);       % MarshMadness, MarshMayhem
nPh = numel(phases);
phAll = zeros(nPh,2); phMar = zeros(nPh,2); phOth = zeros(nPh,2);   % [samples stuck]
for p = 1:nPh
    phAll(p,:) = phaseCounts(TR, Flag, phases(p), true(1,nT));
    phMar(p,:) = phaseCounts(TR, Flag, phases(p), isMarsh);
    phOth(p,:) = phaseCounts(TR, Flag, phases(p), ~isMarsh);
end
totMar = sum(cellfun(@numel, Flag(isMarsh)));
totOth = sum(cellfun(@numel, Flag(~isMarsh)));

fprintf('\n=== By tidal phase ===\n');
fprintf('%%time = share of that group''s samples spent in the phase; %%stuck = share of those samples flagged.\n');
fprintf('%-9s | %-22s | %-22s | %-22s\n','','All tracks','MarshMadness+Mayhem','Other periods');
fprintf('%-9s | %9s %12s | %9s %12s | %9s %12s\n','Phase','%time','%stuck','%time','%stuck','%time','%stuck');
for p = 1:nPh
    fprintf('%-9s | %8.1f%% %11.1f%% | %8.1f%% %11.1f%% | %8.1f%% %11.1f%%\n', phases(p), ...
        100*phAll(p,1)/nAll, 100*phAll(p,2)/max(phAll(p,1),1), ...
        100*phMar(p,1)/totMar, 100*phMar(p,2)/max(phMar(p,1),1), ...
        100*phOth(p,1)/totOth, 100*phOth(p,2)/max(phOth(p,1),1));
end

% Phase at the start of each stuck run (grounding) and just after its end (floating off again)
startCnt = zeros(nPh,1); endCnt = zeros(nPh,1); nRunsEnd = 0; nRunsStart = 0;
for q = 1:nT
    [s,e] = runs(Flag{q});
    n = numel(Flag{q});
    for j = 1:numel(s)
        a = find(phases == TR(q).ph(s(j)), 1);
        if ~isempty(a), startCnt(a) = startCnt(a) + 1; end
        nRunsStart = nRunsStart + 1;
        if e(j) < n
            b = find(phases == TR(q).ph(min(e(j)+1,n)), 1);
            if ~isempty(b), endCnt(b) = endCnt(b) + 1; end
            nRunsEnd = nRunsEnd + 1;
        end
    end
end
fprintf('\n=== Phase when stuck runs start and end (all tracks) ===\n');
fprintf('%-9s %14s %22s\n','Phase','runs starting','runs ending (floating off)');
for p = 1:nPh
    fprintf('%-9s %14d %22d\n', phases(p), startCnt(p), endCnt(p));
end
fprintf('(%d runs in total; %d of them end before the end of the track)\n', nRunsStart, nRunsEnd);

fprintf('\n=== Sensitivity: %% of samples flagged stuck ===\n');
fprintf('rows = box diagonal (m), columns = median Doppler limit (m/s)\n');
fprintf('%8s', ''); fprintf('%10.3f', sens_dop); fprintf('\n');
Sens = nan(numel(sens_radius), numel(sens_dop));
for a = 1:numel(sens_radius)
    for b = 1:numel(sens_dop)
        ns = 0;
        for q = 1:nT
            ns = ns + sum(stuckFlag(TR(q).x, TR(q).y, TR(q).spd, W, sens_radius(a), sens_dop(b), minRun));
        end
        Sens(a,b) = 100*ns/nAll;
    end
    fprintf('%6g m', sens_radius(a)); fprintf('%10.1f', Sens(a,:)); fprintf('\n');
end

%% Figure 1: timeline
[~,ord] = sortrows([[TR.Period]' [TR.Deployment]' [TR.Track]']);
f1 = figure('Units','pixels','Position',[100 50 1300 1000],'Color','w');
ax = axes(f1); hold(ax,'on'); box(ax,'on')
rowOf = zeros(nT,1); rowOf(ord) = 1:nT;
for q = 1:nT
    r = rowOf(q);
    plot(ax, [0 TR(q).t(end)/3600], [r r], '-', 'Color',[0.80 0.80 0.80], 'LineWidth',2);
    [s,e] = runs(Flag{q});
    for j = 1:numel(s)
        plot(ax, [TR(q).t(s(j)) TR(q).t(e(j))]/3600, [r r], '-', 'Color',[0.84 0.15 0.16], 'LineWidth',2.5);
    end
end
mid = nan(1,numel(names)); edges = [];
for i = 1:numel(names)
    rr = rowOf([TR.Period] == i);
    if isempty(rr), continue, end
    mid(i) = mean(rr);
    edges(end+1) = max(rr) + 0.5; %#ok<AGROW>
end
for ed = edges(1:end-1), yline(ax, ed, '-', 'Color',[0.3 0.3 0.3]); end
ok = ~isnan(mid);
yticks(ax, mid(ok)); yticklabels(ax, strrep(names(ok),'_',' '));
set(ax,'YDir','reverse','FontSize',11); ylim(ax,[0 nT+1]);
xlabel(ax,'Hours since release');
title(ax, sprintf('Stuck (red) and moving (gray) time per track: %.1f%% of samples flagged', 100*nStuck/nAll), 'FontSize',13);
saveFig(f1, save_dir, 'StuckDrifters_timeline');

%% Figure 2: by region
f2 = figure('Units','pixels','Position',[100 100 1200 500],'Color','w');
tl = tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(tl); bar(ax1, regs, pctStuck, 'FaceColor',[0.84 0.15 0.16]);
xlabel(ax1,'Region (0 = outside)'); ylabel(ax1,'% of samples flagged stuck'); box(ax1,'on')
title(ax1,'Stuck share by region');
ax2 = nexttile(tl); bb = bar(ax2, regs, [allSp movSp], 'grouped');
bb(1).FaceColor = [0.55 0.55 0.55]; bb(2).FaceColor = [0.00 0.45 0.70];
xlabel(ax2,'Region (0 = outside)'); ylabel(ax2,'Mean 1 Hz Doppler speed (m/s)'); box(ax2,'on')
legend(ax2, {'All samples','Moving samples only'}, 'Location','northeast');
title(ax2,'Mean Doppler speed with and without stuck samples');
saveFig(f2, save_dir, 'StuckDrifters_byRegion');


%% Figure 3: stuck share by tidal phase and run starts/ends
f3 = figure('Units','pixels','Position',[100 100 1200 500],'Color','w');
tl3 = tiledlayout(f3,1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(tl3);
bp = bar(ax1, 1:nPh, [100*phMar(:,2)./max(phMar(:,1),1), 100*phOth(:,2)./max(phOth(:,1),1)], 'grouped');
bp(1).FaceColor = [0.84 0.15 0.16]; bp(2).FaceColor = [0.55 0.55 0.55];
xticks(ax1,1:nPh); xticklabels(ax1,cellstr(phases)); box(ax1,'on')
ylabel(ax1,'% of samples in the phase flagged stuck');
legend(ax1,{'MarshMadness + MarshMayhem','Other periods'},'Location','northeast');
title(ax1,'Stuck share by tidal phase');
ax2 = nexttile(tl3);
bs2 = bar(ax2, 1:nPh, [startCnt endCnt], 'grouped');
bs2(1).FaceColor = [0.84 0.15 0.16]; bs2(2).FaceColor = [0.00 0.45 0.70];
xticks(ax2,1:nPh); xticklabels(ax2,cellstr(phases)); box(ax2,'on')
ylabel(ax2,'Number of stuck runs');
legend(ax2,{'Run starts','Run ends (floats off)'},'Location','northeast');
title(ax2,'Phase when stuck runs start and end');
saveFig(f3, save_dir, 'StuckDrifters_byPhase');

%% Save flags for later use
Meta = rmfield(TR, {'t','x','y','spd','reg','ph','stk'});
Time = {TR.t}';
save(fullfile(save_dir,'StuckDrifters_flags.mat'), 'Flag','Meta','Time','names', ...
    'win_min','radius_m','dop_max','min_run_min','-v7.3');

out.tracks = TT;
out.region = table(regs(:), nSamp, pctStuck, allSp, movSp, ...
    'VariableNames',{'Region','Samples','PctStuck','MeanSpeedAll','MeanSpeedMoving'});
out.sensitivity = Sens;
out.phase = table(phases(:), phAll(:,1), phAll(:,2), phMar(:,1), phMar(:,2), phOth(:,1), phOth(:,2), startCnt, endCnt, ...
    'VariableNames',{'Phase','AllSamples','AllStuck','MarshSamples','MarshStuck','OtherSamples','OtherStuck','RunsStart','RunsEnd'});
fprintf('\nSaved results to %s\n', save_dir);
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


function c = phaseCounts(TR, Flag, phase, sel)
% [samples, stuck samples] in one tidal phase over the selected tracks
c = [0 0];
for q = find(sel)
    m = TR(q).ph == phase;
    c(1) = c(1) + sum(m);
    c(2) = c(2) + sum(m & Flag{q});
end
end


function saveFig(f, save_dir, file_base)
exportgraphics(f, fullfile(save_dir,[file_base '.png']), 'Resolution',300);
savefig(f, fullfile(save_dir,[file_base '.fig']));
end