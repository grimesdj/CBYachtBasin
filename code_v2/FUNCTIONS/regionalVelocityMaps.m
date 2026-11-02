function out = regionalVelocityMaps(cfg)
%% regionalVelocityMaps_v3
% Regional MEAN SPEED and MEAN VELOCITY from drifters, all samples and each
% tidal phase. Same loader and statistics as regionalVelocityMaps_v2 (which is
% left untouched); v3 changes the FIGURES and adds a results table.
%
% Figures saved to THESIS\THESIS_DRAFTING\FIGURES\REGIONAL_ANALYSIS:
%   RegionalVectors_v3_panels_<w>s   one panel per tidal phase (High, Falling, Low, Rising)
%   RegionalVectors_v3_All_<w>s      one large map, all tidal phases
%   RegionalSpeed_v3_bars_<w>s       mean speed by region and phase with 95% intervals
%   RegionalTable_v3_<w>s.csv        results table (also in RegionalStats_v3_<w>s.mat)
%   RegionalVectors_v3_byPeriod_<w>s mean vector per region and observation period (consistency check)
%   RegionalPeriodTable_v3_<w>s.csv  the same, as a table
%
% Settings to record in the decisions log:
%   Uncertainty = 95% bootstrap interval resampling whole deployments
%     (observation period x deployment), only deployments with data in that cell.
%   White arrow = east or north component of the mean has a 95% interval that
%     excludes zero. Grey = cannot rule out zero.
%   Cells with < min_tracks tracks or < min_periods periods are not drawn.
%   Cells with < marginal_tracks tracks are drawn but marked with * (and flagged in the table).
%   exclude_stuck = true drops samples flagged Stuck (Results.Stuck, from flagStuckDrifters and
%     veloAnalysis_v2), widened by 128 s each side (widenStuckMask), the same samples dropped in
%     allVeloAnalysis. Output file names then end in _moving. Run once with false and once with
%     true to compare.
%   Tracks can span several tidal phases, so phase cells are not independent samples.
%   "Rising"/"Falling" are pressure-based phases, not flood/ebb.

close all; rng(1)
if nargin < 1, cfg = []; end %#ok<NASGU>

%% Settings
temporal_res   = 8;      % options: [1 2 4 8 16 32 64 128 256]
min_tracks     = 3;      % cells with fewer tracks are not drawn
min_periods    = 2;      % cells from fewer observation periods are not drawn
marginal_tracks = 5;     % cells with fewer tracks than this are marked *
nBoot          = 2000;
exclude_stuck  = true;   % true = drop samples flagged Stuck (Results.Stuck); false = all samples
velocity_scale = 3000;   % arrow length in metres per 1 m/s
ref_speed      = 0.05;   % reference arrow, m/s
label_dlat     = 0.00045;% region number is placed this far south of the arrow tail (deg)

% % Stationary-deployment mean apparent speed (m/s) by window, for the reference line
% % (from the stationary analysis, Section 3 of the log). Windows [2 4 8 16 32 64 128 256] s.
% stat_windows = [2 4 8 16 32 64 128 256];
% stat_floor   = [1.23 1.18 1.10 0.98 0.80 0.62 0.46 0.34]/100;

regionColors = [0.000 0.447 0.741
                0.850 0.325 0.098
                0.929 0.694 0.125
                0.494 0.184 0.556
                0.466 0.674 0.188
                0.301 0.745 0.933
                0.635 0.078 0.184];

roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', ...
         'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if isempty(k), error('Neither root path contains CB_YachtBasin on this machine.'); end
root_data = roots{k};
fprintf('Using root: %s\n', root_data);

save_dir = fullfile(root_data,'THESIS','THESIS_DRAFTING','FIGURES','REGIONAL_ANALYSIS');
if ~isfolder(save_dir), mkdir(save_dir); end

% Stationary-deployment mean apparent speed (stationary velocity floor) for reference (from stationaryNoiseFloor)
nf = fullfile(root_data,'CB_YachtBasin','StationaryDep','PROCESSED_DATA_V2','DRIFTERS', ...
    'StationaryDep_NoiseFloor.mat');
if ~isfile(nf)
    error('%s not found. Run stationaryNoiseFloor first.', nf);
end
NF = load(nf).Floor;
stat_windows = NF.windows_s;                                % [2 4 8 16 32 64 128 256] s
stat_floor   = mean(NF.Velocity.meanSpd, 1, 'omitnan');     % m/s, mean over the stationary drifters

% disp(stat_floor*100)   % cm/s

sfx = '';
if exclude_stuck, sfx = '_moving'; end
dilate_s = 256;      % mask widened by dilate_s/2 = 128 s each side, same as allVeloAnalysis
nMasked = 0; nTotal = 0;

names   = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};
folders = {'MarshMadness','MarshMayhem','FallFrolic','FallFrolic','NovDep','JuneJamboree','SummerFlood'};
phases  = ["High","Falling","Low","Rising"];   % display order

% Wind (from pullWindData): one table per observation period, times in UTC
WW = load(fullfile(root_data,'CB_YachtBasin','WIND_DATA','Wind_AllPeriods.mat')).Wind;
windNames = {WW.name};
wTotal = zeros(1,numel(names));  wMatched = zeros(1,numel(names));   % for the check printout
tFirst = NaT(1,numel(names),'TimeZone','UTC');  tLast = tFirst;

%% Load every track at the chosen temporal resolution
RD = struct('Observation',{},'Period',{},'Deployment',{},'ID',{},'Cluster',{}, ...
            'V',{},'V_e',{},'V_n',{},'Region',{},'TidePhase',{},'Wu',{},'Wv',{}, ...
            'Vpar',{},'Vperp',{},'Wspd',{}, 'Time', {});

for i = 1:numel(names)

    fn = fullfile(root_data,'CB_YachtBasin',folders{i},'PROCESSED_DATA_V2','DRIFTERS', ...
        [names{i} '_TemporalRes_RMSE.mat']);

    if ~isfile(fn)
        warning('%s: file not found, skipped', names{i});
        continue
    end

    Results = load(fn).Results;
    
    if exclude_stuck && ~isfield(Results, 'Stuck')
        error('%s has no Stuck field. Run main_driver with flagStuckDrifters and veloAnalysis_v2 first.', names{i});
    end
    
    wi = find(strcmp(windNames, names{i}), 1);
    if isempty(wi), error('No wind data for %s. Run pullWindData.', names{i}); end
    Wt = WW(wi).data;
    wx = seconds(Wt.time - Wt.time(1));     % wind times as seconds from the first wind record

    for d = 1:size(Results,1)
        for j = 1:size(Results,2)

            if isempty(Results(d,j).ID), continue, end

            if temporal_res == 1
                V   = Results(d,j).V_1Hz;
                Ve  = Results(d,j).V_e_1Hz;
                Vn  = Results(d,j).V_n_1Hz;
                Reg = Results(d,j).Region;
                Ph  = Results(d,j).TidePhase;
            else
                ri = find([Results(d,j).TemporalRes.Seconds] == temporal_res, 1);
                if isempty(ri), continue, end
                TR  = Results(d,j).TemporalRes(ri);
                V   = TR.V;
                Ve  = TR.V_e;
                Vn  = TR.V_n;
                Reg = TR.Region;
                Ph  = TR.TidePhase;
            end

            if exclude_stuck
                nn = numel(V);
                Sd = widenStuckMask(Results(d,j).Stuck, nn, round(dilate_s/2));
                V = V(:); Ve = Ve(:); Vn = Vn(:);
                V(Sd) = NaN;  Ve(Sd) = NaN;  Vn(Sd) = NaN;
                nMasked = nMasked + sum(Sd);
                nTotal  = nTotal + nn;
            end
            
            % wind at the clock time (UTC) of every sample
            tt = Results(d,j).DateTime(:);
            xq = seconds(tt - Wt.time(1));
            Wu = interp1(wx, Wt.u_to, xq, 'linear', NaN);
            Wv = interp1(wx, Wt.v_to, xq, 'linear', NaN);
            wTotal(i)   = wTotal(i) + numel(Wu);
            wMatched(i) = wMatched(i) + sum(~isnan(Wu));
            tFirst(i) = min([tFirst(i); min(tt)]);  tLast(i) = max([tLast(i); max(tt)]);
            
            % drifter velocity split into along-wind and cross-wind parts
            Wspd  = hypot(Wu, Wv);
            Vpar  = ( Ve(:).*Wu(:) + Vn(:).*Wv(:)) ./ Wspd(:);   % + = moving the way the wind blows
            Vperp = (-Ve(:).*Wv(:) + Vn(:).*Wu(:)) ./ Wspd(:);   % + = to the left of the wind
            Wspd(isnan(Vpar)) = NaN;   % wind speed only where the drifter sample is used
            
            rec = struct('Observation',names{i}, 'Period',i, 'Deployment',d, ...
                'ID',Results(d,j).ID, 'Cluster',i*100 + d, ...
                'V',V(:), 'V_e',Ve(:), 'V_n',Vn(:), ...
                'Region',Reg(:), 'TidePhase',string(Ph(:)), 'Wu',Wu(:), 'Wv',Wv(:), ...
                'Vpar',Vpar(:), 'Vperp',Vperp(:), 'Wspd',Wspd(:), 'Time',tt);

            RD(end+1) = rec; %#ok<AGROW>

        end
    end
end

fprintf('\nLoaded %d tracks at %d s resolution.\n', numel(RD), temporal_res);
if exclude_stuck
    fprintf('Samples excluded as stuck: %d of %d (%.1f%%)\n', nMasked, nTotal, 100*nMasked/max(nTotal,1));
end

fprintf('Along-wind check: mean Vpar = %.3f m/s, mean wind speed = %.2f m/s\n', ...
    mean(vertcat(RD.Vpar),'omitnan'), mean(vertcat(RD.Wspd),'omitnan'));

fprintf('\n=== Wind matched to drifter samples (UTC) ===\n');
for i = 1:numel(names)
    if wTotal(i) == 0, continue, end
    fprintf('%-20s %s to %s   matched %.1f%% of %d samples\n', names{i}, ...
        char(tFirst(i),'yyyy-MM-dd HH:mm'), char(tLast(i),'yyyy-MM-dd HH:mm'), ...
        100*wMatched(i)/wTotal(i), wTotal(i));
end

%% Summarize the wind during each drifter deployment day, each IOP (ONLY WHEN DRIFTERS IN WATER)
fprintf('\n=== Wind while drifters were out, by IOP (UTC) ===\n');
fprintf('%-20s %10s %10s %11s %10s\n', 'period', 'mean m/s', 'max m/s', 'steadiness', 'from deg');
for p = unique([RD.Period])
    k  = [RD.Period] == p;                       % all tracks in this IOP
    tt = vertcat(RD(k).Time);
    Wu = vertcat(RD(k).Wu);
    Wv = vertcat(RD(k).Wv);

    [~, iu] = unique(tt);                        % each clock time counted once
    Wu = Wu(iu);  Wv = Wv(iu);

    spd   = hypot(Wu, Wv);
    meanU = mean(Wu, 'omitnan');
    meanV = mean(Wv, 'omitnan');
    fromDeg = mod(atan2d(meanU, meanV) + 180, 360);

    fprintf('%-20s %10.1f %10.1f %11.2f %10.0f\n', RD(find(k,1)).Observation, ...
        mean(spd,'omitnan'), max(spd), hypot(meanU,meanV)/mean(spd,'omitnan'), fromDeg);
end
%% Statistics: all phases, then each phase
Rg = load('regions.mat');
regions = Rg.regions;
nR = numel(regions);

groups(1).label = 'all tidal phases';
groups(1).tag   = 'All';
for r = 1:nR
    groups(1).S(r) = cellStats(RD, r, "", nBoot);
end
for p = 1:numel(phases)
    groups(p+1).label = [char(phases(p)) ' phase'];
    groups(p+1).tag   = char(phases(p));
    for r = 1:nR
        groups(p+1).S(r) = cellStats(RD, r, phases(p), nBoot);
    end
end

for g = 1:numel(groups)
    fprintf('\n=== %s (%d s) ===\n', groups(g).label, temporal_res);
    printTable(groups(g).S);
end

%% Results table
T = buildTable(groups, min_tracks, min_periods, marginal_tracks);
csvfile = fullfile(save_dir, sprintf('RegionalTable_v3_%ds%s.csv', temporal_res, sfx));
writetable(T, csvfile);
fprintf('\nSaved %s\n', csvfile);

%% Figure 1: one panel per tidal phase
lims = mapLimits(regions);
f1 = figure('Units','pixels','Position',[50 50 1700 900],'Color','w');
tl = tiledlayout(f1, 1, numel(phases), 'TileSpacing','compact', 'Padding','compact');
for p = 1:numel(phases)
    ax = geoaxes(tl);
    ax.Layout.Tile = p;
    drawVectorMap(ax, regions, groups(p+1).S, regionColors, lims, ...
        min_tracks, min_periods, marginal_tracks, velocity_scale, ref_speed, label_dlat);
    title(ax, sprintf('%s', groups(p+1).label), 'FontSize',14);
end
title(tl, sprintf('Mean drifter velocity by region and tidal phase (%d s window)', temporal_res), ...
    'FontSize',15, 'FontWeight','bold');
saveFig(f1, save_dir, sprintf('RegionalVectors_v3_panels_%ds%s', temporal_res, sfx));

%% Figure 2: one large map, all phases
f2 = figure('Units','pixels','Position',[100 100 1000 900],'Color','w');
ax = geoaxes(f2);
[hreg, lab] = drawVectorMap(ax, regions, groups(1).S, regionColors, lims, ...
    min_tracks, min_periods, marginal_tracks, velocity_scale, ref_speed, label_dlat);
legend(ax, hreg, lab, 'Location','southoutside', 'NumColumns',1, 'FontSize',10);
title(ax, sprintf('Mean drifter velocity by region, all tidal phases (%d s window)', temporal_res), 'FontSize',14);
saveFig(f2, save_dir, sprintf('RegionalVectors_v3_All_%ds%s', temporal_res, sfx));

%% Figure 3: mean speed bars with 95% intervals
f3 = figure('Units','pixels','Position',[100 100 1300 600],'Color','w');
ax = axes(f3); hold(ax,'on'); box(ax,'on')
M = nan(nR,5); Lo = M; Up = M;
for g = 1:5
    for r = 1:nR
        s = groups(g).S(r);
        if s.NTracks >= min_tracks && s.NPeriods >= min_periods
            M(r,g) = s.MeanSpeed;  Lo(r,g) = s.SpeedCI(1);  Up(r,g) = s.SpeedCI(2);
        end
    end
end
barColors = [0.55 0.55 0.55; 0.84 0.37 0.00; 0.80 0.47 0.65; 0.00 0.45 0.70; 0.00 0.62 0.45];
b = bar(ax, M, 'grouped');
for g = 1:5
    b(g).FaceColor = barColors(g,:);
    errorbar(ax, b(g).XEndPoints, M(:,g), M(:,g)-Lo(:,g), Up(:,g)-M(:,g), ...
        'k', 'LineStyle','none', 'CapSize',4, 'LineWidth',1);
end
iw = find(stat_windows == temporal_res, 1);
hleg = b;  legTxt = {'All phases','High','Falling','Low','Rising'};
if ~isempty(iw)
    fl = stat_floor(iw);
    hl = yline(ax, fl, '--k', 'LineWidth',1.2);
    hleg = [b, hl];
    legTxt{end+1} = sprintf('Stationary mean apparent speed (%.1f cm/s)', fl*100);
end
xlabel(ax,'Region'); ylabel(ax,'Mean speed (m/s)');
xticks(ax,1:nR); xticklabels(ax, string(1:nR));
legend(ax, hleg, legTxt, 'Location','northeastoutside');
title(ax, sprintf('Mean drifter speed by region and tidal phase (%d s window)', temporal_res), 'FontSize',14);
set(ax,'FontSize',12);
saveFig(f3, save_dir, sprintf('RegionalSpeed_v3_bars_%ds%s', temporal_res, sfx));

%% Diagnostic: mean vector by region and observation period (all tidal phases)
% Question: is the steady northward mean vector in regions 1, 3, 5 and 6 present in
% every observation period, or driven by one or two? Cells with fewer than min_tracks
% tracks are left blank. No intervals here (cells are small); read it as a consistency check.
nP = numel(names);
Pn = zeros(nR,nP);  Pd = zeros(nR,nP);
Pve = nan(nR,nP);   Pvn = nan(nR,nP);
for i = 1:nP
    idx = find([RD.Period] == i);
    for r = 1:nR
        nt = 0; sN = 0; sVe = 0; sVn = 0; cl = [];
        for t = idx
            m = RD(t).Region == r & ~isnan(RD(t).V) & ~isnan(RD(t).V_e) & ~isnan(RD(t).V_n);
            if any(m)
                nt  = nt + 1;
                sN  = sN + sum(m);
                sVe = sVe + sum(RD(t).V_e(m));
                sVn = sVn + sum(RD(t).V_n(m));
                cl(end+1) = RD(t).Cluster; %#ok<AGROW>
            end
        end
        Pn(r,i) = nt;
        Pd(r,i) = numel(unique(cl));
        if sN > 0
            Pve(r,i) = sVe/sN;
            Pvn(r,i) = sVn/sN;
        end
    end
end
Pmag = hypot(Pve, Pvn);
Pdir = mod(atan2d(Pve, Pvn), 360);
Pshow = Pn >= min_tracks;

fprintf('\n=== Mean vector by region and observation period (all phases, %d s) ===\n', temporal_res);
fprintf('Each cell: tracks | |mean V| m/s | direction toward (deg from north). -- = fewer than %d tracks.\n', min_tracks);
fprintf('%-7s', 'Region');
for i = 1:nP, fprintf(' %-22s', strrep(names{i},'_',' ')); end
fprintf('\n');
for r = 1:nR
    fprintf('%-7d', r);
    for i = 1:nP
        if Pshow(r,i)
            fprintf(' %3d | %.3f | %3.0f       ', Pn(r,i), Pmag(r,i), Pdir(r,i));
        else
            fprintf(' %-22s', '--');
        end
    end
    fprintf('\n');
end

[RR, PP] = ndgrid(1:nR, 1:nP);
PeriodName = strrep(string(names(PP(:))), '_', ' ');
PeriodTable = table(RR(:), PeriodName(:), Pn(:), Pd(:), Pve(:), Pvn(:), Pmag(:), Pdir(:), Pshow(:), ...
    'VariableNames', {'Region','Period','Tracks','Deps','Ve','Vn','VecMag','Direction','Drawn'});
pcsv = fullfile(save_dir, sprintf('RegionalPeriodTable_v3_%ds%s.csv', temporal_res, sfx));
writetable(PeriodTable, pcsv);
fprintf('Saved %s\n', pcsv);

%% Wind diagnostic: along-wind and cross-wind drift by region and observation period
% % Same layout idea as the period diagnostic above. Moving samples only. No intervals.
Wsp = nan(nR,nP);  Wpar = nan(nR,nP);  Wperp = nan(nR,nP);
for i = 1:nP
    idx = find([RD.Period] == i);
    for r = 1:nR
        sN = 0; sW = 0; sPar = 0; sPerp = 0;
        for t = idx
            m = RD(t).Region == r & ~isnan(RD(t).Vpar) & ~isnan(RD(t).Vperp) & ~isnan(RD(t).Wspd);
            sN    = sN + sum(m);
            sW    = sW + sum(RD(t).Wspd(m));
            sPar  = sPar + sum(RD(t).Vpar(m));
            sPerp = sPerp + sum(RD(t).Vperp(m));
        end
        if sN > 0
            Wsp(r,i) = sW/sN;  Wpar(r,i) = sPar/sN;  Wperp(r,i) = sPerp/sN;
        end
    end
end

fprintf('\n=== Along-wind and cross-wind drift by region and period (%d s, moving samples) ===\n', temporal_res);
fprintf('along = + moving the way the wind blows; cross = + to the left of the wind; -- = fewer than %d tracks.\n', min_tracks);
for i = 1:nP
    fprintf('\n%s\n', strrep(names{i},'_',' '));
    fprintf('%-7s %7s %9s %11s %11s %11s\n', 'Region','tracks','wind m/s','along cm/s','cross cm/s','along/wind');
    for r = 1:nR
        if Pshow(r,i)
            fprintf('%-7d %7d %9.2f %+11.2f %+11.2f %10.1f%%\n', r, Pn(r,i), Wsp(r,i), ...
                100*Wpar(r,i), 100*Wperp(r,i), 100*Wpar(r,i)/Wsp(r,i));
        else
            fprintf('%-7d %7s\n', r, '--');
        end
    end
end

fprintf('\n=== Pooled along-wind drift by period (%d s, moving samples, all regions) ===\n', temporal_res);
fprintf('%-20s %7s %9s %11s %11s %11s\n', 'period', 'tracks', 'wind m/s', 'along cm/s', 'cross cm/s', 'along/wind');
for i = 1:nP
    idx = find([RD.Period] == i);
    w = vertcat(RD(idx).Wspd);
    a = vertcat(RD(idx).Vpar);
    c = vertcat(RD(idx).Vperp);
    ok = ~isnan(a) & ~isnan(c) & ~isnan(w);
    if ~any(ok), continue, end
    nTr = sum(arrayfun(@(t) any(~isnan(t.Vpar)), RD(idx)));
    fprintf('%-20s %7d %9.2f %+11.2f %+11.2f %10.1f%%\n', names{i}, nTr, ...
        mean(w(ok)), 100*mean(a(ok)), 100*mean(c(ok)), 100*mean(a(ok))/mean(w(ok)));
end

% f4 = figure('Units','pixels','Position',[100 100 1300 650],'Color','w');
% ax = axes(f4); hold(ax,'on'); box(ax,'on')
% diag_ref    = 0.10;                    % reference arrow for this figure, m/s
% arrow_scale = 0.45/diag_ref;           % axis units per m/s (reference arrow = 0.45 units)
% for r = 1:nR
%     for i = 1:nP
%         if ~Pshow(r,i), continue, end
%         quiver(ax, i, r, Pve(r,i)*arrow_scale, Pvn(r,i)*arrow_scale, 0, ...
%             'Color', regionColors(r,:), 'LineWidth',2, 'MaxHeadSize',0.6);
%         plot(ax, i, r, '.', 'Color', regionColors(r,:), 'MarkerSize',10);
%         text(ax, i - 0.42, r + 0.30, sprintf('n=%d', Pn(r,i)), 'FontSize',9, 'Color',[0.35 0.35 0.35]);
%     end
% end
% quiver(ax, nP + 0.30, nR + 0.7, 0.45, 0, 0, 'k', 'LineWidth',2, 'MaxHeadSize',0.6);
% text(ax, nP + 0.30, nR + 0.95, sprintf('%.2f m/s', diag_ref), 'FontSize',10);
% set(ax, 'YDir','reverse', 'FontSize',12, 'XLim',[0.4 nP+1], 'YLim',[0.4 nR+1.2]);
% xticks(ax,1:nP); xticklabels(ax, strrep(names,'_',' '));
% yticks(ax,1:nR); ylabel(ax,'Region');
% title(ax, sprintf('Mean drifter velocity by region and observation period, all tidal phases (%d s window); north is up', temporal_res), 'FontSize',13);
% saveFig(f4, save_dir, sprintf('RegionalVectors_v3_byPeriod_%ds%s', temporal_res, sfx));

%% Save
save(fullfile(save_dir, sprintf('RegionalStats_v3_%ds%s.mat', temporal_res, sfx)), ...
    'groups','T','PeriodTable','temporal_res','min_tracks','min_periods','marginal_tracks','nBoot');
fprintf('Saved RegionalStats_v3_%ds%s.mat to %s\n', temporal_res, sfx, save_dir);

out.groups = groups;
out.table  = T;
out.periodTable = PeriodTable;
end


%% ========================= Local functions =========================

function S = cellStats(RD, r, phase, nBoot)
% Statistics for one region (r) and, if phase ~= "", one tidal phase.
% Bootstrap resamples whole deployments (the Cluster field).

S = struct('N',0,'NTracks',0,'NPeriods',0,'NClusters',0, ...
           'MeanSpeed',NaN,'StdSpeed',NaN,'Ve',NaN,'Vn',NaN,'VecMag',NaN, ...
           'SpeedCI',[NaN NaN],'VeCI',[NaN NaN],'VnCI',[NaN NaN],'Differs',false);

nT  = numel(RD);
n   = zeros(nT,1);  sV = n;  sV2 = n;  sVe = n;  sVn = n;

for t = 1:nT
    m = RD(t).Region == r;
    if phase ~= ""
        m = m & RD(t).TidePhase == phase;
    end
    m = m & ~isnan(RD(t).V) & ~isnan(RD(t).V_e) & ~isnan(RD(t).V_n);

    n(t)   = sum(m);
    sV(t)  = sum(RD(t).V(m));
    sV2(t) = sum(RD(t).V(m).^2);
    sVe(t) = sum(RD(t).V_e(m));
    sVn(t) = sum(RD(t).V_n(m));
end

keep = n > 0;
if ~any(keep), return, end

S.N        = sum(n);
S.NTracks  = sum(keep);
S.NPeriods = numel(unique([RD(keep).Period]));

S.MeanSpeed = sum(sV)/S.N;
S.StdSpeed  = sqrt(max(sum(sV2)/S.N - S.MeanSpeed^2, 0));
S.Ve        = sum(sVe)/S.N;
S.Vn        = sum(sVn)/S.N;
S.VecMag    = hypot(S.Ve, S.Vn);

cl = [RD.Cluster]';
[~,~,g] = unique(cl(keep));
cN  = accumarray(g, n(keep));
cV  = accumarray(g, sV(keep));
cVe = accumarray(g, sVe(keep));
cVn = accumarray(g, sVn(keep));

G = numel(cN);
S.NClusters = G;

if G >= 2
    bs = nan(nBoot,3);
    for b = 1:nBoot
        ix = randi(G,G,1);
        bs(b,:) = [sum(cV(ix)) sum(cVe(ix)) sum(cVn(ix))] / sum(cN(ix));
    end
    S.SpeedCI = ci95(bs(:,1));
    S.VeCI    = ci95(bs(:,2));
    S.VnCI    = ci95(bs(:,3));
    S.Differs = S.VeCI(1) > 0 || S.VeCI(2) < 0 || S.VnCI(1) > 0 || S.VnCI(2) < 0;
end

end


function ci = ci95(x)
x = sort(x(~isnan(x)));
n = numel(x);
ci = x([max(1,round(0.025*n)), min(n,round(0.975*n))])';
end


function printTable(S)
fprintf('%-7s %7s %8s %6s  %-22s %8s %8s %9s  %s\n', 'Region','tracks','periods','deps', ...
    'mean speed [95% CI]','Ve','Vn','|mean V|','mean vector');
yn = {'cannot rule out zero','differs from zero'};
for r = 1:numel(S)
    if S(r).N == 0
        fprintf('%-7d no data\n', r);
        continue
    end
    sp = sprintf('%.3f [%.3f-%.3f]', S(r).MeanSpeed, S(r).SpeedCI(1), S(r).SpeedCI(2));
    fprintf('%-7d %7d %8d %6d  %-22s %+8.3f %+8.3f %9.3f  %s\n', r, ...
        S(r).NTracks, S(r).NPeriods, S(r).NClusters, sp, ...
        S(r).Ve, S(r).Vn, S(r).VecMag, yn{S(r).Differs + 1});
end
end


function T = buildTable(groups, min_tracks, min_periods, marginal_tracks)
% One row per region x phase. Direction = compass bearing the mean vector points TOWARD
% (degrees clockwise from north). Steadiness = |mean V| / mean speed (1 = all tracks
% moving the same way at the same speed; small = no consistent direction).
Phase = strings(0,1); Region = []; Tracks = []; Periods = []; Deps = [];
MeanSpeed = []; SpeedLo = []; SpeedHi = []; Ve = []; Vn = []; VecMag = [];
Direction = []; Steadiness = []; DiffersFromZero = false(0,1); Status = strings(0,1);

for g = 1:numel(groups)
    for r = 1:numel(groups(g).S)
        s = groups(g).S(r);
        Phase(end+1,1) = string(groups(g).label); %#ok<AGROW>
        Region(end+1,1) = r; %#ok<AGROW>
        Tracks(end+1,1) = s.NTracks; %#ok<AGROW>
        Periods(end+1,1) = s.NPeriods; %#ok<AGROW>
        Deps(end+1,1) = s.NClusters; %#ok<AGROW>
        MeanSpeed(end+1,1) = s.MeanSpeed; %#ok<AGROW>
        SpeedLo(end+1,1) = s.SpeedCI(1); %#ok<AGROW>
        SpeedHi(end+1,1) = s.SpeedCI(2); %#ok<AGROW>
        Ve(end+1,1) = s.Ve; %#ok<AGROW>
        Vn(end+1,1) = s.Vn; %#ok<AGROW>
        VecMag(end+1,1) = s.VecMag; %#ok<AGROW>
        Direction(end+1,1) = mod(atan2d(s.Ve, s.Vn), 360); %#ok<AGROW>
        Steadiness(end+1,1) = s.VecMag / s.MeanSpeed; %#ok<AGROW>
        DiffersFromZero(end+1,1) = s.Differs; %#ok<AGROW>
        if s.NTracks < min_tracks || s.NPeriods < min_periods
            Status(end+1,1) = "too few data (not drawn)"; %#ok<AGROW>
        elseif s.NTracks < marginal_tracks
            Status(end+1,1) = "marginal (<" + marginal_tracks + " tracks)"; %#ok<AGROW>
        else
            Status(end+1,1) = "ok"; %#ok<AGROW>
        end
    end
end

T = table(Phase, Region, Tracks, Periods, Deps, MeanSpeed, SpeedLo, SpeedHi, ...
    Ve, Vn, VecMag, Direction, Steadiness, DiffersFromZero, Status);
end


function lims = mapLimits(regions)
allLat = cell2mat(arrayfun(@(x) x.Lat(:), regions(:), 'UniformOutput', false));
allLon = cell2mat(arrayfun(@(x) x.Lon(:), regions(:), 'UniformOutput', false));
lims.lat = [min(allLat,[],'omitnan') - 0.0030, max(allLat,[],'omitnan') + 0.0006];
lims.lon = [min(allLon,[],'omitnan') - 0.0015, max(allLon,[],'omitnan') + 0.0015];
lims.latRef = min(allLat,[],'omitnan') - 0.0007;
lims.lonRef = min(allLon,[],'omitnan') + 0.0002;
end


function [h, lab] = drawVectorMap(ax, regions, S, colors, lims, min_tracks, min_periods, ...
                                  marginal_tracks, vscale, ref_speed, label_dlat)

hold(ax,'on')
geobasemap(ax,'satellite')

nR = numel(regions);
h   = gobjects(nR,1);
lab = strings(nR,1);

for r = 1:nR
    lat = regions(r).Lat(:);
    lon = regions(r).Lon(:);
    c   = colors(mod(r-1,size(colors,1))+1,:);

    h(r) = geoplot(ax, lat, lon, 'Color',c, 'LineWidth',2);

    latc = mean(lat(1:end-1),'omitnan');
    lonc = mean(lon(1:end-1),'omitnan');

    ok = S(r).NTracks >= min_tracks && S(r).NPeriods >= min_periods;
    marginal = ok && S(r).NTracks < marginal_tracks;

    txt = sprintf('%d', r);
    if marginal, txt = [txt '*']; end %#ok<AGROW>

    if ok
        if S(r).Differs, col = 'w'; else, col = [0.65 0.65 0.65]; end
        drawArrow(ax, latc, lonc, S(r).Ve, S(r).Vn, vscale, col);
        lab(r) = sprintf('Region %d%s: |mean V| = %.3f m/s  (%d tracks, %d periods)', ...
            r, repmat('*',1,marginal), S(r).VecMag, S(r).NTracks, S(r).NPeriods);
    else
        lab(r) = sprintf('Region %d: too few data  (%d tracks, %d periods)', ...
            r, S(r).NTracks, S(r).NPeriods);
    end

    text(ax, latc - label_dlat, lonc, txt, 'FontSize',12, 'FontWeight','bold', ...
        'Color','w', 'BackgroundColor',[0.15 0.15 0.15], 'Margin',1.5, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

drawArrow(ax, lims.latRef, lims.lonRef, ref_speed, 0, vscale, 'w');
text(ax, lims.latRef - 0.0004, lims.lonRef, sprintf('%.2f m/s', ref_speed), ...
    'Color','w', 'FontWeight','bold', 'BackgroundColor',[0.15 0.15 0.15], 'Margin',1.5, ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');

geolimits(ax, lims.lat, lims.lon);
end


function drawArrow(ax, lat0, lon0, Ve, Vn, vscale, col)
% Arrow from (lat0,lon0); length = speed * vscale metres. Colored col w a black outline.

mlat = 111320;
mlon = 111320*cosd(lat0);

dx = Ve*vscale;
dy = Vn*vscale;
L  = hypot(dx,dy);

if ~(L > 0), return, end

lat1 = lat0 + dy/mlat;
lon1 = lon0 + dx/mlon;

% shaft: black underlay then white line
geoplot(ax, [lat0 lat1], [lon0 lon1], 'Color','k', 'LineStyle','-', 'LineWidth',5);
geoplot(ax, [lat0 lat1], [lon0 lon1], 'Color',col, 'LineStyle','-', 'LineWidth',2.5);

ux = dx/L;  uy = dy/L;
px = -uy;   py = ux;

hl = 0.25*L;  hl = min(hl, 40);   % head length, metres
hw = 0.5*hl;

x1 = dx - hl*ux + hw*px;   y1 = dy - hl*uy + hw*py;
x2 = dx - hl*ux - hw*px;   y2 = dy - hl*uy - hw*py;

hlat = [lat0 + y1/mlat, lat1, lat0 + y2/mlat];
hlon = [lon0 + x1/mlon, lon1, lon0 + x2/mlon];
geoplot(ax, hlat, hlon, 'Color','k', 'LineStyle','-', 'LineWidth',5);
geoplot(ax, hlat, hlon, 'Color',col, 'LineStyle','-', 'LineWidth',2.5);
end


function saveFig(f, save_dir, file_base)
exportgraphics(f, fullfile(save_dir,[file_base '.png']), 'Resolution',300);
savefig(f, fullfile(save_dir,[file_base '.fig']));
end