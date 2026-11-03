function out = stationaryNoiseFloor(cfg)
%% stationaryNoiseFloor
% Velocity "floor" and position error from NON-MOVING drifters.
% These drifters' true velocity is zero and their true position is fixed, so any
% apparent motion is error. Reads StationaryDep level 1 (1 Hz) directly, prints
% tables for the decisions log / methods, and saves the numbers to
% StationaryDep_NoiseFloor.mat in the same folder.
%
% Velocity for window L seconds = (position(t+L) - position(t)) / L, which is the
% same displacement-over-window definition used in veloAnalysis.
%
% Usage (obperiod = 'ALL' in main_driver):  out = stationaryNoiseFloor(cfg)
% The stationary deployment has its own level 1 file, so cfg is only accepted for
% consistency with the other analysis functions; paths are found automatically.
%
% Caveat: the floor describes THESE drifters (one deployment, one site, one time
% of day, and however they were held). If they floated or swung on a line, part
% of the "floor" is real motion.

% Advisor requests (both done below):
%  1. RMS of each drifter's positions about its own mean (precision) AND the
%     absolute error of its mean position against the RTK truth (accuracy):
%     see "Position error about each drifter's own mean" and the RTK section.
%  2. Decorrelation time of the VELOCITY errors (position-derived and Doppler),
%     compared with the position-error decorrelation: see "Decorrelation of the
%     VELOCITY errors".

close all
if nargin < 1, cfg = []; end %#ok<NASGU>

%% Settings
windows    = [2 4 8 16 32 64 128 256];   % s, same windows as veloAnalysis
minSamples = 600;                        % skip tracks shorter than this (s)
maxLagAcf  = 900;                        % s, for position-error autocorrelation

roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', ...
         'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if isempty(k), error('Neither root path contains CB_YachtBasin on this machine.'); end

fn = fullfile(roots{k},'CB_YachtBasin','StationaryDep','PROCESSED_DATA_V2','DRIFTERS', ...
    'StationaryDep_drifters_level_1.mat');
out_file = fullfile(fileparts(fn),'StationaryDep_NoiseFloor.mat');

L1 = load(fn).drifters_level_1;

%% Collect the usable tracks
T = struct('ID',{},'Deployment',{},'E',{},'N',{},'Dop',{});

for d = 1:size(L1,1)
    for j = 1:size(L1,2)

        D = L1(d,j);
        if isempty(D.East) || isempty(D.North), continue, end

        E = D.East(:);
        N = D.North(:);

        if numel(E) < minSamples
            fprintf('Skipping deployment %d, column %d: only %d samples\n', d, j, numel(E));
            continue
        end

        if isfield(D,'e_time') && ~isempty(D.e_time) && abs(median(diff(D.e_time(:))) - 1) > 1e-6
            warning('Deployment %d, column %d: time step is not 1 s', d, j);
        end

        if isfield(D,'ID') && ~isempty(D.ID), id = D.ID; else, id = j; end

        % GPS Doppler speed, knots -> m/s (level 1 Doppler is interpolated to 1 Hz)
        T(end+1) = struct('ID',id,'Deployment',d,'E',E,'N',N,'Dop',D.Speed(:)*0.51444); %#ok<SAGROW>

    end
end

nT = numel(T);
if nT == 0, error('No usable stationary tracks found.'); end
nW = numel(windows);

fprintf('Stationary tracks used: %d (deployment rows: %s)\n', nT, mat2str(unique([T.Deployment])));

%% Position error about each drifter's own mean position
P = struct('ID',[T.ID]', 'hours',nan(nT,1), 'sigE',nan(nT,1), 'sigN',nan(nT,1), ...
           'rms2D',nan(nT,1), 'cep50',nan(nT,1), 'r95',nan(nT,1), 'trend_mph',nan(nT,1));

for q = 1:nT
    E  = T(q).E;
    N  = T(q).N;
    dE = E - mean(E,'omitnan');
    dN = N - mean(N,'omitnan');
    rr = hypot(dE,dN);
    t  = (0:numel(E)-1)'/3600;
    ok = ~isnan(dE) & ~isnan(dN);
    pE = polyfit(t(ok), dE(ok), 1);
    pN = polyfit(t(ok), dN(ok), 1);

    P.hours(q)     = numel(E)/3600;
    P.sigE(q)      = std(E,'omitnan');
    P.sigN(q)      = std(N,'omitnan');
    P.rms2D(q)     = sqrt(mean(rr.^2,'omitnan'));
    P.cep50(q)     = pct(rr,50);
    P.r95(q)       = pct(rr,95);
    P.trend_mph(q) = hypot(pE(1),pN(1));
end

%% Apparent velocity (the floor) for each window length
F.sigU    = nan(nT,nW);
F.sigV    = nan(nT,nW);
F.biasU   = nan(nT,nW);
F.biasV   = nan(nT,nW);
F.meanSpd = nan(nT,nW);
F.rmsSpd  = nan(nT,nW);
F.p95Spd  = nan(nT,nW);
poolSpd   = cell(1,nW);

maxHalf = max(windows)/2;      % 128 s trimmed from each end, same as veloAnalysis

for q = 1:nT
    for w = 1:nW
        L = windows(w);
        h = L/2;
        c = (1+maxHalf):(numel(T(q).E)-maxHalf);   % same centre samples for every window
        u = (T(q).E(c+h) - T(q).E(c-h)) / L;
        v = (T(q).N(c+h) - T(q).N(c-h)) / L;
        s = hypot(u,v);

        F.sigU(q,w)    = std(u,'omitnan');
        F.sigV(q,w)    = std(v,'omitnan');
        F.biasU(q,w)   = mean(u,'omitnan');
        F.biasV(q,w)   = mean(v,'omitnan');
        F.meanSpd(q,w) = mean(s,'omitnan');
        F.rmsSpd(q,w)  = sqrt(mean(s.^2,'omitnan'));
        F.p95Spd(q,w)  = pct(s,95);
        poolSpd{w}     = [poolSpd{w}; s(~isnan(s))];
    end
end

p95pool = cellfun(@(x) pct(x,95), poolSpd);

%% Position-error autocorrelation (how long errors stay correlated)
% Position errors stay correlated for a long time; the velocity errors are
% examined separately in the next section.
maxLag = min(maxLagAcf, floor(min(arrayfun(@(s) numel(s.E), T))/3));
rE = zeros(maxLag+1,nT);
rN = rE;
efold = nan(nT,1);

for q = 1:nT
    rE(:,q) = acf(T(q).E, maxLag);
    rN(:,q) = acf(T(q).N, maxLag);
    r  = 0.5*(rE(:,q) + rN(:,q));
    k1 = find(r < exp(-1), 1);
    if ~isempty(k1), efold(q) = k1 - 1; end      % seconds
end

%% Decorrelation of the VELOCITY errors (position-derived and Doppler)
% Same idea as above, but for velocity. Compared for NON-MOVING drifters:
%   1 s      position velocity from consecutive 1 Hz positions
%   8 s      position velocity over the 8 s window used in the maps (these windows
%            overlap, so they are correlated for up to 8 s from the method alone)
%   Doppler  GPS Doppler speed (level 1 interpolates it to 1 Hz)
% Each is reported as the first lag where the autocorrelation drops below 1/e and
% the first lag where it reaches zero.
velWin = 8;
rV1  = zeros(maxLag+1,nT);  rV8 = rV1;  rDop = rV1;
tV1  = nan(nT,2);  tV8 = tV1;  tDop = tV1;  tPos = tV1;   % [e-fold, first zero] (s)

for q = 1:nT
    u1 = diff(T(q).E);   v1 = diff(T(q).N);
    u8 = (T(q).E(1+velWin:end) - T(q).E(1:end-velWin))/velWin;
    v8 = (T(q).N(1+velWin:end) - T(q).N(1:end-velWin))/velWin;

    rV1(:,q)  = 0.5*(acf(u1,maxLag) + acf(v1,maxLag));
    rV8(:,q)  = 0.5*(acf(u8,maxLag) + acf(v8,maxLag));
    rDop(:,q) = acf(T(q).Dop, maxLag);

    tV1(q,:)  = lagTimes(rV1(:,q));
    tV8(q,:)  = lagTimes(rV8(:,q));
    tDop(q,:) = lagTimes(rDop(:,q));
    tPos(q,:) = lagTimes(0.5*(rE(:,q) + rN(:,q)));
end

%% Doppler speed at rest
dop = vertcat(T.Dop);
dop = dop(~isnan(dop));

%% Print results
fprintf('\n=== Apparent velocity of NON-MOVING drifters (true velocity = 0), cm/s ===\n');
fprintf('%-8s %8s %8s  %-22s %10s %12s\n','window','sigma_u','sigma_v','mean speed [range]','RMS speed','95th pct');
for w = 1:nW
    ms  = F.meanSpd(:,w)*100;
    txt = sprintf('%.2f [%.2f-%.2f]', mean(ms,'omitnan'), min(ms), max(ms));
    fprintf('%-8s %8.2f %8.2f  %-22s %10.2f %12.2f\n', sprintf('%d s',windows(w)), ...
        100*mean(F.sigU(:,w),'omitnan'), 100*mean(F.sigV(:,w),'omitnan'), txt, ...
        100*mean(F.rmsSpd(:,w),'omitnan'), 100*p95pool(w));
end
fprintf('(sigma_u, sigma_v = std of east/north velocity; 95th pct pools all drifters)\n');

fprintf('\n=== Position error about each drifter''s own mean position (m) ===\n');
fprintf('%-6s %7s %8s %8s %8s %8s %11s %10s\n','ID','hours','sigma_E','sigma_N','RMS 2D','CEP50','95% radius','trend m/h');
for q = 1:nT
    fprintf('%-6d %7.2f %8.2f %8.2f %8.2f %8.2f %11.2f %10.2f\n', P.ID(q), P.hours(q), ...
        P.sigE(q), P.sigN(q), P.rms2D(q), P.cep50(q), P.r95(q), P.trend_mph(q));
end
fprintf('%-6s %7s %8.2f %8.2f %8.2f %8.2f %11.2f %10.2f\n', 'mean', '', ...
    mean(P.sigE), mean(P.sigN), mean(P.rms2D), mean(P.cep50), mean(P.r95), mean(P.trend_mph));

fprintf('\n=== GPS Doppler speed of non-moving drifters (m/s) ===\n');
fprintf('mean %.3f | median %.3f | 95th pct %.3f | share below 0.01 m/s: %.0f%%\n', ...
    mean(dop), pct(dop,50), pct(dop,95), 100*mean(dop < 0.01));

fprintf('\n=== Position-error decorrelation (e-folding time of the autocorrelation) ===\n');
for q = 1:nT
    if isnan(efold(q))
        fprintf('ID %-4d longer than %.0f min\n', P.ID(q), maxLag/60);
    else
        fprintf('ID %-4d %.1f min\n', P.ID(q), efold(q)/60);
    end
end
fprintf('(roughly one independent position sample per 2 x e-folding time)\n');

fprintf('\n=== Decorrelation time of VELOCITY errors (mean over drifters, seconds) ===\n');
fprintf('%-34s %14s %16s\n','quantity','e-fold (s)','first zero (s)');
fprintf('%-34s %14s %16s\n','Position velocity, 1 s',        fmtT(tV1(:,1),maxLag),  fmtT(tV1(:,2),maxLag));
fprintf('%-34s %14s %16s\n','Position velocity, 8 s window',  fmtT(tV8(:,1),maxLag),  fmtT(tV8(:,2),maxLag));
fprintf('%-34s %14s %16s\n','GPS Doppler speed',              fmtT(tDop(:,1),maxLag), fmtT(tDop(:,2),maxLag));
fprintf('%-34s %14s %16s\n','Position error (for comparison)',fmtT(tPos(:,1),maxLag), fmtT(tPos(:,2),maxLag));
fprintf('("> %d" = did not drop that far within the lag window for at least one drifter)\n', maxLag);
fprintf('(1 s differences of white position noise are anti-correlated, so their e-fold is 1 s;\n');
fprintf(' the 8 s window velocities are correlated up to ~8 s just because the windows overlap)\n');

%% Compare with the RTK survey (true positions)
% Survey file: one row per RTK shot (name, lat, lon, height), about 5 shots per
% drifter. Names like d6, d6-2, D10-3 mean drifter 6, 10. Rows that are not
% d<number> (e.g. the reference-station row) are ignored.
% NOTE: RTK networks usually give NAD83(2011) while GPS reports WGS84; these
% differ by roughly a metre in the eastern US, so expect a COMMON offset that is
% the same for every drifter. It does not affect velocities (differences cancel it).

cols = lines(nT);
Sv = [];
% Look in the known RAW folder first, then next to the level 1 file, then search
statDir = fileparts(fileparts(fileparts(fn)));     % ...\StationaryDep
survey_file = fullfile(statDir,'DRIFTERS','DATA','Stationary_Deployment_093026', ...
    'stationary_drifter_survey.csv');
if ~isfile(survey_file)
    survey_file = fullfile(fileparts(fn),'stationary_drifter_survey.csv');
end
if ~isfile(survey_file)
    hit = dir(fullfile(statDir,'**','stationary_drifter_survey.csv'));
    if ~isempty(hit), survey_file = fullfile(hit(1).folder, hit(1).name); end
end

if ~isfile(survey_file)
    fprintf('\nNo survey file at %s - skipping the true-position comparison.\n', survey_file);
else
    Stab = readtable(survey_file, 'ReadVariableNames', false, 'TextType', 'string');
    nm   = string(Stab.Var1);
    tok  = regexp(nm, '^[dD](\d+)(?:-\d+)?$', 'tokens', 'once');
    isD  = ~cellfun(@isempty, tok);
    sidD = cellfun(@(c) str2double(c{1}), tok(isD));

    % Same projection as the drifter pipeline (NC State Plane, EPSG:32119)
    proj = projcrs(32119);
    [xs, ys] = projfwd(proj, Stab.Var2(isD), Stab.Var3(isD));

    ids = unique(sidD)';
    rtkE = nan(numel(ids),1);  rtkN = rtkE;  rtkShots = rtkE;  rtkRep = rtkE;
    for a = 1:numel(ids)
        m = sidD == ids(a);
        rtkE(a)     = mean(xs(m));
        rtkN(a)     = mean(ys(m));
        rtkShots(a) = sum(m);
        rtkRep(a)   = 100*sqrt(std(xs(m))^2 + std(ys(m))^2);     % cm
    end

    errE = cell(nT,1);  errN = cell(nT,1);
    Sv = struct('ID',[T.ID]', 'rtkShots',nan(nT,1), 'rtkRep_cm',nan(nT,1), ...
        'biasE',nan(nT,1), 'biasN',nan(nT,1), 'bias',nan(nT,1), ...
        'rmsErr',nan(nT,1), 'cep50',nan(nT,1), 'r95',nan(nT,1));

    for q = 1:nT
        a = find(ids == T(q).ID, 1);
        if isempty(a)
            fprintf('Track ID %d: no survey point with that drifter number\n', T(q).ID);
            continue
        end
        errE{q} = T(q).E - rtkE(a);
        errN{q} = T(q).N - rtkN(a);
        dd = hypot(errE{q}, errN{q});

        Sv.rtkShots(q)  = rtkShots(a);
        Sv.rtkRep_cm(q) = rtkRep(a);
        Sv.biasE(q)     = mean(errE{q},'omitnan');
        Sv.biasN(q)     = mean(errN{q},'omitnan');
        Sv.bias(q)      = hypot(Sv.biasE(q), Sv.biasN(q));
        Sv.rmsErr(q)    = sqrt(mean(dd.^2,'omitnan'));
        Sv.cep50(q)     = pct(dd,50);
        Sv.r95(q)       = pct(dd,95);
    end

    fprintf('\n=== GPS position vs RTK truth (EPSG:32119), m ===\n');
    fprintf('%-6s %6s %10s %8s %8s %8s %8s %8s %11s\n', 'ID','shots','RTK rep cm', ...
        'bias E','bias N','|bias|','RMS err','CEP50','95% radius');
    for q = 1:nT
        if isnan(Sv.bias(q)), continue, end
        fprintf('%-6d %6d %10.1f %+8.2f %+8.2f %8.2f %8.2f %8.2f %11.2f\n', Sv.ID(q), ...
            Sv.rtkShots(q), Sv.rtkRep_cm(q), Sv.biasE(q), Sv.biasN(q), Sv.bias(q), ...
            Sv.rmsErr(q), Sv.cep50(q), Sv.r95(q));
    end

    ok = ~isnan(Sv.bias);
    cE = mean(Sv.biasE(ok));
    cN = mean(Sv.biasN(ok));
    fprintf('\nCommon offset (mean over drifters): E %+.2f m, N %+.2f m\n', cE, cN);
    fprintf('Spread of bias between drifters (std): E %.2f m, N %.2f m\n', ...
        std(Sv.biasE(ok)), std(Sv.biasN(ok)));
    fprintf('Precision (scatter about own mean, mean RMS 2D): %.2f m\n', mean(P.rms2D(ok)));
    fprintf('Accuracy (RMS error vs RTK truth, mean over drifters): %.2f m\n', mean(Sv.rmsErr(ok)));
    fprintf('\nSummary (mean [min-max] over drifters), m:\n');
    fprintf('  RMS of positions about own mean (precision):        %.2f [%.2f-%.2f]\n', ...
        mean(P.rms2D(ok)), min(P.rms2D(ok)), max(P.rms2D(ok)));
    fprintf('  Error of the mean position vs RTK truth (|bias|):   %.2f [%.2f-%.2f]\n', ...
        mean(Sv.bias(ok)), min(Sv.bias(ok)), max(Sv.bias(ok)));
    fprintf('  RMS of positions about RTK truth (total):           %.2f [%.2f-%.2f]\n', ...
        mean(Sv.rmsErr(ok)), min(Sv.rmsErr(ok)), max(Sv.rmsErr(ok)));
    fprintf('Accuracy after removing the common offset: %.2f m\n', ...
        sqrt(mean((Sv.biasE(ok)-cE).^2 + (Sv.biasN(ok)-cN).^2) + mean(P.rms2D(ok).^2)));

    figure('Name','StationaryDep_RTK_error'); hold on
    for q = 1:nT
        if isempty(errE{q}), continue, end
        plot(errE{q}, errN{q}, '.', 'Color', cols(q,:), 'MarkerSize', 4);
    end
    for q = 1:nT
        if isnan(Sv.bias(q)), continue, end
        plot(Sv.biasE(q), Sv.biasN(q), 'ko', 'MarkerFaceColor', cols(q,:), 'MarkerSize', 9);
        text(Sv.biasE(q), Sv.biasN(q), sprintf('  %d', Sv.ID(q)), 'FontWeight','bold');
    end
    plot(0, 0, 'k+', 'MarkerSize', 16, 'LineWidth', 2);
    axis equal; grid on
    xlabel('East error vs RTK (m)'); ylabel('North error vs RTK (m)')
%     title('GPS error vs RTK truth (circles = mean error per drifter, + = truth)')
end

%% Figures
cols = lines(nT);

figure('Name','StationaryDep_position_scatter')
subplot(1,2,1); hold on
for q = 1:nT
    plot(T(q).E - mean(T(q).E,'omitnan'), T(q).N - mean(T(q).N,'omitnan'), '.', ...
        'Color', cols(q,:), 'MarkerSize', 4);
end
axis equal; grid on
xlabel('East error (m)'); ylabel('North error (m)')
title('Position error about own mean')

subplot(1,2,2); hold on
for q = 1:nT
    rr = hypot(T(q).E - mean(T(q).E,'omitnan'), T(q).N - mean(T(q).N,'omitnan'));
    plot((0:numel(rr)-1)'/3600, rr, 'Color', cols(q,:));
end
grid on
xlabel('Hours since start'); ylabel('Distance from mean position (m)')
title('Look for handling or real motion')

figure('Name','StationaryDep_velocity_floor'); hold on
for q = 1:nT
    loglog(windows, 100*F.meanSpd(q,:), '-', 'Color', [0.7 0.7 0.7]);
end
h1 = loglog(windows, 100*mean(F.meanSpd,1,'omitnan'), 'k-o', 'LineWidth', 2);
h2 = loglog(windows, 100*mean(F.rmsSpd,1,'omitnan'),  'b-o', 'LineWidth', 2);
h3 = loglog(windows, 100*p95pool,                     'r-o', 'LineWidth', 2);
h4 = loglog(windows, 100*mean(F.rmsSpd(:,1),'omitnan')*sqrt(windows(1)./windows), 'k--');       % 1/sqrt(window), anchored to the RMS curve at 2 s
set(gca,'XScale','log','YScale','log'); grid on
xlabel('Window length (s)'); ylabel('Apparent speed of non-moving drifters (cm/s)')
legend([h1 h2 h3 h4], 'Mean speed','RMS speed','95th percentile (pooled)', ...
    '1/\surd(window) scaling', 'Location','southwest')
% title('Velocity noise floor vs window length')

figure('Name','StationaryDep_error_autocorrelation'); hold on
lagMin = (0:maxLag)'/60;
for q = 1:nT
    plot(lagMin, 0.5*(rE(:,q)+rN(:,q)), 'Color', [0.7 0.7 0.7]);
end
plot(lagMin, mean(0.5*(rE+rN),2), 'k', 'LineWidth', 2);
yline(exp(-1),'--'); grid on
xlabel('Lag (min)'); ylabel('Autocorrelation of position error')
title('How long GPS position errors stay correlated')

figure('Name','StationaryDep_velocity_autocorrelation')
lagS = (0:maxLag)';
for sp = 1:2
    subplot(1,2,sp); hold on
    hA = plot(lagS, mean(rV1,2),  'b', 'LineWidth', 2);
%     hB = plot(lagS, mean(rV8,2),  'r', 'LineWidth', 2);
    hC = plot(lagS, mean(rDop,2), 'k', 'LineWidth', 2);
    hD = plot(lagS, mean(0.5*(rE+rN),2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    yline(exp(-1),':'); yline(0,'-','Color',[0.8 0.8 0.8]);
    grid on
    xlabel('Lag (s)'); ylabel('Autocorrelation')
    if sp == 1
        xlim([0 60]); title('First minute')
%         legend([hA hB hC hD], 'Position velocity, 1 s', 'Position velocity, 8 s window', ...
%             'GPS Doppler speed', 'Position error', 'Location','northeast')
        legend([hA hC hD], 'Position velocity, 1 s', ...
            'GPS Doppler speed', 'Position error', 'Location','northeast')
    else
        xlim([0 maxLag]); title('Full lag range')
    end
end

%% Save figures
% Thesis FIGURES folder: use whichever OneDrive root exists on this machine
kf = find(cellfun(@(r) isfolder(fullfile(r,'THESIS','THESIS_DRAFTING','FIGURES')), roots), 1);
if isempty(kf)
    fig_dir = fullfile(pwd,'STATIONARY_DEP');
    warning('THESIS_DRAFTING\\FIGURES not found; saving figures to %s', fig_dir);
else
    fig_dir = fullfile(roots{kf},'THESIS','THESIS_DRAFTING','FIGURES','STATIONARY_DEP');
end
if ~isfolder(fig_dir), mkdir(fig_dir); end
figs = findobj('Type','figure');
for k = 1:numel(figs)
    nm = get(figs(k),'Name');
    if ~startsWith(nm,'StationaryDep_'), continue, end
    set(figs(k),'Color','w');
    exportgraphics(figs(k), fullfile(fig_dir,[nm '.png']), 'Resolution',300);
    savefig(figs(k), fullfile(fig_dir,[nm '.fig']));
    fprintf('Saved figure %s\n', fullfile(fig_dir,[nm '.png']));
end

%% Save
Floor.windows_s   = windows;
Floor.ID          = [T.ID];
Floor.Position    = P;
Floor.Velocity    = F;
Floor.p95Pooled   = p95pool;
Floor.efold_s     = efold;
Floor.maxLag_s    = maxLag;
Floor.DopplerMean = mean(dop);
Floor.DopplerP95  = pct(dop,95);
if ~isempty(Sv)
    Floor.Survey = Sv;
    Floor.RTK    = struct('ID',ids,'E',rtkE,'N',rtkN,'shots',rtkShots,'rep_cm',rtkRep);
end
Floor.VelDecorr = struct('lag_s',lagS, 'r_pos1s',rV1, 'r_pos8s',rV8, 'r_doppler',rDop, ...
    'efold_firstzero_pos1s_s',tV1, 'efold_firstzero_pos8s_s',tV8, ...
    'efold_firstzero_doppler_s',tDop, 'efold_firstzero_poserror_s',tPos);
Floor.created     = datestr(now);
save(out_file,'Floor');
fprintf('\nSaved %s\n', out_file);

out = Floor;

end


%% ========================= Local functions =========================

function p = pct(x, q)
% q-th percentile without the Statistics Toolbox
x = sort(x(~isnan(x)));
if isempty(x), p = NaN; return, end
p = x(max(1, min(numel(x), round(q/100*numel(x)))));
end


function t = lagTimes(r)
% [first lag where r < 1/e, first lag where r <= 0] in samples (= seconds); NaN if never
t = nan(1,2);
k = find(r < exp(-1), 1);   if ~isempty(k), t(1) = k - 1; end
k = find(r <= 0, 1);        if ~isempty(k), t(2) = k - 1; end
end


function s = fmtT(x, maxLag)
% Mean over drifters as text; "> maxLag" if none reached it, trailing + if some did not
if all(isnan(x))
    s = sprintf('> %d', maxLag);
elseif any(isnan(x))
    s = sprintf('>%.1f', mean(x,'omitnan'));
else
    s = sprintf('%.1f', mean(x));
end
end


function r = acf(x, maxLag)
% Autocorrelation (lags 0..maxLag) via FFT, normalised so lag 0 = 1
x = x(:) - mean(x,'omitnan');
x(isnan(x)) = 0;
n = numel(x);
nfft = 2^nextpow2(2*n - 1);
X = fft(x, nfft);
c = real(ifft(abs(X).^2));
r = c(1:maxLag+1) / c(1);
end