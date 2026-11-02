function out = stationaryVsMovingVelocity(cfg)
% STEP 2 (comparison with the stationary deployment):
% how far do the position-derived drifter velocities sit from GPS Doppler speed,
% compared with the same calculation on drifters that are not moving?
%
% Figure: two panels against averaging window.
%   left : RMSE between position velocity and Doppler speed AT THE SAME window
%          (Doppler speed is averaged over the same window)
%   right: mean difference (position minus Doppler), i.e. the bias
% Blue = moving drifters (mean over tracks, 95% interval by resampling deployments).
% Grey = stationary drifters (mean and range), i.e. the noise level.
% The dashed line marks the 8 s window used in the maps.
%
% STUCK DRIFTERS: with excludeStuck = true the moving-drifter statistics use only samples
% not flagged Stuck (Results.Stuck from flagStuckDrifters + veloAnalysis_v2), widened by
% dilate_s/2 each side. The stationary reference is unchanged. Outputs get the suffix
% _moving, and a table of all-samples vs moving-only is printed either way.
%
% Reading the results:
%   - Position velocity is the SIZE of the mean displacement vector; Doppler speed
%     is averaged as a scalar. If a drifter changes direction inside a window the
%     two differ even with perfect GPS, so the difference at long windows is
%     partly real.
%   - At short windows position noise inflates the position velocity (positive bias).

%% Settings
markWin   = 8;                            % window marked in the figure
nBoot     = 1000;

excludeStuck = true;                      % true = drop stuck samples; false = all samples
dilate_s     = 256;                       % mask widened by dilate_s/2 each side (longest window)
minMoving_s  = 300;                       % tracks with less moving data than this are dropped

obsFiles = {'MarshMadness\PROCESSED_DATA_V2\DRIFTERS\MarshMadness_TemporalRes_RMSE.mat'
    'MarshMayhem\PROCESSED_DATA_V2\DRIFTERS\MarshMayhem_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA_V2\DRIFTERS\FallFrolic_100825_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA_V2\DRIFTERS\FallFrolic_101025_TemporalRes_RMSE.mat'
    'NovDep\PROCESSED_DATA_V2\DRIFTERS\NovDep_TemporalRes_RMSE.mat'
    'JuneJamboree\PROCESSED_DATA_V2\DRIFTERS\JuneJamboree_TemporalRes_RMSE.mat'
    'SummerFlood\PROCESSED_DATA_V2\DRIFTERS\SummerFlood_TemporalRes_RMSE.mat'};
names = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};

%% Folders (works on either computer)
roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
kr = find(cellfun(@(r) isfolder(fullfile(r,'THESIS','THESIS_DRAFTING','FIGURES')), roots), 1);
if isempty(kr)
    fig_dir = pwd;
    warning('FIGURES folder not found; saving to %s', fig_dir);
else
    fig_dir = fullfile(roots{kr},'THESIS','THESIS_DRAFTING','FIGURES','STATIONARY_DEP');
    if ~isfolder(fig_dir), mkdir(fig_dir); end
end

%% Load every observation period that exists
data = struct();
for i = 1:numel(obsFiles)
    fn = fullfile(cfg.root_data, obsFiles{i});
    if ~isfile(fn)
        fprintf('%s: file not found, skipped\n', names{i});
        continue
    end
    data.(names{i}) = load(fn).Results;
end
obs_names = fieldnames(data);
rng(1);

%% Per-track statistics (chosen mode, plus the other mode for the comparison)
if ~isfield(data.(obs_names{1}), 'Stuck')
    error(['Results has no Stuck field. Run main_driver with flagStuckDrifters and ' ...
           'veloAnalysis_v2 so the TemporalRes_RMSE files contain it.']);
end
G  = gatherMoving(data, obs_names, excludeStuck,  dilate_s, minMoving_s);
G0 = gatherMoving(data, obs_names, ~excludeStuck, dilate_s, minMoving_s);
if excludeStuck, Gmov = G;  Gall = G0; else, Gmov = G0; Gall = G; end

windows = G.windows;  nW = G.nW;  nTrk = G.nTrk;
rmseDop = G.rmseDop;  biasDop = G.biasDop;  rmseMax = G.rmseMax;
trkCl = G.trkCl;  trkPer = G.trkPer;  trkID = G.trkID;
if nTrk == 0, error('No tracks found in the loaded files.'); end
fprintf('\nTracks used: %d from %d observation periods (%d deployment clusters)\n', ...
    nTrk, numel(unique(trkPer)), numel(unique(trkCl)));
if excludeStuck
    fprintf('Stuck mask ON: %.1f%% of samples removed (mask widened by %g s); %d tracks dropped (< %g s moving data)\n', ...
        100*G.nMasked/max(G.nSamples,1), dilate_s/2, G.nDropped, minMoving_s);
    ttag = 'moving samples only';  sfx = '_moving';
else
    fprintf('Stuck mask OFF: all samples used (%.1f%% of them are flagged stuck)\n', ...
        100*G0.nMasked/max(G0.nSamples,1));
    ttag = 'all samples';  sfx = '';
end

%% Stationary reference (same calculation on the non-moving drifters)
SR = [];
sroot = '';
for q = 1:numel(roots)
    if isfolder(fullfile(roots{q},'CB_YachtBasin')), sroot = roots{q}; break, end
end
if ~isempty(sroot)
    sf = fullfile(sroot,'CB_YachtBasin','StationaryDep','PROCESSED_DATA_V2','DRIFTERS', ...
        'StationaryDep_drifters_level_1.mat');
    if isfile(sf), SR = stationaryReference(sf, windows); end
end
if isempty(SR), fprintf('No stationary reference found; figure 1 will not show it.\n'); end

%% Figure 1: agreement vs window
[mR, loR, hiR] = clusterBoot(rmseDop, trkCl, nBoot);
[mB, loB, hiB] = clusterBoot(biasDop, trkCl, nBoot);
[mM, ~, ~]     = clusterBoot(rmseMax, trkCl, nBoot);

f1 = figure('Name','StationaryComparison_window','Units','pixels','Position',[80 80 1250 520],'Color','w');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% title(tl, sprintf('Position velocity vs GPS Doppler speed, by averaging window (%d tracks)', nTrk), ...
%     'FontWeight','bold');
title(tl, sprintf('Position velocity vs 1 Hz GPS Doppler speed (%d tracks, %s)', nTrk, ttag), ...
    'FontWeight','bold');


for pnl = 1:2
    ax = nexttile; hold on
    if pnl == 1
        mm = mR; lo = loR; hi = hiR;  yl = 'RMSE (m/s)';
        if ~isempty(SR), sm = mean(SR.rmse,1,'omitnan'); smin = min(SR.rmse,[],1); smax = max(SR.rmse,[],1); end
    else
        mm = mB; lo = loB; hi = hiB;  yl = 'Position minus Doppler (m/s)';
        if ~isempty(SR), sm = mean(SR.bias,1,'omitnan'); smin = min(SR.bias,[],1); smax = max(SR.bias,[],1); end
    end
    hS = gobjects(0);
    if ~isempty(SR)
        fill([windows fliplr(windows)], [smin fliplr(smax)], [0.7 0.7 0.7], ...
            'FaceAlpha',0.35, 'EdgeColor','none');
        hS = plot(windows, sm, '-', 'Color',[0.45 0.45 0.45], 'LineWidth',2);
    end
    fill([windows fliplr(windows)], [lo fliplr(hi)], [0.2 0.45 0.8], 'FaceAlpha',0.25, 'EdgeColor','none');
    hM = plot(windows, mm, '-o', 'Color',[0.1 0.35 0.7], 'LineWidth',2, 'MarkerFaceColor',[0.1 0.35 0.7]);
    if pnl == 2, yline(0,'-','Color',[0.5 0.5 0.5]); end
    xline(markWin,'--','Color',[0.3 0.3 0.3]);
    set(ax,'XScale','log','XTick',windows,'XTickLabel',arrayfun(@num2str,windows,'UniformOutput',false));
    xlim([windows(1) windows(end)]);
    grid on; box on
    xlabel('Averaging window (s)'); ylabel(yl);
    if pnl == 1
        ymaxv = max([hi(:); 0.001]);
        if ~isempty(SR), ymaxv = max(ymaxv, max(smax)); end
        ylim([0 ymaxv*1.1]);
        if ~isempty(SR)
            legend([hM hS], 'Moving drifters (mean, 95% interval)', 'Stationary drifters (mean, range)', 'Location','northeast');
        else
            legend(hM, 'Moving drifters (mean, 95% interval)', 'Location','northeast');
        end
        title(sprintf('RMSE (dashed line = %d s window used in maps)', markWin));
    else
        title('Bias');
    end
end
savefig(f1, fullfile(fig_dir,['StationaryComparison_window' sfx '.fig']));
exportgraphics(f1, fullfile(fig_dir,['StationaryComparison_window' sfx '.png']), 'Resolution',300);

%% Printed table
fprintf('\n=== Position velocity vs Doppler speed, by window (%d tracks, %s) ===\n', nTrk, ttag);
fprintf('%-8s %-24s %12s %14s %18s\n', 'window','RMSE vs Doppler [95%]','bias','stationary RMSE','diff from 256 s');
for w = 1:nW
    if isempty(SR), st = '-'; else, st = sprintf('%.4f', mean(SR.rmse(:,w),'omitnan')); end
    fprintf('%-8s %.4f [%.4f-%.4f] %+12.4f %14s %18.4f\n', sprintf('%d s',windows(w)), ...
        mR(w), loR(w), hiR(w), mB(w), st, mM(w));
end
fprintf('(m/s; "diff from 256 s" = RMS difference to the 256 s velocity, how much the velocity changes with window, not an error)\n');

%% Comparison: all samples vs moving samples only
[mRall, ~, ~] = clusterBoot(Gall.rmseDop, Gall.trkCl, nBoot);
[mRmov, ~, ~] = clusterBoot(Gmov.rmseDop, Gmov.trkCl, nBoot);
[mBall, ~, ~] = clusterBoot(Gall.biasDop, Gall.trkCl, nBoot);
[mBmov, ~, ~] = clusterBoot(Gmov.biasDop, Gmov.trkCl, nBoot);
fprintf('\n=== Effect of the stuck mask on moving drifters vs Doppler (m/s) ===\n');
fprintf('all samples: %d tracks;  moving only: %d tracks (%d dropped, < %g s moving data)\n', ...
    Gall.nTrk, Gmov.nTrk, Gmov.nDropped, minMoving_s);
fprintf('%-8s %9s %9s %8s %9s %9s %12s %10s\n', 'window','RMSE all','RMSE mov','change','bias all','bias mov','stationary','mov/stat');
for w = 1:nW
    if isempty(SR), stv = NaN; else, stv = mean(SR.rmse(:,w),'omitnan'); end
    fprintf('%-8s %9.4f %9.4f %7.0f%% %+9.4f %+9.4f %12.4f %9.0f%%\n', sprintf('%d s',windows(w)), ...
        mRall(w), mRmov(w), 100*(mRmov(w)/mRall(w)-1), mBall(w), mBmov(w), stv, 100*stv/mRmov(w));
end
fprintf('(mov/stat = stationary RMSE as %% of the moving-only RMSE)\n');

%% Save
out = struct('windows_s',windows, 'nTracks',nTrk, 'trackCluster',trkCl, 'trackPeriod',trkPer, 'trackID',trkID, ...
    'rmseVsDoppler',rmseDop, 'biasVsDoppler',biasDop, 'rmseVs256s',rmseMax, ...
    'meanRMSE',mR, 'RMSE_CI',[loR; hiR], 'meanBias',mB, 'bias_CI',[loB; hiB], ...
    'stationary',SR, 'excludeStuck',excludeStuck, 'dilate_s',dilate_s, 'minMoving_s',minMoving_s, ...
    'pctSamplesMasked',100*G.nMasked/max(G.nSamples,1), 'nTracksDropped',G.nDropped, ...
    'compare_RMSE_all',mRall, 'compare_RMSE_moving',mRmov, 'compare_bias_all',mBall, 'compare_bias_moving',mBmov, ...
    'created',datestr(now));
save(fullfile(fig_dir,['StationaryComparison_results' sfx '.mat']), 'out');
fprintf('\nSaved figures and results to %s\n', fig_dir);
end


%% ========================= Local functions =========================

function G = gatherMoving(data, obs_names, excludeStuck, dilate_s, minMoving_s)
% Per-track RMSE and bias of position velocity vs 1 Hz Doppler speed, and RMSE vs the 256 s
% velocity. With excludeStuck, samples flagged Stuck (widened by dilate_s/2 each side so no
% averaging window overlaps a stuck period) are set to NaN before the statistics.
windows = [];  nW = 0;
rmseDop = [];  biasDop = [];  rmseMax = [];  trkCl = [];  trkPer = [];  trkID = [];
nTrk = 0;  nSamples = 0;  nMasked = 0;  nDropped = 0;
half = round(dilate_s/2);

for i = 1:numel(obs_names)
    R = data.(obs_names{i});
    for k = 1:numel(R)
        if isempty(R(k).ID) || isempty(R(k).Speed) || isempty(R(k).TemporalRes), continue, end
        [d, ~] = ind2sub(size(R), k);
        TR = R(k).TemporalRes;

        if isempty(windows), windows = [1, TR.Seconds]; nW = numel(windows); end
        Vmax = TR(end).V(:);

        % stuck mask for this track (true = drop), widened to cover the averaging window
        nn = min([numel(R(k).V_1Hz) numel(Vmax) numel(R(k).Speed)]);
        Sd = widenStuckMask(R(k).Stuck, nn, half);
        nSamples = nSamples + nn;
        nMasked  = nMasked + sum(Sd);
        if excludeStuck && sum(~Sd) < minMoving_s
            nDropped = nDropped + 1;
            continue
        end

        nTrk = nTrk + 1;
        rmseDop(nTrk,1:nW) = nan;  biasDop(nTrk,1:nW) = nan;  rmseMax(nTrk,1:nW) = nan;   %#ok<AGROW>
        trkCl(nTrk,1)  = i*100 + d;                                                         %#ok<AGROW>
        trkPer(nTrk,1) = i;                                                                 %#ok<AGROW>
        trkID(nTrk,1)  = R(k).ID;                                                           %#ok<AGROW>

        for w = 1:nW
            if w == 1
                V = R(k).V_1Hz(:);
            else
                V = TR(w-1).V(:);
            end
            S = R(k).Speed(:);          % 1 Hz Doppler speed, as in Section 2
            n = min([numel(V) numel(S) numel(Vmax) numel(R(k).Speed)]);
            V = V(1:n);  S = S(1:n);  Vm = Vmax(1:n);
            if excludeStuck
                V(Sd(1:n)) = NaN;
            end

            dd = V - S;
            rmseDop(nTrk,w) = sqrt(mean(dd.^2,'omitnan'));
            biasDop(nTrk,w) = mean(dd,'omitnan');
            rmseMax(nTrk,w) = sqrt(mean((V - Vm).^2,'omitnan'));
        end
    end
end
G = struct('windows',windows, 'nW',nW, 'nTrk',nTrk, 'rmseDop',rmseDop, 'biasDop',biasDop, ...
    'rmseMax',rmseMax, 'trkCl',trkCl, 'trkPer',trkPer, 'trkID',trkID, ...
    'nSamples',nSamples, 'nMasked',nMasked, 'nDropped',nDropped);
end


function [m, lo, hi] = clusterBoot(vals, cl, nBoot)
% Mean over tracks with a 95% interval from resampling whole clusters (deployments)
m = mean(vals,1,'omitnan');
u = unique(cl);  nU = numel(u);
lo = m;  hi = m;
if nU < 2, return, end
grp = arrayfun(@(x) find(cl == x), u, 'UniformOutput', false);
B = nan(nBoot, size(vals,2));
for b = 1:nBoot
    ix = vertcat(grp{randi(nU, nU, 1)});
    B(b,:) = mean(vals(ix,:),1,'omitnan');
end
B  = sort(B,1);
lo = B(max(1,round(0.025*nBoot)),:);
hi = B(min(nBoot,round(0.975*nBoot)),:);
end


function SR = stationaryReference(sfile, windows)
% Same calculation as for the moving drifters, on drifters that are not moving:
% RMSE and mean difference between position velocity and Doppler speed averaged
% over the same window. Uses the centred-window definition of veloAnalysis.
L1 = load(sfile).drifters_level_1;
res = windows(windows > 1);
maxHalf = max(res)/2;
rm = [];  bs = [];
for d = 1:size(L1,1)
    for j = 1:size(L1,2)
        D = L1(d,j);
        if isempty(D.East) || isempty(D.North) || isempty(D.e_time) || isempty(D.Speed) || isempty(D.V), continue, end
        t = D.e_time(:);  E = D.East(:);  N = D.North(:);
        Spd = D.Speed(:)*0.51444;  V1 = D.V(:);
        if numel(t) <= 2*maxHalf + 10, continue, end
        c = (1+maxHalf):(numel(t)-maxHalf);
        rmRow = nan(1,numel(windows));  bsRow = rmRow;
        for w = 1:numel(windows)
            L = windows(w);
            h  = max(L/2, 0);
            S  = Spd(c);                % 1 Hz Doppler speed
            if L == 1
                V = V1(c);
            else
                dt = t(c+h) - t(c-h);
                V  = hypot(E(c+h) - E(c-h), N(c+h) - N(c-h)) ./ dt;
            end
%             if L == 1
%                 V = V1(c);  S = Spd(c);
%             else
%                 h  = L/2;
%                 dt = t(c+h) - t(c-h);
%                 V  = hypot(E(c+h) - E(c-h), N(c+h) - N(c-h)) ./ dt;
%                 S  = movmean(Spd, L, 'omitnan');
%                 S  = S(c);
%             end
            dd = V - S;
            rmRow(w) = sqrt(mean(dd.^2,'omitnan'));
            bsRow(w) = mean(dd,'omitnan');
        end
        rm(end+1,:) = rmRow;  %#ok<AGROW>
        bs(end+1,:) = bsRow;  %#ok<AGROW>
    end
end
SR = struct('rmse',rm, 'bias',bs);
end
