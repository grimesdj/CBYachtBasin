function [out, all_RMSE] = allVeloAnalysis_v2(cfg)
% Cleaner version of allVeloAnalysis: two figures instead of four, one table.
% allVeloAnalysis.m is left untouched.
%
% Question answered: how well does the position-derived drifter velocity agree
% with GPS Doppler speed, and how does that depend on (1) the averaging window
% and (2) how fast the drifter is moving?
%
% Figure 1  Agreement vs window length (two panels)
%             left : RMSE between position velocity and Doppler speed AT THE SAME
%                    window (Doppler speed is averaged over the same window)
%             right: mean difference (position minus Doppler), i.e. the bias
%           Blue = moving drifters (mean over tracks, 95% interval by resampling
%           deployments). Grey = the same calculation for the STATIONARY drifters,
%           i.e. the noise level. The dashed line marks the 8 s window used in
%           the maps.
% Figure 2  Agreement vs drifter speed (two panels), for four windows
%             left : RMSE, right: RMSE as a percentage of the mean speed in the bin
%
% Printed table (for the thesis log): the numbers behind figure 1, plus the
% difference from the 256 s velocity (a measure of how much the velocity changes
% with window, NOT an error).
%
% Reading the results:
%   - Position velocity is the SIZE of the mean displacement vector; Doppler speed
%     is averaged as a scalar. If a drifter changes direction inside a window the
%     two differ even with perfect GPS, so the difference at long windows is
%     partly real.
%   - At short windows position noise inflates the position velocity (positive bias).
%   - N in the figures is tracks and 1 Hz samples; samples overlap, so tracks are
%     the number to quote.

%% Settings
showWin   = [2 8 32 128];                 % windows drawn in figure 2
markWin   = 8;                            % window marked in figure 1
vel_edges = [0 0.05 0.1 0.3 0.6];         % m/s, speed bins (same as before)
binBy     = 'doppler';                    % 'doppler' = 1 Hz Doppler speed (as before)
                                          % 'vmax'    = 256 s position velocity
nBoot     = 1000;

obsFiles = {'MarshMadness\PROCESSED_DATA\DRIFTERS\MarshMadness_TemporalRes_RMSE.mat'
    'MarshMayhem\PROCESSED_DATA\DRIFTERS\MarshMayhem_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_100825_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_101025_TemporalRes_RMSE.mat'
    'NovDep\PROCESSED_DATA\DRIFTERS\NovDep_TemporalRes_RMSE.mat'
    'JuneJamboree\PROCESSED_DATA\DRIFTERS\JuneJamboree_TemporalRes_RMSE.mat'
    'SummerFlood\PROCESSED_DATA\DRIFTERS\SummerFlood_TemporalRes_RMSE.mat'};
names = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};

%% Folders (works on either computer)
roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
kr = find(cellfun(@(r) isfolder(fullfile(r,'THESIS','THESIS_DRAFTING','FIGURES')), roots), 1);
if isempty(kr)
    fig_dir = pwd;
    warning('FIGURES folder not found; saving to %s', fig_dir);
else
    fig_dir = fullfile(roots{kr},'THESIS','THESIS_DRAFTING','FIGURES','RMSE');
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
all_RMSE = data;
rng(1);

%% Per-track statistics and pooled samples
windows = [];                         % [1 2 4 ... 256]
rmseDop = [];  biasDop = [];  rmseMax = [];  trkCl = [];  trkPer = [];  trkID = [];
poolV = {};  poolS = {};  poolB = {};  poolT = {};
nTrk = 0;

for i = 1:numel(obs_names)
    R = data.(obs_names{i});
    for k = 1:numel(R)
        if isempty(R(k).ID) || isempty(R(k).Speed) || isempty(R(k).TemporalRes), continue, end
        [d, ~] = ind2sub(size(R), k);
        TR = R(k).TemporalRes;

        if isempty(windows), windows = [1, TR.Seconds]; nW = numel(windows); poolV = repmat({{}},nW,1); poolS = poolV; poolB = poolV; poolT = poolV; end
        Vmax = TR(end).V(:);

        nTrk = nTrk + 1;
        rmseDop(nTrk,1:nW) = nan;  biasDop(nTrk,1:nW) = nan;  rmseMax(nTrk,1:nW) = nan;   %#ok<AGROW>
        trkCl(nTrk,1)  = i*100 + d;                                                         %#ok<AGROW>
        trkPer(nTrk,1) = i;                                                                 %#ok<AGROW>
        trkID(nTrk,1)  = R(k).ID;                                                           %#ok<AGROW>

        for w = 1:nW
            if w == 1
                V = R(k).V_1Hz(:);   S = R(k).Speed(:);
            else
                V = TR(w-1).V(:);    S = TR(w-1).Speed(:);
            end
            n = min([numel(V) numel(S) numel(Vmax) numel(R(k).Speed)]);
            V = V(1:n);  S = S(1:n);  Vm = Vmax(1:n);

            dd = V - S;
            rmseDop(nTrk,w) = sqrt(mean(dd.^2,'omitnan'));
            biasDop(nTrk,w) = mean(dd,'omitnan');
            rmseMax(nTrk,w) = sqrt(mean((V - Vm).^2,'omitnan'));

            if strcmp(binBy,'vmax')
                B = Vm;
            else
                B = R(k).Speed(1:n); B = B(:);
            end
            poolV{w}{end+1,1} = V;                      %#ok<AGROW>
            poolS{w}{end+1,1} = S;                      %#ok<AGROW>
            poolB{w}{end+1,1} = B;                      %#ok<AGROW>
            poolT{w}{end+1,1} = nTrk*ones(n,1);         %#ok<AGROW>
        end
    end
end
if nTrk == 0, error('No tracks found in the loaded files.'); end
fprintf('\nTracks used: %d from %d observation periods (%d deployment clusters)\n', ...
    nTrk, numel(unique(trkPer)), numel(unique(trkCl)));

%% Stationary reference (same calculation on the non-moving drifters)
SR = [];
sroot = '';
for q = 1:numel(roots)
    if isfolder(fullfile(roots{q},'CB_YachtBasin')), sroot = roots{q}; break, end
end
if ~isempty(sroot)
    sf = fullfile(sroot,'CB_YachtBasin','StationaryDep','PROCESSED_DATA','DRIFTERS', ...
        'StationaryDep_drifters_level_1.mat');
    if isfile(sf), SR = stationaryReference(sf, windows); end
end
if isempty(SR), fprintf('No stationary reference found; figure 1 will not show it.\n'); end

%% Figure 1: agreement vs window
[mR, loR, hiR] = clusterBoot(rmseDop, trkCl, nBoot);
[mB, loB, hiB] = clusterBoot(biasDop, trkCl, nBoot);
[mM, ~, ~]     = clusterBoot(rmseMax, trkCl, nBoot);

f1 = figure('Name','VeloQAQC_window','Units','pixels','Position',[80 80 1250 520],'Color','w');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Position velocity vs GPS Doppler speed, by averaging window (%d tracks)', nTrk), ...
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
savefig(f1, fullfile(fig_dir,'VeloQAQC_window.fig'));
exportgraphics(f1, fullfile(fig_dir,'VeloQAQC_window.png'), 'Resolution',300);

%% Speed bins (pooled 1 Hz samples, tracks counted separately)
nBins = numel(vel_edges) - 1;
vel_labels = arrayfun(@(a,b) sprintf('%.2f-%.2f',a,b), vel_edges(1:end-1), vel_edges(2:end), 'UniformOutput', false);
rmseBin = nan(nBins,nW);  relBin = nan(nBins,nW);  nSamp = zeros(nBins,nW);  nTrkBin = zeros(nBins,nW);

for w = 1:nW
    Vall = vertcat(poolV{w}{:});  Sall = vertcat(poolS{w}{:});
    Ball = vertcat(poolB{w}{:});  Tall = vertcat(poolT{w}{:});
    bin  = discretize(Ball, vel_edges);
    for b = 1:nBins
        ix = bin == b & ~isnan(Vall) & ~isnan(Sall);
        nSamp(b,w)   = sum(ix);
        nTrkBin(b,w) = numel(unique(Tall(ix)));
        if ~any(ix), continue, end
        rmseBin(b,w) = sqrt(mean((Vall(ix) - Sall(ix)).^2));
        relBin(b,w)  = 100*rmseBin(b,w) / mean(Sall(ix));
    end
end

%% Figure 2: agreement vs speed
wi = find(ismember(windows, showWin));
cols = parula(numel(wi)+1);

f2 = figure('Name','VeloQAQC_speedbins','Units','pixels','Position',[80 80 1250 520],'Color','w');
tl2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
if strcmp(binBy,'vmax'), bl = '256 s position velocity'; else, bl = '1 Hz Doppler speed'; end
title(tl2, sprintf('Position velocity vs Doppler speed, by drifter speed (bins of %s)', bl), 'FontWeight','bold');

xl = cell(1,nBins);
for b = 1:nBins
    xl{b} = sprintf('%s\\newline%d tracks', vel_labels{b}, nTrkBin(b,1));
end
for pnl = 1:2
    ax = nexttile; hold on
    for q = 1:numel(wi)
        if pnl == 1, yv = rmseBin(:,wi(q)); else, yv = relBin(:,wi(q)); end
        plot(1:nBins, yv, '-o', 'Color',cols(q,:), 'LineWidth',2, 'MarkerFaceColor',cols(q,:), ...
            'DisplayName', sprintf('%d s', windows(wi(q))));
    end
    set(ax,'XTick',1:nBins,'XTickLabel',xl);
    xlim([0.7 nBins+0.3]); grid on; box on
    xlabel('Drifter speed (m/s)');
    if pnl == 1
        ylabel('RMSE vs Doppler speed (m/s)'); ylim([0 inf]);
        legend('Location','northwest','Title','Window');
    else
        ylabel('RMSE as % of mean speed in bin'); ylim([0 inf]);
    end
end
savefig(f2, fullfile(fig_dir,'VeloQAQC_speedbins.fig'));
exportgraphics(f2, fullfile(fig_dir,'VeloQAQC_speedbins.png'), 'Resolution',300);

%% Printed table
fprintf('\n=== Position velocity vs Doppler speed, by window (%d tracks) ===\n', nTrk);
fprintf('%-8s %-24s %12s %14s %18s\n', 'window','RMSE vs Doppler [95%]','bias','stationary RMSE','diff from 256 s');
for w = 1:nW
    if isempty(SR), st = '-'; else, st = sprintf('%.4f', mean(SR.rmse(:,w),'omitnan')); end
    fprintf('%-8s %.4f [%.4f-%.4f] %+12.4f %14s %18.4f\n', sprintf('%d s',windows(w)), ...
        mR(w), loR(w), hiR(w), mB(w), st, mM(w));
end
fprintf('(m/s; "diff from 256 s" = RMS difference to the 256 s velocity, how much the velocity changes with window, not an error)\n');

fprintf('\n=== Tracks and samples by speed bin (%s) ===\n', bl);
for b = 1:nBins
    fprintf('%-10s %5d tracks %9d samples\n', vel_labels{b}, nTrkBin(b,1), nSamp(b,1));
end

%% Save
out = struct('windows_s',windows, 'nTracks',nTrk, 'trackCluster',trkCl, 'trackPeriod',trkPer, 'trackID',trkID, ...
    'rmseVsDoppler',rmseDop, 'biasVsDoppler',biasDop, 'rmseVs256s',rmseMax, ...
    'meanRMSE',mR, 'RMSE_CI',[loR; hiR], 'meanBias',mB, 'bias_CI',[loB; hiB], ...
    'speedBinEdges',vel_edges, 'binBy',binBy, 'rmseByBin',rmseBin, 'relRmseByBin',relBin, ...
    'samplesByBin',nSamp, 'tracksByBin',nTrkBin, 'stationary',SR, 'created',datestr(now));
save(fullfile(fig_dir,'VeloQAQC_results.mat'), 'out');
fprintf('\nSaved figures and results to %s\n', fig_dir);
end


%% ========================= Local functions =========================

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
            if L == 1
                V = V1(c);  S = Spd(c);
            else
                h  = L/2;
                dt = t(c+h) - t(c-h);
                V  = hypot(E(c+h) - E(c-h), N(c+h) - N(c-h)) ./ dt;
                S  = movmean(Spd, L, 'omitnan');
                S  = S(c);
            end
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
