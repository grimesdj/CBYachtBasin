function [out, all_RMSE] = allVeloAnalysis_v2(cfg)
% STEP 1 (QAQC of drifter data): how far is the position-derived velocity from
%   (a) the longest-window (256 s) position velocity, and
%   (b) the 1 Hz GPS Doppler speed,
% overall and within speed bins. allVeloAnalysis.m is left untouched.
%
% Figure 1  RMSE vs averaging window (one panel, two lines)
%             blue: vs the 256 s position velocity (the window itself is left out)
%             red : vs the 1 Hz Doppler speed
%           Mean over tracks, shaded = 95% interval from resampling deployments.
%           Dashed line = the 8 s window used in the maps.
% Figure 2  RMSE within speed bins (three panels, one line per window)
%             left  : vs the 256 s position velocity
%             middle: vs the 1 Hz Doppler speed
%             right : the middle panel as a percentage of the mean Doppler speed
%           Bin labels show how many tracks contribute.
%
% Reading the results:
%   - RMSE vs the 256 s velocity is NOT an error. It grows with window length
%     partly because the real flow changes within the longer window.
%   - RMSE vs the 1 Hz Doppler speed mixes position noise (large at short windows),
%     Doppler noise and real direction change (position velocity is the size of a
%     mean displacement; Doppler speed is a scalar speed).
%   - Bins are set by the 1 Hz Doppler speed, the same quantity used as reference,
%     so the lowest and highest bins can be pulled by noise in that speed. Set
%     binBy = 'vmax' to bin by the 256 s position velocity instead and compare.
%   - Samples at 1 Hz overlap heavily; quote tracks, not samples, as N.

%% Settings
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

%% Per-track RMSE and pooled samples
windows = [];                         % [1 2 4 ... 256]
rmseMax = [];  rmseDop = [];  trkCl = [];  trkPer = [];  trkID = [];
poolV = {};  poolM = {};  poolS = {};  poolB = {};  poolT = {};
nTrk = 0;  nW = 0;

for i = 1:numel(obs_names)
    R = data.(obs_names{i});
    for k = 1:numel(R)
        if isempty(R(k).ID) || isempty(R(k).Speed) || isempty(R(k).TemporalRes), continue, end
        [d, ~] = ind2sub(size(R), k);
        TR = R(k).TemporalRes;

        if isempty(windows)
            windows = [1, TR.Seconds];
            nW = numel(windows);
            poolV = repmat({{}},nW,1);  poolM = poolV;  poolS = poolV;  poolB = poolV;  poolT = poolV;
        end
        Vmax = TR(end).V(:);
        Dop  = R(k).Speed(:);                       % 1 Hz Doppler speed, m/s

        nTrk = nTrk + 1;
        rmseMax(nTrk,1:nW) = nan;                                  %#ok<AGROW>
        rmseDop(nTrk,1:nW) = nan;                                  %#ok<AGROW>
        trkCl(nTrk,1)  = i*100 + d;                                %#ok<AGROW>
        trkPer(nTrk,1) = i;                                        %#ok<AGROW>
        trkID(nTrk,1)  = R(k).ID;                                  %#ok<AGROW>

        for w = 1:nW
            if w == 1
                V = R(k).V_1Hz(:);
            else
                V = TR(w-1).V(:);
            end
            n  = min([numel(V) numel(Vmax) numel(Dop)]);
            V  = V(1:n);  Vm = Vmax(1:n);  Sp = Dop(1:n);

            rmseMax(nTrk,w) = sqrt(mean((V - Vm).^2,'omitnan'));
            rmseDop(nTrk,w) = sqrt(mean((V - Sp).^2,'omitnan'));

            if strcmp(binBy,'vmax'), B = Vm; else, B = Sp; end
            poolV{w}{end+1,1} = V;                      %#ok<AGROW>
            poolM{w}{end+1,1} = Vm;                     %#ok<AGROW>
            poolS{w}{end+1,1} = Sp;                     %#ok<AGROW>
            poolB{w}{end+1,1} = B;                      %#ok<AGROW>
            poolT{w}{end+1,1} = nTrk*ones(n,1);         %#ok<AGROW>
        end
    end
end
if nTrk == 0, error('No tracks found in the loaded files.'); end
fprintf('\nTracks used: %d from %d observation periods (%d deployment clusters)\n', ...
    nTrk, numel(unique(trkPer)), numel(unique(trkCl)));

%% Figure 1: RMSE vs window
[mM, loM, hiM] = clusterBoot(rmseMax, trkCl, nBoot);
[mD, loD, hiD] = clusterBoot(rmseDop, trkCl, nBoot);
wM = 1:nW-1;                                  % the 256 s window compared with itself is 0

cBlue = [0.1 0.35 0.7];  cRed = [0.8 0.2 0.15];
f1 = figure('Name','VeloQAQC_RMSE_window','Units','pixels','Position',[80 80 760 560],'Color','w');
ax = axes; hold on
fill([windows(wM) fliplr(windows(wM))], [loM(wM) fliplr(hiM(wM))], cBlue, 'FaceAlpha',0.2, 'EdgeColor','none');
fill([windows fliplr(windows)], [loD fliplr(hiD)], cRed, 'FaceAlpha',0.2, 'EdgeColor','none');
h1 = plot(windows(wM), mM(wM), '-o', 'Color',cBlue, 'LineWidth',2, 'MarkerFaceColor',cBlue);
h2 = plot(windows, mD, '-o', 'Color',cRed, 'LineWidth',2, 'MarkerFaceColor',cRed);
xline(markWin,'--','Color',[0.3 0.3 0.3]);
set(ax,'XScale','log','XTick',windows,'XTickLabel',arrayfun(@num2str,windows,'UniformOutput',false));
xlim([windows(1) windows(end)]);  ylim([0 inf]);
grid on; box on
xlabel('Averaging window (s)'); ylabel('RMSE of position-derived velocity (m/s)');
legend([h1 h2], sprintf('vs %d s position velocity', windows(end)), 'vs 1 Hz Doppler speed', 'Location','northwest');
title(sprintf('Position-derived velocity: RMSE by window (%d tracks; shading = 95%% interval; dashed = %d s)', nTrk, markWin), 'FontSize',11);
savefig(f1, fullfile(fig_dir,'VeloQAQC_RMSE_window.fig'));
exportgraphics(f1, fullfile(fig_dir,'VeloQAQC_RMSE_window.png'), 'Resolution',300);

%% Speed bins (pooled 1 Hz samples; tracks counted separately)
nBins = numel(vel_edges) - 1;
vel_labels = arrayfun(@(a,b) sprintf('%.2f-%.2f',a,b), vel_edges(1:end-1), vel_edges(2:end), 'UniformOutput', false);
rmseBinMax = nan(nBins,nW);  rmseBinDop = nan(nBins,nW);  relBinDop = nan(nBins,nW);
nSamp = zeros(nBins,nW);  nTrkBin = zeros(nBins,nW);

for w = 1:nW
    Vall = vertcat(poolV{w}{:});  Mall = vertcat(poolM{w}{:});
    Sall = vertcat(poolS{w}{:});  Ball = vertcat(poolB{w}{:});  Tall = vertcat(poolT{w}{:});
    bin  = discretize(Ball, vel_edges);
    for b = 1:nBins
        ix = bin == b & ~isnan(Vall);
        nSamp(b,w)   = sum(ix);
        nTrkBin(b,w) = numel(unique(Tall(ix)));
        if ~any(ix), continue, end
        rmseBinMax(b,w) = sqrt(mean((Vall(ix) - Mall(ix)).^2,'omitnan'));
        rmseBinDop(b,w) = sqrt(mean((Vall(ix) - Sall(ix)).^2,'omitnan'));
        relBinDop(b,w)  = 100*rmseBinDop(b,w) / mean(Sall(ix),'omitnan');
    end
end

%% Figure 2: RMSE within speed bins
if strcmp(binBy,'vmax'), bl = '256 s position velocity'; else, bl = '1 Hz Doppler speed'; end
cols = parula(nW+1);
xl = cell(1,nBins);
for b = 1:nBins
    xl{b} = sprintf('%s\\newline%d tracks', vel_labels{b}, nTrkBin(b,1));
end

f2 = figure('Name','VeloQAQC_RMSE_speedbins','Units','pixels','Position',[60 60 1500 520],'Color','w');
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Position-derived velocity: RMSE within speed bins (bins of %s)', bl), 'FontWeight','bold');

panelData  = {rmseBinMax, rmseBinDop, relBinDop};
panelY     = {sprintf('RMSE vs %d s position velocity (m/s)', windows(end)), ...
              'RMSE vs 1 Hz Doppler speed (m/s)', ...
              'RMSE vs Doppler as % of mean Doppler speed'};
panelTitle = {sprintf('vs %d s velocity', windows(end)), 'vs 1 Hz Doppler speed', 'Relative to speed'};
for pnl = 1:3
    ax = nexttile; hold on
    for w = 1:nW
        if pnl == 1 && w == nW, continue, end        % 256 s vs itself is 0
        plot(1:nBins, panelData{pnl}(:,w), '-o', 'Color',cols(w,:), 'LineWidth',1.8, ...
            'MarkerFaceColor',cols(w,:), 'DisplayName', sprintf('%d s', windows(w)));
    end
    set(ax,'XTick',1:nBins,'XTickLabel',xl);
    xlim([0.7 nBins+0.3]); ylim([0 inf]); grid on; box on
    xlabel('Drifter speed (m/s)'); ylabel(panelY{pnl}); title(panelTitle{pnl});
    if pnl == 3
        lg = legend('Location','eastoutside');
        title(lg,'Window');
    end
end
savefig(f2, fullfile(fig_dir,'VeloQAQC_RMSE_speedbins.fig'));
exportgraphics(f2, fullfile(fig_dir,'VeloQAQC_RMSE_speedbins.png'), 'Resolution',300);

%% Printed tables
fprintf('\n=== RMSE of position-derived velocity by window (%d tracks), m/s ===\n', nTrk);
fprintf('%-8s %-28s %-28s\n', 'window', 'vs 256 s velocity [95%]', 'vs 1 Hz Doppler [95%]');
for w = 1:nW
    if w < nW
        a = sprintf('%.4f [%.4f-%.4f]', mM(w), loM(w), hiM(w));
    else
        a = '-';
    end
    fprintf('%-8s %-28s %.4f [%.4f-%.4f]\n', sprintf('%d s',windows(w)), a, mD(w), loD(w), hiD(w));
end

fprintf('\n=== RMSE vs 1 Hz Doppler speed within speed bins (m/s); bins of %s ===\n', bl);
fprintf('%-10s %7s %9s', 'bin (m/s)', 'tracks', 'samples');
wl = arrayfun(@(x) sprintf('%d s',x), windows, 'UniformOutput', false);
fprintf(' %8s', wl{:});
fprintf('\n');
for b = 1:nBins
    fprintf('%-10s %7d %9d', vel_labels{b}, nTrkBin(b,1), nSamp(b,1));
    fprintf(' %8.4f', rmseBinDop(b,:));
    fprintf('\n');
end

fprintf('\n=== RMSE vs %d s position velocity within speed bins (m/s) ===\n', windows(end));
fprintf('%-10s %7s %9s', 'bin (m/s)', 'tracks', 'samples');
fprintf(' %8s', wl{1:end-1});
fprintf('\n');
for b = 1:nBins
    fprintf('%-10s %7d %9d', vel_labels{b}, nTrkBin(b,1), nSamp(b,1));
    fprintf(' %8.4f', rmseBinMax(b,1:end-1));
    fprintf('\n');
end

%% Save
out = struct('windows_s',windows, 'nTracks',nTrk, 'trackCluster',trkCl, 'trackPeriod',trkPer, 'trackID',trkID, ...
    'rmseVs256s',rmseMax, 'rmseVsDoppler',rmseDop, ...
    'meanRMSE_vs256s',mM, 'CI_vs256s',[loM; hiM], 'meanRMSE_vsDoppler',mD, 'CI_vsDoppler',[loD; hiD], ...
    'speedBinEdges',vel_edges, 'binBy',binBy, ...
    'rmseByBin_vs256s',rmseBinMax, 'rmseByBin_vsDoppler',rmseBinDop, 'relRmseByBin_vsDoppler',relBinDop, ...
    'samplesByBin',nSamp, 'tracksByBin',nTrkBin, 'created',datestr(now));
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