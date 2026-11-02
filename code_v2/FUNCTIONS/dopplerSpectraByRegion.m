function Spec = dopplerSpectraByRegion(data, save_dir)
% Doppler speed variance spectra of the MOVING drifters, one per basin region.
%
% Call from allVeloAnalysis (just before "all_RMSE = data;"):
%       dopplerSpectraByRegion(data);
% "data" is the struct allVeloAnalysis already builds (one field per observation
% period, each holding the TemporalRes_RMSE Results array). Needs the fields
% Speed (1 Hz Doppler speed, m/s), Region and Time in every Results entry.
%
% How it works (plain language):
%   A spectrum needs an unbroken stretch of data. A drifter passes through several
%   regions, so each track is cut into short pieces (segLen seconds). A piece is
%   used ONLY if the whole piece lies inside one region, has no gaps and no NaNs.
%   Each piece gets its own spectrum, and all pieces from a region are averaged
%   (Welch's method). The shaded band is a 95% interval from resampling whole
%   tracks (the same idea as the regional maps), so it reflects track-to-track
%   differences, not just the number of pieces.
%   The grey line is the same calculation for the STATIONARY drifters, i.e. what a
%   drifter that is not moving "sees" (GPS noise). Where a region's spectrum sits
%   above it, the variance is real motion.
%
% Limits to keep in mind:
%   - With segLen = 256 s the longest period you can see is 256 s (4 min), and
%     the lowest frequency is 1/256 Hz. Longer segLen resolves slower motion but
%     leaves fewer usable pieces per region.
%   - Doppler speed is a magnitude (never negative), so this is the spectrum of
%     speed, not of velocity components.
%   - Level 1 interpolates Doppler to 1 Hz, so the very highest frequencies are
%     partly an interpolation effect.

if nargin < 2 || isempty(save_dir)
    save_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES\SPECTRA';
end
if ~isfolder(save_dir), mkdir(save_dir); end

%% Settings
fs         = 1;            % Hz
segLen     = 256;          % s, length of each piece (must be even)
overlap    = 0.5;          % fraction of overlap between pieces
regionIDs  = 1:7;
minTracks  = 3;            % need at least this many tracks to plot a region
minSegs    = 5;            % ... and at least this many pieces
nBoot      = 1000;
bandEdges  = [256 32 8 2]; % s, for the variance split printed below

step = round(segLen*(1-overlap));
nf   = segLen/2 + 1;
f    = (0:segLen/2)' * fs/segLen;
df   = fs/segLen;
w    = 0.5*(1 - cos(2*pi*(0:segLen-1)'/(segLen-1)));   % Hann window
U    = sum(w.^2);
rng(1);

%% Gather pieces from the moving drifters, track by track
nR = numel(regionIDs);
trkSum    = cell(nR,1);     % each: nTracks x nf  (sum of piece spectra per track)
trkN      = cell(nR,1);     % each: nTracks x 1   (number of pieces per track)
trkPeriod = cell(nR,1);
for r = 1:nR
    trkSum{r} = zeros(0,nf);  trkN{r} = zeros(0,1);  trkPeriod{r} = zeros(0,1);
end

obs_names = fieldnames(data);
for i = 1:numel(obs_names)
    R = data.(obs_names{i});

    if ~all(isfield(R, {'Speed','Region','Time'}))
        fprintf('%s: Results has no Speed/Region/Time fields - skipped\n', obs_names{i});
        continue
    end

    for d = 1:size(R,1)
        for j = 1:size(R,2)
            x   = R(d,j).Speed;
            reg = R(d,j).Region;
            tim = R(d,j).Time;
            if isempty(x) || isempty(reg) || isempty(tim), continue, end

            x = x(:);  reg = reg(:);
            if numel(x) ~= numel(reg) || numel(x) ~= numel(tim)
                fprintf('%s d%d j%d: Speed/Region/Time lengths differ - skipped\n', obs_names{i}, d, j);
                continue
            end
            t   = toSeconds(tim);
            brk = [false; abs(diff(t) - 1) > 0.05];     % sample not 1 s after previous

            for r = 1:nR
                ok = ~isnan(x) & (reg == regionIDs(r));
                [S, n] = pieceSpectra(x, ok, brk, segLen, step, w, U, fs, nf);
                if n > 0
                    trkSum{r}(end+1,:)    = S;            %#ok<AGROW>
                    trkN{r}(end+1,1)      = n;            %#ok<AGROW>
                    trkPeriod{r}(end+1,1) = i;            %#ok<AGROW>
                end
            end
        end
    end
end

%% Stationary reference spectrum (same method)
statPSD = [];
statN   = 0;
roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if ~isempty(k)
    sf = fullfile(roots{k},'CB_YachtBasin','StationaryDep','PROCESSED_DATA','DRIFTERS', ...
        'StationaryDep_drifters_level_1.mat');
    if isfile(sf)
        L1 = load(sf).drifters_level_1;
        sS = zeros(1,nf);
        for d = 1:size(L1,1)
            for j = 1:size(L1,2)
                D = L1(d,j);
                if isempty(D.Speed) || isempty(D.e_time), continue, end
                x   = D.Speed(:)*0.51444;                  % knots -> m/s
                t   = toSeconds(D.e_time);
                brk = [false; abs(diff(t) - 1) > 0.05];
                [S, n] = pieceSpectra(x, ~isnan(x), brk, segLen, step, w, U, fs, nf);
                sS = sS + S;  statN = statN + n;
            end
        end
        if statN > 0, statPSD = (sS/statN)'; end
    end
end
if isempty(statPSD)
    fprintf('No stationary reference spectrum (file missing or no usable pieces).\n');
end

%% Average, bootstrap, summarise
Reg = struct('Region',num2cell(regionIDs), 'nTracks',0, 'nPeriods',0, 'nPieces',0, ...
    'hours',0, 'PSD',[], 'lo',[], 'hi',[], 'bandFrac',[], 'rmsFluct',NaN);

for r = 1:nR
    nTr = size(trkSum{r},1);
    Reg(r).nTracks  = nTr;
    Reg(r).nPeriods = numel(unique(trkPeriod{r}));
    Reg(r).nPieces  = sum(trkN{r});
    Reg(r).hours    = Reg(r).nPieces * (segLen*(1-overlap)) / 3600;   % approx. independent data
    if nTr == 0, continue, end

    Reg(r).PSD = (sum(trkSum{r},1) / sum(trkN{r}))';

    if nTr >= minTracks
        B = zeros(nBoot,nf);
        for b = 1:nBoot
            ix = randi(nTr, nTr, 1);
            B(b,:) = sum(trkSum{r}(ix,:),1) / sum(trkN{r}(ix));
        end
        B = sort(B,1);
        Reg(r).lo = B(max(1,round(0.025*nBoot)),:)';
        Reg(r).hi = B(min(nBoot,round(0.975*nBoot)),:)';
    end

    P = Reg(r).PSD;
    Reg(r).rmsFluct = sqrt(sum(P(2:end))*df);
    Reg(r).bandFrac = nan(1,numel(bandEdges)-1);
    for b = 1:numel(bandEdges)-1
        loF = 1/bandEdges(b) - 1e-12;
        hiF = 1/bandEdges(b+1) - 1e-12;
        if b == numel(bandEdges)-1, hiF = inf; end       % last band runs to Nyquist
        inBand = f > 0 & f >= loF & f < hiF;             % half-open: no sample counted twice
        Reg(r).bandFrac(b) = sum(P(inBand))*df / (sum(P(2:end))*df);
    end
end

%% Print
fprintf('\n=== Doppler speed spectra by region (segLen %d s, %.0f%% overlap) ===\n', segLen, 100*overlap);
fprintf('%-7s %7s %8s %8s %12s | variance share by period:', 'Region','tracks','periods','pieces','RMS fluct');
for b = 1:numel(bandEdges)-1
    fprintf(' %d-%d s', bandEdges(b), bandEdges(b+1));
    if b < numel(bandEdges)-1, fprintf(' |'); end
end
fprintf('\n');
for r = 1:nR
    if Reg(r).nTracks == 0
        fprintf('%-7d %7d %8d %8d %12s\n', regionIDs(r), 0, 0, 0, '-');
        continue
    end
    fprintf('%-7d %7d %8d %8d %9.2f cm/s |', regionIDs(r), Reg(r).nTracks, Reg(r).nPeriods, ...
        Reg(r).nPieces, 100*Reg(r).rmsFluct);
    fprintf(' %5.0f%%     ', 100*Reg(r).bandFrac);
    fprintf('\n');
end
if ~isempty(statPSD)
    sp = sqrt(sum(statPSD(2:end))*df);
    fprintf('Stationary reference: %d pieces, RMS fluctuation %.2f cm/s (instrument noise level)\n', statN, 100*sp);
end
fprintf('(RMS fluct = square root of the variance in the spectrum, 4 min down to 2 s;\n');
fprintf(' regions with fewer than %d tracks or %d pieces are not drawn)\n', minTracks, minSegs);

%% Figure 1: one panel per region
cols = lines(nR);
f1 = figure('Name','DopplerSpectra_byRegion','Units','pixels','Position',[50 50 1500 800],'Color','w');
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
title(tl, {'Doppler speed variance spectra by region (moving drifters)', ...
    sprintf('%d s pieces inside a single region; shading = 95%% interval (bootstrap by track); grey = stationary drifters', segLen)}, ...
    'FontWeight','bold','FontSize',12);
xlabel(tl,'Frequency (Hz)'); ylabel(tl,'PSD of Doppler speed [(m/s)^2/Hz]');

axs = gobjects(nR+1,1);
allVals = [];
for r = 1:nR
    axs(r) = nexttile; hold on
    set(axs(r),'XScale','log','YScale','log');
    plotted = Reg(r).nTracks >= minTracks && Reg(r).nPieces >= minSegs && ~isempty(Reg(r).lo);
    if ~isempty(statPSD)
        plot(f(2:end), statPSD(2:end), '-', 'Color',[0.6 0.6 0.6], 'LineWidth',1.5);
    end
    if plotted
        fill([f(2:end); flipud(f(2:end))], [Reg(r).lo(2:end); flipud(Reg(r).hi(2:end))], ...
            cols(r,:), 'FaceAlpha',0.25, 'EdgeColor','none');
        plot(f(2:end), Reg(r).PSD(2:end), '-', 'Color',cols(r,:), 'LineWidth',2);
        allVals = [allVals; Reg(r).lo(2:end); Reg(r).hi(2:end)]; %#ok<AGROW>
        title(sprintf('Region %d  (%d tracks, %d pieces)', regionIDs(r), Reg(r).nTracks, Reg(r).nPieces));
    else
        title(sprintf('Region %d  (too few data: %d tracks, %d pieces)', regionIDs(r), Reg(r).nTracks, Reg(r).nPieces));
    end
    grid on; box on
    xlim([f(2) fs/2]);
end

% 8th panel: all regions together
axs(nR+1) = nexttile; hold on
set(axs(nR+1),'XScale','log','YScale','log');
if ~isempty(statPSD)
    hS = plot(f(2:end), statPSD(2:end), '-', 'Color',[0.6 0.6 0.6], 'LineWidth',1.5);
end
hh = gobjects(0);  ll = {};
for r = 1:nR
    if Reg(r).nTracks >= minTracks && Reg(r).nPieces >= minSegs
        hh(end+1) = plot(f(2:end), Reg(r).PSD(2:end), '-', 'Color',cols(r,:), 'LineWidth',1.8); %#ok<AGROW>
        ll{end+1} = sprintf('Region %d', regionIDs(r));                                         %#ok<AGROW>
    end
end
if ~isempty(statPSD), hh(end+1) = hS; ll{end+1} = 'Stationary'; end
if ~isempty(hh), legend(hh, ll, 'Location','southwest'); end
title('All regions'); grid on; box on; xlim([f(2) fs/2]);

% same y-limits everywhere so panels can be compared
yl = [];
if ~isempty(allVals) && ~isempty(statPSD)
    yl = [min([allVals; statPSD(2:end)]) max([allVals; statPSD(2:end)])];
elseif ~isempty(allVals)
    yl = [min(allVals) max(allVals)];
end
if ~isempty(yl), set(axs, 'YLim', yl); end

%% Save
exportgraphics(f1, fullfile(save_dir,'DopplerSpectra_byRegion.png'), 'Resolution',300);
savefig(f1, fullfile(save_dir,'DopplerSpectra_byRegion.fig'));

Spec = struct('f',f, 'segLen_s',segLen, 'overlap',overlap, 'regionIDs',regionIDs, ...
    'Region',Reg, 'StationaryPSD',statPSD, 'StationaryPieces',statN, ...
    'bandEdges_s',bandEdges, 'created',datestr(now));
save(fullfile(save_dir,'DopplerSpectra_byRegion.mat'), 'Spec');
fprintf('\nSaved figure and results to %s\n', save_dir);
end


%% ========================= Local functions =========================

function [S, n] = pieceSpectra(x, ok, brk, segLen, step, w, U, fs, nf)
% Sum of spectra of every piece (segLen samples) that is entirely inside the
% valid set (ok = true), with no time gaps (brk = true marks a gap).
S = zeros(1,nf);
n = 0;
N = numel(x);
if N < segLen, return, end

badCum = cumsum(~ok);            % invalid samples up to and including i
brkCum = cumsum(brk);            % gaps up to and including i

for i0 = 1:step:(N - segLen + 1)
    i1 = i0 + segLen - 1;
    nBad = badCum(i1) - (badCum(i0) - ~ok(i0));
    nBrk = brkCum(i1) - brkCum(i0);          % gap at the first sample itself does not matter
    if nBad > 0 || nBrk > 0, continue, end

    seg = detrend(x(i0:i1), 1);              % remove mean and linear trend
    X   = fft(seg .* w);
    P   = abs(X(1:nf)).^2 / (fs*U);
    P(2:end-1) = 2*P(2:end-1);               % one-sided
    S = S + P(:)';
    n = n + 1;
end
end


function t = toSeconds(T)
% Time vector in seconds from start, whatever the stored format
T = T(:);
if isdatetime(T) || isduration(T)
    t = seconds(T - T(1));
else
    t = double(T);
    if median(diff(t),'omitnan') < 0.01       % datenum (days)
        t = (t - t(1)) * 86400;
    else
        t = t - t(1);
    end
end
end
