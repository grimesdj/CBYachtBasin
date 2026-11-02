clear all;
close all;
%% Drifter trajectories, all observations (left) and moving observations only (right)
% Colour: by tidal phase (D.TidePhase, one label per 1 Hz sample) or by observation period.
% Moving = samples not flagged Stuck (field Stuck in drifters_level_2, from flagStuckDrifters).
% Stuck samples are set to NaN in Lat/Lon, so each track simply has gaps where the drifter
% was stuck. If a level 2 file has no Stuck field yet, it is computed here on the fly
% (needs flagStuckDrifters.m on the path); nothing is saved back to the file.

%% Settings
main_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\';   % where the CB_YachtBasin shared DATA folder is
save_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES';

colorBy  = 'phase';      % 'phase'  = colour by tidal phase of each observation (Rising/High/Falling/Low)
                         % 'period' = colour by observation period
if strcmp(colorBy,'phase')
    fig_name = 'All_Drifter_Trajectories_AllVsMoving_byPhase';
else
    fig_name = 'All_Drifter_Trajectories_AllVsMoving';
end

% tidal phase names and colours (same as the tidal phase plot in findTidalPhaseandRegion)
phaseNames  = {'Rising','High','Falling','Low'};
phaseColors = {'b','r','g','m'};

% display name, folder, file tag, colour
periods = {
    'MarshMadness',            'MarshMadness', 'MarshMadness',      'w'
    'MarshMayhem',             'MarshMayhem',  'MarshMayhem',       'r'
    'FallFrolic - Oct. 08 25', 'FallFrolic',   'FallFrolic_100825', 'b'
    'FallFrolic - Oct. 10 25', 'FallFrolic',   'FallFrolic_101025', 'y'
    'NovDep',                  'NovDep',       'NovDep',            'm'
    'JuneJamboree',            'JuneJamboree', 'JuneJamboree',      'g'
    'SummerFlood',             'SummerFlood',  'SummerFlood',       'k'
    'FallFlood',               'FallFlood',    'FallFlood',         'c'};

%% Get all drifter data (from all IOPs) and store in a structure
data = struct('name',{},'data',{},'color',{});
fprintf('\n%-26s %10s %10s %8s\n', 'period', 'samples', 'stuck', 'stuck %');
for p = 1:size(periods,1)
    fn = fullfile(main_dir,'CB_YachtBasin',periods{p,2},'PROCESSED_DATA','DRIFTERS', ...
        [periods{p,3} '_drifters_level_2.mat']);
    if ~isfile(fn)
        fprintf('%-26s file not found, skipped\n', periods{p,1});
        continue
    end
    S = load(fn);
    drifters = S.drifters_level_2;

    if ~isfield(drifters,'Stuck')
        fprintf('%s: no Stuck field in level 2, computing it now\n', periods{p,1});
        [drifters, ~] = flagStuckDrifters(struct('name',periods{p,3}), drifters);
    end

    nAll = 0;  nStuck = 0;
    for j = 1:numel(drifters)
        if isempty(drifters(j).Lat), continue, end
        st = drifters(j).Stuck;
        nAll = nAll + numel(drifters(j).Lat);
        nStuck = nStuck + sum(st);
    end
    fprintf('%-26s %10d %10d %7.1f%%\n', periods{p,1}, nAll, nStuck, 100*nStuck/max(nAll,1));

    k = numel(data) + 1;
    data(k).name  = periods{p,1};
    data(k).data  = drifters;
    data(k).color = periods{p,4};
end

%% Time to plot it all
fig = figure('Units','pixels','Position',[60 60 1500 700],'Color','w');
tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

gx1 = geoaxes(tl);  gx1.Layout.Tile = 1;
plotTracks(gx1, data, false, colorBy, phaseNames, phaseColors);
title(gx1,'All observations');

gx2 = geoaxes(tl);  gx2.Layout.Tile = 2;
plotTracks(gx2, data, true, colorBy, phaseNames, phaseColors);
title(gx2,'Moving observations (stuck samples removed)');

geobasemap(gx1,'satellite');
geobasemap(gx2,'satellite');

% same map extent in both panels
[latlim, lonlim] = geolimits(gx1);
geolimits(gx2, latlim, lonlim);

legend(gx2, 'Location','best');

% Save as PNG
exportgraphics(fig, fullfile(save_dir, [fig_name '.png']), 'Resolution',300);

% Save as MATLAB figure
savefig(fig, fullfile(save_dir, [fig_name '.fig']));


%% ========================= Local functions =========================

function plotTracks(gx, data, movingOnly, colorBy, phaseNames, phaseColors)
% All tracks are collected into one line per colour group (NaN between tracks), so the map
% draws quickly and the legend has one entry per group.
hold(gx,'on')
if strcmp(colorBy,'phase')
    names = phaseNames;  cols = phaseColors;
else
    names = {data.name};  cols = {data.color};
end
nG  = numel(names);
LAT = repmat({[]},nG,1);  LON = LAT;

for i = 1:numel(data)
    drifters = data(i).data;
    for j = 1:numel(drifters)
        lat = drifters(j).Lat(:);
        lon = drifters(j).Lon(:);
        if isempty(lat), continue, end
        n = min(numel(lat), numel(lon));

        if movingOnly
            st = drifters(j).Stuck(:);
            n  = min(n, numel(st));
            lat(st(1:n)) = NaN;  lon(st(1:n)) = NaN;     % stuck samples left out
        end

        if strcmp(colorBy,'phase')
            ph = string(drifters(j).TidePhase(:));
            if isempty(ph), continue, end
            n = min(n, numel(ph));
            for g = 1:nG
                m = ph(1:n) == phaseNames{g};
                m = m | [m(2:end); false];               % include the next point so segments join
                la = lat(1:n);  lo = lon(1:n);
                la(~m) = NaN;   lo(~m) = NaN;
                LAT{g} = [LAT{g}; la; NaN];              %#ok<AGROW>
                LON{g} = [LON{g}; lo; NaN];              %#ok<AGROW>
            end
        else
            LAT{i} = [LAT{i}; lat(1:n); NaN];            %#ok<AGROW>
            LON{i} = [LON{i}; lon(1:n); NaN];            %#ok<AGROW>
        end
    end
end

for g = 1:nG
    if isempty(LAT{g}), LAT{g} = NaN;  LON{g} = NaN; end
    geoplot(gx, LAT{g}, LON{g}, 'Color',cols{g}, 'LineWidth',1.5, 'DisplayName',names{g});
end
end
