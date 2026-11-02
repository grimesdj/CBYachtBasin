%% regionalVelocityMaps_v2
% Regional MEAN SPEED (scalar) and MEAN VELOCITY (vector) maps from drifters,
% for all samples and for each tidal phase, from one loader and one set of
% statistics. Does NOT replace regionalVeloAnalysis.m or RegionalVeloVectors.m.
%
% Choices to record in the decisions log:
%   temporal_res = 8 s. Stationary-drifter noise floor at 8 s = 0.011 m/s (9 drifters).
%   Uncertainty = 95% bootstrap interval that resamples whole deployments
%     (observation period x deployment), using only deployments that have data
%     in that region/phase.
%   Arrows: solid  = east or north component of the mean has a 95% interval that
%                    excludes zero (clearly different from zero)
%           dashed = cannot rule out zero mean flow
%   Cells with fewer than min_tracks tracks or min_periods periods are not drawn.
%   "Rising"/"Falling" are pressure-based phases, not flood/ebb.

clear; close all; rng(1)

%% Settings
temporal_res   = 8;      % options: [1 2 4 8 16 32 64 128 256]
min_tracks     = 3;      % cells with fewer tracks are not drawn
min_periods    = 2;      % cells from fewer observation periods are not drawn
nBoot          = 2000;
velocity_scale = 3000;   % arrow length in metres per 1 m/s
ref_speed      = 0.05;   % reference arrow, m/s

roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', ...
         'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if isempty(k), error('Neither root path contains CB_YachtBasin on this machine.'); end
root_data = roots{k};
fprintf('Using root: %s\n', root_data);

save_dir = fullfile(root_data,'THESIS','THESIS_DRAFTING','FIGURES','REGIONAL_ANALYSIS');
if ~isfolder(save_dir), mkdir(save_dir); end

names   = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};
folders = {'MarshMadness','MarshMayhem','FallFrolic','FallFrolic','NovDep','JuneJamboree','SummerFlood'};
phases  = ["High","Low","Rising","Falling"];

%% Load every track at the chosen temporal resolution
RD = struct('Observation',{},'Period',{},'Deployment',{},'ID',{},'Cluster',{}, ...
            'V',{},'V_e',{},'V_n',{},'Region',{},'TidePhase',{});

for i = 1:numel(names)

    fn = fullfile(root_data,'CB_YachtBasin',folders{i},'PROCESSED_DATA','DRIFTERS', ...
        [names{i} '_TemporalRes_RMSE.mat']);

    if ~isfile(fn)
        warning('%s: file not found, skipped', names{i});
        continue
    end

    Results = load(fn).Results;

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

            rec = struct('Observation',names{i}, 'Period',i, 'Deployment',d, ...
                'ID',Results(d,j).ID, 'Cluster',i*100 + d, ...
                'V',V(:), 'V_e',Ve(:), 'V_n',Vn(:), ...
                'Region',Reg(:), 'TidePhase',string(Ph(:)));

            RD(end+1) = rec; %#ok<SAGROW>

        end
    end
end

fprintf('\nLoaded %d tracks at %d s resolution.\n', numel(RD), temporal_res);

%% Statistics for every region, all phases and each phase
load('regions.mat');
nR = numel(regions);

for r = 1:nR
    All(r) = cellStats(RD, r, "", nBoot); %#ok<SAGROW>
end

groups(1).label = 'all tidal phases';
groups(1).tag   = 'All';
groups(1).S     = All;

for p = 1:numel(phases)
    ph = char(phases(p));
    for r = 1:nR
        PS.(ph)(r) = cellStats(RD, r, phases(p), nBoot);
    end
    groups(p+1).label = [ph ' phase'];
    groups(p+1).tag   = ph;
    groups(p+1).S     = PS.(ph);
end

%% Print tables
for g = 1:numel(groups)
    fprintf('\n=== %s (%d s) ===\n', groups(g).label, temporal_res);
    printTable(groups(g).S);
end

%% Maps
for g = 1:numel(groups)

    drawMap(regions, groups(g).S, 'speed', ...
        sprintf('Mean drifter speed by region, %s (%d s)', groups(g).label, temporal_res), ...
        sprintf('RegionalMeanScalarVelo_v2_%s_%ds', groups(g).tag, temporal_res), ...
        save_dir, min_tracks, min_periods, velocity_scale, ref_speed);

    drawMap(regions, groups(g).S, 'vector', ...
        sprintf('Mean drifter velocity by region, %s (%d s)', groups(g).label, temporal_res), ...
        sprintf('RegionalMeanVelocityVectors_v2_%s_%ds', groups(g).tag, temporal_res), ...
        save_dir, min_tracks, min_periods, velocity_scale, ref_speed);

end

%% Save the statistics so figures/tables do not need a rerun
save(fullfile(save_dir, sprintf('RegionalStats_v2_%ds.mat', temporal_res)), ...
    'groups','temporal_res','min_tracks','min_periods','nBoot');
fprintf('\nSaved RegionalStats_v2_%ds.mat to %s\n', temporal_res, save_dir);


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

% Bootstrap over deployments
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


function drawMap(regions, S, mode, ttl, file_base, save_dir, min_tracks, min_periods, vscale, ref_speed)

nR = numel(regions);

figure('Units','pixels','Position',[100 100 1400 850], 'Color','w');
gx = geoaxes;
hold(gx,'on')
geobasemap(gx,'satellite')

h   = gobjects(nR,1);
lab = strings(nR,1);

for r = 1:nR

    lat = regions(r).Lat(:);
    lon = regions(r).Lon(:);

    h(r) = geoplot(gx, lat, lon, 'LineWidth', 2);

    latc = mean(lat(1:end-1),'omitnan');
    lonc = mean(lon(1:end-1),'omitnan');

    ok = S(r).NTracks >= min_tracks && S(r).NPeriods >= min_periods;

    if ~ok
        lab(r) = sprintf('Region %d: too few data  (%d tracks, %d periods)', ...
            r, S(r).NTracks, S(r).NPeriods);
    elseif strcmp(mode,'speed')
        lab(r) = sprintf('Region %d: %.3f m/s [%.3f-%.3f]  (%d tracks, %d periods)', ...
            r, S(r).MeanSpeed, S(r).SpeedCI(1), S(r).SpeedCI(2), S(r).NTracks, S(r).NPeriods);
    else
        if S(r).Differs, ls = '-'; else, ls = '--'; end
        drawArrow(gx, latc, lonc, S(r).Ve, S(r).Vn, vscale, ls);
        lab(r) = sprintf('Region %d: |mean V| = %.3f m/s  (%d tracks, %d periods)', ...
            r, S(r).VecMag, S(r).NTracks, S(r).NPeriods);
    end

    text(gx, latc, lonc, sprintf('%d',r), ...
        'FontSize',14, 'FontWeight','bold', 'Color','w', ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom');

end

if strcmp(mode,'vector')
    allLat = cell2mat(arrayfun(@(x) x.Lat(:), regions(:), 'UniformOutput', false));
    allLon = cell2mat(arrayfun(@(x) x.Lon(:), regions(:), 'UniformOutput', false));
    lat_ref = min(allLat) - 0.0015;
    lon_ref = min(allLon);
    drawArrow(gx, lat_ref, lon_ref, ref_speed, 0, vscale, '-');
    text(gx, lat_ref - 0.0004, lon_ref, sprintf('%.2f m/s', ref_speed), ...
        'Color','w', 'FontWeight','bold', 'HorizontalAlignment','left');
    sub = {'solid arrow = clearly different from zero;', 'dashed = cannot rule out zero (95%, bootstrap by deployment)'};
else
    sub = {'brackets = 95% interval, bootstrap by deployment'};
end

legend(gx, h, lab, 'Location','bestoutside');
sgtitle([{ttl}; sub(:)], 'FontSize',13, 'FontWeight','bold');

exportgraphics(gcf, fullfile(save_dir,[file_base '.png']), 'Resolution',300);
savefig(gcf, fullfile(save_dir,[file_base '.fig']));

end


function drawArrow(gx, lat0, lon0, Ve, Vn, vscale, ls)
% Arrow from (lat0,lon0); length = speed * vscale metres.

mlat = 111320;
mlon = 111320*cosd(lat0);

dx = Ve*vscale;
dy = Vn*vscale;
L  = hypot(dx,dy);

if ~(L > 0), return, end

lat1 = lat0 + dy/mlat;
lon1 = lon0 + dx/mlon;

geoplot(gx, [lat0 lat1], [lon0 lon1], 'Color','w', 'LineStyle',ls, 'LineWidth',2.5);

ux = dx/L;  uy = dy/L;
px = -uy;   py = ux;

hl = 0.20*L;
hw = 0.10*L;

x1 = dx - hl*ux + hw*px;   y1 = dy - hl*uy + hw*py;
x2 = dx - hl*ux - hw*px;   y2 = dy - hl*uy - hw*py;

geoplot(gx, [lat0 + y1/mlat, lat1, lat0 + y2/mlat], ...
            [lon0 + x1/mlon, lon1, lon0 + x2/mlon], ...
        'Color','w', 'LineStyle','-', 'LineWidth',2.5);

end
