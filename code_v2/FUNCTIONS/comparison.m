function [ROI_ADCP, ROI_EXTRAP, ROI_EOF, COMP] = comparison(cfg, adcpLoc)
%COMPARISON  Level 3: drifter vs ADCP near-surface velocity inside the ROI.
%
% Revised 2026-10-13. Changes from the previous version are tagged [FIX n]:
%   [FIX 1] ROI is the North ADCP only; non-North runs are skipped
%   [FIX 2] Products now mean:
%             ADCP   = raw QC'd bins, highest reliable bin
%             EXTRAP = extrapolated observed profile (Level 2 U_mapped),
%                      which includes the surface value
%             EOF    = EOF reconstruction, mean + leading modes (U_recon)
%           Previously "EOF" was actually the mapped observations.
%   [FIX 3] Target depth uses Level 2 water depth h, not raw pressure
%   [FIX 4] Speed = magnitude of the mean velocity, for drifter and ADCP
%           alike (mean of 1 s GPS speeds is biased high by noise)
%   [FIX 5] Statistics per product, NaN-safe, paired by crossing
%           (previously all three products were pooled into one r/RMSE)
%   [FIX 6] ADCP ensembles averaged over a fixed-length window centred on
%           the drifter transit midpoint (cfg.adcp_window_min), widened
%           only if the transit itself is longer; max time-offset check
%   [FIX 7] Raw ADCP top bin = highest bin below the surface that is
%           valid in >= 50% of the transit (was: first sample only)
%   [FIX 8] "No drifters in ROI" check works when some slots are empty
%   [FIX 9] Figure file names include obsTag (no more overwriting)
%
% COMP layout (changed):
%   COMP.ADCP / COMP.EXTRAP / COMP.EOF  -> .V, .VE, .VN each with
%       n, r, bias (drifter - ADCP), rmse, mean_abs_drifter, mean_abs_adcp
%   COMP.pairs -> paired values per crossing (j, i, drifter_*, ADCP_*, ...)
%
% Optional cfg fields:
%   cfg.adcp_latlon          [lat lon] of ROI centre (default: North ADCP)
%   cfg.roi_radius           m        (default 50)
%   cfg.drifter_depth        m below surface compared to drifter (default 0.5)
%   cfg.max_time_offset_min  min      (default 15)
%   cfg.adcp_window_min      total length of the ADCP averaging window,
%                            centred on the transit midpoint (default 20)
%   cfg.adcp_min_valid_frac  fraction of transit the raw top bin must be
%                            valid in (default 0.5; use eps for "any")

ROI_ADCP = []; ROI_EXTRAP = []; ROI_EOF = []; COMP = [];

%% [FIX 1] North ADCP only for now
if ~contains(cfg.obsTag, 'NorthADCP')
    warning('comparison:northOnly', ...
        'ROI is defined for the North ADCP only; skipping comparison for %s.', cfg.obsTag);
    return
end

%% Settings
adcp_latlon  = getcfg(cfg, 'adcp_latlon', [34.05673145 -77.88949171]);   % North ADCP
thresh       = getcfg(cfg, 'roi_radius', 50);
z_target     = getcfg(cfg, 'drifter_depth', 0.5);
maxOffset    = getcfg(cfg, 'max_time_offset_min', 15)/1440;              % days
minValidFrac = getcfg(cfg, 'adcp_min_valid_frac', 0.5);   % top bin must be valid in this fraction of the transit
winHalf      = getcfg(cfg, 'adcp_window_min', 20)/2/1440;  % half-window, days

%% Load
base_name = regexprep(cfg.name, '_sheet$', '');
% L1 = load(fullfile(cfg.out.adcp_data, [cfg.obsTag '_adcp_level_1.mat']), 'adcp_level_1');
L1 = load(fullfile(cfg.out.adcp_data, [base_name '_' adcpLoc '_adcp_level_1.mat']));

adcp_level_1 = L1.adcp_level_1;
% L2 = load(fullfile(cfg.out.adcp_data, [cfg.obsTag '_adcp_level_2.mat']), 'adcp_level_2');
L2 = load(fullfile(cfg.out.adcp_data, [base_name '_' adcpLoc '_adcp_level_2.mat']));
adcp_level_2 = L2.adcp_level_2;

S = load(fullfile(cfg.out.drifters_data, [cfg.name '_drifters_level_2.mat']));
if isfield(S, 'drifters_level_2')
    drifters_level_2 = S.drifters_level_2;
elseif isfield(S, 'trimmed')
    drifters_level_2 = S.trimmed;
else
    error('File does not contain drifters_level_2 or trimmed')
end

if ~isfield(adcp_level_2, 'U_recon') || ~isfield(adcp_level_2, 'h')
    error(['comparison: Level 2 for %s was made by the old adcp_post_processing ' ...
           '(no U_recon/h). Re-run Level 2 first.'], cfg.obsTag);
end

proj   = projcrs(32119);
time   = adcp_level_1.time(:);
h      = adcp_level_2.h(:);
qcFlag = adcp_level_1.qcFlag';          % bins x time
Rm     = adcp_level_2.rot_mat;

%% Products in east/north, bins x time  [FIX 2]
[ve_x, vn_x] = unrotate(adcp_level_2.U_mapped, adcp_level_2.V_mapped, Rm);  % z/h x time
[ve_e, vn_e] = unrotate(adcp_level_2.U_recon,  adcp_level_2.V_recon,  Rm);  % time x z/h

Products(1) = struct('name', 'ADCP',   'VE', adcp_level_1.V_E, 'VN', adcp_level_1.V_N, ...
                     'bins', adcp_level_1.bins(:)', 'vcoord', 'mab');
Products(2) = struct('name', 'EXTRAP', 'VE', ve_x,   'VN', vn_x, ...
                     'bins', adcp_level_2.zonh_grid(:)', 'vcoord', 'zonh');
Products(3) = struct('name', 'EOF',    'VE', ve_e.', 'VN', vn_e.', ...
                     'bins', adcp_level_2.zonh_grid(:)', 'vcoord', 'zonh');

%% Output template
R_template = struct( ...
    'inds', [], 'inde', [], ...
    'East', [], 'North', [], 'Lat', [], 'Lon', [], ...
    'v_e', [], 'v_n', [], 'V', [], ...
    'avg_V', [], 'avg_ve', [], 'avg_vn', [], 'avg_speed_1s', [], ...
    'drifter_Time_dattime', [], 'drifter_Time', [], 'transit_min', [], ...
    'adcp_ve', [], 'adcp_vn', [], 'adcp_V', [], ...
    'last_idx', [], 'V_topbin', [], 'VE_topbin', [], 'VN_topbin', [], ...
    'adcp_P', [], 'adcp_h', [], 'qcflag', [], 'adcp_bin', [], ...
    'adcp_Time_datnum', [], 'adcp_Time_dattime', [], 'adcp_dbins', []);

[numDep, numDr] = size(drifters_level_2);
ROIs = repmat({repmat(R_template, numDep, numDr)}, 1, numel(Products));

[E_loc, N_loc] = projfwd(proj, adcp_latlon(1), adcp_latlon(2));

%% Main loop
for j = 1:numDep
    for i = 1:numDr

        D = drifters_level_2(j,i);
        if isempty(D.East) || isempty(D.North)
            continue
        end

        dist = hypot(D.East - E_loc, D.North - N_loc);
        in   = dist <= thresh;
        if ~any(in)
            continue
        end

        T  = D.Date_and_time(in);
        tD = datenum(T);  tD = tD(:);
        t0 = tD(1);  t1 = tD(end);

        % [FIX 6] Fixed-length window centred on the transit midpoint,
        % widened to cover the transit if the transit is longer
        tMid = 0.5*(t0 + t1);
        hw   = max(winHalf, 0.5*(t1 - t0));
        k = find(time >= tMid - hw & time <= tMid + hw);
        if isempty(k)
            [off, k] = min(abs(time - mean([t0 t1])));
            if off > maxOffset
                fprintf('Dep %d drifter %d: nearest ADCP ensemble is %.0f min from transit; skipped\n', ...
                    j, i, off*1440);
                continue
            end
        end
        inds = k(1);  inde = k(end);

        transit_min = (t1 - t0)*1440;
        if transit_min > 60
            fprintf(['Dep %d drifter %d: ROI points span %.0f min ' ...
                     '(multiple passes?) - averaged together\n'], j, i, transit_min);
        end

        hT = mean(h(inds:inde), 'omitnan');
        zTarget = hT - z_target;                 % height above bed

        % [FIX 4] drifter: magnitude of the mean velocity
        ve_d = D.v_e(in);  vn_d = D.v_n(in);
        avg_ve = mean(ve_d, 'omitnan');
        avg_vn = mean(vn_d, 'omitnan');
        avg_V  = hypot(avg_ve, avg_vn);

        for p = 1:numel(Products)

            R  = R_template;
            VE = Products(p).VE(:, inds:inde);
            VN = Products(p).VN(:, inds:inde);
            VEm = mean(VE, 2, 'omitnan');
            VNm = mean(VN, 2, 'omitnan');
            bins = Products(p).bins(:);

            VE_top = NaN;  VN_top = NaN;  adcp_bin = NaN;  last_idx = NaN;

            switch Products(p).name
                case 'ADCP'    % [FIX 7]
                    frac = mean(~isnan(VE) & ~isnan(VN), 2);
                    cand = find(frac >= minValidFrac & bins < hT);
                    if ~isempty(cand)
                        last_idx = cand(end);
                        VE_top   = VEm(last_idx);
                        VN_top   = VNm(last_idx);
                        adcp_bin = bins(last_idx);
                    end

                otherwise      % EXTRAP, EOF on z/h  [FIX 3]
                    zab = bins*hT;
                    ok  = ~isnan(VEm) & ~isnan(VNm);
                    if nnz(ok) >= 2
                        VE_top = interp1(zab(ok), VEm(ok), zTarget);
                        VN_top = interp1(zab(ok), VNm(ok), zTarget);
                    end
            end

            R.inds = inds;  R.inde = inde;

            R.adcp_ve   = VE;
            R.adcp_vn   = VN;
            R.adcp_V    = hypot(VE, VN);
            R.VE_topbin = VE_top;
            R.VN_topbin = VN_top;
            R.V_topbin  = hypot(VE_top, VN_top);      % [FIX 4]
            R.last_idx  = last_idx;
            R.adcp_bin  = adcp_bin;
            R.adcp_P    = adcp_level_1.pressure(inds:inde);
            R.adcp_h    = h(inds:inde);
            R.qcflag    = qcFlag(:, inds:inde);
            R.adcp_Time_datnum  = time(inds:inde);
            R.adcp_Time_dattime = datetime(time(inds:inde), 'ConvertFrom', 'datenum');
            R.adcp_dbins = Products(p).bins;

            R.East  = D.East(in);   R.North = D.North(in);
            R.Lat   = D.Lat(in);    R.Lon   = D.Lon(in);
            R.v_e   = ve_d;         R.v_n   = vn_d;   R.V = D.V(in);
            R.avg_ve = avg_ve;      R.avg_vn = avg_vn;  R.avg_V = avg_V;
            R.avg_speed_1s = mean(D.V(in), 'omitnan');  % old definition, for reference
            R.drifter_Time_dattime = T;
            R.drifter_Time = tD;
            R.transit_min  = transit_min;

            ROIs{p}(j,i) = R;
        end
    end
end

ROI_ADCP   = ROIs{1};
ROI_EXTRAP = ROIs{2};
ROI_EOF    = ROIs{3};

%% [FIX 8] Any crossings at all?
hasCross = arrayfun(@(r) ~isempty(r.inds), ROI_ADCP);
if ~any(hasCross(:))
    disp('No events where drifters pass within ROI');
    ROI_ADCP = []; ROI_EXTRAP = []; ROI_EOF = []; COMP = [];
    return
end

%% [FIX 5] Paired values and per-product statistics
[cj, ci] = find(hasCross);
pairs.j = cj;
pairs.i = ci;
pairs.drifter_V  = [ROI_ADCP(hasCross).avg_V]';
pairs.drifter_VE = [ROI_ADCP(hasCross).avg_ve]';
pairs.drifter_VN = [ROI_ADCP(hasCross).avg_vn]';

comps = {'V', 'VE', 'VN'};
fprintf('\n%s: %d ROI crossings (ROI radius %g m, %.2f m below surface)\n', ...
    cfg.obsTag, numel(cj), thresh, z_target);
fprintf('%-7s %-3s %4s %7s %9s %8s\n', 'Product', 'Var', 'n', 'r', 'bias', 'RMSE');

for p = 1:numel(Products)
    nm = Products(p).name;
    Rp = ROIs{p}(hasCross);
    pairs.([nm '_V'])  = [Rp.V_topbin]';
    pairs.([nm '_VE']) = [Rp.VE_topbin]';
    pairs.([nm '_VN']) = [Rp.VN_topbin]';
    for c = 1:numel(comps)
        st = compStats(pairs.(['drifter_' comps{c}]), pairs.([nm '_' comps{c}]));
        COMP.(nm).(comps{c}) = st;
        fprintf('%-7s %-3s %4d %7.2f %9.3f %8.3f\n', nm, comps{c}, st.n, st.r, st.bias, st.rmse);
    end
end
fprintf('(bias = drifter - ADCP, m/s)\n\n');

COMP.pairs    = pairs;
COMP.settings = struct('adcp_latlon', adcp_latlon, 'roi_radius', thresh, ...
    'drifter_depth', z_target, 'max_time_offset_min', maxOffset*1440, ...
    'min_valid_frac', minValidFrac, 'adcp_window_min', 2*winHalf*1440);

%% Figures  [FIX 9]
allVals = [pairs.drifter_V; pairs.drifter_VE; pairs.drifter_VN];
for p = 1:numel(Products)
    nm = Products(p).name;
    allVals = [allVals; pairs.([nm '_V']); pairs.([nm '_VE']); pairs.([nm '_VN'])]; %#ok<AGROW>
end
L = ceil(10*max(abs(allVals), [], 'omitnan'))/10;
if isempty(L) || isnan(L) || L == 0, L = 0.4; end

titles = {'ADCP OBSERVATIONS (top reliable bin)', ...
          'EXTRAPOLATED OBSERVATIONS', ...
          sprintf('EOF RECONSTRUCTION (mean + %d modes)', getEOFmodes(adcp_level_2))};
files  = {'adcp', 'extrap', 'eof'};
mk     = {'o', 's', '^'};
cl     = {'r', 'b', 'g'};

for p = 1:numel(Products)
    nm = Products(p).name;
    f = figure; hold on
    hs = gobjects(3,1);
    for c = 1:3
        st = COMP.(nm).(comps{c});
        hs(c) = scatter(pairs.(['drifter_' comps{c}]), pairs.([nm '_' comps{c}]), ...
            36, cl{c}, mk{c}, 'filled', ...
            'DisplayName', sprintf('%s  r=%.2f  RMSE=%.3f  n=%d', ...
            compLabel(comps{c}), st.r, st.rmse, st.n));
    end
    plot([-L L], [-L L], 'k', 'DisplayName', '1:1 line');
    hold off
    axis equal; xlim([-L L]); ylim([-L L]); grid on
    legend('Location', 'best')
    xlabel('Drifter Velocities (m/s)'); ylabel('ADCP Velocities (m/s)');
    title(titles{p});
    exportgraphics(f, fullfile(cfg.out.comp_figures, [cfg.obsTag '_' files{p} '_drifter_comp.pdf']), 'Resolution', 300);
    exportgraphics(f, fullfile(cfg.out.comp_figures, [cfg.obsTag '_' files{p} '_drifter_comp.png']), 'Resolution', 300);
end

% All products, one panel per component
All = figure;
for c = 1:3
    ax = subplot(1,3,c); hold(ax, 'on')
    for p = 1:numel(Products)
        nm = Products(p).name;
        scatter(ax, pairs.(['drifter_' comps{c}]), pairs.([nm '_' comps{c}]), ...
            30, cl{p}, mk{p}, 'filled', 'DisplayName', nm);
    end
    plot(ax, [-L L], [-L L], 'k', 'DisplayName', '1:1 line');
    hold(ax, 'off'); axis(ax, 'equal'); xlim(ax, [-L L]); ylim(ax, [-L L]); grid(ax, 'on')
    title(ax, compLabel(comps{c}));
    xlabel(ax, 'Drifter (m/s)');
    if c == 1, ylabel(ax, 'ADCP product (m/s)'); legend(ax, 'Location', 'best'); end
end
sgtitle(strrep(cfg.obsTag, '_', '\_'));
exportgraphics(All, fullfile(cfg.out.comp_figures, [cfg.obsTag '_all_drifter_comp.pdf']), 'Resolution', 300);
exportgraphics(All, fullfile(cfg.out.comp_figures, [cfg.obsTag '_all_drifter_comp.png']), 'Resolution', 300);

end

%% ------------------------------------------------------------------------
function [ve, vn] = unrotate(u, v, R)
% Rotated (major, minor) -> (east, north); works for any array shape
sz = size(u);
uv = [u(:) v(:)]*R';
ve = reshape(uv(:,1), sz);
vn = reshape(uv(:,2), sz);
end

function st = compStats(x, y)
% Paired statistics; bias and RMSE are x - y (drifter - ADCP),
% matching pairStats.m used by the all-period scripts
ok = ~isnan(x) & ~isnan(y);
st.n = nnz(ok);
if st.n >= 3
    c = corrcoef(x(ok), y(ok));
    st.r = c(1,2);
else
    st.r = NaN;
end
d = x(ok) - y(ok);
if st.n > 0
    st.bias = mean(d);
    st.rmse = sqrt(mean(d.^2));
    st.mean_abs_drifter = mean(abs(x(ok)));
    st.mean_abs_adcp    = mean(abs(y(ok)));
else
    st.bias = NaN; st.rmse = NaN;
    st.mean_abs_drifter = NaN; st.mean_abs_adcp = NaN;
end
end

function s = compLabel(c)
switch c
    case 'V',  s = '|V|';
    case 'VE', s = 'u (east)';
    case 'VN', s = 'v (north)';
end
end

function n = getEOFmodes(L2)
if isfield(L2, 'settings') && isfield(L2.settings, 'eof_nmodes')
    n = L2.settings.eof_nmodes;
else
    n = 3;
end
end

function v = getcfg(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end
