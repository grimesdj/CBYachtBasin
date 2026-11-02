function [adcp_level_2] = adcp_post_processing(cfg)
%ADCP_POST_PROCESSING  ADCP Level 2: principal axes, surface/bed
% extrapolation, z/h mapping and complex EOF analysis.
%
% Revised 2026-10-13. Changes from the previous version are tagged [FIX n]:
%   [FIX 1] PCA covariances use valid (non-NaN) ensembles only
%   [FIX 2] Major-axis sign set by the tide: +major = flood (rising eta)
%   [FIX 3] Smoothing for the figures no longer spreads NaNs
%   [FIX 4] Water depth h = pressure (dbar -> m) + pressure-sensor height
%   [FIX 5] Extrapolation rewritten per profile: separate u and v slopes,
%           no-slip at z = 0, nothing extrapolated above the surface
%   [FIX 6] z/h mapping per profile with interp1 (replaces griddata)
%   [FIX 7] Only the time mean is removed before the EOF (was linear detrend)
%   [FIX 8] Complex EOF phase fixed so mode-shape plots are reproducible
%   [FIX 9] Built-ins norm/eig no longer shadowed
%
% Output fields read downstream keep their names and shapes:
%   rot_mat, u2_extrap (time x bins+1), v2_extrap, dbins_extrap,
%   U_mapped (z/h x time), V_mapped, zonh_grid, EC_cmplx, EOFs_cmplx,
%   trend_cmplx, stdMode_cmplx, scaledEOFs_cmplx, frac_of_var,
%   U_predicted, V_predicted (anomalies, time x z/h)
% New fields: h, flood_corr, pca_var, time_valid, U_recon, V_recon,
%   extrap, settings
%
% Optional cfg fields (set in get_obs_config):
%   cfg.p_sensor_hab   pressure-sensor height above bed, m   (default 0 + warning)
%   cfg.rho            seawater density, kg/m^3              (default 1025)
%   cfg.extrap_method  'linear' or 'constant'                (default 'linear')
%   cfg.extrap_nfit    top bins used for the surface slope   (default 3)
%   cfg.eof_nmodes     modes kept in the reconstruction      (default 3)

fname = fullfile(cfg.out.adcp_data, [cfg.obsTag '_adcp_level_1.mat']);
S = load(fname, 'adcp_level_1');
A = S.adcp_level_1;

%% Settings
p_hab        = getcfg(cfg, 'p_sensor_hab', NaN);
rho          = getcfg(cfg, 'rho', 1025);
extrapMethod = getcfg(cfg, 'extrap_method', 'linear');
nfit         = getcfg(cfg, 'extrap_nfit', 3);
nModesPred   = getcfg(cfg, 'eof_nmodes', 3);
zonh_grid    = 0:0.05:1;

if isnan(p_hab)
    warning('adcp_post_processing:noSensorHeight', ...
        ['cfg.p_sensor_hab not set for %s: assuming the pressure sensor ' ...
         'is at the bed (h = pressure).'], cfg.obsTag);
    p_hab = 0;
end

%% Extract Level 1
V1        = A.V_E;              % bins x time, NaN where QC failed
V2        = A.V_N;
dbins     = A.bins(:)';         % 1 x nb, m above bed
East_avg  = A.East_avg(:);
North_avg = A.North_avg(:);
pressure  = A.pressure(:);
eta       = A.eta(:);
east      = A.east;             % time x bins, NaN where QC failed
north     = A.north;
f1        = A.f1;
f2        = A.f2;
qcFlag    = A.qcFlag;           % time x bins
time      = A.time(:);
uv_caxis  = A.uv_caxis;

nt = numel(time);
nb = numel(dbins);
assert(isequal(size(east), [nt nb]), 'Expected east/north as time x bins.');

EOF = struct();
EOF.qcFlag = qcFlag;

%% [FIX 4] Total water depth
dbar2m = 1e4/(rho*9.81);              % ~0.995 m per dbar
h = pressure*dbar2m + p_hab;          % nt x 1, m (surface height above bed)
EOF.h = h;

%% Principal axes of the depth-averaged current  [FIX 1, 2, 9]
ok  = ~isnan(East_avg) & ~isnan(North_avg);
Nok = nnz(ok);
ue  = East_avg(ok)  - mean(East_avg(ok));
vn  = North_avg(ok) - mean(North_avg(ok));

u2avg = (ue'*ue)/(Nok-1);
v2avg = (vn'*vn)/(Nok-1);
uvavg = (ue'*vn)/(Nok-1);
disc    = sqrt((u2avg - v2avg)^2 + 4*uvavg^2);
sig2pos = 0.5*((u2avg + v2avg) + disc);
sig2neg = 0.5*((u2avg + v2avg) - disc);

if abs(uvavg) < eps
    % Already aligned with east/north
    if u2avg >= v2avg
        Vrot = eye(2);
    else
        Vrot = [0 -1; 1 0];
    end
    lam = [sig2pos sig2neg];
else
    Vrot = [uvavg,          sig2neg - v2avg;
            sig2pos - u2avg, uvavg];
    Vrot = Vrot*diag(1./sqrt(diag(Vrot'*Vrot)));   % unit eigenvectors
    lam  = [sig2pos sig2neg];
    [lam, srt] = sort(lam, 'descend');
    Vrot = Vrot(:, srt);
end

% Right-handed coordinate system
if det(Vrot) < 0
    Vrot(:,2) = -Vrot(:,2);
end

% Sign convention: positive major axis = flood (current correlated with
% rising water level). A 180-degree flip of both axes keeps det > 0.
uv1  = [ue vn]*Vrot;
deta = gradient(eta(ok), time(ok));
cc   = corrcoef(uv1(:,1), deta, 'Rows', 'complete');
floodCorr = cc(1,2);
if floodCorr < 0
    Vrot = -Vrot;
    uv1  = -uv1;
    floodCorr = -floodCorr;
end

EOF.rot_mat    = Vrot;
EOF.flood_corr = floodCorr;   % near 0 = sign choice is weak; check by eye
EOF.pca_var    = lam;

fig1 = figure;
p1 = plot(East_avg, North_avg, '.r'); hold on
p2 = plot(uv1(:,1), uv1(:,2), '.c');
p3 = plot([0 sqrt(lam(1))*Vrot(1,1)], [0 sqrt(lam(1))*Vrot(2,1)], '--b', ...
          [0 sqrt(lam(2))*Vrot(1,2)], [0 sqrt(lam(2))*Vrot(2,2)], '--b');
hold off
axis equal; grid on
legend([p1(1) p2(1) p3(1)], 'input data', 'rotated data', 'eigen vectors (major = flood)')
xlabel('rotated-major-axis'); ylabel('rotated-minor-axis')
title(sprintf('%s   corr(u_{major}, d\\eta/dt) = %.2f', ...
    strrep(cfg.obsTag, '_', '\_'), floodCorr))
exportgraphics(fig1, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_adcp_PCA.png']), 'Resolution', 300);

%% Rotate the depth-varying currents
uv2 = [east(:), north(:)]*Vrot;
u2  = reshape(uv2(:,1), nt, nb);      % time x bins, major axis
v2  = reshape(uv2(:,2), nt, nb);      % time x bins, minor axis

%% Smoothed fields for the figures only  [FIX 3]
q   = double(qcFlag ~= 0 & ~isnan(u2) & ~isnan(v2));
u2z = u2;  u2z(q == 0) = 0;
v2z = v2;  v2z(q == 0) = 0;
on  = conv2(f1, f2, q', 'same');      % bins x time
V_1 = conv2(f1, f2, u2z', 'same')./on;
V_2 = conv2(f1, f2, v2z', 'same')./on;
qplot = q';
qplot(on <= 0.95) = 0;

ylims = [0 6];

fig2 = figure;
ax1 = subplot(2,1,1);
imagesc(time, dbins, V_1, 'AlphaData', qplot); caxis(uv_caxis)
colormap(ax1, cmocean('balance')); colorbar
hold on; plot(time, h, 'k'); hold off
title('Principal axis (+ = flood)')
set(ax1, 'ydir', 'normal', 'ticklabelinterpreter', 'latex', 'ylim', ylims, 'XTickLabel', [])

ax2 = subplot(2,1,2);
imagesc(time, dbins, V_2, 'AlphaData', qplot); caxis([-0.3 0.3])
colormap(ax2, cmocean('balance')); cb = colorbar; set(cb, 'ytick', -0.3:0.1:0.3)
hold on; plot(time, h, 'k'); hold off
title('Minor axis')
set(ax2, 'ydir', 'normal', 'ticklabelinterpreter', 'latex', 'ylim', ylims)
datetick(ax2, 'x', 'mm/dd HH:MM', 'keeplimits')
ylabel('mab', 'interpreter', 'latex')
linkaxes([ax1 ax2], 'x')
exportgraphics(fig2, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_depth_varying_PCA.png']), 'Resolution', 300);

fig3 = figure;
ax1 = subplot(2,1,1);
him = imagesc(time, dbins, V1); set(him, 'AlphaData', ~isnan(V1)); caxis(uv_caxis)
colormap(ax1, cmocean('balance')); colorbar
hold on; plot(time, h, 'k'); hold off
title('East')
set(ax1, 'ydir', 'normal', 'ticklabelinterpreter', 'latex', 'ylim', ylims, 'XTickLabel', [])

ax2 = subplot(2,1,2);
him = imagesc(time, dbins, V2); set(him, 'AlphaData', ~isnan(V2)); caxis(uv_caxis)
colormap(ax2, cmocean('balance')); colorbar
hold on; plot(time, h, 'k'); hold off
title('North')
set(ax2, 'ydir', 'normal', 'ticklabelinterpreter', 'latex', 'ylim', ylims)
datetick(ax2, 'x', 'mm/dd', 'keeplimits')
ylabel('mab', 'interpreter', 'latex')
linkaxes([ax1 ax2], 'x')
exportgraphics(fig3, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_depth_varying_raw.png']), 'Resolution', 300);

%% Extrapolate each profile to bed and surface, and map to z/h  [FIX 5, 6]
% Each profile is built on its own vertical coordinate:
%   z = 0 (bed, no-slip u = v = 0) -> valid bins -> z = h (surface)
% Interior gaps are filled linearly in z. Bins above the surface stay NaN.
dbins_extrap = [0 dbins];
u2_extrap = nan(nt, nb+1);
v2_extrap = nan(nt, nb+1);
U_mapped  = nan(numel(zonh_grid), nt);
V_mapped  = nan(numel(zonh_grid), nt);

ex.z_top   = nan(nt,1);     % height of highest valid bin, m
ex.u_top   = nan(nt,1);  ex.v_top  = nan(nt,1);
ex.u_surf  = nan(nt,1);  ex.v_surf = nan(nt,1);
ex.n_valid = zeros(nt,1);
minBins = 2;

for i = 1:nt
    if isnan(h(i)), continue; end

    good = ~isnan(u2(i,:)) & ~isnan(v2(i,:)) & dbins > 0 & dbins < h(i);
    ex.n_valid(i) = nnz(good);
    if nnz(good) < minBins, continue; end

    zg = dbins(good);
    ug = u2(i,good);
    vg = v2(i,good);

    switch extrapMethod
        case 'linear'
            % Least-squares shear over the top nfit valid bins, applied
            % separately to u and v, anchored at the top valid bin
            k  = max(1, numel(zg) - nfit + 1):numel(zg);
            pu = polyfit(zg(k), ug(k), 1);
            pv = polyfit(zg(k), vg(k), 1);
            us = ug(end) + pu(1)*(h(i) - zg(end));
            vs = vg(end) + pv(1)*(h(i) - zg(end));
        case 'constant'
            us = ug(end);
            vs = vg(end);
        otherwise
            error('Unknown cfg.extrap_method: %s', extrapMethod)
    end

    zp = [0 zg h(i)];
    up = [0 ug us];
    vp = [0 vg vs];

    inWater = dbins_extrap <= h(i);
    u2_extrap(i,inWater) = interp1(zp, up, dbins_extrap(inWater));
    v2_extrap(i,inWater) = interp1(zp, vp, dbins_extrap(inWater));

    U_mapped(:,i) = interp1(zp, up, zonh_grid*h(i));
    V_mapped(:,i) = interp1(zp, vp, zonh_grid*h(i));

    ex.z_top(i)  = zg(end);
    ex.u_top(i)  = ug(end);  ex.v_top(i)  = vg(end);
    ex.u_surf(i) = us;       ex.v_surf(i) = vs;
end

EOF.u2_extrap    = u2_extrap;
EOF.v2_extrap    = v2_extrap;
EOF.dbins_extrap = dbins_extrap;
EOF.U_mapped     = U_mapped;
EOF.V_mapped     = V_mapped;
EOF.zonh_grid    = zonh_grid;
EOF.extrap       = ex;

% QC: how far, and by how much, each profile was extrapolated
figX = figure;
axA = subplot(2,1,1);
plot(time, ex.u_top, 'k', time, ex.u_surf, 'r')
legend('top valid bin', 'extrapolated surface', 'Location', 'best')
ylabel('u_{major} (m/s)'); grid on
title(sprintf('Surface extrapolation (%s)', extrapMethod))
set(axA, 'XTickLabel', [])
axB = subplot(2,1,2);
plot(time, h - ex.z_top, 'k')
ylabel('top bin to surface (m)'); grid on
datetick(axB, 'x', 'mm/dd HH:MM', 'keeplimits')
linkaxes([axA axB], 'x')
exportgraphics(figX, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_extrap_QC.png']), 'Resolution', 300);

%% Complex EOF  [FIX 7, 8]
d0     = (U_mapped + 1i*V_mapped).';     % time x z/h
tvalid = all(~isnan(d0), 2);             % profiles with a complete z/h column
Nv     = nnz(tvalid);
if Nv < 2
    error('adcp_post_processing: fewer than 2 complete profiles for the EOF (%s).', cfg.obsTag);
end

trend = mean(d0(tvalid,:), 1);           % time-mean profile
d_complex = nan(size(d0));
d_complex(tvalid,:) = d0(tvalid,:) - trend;

[E, S_pc, EOFs] = svd(d_complex(tvalid,:), 'econ');
ECv = E*S_pc;

% Fix each mode's arbitrary complex phase: the mode is real and positive
% where its magnitude peaks. Reconstruction EC*EOFs' is unchanged.
for k = 1:size(EOFs,2)
    [~, kmax] = max(abs(EOFs(:,k)));
    ph = exp(-1i*angle(EOFs(kmax,k)));
    EOFs(:,k) = EOFs(:,k)*ph;
    ECv(:,k)  = ECv(:,k)*ph;
end

EC = nan(nt, size(EOFs,2));              % NaN rows for incomplete profiles
EC(tvalid,:) = ECv;

stdMode    = diag(S_pc)/sqrt(Nv);
scaledEOFs = EOFs.*stdMode.';
fraction_of_variance = stdMode.^2./sum(stdMode.^2);

EOF.EC_cmplx         = EC;
EOF.stdMode_cmplx    = stdMode;
EOF.scaledEOFs_cmplx = scaledEOFs;
EOF.EOFs_cmplx       = EOFs;
EOF.trend_cmplx      = trend;
EOF.frac_of_var      = fraction_of_variance;
EOF.time_valid       = tvalid;

fig4 = figure;
nShow = min(10, numel(fraction_of_variance));
plot(1:nShow, fraction_of_variance(1:nShow), 'o-'); grid on
xlabel('Mode'); ylabel('Fraction of variance')
exportgraphics(fig4, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_frac_of_var.png']), 'Resolution', 300);

%% Mode shapes
nPlot = min(3, size(EOFs,2));

fig5 = figure;
plot(real(trend), zonh_grid, 'k', real(scaledEOFs(:,1:nPlot)), zonh_grid, 'linewidth', 2);
xline(0); grid on
xlabel('Amplitude (m/s)'); ylabel('z/h')
legend(['Time Averaged Velocity', arrayfun(@(k) sprintf('Mode %d', k), 1:nPlot, 'UniformOutput', false)], ...
    'location', 'best');
exportgraphics(fig5, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_3_modes_w_trend.png']), 'Resolution', 300);

fig6 = figure;
for k = 1:nPlot
    subplot(1, nPlot, k)
    plot(real(trend), zonh_grid, 'k', imag(trend), zonh_grid, '--k', ...
         real(scaledEOFs(:,k)), zonh_grid, 'r', imag(scaledEOFs(:,k)), zonh_grid, '--r')
    title(sprintf('Mode %d (%.0f%%)', k, 100*fraction_of_variance(k)))
    xlabel('m/s'); grid on
    if k == 1
        ylabel('z/h')
        legend('mean, major', 'mean, minor', 'mode, major', 'mode, minor', 'Location', 'best')
    end
end
exportgraphics(fig6, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_sep_modes.png']), 'Resolution', 300);

%% Reconstruction from the leading modes
nm   = min(nModesPred, size(EOFs,2));
pred = EC(:,1:nm)*EOFs(:,1:nm)';          % time x z/h, anomaly

U_pred_cmplx_EOF = real(pred);
V_pred_cmplx_EOF = imag(pred);

EOF.U_predicted = U_pred_cmplx_EOF;                       % anomaly (as before)
EOF.V_predicted = V_pred_cmplx_EOF;
EOF.U_recon     = U_pred_cmplx_EOF + real(trend);         % mean + modes
EOF.V_recon     = V_pred_cmplx_EOF + imag(trend);

fig8 = figure;
plotField(subplot(3,1,1), time, zonh_grid, EOF.U_recon', [-0.5 0.5], ...
    sprintf('U reconstructed (mean + %d modes)', nm));
plotField(subplot(3,1,2), time, zonh_grid, U_mapped, [-0.5 0.5], 'U observed (mapped)');
plotField(subplot(3,1,3), time, zonh_grid, real(d_complex)', [-0.5 0.5], 'U anomaly (mean removed)');
exportgraphics(fig8, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_U_pred_from_EOF.png']), 'Resolution', 300);

fig9 = figure;
plotField(subplot(2,1,1), time, zonh_grid, EOF.V_recon', [-0.2 0.2], ...
    sprintf('V reconstructed (mean + %d modes)', nm));
plotField(subplot(2,1,2), time, zonh_grid, V_mapped, [-0.2 0.2], 'V observed (mapped)');
exportgraphics(fig9, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_V_pred_from_EOF.png']), 'Resolution', 300);

%% Residuals (anomaly minus reconstructed anomaly)
residuals_U = real(d_complex) - U_pred_cmplx_EOF;
residuals_V = imag(d_complex) - V_pred_cmplx_EOF;

fig10 = figure;
plotField(subplot(2,1,1), time, zonh_grid, residuals_U', [-0.2 0.2], 'Residuals from EOF; U');
plotField(subplot(2,1,2), time, zonh_grid, residuals_V', [-0.1 0.1], 'Residuals from EOF; V');
exportgraphics(fig10, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_EOF_residuals.png']), 'Resolution', 300);

%% Record settings and return
EOF.settings = struct('p_sensor_hab', p_hab, 'rho', rho, 'dbar2m', dbar2m, ...
    'extrap_method', extrapMethod, 'extrap_nfit', nfit, 'eof_nmodes', nm, ...
    'min_bins', minBins);

adcp_level_2 = EOF;

end

%% ------------------------------------------------------------------------
function v = getcfg(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end

function plotField(ax, time, z, C, clim, ttl)
him = imagesc(ax, time, z, C);
set(him, 'AlphaData', ~isnan(C));
set(ax, 'YDir', 'normal');
caxis(ax, clim);
colormap(ax, cmocean('balance'));
colorbar(ax);
title(ax, ttl);
ylabel(ax, 'z/h');
datetick(ax, 'x', 'mm/dd HH:MM', 'keeplimits');
end
