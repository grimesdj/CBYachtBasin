function ROI(cfg, adcpLoc)
%ROI  Level 4: vertical profiles at each ROI crossing (ADCP, EOF, drifter).
%
% Revised 2026-10-13:
%   - EOF profile is Level 2 U_recon/V_recon (mean + leading modes)
%   - Averages over the transit with 'omitnan' (Level 2 now leaves NaN
%     rows for profiles that could not be extrapolated)
%   - Raw ADCP profile averaged over the transit (was first sample only)
%   - Water depth from Level 2 h; y-axis is height above bed
%   - Returns quietly when comparison was skipped (e.g. SouthADCP)

base_name = regexprep(cfg.name, '_sheet$', '');

L2 = load(fullfile(cfg.out.adcp_data, [base_name '_' adcpLoc '_adcp_level_2.mat']));
adcp_level_2 = L2.adcp_level_2;

fname_ROI = fullfile(cfg.out.comp_data, [cfg.obsTag '_ROI_ADCP.mat']);
S = load(fname_ROI, 'ROI_ADCP');
ROI_ADCP = S.ROI_ADCP;

if isempty(ROI_ADCP)
    fprintf('ROI: no comparison results for %s; nothing to plot.\n', cfg.obsTag);
    return
end
if ~isfield(adcp_level_2, 'U_recon')
    error('ROI: Level 2 for %s has no U_recon. Re-run Level 2 first.', cfg.obsTag);
end

z_target  = getcfg(cfg, 'drifter_depth', 0.5);
zonh_grid = adcp_level_2.zonh_grid;

%% EOF reconstruction back to east/north (time x z/h)
U = adcp_level_2.U_recon;
V = adcp_level_2.V_recon;
uv = [U(:) V(:)]*adcp_level_2.rot_mat';
ve_eof = reshape(uv(:,1), size(U));
vn_eof = reshape(uv(:,2), size(U));

%% One figure per crossing
for j = 1:size(ROI_ADCP,1)
    for i = 1:size(ROI_ADCP,2)

        R = ROI_ADCP(j,i);
        if isempty(R.inds)
            continue
        end

        inds = R.inds;  inde = R.inde;
        hT = mean(R.adcp_h, 'omitnan');

        u_eof = mean(ve_eof(inds:inde,:), 1, 'omitnan');
        v_eof = mean(vn_eof(inds:inde,:), 1, 'omitnan');

        adcp_ve = mean(R.adcp_ve, 2, 'omitnan');
        adcp_vn = mean(R.adcp_vn, 2, 'omitnan');

        fig = figure;

        subplot(1,2,1)
        plot(adcp_vn, R.adcp_dbins, 'k', 'LineWidth', 2); hold on
        plot(v_eof, zonh_grid*hT, 'r', 'LineWidth', 2)
        plot(R.avg_vn, hT - z_target, 'ro', 'LineWidth', 2)
        yline(hT, '--b');
        hold off
        xlabel('m/s'); ylabel('Height above bed (m)');
        title('v (north)')
        xlim([-0.7 0.7]); ylim([0 6]); grid on

        subplot(1,2,2)
        plot(adcp_ve, R.adcp_dbins, 'k', 'LineWidth', 2); hold on
        plot(u_eof, zonh_grid*hT, 'r', 'LineWidth', 2)
        plot(R.avg_ve, hT - z_target, 'ro', 'LineWidth', 2)
        yline(hT, '--b');
        hold off
        xlabel('m/s'); ylabel('Height above bed (m)');
        title('u (east)')
        legend('ADCP (transit mean)', 'EOF (mean + modes)', 'Drifter', 'Surface', ...
            'Location', 'northeast');
        xlim([-0.4 0.4]); ylim([0 6]); grid on

        sgtitle(sprintf('Dep %d, drifter %d: %s', j, i, ...
            datestr(R.adcp_Time_datnum(1), 'mm/dd HH:MM')));

        fname   = sprintf('_ROI_Profile_raw_j%02d_i%02d.png', j, i);
        fname_1 = sprintf('_ROI_Profile_raw_j%02d_i%02d.fig', j, i);
        exportgraphics(fig, fullfile(cfg.out.comp_figures, [cfg.obsTag fname]), 'Resolution', 300);
        savefig(fig, fullfile(cfg.out.comp_figures, [cfg.obsTag fname_1]));
%         close(fig)
    end
end

end

%% ------------------------------------------------------------------------
function v = getcfg(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end
