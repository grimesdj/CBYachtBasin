function ALL = comparisonAll(cfg)
%COMPARISONALL  Drifter vs ADCP agreement as a function of the drifter
% differencing interval, all observation periods combined.
%
% Revised 2026-10-13:
%   - Period list comes from allObsPeriods.m. (Previously the Results-file
%     list had 5 entries for 6 periods - SummerFlood was missing, so the
%     loop would index past the end of resultsFiles.)
%   - ADCP values, ROI centre and ROI radius come from each period's COMP
%     (written by the revised comparison.m), so they match Level 3 exactly
%   - Two drifter speed definitions at every differencing interval:
%       SpeedOfMean = |(mean v_e, mean v_n)|  - same definition as Level 3
%       MeanSpeed   = mean(|V|)               - the old definition
%     plus the u and v components. The ADCP side is always the magnitude
%     of the mean velocity.
%   - Statistics per product only; bias = drifter - ADCP
%   - Stops with a clear message if a period still has old Level 3 output
%
% Note: drifter ROI points here come from the veloAnalysis Results, which
% drop the first/last 128 s of each track, so a crossing near the start or
% end of a track can have fewer points than in Level 3, or none.
%
% Run with cfg from get_obs_config('ALL', 'NorthADCP', ...).

P       = allObsPeriods();
res_sec = [1 2 4 8 16 32 64 128 256];
nRes    = numel(res_sec);
prods   = {'ADCP', 'EXTRAP', 'EOF'};
proj    = projcrs(32119);

% Drifter values: crossings x resolutions
Dsom = zeros(0, nRes);  Dmsp = Dsom;  Dve = Dsom;  Dvn = Dsom;
% ADCP values: crossings x 1
for p = 1:3
    AV.(prods{p}) = struct('V', zeros(0,1), 'VE', zeros(0,1), 'VN', zeros(0,1));
end
period = strings(0,1);
nNoResults = 0;

%% Collect drifter values at each interval, paired with Level 3 ADCP values
for n = 1:numel(P)

    C = loadPeriodComp(cfg.root_data, P(n));
    if isempty(C), continue, end
    pr = C.pairs;

    fR = fullfile(cfg.root_data, P(n).folder, 'PROCESSED_DATA', 'DRIFTERS', ...
        [P(n).name '_TemporalRes_RMSE.mat']);
    if ~isfile(fR)
        warning('comparisonAll:noResults', 'Missing %s; skipping %s.', fR, P(n).obsTag);
        continue
    end
    S = load(fR, 'Results');
    Results = S.Results;

    [E0, N0] = projfwd(proj, C.settings.adcp_latlon(1), C.settings.adcp_latlon(2));
    Rroi = C.settings.roi_radius;

    fprintf('%s: %d crossings\n', P(n).obsTag, numel(pr.j));

    for k = 1:numel(pr.j)
        j = pr.j(k);
        i = pr.i(k);
        row = nan(4, nRes);

        haveR = j <= size(Results,1) && i <= size(Results,2) && ~isempty(Results(j,i).RMSE);
        if haveR
            R = Results(j,i);
            for r = 1:nRes
                [E, N, ve, vn, sp] = seriesAt(R, res_sec(r));
                if isempty(E), continue, end
                in = hypot(E - E0, N - N0) <= Rroi;
                if ~any(in), continue, end
                mve = mean(ve(in), 'omitnan');
                mvn = mean(vn(in), 'omitnan');
                row(:,r) = [hypot(mve, mvn); mean(sp(in), 'omitnan'); mve; mvn];
            end
        else
            nNoResults = nNoResults + 1;   % track too short for veloAnalysis
        end

        Dsom(end+1,:) = row(1,:); %#ok<AGROW>
        Dmsp(end+1,:) = row(2,:); %#ok<AGROW>
        Dve(end+1,:)  = row(3,:); %#ok<AGROW>
        Dvn(end+1,:)  = row(4,:); %#ok<AGROW>
        for p = 1:3
            AV.(prods{p}).V(end+1,1)  = pr.([prods{p} '_V'])(k);
            AV.(prods{p}).VE(end+1,1) = pr.([prods{p} '_VE'])(k);
            AV.(prods{p}).VN(end+1,1) = pr.([prods{p} '_VN'])(k);
        end
        period(end+1,1) = string(P(n).obsTag); %#ok<AGROW>
    end
end

if isempty(period)
    error('comparisonAll: no ROI crossings found in any period.');
end
if nNoResults > 0
    fprintf(['%d crossing(s) have no veloAnalysis Results (track shorter than ' ...
             '2 x 128 s) and are NaN at every interval.\n'], nNoResults);
end

%% Statistics at each interval, per product
metrics  = {'SpeedOfMean', 'MeanSpeed', 'u', 'v'};
Dm       = {Dsom, Dmsp, Dve, Dvn};
adcpComp = {'V', 'V', 'VE', 'VN'};

for m = 1:numel(metrics)
    RMSE = nan(nRes,3);  Bias = RMSE;  Rr = RMSE;  Nn = RMSE;
    for p = 1:3
        a = AV.(prods{p}).(adcpComp{m});
        for r = 1:nRes
            st = pairStats(Dm{m}(:,r), a);
            RMSE(r,p) = st.rmse;
            Bias(r,p) = st.bias;
            Rr(r,p)   = st.r;
            Nn(r,p)   = st.n;
        end
    end
    ALL.(metrics{m}).RMSE = resTable(res_sec, RMSE, prods);
    ALL.(metrics{m}).Bias = resTable(res_sec, Bias, prods);
    ALL.(metrics{m}).r    = resTable(res_sec, Rr,   prods);
    ALL.(metrics{m}).N    = resTable(res_sec, Nn,   prods);
end

disp(' '); disp('RMSE BY DIFFERENCING INTERVAL - drifter |mean velocity| (Level 3 definition)')
disp(ALL.SpeedOfMean.RMSE)
disp('BIAS (drifter - ADCP) - drifter |mean velocity|')
disp(ALL.SpeedOfMean.Bias)
disp('RMSE BY DIFFERENCING INTERVAL - drifter mean |V| (old definition)')
disp(ALL.MeanSpeed.RMSE)
disp('BIAS (drifter - ADCP) - drifter mean |V|')
disp(ALL.MeanSpeed.Bias)

%% Store paired data
ALL.res_sec  = res_sec;
ALL.period   = period;
ALL.drifter  = struct('SpeedOfMean', Dsom, 'MeanSpeed', Dmsp, 'u', Dve, 'v', Dvn);
ALL.adcp     = AV;

%% Figures
cmap = lines(3);

f = figure;
ax1 = subplot(1,2,1); hold(ax1, 'on')
ax2 = subplot(1,2,2); hold(ax2, 'on')
for p = 1:3
    semilogx(ax1, res_sec, ALL.SpeedOfMean.RMSE.(prods{p}), '-o', 'Color', cmap(p,:), ...
        'LineWidth', 1.5, 'DisplayName', [prods{p} ', |mean v|']);
    semilogx(ax1, res_sec, ALL.MeanSpeed.RMSE.(prods{p}), '--', 'Color', cmap(p,:), ...
        'LineWidth', 1.2, 'DisplayName', [prods{p} ', mean |v|']);
    semilogx(ax2, res_sec, ALL.SpeedOfMean.Bias.(prods{p}), '-o', 'Color', cmap(p,:), ...
        'LineWidth', 1.5, 'DisplayName', [prods{p} ', |mean v|']);
    semilogx(ax2, res_sec, ALL.MeanSpeed.Bias.(prods{p}), '--', 'Color', cmap(p,:), ...
        'LineWidth', 1.2, 'DisplayName', [prods{p} ', mean |v|']);
end
yline(ax2, 0, 'k--', 'HandleVisibility', 'off');
set([ax1 ax2], 'XScale', 'log');
grid(ax1, 'on'); grid(ax2, 'on')
xlabel(ax1, 'Differencing Interval (s)'); ylabel(ax1, 'RMSE (m/s)');
xlabel(ax2, 'Differencing Interval (s)'); ylabel(ax2, 'Bias: Drifter - ADCP (m/s)');
title(ax1, 'Speed RMSE'); title(ax2, 'Speed bias');
legend(ax1, 'Location', 'best');
sgtitle('Drifter vs ADCP speed vs differencing interval (all periods)');
saveFig(f, cfg, 'ALL_resolution_speed');

f = figure;
comp = {'u', 'v'};
for c = 1:2
    ax = subplot(1,2,c); hold(ax, 'on')
    for p = 1:3
        semilogx(ax, res_sec, ALL.(comp{c}).RMSE.(prods{p}), '-o', 'Color', cmap(p,:), ...
            'LineWidth', 1.5, 'DisplayName', prods{p});
    end
    set(ax, 'XScale', 'log'); grid(ax, 'on')
    xlabel(ax, 'Differencing Interval (s)'); ylabel(ax, 'RMSE (m/s)');
    title(ax, [comp{c} ' RMSE']);
    if c == 1, legend(ax, 'Location', 'best'); end
end
sgtitle('Drifter vs ADCP components vs differencing interval (all periods)');
saveFig(f, cfg, 'ALL_resolution_components');

end

%% ------------------------------------------------------------------------
function [E, N, ve, vn, sp] = seriesAt(R, res)
% Positions and velocities of one trajectory at one differencing interval
if res == 1
    E = R.Easting;  N = R.Northing;
    ve = R.V_e_1Hz; vn = R.V_n_1Hz; sp = R.V_1Hz;
else
    k = find([R.TemporalRes.Seconds] == res, 1);
    if isempty(k)
        E = []; N = []; ve = []; vn = []; sp = [];
        return
    end
    TR = R.TemporalRes(k);
    E = TR.Easting;  N = TR.Northing;
    ve = TR.V_e;     vn = TR.V_n;     sp = TR.V;
end
end

function T = resTable(res_sec, M, prods)
T = array2table([res_sec(:) M], 'VariableNames', [{'Seconds'}, prods]);
end

function saveFig(f, cfg, base)
exportgraphics(f, fullfile(cfg.out.comp_figures, [base '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(cfg.out.comp_figures, [base '.pdf']), 'Resolution', 300);
end
