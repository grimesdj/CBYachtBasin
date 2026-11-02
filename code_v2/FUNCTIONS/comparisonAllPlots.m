function ALL = comparisonAllPlots(cfg)
%COMPARISONALLPLOTS  Drifter vs ADCP scatter plots and statistics,
% all observation periods combined.
%
% Revised 2026-10-13:
%   - Reads the paired values from each period's COMP (written by the
%     revised comparison.m), so drifter and ADCP values are matched
%     crossing by crossing
%   - Period list comes from allObsPeriods.m
%   - Statistics (n, r, bias, RMSE) per product; the pooled
%     "ALL_PRODUCTS" row is gone (it counted each crossing three times)
%   - Per-period breakdown table for |V|
%   - Bias = drifter - ADCP throughout
%   - Figures saved to cfg.out.comp_figures with an ALL_ prefix
%
% Run with cfg from get_obs_config('ALL', 'NorthADCP', ...).

P      = allObsPeriods();
prods  = {'ADCP', 'EXTRAP', 'EOF'};
comps  = {'V', 'VE', 'VN'};
titles = {'ADCP OBSERVATIONS (top reliable bin)', ...
          'EXTRAPOLATED OBSERVATIONS', ...
          'EOF RECONSTRUCTION'};
files  = {'adcp', 'extrap', 'eof'};

%% Collect paired values from every period
for c = 1:3
    Dr.(comps{c}) = zeros(0,1);
    for p = 1:3
        Ad.(prods{p}).(comps{c}) = zeros(0,1);
    end
end
period = strings(0,1);

for n = 1:numel(P)
    C = loadPeriodComp(cfg.root_data, P(n));
    if isempty(C), continue, end
    pr = C.pairs;
    m  = numel(pr.j);
    fprintf('%s: %d crossings\n', P(n).obsTag, m);

    for c = 1:3
        Dr.(comps{c}) = [Dr.(comps{c}); pr.(['drifter_' comps{c}])(:)];
        for p = 1:3
            Ad.(prods{p}).(comps{c}) = [Ad.(prods{p}).(comps{c}); pr.([prods{p} '_' comps{c}])(:)];
        end
    end
    period = [period; repmat(string(P(n).obsTag), m, 1)]; %#ok<AGROW>
end

if isempty(period)
    error('comparisonAllPlots: no ROI crossings found in any period.');
end

%% Overall statistics, per product
fprintf('\nALL PERIODS: %d crossings\n', numel(period));
fprintf('%-7s %-3s %4s %7s %9s %8s\n', 'Product', 'Var', 'n', 'r', 'bias', 'RMSE');
for p = 1:3
    for c = 1:3
        st = pairStats(Dr.(comps{c}), Ad.(prods{p}).(comps{c}));
        ALL.stats.(prods{p}).(comps{c}) = st;
        fprintf('%-7s %-3s %4d %7.2f %9.3f %8.3f\n', prods{p}, comps{c}, st.n, st.r, st.bias, st.rmse);
    end
end
fprintf('(bias = drifter - ADCP, m/s)\n\n');

% Summary table (RMSE) in the same layout as before, minus the pooled row
rmseVals = nan(3,3);
for p = 1:3
    for c = 1:3
        rmseVals(p,c) = ALL.stats.(prods{p}).(comps{c}).rmse;
    end
end
ALL.RMSE_Table = array2table(rmseVals, 'VariableNames', comps, 'RowNames', prods);

%% Per-period breakdown for |V|
tags = unique(period, 'stable');
vals = nan(numel(tags), 9);
vnames = {};
for p = 1:3
    vnames = [vnames, {[prods{p} '_n'], [prods{p} '_bias'], [prods{p} '_RMSE']}]; %#ok<AGROW>
end
for t = 1:numel(tags)
    idx = period == tags(t);
    for p = 1:3
        st = pairStats(Dr.V(idx), Ad.(prods{p}).V(idx));
        vals(t, (p-1)*3 + (1:3)) = [st.n st.bias st.rmse];
    end
end
ALL.byPeriod_V = array2table(vals, 'VariableNames', vnames, 'RowNames', cellstr(tags));
disp('|V| BY OBSERVATION PERIOD (bias = drifter - ADCP, m/s)')
disp(ALL.byPeriod_V)

%% Store paired data
ALL.period  = period;
ALL.drifter = Dr;
for p = 1:3
    ALL.(prods{p}) = Ad.(prods{p});
end

%% Axis limit shared by all figures
allVals = [Dr.V; Dr.VE; Dr.VN];
for p = 1:3
    allVals = [allVals; Ad.(prods{p}).V; Ad.(prods{p}).VE; Ad.(prods{p}).VN]; %#ok<AGROW>
end
L = ceil(10*max(abs(allVals), [], 'omitnan'))/10;
if isempty(L) || isnan(L) || L == 0, L = 0.4; end

mk = {'o', 's', '^'};
cl = {'r', 'b', 'g'};
compName = {'|V|', 'u (east)', 'v (north)'};

%% One figure per product (markers = component)
for p = 1:3
    f = figure; hold on
    for c = 1:3
        st = ALL.stats.(prods{p}).(comps{c});
        scatter(Dr.(comps{c}), Ad.(prods{p}).(comps{c}), 36, cl{c}, mk{c}, 'filled', ...
            'DisplayName', sprintf('%s  r=%.2f  bias=%.3f  RMSE=%.3f  n=%d', ...
            compName{c}, st.r, st.bias, st.rmse, st.n));
    end
    plot([-L L], [-L L], 'k', 'DisplayName', '1:1 line');
    hold off
    axis equal; xlim([-L L]); ylim([-L L]); grid on
    legend('Location', 'best')
    xlabel('Drifter Velocities (m/s)'); ylabel('ADCP Velocities (m/s)');
    title([titles{p} ' - ALL OBSERVATION PERIODS']);
    saveFig(f, cfg, ['ALL_' files{p} '_all_obs']);
end

%% Products side by side, one panel per component
f = figure;
for c = 1:3
    ax = subplot(1,3,c); hold(ax, 'on')
    for p = 1:3
        scatter(ax, Dr.(comps{c}), Ad.(prods{p}).(comps{c}), 30, cl{p}, mk{p}, 'filled', ...
            'DisplayName', prods{p});
    end
    plot(ax, [-L L], [-L L], 'k', 'DisplayName', '1:1 line');
    hold(ax, 'off'); axis(ax, 'equal'); xlim(ax, [-L L]); ylim(ax, [-L L]); grid(ax, 'on')
    title(ax, compName{c}); xlabel(ax, 'Drifter (m/s)');
    if c == 1, ylabel(ax, 'ADCP product (m/s)'); legend(ax, 'Location', 'best'); end
end
sgtitle('ALL OBSERVATION PERIODS - products compared');
saveFig(f, cfg, 'ALL_products_by_component');

%% |V| coloured by observation period, one panel per product
f = figure;
cmap = lines(numel(tags));
for p = 1:3
    ax = subplot(1,3,p); hold(ax, 'on')
    for t = 1:numel(tags)
        idx = period == tags(t);
        scatter(ax, Dr.V(idx), Ad.(prods{p}).V(idx), 30, cmap(t,:), 'filled', ...
            'DisplayName', strrep(char(tags(t)), '_', '\_'));
    end
    plot(ax, [0 L], [0 L], 'k', 'DisplayName', '1:1 line');
    hold(ax, 'off'); axis(ax, 'equal'); xlim(ax, [0 L]); ylim(ax, [0 L]); grid(ax, 'on')
    title(ax, [prods{p} ' |V|']); xlabel(ax, 'Drifter |V| (m/s)');
    if p == 1, ylabel(ax, 'ADCP |V| (m/s)'); end
    if p == 3, legend(ax, 'Location', 'bestoutside'); end
end
sgtitle('|V| by observation period');
saveFig(f, cfg, 'ALL_speed_by_period');

end

%% ------------------------------------------------------------------------
function saveFig(f, cfg, base)
exportgraphics(f, fullfile(cfg.out.comp_figures, [base '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(cfg.out.comp_figures, [base '.pdf']), 'Resolution', 300);
end
