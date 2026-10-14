function ALL = comparisonAllPlots(cfg)

%% Observation periods / tags

obsTags = { ...
    'MarshMadness_NorthADCP'
    'MarshMayhem_NorthADCP'
    'FallFrolic_100825_NorthADCP'
    'FallFrolic_101025_NorthADCP'
%     'NovDep_NorthADCP'
    'JuneJamboree_NorthADCP'};

names = { ...
    'MarshMadness'
    'MarshMayhem'
    'FallFrolic'
    'FallFrolic'
%     'NovDep'
    'JuneJamboree'};


%% Initialize combined arrays

drifter_V  = [];
drifter_ve = [];
drifter_vn = [];

ADCP_V  = [];
ADCP_VE = [];
ADCP_VN = [];

EXTRAP_V  = [];
EXTRAP_VE = [];
EXTRAP_VN = [];

EOF_V  = [];
EOF_VE = [];
EOF_VN = [];

ObsName = strings(0,1);


%% Loop through observation periods

for n = 1:numel(obsTags)

    obsTag = obsTags{n};

    fprintf('Loading %s\n',obsTag);
    
    name = names{n};


    %% Load ROI structures

    A = load(fullfile( ...
        cfg.root_data, ...
        [name '\PROCESSED_DATA\COMP\' obsTag '_ROI_ADCP.mat']));

    E = load(fullfile( ...
        cfg.root_data, ...
        [name '\PROCESSED_DATA\COMP\' obsTag '_ROI_EXTRAP.mat']));

    F = load(fullfile( ...
        cfg.root_data, ...
        [name '\PROCESSED_DATA\COMP\' obsTag '_ROI_EOF.mat']));


    ROI_ADCP = A.ROI_ADCP;
    ROI_EXTRAP = E.ROI_EXTRAP;
    ROI_EOF = F.ROI_EOF;


    %% Find valid ROI events

    valid = ~arrayfun( ...
        @(x) isempty(x.avg_V), ...
        ROI_ADCP);


    nValid = sum(valid(:));

    if nValid == 0
        continue
    end


    %% Drifter data

    drifter_V = [ ...
        drifter_V
        [ROI_ADCP(valid).avg_V]' ];

    drifter_ve = [ ...
        drifter_ve
        [ROI_ADCP(valid).avg_ve]' ];

    drifter_vn = [ ...
        drifter_vn
        [ROI_ADCP(valid).avg_vn]' ];


    %% Raw ADCP

    ADCP_V = [ ...
        ADCP_V
        [ROI_ADCP(valid).V_topbin]' ];

    ADCP_VE = [ ...
        ADCP_VE
        [ROI_ADCP(valid).VE_topbin]' ];

    ADCP_VN = [ ...
        ADCP_VN
        [ROI_ADCP(valid).VN_topbin]' ];


    %% Extrapolated

    EXTRAP_V = [ ...
        EXTRAP_V
        [ROI_EXTRAP(valid).V_topbin]' ];

    EXTRAP_VE = [ ...
        EXTRAP_VE
        [ROI_EXTRAP(valid).VE_topbin]' ];

    EXTRAP_VN = [ ...
        EXTRAP_VN
        [ROI_EXTRAP(valid).VN_topbin]' ];


    %% EOF

    EOF_V = [ ...
        EOF_V
        [ROI_EOF(valid).V_topbin]' ];

    EOF_VE = [ ...
        EOF_VE
        [ROI_EOF(valid).VE_topbin]' ];

    EOF_VN = [ ...
        EOF_VN
        [ROI_EOF(valid).VN_topbin]' ];


    %% Keep observation-period identity

    ObsName = [ ...
        ObsName
        repmat(string(names{n}),nValid,1) ];

end


%% Store combined output

ALL.drifter_V  = drifter_V;
ALL.drifter_ve = drifter_ve;
ALL.drifter_vn = drifter_vn;

ALL.ADCP.V  = ADCP_V;
ALL.ADCP.VE = ADCP_VE;
ALL.ADCP.VN = ADCP_VN;

ALL.EXTRAP.V  = EXTRAP_V;
ALL.EXTRAP.VE = EXTRAP_VE;
ALL.EXTRAP.VN = EXTRAP_VN;

ALL.EOF.V  = EOF_V;
ALL.EOF.VE = EOF_VE;
ALL.EOF.VN = EOF_VN;

ALL.ObsName = ObsName;


%% Put products into plotting structure

PlotData(1).name = 'ADCP OBSERVATIONS';
PlotData(1).file = 'adcp_all_obs';
PlotData(1).V  = ADCP_V;
PlotData(1).VE = ADCP_VE;
PlotData(1).VN = ADCP_VN;

PlotData(2).name = 'EXTRAPOLATED OBSERVATIONS';
PlotData(2).file = 'extrap_all_obs';
PlotData(2).V  = EXTRAP_V;
PlotData(2).VE = EXTRAP_VE;
PlotData(2).VN = EXTRAP_VN;

PlotData(3).name = 'EOF';
PlotData(3).file = 'eof_all_obs';
PlotData(3).V  = EOF_V;
PlotData(3).VE = EOF_VE;
PlotData(3).VN = EOF_VN;


%% Plot each product

for p = 1:numel(PlotData)

    f = figure;

    scatter( ...
        drifter_V, ...
        PlotData(p).V, ...
        36,'r','o','filled');

    hold on

    scatter( ...
        drifter_ve, ...
        PlotData(p).VE, ...
        36,'b','s','filled');

    scatter( ...
        drifter_vn, ...
        PlotData(p).VN, ...
        36,'g','^','filled');

    plot([-0.4 0.4],[-0.4 0.4],'k');

    grid on

    legend( ...
        'V','u','v','1:1 line', ...
        'Location','best');

    xlabel('Drifter Velocities (m/s)');
    ylabel('ADCP Velocities (m/s)');
    
    %% Bias
    
    Bias_V = mean(drifter_V - PlotData(p).V,'omitnan');
    
    Bias_VE = mean(drifter_ve - PlotData(p).VE,'omitnan');
    
    Bias_VN = mean(drifter_vn - PlotData(p).VN,'omitnan');
    
    txt = sprintf('Bias V = %.3f m/s \nBias u = %.3f m/s \nBias v = %.3f m/s', ...
        Bias_V, Bias_VE, Bias_VN);
    
    text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment','top','FontSize',10,'BackgroundColor','w','EdgeColor','k');

    title([PlotData(p).name ' - ALL OBSERVATION PERIODS']);


%     exportgraphics( ...
%         f, ...
%         fullfile( ...
%         cfg.out.comp_figures, ...
%         ['_' PlotData(p).file '.pdf']), ...
%         'Resolution',300);
% 
%     exportgraphics( ...
%         f, ...
%         fullfile( ...
%         cfg.out.comp_figures, ...
%         ['_' PlotData(p).file '.png']), ...
%         'Resolution',300);

end


%% All products together

all_drifter_V = repmat(drifter_V,3,1);
all_drifter_ve = repmat(drifter_ve,3,1);
all_drifter_vn = repmat(drifter_vn,3,1);

all_adcp_V = [ ...
    ADCP_V
    EXTRAP_V
    EOF_V ];

all_adcp_VE = [ ...
    ADCP_VE
    EXTRAP_VE
    EOF_VE ];

all_adcp_VN = [ ...
    ADCP_VN
    EXTRAP_VN
    EOF_VN ];

all_Bias_V = mean(all_drifter_V - all_adcp_V);

all_Bias_VE = mean(all_drifter_ve - all_adcp_VE);

all_Bias_VN = mean(all_drifter_vn - all_adcp_VN);


f = figure;

scatter( ...
    all_drifter_V, ...
    all_adcp_V, ...
    36,'r','o','filled');

hold on

scatter( ...
    all_drifter_ve, ...
    all_adcp_VE, ...
    36,'b','s','filled');

scatter( ...
    all_drifter_vn, ...
    all_adcp_VN, ...
    36,'g','^','filled');

plot([-0.4 0.4],[-0.4 0.4],'k');

grid on

legend( ...
    'V','u','v','1:1 line', ...
    'Location','best');

xlabel('Drifter Velocities (m/s)');
ylabel('ADCP Velocities (m/s)');

txt = sprintf('Bias V = %.3f m/s \nBias u = %.3f m/s \nBias v = %.3f m/s', ...
        all_Bias_V, all_Bias_VE, all_Bias_VN);
    
text(0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment','top','FontSize',10,'BackgroundColor','w','EdgeColor','k');

title('ALL OBSERVATION PERIODS - ALL PRODUCTS');


% exportgraphics( ...
%     f, ...
%     fullfile( ...
%     cfg.out.comp_figures, ...
%     '_all_obs_all_products.pdf'), ...
%     'Resolution',300);
% 
% exportgraphics( ...
%     f, ...
%     fullfile( ...
%     cfg.out.comp_figures, ...
%     '_all_obs_all_products.png'), ...
%     'Resolution',300);


%% Combined RMSE

ALL.RMSE.ADCP.V = ...
    rms(drifter_V - ADCP_V,'omitnan');

ALL.RMSE.ADCP.VE = ...
    rms(drifter_ve - ADCP_VE,'omitnan');

ALL.RMSE.ADCP.VN = ...
    rms(drifter_vn - ADCP_VN,'omitnan');


ALL.RMSE.EXTRAP.V = ...
    rms(drifter_V - EXTRAP_V,'omitnan');

ALL.RMSE.EXTRAP.VE = ...
    rms(drifter_ve - EXTRAP_VE,'omitnan');

ALL.RMSE.EXTRAP.VN = ...
    rms(drifter_vn - EXTRAP_VN,'omitnan');


ALL.RMSE.EOF.V = ...
    rms(drifter_V - EOF_V,'omitnan');

ALL.RMSE.EOF.VE = ...
    rms(drifter_ve - EOF_VE,'omitnan');

ALL.RMSE.EOF.VN = ...
    rms(drifter_vn - EOF_VN,'omitnan');


ALL.RMSE.ALL.V = ...
    rms(all_drifter_V - all_adcp_V,'omitnan');

ALL.RMSE.ALL.VE = ...
    rms(all_drifter_ve - all_adcp_VE,'omitnan');

ALL.RMSE.ALL.VN = ...
    rms(all_drifter_vn - all_adcp_VN,'omitnan');

%% RMSE table
RMSE_values = [ ...
    ALL.RMSE.ADCP.V,   ALL.RMSE.ADCP.VE,   ALL.RMSE.ADCP.VN;
    ALL.RMSE.EXTRAP.V, ALL.RMSE.EXTRAP.VE, ALL.RMSE.EXTRAP.VN;
    ALL.RMSE.EOF.V,    ALL.RMSE.EOF.VE,    ALL.RMSE.EOF.VN;
    ALL.RMSE.ALL.V,    ALL.RMSE.ALL.VE,    ALL.RMSE.ALL.VN];

RMSE_Table = array2table(RMSE_values);

RMSE_Table.Properties.VariableNames = {'V','VE','VN'};

RMSE_Table.Properties.RowNames = ...
    {'ADCP','EXTRAP','EOF','ALL_PRODUCTS'};

disp(RMSE_Table)

end