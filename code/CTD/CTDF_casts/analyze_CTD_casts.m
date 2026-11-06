clear all
close all
% -------------------------------------------------
%% Analyze casts along navigation channel
% -------------------------------------------------
%% 0) list of input L0 files: ADCP (pressure) & CTD
rootDIR ='~/OneDriveUNCW/DATA/CB_YachtBasin/'
CTDFfiles= {fullfile(rootDIR,'JuneJamboree','PROCESSED_DATA','CTD','230442_CTDF_L0.mat');...
            fullfile(rootDIR,'SummerFlood' ,'PROCESSED_DATA','CTD','230442_CTDF_L0.mat');...
            fullfile(rootDIR,'FallFlood' ,'PROCESSED_DATA','RBR','230442_20260930_CTDF_L0.mat');...
            fullfile(rootDIR,'FallFlood' ,'PROCESSED_DATA','RBR','230442_20261028_CTDF_L0.mat');...
            fullfile(rootDIR,'FallFlood' ,'PROCESSED_DATA','RBR','230442_20261104_CTDF_L0.mat')};
ADCPfiles= {fullfile(rootDIR,'JuneJamboree','PROCESSED_DATA_v2','ADCP','JuneJamboree_NorthADCP_adcp_level_1.mat');...
            fullfile(rootDIR,'SummerFlood','PROCESSED_DATA_v2','ADCP','SummerFlood_NorthADCP_adcp_level_1.mat');...
            fullfile(rootDIR,'FallFlood','PROCESSED_DATA_v2','ADCP','FallFlood_NorthADCP_adcp_level_1.mat');...            
            fullfile(rootDIR,'FallFlood','PROCESSED_DATA_v2','ADCP','FallFlood_NorthADCP_adcp_level_1.mat');...
            fullfile(rootDIR,'FallFlood','PROCESSED_DATA_v2','ADCP','FallFlood_NorthADCP_adcp_level_1.mat')};
% 
% ---------------------------------
%% 1) define ROI for casts:
% ---------------------------------
[perim,Lat0,Lon0,LatInterp,LonInterp,East0,North0,EastInterp,NorthInterp,distInterp] = define_ROI_for_casts(CTDFfiles);
%
% 1.1) pre-allocate regular along-channel and depth grid
dInterp = [0.1:0.1:6]';
[xx0,dd0] = meshgrid(distInterp,dInterp);
%
%
% ---------------------------------
%% 2) load indivual dates and
%     interpolate to 
% ---------------------------------
% for ii=1:length(CTDFfiles)
    ii=2;
    CTDFfile = CTDFfiles{ii};
    ADCPfile = ADCPfiles{ii};
    %
    % 2.0) load CTDF & ADCP
    CTD = load(CTDFfile);
    load(ADCPfile);
    ADCP= adcp_level_1;
    % 2.1) get the time/pressure grid for ADCP
    %
    % 2.2) get casts in ROI
    inROI = inpolygon(CTD.Latitude,CTD.Longitude,perim(:,1),perim(:,2));    
    %
    fig0 = figure; geoplot(perim(:,1),perim(:,2))
    hold on,geoplot(LatInterp, LonInterp,'-m')
    hold on,geoplot(CTD.Latitude, CTD.Longitude,'*r')
    hold on,geoplot(CTD.Latitude(inROI), CTD.Longitude(inROI),'*g')
    %
    % 2.3) map to nearest point on alongchannel curve (LatInterp,LonInterp)
    Latitude = CTD.Latitude(inROI)';
    Longitude= CTD.Longitude(inROI)';
    [Pnear, distNear, iseg, tNear, ivert] = nearest_point_on_curve([Latitude Longitude], [LatInterp, LonInterp]);
    %
    dumLat = [Pnear(:,1)'; CTD.Latitude(inROI); nan*iseg'];
    dumLon = [Pnear(:,2)'; CTD.Longitude(inROI); nan*iseg'];
% $$$     dumLat = [LatInterp(ivert)'; CTD.Latitude(inROI); nan*iseg'];
% $$$     dumLon = [LonInterp(ivert)'; CTD.Longitude(inROI); nan*iseg'];
    hold on,geoplot(dumLat(:), dumLon(:) ,'-g')    

    % 2.4) interpolate to grid
    % define projection:
    proj = projcrs(32119);
    % Convert to ENU:
    [East,North] = projfwd(proj,Pnear(:,1),Pnear(:,2));
    Distance = sqrt( (East-East0).^2 + (North-North0).^2);
    %
    % 2.5) need to segment the data to individual transects...
    Time        = CTD.Time       (:,inROI);
    Pressure    = CTD.Pressure;
    Temperature = CTD.Temperature(:,inROI);
    Salinity    = CTD.Salinity   (:,inROI);    
    Density     = CTD.Density    (:,inROI);
    CHLa        = CTD.CHLa       (:,inROI);
    NTU         = CTD.NTU        (:,inROI);
    FDOM        = CTD.FDOM       (:,inROI);
    %
    % 2.6) remove surface and bottom bins
    ivalid = Pressure>0.5 & ~isnan(Density);
    rows   = repmat([1:size(ivalid,1)]',1,size(ivalid,2));
    %    rows(ivalid) = 9999;
    last_valid_row = max(rows-~ivalid*9999);
    last_valid_row(last_valid_row==9999) = 3;
    ivalid = Pressure>0.5 & rows<last_valid_row-2;
    %
    Temperature(~ivalid) = nan;
    Salinity   (~ivalid) = nan;    
    Density    (~ivalid) = nan;
    CHLa       (~ivalid) = nan;
    NTU        (~ivalid) = nan;
    FDOM       (~ivalid) = nan;
    %
    fig1 = figure;
    ax1  = subplot(2,1,1);
    plot(datetime(ADCP.time,'convertFrom','datenum'),ADCP.pressure,'*')
    ax2  = subplot(2,1,2);
    plot(datetime(Time,'convertFrom','datenum'),Distance,'*')
    set(ax1,'xlim',get(ax2,'xlim'))
    switch ii
      case 1
        TransectTimes = [datenum('14-Jun-2026 08:15:00'), datenum('14-Jun-2026 11:15:00');...
        %                         datenum('14-Jun-2026 09:15:00'), datenum('14-Jun-2026 10:15:00');...
                         datenum('14-Jun-2026 11:15:00'), datenum('14-Jun-2026 14:15:00');...
% $$$                          datenum('14-Jun-2026 13:15:00'), datenum('14-Jun-2026 14:15:00');...
                         %                         datenum('14-Jun-2026 14:15:00'), datenum('14-Jun-2026 15:15:00');...
                         datenum('14-Jun-2026 14:15:00'), datenum('14-Jun-2026 15:30:00');...
                         datenum('14-Jun-2026 15:30:00'), datenum('14-Jun-2026 17:00:00');...                         
                         %                         datenum('14-Jun-2026 16:30:00'), datenum('14-Jun-2026 17:30:00');...
                         datenum('14-Jun-2026 17:00:00'), datenum('14-Jun-2026 20:30:00');...
                         %                         datenum('14-Jun-2026 18:45:00'), datenum('14-Jun-2026 20:30:00')...
                        ];
        flag = [61 90];
      case 2
        TransectTimes = [datenum('12-Aug-2026 15:00:00'), datenum('12-Aug-2026 16:30:00');...
        %                         datenum('14-Jun-2026 09:15:00'), datenum('14-Jun-2026 10:15:00');...
                         datenum('12-Aug-2026 16:30:00'), datenum('12-Aug-2026 18:00:00');...
% $$$                          datenum('14-Jun-2026 13:15:00'), datenum('14-Jun-2026 14:15:00');...
                         %                         datenum('14-Jun-2026 14:15:00'), datenum('14-Jun-2026 15:15:00');...
                         datenum('14-Jun-2026 14:15:00'), datenum('14-Jun-2026 15:30:00');...
                         datenum('14-Jun-2026 15:30:00'), datenum('14-Jun-2026 17:00:00');...                         
                         %                         datenum('14-Jun-2026 16:30:00'), datenum('14-Jun-2026 17:30:00');...
                         datenum('14-Jun-2026 17:00:00'), datenum('14-Jun-2026 20:30:00');...
                         %                         datenum('14-Jun-2026 18:45:00'), datenum('14-Jun-2026 20:30:00')...
                        ];
      case 3
        TransectTimes = [datenum('30-Sep-2026 13:00:00'), datenum('30-Sep-2026 13:45:00');...
                         datenum('30-Sep-2026 13:45:00'), datenum('30-Sep-2026 14:45:00');...
                         datenum('30-Sep-2026 14:45:00'), datenum('30-Sep-2026 15:30:00')...
                        ];
        %% Manually flagged casts:
        flag = 6;
    end
    %
    if ~isempty(flag)
        Time         (flag)  = [];
        Latitude     (flag)  = [];
        Longitude    (flag)  = [];
        East         (flag)  = [];
        North        (flag)  = [];
        Distance     (flag)  = [];
        Pnear      (flag,:)  = [];
        Temperature(:,flag)  = [];
        Salinity   (:,flag)  = [];
        Density    (:,flag)  = [];
        CHLa       (:,flag)  = [];
        NTU        (:,flag)  = [];
        FDOM       (:,flag)  = [];
    end
    axes(ax1)
    hold on,
    xline(datetime(TransectTimes(:,2),'convertFrom','datenum'),'-r','linewidth',2)
    xline(datetime(TransectTimes(:,1),'convertFrom','datenum'),'--g','linewidth',2)
    axes(ax2)
    hold on,
    xline(datetime(TransectTimes(:,2),'convertFrom','datenum'),'-r','linewidth',2)
    xline(datetime(TransectTimes(:,1),'convertFrom','datenum'),'--g','linewidth',2)
    %
    %
% $$$     tmp    = Density;
% $$$     ivalid = ~isnan(tmp);
% $$$     % create grids, sort distance S-N
% $$$     [Distance_tmp, srt] = sort(Distance');
% $$$     tmp = tmp(:,srt);
% $$$     ivalid = ivalid(:,srt);
% $$$     [xx,dd] = meshgrid(dist_tmp,CTD.Pressure);
% $$$         F = scatteredInterpolant(xx(ivalid),dd(ivalid),tmp(ivalid));
% $$$         DensityInterpAll = F(xx0,dd0);
% $$$     %
    DensityInterp=nan([size(dd0), size(TransectTimes,1)]);
    fig2 = figure;
    NT   = size(TransectTimes,1);
    %% get density range from all casts
    dlims = mean(Density(:),'omitnan') - 1e3 + 1.5*std(Density(:),'omitnan')*[-1 1];
    for jj=1:NT
        idx    = find(Time>=TransectTimes(jj,1) & Time<=TransectTimes(jj,2));
        tmp    = Density(:,idx);
        ivalid = ~isnan(tmp);
        %
        figure(fig0)
        geoscatter(Latitude(idx),Longitude(idx),25,nanmean(Density(:,idx),1),'filled'),colormap(cmocean('thermal')),caxis([1010 1025])        
% $$$         geoscatter(Latitude(idx),Longitude(idx),25,nanmean(Temperature(:,idx),1),'filled'),colormap(cmocean('thermal')),caxis([15 24])        
        %
        if sum(ivalid,'all')/numel(ivalid)<0.25
            continue
        end              
        % create grids, sort distance S-N
        [dist_tmp, srt] = sort(Distance(idx)');
        tmp = tmp(:,srt);
        ivalid = ivalid(:,srt);
        %
        %
        [xx,dd] = meshgrid(dist_tmp,CTD.Pressure);
        F = scatteredInterpolant(xx(ivalid),dd(ivalid),tmp(ivalid),'linear','none');
        DensityInterp(:,:,jj) = F(xx0,dd0);
        %
        %
        figure(fig2)
        ax = subplot(NT,1,jj);
        imagesc(xx0(1,:),dd0(:,1),DensityInterp(:,:,jj)-1e3,'alphadata',~isnan(DensityInterp(:,:,jj))),colormap(cmocean('dense'))
        caxis(dlims)
        cb = colorbar;
        ylabel(cb,'$\sigma_\theta$ [kg/m$^3$]','interpreter','latex','fontsize',10)
        set(cb,'ticklabelinterpreter','latex')
        if jj<NT
            set(ax,'xticklabel',[],'ticklabelinterpreter','latex','tickdir','out')
            ylabel(ax,'depth [m]','interpreter','latex')
        end
    end
    set(ax,'ticklabelinterpreter','latex','tickdir','out')
    xlabel(ax,'Meters North of Municipal Dock','interpreter','latex')
    ylabel(ax,'depth [m]','interpreter','latex')    
% __________
%% 3) Plot individual casts:
%
% 3.1) Temperature (y,z)
%
% 3.2) Salinity (y,z)
%
% 3.3) Density (y,z)
%
% 3.4) NTU (y,z)
%
% 3.5) CDOM (y,z)
%
% 3.7) Chl-a (y,z)
%
% $$$ end;
% ---------
%% 4) Plot compendium by tidal phase:
