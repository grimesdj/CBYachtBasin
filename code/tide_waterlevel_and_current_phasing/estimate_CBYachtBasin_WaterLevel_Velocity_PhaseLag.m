function [depth_avg_lag, surf_avg_lag] = estimate_CBYachtBasin_WaterLevel_Velocity_PhaseLag(adcpFile)

figDir  = './figures';
adcp    = load(adcpFile);
%
N = length(adcp.time);
tmp = adcp.VelNorth;
velDepthAvg = mean(tmp,1,'omitnan');

bin_depth = (adcp.dbins-adcp.pressure)';
surf_bins = bin_depth<0 & bin_depth./adcp.pressure'>-1/3;
tmp(~surf_bins)=nan;
velSurfAvg = mean(tmp,1,'omitnan');

basinWaterLevel = adcp.pressure-mean(adcp.pressure);

% low-pass w/ 3hr filter:
dt = (adcp.time(2)-adcp.time(1))*86400;
Nf = round(10800/dt); if iseven(Nf), Nf=Nf+1; end
flt= hanning(Nf); flt = flt/sum(flt);

basinWaterLevel = conv(basinWaterLevel,flt,'same');
velDepthAvg     = conv(velDepthAvg    ,flt,'same');
velSurfAvg      = conv(velSurfAvg     ,flt,'same');


% 6 hour max lag
s12 = std(velDepthAvg).*std(basinWaterLevel);
maxLag = 36;
[R12,lags] = xcorr(basinWaterLevel,velDepthAvg,maxLag);
R12 = R12/(N-1)/s12;
%
[maxR,lag] = max(R12);
depth_avg_lag   = lags(lag)*dt/3600;% minutes

s12 = std(velSurfAvg).*std(basinWaterLevel);
maxLag = 36;
[R12,lags] = xcorr(basinWaterLevel,velSurfAvg,maxLag);
R12 = R12/(N-1)/s12;
%
[maxR,lag] = max(R12);
surf_avg_lag   = lags(lag)*dt/60;% minutes

fig1 = figure;
ax1  = axes;
plot(datetime(adcp.time,'convertfrom','datenum'),basinWaterLevel,'-k','linewidth',2)
ylabel('$\eta$ [m]','interpreter','latex')
set(ax1,'fontsize',15,'tickdir','out','ticklabelinterpreter','latex','box','off')
% $$$ title( sprintf('$\\tau=%2.2f~[h]$',time_lag),'interpreter','latex')
pos1 = get(ax1,'position');
pos1(3) = 0.95*pos1(3)
set(ax1,'position',pos1,'plotboxaspectratio',[1 0.5 1],'ylim',[-1.0 1.5])
%
ax2 = axes;
plot(datetime(adcp.time,'convertfrom','datenum'),velDepthAvg,'-b',datetime(adcp.time,'convertfrom','datenum'),velSurfAvg,'-r','linewidth',2)
legend('Depth Avg.', 'Surface')
% hold on,
% plot(datetime(adcp.time+time_lag/24,'convertfrom','datenum'),adcp.VelNorth(inds),'--g','linewidth',1.5)
ylabel('$v$ [m/s]','interpreter','latex')
set(ax2,'fontsize',16,'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','color','none','ycolor','r','box','off','xtick',[])
set(ax2,'position',pos1,'plotboxaspectratio',[1 0.5 1]);
linkaxes([ax1, ax2],'x')
time = datetime(adcp.time,'convertfrom','datenum');
if month(mean(time))==5
    lims = [datetime('26-May-2025 00:00:00') datetime('28-May-2025 00:00:00')];
else
    lims = mean(time) + hours([-24 24]);
end
set(ax2,'xlim',lims)
figName = [figDir,filesep,'surface_velocity_lag.png'];
exportgraphics(fig1,figName)

title(sprintf('$\\tau_\\mathrm{avg}=%2.1f$~hr, $\\tau_\\mathrm{surf}=%2.1f$~min.',abs(depth_avg_lag),abs(surf_avg_lag)))
str = split(adcpFile,filesep);
figName = [figDir,filesep,'velocity_phase_lag_',str{end-2},'.png'];
exportgraphics(fig1,figName)

end

