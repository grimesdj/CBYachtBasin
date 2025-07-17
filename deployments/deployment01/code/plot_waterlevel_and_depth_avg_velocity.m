clear all
close all
% Enter input /directory/ and fileName root without file extension
inputDir  = '/Users/derekgrimes/OneDriveUNCW/Data/CB_YachtBasin/MarshMadness/AQD5459/';
inputFile = 'CB_N101';
fileName  = [inputDir,filesep,inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = inputDir;
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = inputDir;
L0Name  = [inputFile,'_L0'];
figDir  = [inputDir,filesep,'figures',filesep];
adcp    = load([L0Dir,filesep,L0Name,'.mat']);
%
%
N = length(adcp.time);
% $$$ % this finds the nearest bin to surface
% $$$ rows    = repmat([1:length(adcp.dbins)]',1,N);
% $$$ rows(isnan(adcp.VelNorth)) = 0;
% $$$ inwater = max(rows,[],1);
% $$$ inds    = sub2ind(size(adcp.VelNorth),inwater,1:N);
velocityDepthAverage = mean(adcp.VelNorth,1,'omitnan');
%
% 6 hour max lag
dt = adcp.time(2)-adcp.time(1);
basinWaterLevel = adcp.pressure-mean(adcp.pressure);
s12 = std(velocityDepthAverage).*std(basinWaterLevel);
maxLag = 600;
[R12,lags] = xcorr(basinWaterLevel,velocityDepthAverage,maxLag);
R12 = R12/(N-1)/s12;
%
%
[maxR,lag] = max(R12);
time_lag   = lags(lag)*dt*24;
%
fig1 = figure;
ax1 = axes;
plot(datetime(adcp.time,'convertfrom','datenum'),adcp.pressure-nanmean(adcp.pressure),'-k','linewidth',2)
ylabel('$\eta$ [m]','interpreter','latex')
set(ax1,'fontsize',15,'tickdir','out','ticklabelinterpreter','latex','box','off')
title( sprintf('$\\tau=%2.2f~[h]$',time_lag),'interpreter','latex')
pos1 = get(ax1,'position');
pos1(3) = 0.95*pos1(3)
set(ax1,'position',pos1)
%
ax2 = axes;
plot(datetime(adcp.time,'convertfrom','datenum'),velocityDepthAverage,'-r','linewidth',2)
hold on,
plot(datetime(adcp.time+time_lag/24,'convertfrom','datenum'),velocityDepthAverage,'--g','linewidth',1.5)
ylabel('$v$ [m/s]','interpreter','latex')
set(ax2,'fontsize',16,'tickdir','out','ticklabelinterpreter','latex','yaxislocation','right','color','none','ycolor','r','box','off','xtick',[])
set(ax2,'position',pos1);
linkaxes([ax1, ax2],'x')
time = datetime(adcp.time,'convertfrom','datenum');
set(ax2,'xlim',mean(time) + hours([-24 24]))
figName = [figDir,filesep,'depth_averaged_velocity_lag.png'];
exportgraphics(fig1,figName)
% $$$ 
% $$$ 
% $$$ fig1 = figure;
% $$$ ax1 = subplot(2,1,1);
% $$$ plot(datetime(adcp.time,'convertfrom','datenum'),adcp.pressure-nanmean(adcp.pressure),'-k','linewidth',2)
% $$$ ylabel('$\eta$ [m]','interpreter','latex')
% $$$ set(gca,'fontsize',15,'tickdir','out','ticklabelinterpreter','latex')
% $$$ ax2 = subplot(2,1,2);
% $$$ plot(datetime(adcp.time,'convertfrom','datenum'),adcp.VelNorth(inds),'-k','linewidth',2)
% $$$ ylabel('$v$ [m/s]','interpreter','latex')
% $$$ set(gca,'fontsize',15,'tickdir','out','ticklabelinterpreter','latex')
% $$$ figName = [figDir,filesep,'surface_velocity_lag.png'];
% $$$ exportgraphics(fig1,figName)
