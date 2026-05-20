function [adcp_level_1] = load_and_save_adcp(inputDir, inputFile, cfg);


%% Load the ADCP data

fileName  = [inputDir,filesep,inputFile];
% Enter raw output /directory/ and fileName without .mat
outputDir = inputDir;
outputName= [inputFile,'_raw'];
% Enter processed output fileName without .mat
L0Dir   = inputDir;
L0Name  = [inputFile,'_L0'];
load([L0Dir,'/',L0Name,'.mat'])
% % Enter time when instrument was in air for pressure offset
% atmTime = [datenum('03/6/2025 15:00:00'), datenum('03/6/2025 16:45:00')];
% depTime = [datenum('03/6/2025 19:30:00'), datenum('04/1/2025 11:30:00')];
% % Enter path to save figures
figDir = [inputDir,filesep,'figures',filesep];

%% quick convolution running mean filter
np1 = round(0.3/config.binSize);% 30 cm vertical 
np2 = 30*60/config.dt;% 15 data point averaging  % try limiting to 1 - 3 hours
f1 = hamming(np1);f1 = f1./sum(f1);
f2 = hamming(np2);f2 = f2./sum(f2);
% do a nan-mean filter, keep track of normalization
on = conv2(f1,f2,qcFlag','same');
%
A1 = conv2(f1,f2,(a1.*qcFlag)','same')./on;
A2 = conv2(f1,f2,(a2.*qcFlag)','same')./on;
A3 = conv2(f1,f2,(a3.*qcFlag)','same')./on;

%% Get the currents
% V1 = conv2(f1,f2,(east.*qcFlag)','same')./on;
% % V1 = (east.*qcFlag)'
east(~qcFlag) = NaN;
V1 = east';
% V2 = conv2(f1,f2,(north.*qcFlag)','same')./on;
% % V2 = (north.*qcFlag)'
north(~qcFlag) = NaN;
V2 = north';
V3 = conv2(f1,f2,(up.*qcFlag)','same')./on;
% %
uv_caxis = 5*[-1 1]*max(std(V1(:),'omitnan'),std(V2(:),'omitnan'));
up_caxis = [-0.1 0.1];

East_avg  = sum(east.*qcFlag,2,'omitnan')./sum(qcFlag,2);
North_avg = sum(north.*qcFlag,2,'omitnan')./sum(qcFlag,2);

N = numel(East_avg);

%% Plot the data

ylims = [0 6];

fig1 = figure;

ax1 = subplot(2,1,1);
imagesc(time,dbins',V1.*qcFlag','AlphaData',qcFlag'),caxis(uv_caxis),colormap(cmocean('balance')),colorbar
% xlim([datenum(2025,11,07) datenum(2025,11,08)]);
hold on
% xline([datenum(2025,11,07,13,04,0)],'r','Linewidth',2);
% xline([datenum(2025,11,07,13,47,0)],'k','Linewidth',2);
% xline([datenum(2025,11,07,14,40,0)],'g','Linewidth',2);
% xline([datenum(2025,11,07,15,56,0)],'m','Linewidth',2);
plot(time,pressure)
hold off
datetick('x','mm/dd HH:MM','keeplimits','keepticks');
% datetick('x','HH:MM','keeplimits')
% text(time(655),ylims(2)+0.05*range(ylims),'East')
title('East');
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)

ax2 = subplot(2,1,2);
imagesc(time,dbins',V2.*qcFlag','AlphaData',qcFlag'),caxis(uv_caxis),colormap(cmocean('balance')),colorbar
% xlim([datenum(2025,11,07) datenum(2025,11,08)]);
hold on
% xline([datenum(2025,11,07,13,04,0)],'r','Linewidth',2);
% xline([datenum(2025,11,07,13,47,0)],'k','Linewidth',2);
% xline([datenum(2025,11,07,14,40,0)],'g','Linewidth',2);
% xline([datenum(2025,11,07,15,56,0)],'m','Linewidth',2);
plot(time,pressure)
hold off
datetick('x','mm/dd HH:MM','keeplimits','keepticks');
% datetick('x','HH:MM','keeplimits')
title('North');
% text(time(1),ylims(2)+0.05*range(ylims),'North')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
ylabel('mab','interpreter','latex')

%% Save the ADCP data for post processing

% call this L1

ADCP_comp.time = time;
ADCP_comp.V_E = V1;
% ADCP_comp.V_E_raw = V1_1;
ADCP_comp.V_N = V2;
% ADCP_comp.V_N_raw = V2_1;
ADCP_comp.bins = dbins;
ADCP_comp.East_avg = East_avg;
ADCP_comp.North_avg = North_avg;
ADCP_comp.qcFlag = qcFlag;
ADCP_comp.pressure = pressure;
ADCP_comp.east = east;
ADCP_comp.north = north;
ADCP_comp.f1 = f1;
ADCP_comp.f2 = f2;
ADCP_comp.time = time;
ADCP_comp.on = on;
ADCP_comp.uv_caxis = uv_caxis;
ADCP_comp.up_caxis = up_caxis;

%%
adcp_level_1 = ADCP_comp;

exportgraphics(fig1, fullfile(cfg.out.adcp_figures, [cfg.obsTag '_adcp_level_1.png']), 'Resolution', 300);

end