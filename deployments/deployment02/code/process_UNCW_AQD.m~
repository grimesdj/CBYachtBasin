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
% Enter time when instrument was in air for pressure offset
atmTime = [datenum('03/6/2025 15:00:00'), datenum('03/6/2025 16:45:00')];
depTime = [datenum('03/6/2025 19:30:00'), datenum('04/1/2025 11:30:00')];
% Enter path to save figures
figDir = [inputDir,filesep,'figures',filesep];
if ~exist(figDir,'dir'), eval(['!mkdir -p ',figDir]), end
%
% Enter time-offset (UTC->EDT) tos = -4 hours
tos = 0;
%
% returns structure A with all aquadopp data
load_AQD_data
save([outputDir,'/',outputName,'.mat'],'-struct','A')
%
%
% plot some stuff
if ~exist('atmTime','var')
    disp('pick 2 points bounding when out of water for ATM pressure offset')
    plot(pressure)
    l = ginput(2);
    l = round(l(:,1));
    atmTime = [A.time(l(1)), A.time(l(2))];
    fprintf('atmTime = \n')
    fprintf('%s --- %s', datestr(atmTime(1)), datestr(atmTime(2)));
else
    l = find(A.time>=atmTime(1) & A.time<=atmTime(2));
end
A.pressureOffset = mean(pressure(l(1):l(2)));
%
if ~exist('depTime','var')
    % now trim the data to when it was in the water
    disp('pick start/end points of deployment')
    l = ginput(2);
    l = round(l(:,1));
    depTime = [A.time(l(1)), A.time(l(2))];
    fprintf('depTime = \n')
    fprintf('%s --- %s', datestr(depTime(1)), datestr(depTime(2)));
else
    valid = find(A.time>=depTime(1) & A.time<=depTime(2));
end
%
if corrFlag
    vars  = {'time','volt','seconds','sspeed','heading','pitch','roll','pressure','temperature','a1','a2','a3','v1','v2','v3','c1','c2','c3','b1','b2','b3','east','north','up'};
else
    vars  = {'time','volt','seconds','sspeed','heading','pitch','roll','pressure','temperature','a1','a2','a3','v1','v2','v3','b1','b2','b3','east','north','up'};
end
for jj = 1:length(vars)
    eval(['A.',vars{jj},' = A.',vars{jj},'(valid,:);'])
end
nsamples = length(valid);
%
A.maxRange = (A.pressure-A.pressureOffset).*cosd(20)-1*binsize;
ylims      = [0 max(min(max(A.maxRange),max(A.dbins)),1)];
dum1       = A.maxRange.*ones(1,nbins);
dum2       = ones(nsamples,1)*A.dbins;
qcFlag0    =  (dum2<=dum1);
if corrFlag
    A.qcFlag   =  double( (dum2<=dum1) & min(A.a1,min(A.a2,A.a3))>75 & min(A.c1,min(A.c2,A.c3))>30 );
else
    A.qcFlag   =  double( (dum2<=dum1) & min(A.a1,min(A.a2,A.a3))>75 );
end
%
time = datetime(A.time,'convertFrom','datenum');
%
% use acceleration and jolt to filter bad data
u1   = A.b1;
d1   = gradientDG(u1)/dt;
dd1  = gradientDG(d1)/dt;
u2   = A.b2;
d2   = gradientDG(u2)/dt;
dd2  = gradientDG(d2)/dt;
u3   = A.b3;
d3   = gradientDG(u3)/dt;
dd3  = gradientDG(d3)/dt;
%
r01  =  nanstd(u1(:));
r02  =  nanstd(u2(:));
r03  =  nanstd(u3(:));
R0   = (u1./r01).^2 + (u2./r02).^2 + (u3./r03).^2;
%
r11  = nanstd(d1(:));
r12  = nanstd(d2(:));
r13  = nanstd(d3(:));
R1   = (d1./r11).^2 + (d2./r12).^2 + (d3./r13).^2;
%
r21  = nanstd(dd1(:));
r22  = nanstd(dd2(:));
r23  = nanstd(dd3(:));
R2   = (dd1./r21).^2 + (dd2./r22).^2 + (dd3./r23).^2;
%
valid = (R0<15);% & (R1<5) & (R2<10);
A.qcFlag = A.qcFlag & valid;
%
% make a few quick plots
fig0 = figure;
ax1 = subplot(2,1,1);
plot(time,A.temperature)
ylabel(ax1,'$T$ [$^\circ$]','interpreter','latex')
set(ax1,'xticklabel','','ticklabelinterpreter','latex','tickdir','out')
ax2 = subplot(2,1,2);
plot(time,A.pressure)
ylabel(ax2,'$P$ [m]','interpreter','latex')
xlabel(ax2,'time [s]','interpreter','latex')
set(ax2,'ticklabelinterpreter','latex','tickdir','out')
figName = [figDir,'/',inputFile,'_temperature_pressure.png'];
exportgraphics(fig0,figName)
%
%
% $$$ % make a few quick plots
% $$$ fig0 = figure;
% $$$ plot(time,A.volt)
% $$$ ylabel('[volts]','interpreter','latex')
% $$$ xlabel('time [s]','interpreter','latex')
% $$$ set(gca,'ticklabelinterpreter','latex','tickdir','out')
% $$$ figName = [figDir,'/',inputFile,'_battery_voltage.png'];
% $$$ exportgraphics(fig0,figName)
%
%
% quick convolution running mean filter
np1 = round(0.3/A.config.binSize);% 10 cm vertical 
np2 = 31;% 5min for 1Hz data
f1 = hamming(np1);f1 = f1./sum(f1);
f2 = hamming(np2);f2 = f2./sum(f2);
% do a nan-mean filter, keep track of normalization
on = conv2(f1,f2,A.qcFlag','same');
%
A1 = conv2(f1,f2,(A.a1.*A.qcFlag)','same')./on;
A2 = conv2(f1,f2,(A.a2.*A.qcFlag)','same')./on;
A3 = conv2(f1,f2,(A.a3.*A.qcFlag)','same')./on;
%
fig1 = figure;
ax1 = subplot(3,1,1);
imagesc(time,A.dbins',A1.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'X')
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
title(ax1,'Amplitude')
%
ax2 = subplot(3,1,2);
imagesc(time,A.dbins',A2.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Y')
ylabel('mab','interpreter','latex')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
%
ax3 = subplot(3,1,3);
imagesc(time,A.dbins',A3.*qcFlag0'),caxis([100 180]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir,'/',inputFile,'_amplitude.png'];
exportgraphics(fig1,figName)
%
%
if corrFlag
C1 = conv2(f1,f2,(A.c1.*A.qcFlag)','same')./on;
C2 = conv2(f1,f2,(A.c2.*A.qcFlag)','same')./on;
C3 = conv2(f1,f2,(A.c3.*A.qcFlag)','same')./on;
%
fig1 = figure;
ax1 = subplot(3,1,1);
imagesc(time,A.dbins',C1.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'X')
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
title(ax1,'Correlation')
%
ax2 = subplot(3,1,2);
imagesc(time,A.dbins',C2.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Y')
ylabel('mab','interpreter','latex')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
%
ax3 = subplot(3,1,3);
imagesc(time,A.dbins',C3.*qcFlag0'),caxis([0 100]),colormap(cmocean('thermal')),colorbar
text(time(1),ylims(2),'Z')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir,'/',inputFile,'_correlation.png'];
exportgraphics(fig1,figName)
end
%
%
% Now plot currents (ENU)
V1 = conv2(f1,f2,(A.east.*A.qcFlag)','same')./on;
V2 = conv2(f1,f2,(A.north.*A.qcFlag)','same')./on;
V3 = conv2(f1,f2,(A.up.*A.qcFlag)','same')./on;
%
uv_caxis = 5*[-1 1]*max(std(V1(:),'omitnan'),std(V2(:),'omitnan'));
up_caxis = [-0.1 0.1];
fig2 = figure;
ax1 = subplot(3,1,1);
imagesc(time,A.dbins',V1.*qcFlag0'),caxis(uv_caxis),colormap(cmocean('balance')),colorbar
text(time(1),ylims(2)+0.05*range(ylims),'East')
%
set(ax1,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
ax2 = subplot(3,1,2);
imagesc(time,A.dbins',V2.*qcFlag0'),caxis(uv_caxis),colormap(cmocean('balance')),colorbar
text(time(1),ylims(2)+0.05*range(ylims),'North')
set(ax2,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
ylabel('mab','interpreter','latex')
%
ax3 = subplot(3,1,3);
imagesc(time,A.dbins',V3.*qcFlag0'),caxis(up_caxis),colormap(cmocean('balance')),colorbar
text(time(1),ylims(2)+0.05*range(ylims),'Up')
set(ax3,'ydir','normal','ticklabelinterpreter','latex','ylim',ylims)
xlabel('time [s]','interpreter','latex')
figName = [figDir,'/',inputFile,'_velocity.png'];
exportgraphics(fig2,figName)
%
% get the time-averaged current.
% Note, this is not in depth normalized (sigma) coordinates.
% first, nan any values that don't pass QC.
V1(~A.qcFlag')=nan;
V2(~A.qcFlag')=nan;
V3(~A.qcFlag')=nan;
flag = sum(A.qcFlag,1) > 0.50*nsamples;
% Uz = time averegd; U = depth & time averaged
Uz = nanmean(V1,2); U = nanmean(Uz); Uz(~flag)=nan;
Vz = nanmean(V2,2); V = nanmean(Vz); Vz(~flag)=nan;
Wz = nanmean(V3,2); W = nanmean(Wz); Wz(~flag)=nan;
%
fig3 = figure;
plot(Uz,A.dbins,'.-r',Vz,A.dbins,'.-b',Wz,A.dbins,'.-k')
hh = legend('East','North','Up');
set(hh,'fontsize',9)
xlabel('m/s','interpreter','latex')
ylabel('mab','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out')
figName = [figDir,'/',inputFile,'_mean_velocity.png'];
exportgraphics(fig3,figName)
%
fig4 = figure;
pos  = get(fig4,'position');
pos(3:4) = pos(3:4).*[0.5, 1];
set(fig4, 'position',pos)
plot(Uz,A.dbins,'.-k','linewidth',2)
set(hh,'fontsize',9)
xlabel('Avg. East Velocity [m/s]','interpreter','latex')
ylabel('meters-above-bottom','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 1.5],'xlim',[0 1])
figName = [figDir,'/',inputFile,'_mean_alongshore_velocity.png'];
exportgraphics(fig4,figName)
%
A.VelEast = V1;
A.VelNorth= V2;
A.VelZ   = V3;
%
A.VelXAvg = Uz;
A.VelYAvg= Vz;
A.VelZAvg   = Wz;
%
save([L0Dir,'/',L0Name,'.mat'],'-struct','A')
%
% add the config info to the structure A to quick save as netcdf4
fieldNames = fields(A.config);
originalFields = fields(A);
%
for j = 1:length(fieldNames)
 A.(fieldNames{j}) = A.config.(fieldNames{j});
end
A = orderfields(A,cat(1,fieldNames,originalFields));
ncfile = [L0Dir,'/',L0Name,'.nc'];
if exist(ncfile,'file')
    eval(['!rm ',ncfile])
end
struct2nc(A,ncfile,'NETCDF4');
%
% $$$ %
% $$$ v1bar = nansum( A.v1.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ v2bar = nansum( A.v2.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ v3bar = nansum( A.v3.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ % rotate 60 degrees:
% $$$ vrot = [cosd(60) sind(60);-sind(60) cosd(60)]*[v1bar';v2bar'];
% $$$ v1rot = vrot(1,:)';
% $$$ v2rot = vrot(2,:)';
% $$$ %
% $$$ %
% $$$ dv    = 0.01;
% $$$ vbins = [-1:dv:1];
% $$$ v1H   = hist(v1rot,vbins);
% $$$ v2H   = hist(v2rot,vbins);
% $$$ v3H   = hist(v3bar,vbins);
% $$$ %
% $$$ figure, plot(v1rot,v2rot,'.')
% $$$ hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
% $$$ xlabel('$\bar{v}_X$','interpreter','latex')
% $$$ ylabel('$\bar{v}_Y$','interpreter','latex')
% $$$ figName = [figDir,'/',inputFile,'_depth_avg_XYvel_checker_pattern.png'];
% $$$ exportgraphics(gcf,figName)
% $$$ %
% $$$ %
% $$$ figure,
% $$$ subplot(3,1,1)
% $$$ bar(vbins, v1H./length(v1bar))
% $$$ ylabel('$p(\bar{v}_{X''})$','interpreter','latex')
% $$$ subplot(3,1,2)
% $$$ bar(vbins, v2H./length(v1bar))
% $$$ ylabel('$p(\bar{v}''_{Y''})$','interpreter','latex')
% $$$ subplot(3,1,3)
% $$$ bar(vbins, v3H./length(v1bar))
% $$$ ylabel('$p(\bar{v}_Z)$','interpreter','latex')
% $$$ xlabel('velocity')
% $$$ figName = [figDir,'/',inputFile,'_depth_avg_vel_XYrot60deg_histograms.png'];
% $$$ exportgraphics(gcf,figName)
% $$$ %
% $$$ %
% $$$ %
% $$$ %
% $$$ b1bar = nansum( A.b1.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ b2bar = nansum( A.b2.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ b3bar = nansum( A.b3.*A.qcFlag ,2)./sum(A.qcFlag,2);
% $$$ figure, plot(b1bar,b2bar,'.')
% $$$ hold on, plot(1.5*[-cosd(60) cosd(60)],1.5*[-sind(60) sind(60)],'--r')
% $$$ xlabel('$\bar{v}_X$','interpreter','latex')
% $$$ ylabel('$\bar{v}_Y$','interpreter','latex')
% $$$ figName = [figDir,'/',inputFile,'_depth_avg_Beam1Beam2_checker_pattern.png'];
% $$$ exportgraphics(gcf,figName)
