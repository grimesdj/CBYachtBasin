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
% get time of first/last observation...
tlimits = [adcp.time(1) adcp.time(end)];
%
%
% load WB tide data
% we need a few things to retreive waterlevel data:
% 1) the tide gauge's identifier
gaugeID   = '8658163';
% 2) the start/end date for the data-set
startDate = datestr(tlimits(1),'yyyymmdd');%'20230701';%'20220701';%
endDate   = datestr(tlimits(2),'yyyymmdd');%'20230801';%'20230701';%
% 3) the desired variable/interval/units/formatting,
%    see string after "product" below. See also:
%    https://api.tidesandcurrents.noaa.gov/api/prod/
% 4) the URL syntax to queiry data using their API (application
%    programming interface)--a method for two computers to communicate,
%    here it's your computer talking to the file-system housing the data.
%    I'm using either product=water_level or product=hourly_height
url = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=',startDate,'&end_date=',endDate,'&station=',gaugeID,'&product=water_level&interval=6&datum=NAVD&time_zone=gmt&units=metric&format=csv'];
%
%
% we'll look at two methods to read in the ascii formatted data. We'll time each using tic/toc.
% create a filename to save the data on hard-drive:
fileName = [inputDir,filesep,'water_levels_',gaugeID,'_',startDate,'_',endDate,'.csv'];
websave(fileName,url);
%
%
data = readtable(fileName);
time       = table2array(data(:,1));
waterLevel = table2array(data(:,2));
waterSigma = table2array(data(:,3));
nSigma     = 6;
%
% trim WL-data to ADCP deployment period:
idx       = find(time>datetime(adcp.time(1),'convertFrom','datenum') & time<datetime(adcp.time(end),'convertFrom','datenum'));
time      = time(idx);
waterLevel=waterLevel(idx);
waterSigma=waterSigma(idx);
N         = length(time);
dt        = seconds(time(2)-time(1));
%
%
% map to WL-time vector:
basinWaterLevel = interp1(adcp.time, adcp.pressure-mean(adcp.pressure), datenum(time));
%
% fill nans:
flag = isnan(basinWaterLevel);
basinWaterLevel(flag)=interp1(time(~flag),basinWaterLevel(~flag),time(flag),'linear','extrap');
%
%
% apply a 30-min filter:
Nf = 7200/dt; if iseven(Nf), Nf=Nf+1; end
flt= hanning(Nf); flt = flt/sum(flt);
basinWaterLevel = conv(basinWaterLevel,flt,'same');
waterLevel      = conv(waterLevel,flt,'same');
%
fig1 = figure; plot(time,waterLevel,'-k',time,basinWaterLevel,'-r','linewidth',2)
%
% 6 hour max lag
s12 = std(waterLevel).*std(basinWaterLevel);
maxLag = 600;
[R12,lags] = xcorr(waterLevel,basinWaterLevel,maxLag);
R12 = R12/(N-1)/s12;
%
%
[maxR,lag] = max(R12);
time_lag   = lags(lag)*dt
%
fig2 = figure; plot(lags*dt,R12)
xlabel('$\Delta t$ [s]','interpreter','latex')
ylabel('$R_{1,2}$ [n/a]','interpreter','latex')
%
%
figure(fig1)
hold on,plot(time+seconds(time_lag),basinWaterLevel,'--g','linewidth',1.5)
title(sprintf('$\\tau=%d$',abs(time_lag)))
set(gca,'xlim',mean(time) + hours([-24 24]),'tickdir','out','ticklabelinterpreter','latex')
leg = legend('WB','CB','CB($t-\tau$)','AutoUpdate','off');
set(leg,'interpreter','latex')
figName = [figDir,filesep,'tide_phase_lag.png'];
exportgraphics(fig1,figName)
%
%
%
% use min_max to estimate lag between high/low tides:
[idx_WB_max,idx_WB_min,cnt_WB_max,cnt_WB_min] = min_max(waterLevel);
[idx_CB_max,idx_CB_min,cnt_CB_max,cnt_CB_min] = min_max(basinWaterLevel);
% 
plot(time(idx_WB_max),waterLevel(idx_WB_max),'ok','markerfacecolor','c')
plot(time(idx_CB_max),basinWaterLevel(idx_CB_max),'or','markerfacecolor','c')

plot(time(idx_WB_min),waterLevel(idx_WB_min),'ok','markerfacecolor','y')
plot(time(idx_CB_min),basinWaterLevel(idx_CB_min),'or','markerfacecolor','y')

% get closest WB high tide to each CB high-tide...
dt_max = idx_WB_max'-idx_CB_max;
[dt_max_WB_CB,idx_CB_to_WB] = min(abs( dt_max ),[],2);
d_idx_max = -(idx_WB_max-idx_CB_max(idx_CB_to_WB));
dt_max = dt*d_idx_max/60;
% remove outliers
dt_max(abs(dt_max)>120)=[];
dt_high_tide = mean(dt_max)

% get closest WB low tide to each low high-tide...
dt_min = idx_WB_min'-idx_CB_min;
[dt_min_WB_CB,idx_CB_to_WB] = min(abs( dt_min ),[],2);
d_idx_min = -(idx_WB_min-idx_CB_min(idx_CB_to_WB));
dt_min = dt*d_idx_min/60;
% remove outliers
dt_min(abs(dt_min)>120)=[];
dt_low_tide = mean(dt_min)


function [idx_max,idx_min,cnt_max,cnt_min] = min_max(y)
N   = length(y);
idx_max = [];
idx_min = [];
cnt_max = 0;
cnt_min = 0;
for ii = 2:N-1;
    % local minimum:
    if y(ii-1)>y(ii) & y(ii+1)>y(ii)
        cnt_min = cnt_min+1;
        idx_min(cnt_min) = ii;
    % local maximum:
    elseif y(ii-1)<y(ii) & y(ii+1)<y(ii)
        cnt_max = cnt_max+1;
        idx_max(cnt_max) = ii;
    end
end

end
