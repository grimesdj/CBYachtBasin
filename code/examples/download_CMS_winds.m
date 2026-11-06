clear all
close all
% $$$ dataPath = dir('~/PHY477/data/REUdataJune2022.txt');
% $$$ fid      = fopen([dataPath(1).folder,filesep,dataPath(1).name]);
% $$$ nHeader  = 6;
% $$$ tEmpty   = '"NAN"';
% $$$ format   = ['%s %*u %*f %f ',repmat('%*s ',1,3),'%f %*s %f %*s %f ', repmat('%*s ',1,5), '%f %f ',repmat('%*s ',1,6),'%f ',repmat('%*s ',1,3)];
% $$$ % keeping time, temp, humidity, pressure, wind spd U, wind spd V, solar,
% $$$ data     = textscan(fid,format,'HeaderLines',nHeader,'Delimiter',',','treatasempty',tEmpty);
% $$$ 
% $$$ time     = datenum(char(data{1}));
% $$$ dtime    = datetime(time,'ConvertFrom','datenum');
% $$$ temp     = data{2};
% $$$ humid    = data{3};
% $$$ pres     = data{4};
% $$$ wSpeed   = data{5};
% $$$ uWind    = data{6};
% $$$ vWind    = data{7};
% $$$ solar    = data{8};
% $$$
tStart = '2026-09-01T00:00:00';
tEnd   = '2026-09-30T23:59:59';
url = ['https://cormp.org/services/download.php?time=',tStart,'-04:00/',tEnd,'-04:00&tz=UTC&standard=false&output=csv&pretty=true&parameters[]=CMSDOCK+Air+temperature&parameters[]=CMSDOCK+Solar+radiation&parameters[]=CMSDOCK+Wind+speed&parameters[]=CMSDOCK+Wind+from+direction'];
str = urlread(url);
%
nHeader = 10;
format  = '%s %f %*s %f %*s %f %*s %f %*s';
data    = textscan(str,format,'HeaderLines',nHeader,'Delimiter',',');
time    = datetime(data{1});
temp    = (data{2}-32)*5/9;% Celsius
solar   = data{3};
wSpd    = data{4}/2.237;% mps
wDir    = data{5};
uWind = wSpd.*cosd(wDir);
N     = length(time);
