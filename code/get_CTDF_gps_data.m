function [gps,rawFiles] = get_CTDF_gps_data(releaseNumber,SN)
%
% USAGE: [gps,rawFiles] = get_CTDF_gps_data(releaseNumber,SN,)
%
% get the gps data from multiple devices for pairing with CTDF casts

% get gps info
gpsInfoFile = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/info/release%d_ctdf_cast_gps/gps_ctd_serial_number_info.csv',releaseNumber);
gpsFID = fopen(gpsInfoFile,'r');
gpsDAT = textscan(gpsFID,'%s %d %s','Headerlines',1,'delimiter',',');
%
rootDir  = sprintf('/Users/derekgrimes/Library/CloudStorage/OneDrive-UNC-Wilmington/KELP-vingilote/data/Release%d/',releaseNumber);
rawDir    = [rootDir,filesep,'raw'];
%
% Pair CTDF and GPS data
rawFiles = {};
isPair = ismember(gpsDAT{2},SN);
indPair= find(isPair);
% for each pair, get the GPS data file and read time/lat/lon/etc.
gps = struct([]);
for jj=1:sum(isPair)
    string  = char(gpsDAT{3}(indPair(jj)));
    gpsRoot = sprintf([rawDir,filesep,'%s','*.TXT'],string);
    gpsFiles= dir(gpsRoot);
    gpsPath = [gpsFiles.folder,filesep,gpsFiles.name];
    rawFiles{jj}=gpsPath;
    % get date from filename
    gpsStr = split(gpsFiles.name,{'_','.'});
    gpsDate= gpsStr{3}();
    gpsYear= str2num(gpsDate(1:4)); gpsMon = str2num(gpsDate(5:6)); gpsDay = str2num(gpsDate(7:8));
    gpsTime= gpsStr{4}();
    gpsHH  = str2num(gpsTime(1:2)); gpsMM = str2num(gpsTime(3:4)); gpsSS = str2num(gpsTime(5:6));
    %
    [gps_hour, gps_minute, gps_second, velocity, numsats, height, lat, lon, quality, geoidsep, pdop, hdop, vdop] = readGPS(gpsPath);
    N_gps = length(gps_hour);
    time_gps = datenum([repmat(gpsYear,N_gps,1),repmat(gpsMon,N_gps,1),repmat(gpsDay,N_gps,1),gps_hour,gps_minute,gps_second]);
    % add a day if necessary
    dTime = diff(time_gps);
    indDayChng = find(dTime<-86398/86400)+1;
    for kk = 1:length(indDayChng)
        time_gps(indDayChng(kk):end)=time_gps(indDayChng(kk):end)+1;
    end
    valid = quality==1;% (lat~=0 & lon~=0);
    gps(jj).file=gpsPath;
    gps(jj).time=time_gps(valid);
    gps(jj).velocity=velocity(valid);
    gps(jj).latitude=lat(valid);
    gps(jj).longitude=lon(valid);
    gps(jj).height = height(valid);
    gps(jj).quality=quality(valid);
    gps(jj).geoidsep=geoidsep(valid);
    gps(jj).pdop=pdop(valid);
    gps(jj).hdop=hdop(valid);
    gps(jj).vdop=vdop(valid);
end
