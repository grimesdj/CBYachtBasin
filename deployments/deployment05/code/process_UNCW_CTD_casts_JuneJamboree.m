clear all
close all
%% Where is the data?
rootDir = '~/OneDriveUNCW/DATA/CB_YachtBasin/JuneJamboree/'
L0Dir  = [rootDir,'PROCESSED_DATA',filesep,'CTD',filesep];
%% RBR-CTD+F:
ctdSerialNumber = 230442;
rbrFile = dir(sprintf('%sRBR/%d_*.rsk',rootDir,ctdSerialNumber));
%% Locsys GT-31 GPS:
gpsFile = [rootDir,'GPS/GPSUSER05_133200472_20260614_110619.TXT'];

%% Output file:
fout = [L0Dir,filesep,num2str(ctdSerialNumber),'_casts_L0.mat'];



%% get the GPS data
gps_dir = dir(gpsFile);
gpsStr = split(gps_dir.name,{'_','.'});
gpsDate= gpsStr{3}();
gpsYear= str2num(gpsDate(1:4)); gpsMon = str2num(gpsDate(5:6)); gpsDay = str2num(gpsDate(7:8));
gpsTime= gpsStr{4}();
gpsHH  = str2num(gpsTime(1:2)); gpsMM = str2num(gpsTime(3:4)); gpsSS = str2num(gpsTime(5:6));
[gps_hour, gps_minute, gps_second, velocity, numsats, height, lat, lon, quality, geoidsep, pdop, hdop, vdop] = readGPS(gpsFile);
N_gps = length(gps_hour);
time_gps = datenum([repmat(gpsYear,N_gps,1),repmat(gpsMon,N_gps,1),repmat(gpsDay,N_gps,1),gps_hour,gps_minute,gps_second]);
% add a day if necessary
dTime = diff(time_gps);
indDayChng = find(dTime<-86398/86400)+1;
for kk = 1:length(indDayChng)
    time_gps(indDayChng(kk):end)=time_gps(indDayChng(kk):end)+1;
end
valid = quality==1;% (lat~=0 & lon~=0);
gps.file=gpsFile;
gps.time_utc=time_gps(valid);
gps.velocity=velocity(valid);
gps.latitude=lat(valid);
gps.longitude=lon(valid);
gps.height = height(valid);
gps.quality=quality(valid);
gps.geoidsep=geoidsep(valid);
gps.pdop=pdop(valid);
gps.hdop=hdop(valid);
gps.vdop=vdop(valid);
%% convert to EDT
gps.time = gps.time_utc-4/24;


%% open the RBR file
rsk = RSKopen([rbrFile.folder,filesep,rbrFile.name]);
%
rsk = RSKreaddata(rsk);
serialNum = rsk.instruments.serialID;
time = rsk.data.tstamp;
data = rsk.data.values;
%
%
cond     = data(:,1);
temp     = data(:,2);
pres_abs = data(:,3);
chla_raw  = data(:,4);
fdom_raw  = data(:,5);
turb_raw  = data(:,6);
salt      = data(:,9);

%% Atmospheric Pressure offset
ATM_Time = [datenum('14-June-2026 12:50:00'), datenum('14-June-2026 13:00:00')];
inAir  = find(time>=ATM_Time(1) & time<=ATM_Time(2));
pres0 = mean(pres_abs(inAir));
pres  = pres_abs-pres0;

%% interpolate downcasts to common depth grid:
% break record into individual casts
% 1) first smooth pressure record to find zero-up crossings
tau = 3;% smooth timescale in seconds
dt  = 1/16;% sample rate Hz
N   = tau/dt;
f   = hanning(N); f = f./sum(f);
pres_avg = conv(pres,f,'same');
dpres_avg = gradient(pres_avg);

%% locate the 'casts'
cast          = pres_avg>=0.15;
downcast      = find(cast & dpres_avg>0);
upcast        = find(cast & dpres_avg<0);
cast          = find(cast);
%
[startINDs,endINDs] = Segment(cast,(1/dt));
cast_start     = cast(startINDs);
cast_end       = cast(endINDs);
valid          = (cast_end-cast_start)>10;
cast_start(~valid)=[];
cast_end(~valid)  =[];
Ncasts         = length(cast_start);

%% Create a uniform pressure grid with 10cm resolution to the maximum depth (6 m)
pres_grid = [0:0.1:6]';
Npres     = length(pres_grid);
time_grid = nan(1,Ncasts);
latitude_grid  = nan(1,Ncasts);
longitude_grid = nan(1,Ncasts);
temp_grid = nan(Npres,Ncasts);
fdom_grid  = nan(Npres,Ncasts);
turb_grid  = nan(Npres,Ncasts);
chla_grid  = nan(Npres,Ncasts);
cond_grid = nan(Npres,Ncasts);
salt_grid  = nan(Npres,Ncasts);
for jj = 1:Ncasts
    this_cast = find(downcast>=cast_start(jj) & downcast<=cast_end(jj));
    [~,inds]  = unique(pres(downcast(this_cast)));
    time_grid(1,jj) = time(downcast(this_cast(inds(1))));
    try latitude_grid(1,jj)  = interp1(gps(1).time,gps(1).latitude,time_grid(1,jj));
    catch
        latitude_grid(1,jj)  = nan;
    end
    try longitude_grid(1,jj) = interp1(gps(1).time,gps(1).longitude,time_grid(1,jj));
    catch
        longitude_grid(1,jj) = nan;
    end
    temp_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),temp(downcast(this_cast(inds))),pres_grid);
    try cond_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),cond(downcast(this_cast(inds))),pres_grid);
    catch
        cond_grid(:,jj) = nan;
    end
    try salt_grid(:,jj) = interp1(pres(downcast(this_cast(inds))),salt(downcast(this_cast(inds))),pres_grid);
    catch
        salt_grid(:,jj) = nan;
    end
    try fdom_grid(:,jj)  = interp1(pres(downcast(this_cast(inds))),fdom_raw(downcast(this_cast(inds))),pres_grid);
    catch
        fdom_grid(:,jj)  = nan;
    end
    try turb_grid(:,jj)  = interp1(pres(downcast(this_cast(inds))),turb_raw(downcast(this_cast(inds))),pres_grid);
    catch
        turb_grid(:,jj)  = nan;
    end
    try chla_grid(:,jj)  = interp1(pres(downcast(this_cast(inds))),chla_raw(downcast(this_cast(inds))),pres_grid);
    catch
        chla_grid(:,jj)  = nan;
    end
end
%
% save gridded data to netcdf files
CTD = struct('Time',time_grid,'Latitude',latitude_grid,'Longitude',longitude_grid,'Pressure',pres_grid,'CHLa',chla_grid,'FDOM',fdom_grid,'NTU',turb_grid,'Temperature',temp_grid,'Conductivity',cond_grid,'Salinity',salt_grid)
%
% archive raw and gridded data
save(fout,'time','cond','temp','pres','salt','chla_raw','fdom_raw','turb_raw','gps','time_grid','pres_grid','temp_grid','cond_grid','salt_grid','chla_grid','turb_grid','fdom_grid','latitude_grid','longitude_grid','pres0')
