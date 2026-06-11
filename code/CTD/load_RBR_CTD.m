% directory where all data is stored:
data_dir = '/Users/derekgrimes/OneDriveUNCW/DATA/CB_YachtBasin/NovDep/';
% RBR-CTD is in "RBR" sub-directory:
ctd_dir = [data_dir,filesep,'CTD',filesep];
% example serial number:
SN = 210866;

% build full path to data file
rbrFileStr = sprintf('%s%04d*.rsk',ctd_dir,SN);
rbrFile = dir(rbrFileStr);
fin     = [rbrFile.folder,filesep,rbrFile.name];


% open the file, then read data
try rsk = RSKopen(fin);
catch
    disp(['missing data file: ', rbrFileStr])
    return
end
rsk = RSKreaddata(rsk);
% parse data from structure:
rbr_time = rsk.data.tstamp;
rbr_data = rsk.data.values;
rbr_cond = rbr_data(:,1);
rbr_temp = rbr_data(:,2);
rbr_pres = rbr_data(:,3);
% need to offset pressure, this is a quick and dirty method
rbr_pres_offset = mean(rbr_pres(rbr_pres<10.5 & rbr_pres>9.5));

% convert conductivity to salinity
rbr_salt = gsw_SP_from_C(rbr_cond,rbr_temp,rbr_pres-rbr_pres_offset);


figure, subplot(3,1,1) , plot(datetime(rbr_time,'convertFrom','datenum'), rbr_pres-rbr_pres_offset),

subplot(3,1,2), plot(datetime(rbr_time,'convertFrom','datenum'), rbr_temp),
subplot(3,1,3), plot(datetime(rbr_time,'convertFrom','datenum'), rbr_salt),