% directory where all data is stored:
data_dir = '/Users/derekgrimes/OneDriveUNCW/DATA/CB_YachtBasin/MarshMadness/';
% RBR-CTD is in "RBR" sub-directory:
ctd_dir = [data_dir,filesep,'RBR',filesep];
% example serial number:
SN = 210864;
rbrFileStr = sprintf('%s%04d*.rsk',ctd_dir,SN);
rbrFile = dir(rbrFileStr);
fin     = [rbrFile.folder,filesep,rbrFile.name];
% open the file, then read data
try rsk = RSKopen(fin);
catch
    disp(['missing data file: ', rbrFileStr])
    continue
end
rsk = RSKreaddata(rsk);
rbr_time = rsk.data.tstamp;
rbr_data = rsk.data.values;
rbr_cond = rbr_data(:,1);
rbr_temp = rbr_data(:,2);
rbr_pres = rbr_data(:,3);
% need to offset pressure,
rbr_pres_offset = mean(rbr_pres(rbr_pres<10.5 & rbr_pres>9.5));
rbr_salt = gsw_SP_from_C(rbr_cond,rbr_temp,rbr_pres-rbr_pres_offset);
% put info into structure array
