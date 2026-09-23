% directory where all data is stored:
% % data_dir = '/Users/derekgrimes/OneDriveUNCW/DATA/CB_YachtBasin/MarshMadness/';
data_dir = '/Users/derekgrimes/OneDriveUNCW/DATA/CB_YachtBasin/SummerFlood/';
% RBR-CTD is in "RBR" sub-directory:
ctd_dir = [data_dir,filesep,'RBR',filesep];
% example serial number:
% SN = 210864;
SN = [210864,210866,242452]; %210864;

% need to offset pressure, this is a quick and dirty method
atmTime = [datenum('31-Jul-2026 16:00:00') datenum('31-Jul-2026 17:37:00')];

fig0 = figure;
for jj=1:length(SNs)
    SN = SN(jj);

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

idx = find(rbr_time>=atmTime(1) & rbr_time<=atmTime(2));
rbr_pres_offset = median(rbr_pres(idx),'omitnan');

% convert conductivity to salinity
rbr_salt = gsw_SP_from_C(rbr_cond,rbr_temp,rbr_pres-rbr_pres_offset);


% $$$ figure, subplot(3,1,1) , plot(datetime(rbr_time,'convertFrom','datenum'), rbr_pres-rbr_pres_offset), ylim([-1 2])
% $$$ 
% $$$ subplot(3,1,2), plot(datetime(rbr_time,'convertFrom','datenum'), rbr_temp), ylim([25 30])
% $$$ subplot(3,1,3), plot(datetime(rbr_time,'convertFrom','datenum'), rbr_salt), ylim([0 20])

% $$$ 
hold on,
plot(datetime(rbr_time,'convertFrom','datenum'),rbr_salt,'-');

end
ylabel('Salinity [psu]','interpreter','latex')
set(gca,'ticklabelinterpreter','latex','tickdir','out','ylim',[0 36], 'xlim', [datetime('31-Jul-2026 18:00:00') datetime('02-Sep-2026 09:00:00')])

legend({'South-Bottom','North-Bottom','North-Surface'},'interpreter','latex','location','southeast')
exportgraphics(gcf,[data_dir,'salinity_vs_time.pdf'])