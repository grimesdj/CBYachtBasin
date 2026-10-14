function [drifters_QAQC] = veloQAQC(cfg)



proj = projcrs(32119);

% load the drifter data
filename = fullfile(cfg.out.drifters_data, [cfg.name '_raw_drifters.mat']);
load(filename);

% Get the deployment times
dep_fin = readtable(cfg.drifters_dep_times);

for i = 1:numel(raw_drifters)

   idx = raw_drifters(i).ID;
   idy = find(idx == table2array(dep_fin(:,1)));

   tstarts = table2array(dep_fin(idy,3));
   tstarts = datetime(tstarts,'ConvertFrom','excel');
   tstarts = timeofday(tstarts);
   tends = table2array(dep_fin(idy,5));
   tends = datetime(tends,'ConvertFrom','excel');
   tends = timeofday(tends);
   dates_end = table2array(dep_fin(idy,6));
   dates_start = table2array(dep_fin(idy,4));

   tstarts = tstarts + dates_start;
   tends = tends + dates_end;





    
   figure;
   plot(raw_drifters(i).Date_and_time,raw_drifters(i).Speed)
%    ylim([0 15]);
   xlim([min(tstarts) - hours(1), max(tends) + hours(1)])
   title(raw_drifters(i).ID)
   hold on;
   yline(2,'k--','LineWidth',2);
   xline(tstarts, 'r', 'LineWidth', 2);
   xline(tends, 'k', 'LineWidth', 2);
   hold off;



    
    
end

drifters_QAQC = raw_drifters;





end