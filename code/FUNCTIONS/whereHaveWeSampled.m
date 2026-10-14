clear all;
close all;
%%
main_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'; % where is the CB_YachtBasin shared DATA folder on your machine

% Where to save the 'data' struct
out_dir = ([main_dir 'CB_YachtBasin\ALL_IOP_PROCESSING\DRIFTERS\']);

%% Get alltime drifter data (from all IOPs) and store them in a structure

march = load([main_dir,'CB_YachtBasin\MarshMadness\PROCESSED_DATA\DRIFTERS\MarshMadness_drifters_level_1.mat']);

data(1).name = 'MarchMadness';
data(1).data = march.drifters_level_1
data(1).color = 'k';

may = load([main_dir,'CB_YachtBasin\MarshMayhem\PROCESSED_DATA\DRIFTERS\MarshMayhem_drifters_level_1.mat']);

data(2).name = 'MarchMayhem';
data(2).data = may.drifters_level_1;
data(2).color = 'r';

oct_1008 = load([main_dir,'CB_YachtBasin\FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_100825_drifters_level_1.mat']);

data(3).name = 'FallFrolic100825';
data(3).data = oct_1008.drifters_level_1;
data(3).color = 'b';

oct_1010 = load([main_dir,'CB_YachtBasin\FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_101025_drifters_level_1.mat']);

data(4).name = 'FallFrolic101025';
data(4).data = oct_1010.drifters_level_1;
data(4).color = 'y';

nov = load([main_dir,'CB_YachtBasin\NovDep\PROCESSED_DATA\DRIFTERS\NovDep_drifters_level_1.mat']);

data(5).name = 'NovDep';
data(5).data = nov.drifters_level_1;
data(5).color = 'm';

june = load([main_dir,'CB_YachtBasin\JuneJamboree\PROCESSED_DATA\DRIFTERS\JuneJamboree_drifters_level_1.mat']);

data(6).name = 'JuneJamboree';
data(6).data = june.drifters_level_1;
data(6).color = 'g';

%% Save the 'data' struct
save(fullfile(out_dir, ['all_drifter_data.mat']), 'data');


%% Time to plot it all

figure;
gx = geoaxes;
geobasemap = 'satellite';
hold(gx,'on')

legendHandles = gobjects(numel(data),1);

for i = 1:numel(data)
    
  dep_name = data(i).name;
   
  drifters = data(i).data;
  
  color = data(i).color;
  
  % Create one proxy line for this deployment's legend entry
  legendHandles(i) = geoplot(gx, NaN, NaN, ...
        'Color',color, ...
        'LineWidth',1.5, ...
        'DisplayName',data(i).name);
  

  
  for j = 1: numel(drifters)
      

      h = geoplot(gx, drifters(j).Lat,drifters(j).Lon,'Color',color,'LineWidth',1.5);
      
% Store for legend

 

      
  end
    
    
    
    
end

legend(gx, legendHandles, {data.name},'Location','best');