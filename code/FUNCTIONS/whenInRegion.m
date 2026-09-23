clear all;
close all;
%%
proj = projcrs(32119);

main_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'; % where is the CB_YachtBasin shared DATA folder on your machine

%% Load the Drifter data
load(fullfile([main_dir 'CB_YachtBasin\ALL_IOP_PROCESSING\DRIFTERS\all_drifter_data.mat']));

% %% Create the regions
% regions = struct([]);
% 
% % Shared latitude boundaries
% latNorth = 34.0577;
% latMid1  = 34.0502;
% latMid2  = 34.0444;
% latSouth = 34.0355;
% 
% % Region 1: northwest
% regions(1).Name = 'Region 1';
% 
% regions(1).Lat = [
%     latNorth
%     latNorth
%     latMid1
%     latMid1
%     latNorth
% ];
% 
% regions(1).Lon = [
%    -77.8934
%    -77.8875
%    -77.8875
%    -77.8934
%    -77.8934
% ];
% 
% % Region 2: northeast
% regions(2).Name = 'Region 2';
% 
% regions(2).Lat = [
%     latNorth
%     latNorth
%     latMid1
%     latMid1
%     latNorth
% ];
% 
% regions(2).Lon = [
%    -77.8899
%    -77.8845
%    -77.8845
%    -77.8899
%    -77.8899
% ];
% 
% % Shared boundary between Regions 3 and 4
% lonSplit34 = -77.88945;
% 
% % Region 3: middle-west
% regions(3).Name = 'Region 3';
% 
% regions(3).Lat = [
%     latMid1
%     latMid1
%     latMid2
%     latMid2
%     latMid1
% ];
% 
% regions(3).Lon = [
%    -77.8945
%     lonSplit34
%     lonSplit34
%    -77.8945
%    -77.8945
% ];
% 
% % Region 4: middle-east
% regions(4).Name = 'Region 4';
% 
% regions(4).Lat = [
%     latMid1
%     latMid1
%     latMid2
%     latMid2
%     latMid1
% ];
% 
% regions(4).Lon = [
%     lonSplit34
%    -77.8865
%    -77.8865
%     lonSplit34
%     lonSplit34
% ];
% 
% % Region 5: southern basin
% regions(5).Name = 'Region 5';
% 
% regions(5).Lat = [
%     latMid2
%     latMid2
%     latSouth
%     latSouth
%     latMid2
% ];
% 
% regions(5).Lon = [
%    -77.8937
%    -77.8883
%    -77.8883
%    -77.8937
%    -77.8937
% ];

%% Define 7 geographic regions

% -------------------------
% Region 1
% -------------------------
regions(1).Lon = [ ...
    -77.8912186 ... % R1_1
    -77.8914118 ... % R1_2
    -77.8878713 ... % R1_3
    -77.8889120 ... % R1_4
    -77.8912186];   % close polygon

regions(1).Lat = [ ...
     34.0571503 ...
     34.0537281 ...
     34.0540570 ...
     34.0576570 ...
     34.0571503];


% -------------------------
% Region 2
% -------------------------
regions(2).Lon = [ ...
    -77.8914118 ... % R2_1
    -77.8912502 ... % R2_2
    -77.8888898 ... % R2_3
    -77.8878713 ... % R2_4
    -77.8914118];

regions(2).Lat = [ ...
     34.0537281 ...
     34.0508411 ...
     34.0515167 ...
     34.0540570 ...
     34.0537281];


% -------------------------
% Region 3
% -------------------------
regions(3).Lon = [ ...
    -77.8878713 ... % R3_1
    -77.8888898 ... % R3_2
    -77.8872221 ... % R3_3
    -77.8860741 ... % R3_4
    -77.8851836 ... % R3_5
    -77.8878713];

regions(3).Lat = [ ...
     34.0540570 ...
     34.0515167 ...
     34.0512945 ...
     34.0508627 ...
     34.0531739 ...
     34.0540570];


% -------------------------
% Region 4
% -------------------------
regions(4).Lon = [ ...
    -77.8912502 ... % R4_1
    -77.8912561 ... % R4_2
    -77.8897005 ... % R4_3
    -77.8888898 ... % R4_4
    -77.8912502];

regions(4).Lat = [ ...
     34.0508411 ...
     34.0461957 ...
     34.0461690 ...
     34.0515167 ...
     34.0508411];


% -------------------------
% Region 5
% -------------------------
regions(5).Lon = [ ...
    -77.8888898 ... % R5_1
    -77.8897005 ... % R5_2
    -77.8872221 ... % R5_3
    -77.8888898];

regions(5).Lat = [ ...
     34.0515167 ...
     34.0461690 ...
     34.0512945 ...
     34.0515167];


% -------------------------
% Region 6
% -------------------------
regions(6).Lon = [ ...
    -77.8872221 ... % R6_1
    -77.8897005 ... % R6_2
    -77.8884393 ... % R6_3
    -77.8860741 ... % R6_4
    -77.8872221];

regions(6).Lat = [ ...
     34.0512945 ...
     34.0461690 ...
     34.0459614 ...
     34.0508627 ...
     34.0512945];


% -------------------------
% Region 7
% -------------------------
regions(7).Lon = [ ...
    -77.8912561 ... % R7_1
    -77.8936429 ... % R7_2
    -77.8917117 ... % R7_3
    -77.8884393 ... % R7_4
    -77.8897005 ... % R7_5
    -77.8912561];

regions(7).Lat = [ ...
     34.0461957 ...
     34.0369647 ...
     34.0365913 ...
     34.0459614 ...
     34.0461690 ...
     34.0461957];
%% Plot the regions
reg = figure;
gx = geoaxes;
hold(gx,'on')

geobasemap(gx, 'satellite')

colors = lines(numel(regions));

for i = 1:numel(regions)
    
   geoplot(regions(i).Lat,regions(i).Lon,'Color',colors(i,:),'LineWidth',2) 
    
end

%% Convert the regions from lat lon to ENU
for k = 1:numel(regions)

    [regions(k).Easting, regions(k).Northing] = projfwd(proj, regions(k).Lat, regions(k).Lon);

end

%% Asign region ID to each drifter observation


%% Make figure for a specific region for ONE OBSERVATION PERIOD



