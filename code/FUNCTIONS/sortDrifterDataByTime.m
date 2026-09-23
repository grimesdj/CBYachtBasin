function time_chunks = sortDrifterDataByTime(trimmed)

% Flatten the structure temporarily
T = trimmed(:);

% Find struct elements with an empty ID
empty_idx = arrayfun(@(x) isempty(x.ID), T);

% Remove those entire elements
T(empty_idx) = [];

% Get all start times
start_times = [T.Start_time];

% Find overall time range
t0 = min(start_times);
t1 = max(start_times) + hours(1.5);

% Define 1.5-hour bins
binWidth = hours(1.5);

% Create bin edges
edges = t0:binWidth:(t1 + binWidth);

% Assign each drifter record to a bin based on Start_time
bin_idx = discretize(start_times, edges);

% Number of 1.5-hour periods
nBins = length(edges) - 1;

% Create output structure
time_chunks = struct([]);

for i = 1:nBins

    % Find records starting within this 1.5-hour interval
    idx = find(bin_idx == i);

    time_chunks(i).Start_time = edges(i);
    time_chunks(i).End_time   = edges(i+1);

    % Store corresponding drifter structures
    time_chunks(i).Drifters = T(idx);

end

end

% How do i divide up the domain -> define polygons -> is it in polygon
% find all drifters and all times when in each region/polygon 
%
% Plot all (all time) the drifter trajectories on one plot; what portion of the basin
% have we sampled?