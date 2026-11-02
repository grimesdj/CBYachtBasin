function [phase, high_idx, low_idx] = assignTidalPhase(time, eta, blank_min, minSep_hr, smooth_min)
% Same logic as findTidalPhaseandRegion, with the sample interval measured
% from the data. time = datenum (days); blank_min = minutes either side of
% each high/low labelled "High"/"Low"; minSep_hr = minimum hours between highs.

sz   = size(eta);
time = time(:);  eta = eta(:);

dt        = median(diff(time));
minDist   = round((minSep_hr/24)/dt);
peakBlank = blank_min/(24*60);

if nargin < 5, smooth_min = 0; end
nS   = max(1, round(smooth_min/(dt*1440)));
nS   = 2*floor(nS/2) + 1;                 % odd window so no time shift
etaS = movmean(eta, nS);
[~,high_idx] = findpeaks( etaS, 'MinPeakDistance', minDist);
[~,low_idx]  = findpeaks(-etaS, 'MinPeakDistance', minDist);

phase = strings(size(eta));

ext_idx  = [high_idx(:); low_idx(:)];
ext_type = [repmat("High",numel(high_idx),1); repmat("Low",numel(low_idx),1)];
[ext_idx,order] = sort(ext_idx);
ext_type = ext_type(order);

if any(ext_type(1:end-1) == ext_type(2:end))
    warning('Two highs or two lows in a row: an extremum was missed or doubled.');
end

for i = 1:numel(ext_idx)-1
    if ext_type(i)=="High" && ext_type(i+1)=="Low"
        phase(ext_idx(i):ext_idx(i+1)) = "Falling";
    elseif ext_type(i)=="Low" && ext_type(i+1)=="High"
        phase(ext_idx(i):ext_idx(i+1)) = "Rising";
    end
end

for i = 1:numel(high_idx)
    phase(abs(time - time(high_idx(i))) <= peakBlank) = "High";
end
for i = 1:numel(low_idx)
    phase(abs(time - time(low_idx(i))) <= peakBlank) = "Low";
end

phase = reshape(phase, sz);
end