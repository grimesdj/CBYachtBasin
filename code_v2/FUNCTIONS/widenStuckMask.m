function Sd = widenStuckMask(St, nn, half)
% Stuck flag -> mask of samples to drop, widened by half samples on each side.
%   St   = R(k).Stuck for one track (may be empty or the wrong length)
%   nn   = number of samples the track has in the analysis
%   half = margin in seconds each side (dilate_s/2 = 128)
St = St(:);
if isempty(St), St = false(nn,1); end
if numel(St) < nn, St(end+1:nn,1) = false; end
St = St(1:nn);
if half > 0 && any(St)
    Sd = movmax(double(St), 2*half+1, 'Endpoints','shrink') > 0;
else
    Sd = St;
end
end