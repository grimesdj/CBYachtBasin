lat_lon = [34.050288,-77.88789331];
lat_lon = mean(lat_lon,1)
HHMM = abs(rem(lat_lon,1));
MM = HHMM*60
SS = rem(HHMM*60,1)*60