function [lat_new, lon_new] = convertToDecimalDegrees(lats_raw, lons_raw)
    % Converts latitude and longitude from DDMM.MMMM format to decimal degrees
    
    % Convert latitude
    lat_dd = floor(lats_raw / 100);      % Extract degrees
    lat_mm = lats_raw - (lat_dd * 100);  % Extract minutes
    lat_new = lat_dd + (lat_mm / 60);    % Convert to decimal degrees

    % Convert longitude
    lon_dd = floor(lons_raw / 100);      % Extract degrees
    lon_mm = lons_raw - (lon_dd * 100);  % Extract minutes
    lon_new = lon_dd + (lon_mm / 60);    % Convert to decimal degrees
    lon_new = -lon_new;
end