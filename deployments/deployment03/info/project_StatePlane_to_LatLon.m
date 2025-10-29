clear all
close all
%
proj = projcrs(32119);%(2264)
%
gps_file = '~/Grants/CBFLOOD/deployments/2025_10/NSF-CB.csv';
format   = '%s %f %f %f';
fid      = fopen(gps_file);
data     = textscan(fid,format,'delimiter',',');
%
Name     = data{1};
Easting  = data{3};
Northing = data{2};
Elevation = data{4};
% $$$ 
% $$$ 
% $$$ Northing = [3.348013000000000e+04, 3.460218233333333e+04];
% $$$ Easting  = [7.122062629999999e+05, 7.121378996666666e+05];

[Latitude,Longitude] = projinv(proj,Easting,Northing);

figure, geoplot(Latitude,Longitude,'*r')
geobasemap('satellite')

header = {'Site ID', 'Northing', 'Easting', 'Elevation', 'Latitude','Longitude'};
N      = length(Northing);
data   = cat(2,Name,mat2cell(Northing,ones(N,1)), mat2cell(Easting,ones(N,1)),mat2cell(Elevation,ones(N,1)), ...
                    mat2cell(Latitude,ones(N,1)), mat2cell(Longitude,ones(N,1)));
out = cat(1,header, data);

file_out = '/Users/derekgrimes/git/cbflood/deployments/deployment03/info/GPS_points_20251003.csv';
writecell(out,file_out)