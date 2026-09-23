function [hour, minute, second, velocity, numsats, height, lat, long, quality, geoidsep, pdop, hdop, vdop] = readGPS(Filename)
%% Ließt NMEA GPS Daten ein
% GPSdaten=readGPS(filename)
% GSPdaten ist eine Matrix mit folgendem Inhalt:
% 
% Definition: http://de.wikipedia.org/wiki/NMEA_0183

fid = fopen(Filename,'rt');

%nach den Zeilenvorschüben suchen, 10 ist der ASCII-Code für einen
%Zeilenvorschub
dateilange = nnz(fread(fid)==10);
disp(sprintf('Approximately %d lines are read...',dateilange));

%Pointer wieder auf Beginn setzen
fseek(fid,0,'bof');

%Preallocation der Daten
zeile = 1;
GPS_data=zeros(dateilange,13);
hour = 0;
minute = 0;
second = 0;
UTC = 0;
lat_degree = 0;
lat_decimal = 0;
long_decimal = 0;
lat = 0;
long = 0;
lat_A = '0';
long_A = '0';

velocity = 0;
course = 0;
numsats = 0;
height=0;
lat_decimal=0;
long_decimal=0;
quality = 0;
geoidsep = 0;

pdop = 0;
hdop = 0;
vdop = 0;

count = 0;
while ~feof(fid) % so lange Dateiende nicht erreicht ist
    count = count+1;
    line = fgetl(fid); % Zeile lesen
    if isempty(line) || ~isequal(line(1),'$') % wenn leere Zeile
        continue     % überspringen
    elseif strncmp(line,'$GPGGA',6) % wenn Zeile mit $GPRMC
      data = textscan(line,'%s%f%f%c%f%c%f%f%f%f%c%f',1,'delimiter',','); 
        % compute UTC(HHMMSS.SSS), Universal Time Coordinated 
        hour = floor(data{2}/10000); 
        minute = floor((data{2}-hour*10000)/100); 
        second = round(data{2}-floor(data{2}/100)*100); 
        UTC = strcat(num2str(hour),':',num2str(minute),':',num2str(second)); 

        %compute latitude(DDMM.MMM) and longitude(DDDMM.MMMM) 
        lat_degree = floor(data{3}/100);
        lat_decimal = round((data{3}-lat_degree*100)/60*10^8)/10^8; 
        lat = lat_degree+lat_decimal;
        lat_A = strcat(num2str(lat_degree+lat_decimal),data{4});

        long_degree= floor(data{5}/100); 
        long_decimal = round((data{5}-long_degree*100)/60*10^8)/10^8; 
        long = long_degree+long_decimal;
        long_A = strcat(num2str(long_degree+long_decimal),data{6});

        % switch the sign for N/S or E/W
        if data{4}=='S'
            lat=-lat;
        end
        if data{6}=='W'
            long=-long;
        end

        numsats = data{8};

        %GPS-Qualität:
        % 0 für ungültig
        % 1 für GPS fix
        % 2 für DGPS fix
        % 6 für geschätzt (nur bei NMEA-0183 ab Version 2.3)
        quality = data{7};

        hdop = data{9};
        
        %Höhe der Antenne über Geoid oder MSL (mean sea level)
        height = data{10};
        
        % get the geoid separation used
        geoidsep = data{12};
        
    elseif strncmp(line,'$GPVTG',6) % wenn Zeile mit $GPRMC
        GPVTGdata = textscan(line,'%s%f%c%f%c%f%c%f%c',1,'delimiter',',');
        % compute velocity(km/h) and course 
        velocity = GPVTGdata{8}; 
        course = GPVTGdata{2}; 
        
    elseif strncmp(line,'$GPGSA',6)
      GPGSAdata = textscan(line,'%s%c%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',1,'delimiter',',');  
      pdop = GPGSAdata{16};
      hdop = GPGSAdata{17};
      vdop = GPGSAdata{18};

    end
    tmp = [hour, minute, second, velocity, numsats, height, lat, long, quality, geoidsep, pdop, hdop, vdop];
    if length(tmp)==13
    GPS_data(zeile,:)=tmp;
    zeile = zeile + 1; % eine Zeile weiter
    end
end
fclose(fid); %schließen

hour = GPS_data(:,1);
minute = GPS_data(:,2);
second = GPS_data(:,3);
velocity = GPS_data(:,4);
numsats = GPS_data(:,5);
height = GPS_data(:,6);
lat = GPS_data(:,7);
long = GPS_data(:,8);
quality = GPS_data(:,9);
geoidsep = GPS_data(:,10);
pdop = GPS_data(:,11);
hdop = GPS_data(:,12);
vdop = GPS_data(:,13);
