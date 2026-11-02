function [RSK] = RSKderiveAPT(RSK, varargin)

% RSKderiveAPT - convert paroscientific accelerations period
% and temperature period into accelerations and temperature
%
% Syntax:  [RSK] = RSKderiveAPT(RSK, 'AlignmentCoefficients',Alignment_coefficients,
% 'TemperatureCoefficients',Temperature_coefficients,'AccelerationCoefficients',Acceleration_coefficients)
% 
% The RBRquartz³ APT is a combined triaxial quartz accelerometer and a bottom pressure recorder. 
% It is equipped with a Paroscientific, Inc. triaxial accelerometer and records the acceleration 
% output periods from the accelerometer. RSK files contain only periods for acceleration and temperature.
% 
% This function derives 3-axis accelerations from the accelerometer frequency channels for RSK files. 
% It implements the calibration equations developed by Paroscientific, Inc. to derive accelerations. 
% It requires users to input the alignment and acceleration coefficients. The coefficients are available 
% on the Paroscientific, Inc. triax accelerometer instrument configuration sheet, which is shipped along 
% with the logger. Derived accelerations and temperature are added to the RSK data structure. The channel
% list is updated.
%
%
% Inputs: 
%    [Required] - RSK - Structure containing the logger metadata and data
%
%                 AlignmentCoefficients - a matrix of alignment coefficients on the
%                 Paroscientific, Inc. triax accelerometer instrument configuration sheet
%
%                 TemperatureCoefficients - a matrix of temperature coefficients on the
%                 Paroscientific, Inc. triax accelerometer instrument configuration sheet
%
%                 AccelerationCoefficients - a matrix of acceleration coefficients on the
%                 Paroscientific, Inc. triax accelerometer instrument configuration sheet
%
%
% Outputs:
%    RSK - Updated structure containing the derived APT accelerations and
%    temperature channels
%
% See also: RSKderiveBPR, RSKderiveseapressure, RSKderivedepth.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2022-09-07


p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'AlignmentCoefficients',[], @isnumeric);
addParameter(p, 'TemperatureCoefficients',[], @isnumeric);
addParameter(p, 'AccelerationCoefficients',[], @isnumeric);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
AlignmentCoefficients = p.Results.AlignmentCoefficients;
TemperatureCoefficients = p.Results.TemperatureCoefficients;
AccelerationCoefficients = p.Results.AccelerationCoefficients;

checkDataField(RSK)

if isempty(AlignmentCoefficients) || isempty(TemperatureCoefficients) || isempty(AccelerationCoefficients)
    RSKwarning('Please input the required coefficients from the Paroscientific, Inc. triax accelerometer instrument configuration sheet')
    return
end

% Ensure coefficients from the Paro configuration sheet are input correctly
if all((size(AlignmentCoefficients) == [3,3]) .* (size(TemperatureCoefficients) == [4,3]) .* (size(AccelerationCoefficients) == [10,3]))
else
    RSKwarning('Input coefficient dimensions wrong. Please review the required coefficients from the Paroscientific, Inc. triax accelerometer instrument configuration sheet')
    return 
end

% Find APT by looking for "SACC" in the partNumber or 7-digit partNumber list for APT
isAPT0 = all(strfind(RSK.instruments.partNumber,'SACC'));
isAPT1 = any(strcmp(RSK.instruments.partNumber,{'0004274','0006501','0010311'}));
isAPT = any([isAPT0,isAPT1],'all');

% Find the 4 channel indices for APT acceleration and temperature periods (3 * peri00 + 1 * peri01)
if isAPT
    if length(RSK.channels) < 4
        RSKerror('RBRquartzAPT should contain at least 4 channels. Please apply RSK.deriveAPT to RSK files with triax accelerometer measurements.')
    else    
        for i = 1: length(RSK.channels)
            if i > length(RSK.channels)-3
                RSKerror('No 3-axis acceleration period measurements detected.');
            end
            if all([strcmp(RSK.channels(i).shortName,'peri00') , strcmp(RSK.channels(i+1).shortName,'peri00') , strcmp(RSK.channels(i+2).shortName,'peri00') , strcmp(RSK.channels(i+3).shortName,'peri01')])
                APTchannelx = i;
                APTchannely = i + 1;
                APTchannelz = i + 2;
                APTchanneltemp = i + 3;
                break
            end
        end
    end
else
    RSKerror('No RBRquartzAPT detected. Please apply RSK.deriveAPT to RSK files with APT measurements.');
end


[ accX,tempX ] = deriveParos( RSK.data.values(:,APTchannelx)*1e-12,RSK.data.values(:,APTchanneltemp)*1e-12,TemperatureCoefficients(:,1) , AccelerationCoefficients(:,1));
[ accY,tempY ] = deriveParos( RSK.data.values(:,APTchannely)*1e-12,RSK.data.values(:,APTchanneltemp)*1e-12,TemperatureCoefficients(:,2) , AccelerationCoefficients(:,2));
[ accZ,tempZ ] = deriveParos( RSK.data.values(:,APTchannelz)*1e-12,RSK.data.values(:,APTchanneltemp)*1e-12,TemperatureCoefficients(:,3) , AccelerationCoefficients(:,3));
[ accX,accY,accZ ] = alignParos( accX,accY,accZ, AlignmentCoefficients );

RSK = addchannelmetadata(RSK, 'shortName','accx00', 'longName','X axis acceleration', 'units','m/s²', 'unitsPlainText','m.s-2');
RSK = addchannelmetadata(RSK, 'shortName','accy00', 'longName','Y axis acceleration', 'units','m/s²', 'unitsPlainText','m.s-2');
RSK = addchannelmetadata(RSK, 'shortName','accz00', 'longName','Z axis acceleration', 'units','m/s²', 'unitsPlainText','m.s-2');
RSK = addchannelmetadata(RSK, 'shortName','temp41', 'longName','Accelerometer temperature', 'units','°C', 'unitsPlainText','Degrees_C');
[accXcol,accYcol,accZcol, accTcol] = getchannelindex(RSK, {'X axis acceleration','Y axis acceleration','Z axis acceleration','Accelerometer temperature'});

RSK.data.values(:,accXcol) = accX;
RSK.data.values(:,accYcol) = accY;
RSK.data.values(:,accZcol) = accZ;
RSK.data.values(:,accTcol) = (tempX + tempY + tempZ)/3;

logentry = ('APT temperature and accelerations derived from period data.');
RSK = RSKappendtolog(RSK, logentry);


    %% Nested Functions
    function [ acceleration,temperature ] = deriveParos( accPeriods,tempPeriods,tempCoefficients,accCoefficients )
    % Returns the acceleration for each axis and the temperature
    % periods in seconds
    % output acceleration in m/s^-2
    % output temperature in °C

    u0 = tempCoefficients(1);
    y1 = tempCoefficients(2);
    y2 = tempCoefficients(3);
    y3 = tempCoefficients(4);

    c1 = accCoefficients(1);
    c2 = accCoefficients(2);
    c3 = accCoefficients(3);
    d1 = accCoefficients(4);
    d2 = accCoefficients(5);
    t1 = accCoefficients(6);
    t2 = accCoefficients(7);
    t3 = accCoefficients(8);
    t4 = accCoefficients(9);
    t5 = accCoefficients(10);

    U = (tempPeriods/(1e-6)) - u0;
    temperature = y1 .* U + y2 .* U .* U + y3 .* U .* U .* U;
    Tsquare = (accPeriods/(1e-6)) .* (accPeriods/(1e-6));

    C = c1 + c2 .* U + c3 .* U .^ 2;
    D = d1 + d2 .* U;
    T0 = t1 + t2 .* U + t3 .* U .^ 2 + t4 .* U .^ 3 + t5 .* U .^ 4;

    acceleration = C .* (1 - T0 .^ 2 ./ Tsquare) .* (1 - D .* (1 - T0 .^ 2 ./ Tsquare));
    end


    function [ gx,gy,gz ] = alignParos( accX,accY,accZ, alignCoefficients)
    % Aligns x, y, z to standard coordinate
    % input acceleration in m/s^-2
    % output acceleration in m/s^-2
        gx=alignCoefficients(1,1).*accX+alignCoefficients(1,2).*accY+alignCoefficients(1,3).*accZ;
        gy=alignCoefficients(2,1).*accX+alignCoefficients(2,2).*accY+alignCoefficients(2,3).*accZ;
        gz=alignCoefficients(3,1).*accX+alignCoefficients(3,2).*accY+alignCoefficients(3,3).*accZ;
    end

end
