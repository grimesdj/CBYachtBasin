function RSK = addchannelmetadata(RSK, varargin)

% ADDCHANNELMETADATA - Add the metadata for a new channel.
%
% Syntax:  [RSK] = ADDCHANNELMETADATA(RSK, 'shortName', 'longName', 'units', [OPTIONS])
% 
% Adds all the metadata associated with a new channel in the channel field
% of the RSK structure.
%
% Inputs:
%   [Required] - RSK - Input RSK structure
%
%              - shortName - Short name of the new channel
%
%              - longName - Full name of the new channel
%            
%              - units - Units of the new channel
% 
%   [Optional] - unitsPlainText - Units of the new channel in plain text. 
%
% Outputs:
%    RSK - RSK structure containing new channel metadata.
%
% See also: RSKderivedepth, RSKderiveseapressure, RSKderivesalinity.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2018-11-09

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p,'shortName','', @ischar);
addParameter(p,'longName','', @ischar);
addParameter(p,'units','', @ischar);
addParameter(p,'unitsPlainText','',@ischar);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
shortName = p.Results.shortName;
longName = p.Results.longName;
units = p.Results.units;
unitsPlainText = p.Results.unitsPlainText;

if isempty(shortName) || isempty(longName) || isempty(units)
    RSKwarning('Please specify channel shortName, channel longName, and units.')
    return
end

hasChan = any(strcmpi({RSK.channels.longName}, longName));

if ~hasChan
    nchannels = length(RSK.channels);
    RSK.channels(nchannels+1).shortName = shortName;
    RSK.channels(nchannels+1).longName = longName;
    RSK.channels(nchannels+1).units = units;
    if ~isempty(unitsPlainText)
        RSK.channels(nchannels+1).unitsPlainText = unitsPlainText;
    end
    if ~isempty(RSK.channels(nchannels).channelID)
        RSK.channels(nchannels+1).channelID = max(nchannels+1,RSK.channels(nchannels).channelID+1);
    end
end

end