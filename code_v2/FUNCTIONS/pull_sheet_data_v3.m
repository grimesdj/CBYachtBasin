%% pull_sheet_data_v3
% Downloads each drifter's tab from the published Google Sheet (same links as
% process_google_sheet_more_drifters_V3) and writes one data_XX_v3.xlsx per drifter,
% in the same 12-column layout as data_XX_v2.xlsx, so drifter_post_processing_sheet
% can read it with cfg.sheet_file_suffix = '_v3'.
%
% Layout written (no header row):
%   date | time as text hh:mm:ss.SSS | counter | lat ddmm.mmmm | N | lon dddmm.mmmm | W |
%   altitude | speed (knots) | angle | volt | antenna
%
% Nothing is skipped silently: every drifter and every problem is printed.

clear; close all

%% Settings
drifter_ids = 1:17;      % drifters to download
out_dir = 'C:\Users\bcm3620\OneDrive - UNC-Wilmington\CB_YachtBasin\Archived_Sheet_Data\RAW';
suffix  = '_v3';         % written file names: data_XX_v3.xlsx
dateFmt = 'M/d/yyyy';    % date format of the date field in the message (as in findLastFix)

%% Sheet information (copied from process_google_sheet_more_drifters_V3)
deviceIDs = {
    'e00fce682357c1951e4cdf1a';   % 01
    'e00fce68fa48f531becc05c6';   % 02
    'e00fce6885671f4c89fdb2f7';   % 03
    'e00fce68b7c9e4800d61f9d2';   % 04
    'e00fce68b333416c10babe5c';   % 05
    'e00fce68a00a5d36679a786b';   % 06
    'e00fce6859db454acd21b18e';   % 07
    'e00fce68d33a4d868b70651a';   % 08
    'e00fce6877578761c67a3f72';   % 09
    'e00fce68d25dd33847451b43';   % 10
    'e00fce68d9662407c784127e';   % 11
    'e00fce682cdba2c74331e417';   % 12
    'e00fce68539c5e2c221cdf73';   % 13
    'e00fce68985dad7ee73f458e';   % 14
    'e00fce683fd35d2e0db6c0cf';   % 15
    'e00fce68d03baacd17a3a4ec';   % 16
    'e00fce683f96ded01a5c5b7e'    % 17
    };

gid = [
    1076537105      % 01
    1755712099      % 02
    1467086919      % 03
    922568942       % 04
    696576143       % 05
    1694335875      % 06
    567051619       % 07
    1927134163      % 08
    333865808       % 09
    2136302423      % 10
    1214851779      % 11
    1214851779      % 12   <-- NOTE: same GID as 11, check this
    903604567       % 13
    1342972763      % 14
    410768781       % 15
    291636615       % 16
    2024193406      % 17
    ];

url_fmt = ['https://docs.google.com/spreadsheets/d/e/' ...
           '2PACX-1vSsNKGE-16ACQ_fSPChaLZRRswo0tPqe_' ...
           'dghihqE9XDRUy0CF9-zvRThc0MqyteMN-T0IJ_YUeodkuC/' ...
           'pub?gid=%d&single=true&output=csv'];

if ~isfolder(out_dir), mkdir(out_dir); end
opts = weboptions('ContentType','text', 'Timeout',120);

%% One drifter at a time
for i = drifter_ids

    % Two drifters with the same GID would download the same data twice
    prev = drifter_ids(drifter_ids < i);
    dup  = prev(gid(prev) == gid(i));
    if ~isempty(dup)
        fprintf('SKIPPED drifter %02d: same GID as drifter %02d. Fix gid(%d) in this script.\n\n', i, dup(1), i);
        continue
    end

    % Download the whole tab as text
    try
        raw = webread(sprintf(url_fmt, gid(i)), opts);
    catch ME
        fprintf('FAILED drifter %02d: download failed (%s)\n\n', i, ME.message);
        continue
    end
    lines = splitlines(string(raw));
    lines = lines(strlength(strip(lines)) > 0);
    nL = numel(lines);

    % Split each line: column 3 holds the comma-separated message of 12 fields
    F   = strings(nL,12);
    dev = strings(nL,1);
    ok  = false(nL,1);
    for k = 1:nL
        fl = parseCsvLine(lines(k));
        if numel(fl) < 3, continue, end
        p = strip(split(fl(3), ','));
        if numel(p) >= 12
            F(k,:) = p(1:12)';
            dev(k) = fl(1);
            ok(k)  = true;
        end
    end
    if ~any(ok)
        fprintf('SKIPPED drifter %02d: no line has 12 comma-separated fields in column 3 (%d lines downloaded)\n\n', i, nL);
        continue
    end
    devFound = unique(dev(ok));
    F = F(ok,:);

    % Device ID check: the ID in column 1 should be this drifter's
    if all(devFound == string(deviceIDs{i}))
        devMsg = 'device ID matches';
    else
        devMsg = ['DEVICE ID MISMATCH, found ' char(strjoin(devFound, ', ')) ', expected ' deviceIDs{i}];
    end

    % Date: whole column first, row by row if any row does not match the format
    try
        dts = datetime(F(:,1), 'InputFormat', dateFmt);
    catch
        dts = NaT(size(F,1),1);
        for k = 1:size(F,1)
            try
                dts(k) = datetime(F(k,1), 'InputFormat', dateFmt);
            catch
            end
        end
    end

    % Time of day: hh:mm:ss or hh:mm:ss.SSS -> seconds
    secs = nan(size(F,1),1);
    for k = 1:size(F,1)
        tk = regexp(F(k,2), '^(\d+):(\d+):(\d+(?:\.\d+)?)$', 'tokens', 'once');
        if ~isempty(tk)
            secs(k) = str2double(tk{1})*3600 + str2double(tk{2})*60 + str2double(tk{3});
        end
    end

    good = ~isnat(dts) & ~isnan(secs);
    fprintf('Drifter %02d: %d lines, %d with 12 fields, %d with valid date and time | %s\n', ...
        i, nL, size(F,1), nnz(good), devMsg);
    if ~any(good)
        fprintf('SKIPPED drifter %02d: no valid date and time. Check dateFmt.\n\n', i);
        continue
    end
    dts = dts(good);  secs = secs(good);  F = F(good,:);

    % Sort in time and remove repeated times
    [~, ia] = unique(datenum(dts) + secs/86400);
    dts = dts(ia);  secs = secs(ia);  F = F(ia,:);

    % Time as text hh:mm:ss.SSS
    hh = floor(secs/3600);  mm = floor(mod(secs,3600)/60);  ss = mod(secs,60);
    timeTxt = string(arrayfun(@(a,b,c) sprintf('%02d:%02d:%06.3f', a, b, c), hh, mm, ss, 'UniformOutput', false));

    T = table(dts, timeTxt, str2double(F(:,3)), str2double(F(:,4)), F(:,5), ...
              str2double(F(:,6)), F(:,7), str2double(F(:,8)), str2double(F(:,9)), ...
              str2double(F(:,10)), str2double(F(:,11)), str2double(F(:,12)));

    outfile = fullfile(out_dir, sprintf('data_%02d%s.xlsx', i, suffix));
    if isfile(outfile), delete(outfile); end
    writetable(T, outfile, 'WriteVariableNames', false);

    fprintf('   saved %s | %d rows | %s to %s\n\n', outfile, height(T), ...
        char(dts(1) + seconds(secs(1)), 'dd-MMM-yyyy HH:mm:ss'), ...
        char(dts(end) + seconds(secs(end)), 'dd-MMM-yyyy HH:mm:ss'));
end

fprintf('Done. Set cfg.sheet_file_suffix = ''%s'' in get_obs_config for the periods that use these files.\n', suffix);

%% ======================== Local functions ==============================

function fields = parseCsvLine(line)
% Split one CSV line into fields, respecting "quoted, fields".

    tok = regexp(char(line), '(?:^|,)("(?:[^"]|"")*"|[^,]*)', 'tokens');
    fields = strings(1, numel(tok));

    for n = 1:numel(tok)
        c = tok{n}{1};
        if numel(c) >= 2 && c(1) == '"' && c(end) == '"'
            c = strrep(c(2:end-1), '""', '"');
        end
        fields(n) = string(c);
    end
end
