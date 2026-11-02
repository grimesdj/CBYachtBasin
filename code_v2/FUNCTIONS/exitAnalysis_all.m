%% exitAnalysis_v2  -  mouth-exit probability from UNTRIMMED level 2 drifter files
% Outcome per track:
%   OutsideStart : first sample is in Region 0 (excluded from all denominators)
%   MouthExit    : first Region 1 -> Region 0 transition anywhere in the track
%   NoMouthExit  : everything else (lateral marsh/creek exits count here, as do
%                  tracks still inside when the drifter was pulled = censored)
% Results saved to ExitAnalysis_results.mat
clear; close all; rng(1)

%% Settings
roots = {'C:\Users\Ben\OneDrive - UNC-Wilmington\', ...
         'C:\Users\bcm3620\OneDrive - UNC-Wilmington\'};
k = find(cellfun(@(r) isfolder(fullfile(r,'CB_YachtBasin')), roots), 1);
if isempty(k), error('Neither root path contains CB_YachtBasin on this machine.'); end
root_data = roots{k};
fprintf('Using root: %s\n', root_data);
% Do you want to use the SHEEET (5 min) or onboard (1hz w various
% differencing window options) for the FallFrolic drifters??
% names   = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};
names   = {'MarshMadness','MarshMayhem','FallFrolic_100825_sheet','FallFrolic_101025_sheet','NovDep','JuneJamboree','SummerFlood'};
folders = {'MarshMadness','MarshMayhem','FallFrolic','FallFrolic','NovDep','JuneJamboree','SummerFlood'};
% Add when FallFlood level 2 exists:
% names{end+1} = 'FallFlood';  folders{end+1} = 'FallFlood';

horizons_hr = [0.5 1 2 3 6.2];
phases      = ["High","Low","Rising","Falling"];
nBoot       = 2000;
mouthRegion = 1;

%% Build one row per track
rows = struct('Period',{},'Deployment',{},'Track',{},'ID',{},'ReleaseTime',{}, ...
              'StartRegion',{},'StartPhase',{},'Duration',{},'Outcome',{}, ...
              'ExitTime',{},'ExitPhase',{},'AnyStuck',{},'StuckBeforeExit',{});
for i = 1:numel(names)
    fn = fullfile(root_data,'CB_YachtBasin',folders{i},'PROCESSED_DATA_V2','DRIFTERS',[names{i} '_drifters_level_2.mat']);
    if ~isfile(fn), fprintf('%s: file not found, skipped\n',names{i}); continue, end
    L2 = load(fn).drifters_level_2;
    [nD,nC] = size(L2);
    for d = 1:nD
        for k = 1:nC
            D = L2(d,k);
            if isempty(D.Regions) || isempty(D.Date_and_time), continue, end
            R  = D.Regions(:);
            t  = seconds(D.Date_and_time(:) - D.Date_and_time(1));
            ph = string(D.TidePhase(:));
            r = struct();
            r.Period = string(names{i});  r.Deployment = d;  r.Track = k;
            if isfield(D,'ID') && ~isempty(D.ID), r.ID = D.ID; else, r.ID = NaN; end
            r.ReleaseTime = D.Date_and_time(1);
            r.StartRegion = R(1);  r.StartPhase = ph(1);  r.Duration = t(end);
            r.ExitTime = NaN;      r.ExitPhase = "";
            e = find(R(1:end-1)==mouthRegion & R(2:end)==0, 1, 'first');
            if R(1) <= 0
                r.Outcome = "OutsideStart";
            elseif ~isempty(e)
                r.Outcome = "MouthExit";  r.ExitTime = t(e+1);  r.ExitPhase = ph(e+1);
                        else
                r.Outcome = "NoMouthExit";
            end
            if ~isfield(D,'Stuck') || numel(D.Stuck) ~= numel(R)
                error('%s has no usable Stuck field in level 2. Run flagStuckDrifters first.', names{i});
            end
            sk = logical(D.Stuck(:));                              % true where this track was flagged stuck
            r.AnyStuck = any(sk);                                  % true if the track has any stuck sample
            r.StuckBeforeExit = ~isempty(e) && any(sk(1:e+1));     % true if a stuck sample comes before the exit
            rows(end+1,1) = r; %#ok<SAGROW>
        end
    end
end
tbl = struct2table(rows);
tbl.Cluster = findgroups(tbl.Period, tbl.Deployment);   % bootstrap unit = period x deployment

fprintf('\nTrack outcomes (all periods):\n');
disp(groupsummary(tbl,'Outcome'))
ins = tbl.StartRegion > 0;
exited = tbl.Outcome == "MouthExit";

%% Check: stuck samples in the exit analysis
fprintf('\n=== Stuck samples in the exit analysis ===\n');
fprintf('Tracks that started inside: %d, with any stuck sample: %d\n', sum(ins), sum(ins & tbl.AnyStuck));
fprintf('Mouth exits: %d, with a stuck sample before the exit: %d\n', sum(ins & exited), sum(ins & exited & tbl.StuckBeforeExit));
fprintf('Tracks that never exited: %d, with any stuck sample: %d\n', sum(ins & ~exited), sum(ins & ~exited & tbl.AnyStuck));

%% Exit percentages: whole basin, by starting region, by starting tidal phase
fprintf('\n=== Share of tracks that left through the mouth ===\n');
fprintf('(x of n = tracks that left through the mouth, out of tracks that started inside the basin)\n');

x = sum(ins & exited);
n = sum(ins);
fprintf('Whole basin: %d of %d = %.1f%%\n', x, n, 100*x/n);

fprintf('\nBy starting region:\n');
for r = unique(tbl.StartRegion(ins))'
    sel = ins & tbl.StartRegion == r;
    x = sum(sel & exited);
    n = sum(sel);
    fprintf('Region %d: %d of %d = %.1f%%\n', r, x, n, 100*x/n);
end

fprintf('\nBy starting tidal phase:\n');
for p = phases
    sel = ins & tbl.StartPhase == p;
    x = sum(sel & exited);
    n = sum(sel);
    fprintf('%s: %d of %d = %.1f%%\n', p, x, n, 100*x/n);
end

fprintf('\nBy starting region and starting tidal phase (exits / tracks):\n');
fprintf('%-8s', 'Region');
for p = phases
    fprintf(' %-14s', p);
end
fprintf('\n');
for r = unique(tbl.StartRegion(ins))'
    fprintf('%-8d', r);
    for p = phases
        sel = ins & tbl.StartRegion == r & tbl.StartPhase == p;
        x = sum(sel & exited);
        n = sum(sel);
        if n == 0
            fprintf(' %-14s', '--');
        else
            fprintf(' %-14s', sprintf('%d/%d = %.0f%%', x, n, 100*x/n));
        end
    end
    fprintf('\n');
end

%% Surface current U at the moment of each mouth exit
smooth_samples = 2;   % number of ADCP samples averaged to smooth U (about 20 min)
tbl.ExitDatenum = datenum(tbl.ReleaseTime + seconds(tbl.ExitTime));   % clock time of each exit, as a number
tbl.ExitU = nan(height(tbl),1);                                       % U at each exit (filled in below)

for i = 1:numel(names)
    sel = find(exited & tbl.Period == string(names{i}));   % mouth exits in this observation period
    if isempty(sel), continue, end
    adcpDir = fullfile(root_data,'CB_YachtBasin',folders{i},'PROCESSED_DATA_V2','ADCP');
    adcpName = erase(names{i},'_sheet');   % ADCP files do not have _sheet in the name
    A1 = load(fullfile(adcpDir,[adcpName '_NorthADCP_adcp_level_1.mat'])).adcp_level_1;
    A2 = load(fullfile(adcpDir,[adcpName '_NorthADCP_adcp_level_2.mat'])).adcp_level_2;
    tADCP = A1.time(:);                                                    % times of the ADCP profiles
    Usurf = movmean(A2.U_mapped(end,:)', smooth_samples, 'omitnan');      % smoothed surface U at those times
    tbl.ExitU(sel) = interp1(tADCP, Usurf, tbl.ExitDatenum(sel));         % U at the exit times
end

fprintf('\nSurface current U at each mouth exit (negative = ebb, toward the mouth):\n');
for k = find(exited)'
    fprintf('%-18s | Drifter %2d | %s | %-7s | U = %+.3f m/s\n', ...
        tbl.Period(k), tbl.ID(k), datestr(tbl.ExitDatenum(k),'mm/dd HH:MM'), tbl.ExitPhase(k), tbl.ExitU(k));
end

fprintf('\nSurface current U at mouth exits, by tidal phase at exit (negative = ebb, toward the mouth):\n');
for p = phases
    sel = exited & tbl.ExitPhase == p & ~isnan(tbl.ExitU);
    if ~any(sel)
        fprintf('%-8s no exits\n', p);
        continue
    end
    u = tbl.ExitU(sel);
    fprintf('%-8s n = %2d | mean %+.3f | min %+.3f | max %+.3f m/s\n', p, sum(sel), mean(u), min(u), max(u));
end
u = tbl.ExitU(exited & ~isnan(tbl.ExitU));
fprintf('All exits: n = %d | mean %+.3f | min %+.3f | max %+.3f m/s | ebb (negative) %d of %d\n', ...
    numel(u), mean(u), min(u), max(u), sum(u < 0), numel(u));

%% Fixed-horizon exit probability (OFF: biased for this dataset because tracks
%  end soon after exiting; use the Kaplan-Meier sections below instead)
showFixedHorizon = false;
Hres = [];
if showFixedHorizon
fprintf('\n=== P(mouth exit within T of release) ===\n');
fprintf('Only tracks observed for at least T are counted.\n');
Hres = struct('H',{},'x',{},'n',{},'p',{},'ci',{},'bootci',{});
fprintf('\n%-6s %8s %8s %7s   %-16s   %-16s\n','T(h)','eligible','exits','p','Wilson 95%','Bootstrap(deployment) 95%');
for h = 1:numel(horizons_hr)
    H = horizons_hr(h)*3600;
    elig = ins & (tbl.Duration >= H);
    ev   = elig & exited & tbl.ExitTime <= H;
    x = sum(ev); n = sum(elig);
    ci = wilson(x,n);
    % bootstrap by cluster
    cx = accumarray(tbl.Cluster, ev);  cn = accumarray(tbl.Cluster, elig);
    G = numel(cx);  pb = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(G,G,1);
        if sum(cn(idx)) > 0, pb(b) = sum(cx(idx))/sum(cn(idx)); end
    end
    pb = sort(pb(~isnan(pb)));
    if numel(pb) > 20
        bci = pb(max(1,round([0.025 0.975]*numel(pb))))';
    else
        bci = [NaN NaN];
    end
    Hres(h) = struct('H',horizons_hr(h),'x',x,'n',n,'p',x/max(n,1),'ci',ci,'bootci',bci);
    fprintf('%-6.1f %8d %8d %7.3f   [%5.3f %5.3f]      [%5.3f %5.3f]\n', ...
        horizons_hr(h), n, x, x/max(n,1), ci(1), ci(2), bci(1), bci(2));
end

%% By starting tidal phase and starting region, for each horizon
fprintf('\n=== By STARTING tidal phase ===\n');
for h = 1:numel(horizons_hr)
    H = horizons_hr(h)*3600;
    fprintf('\nT = %.1f h\n',horizons_hr(h));
    for p = phases
        sel  = ins & tbl.StartPhase == p;
        elig = sel & (tbl.Duration >= H);
        ev   = elig & exited & tbl.ExitTime <= H;
        printRow(char(p), sum(ev), sum(elig));
    end
end
fprintf('\n=== By STARTING region ===\n');
for h = 1:numel(horizons_hr)
    H = horizons_hr(h)*3600;
    fprintf('\nT = %.1f h\n',horizons_hr(h));
    for r = unique(tbl.StartRegion(ins))'
        sel  = ins & tbl.StartRegion == r;
        elig = sel & (tbl.Duration >= H);
        ev   = elig & exited & tbl.ExitTime <= H;
        printRow(sprintf('Region %d',r), sum(ev), sum(elig));
    end
end

end   % showFixedHorizon

%% Phase at the moment of exit (descriptive only, no denominator)
fprintf('\n=== Tidal phase at the moment of mouth exit ===\n');
for p = phases
    fprintf('%-8s %d exits\n', p, sum(exited & tbl.ExitPhase == p));
end

%% Cumulative exit probability, no horizon needed (Kaplan-Meier)
obsT = tbl.Duration(ins);
evK  = exited(ins);
obsT(evK) = tbl.ExitTime(ins & exited);
tu = unique(obsT(evK));
S = 1; Scum = nan(numel(tu),1);
for j = 1:numel(tu)
    ni = sum(obsT >= tu(j));
    di = sum(obsT == tu(j) & evK);
    S = S*(1 - di/ni);
    Scum(j) = 1 - S;
end
fprintf('\n=== Cumulative mouth-exit probability (Kaplan-Meier, pooled) ===\n');
fprintf('%-10s %10s %10s\n','minutes','P(exit)','at risk');
for q = [15 30 60 90 120 180]
    j = find(tu <= q*60, 1, 'last');
    if isempty(j), pq = 0; else, pq = Scum(j); end
    fprintf('%-10d %10.3f %10d\n', q, pq, sum(obsT >= q*60));
end

figure
subplot(2,1,1)
stairs([0; tu/60], [0; Scum], 'LineWidth',2); grid on
ylabel('P(exit via mouth)'); title('Cumulative exit probability since release (pooled)')
subplot(2,1,2)
tg = 0:5:max(obsT)/60;
plot(tg, arrayfun(@(q) sum(obsT >= q*60), tg), 'LineWidth',2); grid on
xlabel('Minutes since release'); ylabel('Tracks still observed')

exportgraphics(gcf, fullfile('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES\EXITS\','KM_pooled.png'), 'Resolution', 300);

%% Kaplan-Meier by starting region and phase (uses every track)
fprintf('\n=== Kaplan-Meier P(mouth exit) by START region / phase ===\n');
obsAll = tbl.Duration;  obsAll(exited) = tbl.ExitTime(exited);
for q = [60 120 180]*60
    fprintf('\nWithin %d min\n', q/60);
    for r = unique(tbl.StartRegion(ins))'
        s = ins & tbl.StartRegion == r;
        [P,nr] = kmAt(obsAll(s), exited(s), q);
        fprintf('  Region %d  P = %.3f  at risk %2d  (n=%d, exits=%d)\n', r, P, nr, sum(s), sum(s & exited));
    end
    for p = phases
        s = ins & tbl.StartPhase == p;
        if ~any(s), fprintf('  %-8s no data\n', p); continue, end
        [P,nr] = kmAt(obsAll(s), exited(s), q);
                fprintf('  %-8s  P = %.3f  at risk %2d  (n=%d, exits=%d)\n', p, P, nr, sum(s), sum(s & exited));
    end
    fprintf('  By starting region and phase:\n');
    for r = unique(tbl.StartRegion(ins))'
        for p = phases
            s = ins & tbl.StartRegion == r & tbl.StartPhase == p;
            if ~any(s), continue, end
            [P,nr] = kmAt(obsAll(s), exited(s), q);
            fprintf('    Region %d %-8s  P = %.3f  at risk %2d  (n=%d, exits=%d)\n', r, p, P, nr, sum(s), sum(s & exited));
        end
    end
end

%% Save
save(fullfile('C:\Users\bcm3620\OneDrive - UNC-Wilmington\THESIS\THESIS_DRAFTING\FIGURES\EXITS\','ExitAnalysis_results.mat'), 'tbl','Hres','tu','Scum','horizons_hr');
fprintf('\nSaved ExitAnalysis_results.mat\n');

%% Local functions
function ci = wilson(x,n)
    if n == 0, ci = [NaN NaN]; return, end
    z = 1.96; p = x/n; den = 1 + z^2/n;
    c = (p + z^2/(2*n))/den;
    h = z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den;
    ci = [max(0,c-h) min(1,c+h)];
end

function [P,nrisk] = kmAt(obsT, evK, q)
    tu = unique(obsT(evK)); S = 1;
    for j = 1:numel(tu)
        if tu(j) > q, break, end
        S = S*(1 - sum(obsT==tu(j) & evK)/sum(obsT >= tu(j)));
    end
    P = 1 - S;  nrisk = sum(obsT >= q);
end

function printRow(label,x,n)
    if n == 0
        fprintf('  %-10s no data\n',label);
    else
        ci = wilson(x,n);
        fprintf('  %-10s %3d of %3d = %5.3f   Wilson [%5.3f %5.3f]\n',label,x,n,x/n,ci(1),ci(2));
    end
end