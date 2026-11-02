function C = loadPeriodComp(root_data, Pn)
%LOADPERIODCOMP  Load one period's Level 3 COMP structure for the
% all-period analyses. Returns [] if the period has no output or no ROI
% crossings; errors if the output is from the old comparison.m.
%
%   Pn is one element of allObsPeriods().

f = fullfile(root_data, Pn.folder, 'PROCESSED_DATA_V2', 'COMP', [Pn.obsTag '_COMP.mat']);

if ~isfile(f)
    warning('loadPeriodComp:missing', 'No Level 3 output for %s (%s); skipped.', Pn.obsTag, f);
    C = [];
    return
end

S = load(f, 'COMP');
C = S.COMP;

if isempty(C)
    fprintf('%s: no ROI crossings.\n', Pn.obsTag);
    return
end

if ~isfield(C, 'pairs')
    error('loadPeriodComp:oldFormat', ...
        ['%s was written by the old comparison.m.\n' ...
         'Re-run Levels 2-4 (adcp_post_processing, comparison, ROI) for %s first.'], ...
        f, Pn.obsTag);
end

end
