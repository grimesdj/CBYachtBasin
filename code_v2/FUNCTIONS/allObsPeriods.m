function P = allObsPeriods()
%ALLOBSPERIODS  The one list of observation periods used by the
% all-period analyses (comparisonAll, comparisonAllPlots).
% Add or comment out a row here and both scripts follow.
%
% Columns: data folder | cfg.name | cfg.obsTag

L = { ...
    'MarshMadness',  'MarshMadness',       'MarshMadness_NorthADCP'
    'MarshMayhem',   'MarshMayhem',        'MarshMayhem_NorthADCP'
    'FallFrolic',    'FallFrolic_100825',  'FallFrolic_100825_NorthADCP'
    'FallFrolic',    'FallFrolic_101025',  'FallFrolic_101025_NorthADCP'
%   'NovDep',        'NovDep',             'NovDep_NorthADCP'
    'JuneJamboree',  'JuneJamboree',       'JuneJamboree_NorthADCP'
    'SummerFlood',   'SummerFlood',        'SummerFlood_NorthADCP'
%   'FallFlood',     'FallFlood',          'FallFlood_NorthADCP'
    };

P = struct('folder', L(:,1), 'name', L(:,2), 'obsTag', L(:,3));

end
