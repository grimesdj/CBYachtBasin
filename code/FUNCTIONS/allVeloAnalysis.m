function [all_RMSE] = allVeloAnalysis(cfg)
%% Load all the Temporal res data from all obs periods
obsFiles = {'MarshMadness\PROCESSED_DATA\DRIFTERS\MarshMadness_TemporalRes_RMSE.mat'
    'MarshMayhem\PROCESSED_DATA\DRIFTERS\MarshMayhem_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_100825_TemporalRes_RMSE.mat'
    'FallFrolic\PROCESSED_DATA\DRIFTERS\FallFrolic_101025_TemporalRes_RMSE.mat'
    'NovDep\PROCESSED_DATA\DRIFTERS\NovDep_TemporalRes_RMSE.mat'
    'JuneJamboree\PROCESSED_DATA\DRIFTERS\JuneJamboree_TemporalRes_RMSE.mat'
    'SummerFlood\PROCESSED_DATA\DRIFTERS\SummerFlood_TemporalRes_RMSE.mat'};

names = {'MarshMadness','MarshMayhem','FallFrolic_100825','FallFrolic_101025','NovDep','JuneJamboree','SummerFlood'};

for i = 1:numel(obsFiles)

    filename = fullfile(cfg.root_data, obsFiles{i});

    S = load(filename);

    data.(names{i}) = S.Results;

end

obs_names = fieldnames(data);

all_vals_RMSE = [];

all_vals_RMSE_v_Speed = [];

all_vals_RMSE_v_smoothed_Speed = [];

figure;
hold on

for i = 1:numel(obs_names)
    
   obs = obs_names{i};
   
   R = data.(obs);
   
   obs_vals_RMSE = vertcat(R.RMSE);
   
   obs_vals_RMSE_v_Speed = vertcat(R.RMSE_v_Speed);
   
   obs_vals_RMSE_v_smoothed_Speed = vertcat(R.RMSE_v_Speed_Smoothed);
   
   all_vals_RMSE = [all_vals_RMSE; obs_vals_RMSE];
   
   all_vals_RMSE_v_Speed = [all_vals_RMSE_v_Speed; obs_vals_RMSE_v_Speed];
   
   all_vals_RMSE_v_smoothed_Speed = [all_vals_RMSE_v_smoothed_Speed; obs_vals_RMSE_v_smoothed_Speed];
   
   for k = 1:numel(R)
       
      scatter(R(k).Resolution_sec,R(k).RMSE,20,'b','filled','HandleVisibility','off')
      
      scatter(R(k).Resolution_sec,R(k).RMSE_v_Speed,20,'r','filled','HandleVisibility','off')
      
      scatter(R(k).Resolution_sec,R(k).RMSE_v_Speed_Smoothed,20,'g','filled','HandleVisibility','off')
      
      
       
   end
    
    
end

mean_RMSE_all = mean(all_vals_RMSE,1,'omitnan')
mean_RMSE_v_Speed_all = mean(all_vals_RMSE_v_Speed,1,'omitnan')
mean_RMSE_v_smoothed_Speed_all = mean(all_vals_RMSE_v_smoothed_Speed,1,'omitnan')
Resolutions = R(1).Resolution_sec;

 % Dummy points just for legend
h1 = scatter(nan,nan,20,'b','filled');
h2 = scatter(nan,nan,20,'r','filled');
h3 = scatter(nan,nan,20,'g','filled');
h4 = plot(Resolutions,mean_RMSE_all,'b-','LineWidth',2);
h5 = plot(Resolutions,mean_RMSE_v_Speed_all,'r-','LineWidth',2);
h6 = plot(Resolutions,mean_RMSE_v_smoothed_Speed_all,'g-','LineWidth',2);

xlabel('Temporal Resolution (s)')
ylabel('RMSE (m/s)')

legend([h1 h2 h3 h4 h5 h6], ...
  'RMSE vs Max Lag', ...
  'RMSE vs Doppler Speed', ...
  'RMSE vs Smoothed Doppler Speed', ...
  'Mean RMSE v maxLag', ...
  'Mean RMSE v Doppler Speed',...
  'Mean RMSE v Smoothed Doppler Speed',...
  'Location','best')

title(sprintf('RMSE Across All Trajectories - %s', cfg.name))
grid on

%% RMSE as a function of drifter velocity across all observation periods

% Velocity bins
vel_edges = [0 0.05 0.1 0.3 0.6];

vel_centers = 0.5 * ...
    (vel_edges(1:end-1) + vel_edges(2:end));

vel_labels = { ...
    '0-0.05', ...
    '0.05-0.1', ...
    '0.1-0.3', ...
    '0.3-0.6'};

nBins = length(vel_edges) - 1;
nRes  = length(Resolutions);

% Containers
RMSE_byVelocity_MaxLag = nan(nBins,nRes);
RMSE_byVelocity_Speed_1Hz  = nan(nBins,nRes);
RMSE_byVelocity_Speed_smoothed  = nan(nBins,nRes);

N_byVelocity = zeros(nBins,nRes);
meanSpeed_byVelocity = nan(nBins,nRes);
mean_smoothed_Speed_byVelocity = nan(nBins,nRes);


%% Loop over each temporal resolution

for rr = 1:nRes

    all_Vcurrent     = [];
    all_Speedcurrent = [];
    all_Vmax         = [];
    all_Speed        = [];


    %% Gather all observations from all observation periods

    for i = 1:numel(obs_names)

        obs = obs_names{i};

        R = data.(obs);

        for k = 1:numel(R)

            % Skip empty trajectory entries
            if isempty(R(k).Resolution_sec) || ...
               isempty(R(k).Speed)
                continue
            end


            %% Current temporal-resolution velocity

            if rr == 1

                % 1 Hz position-derived velocity
                V_current = R(k).V_1Hz;
                
                % At 1Hz the corresponding Doppler speed is just the 
                % original centered 1Hz Doppler Speed
                Speed_current = R(k).Speed;

            else

                % TemporalRes(1) = 2 s
                % TemporalRes(2) = 4 s
                % etc.
                V_current = ...
                    R(k).TemporalRes(rr-1).V;
                
                Speed_current = ...
                    R(k).TemporalRes(rr-1).Speed;

            end


            %% Maximum-lag reference velocity

            V_max = ...
                R(k).TemporalRes(end).V;


            %% 1 Hz Doppler speed

            Speed = R(k).Speed;


            % Force column vectors
            V_current     = V_current(:);
            Speed_current = Speed_current(:);
            V_max         = V_max(:);
            Speed         = Speed(:);


            %% Make sure vectors are same length

            n = min([ ...
                length(V_current), ...
                length(Speed_current), ...
                length(V_max), ...
                length(Speed)]);

            V_current = V_current(1:n);
            Speed_current = Speed_current(1:n);
            V_max     = V_max(1:n);
            Speed     = Speed(1:n);


            %% Append

            all_Vcurrent = ...
                [all_Vcurrent; V_current];
            
            all_Speedcurrent = ...
                [all_Speedcurrent; Speed_current];

            all_Vmax = ...
                [all_Vmax; V_max];

            all_Speed = ...
                [all_Speed; Speed];

        end
        
        

    end
    
        fprintf('\nResolution = %d s\n',Resolutions(rr))

        fprintf('Mean abs difference Speed vs smoothed Speed = %.8f m/s\n', ...
            mean(abs(all_Speed - all_Speedcurrent),'omitnan'));

        fprintf('Max abs difference Speed vs smoothed Speed = %.8f m/s\n', ...
            max(abs(all_Speed - all_Speedcurrent),[],'omitnan'));


    %% Bin observations according to Doppler speed

    vel_bin = discretize(all_Speed,vel_edges);


    %% Calculate RMSE within each velocity bin

    for b = 1:nBins

        idx = vel_bin == b;

        meanSpeed_byVelocity(b,rr) = ...
            mean(all_Speed(idx),'omitnan');
        
        mean_smoothed_Speed_byVelocity(b,rr) = ...
            mean(all_Speedcurrent(idx),'omitnan');
        
        N_byVelocity(b,rr) = sum(idx);


        if ~any(idx)
            continue
        end


        %% RMSE relative to maximum lag

        diff_max = ...
            all_Vcurrent(idx) - all_Vmax(idx);

        RMSE_byVelocity_MaxLag(b,rr) = ...
            sqrt(mean(diff_max.^2,'omitnan'));


        %% RMSE relative to Doppler Speed

        diff_speed = ...
            all_Vcurrent(idx) - all_Speed(idx);

        RMSE_byVelocity_Speed_1Hz(b,rr) = ...
            sqrt(mean(diff_speed.^2,'omitnan'));
        
        %% RMSE relative to Smoothed Doppler Speed

        diff_speed_smoothed = ...
            all_Vcurrent(idx) - all_Speedcurrent(idx);

        RMSE_byVelocity_Speed_smoothed(b,rr) = ...
            sqrt(mean(diff_speed_smoothed.^2,'omitnan'));

    end

end

%% Normalized RMSE

NRMSE_byVelocity_Speed = ...
    RMSE_byVelocity_Speed_1Hz ./ meanSpeed_byVelocity;

% Convert to percent
NRMSE_percent = ...
    NRMSE_byVelocity_Speed * 100;

%% Bin counts
N_bins = N_byVelocity(:,1);

labels_with_N = cell(nBins,1);

for b = 1:nBins
    labels_with_N{b} = sprintf('%s (N=%d)', ...
        vel_labels{b}, ...
        N_bins(b));
end

xpos = 1:nBins;



%% Figure: RMSE vs drifter velocity relative to Doppler Speed 

figure
hold on

for rr = 1:nRes

    plot( ...
        xpos, ...
        RMSE_byVelocity_Speed_1Hz(:,rr), ...
        '-o', ...
        'LineWidth',1.5, ...
        'DisplayName', ...
        sprintf('%d s',Resolutions(rr)));
    
    

end

  plot( ...  
      xpos, ...
      RMSE_byVelocity_Speed_smoothed(:,5), ...
      '--*', ...
      'LineWidth',1.5, ...
      'DisplayName', ...
      sprintf('%d s',Resolutions(5)));

xlabel('Doppler Drifter Velocity (m/s)')
ylabel('RMSE Relative to Doppler Speed (m/s)')

xticks(xpos)
xticklabels(labels_with_N)

legend('Location','bestoutside')

title(sprintf( ...
    'RMSE vs Drifter Velocity - %s', ...
    cfg.name))

grid on


%% Figure: RMSE vs drifter velocity relative to maximum lag

figure
hold on

for rr = 1:nRes

    plot( ...
        xpos, ...
        RMSE_byVelocity_MaxLag(:,rr), ...
        '-o', ...
        'LineWidth',1.5, ...
        'DisplayName', ...
        sprintf('%d s',Resolutions(rr)));

end

xlabel('Doppler Drifter Velocity (m/s)')
ylabel('RMSE Relative to Max Lag (m/s)')

xticks(xpos)
xticklabels(labels_with_N)

legend('Location','bestoutside')

title(sprintf( ...
    'RMSE vs Drifter Velocity - Max Lag Reference - %s', ...
    cfg.name))

grid on

%% NRMSE figure

figure
hold on



for rr = 1:nRes

    plot( ...
        xpos, ...
        NRMSE_percent(:,rr), ...
        '-o', ...
        'LineWidth',1.5, ...
        'DisplayName', ...
        sprintf('%d s',Resolutions(rr)));

end

xlabel('Doppler Drifter Velocity (m/s)')
ylabel('Normalized RMSE (%)')

xticks(xpos)
xticklabels(labels_with_N)

legend('Location','bestoutside')

title(sprintf( ...
    'Normalized RMSE vs Drifter Velocity - %s', ...
    cfg.name))

grid on

all_RMSE = data;


end
