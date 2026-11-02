 function trimmed_drifters = trimDrifterDataByDeployment(drifters, deployments)
    % Trim each drifter's track based on its deployment time
    % Inputs:
    % - drifters: struct array with fields lat, lon, time
    % - deployments: table with fields drifter_id and time
    % Output:
    % - trimmed_drifters: struct array with trimmed lat/lon/time
    
    drifterIDs = cell2mat({drifters.ID}');
    trimmed_drifters = struct([]);
    rowx=1;
    for i = 1:numel(deployments.Starttime)
        idx = deployments.drifter_id(i);      % index of drifter
        if size(trimmed_drifters,2)<idx
            % don't do anything to rowx, until we fill first column
        elseif ~isempty(trimmed_drifters(rowx,idx).Lat)

            rowx = rowx+1;
        end
%         deploy_time = deployments.Startdattime(i);    % deployment timestamp
%         pull_time = deployments.Pulldattime(i);
        %% pair the current idx to appropriate drifter ID in structure
        idy = find(drifterIDs == idx);
        
        if isempty(idy)
            disp(['No drifter ID matching deployment log ID: ', num2str(idx)])
            continue
        end
        
        disp(['current notebook ID: ', num2str( idx)])
        disp(['current drifter ID: ', num2str(drifterIDs(idy))]);
        disp(numel(drifters))

        % Get drifter data
        drifterID = drifters(idy).ID;
        if drifterID~=idx
            disp('ERROR: drifter number does not match currrent deployment notes')
            keyboard
        end
        lat = drifters(idy).Lat;
        lon = drifters(idy).Lon;
        E = drifters(idy).East;
        N = drifters(idy).North;
        t = drifters(idy).Date_and_time;
        e_time = drifters(idy).etime;
        %t = datetime(t,'ConvertFrom','datenum');
%         A = drifters(idy).Alt;
        Velo = drifters(idy).Speed;
%         Ang = drifters(idy).Angle;
%         Volt = drifters(idy).Volt;
        v_e = drifters(idy).v_e;
        v_n = drifters(idy).v_n;
        V = drifters(idy).V;
        
        if isempty(t)
            trimmed_drifters(rowx,idx).Start_time = NaT;
            trimmed_drifters(rowx,idx).End_time = NaT;   
            continue
        end
%         
%         if idx == 5 || idx == 6 || idx == 10
% 
%         fprintf('\nDrifter %d BEFORE GPS masking\n',idx)
%         fprintf('Last raw timestamp: %s\n',string(t(end)))
%         fprintf('Number of observations: %d\n',numel(t))
% 
%         end
        
        
        latThresh = 0.0002;
        lonThresh = 0.0002;

        dLat = abs(diff(lat));
        dLon = abs(diff(lon));

        badIdx = [false; (dLat > latThresh) | (dLon > lonThresh)];

        % Remove the same observations from EVERY field
        lat(badIdx) = [];
        lon(badIdx) = [];
        E(badIdx) = [];
        N(badIdx) = [];
        t(badIdx) = [];
        e_time(badIdx) = [];
        Velo(badIdx) = [];
        v_e(badIdx) = [];
        v_n(badIdx) = [];
        V(badIdx) = [];
        
        buffer = seconds(30);

        target_start = deployments.Startdattime(i);  % Example: pick one from col1
        %[~, idy_start(i)] = min(abs(drifters(i).Date_and_time - target_start));

        target_pull = deployments.Pulldattime(i);  % Example: pick one from col1
        %[~, idy_pull(i)] = min(abs(drifters(1).Date_and_time - target_pull));
%         
%         if idx == 9
%             fprintf('\n================ DRIFTER 9 DEBUG ================\n')
%             fprintf('Structure index idy: %d\n', idy)
%             fprintf('Raw drifter ID: %d\n', drifters(idy).ID)
% 
%             fprintf('Target start: %s\n', string(target_start))
%             fprintf('Target pull : %s\n', string(target_pull))
% 
%             fprintf('First available timestamp: %s\n', string(t(1)))
%             fprintf('Last available timestamp : %s\n', string(t(end)))
%             fprintf('Number of timestamps: %d\n', numel(t))
%         end

        % trim the pull index to account for pull time being next deployment time!!
        % need to keep better track of pull times
        %idy_pull = idy_pull-12;
        idt = find( t >= target_start + buffer &...
                    t <= target_pull - buffer);
                
%         if idx == 9
%             fprintf('D9 points inside deployment window: %d\n', numel(idt))
% 
%             if isempty(idt)
%                 fprintf('*** D9 IS GETTING DROPPED HERE ***\n')
%             else
%                 fprintf('Trimmed start: %s\n', string(t(idt(1))))
%                 fprintf('Trimmed end  : %s\n', string(t(idt(end))))
%             end
% 
%             fprintf('=================================================\n\n')
%         end
                
        if isempty(idt)
            trimmed_drifters(rowx,idx).Start_time = NaT;
            trimmed_drifters(rowx,idx).End_time = NaT;
            
           continue 
        end

        % Trim data
        trimmed_drifters(rowx,idx).ID = drifters(idy).ID;
        trimmed_drifters(rowx,idx).Lat = lat(idt);
        trimmed_drifters(rowx,idx).Lon = lon(idt);
        trimmed_drifters(rowx,idx).East = E(idt);
        trimmed_drifters(rowx,idx).North = N(idt);
        trimmed_drifters(rowx,idx).Date_and_time = t(idt);
        trimmed_drifters(rowx,idx).Start_time = t(idt(1));
        trimmed_drifters(rowx,idx).End_time = t(idt(end));
        trimmed_drifters(rowx,idx).e_time = e_time(idt);
%         trimmed_drifters(rowx,idx).Alt = A(idt);
        trimmed_drifters(rowx,idx).Speed = Velo(idt);
%         trimmed_drifters(rowx,idx).Angle = Ang(idt);
%         trimmed_drifters(rowx,idx).Volt = Volt(idt);
        trimmed_drifters(rowx,idx).v_e = v_e(idt);
        trimmed_drifters(rowx,idx).v_n = v_n(idt);
        trimmed_drifters(rowx,idx).V = V(idt);

    end
end

