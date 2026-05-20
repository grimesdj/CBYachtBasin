 function trimmed_drifters = trimDrifterDataByDeployment_SHEET(drifters, deployments)
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
        %t = datetime(t,'ConvertFrom','datenum');
%         A = drifters(idy).Alt;
%         Velo = drifters(idy).Speed;
%         Ang = drifters(idy).Angle;
%         Volt = drifters(idy).Volt;
        v_e = drifters(idy).v_e;
        v_n = drifters(idy).v_n;
        V = drifters(idy).V;
        
        

        target_start = deployments.Startdattime(i);  % Example: pick one from col1
        %[~, idy_start(i)] = min(abs(drifters(i).Date_and_time - target_start));

        target_pull = deployments.Pulldattime(i);  % Example: pick one from col1
        %[~, idy_pull(i)] = min(abs(drifters(1).Date_and_time - target_pull));

        % trim the pull index to account for pull time being next deployment time!!
        % need to keep better track of pull times
        %idy_pull = idy_pull-12;
        idt = find( t >= target_start &...
                    t <= target_pull);
                
        if isempty(idt)
            trimmed_drifters(rowx,idx).Start_time = NaT;
            trimmed_drifters(rowx,idx).End_time = NaT;
            
           continue 
        end

        % Trim data
        trimmed_drifters(rowx,idx).Lat = lat(idt);
        trimmed_drifters(rowx,idx).Lon = lon(idt);
        trimmed_drifters(rowx,idx).East = E(idt);
        trimmed_drifters(rowx,idx).North = N(idt);
        trimmed_drifters(rowx,idx).Date_and_time = t(idt);
        trimmed_drifters(rowx,idx).Start_time = t(idt(1));
        trimmed_drifters(rowx,idx).End_time = t(idt(end));
%         trimmed_drifters(rowx,idx).Alt = A(idt);
%         trimmed_drifters(rowx,idx).Speed = Velo(idt);
%         trimmed_drifters(rowx,idx).Angle = Ang(idt);
%         trimmed_drifters(rowx,idx).Volt = Volt(idt);
        trimmed_drifters(rowx,idx).v_e = v_e(idt);
        trimmed_drifters(rowx,idx).v_n = v_n(idt);
        trimmed_drifters(rowx,idx).V = V(idt);

    end
end

