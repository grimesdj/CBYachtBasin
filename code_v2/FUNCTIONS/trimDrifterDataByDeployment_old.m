 function trimmed_drifters = trimDrifterDataByDeployment(drifters, deployments)
    % Trim each drifter's track based on its deployment time
    % Inputs:
    % - drifters: struct array with fields lat, lon, time
    % - deployments: table with fields drifter_id and time
    % Output:
    % - trimmed_drifters: struct array with trimmed lat/lon/time

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

        disp(idx)
        disp(numel(drifters))

        % Get drifter data
        lat = drifters(idx).Lat;
        lon = drifters(idx).Lon;
        E = drifters(idx).East;
        N = drifters(idx).North;
        t = drifters(idx).Date_and_time;
        t = datetime(t,'ConvertFrom','datenum');
        A = drifters(idx).Alt;
        Velo = drifters(idx).Speed;
        Ang = drifters(idx).Angle;
        Volt = drifters(idx).Volt;
        v_e = drifters(idx).v_e;
        v_n = drifters(idx).v_n;
        V = drifters(idx).V;
        
        

        target_start = deployments.Startdattime(i);  % Example: pick one from col1
        %[~, idx_start(i)] = min(abs(drifters(i).Date_and_time - target_start));

        target_pull = deployments.Pulldattime(i);  % Example: pick one from col1
        %[~, idx_pull(i)] = min(abs(drifters(1).Date_and_time - target_pull));

        % trim the pull index to account for pull time being next deployment time!!
        % need to keep better track of pull times
        %idx_pull = idx_pull-12;
        idt = find( t >= target_start &...
                    t <= target_pull);

        % Trim data
        trimmed_drifters(rowx,idx).Lat = lat(idt);
        trimmed_drifters(rowx,idx).Lon = lon(idt);
        trimmed_drifters(rowx,idx).East = E(idt);
        trimmed_drifters(rowx,idx).North = N(idt);
        trimmed_drifters(rowx,idx).Date_and_time = t(idt);
        trimmed_drifters(rowx,idx).Alt = A(idt);
        trimmed_drifters(rowx,idx).Speed = Velo(idt);
        trimmed_drifters(rowx,idx).Angle = Ang(idt);
        trimmed_drifters(rowx,idx).Volt = Volt(idt);
        trimmed_drifters(rowx,idx).v_e = v_e(idt);
        trimmed_drifters(rowx,idx).v_n = v_n(idt);
        trimmed_drifters(rowx,idx).V = V(idt);

    end
end

