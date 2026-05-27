function ROI(adcp_level_2, ROI_ADCP, cfg)

%% U and V predicted from the complex EOF
U_pred_cmplx = real(adcp_level_2.EC_cmplx(:,1:3)*adcp_level_2.EOFs_cmplx(:,1:3)')+real(adcp_level_2.trend_cmplx);
V_pred_cmplx = imag(adcp_level_2.EC_cmplx(:,1:3)*adcp_level_2.EOFs_cmplx(:,1:3)')+imag(adcp_level_2.trend_cmplx);

% U_pred_cmplx = real(adcp_level_2.EC_cmplx*adcp_level_2.EOFs_cmplx');
% V_pred_cmplx = imag(adcp_level_2.EC_cmplx*adcp_level_2.EOFs_cmplx');
%% Unrotate the EOF predictions
[M,N] = size(U_pred_cmplx);

U_col = reshape(U_pred_cmplx,[],1);
V_col = reshape(V_pred_cmplx,[],1);

uv2 = [U_col,V_col]*adcp_level_2.rot_mat';
u2 = uv2(:,1);
v2 = uv2(:,2);

u2 = reshape(u2,M,N);
v2 = reshape(v2,M,N);

%%
for j = 1:size(ROI_ADCP,1);
    for i = 1:size(ROI_ADCP,2);

        if isempty(ROI_ADCP(j,i).inds)
            continue;
        end

        inds = ROI_ADCP(j,i).inds;
        inde = ROI_ADCP(j,i).inde;
        u2_ROI = mean(u2(inds:inde,:),1);
        v2_ROI = mean(v2(inds:inde,:),1);

        
        adcp_vn = ROI_ADCP(j,i).adcp_vn(:,1);
        adcp_dbins = ROI_ADCP(j,i).adcp_dbins;
        P = ROI_ADCP(j,i).adcp_P(1);
        adcp_ve = ROI_ADCP(j,i).adcp_ve(:,1);

        avg_vn = ROI_ADCP(j,i).avg_vn;
        avg_ve = ROI_ADCP(j,i).avg_ve;

        fig = figure;
        subplot(1,2,1)
        plot(adcp_vn,adcp_dbins,'k','LineWidth',2)
        hold on
        plot(v2_ROI(1,:),adcp_level_2.zonh_grid*P,'r','LineWidth',2)
        plot(avg_vn,P,'ro','LineWidth',2)
        hold off
        xlabel('m/s');
        ylabel('Depth (m)');
        title('v')
        % legend('ADCP','EOF','Drifter','Location','best');
        xlim([-0.7 0.7])
        ylim([0 6])
        grid on

        subplot(1,2,2)
        plot(adcp_ve,adcp_dbins,'k','LineWidth',2)
        hold on
        plot(u2_ROI(1,:),adcp_level_2.zonh_grid*P,'r','LineWidth',2)
        plot(avg_ve,P,'ro','LineWidth',2)
        hold off
        xlabel('m/s');
        ylabel('Depth (m)');
        title('u')
        legend('ADCP','EOF','Drifter','Location','northeast');
        xlim([-0.4 0.4])
        ylim([0 6])
        grid on
        
        
        fname = sprintf('_ROI_Profile_raw_j%02d_i%02d.png', j, i);
        fname_1 = sprintf('ROI_Profile_raw_j%02d_i%02d.fig', j, i);
%         exportgraphics(fig,fname, 'Resolution', 300);
%         savefig(fig,fname_1);


        exportgraphics(fig, fullfile(cfg.out.comp_figures, [cfg.obsTag fname]), 'Resolution', 300);







    end
end


end