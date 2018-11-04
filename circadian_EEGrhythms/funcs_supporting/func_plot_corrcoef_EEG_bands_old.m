


% To run, drop this into the appropriate sectoin of run_all_ergodic.m        
%         if analyze_only_seizing == 1
%             [T X] = rat_2_ts (r_seiz,prepostchronic, smoothmode_ts);
%         elseif analyze_only_seizing == 2
%             [T X] = rat_2_ts (r_ctrl,prepostchronic, smoothmode_ts);
%         else
%             [T X] = rat_2_ts (r_ns,prepostchronic, smoothmode_ts);
%         end
%         
%         func_plot_corrcoef_EEG_bands_old(T,X,phi,prepostchronic)



function func_plot_corrcoef_EEG_bands_old(T,X,phi,prepostchronic)
    plot_individual_rats = 1;
        plot_vs_phi = 0; % if 1, have the x and y axes of the colour plots be phase. Otherwise, they're freq band.
        
    ratrange = 1:length(T);
    corrcoef_arr = [];
    
    for j = ratrange
        [coef, score, latent] = princomp(X{j});
        corrcoef_PC = (coef)*diag(latent)*inv(coef);
        if plot_individual_rats == 1
            
            if ~plot_vs_phi
                figure;
                colourmapval = [summer; flipud(autumn)]; colormap(colourmapval);
                %colormap winter;
                imagesc(corrcoef_PC,[-1 1]); colorbar;set(gca,'YDir','normal')
            else
               j
               phi_curr = phi(:,j,prepostchronic);
               plot_corr_surf(corrcoef_PC,phi_curr)
                
            end
        end
        corrcoef_arr = cat(3,corrcoef_arr, corrcoef_PC);
    end

    if prepostchronic == 0
        phi_stage = mean(phi,3);
    else
        phi_stage = phi(:,:,prepostchronic);
    end
    phi_stage1 = permute(phi_stage,[3 1 2]); phi_stage1 = repmat(phi_stage1,size(phi,1),1);
    phi_stage2 = permute(phi_stage,[1 3 2]); phi_stage2 = repmat(phi_stage2,1,size(phi,1));
    phidiffs = phi_stage1 - phi_stage2;     % Position i,j corresponds to freqband i - freqband j

    ut = logical(repmat(triu(ones(size(phi,1),size(phi,1)),1),[1 1 size(phi,2)]));  % Generate indices to pull out upper triangle
    phidiffs = phidiffs(ut);
    corrcoef_arr = corrcoef_arr(ut);
    figure;plot(phidiffs,corrcoef_arr,'.'); xlabel('Phase difference (days)'); ylabel('Corr coef');

end



function plot_corr_surf(corrcoef_PC,phi_curr)

    index = ~isnan(phi_curr);
    phi_curr = phi_curr(index);
    corrcoef_PC = corrcoef_PC(index,:);
    corrcoef_PC = corrcoef_PC(:,index);

    [phi1 phi2] = meshgrid(phi_curr,phi_curr);
    phi1=phi1(:);phi2=phi2(:);corrcoef_PC=corrcoef_PC(:);


    min_phi=0.2; step_phi = 0.01; max_phi = 1.2;
    min_phi=min(phi_curr); step_phi = 0.01; max_phi = max(phi_curr);
    [phi1I phi2I] = meshgrid(min_phi:step_phi:max_phi,min_phi:step_phi:max_phi);
    corrcoef_PC_i = griddata(phi1,phi2,corrcoef_PC,phi1I,phi2I);
    
    % 3D mesh plot
%     figure; mesh(phi1I,phi2I,corrcoef_PC_i);
%     hold on; plot3(phi1,phi2,corrcoef_PC,'.')
    
    % 2D colour plot
    figure; imagesc([min_phi max_phi],[min_phi max_phi],corrcoef_PC_i,[-1 1]); colorbar;set(gca,'YDir','normal')

    % 2D plot with data points for each phi pair, colour representing
    % corrrelation coefficient
%     figure;
%     colourarr = autumn;
%     corrcoef_PC_2_colorindices = round((corrcoef_PC+1)*(size(colourarr,1)/2 - 0.5) + 1);
%     for i = 1:(size(corrcoef_PC,1)*size(corrcoef_PC,2))
%     hold on; plot(phi1(i),phi2(i),'x','Color',colourarr(corrcoef_PC_2_colorindices(i),:),'MarkerSize',20,'LineWidth',2)
%     end

end