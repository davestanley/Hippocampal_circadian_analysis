
function corrcoef_3D = func_plot_corrcoef_EEG_bands(r_seiz,r_ns,r_ctrl,analyze_only_seizing,prepostchronic,smoothmode_ts,phi,rat2plot)

    plot_individual_rats = 1;
    plot_on = 0;
    show_significance = 1;

    if analyze_only_seizing == 1
        [T X] = rat_2_ts (r_seiz(rat2plot),prepostchronic, smoothmode_ts);
    elseif analyze_only_seizing == 2
        [T X] = rat_2_ts (r_ctrl,prepostchronic, smoothmode_ts);
    else
        [T X] = rat_2_ts (r_ns,prepostchronic, smoothmode_ts);
    end
    
    if prepostchronic == 0
        phi_stage = mean(phi,3);
    else
        phi_stage = phi(:,:,prepostchronic);
    end
    
    
    corrcoef_EEG_bands(T,X,phi_stage,plot_on,plot_individual_rats, show_significance)
    
end