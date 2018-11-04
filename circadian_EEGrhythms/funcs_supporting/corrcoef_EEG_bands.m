
function [corrcoef_3D, hi] = corrcoef_EEG_bands(T,X,phi_stage,plot_on,plot_individual_rats, show_significance)

    %plot_individual_rats = 0;
    %plot_on = 1;
    
    % phi is not necessary if not doing any plotting...
    
        plot_vs_phi = 0; % if 1, have the x and y axes of the colour plots be phase. Otherwise, they're freq band.
        
    ratrange = 1:length(X);
    corrcoef_arr = [];
    P_arr = [];
    
    for j = ratrange
        %corrcoef_PC = corrcoef_mat(X{j});
        [corrcoef_PC,P] = corrcoef(X{j});
        corrcoef_arr = cat(3,corrcoef_arr, corrcoef_PC);
        P_arr = cat(3,P_arr,P);
    end

    if plot_individual_rats
        for j = ratrange
            phi_curr = phi_stage(:,j);
            hi(j) = plot_corrcoefs(corrcoef_arr(:,:,j),P_arr(:,:,j),plot_vs_phi,phi_curr,show_significance);
            xlabel('Freq band'); ylabel('Freq band');
        end
    end
    
    corrcoef_3D = corrcoef_arr;
       
    
    
    
    if plot_on
        
        Nbands = size(corrcoef_3D,1);
        Nrats = length(X);
        phi_stage1 = permute(phi_stage,[3 1 2]); phi_stage1 = repmat(phi_stage1,Nbands,1);
        phi_stage2 = permute(phi_stage,[1 3 2]); phi_stage2 = repmat(phi_stage2,1,Nbands);
        phidiffs = phi_stage1 - phi_stage2;     % Position i,j corresponds to freqband i - freqband j

        ut = logical(repmat(triu(ones(Nbands,Nbands),1),[1 1 Nrats]));  % Generate indices to pull out upper triangle
        phidiffs = phidiffs(ut);
        corrcoef_arr = corrcoef_arr(ut);
    
        figure;plot(phidiffs,corrcoef_arr,'.'); xlabel('Phase difference (days)'); ylabel('Corr coef');
    end

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


function hi = plot_corrcoefs(mycorrcoef,myP,plot_vs_phi,phi_curr,show_significance)
    %show_significance = 1;

    if ~plot_vs_phi
        figure;
        colourmapval = [summer; flipud(autumn)]; colormap(colourmapval);
        %colormap winter;
        imagesc(mycorrcoef,[-1 1]); colorbar;set(gca,'YDir','normal')
        hi = gca;
        
        if show_significance
            myh = myP < 0.05;
            myh(logical(eye(size(myh)))) = 1;
            h = findall(gca,'Type','image');
            set (h(1), 'AlphaData', myh);                 % Set transparency to account for error.
        end
        
    else
       plot_corr_surf(mycorrcoef,phi_curr)
    end

end

