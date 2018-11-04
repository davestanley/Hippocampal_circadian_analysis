
function func_plot_correlation_basic(phi,Amp)
        corr_within_rats = 1;
        % Should have amp normlization turned on    
        plot_individual_correlations = 0;
        
        if corr_within_rats
            ts_latent = [];
            ts_SS = [];
            corr_all = zeros(size(Amp,2),2);
            pcorr_all = zeros(size(Amp,2),2);
            corrLB = zeros(size(Amp,2),2);
            corrUB = zeros(size(Amp,2),2);
            
            for j = 1:size(Amp,2)
                temp2 = Amp(:,j,2);
                temp3 = phi(:,j,2);
                index = ~isnan(temp2); temp2 = temp2(index); temp3 = temp3(index);      % Remove NaNs
                [temp ptemp lbtemp ubtemp] =  corrcoef(temp2,temp3);
                corr_all(j,1) = temp(1,2); pcorr_all(j,1) = ptemp(1,2); corrLB(j,1) = lbtemp(1,2); corrUB(j,1) = ubtemp(1,2);
                ts_latent = [ts_latent; [temp2 temp3]];
                if plot_individual_correlations
                    figure; subplot(211); hold on;
%                     plot(temp2,temp3,'.');
                    plot(temp2,'b'); plot(temp3,'r');
                end

                temp2 = Amp(:,j,3);
                temp3 = phi(:,j,3);
                index = ~isnan(temp2); temp2 = temp2(index); temp3 = temp3(index);      % Remove NaNs
                [temp ptemp ] =  corrcoef(temp2,temp3);
                corr_all(j,2) = temp(1,2); pcorr_all(j,2) = ptemp(1,2); corrLB(j,2) = lbtemp(1,2); corrUB(j,2) = ubtemp(1,2);
                
                ts_SS = [ts_SS; [temp2(:) temp3(:)]];
                if plot_individual_correlations
                    subplot(212); hold on;
    %                 plot(temp2,temp3,'.')
                    plot(temp2,'b'); plot(temp3,'r');
                end
            end
            
        else
            ts_latent = [];
            ts_SS = [];
            corr_all = zeros(size(Amp,1),2);
            for i = 1:size(Amp,1)

                temp2 = Amp(i,:,2);
                temp3 = phi(i,:,2);
                index = ~isnan(temp2); temp2 = temp2(index); temp3 = temp3(index);      % Remove NaNs
                temp =  corrcoef(temp2,temp3);
                corr_all(i,1) = temp(1,2);
                ts_latent = [ts_latent; [temp2(:) temp3(:)]];

                temp2 = Amp(i,:,3);
                temp3 = phi(i,:,3);
                index = ~isnan(temp2); temp2 = temp2(index); temp3 = temp3(index);      % Remove NaNs
                temp =  corrcoef(temp2,temp3);
                corr_all(i,2) = temp(1,2);
                ts_SS = [ts_SS; [temp2(:) temp3(:)]];
            end
            
        end
        
        figure('Color','w'); bar(mean(corr_all),'EdgeColor',[0.5 0.5 0.5]); colormap ([0.5 0.5 0.5])
        hold on; errorbar(1:2,mean(corr_all),confidence(corr_all,0.05),'k','LineStyle','none','LineWidth',2);
        hold on; plot(repmat([1 2], size(corr_all,1),1),corr_all,'kx','LineWidth',2,'MarkerSize',20)
        set(gca,'FontSize',30,'XTickLabel',{},'Box','off')
        
end